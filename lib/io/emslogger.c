/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/io/emslogger.c
 *
 *  \brief Message logging with timestamp
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: emslogger.c 5833 2018-06-27 00:21:35Z riz008 $
 */

#include <stdlib.h>
#include <stdio.h>
#if defined(__sparc)
#include <strings.h>
#else
#include <string.h>
#endif
#include <ctype.h>
#include <time.h>
#include "ems.h"

#ifndef DEFAULTTAGS
#define DEFAULTTAGS
#endif

#define MAXNUMMAND  10
#define MAXNUMGRIDS 20
#define NUMTAGS   7
#define MAXCUSTOMTAGS   10

#define ZEROPAD 1
#define SPECIAL  32          /* 0x */
#define SPACE        8       /* space if plus */
#define PLUS 4               /* show plus */
#define LEFT 16              /* left justified */
#define LARGE        64      /* use 'ABCDEF' instead of 'abcdef' */
#define SIGN 2               /* unsigned/signed long */

/* # define ALLOW_CUSTOMLOG
*/

/* Local variables.
 */

typedef struct {
  char* tag;                    /* The debug tag name */
  char* mandatory[MAXNUMMAND];  /* The list of mandatory tags */
  int  id;
  char* nice;
} log_tag_map_t;


static int ndtags = 0;


static char *dtags[NUMTAGS];
static int *itags[NUMTAGS];
#if defined(ALLOW_CUSTOMLOG)
static int cndtags;
static log_tag_map_t *custom_map = NULL;
#endif

FILE *log_fp = NULL;
long  log_init_offset = 0L; // This is used as the end of init offset
			    // for the runlog wrap

static log_tag_map_t dtags_map[NUMTAGS] = {
  {"main", {NULL},LMAIN,"main  "},
  {"debug", {NULL},LDEBUG,"debug "},
  {"metric", {NULL},LMETRIC,"metric"},
  {"trace", {NULL},LTRACE,"trace "},
  {"info", {NULL},LINFO,"info  "},
  {"error",{NULL},LERROR,"error "},
  {"warn",{NULL},LWARN,"warn  "}
};

/* Run though the setting command and split out the tags into
 * a large NULL terminated array.
 */
static char **split_tags(const char *command, char delim)
{
  char **tags = NULL;
  int i;
  int n = 1;                    /* Number of commands (one more than
                                   delim) */
  int start_pos = 0;
  char buffer[MAXSTRLEN];

/* Count the number of delimitars to determine how long the return array
 * will need to be.
 */
  for (i = 0; i < (int)strlen(command); ++i)
    if (command[i] == delim)
      ++n;

  tags = (char **)malloc(sizeof(char *) * (n + 1));
  memset(tags, 0, sizeof(char *) * (n + 1));

/* Extract the tags.
 */
  n = 0;
  for (i = 0; i <= (int)strlen(command); ++i) {
    if ((command[i] == delim) || (command[i] == '\000')) {
      if (start_pos != i) {
        strncpy(buffer, &command[start_pos], i - start_pos);
        buffer[i - start_pos] = 0;
        tags[n] = strdup(buffer);
        start_pos = i + 1;
        ++n;
      }
    }
  }
  return tags;
}


/* Insert a tag in the dtags array by checking the mapping array.
 * If the tag is located then insert it in the dtags table, and
 * any of it's mandatory tags.
 */
void enable_loglevel(char* tag)
{
  int i;
  log_tag_map_t *map = NULL;

  if(log_if(tag))
    return;

/* Locate the tag in the mapping list.
 */
  for (i = 0; i < NUMTAGS; ++i) {
    if (strcmp(tag, dtags_map[i].tag) == 0) {
      map = &dtags_map[i];
      break;
    }
  }

  if (map == NULL) {
#if defined(ALLOW_CUSTOMLOG)
  int ll = atoi(tag);
  if(ll >= LNONE)
    fprintf(stderr, "Custom log levels must be smaller %d \n", LNONE);
  else
  {
     if(custom_map == NULL)
     {
      custom_map = (log_tag_map_t *)calloc(MAXCUSTOMTAGS,sizeof(log_tag_map_t));
      cndtags =0;
     }

     //map = malloc(log_tag_map_t);
     sprintf(custom_map[cndtags].tag,"CUSTOM%d",ll);
     custom_map[cndtags].id=ll;
     cndtags++;
     return;
  }
#else
    fprintf(stderr, "The log tag '%s' is not supported.\n", tag);
    /*log_usage();*/
#endif
  }

/* Now insert the tag (if not already in the list), and all
 * of it's mandatory tags.
 */
  dtags[ndtags] = map->tag;
  itags[ndtags] = &map->id;
  ++ndtags;
  for (i = 0; map->mandatory[i] != NULL; ++i)
    enable_loglevel(map->mandatory[i]);
}


void disable_loglevel(char* tag)
{
  int i;
  log_tag_map_t *map = NULL;

  if(!log_if(tag))
    return;

/* Locate the index of the tag in the mapping list.
 */
  for (i = 0; i < NUMTAGS; ++i) {
    if (strcmp(tag, dtags_map[i].tag) == 0) {
      map = &dtags_map[i];
      break;
    }
  }

  if (map == NULL)
    return;
#if defined(ALLOW_CUSTOMLOG)
  if(custom_map != NULL)
  {
    i=0;
    for(;i<cndtags;i++)
    {
      if(strcmp(tag, custom_map[i].tag) == 0)
      {
        for (; i < cndtags-1; ++i)
          custom_map[i] = custom_map[i+1];
        cndtags--;
        return;
      }
    }
  }
#endif



  /* Now move the tags forward
   */
  for (i = 0; i < NUMTAGS-1; ++i)
  {
    dtags[i] = dtags[i+1];
    itags[i] = itags[i+1];
  }
  dtags[ndtags] = NULL;
  itags[ndtags] = NULL;
  ndtags--;
}



/* Set the log flag, and parses the log settings
 * to determine which log flags will be set.
 */
void log_init(const char *dcommand)
{
  int i, j;
  char **tags;

  if (log_fp == NULL)
    log_fp = stdout;

  if(dcommand == NULL)
    dcommand="warn"; /* ,info,main" */

  tags = split_tags(dcommand, ',');

  for (i = 0; tags[i] != NULL; ++i) {

    if (strcmp(tags[i], "all") == 0) {
      ndtags = NUMTAGS;
      for (j = 0; j < NUMTAGS; ++j)
      {
        dtags[j] = dtags_map[j].tag;
        itags[j] = &dtags_map[j].id;
      }
      break;
    } else
      enable_loglevel(tags[i]);

    free(tags[i]);            /* Release the current tags memory */
  }

  free(tags);                 /* Release the tags array */
#if defined(DEFAULTTAGS)
/*
  enable_loglevel("main");
  enable_loglevel("info");*/
  enable_loglevel("warn");

#endif
}

void log_set_output(const char *filename)
{

  if (filename == NULL || strcmp(filename, "-") == 0 )
    log_fp = stdout;
  else
    log_fp = fopen(filename, "w");
}



/* Cleanup the log settings.
 */
void log_end(void)
{
  if(log_fp != NULL)
  {
    fflush(log_fp);
    if (log_fp != stderr && log_fp != stdout)
      fclose(log_fp);

    log_fp = NULL;
  }
}


/* If tag is an enabled command the return true (1) else false (0).
 */
int log_if(const char *tag)
{
 /* Run though the list looking for a match.
  */
  int i;
  for (i = 0; i < ndtags; ++i) {
    if (strcmp(dtags[i], tag) == 0)
      return 1;
  }
  return 0;
}


int is_log_fatal(int level)
{
  if(level < 0)
    return 1;
  return 0;
}


/* If log level is an enabled command the return true (1) else false (0).
 */
int is_log_enabled(int level)
{
  int i=0;
  if(is_log_fatal(level))
    return 1;

  for(;i<ndtags;i++)
  {
    if(*itags[i] == level)
      return 1;
  }
  return 0;
}


/* Flush the logging output stream
 */
void emslog_flush()
{
  if(log_fp != NULL)
    fflush(log_fp);
}


/* Diagnostic usage.
 */
void log_usage(void)
{
  int i, j;

  fprintf(stderr, "   Supported log tags are:\n");
  fprintf(stderr, "      - all\n");
  for (i = 0; i < NUMTAGS; ++i) {
    fprintf(stderr, "      - %s", dtags_map[i].tag);
    if (dtags_map[i].mandatory[0] != NULL) {
      int pad = 20 - strlen(dtags_map[i].tag);
      for (j = 0; j < pad; ++j)
        fprintf(stderr, " ");
      fprintf(stderr, " [");
      for (j = 0; dtags_map[i].mandatory[j] != NULL; ++j)
        fprintf(stderr, " %s", dtags_map[i].mandatory[j]);
      fprintf(stderr, " ]");
    }
    fprintf(stderr, "\n");
  }
}


FILE* emslog_fp(int level)
{
  if(is_log_enabled(level))
    return  log_fp;
  return NULL;
}

/*
 * Prints a marker and then stores the offset so that we wrap back to
 * this point in the log file
 */
void emslog_set_log_offset(const char str[])
{
  if(log_fp != NULL) {
    fprintf(log_fp, "<----------  %s  ----------->\n", str);
    fflush(log_fp);
    log_init_offset = ftell(log_fp);
  }
}

char* level2tag(int t)
{
  int i=0;
  if(t == LERROR)
    return "ERROR";
  if(t == LFATAL)
    return "FATAL";
  if(t == LPANIC)
    return "PANIC";

  for(;i< NUMTAGS;i++)
  {
    if(dtags_map[i].id == t)
      return dtags_map[i].tag;
  }
#if defined(ALLOW_CUSTOMLOG)
  if(custom_map != NULL)
  {
    i=0;
    for(;i<cndtags;i++)
    {
      if(custom_map[i].id == t)
        return  custom_map[i].tag;
    }
  }
#endif
  return "NA";
}

/* Nice converts tags in string s of equal length
 * currently 6 characters are expected
 *
 */
char* level2nice(int t)
{
  int i=0;
  if(t == LERROR)
    return "ERROR ";
  if(t == LFATAL)
    return "FATAL ";
  if(t == LPANIC)
    return "PANIC ";

  for(;i< NUMTAGS;i++)
  {
    if(dtags_map[i].id == t)
      return dtags_map[i].nice;
  }
#if defined(ALLOW_CUSTOMLOG)
  if(custom_map != NULL)
  {
    i=0;
    for(;i<cndtags;i++)
    {
      if(custom_map[i].id == t)
        return  custom_map[i].nice;
    }
  }
#endif
  return "NA";
}





/* Write the log output to the standard out or a log file (if requested).
 *
 * param: level - the log level
 * param: tag   - the origin of the log request
 * param: s     - the log message
 * param: args  - the initialised variable parametere list
 */
static int do_writelog(int level, const char *tag, const char *s,va_list args)
{
    char buffer[MAXSTRLEN];
    time_t now = time(NULL);
    struct tm *t = localtime(&now);
    /*
     * If we're over the file size limit, then wrap
     */
    if (log_fp != NULL && log_init_offset > 0L) {
      // hd_init sets a marker so wait until then
      char *flim_s = getenv("EMS_RUNLOG_FILE_LIMIT");
      int   flim   = (flim_s != NULL ? atoi(flim_s) : 10e6); // 10 MB
                                                             // as default
      // This the size in bytes
      long fs = ftell(log_fp);
      if (fs > flim)
	fseek(log_fp, log_init_offset, SEEK_SET);
    }

    vsprintf(buffer, s, args);
    strip(buffer, "\n");
    if (strlen(buffer) > 0) {
      if(log_fp == NULL)
      {
          fprintf(stderr, "PANIC-OUT [%d/%2.2d/%2.2d %2.2d:%2.2d:%2.2d]-[%s](%s) %s\n",
              1900 + t->tm_year, t->tm_mon, t->tm_mday,
              t->tm_hour, t->tm_min, t->tm_sec,
              level2nice(level),tag, buffer);
      } else {
        fprintf(log_fp, "[%d/%2.2d/%2.2d %2.2d:%2.2d:%2.2d]-[%s](%s) %s\n",
              1900 + t->tm_year, t->tm_mon, t->tm_mday,
              t->tm_hour, t->tm_min, t->tm_sec,
              level2nice(level),tag, buffer);
        emslog_flush();

        /* if this is a severe log entry, duplicate it to stderr */
        if(is_log_fatal(level))
        {
          fprintf(stderr, "[%d/%2.2d/%2.2d %2.2d:%2.2d:%2.2d]-[%s](%s) %s\n",
              1900 + t->tm_year, t->tm_mon, t->tm_mday,
              t->tm_hour, t->tm_min, t->tm_sec,
              level2nice(level),tag, buffer);
        }
      }
    }
    return 1;
}

/* Write the log output to the standard out or a log file (if requested).
 *
 * param: level - the log level
 * param: tag   - the origin of the log request
 * param: s     - the log message
 * param: args  - the initialised variable parametere list
 */
int writelog(int level, const char *tag, const char *s,va_list args)
{
  /*provide a 'write direct' fo module depending log level checking
  if (!is_log_enabled(level))
    return (0);*/
  return do_writelog(level,tag,s,args);
}

/* Write the log output to the standard out or a log file (if requested).
 *
 * param: level - the log level
 * param: format- the log message and formatting
 * param: ...  - the initialised variable parametere list
 */
int errfn_log_default(int level, char *format, ...)
{
  int i = 0;
  va_list args;
    /* Check that the log flag is set here to not waste extra cycles.*/
  if (!is_log_enabled(level))
    return (0);

  va_start(args, format);
  i =do_writelog(level,"",format,args);
  va_end(args);
  return i;
}


/* Write the log output to the standard out or a log file (if requested).
 *
 * param: level - the log level
 * param: tag   - the origin of the log request
 * param: format- the log message and formatting
 * param: ...  - the initialised variable parametere list
 */
int errfn_tag_default(int level, char *tag, char *format, ...)
{
  int i = 0;
  va_list args;
    /* Check that the log flag is set here to not waste extra cycles.*/
  if (!is_log_enabled(level))
    return (0);

  va_start(args, format);
  i =do_writelog(level,tag,format,args);
  va_end(args);
  return i;
}

logfn   emslog    = errfn_log_default;
logtag  emstag    = errfn_tag_default;

/*
void errfn_trace_default(char *format, ...)
{
  va_list args;
  va_start(args, format);
  dlog(LTRACE,"",format,args);
  va_end(args);
}


void errfn_metric_default(char *format, ...)
{
  va_list args;
  va_start(args, format);
  dlog(LMETRIC,"",format,args);
  va_end(args);
}


void errfn_main_default(char *format, ...)
{
  va_list args;
  va_start(args, format);
  dlog(LMAIN,"",format,args);
  va_end(args);
}


void errfn_info_default(char *format, ...)
{
  va_list args;
  va_start(args, format);
  dlog(LINFO,"",format,args);
  va_end(args);
}


void errfn_debug_default(char *format, ...)
{
  va_list args;
  va_start(args, format);
  dlog(LDEBUG,"",format,args);
  va_end(args);
}
*/


