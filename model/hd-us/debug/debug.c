/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/debug/debug.c
 *  
 *  Description:
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: debug.c 5873 2018-07-06 07:23:48Z riz008 $
 *
 */

#include <stdarg.h>
#include <stdlib.h>
#if defined(__sparc)
#include <strings.h>
#else
#include <string.h>
#endif
#include <time.h>
#include "hd.h"

#define MAXNUMMAND	10
#define MAXNUMGRIDS	20
#define NUMTAGS 	8

/* List of debug tags and their associated tags.
 */
typedef struct {
  char *tag;                    /* The debug tag name */
  char *mandatory[MAXNUMMAND];  /* The list of mandatory tags */
} debug_tag_map_t;

static debug_tag_map_t dtags_map[NUMTAGS] = {
  {"conversions", {NULL}},
  {"time", {NULL}},
  {"dump", {"time", NULL}},
  {"particles", {"time", NULL}},
  {"init_m", {NULL}},
  {"init_w", {NULL}},
  {"ecology", {NULL}},
  {"sediments", {NULL}}
};


/* Local variables.
 */
static int ndtags;
static char *dtags[NUMTAGS];
static FILE *debug_fp = NULL;


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
      } else
        hd_quit("split_tags: Malformed diagnostic tag.\n");
    }
  }

  return tags;
}


/* Insert a tag in the dtags array by checking the mapping array.
 * If the tag is located then insert it in the dtags table, and
 * any of it's mandatory tags.
 */
static void insert_tag(char *tag)
{
  int i;
  debug_tag_map_t *map = NULL;

/* Locate the tag in the mapping list.
 */
  for (i = 0; i < NUMTAGS; ++i) {
    if (strcmp(tag, dtags_map[i].tag) == 0) {
      map = &dtags_map[i];
      break;
    }
  }

  if (map == NULL) {
    fprintf(stderr, "The debug tag '%s' is not supported.\n", tag);
    debug_usage();
    hd_quit("\n");
  }

/* Now insert the tag (if not already in the list), and all
 * of it's mandatory tags.
 */
  dtags[ndtags] = map->tag;
  ++ndtags;
  for (i = 0; map->mandatory[i] != NULL; ++i)
    insert_tag(map->mandatory[i]);
}


/* Set the debug flag, and parses the debug settings
 * to determine which debug flags will be set.
 */
void debug_init(const char *dcommand)
{
  int i, j;

  debug = 1;
  if (debug_fp == NULL)
    debug_fp = stdout;
  ndtags = 0;
  if (debug) {
    char **tags = split_tags(dcommand, ',');

    for (i = 0; tags[i] != NULL; ++i) {

      if (strcmp(tags[i], "all") == 0) {
        ndtags = NUMTAGS;
        for (j = 0; j < NUMTAGS; ++j)
          dtags[j] = dtags_map[j].tag;
      } else
        insert_tag(tags[i]);

      free(tags[i]);            /* Release the current tags memory */
    }

    free(tags);                 /* Release the tags array */

    if (ndtags <= 0) {
      hd_warn
        ("debug_init: No tags were located, turning debug mode off.\n");
      debug = 0;
    }else
      enable_loglevel("debug");
  }
}

void debug_set_output(const char *filename)
{
/*
  if (strcmp(filename, "-") == 0)
    debug_fp = stdout;
  else
    debug_fp = fopen(filename, "w");
*/
}



/* Cleanup the debug settings.
 */
void debug_end(void)
{
  /* UR taken care of in lib
  if (debug) {
    fflush(debug_fp);
    if (debug_fp != stderr)
      fclose(debug_fp);
  }
  debug = 0;
  debug_fp = NULL;
  */
}


/* If dtag is an enabled command the return true (1) else false (0).
 */
int debug_if(const char *dtag)
{

  if (debug) {
/* Run though the list looking for a match.
 */
    int i;
    for (i = 0; i < ndtags; ++i) {
      if (strcmp(dtags[i], dtag) == 0)
        return 1;
    }
  }

  return 0;
}


/* Diagnostic usage.
 */
void debug_usage(void)
{
  int i, j;

  fprintf(stderr, "   Supported debug tags are:\n");
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


/* Write the debug output to the standard error or a log file (if requested).
 */
void dlog(const char *tag, const char *s, ...)
{
  va_list args;

/* Check that the debug flag is set.
 * /
  if (debug) {
    char buffer[MAXSTRLEN];
   time_t now = time(NULL);
    struct tm *t = localtime(&now);
*/
    va_start(args, s);
/*
    vsprintf(buffer, s, args);
  
    strip(buffer, "\n");
    if (strlen(buffer) > 0) {
      fprintf(debug_fp, "%s|%d/%2.2d/%2.2d|%2.2d:%2.2d:%2.2d|%s\n", tag,
              1900 + t->tm_year, t->tm_mon, t->tm_mday,
              t->tm_hour, t->tm_min, t->tm_sec, buffer);
      fflush(debug_fp);
    }
    
    UR doubled .... 
    emstag(LDEBUG,tag,buffer);
    */
    writelog(LDEBUG,tag,s,args);
    fflush(debug_fp);

    va_end(args);
    /*
  } else
    hd_quit
      ("dlog: Request to output debug information without debug being set.\n");
  */
}
