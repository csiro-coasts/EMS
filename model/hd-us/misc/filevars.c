/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/misc/filevars.c
 *  
 *  Description:
 *  Function to help in extracting the variable
 *  names appended to the filename.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: filevars.c 5873 2018-07-06 07:23:48Z riz008 $
 *
 */

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "hd.h"

static int decode_args(char *oprm, char *args[], int maxn);
static int find_argument(char *tag, int nargs, char *args[]);
static char *get_argument_value(char *tag, int nargs, char *args[],
                                int req);

struct {
  char *tag;
  char *name;
} fv_default_names[] = {
  {
  "air_temp", "air_temp"}, {
  "cloud", "cloud"}, {
  "direction", "direction"}, {
  "eta", "eta"}, {
  "evaporation", "evaporation"}, {
  "humidity", "humidity"}, {
  "period", "period"}, {
  "pressure", "pressure"}, {
  "precipitation", "precipitation"}, {
  "swr", "swr"}, {
  "ub", "ub"}, {
  "wind_u", "u"}, {
  "wind_v", "v"}, {
  "u", "u"}, {
  "u1", "u1"}, {
  "v", "v"}, {
  "u2", "u2"}, {
NULL, NULL}};

char *fv_get_filename(char *oprm, char *buf)
{
  char *args[256];
  int nargs;
  strcpy(buf, oprm);
  nargs = decode_args(buf, args, 256);
  return args[0];
}

char *fv_get_varname(char *oprm, char *tag, char *buf)
{
  char *args[256];
  int nargs;
  int i = 0;
  strcpy(buf, oprm);
  nargs = decode_args(buf, args, 256);

  if (nargs > 1) {
    char *v = NULL;
    args[0][0] = '\000';        /* Blank the filename */
    v = get_argument_value(tag, nargs, args, 0);
    if (v != NULL)
      return v;
  }

  /* Check in default list */
  while (fv_default_names[i].tag != NULL) {
    if (strcasecmp(fv_default_names[i].tag, tag) == 0)
      return fv_default_names[i].name;
    ++i;
  }

  return tag;
}


static int decode_args(char *p, char *args[], int maxn)
{
  int len = strlen(p);
  int i;

  for (i = 0; i < len; ++i) {
    if (p[i] == '(' || p[i] == ')' || p[i] == ',')
      p[i] = ' ';
  }

  return parseline(p, args, maxn);
}

/* Search for the specified argument within the list, and return
 * the value as a string. NULL if it cannot be found.
 */
static char *get_argument_value(char *tag, int nargs, char *args[],
                                int req)
{
  int index = find_argument(tag, nargs, args);

  if (index >= 0) {
    char *arg = args[index];
    while (*arg) {
      if (*arg++ == '=')
        return arg;
    }
  }

  if (req)
    quit("Unable to locate argument '%s'.\n", tag);

  return NULL;
}

/* Locate an argument in the list by name.
 * @return index of tag in list or -1.
 */
static int find_argument(char *tag, int nargs, char *args[])
{
  int i;

  for (i = 0; i < nargs; ++i) {
    if (strncasecmp(tag, args[i], strlen(tag)) == 0)
      return i;
  }

  return -1;
}
