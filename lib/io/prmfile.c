/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/io/prmfile.c
 *
 *  \brief Parse ASCII file using keywords
 *
 *  Routines to handle ascii keyword files. These files are commonly
 *  used read parameter "key" or "token" pairs in various programs
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: prmfile.c 5833 2018-06-27 00:21:35Z riz008 $
 */

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include "ems.h"
#include "errno.h"
#include "ems.h"


errfn keyprm_errfn = errfn_quit_default;

int keyprm_verbose = 0;
static int keyprm_case_sensitive = 1;

#define STRCMP(a,b) ((keyprm_case_sensitive) ? strcmp(a,b) : strcasecmp(a,b))
#define STRNCMP(a,b,l) ((keyprm_case_sensitive) ? strncmp(a,b,l) : strncasecmp(a,b,l))

/** Set the behavior of the key finding routines
  * if a key is not found.
  *
  * @param fn pointer to error handling function.
  */
void prm_set_errfn(void (*fn) (const char *, ...))
{
  if (fn)
    keyprm_errfn = fn;
  else
    keyprm_errfn = errfn_quiet;
}

/** Set the case sensitivity for key finding routines.
  *
  * @param c 1 if case sensitive, else 0.
  */
void prm_set_case(int c)
{
  keyprm_case_sensitive = c;
}


/** Skip forward from the current file position to
  * the start of the next line beginning with key.
  *
  * @param fp pointer to stdio FILE structure.
  * @param key keyname to locate in file.
  * @return non-zero if successful.
  */
int prm_skip_to_start_of_key(FILE * fp, char *key)
{
  char buf[MAXLINELEN];
  int len = strlen(key);
  long fpos;
  char *s;
  int rewound = 0;
  /*UR 22/05/2008
   * reduced failure of finding a key to
   * 'TRACE'
   * since it is not vital and allowed.
   */
  do {
    fpos = ftell(fp);
    s = fgets(buf, MAXLINELEN, fp);
    if (s == NULL) {
      if (!rewound) {
        fpos = 0L;
        if (fseek(fp, fpos, 0))
          quit("prm_skip_to_start_of_key: %s\n", strerror(errno));
        s = fgets(buf, MAXLINELEN, fp);
        rewound = 1;
      } else {
        emstag(LTRACE,"lib:prmfile:prm_skip_to_start_of_key","key %s not found",key);
       /* (*keyprm_errfn) ("prm_skip_to_start_of_key: key %s not found\n",
                         key);*/
        return (0);
      }
    }

    while (is_blank(s[0]))
      s++;

    /* Truncate the string at the first space after the key length. */
    if (is_blank(s[len]))
      s[len] = '\000';
  } while (STRCMP(key, s) != 0);

  if (s == NULL) {
    emstag(LTRACE,"lib:prmfile:prm_skip_to_start_of_key","key %s not found", key);
    /*(*keyprm_errfn) ("prm_skip_to_start_of_key: key %s not found\n", key);*/
    return (0);
  }

  /* seek to start of line */
  if (fseek(fp, fpos + (s - buf), 0))
    quit("prm_skip_to_start_of_key: %s\n", strerror(errno));

  return (1);
}



/** Skip forward from the current file position to
  * the next line beginning with key, positioned at the character
  * immediately after the key.
  *
  * @param fp pointer to stdio FILE structure.
  * @param key keyname to locate in file.
  * @return non-zero if successful.
  */
int prm_skip_to_end_of_key(FILE * fp, char *key)
{
  char buf[MAXLINELEN];
  int len = strlen(key);
  long fpos;
  char *s;
  int rewound = 0;
  /*UR 10/06/2005
   * reduced failure of finding a key to 
   * 'TRACE'
   * since it is not vital and allowed. 
   */
  do {
    fpos = ftell(fp);
    s = fgets(buf, MAXLINELEN, fp);
    if (s == NULL) {
      if (!rewound) {
        fpos = 0L;
        if (fseek(fp, fpos, 0))
          quit("prm_skip_to_end_of_key: %s\n", strerror(errno));
        s = fgets(buf, MAXLINELEN, fp);
        rewound = 1;
      } else {
        emstag(LTRACE,"lib:prmfile:prm_skip_to_end_of_key"," key %s not found\n",
                         key);
        /*
        (*keyprm_errfn) ("prm_skip_to_end_of_key: key %s not found\n",
                         key);*/
        return (0);
      }
    }

    if (s == NULL)
      break;

    while (is_blank(s[0]))
      s++;

    /* Truncate the string at the first space after the key length. */
    if (strlen(s) > len && is_blank(s[len]))/*UR added length check to prevent invalid write*/
      s[len] = '\000';
  } while (STRCMP(key, s) != 0);

  if (s == NULL) {
    emstag(LTRACE,"lib:prmfile:prm_skip_to_end_of_key","key %s not found\n", key);
    /*(*keyprm_errfn) ("prm_skip_to_end_of_key: key %s not found\n", key);*/
    return (0);
  }

  /* seek to character after key */
  if (fseek(fp, fpos + (s - buf) + len, 0))
    quit("prm_skip_to_end_of_key: %s\n", strerror(errno));

  return (1);
}


/** Read an integer value following a key.
  *
  * @param fp pointer to stdio FILE structure.
  * @param key keyname to locate in file.
  * @param p pointer to returned integer value.
  * @return non-zero if successful.
  */
int prm_read_int(FILE * fp, char *key, int *p)
{
  /* Read parameter */
  if (!prm_skip_to_end_of_key(fp, key))
    return (0);

  if (fscanf(fp, "%d", p) != 1)
    quit("prm_read_int: Can't read %s\n", key);

  /* Echo to stdout if verbose enough */
  emstag(LTRACE,"lib:prmfile:prm_read_int","%s  %d", key, *p);

  return (1);
}

/** Read an double value following a key.
  *
  * @param fp pointer to stdio FILE structure.
  * @param key keyname to locate in file.
  * @param p pointer to returned double value.
  * @return non-zero if successful.
  */
int prm_read_double(FILE * fp, char *key, double *p)
{
  /* Read parameter */
  if (!prm_skip_to_end_of_key(fp, key))
    return (0);

  if (fscanf(fp, "%lf", p) != 1)
    quit("prm_read_double: Can't read %s\n", key);

  /* Echo to stdout if verbose enough */
  emstag(LTRACE,"lib:prmfile:prm_read_double","%s  %f", key, *p);

  return (1);
}


/** Read an array of double value following a key.
  *
  * @param fp pointer to stdio FILE structure.
  * @param key keyname to locate in file.
  * @param p pointer to returned array of double value.
  * @param size pointer to returned array size.
  * @return non-zero if successful.
  */
int prm_read_darray(FILE * fp, char *key, double **p, int *size)
{
  int i = 0;

  /* Read parameter */
  if (!prm_skip_to_end_of_key(fp, key))
    return (0);

  if (fscanf(fp, "%d", size) != 1 || *size < 0) {
    emstag(LWARN,"lib:prmfile:prm_read_darray","Parameter %s array size format bad", key);
    /*keyprm_errfn("prm_read_darray: Parameter %s array size format bad\n",
                 key);*/
    return -1;
  }

  if (*size > 0) {
    double val = 0.0;
    double lastval = 0.0;

/* allocate memory for the array if required */

    if (*p == NULL &&
        (*p = (double *)malloc((*size) * sizeof(double))) == NULL) {
      emstag(LWARN,"lib:prmfile:prm_read_darray","%s: Can't read array value",key);
      /*keyprm_errfn("prm_read_darray: %s: Can't read array value\n", key);*/
      return (-1);
    }

/**** Now populate the array, if we run out of values, then pad
 **** to end with the last value
 ****/
    for (i = 0; i < (*size); ++i) {
      if (fscanf(fp, "%lf", &val) != 1) {
        if (i == 0) {
          emstag(LWARN,"lib:prmfile:prm_read_darray","%s: Can't read any array values\n", key);
          /*keyprm_errfn
            ("prm_read_darray: %s: Can't read any array values\n", key);*/
          return -1;
        } else {
/* not enough values - fill remainder with last value */
          /* warn("prm_read_darray: %s: Filling to end of array with last
             value, %f\n",key, lastval); */
          for (; i < (*size); ++i)
            (*p)[i] = lastval;
          break;
        }
      } else {
/* A single value was located, so populate the array */
        (*p)[i] = val;
        lastval = val;
      }
    }
  }

  emstag(LTRACE,"lib:prmfile:prm_read_darray","%s  %d\n", key, *size);
  /* Echo to stdout if verbose enough */
  if (keyprm_verbose > 1 || is_log_enabled(LMETRIC))
  {
    for (i = 0; i < (*size); ++i)
      emslog(LMETRIC,"%f\n", (*p)[i]);
  }

  return 1;
}

/** Read a character string following a key.
  *
  * @param fp pointer to stdio FILE structure.
  * @param key keyname to locate in file.
  * @param p pointer to returned character string.
  * @return non-zero if successful.
  */
int prm_read_char(FILE * fp, char *key, char *p)
{
  char line[MAXLINELEN];
  char *s;
  char *r;

  /* Read parameter */
  if (!prm_skip_to_end_of_key(fp, key))
    return (0);

  if (fgets(line, MAXLINELEN, fp) != line)
    quit("prm_read_char: Can't read %s\n", key);

  /* Strip leading space */
  for (s = line; *s && is_blank(*s); s++) /* loop */
    ;

  /* Copy out result */
  for (r = p; *s && (*s != '\n'); *r++ = *s++)  /* loop */
    ;
  *r = 0;


  emstag(LTRACE,"lib:prmfile:prm_read_char","%s  %s", key, p);

  return (1);
}

/** Read a time value and convert to seconds.
  *
  * @param fp pointer to stdio FILE structure.
  * @param key keyname to locate in file.
  * @param val pointer to returned time value.
  * @return non-zero if successful.
  */
int prm_get_time_in_secs(FILE * fp, char *key, double *val)
{
  char buf[MAXLINELEN];

  if (!prm_skip_to_end_of_key(fp, key))
    return (0);

  if (fgets(buf, MAXLINELEN, fp) == NULL || !tm_scale_to_secs(buf, val))
    quit("timeparam: Can't read %s\n", key);

  return (1);
}

#if 0
/* scanf-like routine which matches a key at the start of
  * the line * before doing the read.
  *
  * @param fp pointer to stdio FILE structure.
  * @param key keyname to locate in file.
  * @param format scanf format character string.
  * @param ellipsis remaining arguments as per scanf.
  * @return non-zero if successful.
  *
  * @see scanf.
  */
#include <stdarg.h>

int keyfscanf(FILE * fp, char *key, char *format, ...)
{
  va_list args;
  int n;

  va_start(args, format);
  prm_skip_to_end_of_key(fp, key);
  n = vfscanf(fp, format, args);
  va_end(args);
  return (n);
}
#endif

/** Read the next non-blank non-comment line
  * from an ascii file. Uses fgets() to read
  * the line, so the terminating newline
  * character is read.
  * Comment lines are those lines which have
  * '#' as their first character.
  * 
  * @param line pointer to storage for line read
  * @param n max length of line to read
  * @param fp stream from which to read
  * 
  * @return non-zero if line was read successfilly.
  */
int prm_next_line(char *line, long int n, FILE * fp)
{
  char *s;

  do {
    s = fgets(line, n, fp);
  } while (s && (prm_is_comment_line(line) || prm_is_blank_line(line)));

  if (s == NULL)
    return (0);

  return (1);
}

int prm_is_comment_line(char *line)
{
  return (line[0] == '#');
}

int prm_is_blank_line(char *line)
{
  while (*line && isspace((int)*line))
    line++;
  if (*line)
    return (0);
  return (1);
}

/** Read stream to end of line.
  *
  * @param fp stream from which to read
  */
void prm_flush_line(FILE * fp)
{
  int c;

  while ((c = fgetc(fp)) != EOF && c != '\n')
    /* loop */ ;
}

/** Compiles a key string.
 * E.g., `prm_get_key(keybuf, "GRID0", "TRACER%1.1d.%s, 11, "long_name")'
 * results in `keybuf' equal to "GRID0.TRACER11.long_name"
 * while `prm_get_key(keybuf, NULL, "TRACER%1.1d.%s, 11, "long_name")'
 * prodices "TRACER11.long_name".
 * @param keybuf output
 * @param prefix input prefix like: "GRID0" or NULL
 * @param suffix input suffix
 * @param ... other input
 * @return key string
 */
char *prm_get_key(char *keybuf, const char *prefix, const char *suffix,
                  ...)
{
  char buffer[MAXLINELEN];
  va_list args;

  va_start(args, suffix);
  vsprintf(buffer, suffix, args);
  va_end(args);

  if (prefix != NULL)
    sprintf(keybuf, "%s.%s", prefix, buffer);
  else
    sprintf(keybuf, "%s", buffer);

  return keybuf;
}

