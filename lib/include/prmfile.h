/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/prmfile.h
 *
 *  \brief Parameter file reader prototypes
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: prmfile.h 6866 2021-07-23 05:23:13Z her127 $
 */


#ifndef	_PRMFILE_H
#define _PRMFILE_H

#include <stdio.h>
#include "errfn.h"

/* Prototypes */
void prm_flush_line(FILE * fp);
int prm_next_line(char *line, long n, FILE * fp);
int prm_is_comment_line(char *line);
int prm_is_blank_line(char *line);
char *prm_get_key(char *keybuf, const char *prefix, const char *suffix,
                  ...);
           
void prm_set_errfn(errfn fn);
void errfn_set_warn(errfn fn);
void prm_set_case(int c);
int prm_skip_to_start_of_key(FILE * fp, char *key);
int prm_skip_to_end_of_key(FILE * fp, char *key);
int prm_read_int(FILE * fp, char *key, int *p);
int prm_read_double(FILE * fp, char *key, double *p);
int prm_read_darray(FILE * fp, char *key, double **p, int *size);
int prm_read_char(FILE * fp, char *key, char *p);
int keyfscanf(FILE * fp, char *key, char *format, ...);
int prm_get_time_in_secs(FILE * fp, char *key, double *v);

/* External variables */
extern int keyprm_verbose;

#endif                          /* _PRMFILE_H */
