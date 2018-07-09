/**
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/eqn_parser.h
 *
 *  \brief Prototypes for the equation parser
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: eqn_parser.h 5834 2018-06-27 00:55:37Z riz008 $
 */

#ifndef _EQN_PARSER_H
#define _EQN_PARSER_H

typedef double *(*eqn_custom_fcn)(const char *str, void *data);

void *EqnCreateParser(const char *str, eqn_custom_fcn fcn, void *data, char *err);
char *EqnDisplayStr(void *e);
double EqnGetValue(void *e);
void EqnFree(void *e);

#endif

/* end of eqn_parser.h */
