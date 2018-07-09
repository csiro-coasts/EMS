/**
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/istach.h
 *
 *  \brief Header file for istack
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: istack.h 5834 2018-06-27 00:55:37Z riz008 $
 */

#if !defined(_ISTACK_H)
#define _ISTACK_H

#if !defined(_ISTACK_STRUCT)
#define _ISTACK_STRUCT
struct istack;
typedef struct istack istack;
#endif

struct istack {
    int n;
    int nallocated;
    int* v;
};

istack* istack_create(void);
void istack_destroy(istack* s);
void istack_push(istack* s, int v);
int istack_pop(istack* s);
int istack_contains(istack* s, int v);
void istack_reset(istack* s);

#endif
