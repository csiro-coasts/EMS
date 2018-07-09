/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/misc/istack.c
 *
 *  \brief Integer stack abstract datatype
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: istack.c 5831 2018-06-26 23:48:06Z riz008 $
 */

#define STACK_NSTART 50
#define STACK_NINC 50

#include <stdlib.h>
#include <string.h>
#include "istack.h"

istack* istack_create(void)
{
    istack* s = malloc(sizeof(istack));

    s->n = 0;
    s->nallocated = STACK_NSTART;
    s->v = malloc(STACK_NSTART * sizeof(int));
    return s;
}

void istack_destroy(istack* s)
{
    if (s != NULL) {
        free(s->v);
        free(s);
    }
}

void istack_reset(istack* s)
{
    s->n = 0;
}

int istack_contains(istack* s, int v)
{
    int i;

    for (i = 0; i < s->n; ++i)
        if (s->v[i] == v)
            return 1;
    return 0;
}

void istack_push(istack* s, int v)
{
    if (s->n == s->nallocated) {
        s->nallocated *= 2;
        s->v = realloc(s->v, s->nallocated * sizeof(int));
    }

    s->v[s->n] = v;
    s->n++;
}

int istack_pop(istack* s)
{
    s->n--;
    return s->v[s->n];
}

int istack_getnentries(istack* s)
{
    return s->n;
}

int* istack_getentries(istack* s)
{
    return s->v;
}
