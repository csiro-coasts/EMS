/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file misc/stringtable.c
 *
 *  \brief Stringtable abstract datatype
 *
 * Stringtable is an expandable array of strings and associated indices.
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: stringtable.c 5831 2018-06-26 23:48:06Z riz008 $
 */

#define STRINGTABLE_NSTART 50
#define STRINGTABLE_NINC 50

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "ems.h"

/** Compare procedure for qsort(). Non-case-sensitive.
 * @param se1 Pointer to stringentry
 * @param se2 Pointer to stringentry
 * @return 0 if se1 and se2 are identical, -1 or 1 otherwise.
 */
static int cmpse(const void* se1, const void* se2)
{
    return strcasecmp((*(stringentry**) se1)->s, (*(stringentry**) se2)->s);
}

/** Stringtable constructor.
 * @param name Stringtable name
 * @return Stringtable
 */
stringtable* stringtable_create(char* name)
{
    stringtable* st = malloc(sizeof(stringtable));

    if (name != NULL)
        st->name = strdup(name);
    else
        st->name = strdup("stringtable");
    st->unique = 1;
    st->n = 0;
    st->nallocated = STRINGTABLE_NSTART;
    st->sorted = 0;
    st->se = calloc(STRINGTABLE_NSTART, sizeof(void*));

    return st;
}

/** Stringtable destructor.
 * @param st Stringtable
 */
void stringtable_destroy(stringtable* st)
{
    int i;

    if (st == NULL)
        return;

    for (i = 0; i < st->n; ++i) {
        free(st->se[i]->s);
        free(st->se[i]);
    }
    free(st->name);
    free(st->se);
    free(st);
}

/** Resets contents of a stringtable.
 * @param st Stringtable
 */
void stringtable_reset(stringtable* st)
{
    st->n = 0;
}

/** Adds an entry to a stringtable.
 * @param st Stringtable
 * @param s String to be added
 * @param index External index of the added string to be stored
 */
void stringtable_add(stringtable* st, char* s, int index)
{
    stringentry* se = malloc(sizeof(stringentry));

    se->s = strdup(s);
    se->index = (index >= 0) ? index : st->n;
    se->naccess = 0;
    se->n = 0;

    if (st->unique)
        if (stringtable_findstringindex(st, s) >= 0) {
            /*
             * (there is no error quit function for stringtable for now)
             */
            emslog(LPANIC, "error: stringtable \"%s\": entry \"%s\" duplicated\n", st->name, s);
            exit(1);
        }

    if (st->n == st->nallocated) {
        st->se = realloc(st->se, (st->nallocated + STRINGTABLE_NINC) * sizeof(void*));
        st->nallocated += STRINGTABLE_NINC;
    }

    st->se[st->n] = se;
    st->n++;
    st->sorted = 0;
}

/** Adds an entry to a stringtable if it is not already there.
 * @param st Stringtable
 * @param s String to be added
 * @param index External index of the added string to be stored
 */
void stringtable_add_ifabscent(stringtable* st, char* s, int index)
{
    stringentry* se = malloc(sizeof(stringentry));

    se->s = strdup(s);
    se->index = (index >= 0) ? index : st->n;
    se->naccess = 0;
    se->n = 0;

    if (stringtable_findstringindex(st, s) >= 0) {
        free(se->s);
        free(se);
        return;
    }

    if (st->n == st->nallocated) {
        st->se = realloc(st->se, (st->nallocated + STRINGTABLE_NINC) * sizeof(void*));
        st->nallocated += STRINGTABLE_NINC;
    }

    st->se[st->n] = se;
    st->n++;
    st->sorted = 0;
}

/** Finds index associated with a string in a stringtable. Uses
 * binary search for a sorted table; cycles trhough all entries otherwise.
 * @param st Stringtable
 * @param s String
 * @return Index of the string if found; -1 otherwise
 */
int stringtable_findstringindex(stringtable* st, char* s)
{
    int index = -1;

    if (st->sorted) {
        stringentry tmp;

        tmp.s = s;
        tmp.index = -1;
	tmp.n = 0;

        {
            stringentry* se = &tmp;
            stringentry** ans = (stringentry**) bsearch(&se, st->se, st->n, sizeof(void*), cmpse);

            if (ans != NULL) {
                index = (*ans)->index;
                (*ans)->naccess++;
            }
        }
    } else {
        int n = st->n;
        int i;

        for (i = 0; i < n; ++i)
            if (strcasecmp(st->se[i]->s, s) == 0) {
                index = st->se[i]->index;
                st->se[i]->naccess++;
                break;
            }
    }

    return index;
}

/** Finds string associated with an index in a stringtable. Uses
 * linear search.
 * @param st Stringtable
 * @param index Index
 * @return String if found; NULL otherwise
 */
char* stringtable_findindexstring(stringtable* st, int index)
{
    int i;

    for (i = 0; i < st->n; ++i)
        if (st->se[i]->index == index)
            return st->se[i]->s;

    return NULL;
}


void stringtable_parse(stringtable* st, char* line)
{
	int i,n;
	char* pbuf[1000];
	n = parseline(line,pbuf,100);
	for(i=0; i < n ;i++)
		stringtable_add(st,pbuf[i],-1);
}


int stringtable_contains(stringtable* st, char* s)
{
	return (stringtable_findstringindex(st,s) > -1 );
}

/** Gets the value of n for the stringtable entry specified by its
 * index
 * @param st Stringtable
 * @param i Index
 * @return value of n
 */
int stringtable_entry_get_n(stringtable *st, int i)
{
  if (i < st->n)
    return(st->se[i]->n);
  else
    emstag(LERROR,"eco:stringtable:stringtable_entry_get_n",
	   "Given index (%d) outside range (0->%d)", i, st->n);
    
  return(0);
}


/** Sorts a stringtable using qsort().
 * @param st pointer to stringtable
 */
void stringtable_sort(stringtable* st)
{
    if (st->sorted)
        return;
    qsort(st->se, st->n, sizeof(void*), cmpse);
    st->sorted = 1;
}

/** Prints contents of a stringtable to standard error.
 * @param st Stringtable
 */
void stringtable_print(stringtable* st)
{
    /*UR consolidated */
    stringtable_print_fp(st,stderr);
    /*
    int i;
    fprintf(stderr, "%s:\n", st->name);
    for (i = 0; i < st->n; ++i)
        fprintf(stderr, "  %s(%d)\n", st->se[i]->s, st->se[i]->index);
    fflush(stderr);
    */
}


/** Prints contents of a stringtable to file stream.
 * @param st Stringtable
 * @param fp open file pointer
 */
void stringtable_print_fp(stringtable* st, FILE* fp)
{
    int i;
    if(fp == NULL)
      return;

    fprintf(fp, "%s:\n", st->name);
    for (i = 0; i < st->n; ++i)
        fprintf(fp, "  %s(%d)\n", st->se[i]->s, st->se[i]->index);
    fflush(fp);
}


/** Prints stringtable entries with specified separator to standard error.
 * @param st Stringtable
 * @param sep Separator
 */
void stringtable_printentries(stringtable* st, char* sep)
{
    int n = st->n;
    int i;

    assert(sep != NULL);

    if (n < 1)
        return;

    fprintf(stderr, "%s", st->se[0]->s);
    for (i = 1; i < n; ++i)
        fprintf(stderr, "%s%s", sep, st->se[i]->s);
    fprintf(stderr, "%s", sep);
    fflush(stderr);
}


/** Prints specified stringtable entry;
 * @param st Stringtable
 * @param i Index
 * @param fp open file pointer
 */
void stringtable_printentry_fp(stringtable* st, int i, FILE* fp)
{
    if(fp == NULL)
      return;
    if(st == NULL)
      emstag(LWARN,"eco:stringtable:stringtable_printentry_fp","StringTbale content NULL");
    else{
      fprintf(fp, "  %s: %s(%d)\n", st->name, st->se[i]->s, st->se[i]->index);
      fflush(fp);
    }
}


/** Prints specified stringtable entry;
 * @param st Stringtable
 * @param i Index
 */
void stringtable_printentry(stringtable* st, int i)
{
    stringtable_printentry_fp(st,i,stderr);
}


/** Return concatenate contents of a stringtable .
 *
 * @param st Stringtable
 * @param buf string buffer
 * @return string comma separated list of content
 */
char* stringtable_to_string(stringtable* st, char* buf)
{
    int i;
    if (st->n > 0)
        strcat(buf, st->se[0]->s);
    for (i = 1; i < st->n; ++i) {
        strcat(buf, ", ");
        strcat(buf, st->se[i]->s);
    }
    return buf;
}
