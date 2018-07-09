/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/stringtable.h
 *
 *  \brief Stringtable is an expandable array of strings and associated indices
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: stringtable.h 5834 2018-06-27 00:55:37Z riz008 $
 */

#if !defined(_STRINGTABLE_H)

#include <stdio.h>

typedef struct {
  char* s;
  int index;   /* 
		* An index. Is set during the string entry.
		* If "-1" entered, is set to index of this
		* entry in the stringtable.
		*/
  int naccess; /* Number of times looked for. */
  int n;       /* 
		* Size associated with this entry - used in ecology 
		* NOTE: Instead of this specific integer, another way
		*       would to create a void pointer for a user struct
		*/
} stringentry;

#if !defined(STRINGTABLE_STRUCT)
struct stringtable;
typedef struct stringtable stringtable;

#define STRINGTABLE_STRUCT
#endif

struct stringtable {
    char* name;
    int unique;                 /* flag: whether all entries must be unique;
                                 * 1 by default */
    int n;
    int nallocated;
    int sorted;                 /* flag */
    stringentry** se;
};

/** Stringtable constructor.
 * @param name Stringtable name
 * @return Stringtable
 */
stringtable* stringtable_create(char* name);

/** Stringtable destructor.
 * @param st Stringtable
 */
void stringtable_destroy(stringtable* st);

/** Adds an entry to a stringtable.
 * @param st Stringtable
 * @param s String to be added
 * @param index External index of the added string to be stored
 */
void stringtable_add(stringtable* st, char* s, int index);

/** Adds an entry to a stringtable if it is not already there.
 * @param st Stringtable
 * @param s String to be added
 * @param index External index of the added string to be stored
 */
void stringtable_add_ifabscent(stringtable* st, char* s, int index);

/** Finds index associated with a string in a stringtable. Uses
 * binary search for a sorted table; cycles trhough all entries otherwise.
 * @param st Stringtable
 * @param s String
 * @return Index of the string if found; -1 otherwise
 */
int stringtable_findstringindex(stringtable* st, char* s);

/** Finds string associated with an index in a stringtable. Uses
 * linear search.
 * @param st Stringtable
 * @param i Index
 * @return String if found; NULL otherwise
 */
char* stringtable_findindexstring(stringtable* st, int i);

/** Resets contents of a stringtable.
 * @param st Stringtable
 */
void stringtable_reset(stringtable* st);

/** Sorts a stringtable using qsort().
 * @param st Stringtable
 */
void stringtable_sort(stringtable* st);


/** test if the s tring table contains this entry
 * @param st the stringtable to test for s
 * @param s - the string totest for
 * @return 1 if st contains s, 0 otherwise
 */
int stringtable_contains(stringtable* st, char* s);


/** Split line l by spaces and add the content to string table st
 *
 * @param st - the string table to add to
 * @param l - the line to parse
 */
void stringtable_parse(stringtable* st, char* l);


/** Prints contents of a stringtable to standard error.
 * @param st Stringtable
 */
void stringtable_print(stringtable* st);


/** Return concatenate contents of a stringtable .
 *
 * @param st Stringtable
 * @return string comma separated list of content
 */
char* stringtable_to_string(stringtable* st, char* buf);


/*UR ADDED */
/** Prints contents of a stringtable to a file stream.
 * @param st Stringtable
 * @param fp the file stream
 */
void stringtable_print_fp(stringtable* st, FILE* fp);


/** Prints stringtable entries with specified separator to standard error.
 * @param st Stringtable
 * @param sep Separator
 */
void stringtable_printentries(stringtable* st, char* sep);

/** Prints specified stringtable entry;
 * @param st Stringtable
 * @param i Index
 */
void stringtable_printentry(stringtable* st, int i);


/** Prints specified stringtable entry;
 * @param st Stringtable
 * @param i Index
 * @param fp outputstream
 */
void stringtable_printentry_fp(stringtable* st, int i, FILE* fp);


/** Gets the value of n for the stringtable entry specified by its
 * index
 * @param st Stringtable
 * @param i Index
 * @return value of n
 */
int stringtable_entry_get_n(stringtable *st, int i);

#define _STRINGTABLE_H
#endif
