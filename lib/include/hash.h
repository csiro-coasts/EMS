/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/hash.h
 *
 *  \brief Creating and managing hash tables
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: hash.h 5834 2018-06-27 00:55:37Z riz008 $
 */

#ifndef _HASH_H
#define _HASH_H

#include "errfn.h"
#include "emslogger.h"

typedef struct hash_element hash_element_t;

/*
 * This is the hashing function itself
 */
typedef size_t (*ht_hash) (void *);

/*
 * Compares two keys
 */
typedef int (*ht_compare) (void *, void *);

/*
 * Generic copy function for keys
 */
typedef void* (*ht_copy) (void *);				

/*
 * Generic free function
 */
typedef void (*ht_free) (void *,void *);

struct hash_element {
  hash_element_t *child;
  void *key;
  void *data;
  ht_free free;
};

typedef struct {
  size_t nelems;
  size_t noccupied;
  size_t nentries;
  hash_element_t **elems;
  ht_hash    makehash;
  ht_compare compare;
  ht_copy    copy;
  ht_free    free;
  hash_element_t *lastelem;
} hash_table_t;


extern hash_table_t *ht_create(size_t nelems,
                               ht_hash makehash,
                               ht_compare compare);
extern hash_table_t *ht_create_complex(size_t nelems,
				       ht_hash makehash,
				       ht_compare compare,
				       ht_copy    copy,
				       ht_free    free);

extern size_t ht_string_hash(void *s);
extern int ht_string_compare(void *s1, void *s2);
extern void ht_add(hash_table_t *, void *, void *);
extern void ht_add_complex(hash_table_t *, void *,
                           void *,
			   ht_free free);
extern void *ht_find(hash_table_t *, void *);
extern void ht_print_stats(hash_table_t *ht);
/*UR-ADDED */
extern void ht_destroy(hash_table_t *ht);
extern void ht_void_destroy(void* ht);

/* FR added */
hash_table_t* ht_create_d1(int size);
hash_table_t* ht_create_d2(int size);
hash_table_t* ht_create_d3(int size);
hash_table_t* ht_create_i1(int size);
hash_table_t* ht_create_i2(int size);

/* Deletes an entry from the table.  Returns a pointer to the data that
 * was associated with the key so that the calling code can dispose it
 * properly.
 *
 * @param table The hash table
 * @param key The key
 * @return The associated data or NULL
 */
void* ht_delete(hash_table_t* ht, void* key);

/* For each entry, calls a specified function with corresponding data as a
 * parameter.
 *
 * @param table The hash table
 * @param func The action function
 */
void ht_process(hash_table_t* table, void (*func) (void*));
int ht_getnentries(hash_table_t* ht);
int ht_getnelems(hash_table_t* ht);
int ht_getnoccupied(hash_table_t* ht);

#endif
