/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/misc/hash.c
 *
 *  \brief Hashh infrastructure
 *
 *  Create, populate and locate items within a hash table
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: hash.c 5833 2018-06-27 00:21:35Z riz008 $
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include "hash.h"

// FR : This holds true for most chips but need to make sure. This
// should also come as part of configure not hard-coded as such
#define INT_PER_DOUBLE 2
#define INT_PER_POINTER 2 // Need to solve this for anything other than 64 bit!
#define BYTE_PER_INT 4

// local prototypes
static void ht_free_pair(void *key, void *data);
static void ht_free_destroy(void *key, void *data);

size_t ht_string_hash(void *s)
{
  size_t i, n;
  size_t hashval = 0;
  size_t ival = 0;

  n = sizeof(size_t);
  for (i = 0; i < strlen((char *)s); i += n) {
    memcpy(&ival, &((char *)s)[i], n);
    hashval += ival;
  }

  return hashval;
}


hash_table_t *ht_create(size_t nelems, ht_hash makehash, ht_compare compare)
{
  hash_table_t *ht = (hash_table_t *)malloc(sizeof(hash_table_t));
  emstag(LTRACE,"lib:hash.c:ht_create","Creating Hashtable with size %d ",nelems);
  memset(ht, 0, sizeof(hash_table_t));
  ht->elems = (hash_element_t **)malloc(sizeof(hash_element_t *) * nelems);
  memset(ht->elems, 0, sizeof(hash_element_t *) * nelems);
  ht->nelems = nelems;
  ht->noccupied = 0;
  ht->makehash = makehash;
  ht->compare = compare;
  ht->free = NULL;
  return ht;
}

hash_table_t *ht_create_complex(size_t nelems,
				ht_hash makehash,
				ht_compare compare,
				ht_copy    copy,
				ht_free    free)
{
  hash_table_t *ht = ht_create(nelems,makehash,compare);
  emstag(LTRACE,"lib:hash.c:ht_create","Creating Hashtable with size %d ",nelems);
  ht->copy = copy;
  ht->free = free;
  return ht;
}


int ht_string_compare(void *s1, void *s2)
{
  return !strcmp((char *)s1, (char *)s2);
}

/*
 * Generic add function.
 */
void ht_add(hash_table_t *ht, void *key, void *data)
{
  hash_element_t *newhash;
  hash_element_t *curhash;
  size_t hashval;

  newhash = (hash_element_t *)(malloc(sizeof(hash_element_t)));
  memset(newhash, 0, sizeof(hash_element_t));
  // Use the copy function if available.  Needed for multidimensional keys
  if (ht->copy != NULL) {
    newhash->key = ht->copy(key);
  } else {
    newhash->key = key;
  }
  newhash->data = data;
  newhash->free = NULL;

  hashval = ht->makehash(key) % ht->nelems;

  if (ht->elems[hashval] == NULL) {
    ht->elems[hashval] = newhash;
    ht->elems[hashval]->child = NULL;
    ++ht->noccupied;
  } else {
    curhash = ht->elems[hashval];
    while (curhash->child != NULL)
      curhash = curhash->child;

    curhash->child = newhash;
    newhash->child = NULL;
  }
  ht->nentries++;
}

void ht_add_complex(hash_table_t *ht, void *key, void *data, 
		    ht_free free)
{
  hash_element_t *newhash;
  hash_element_t *curhash;
  size_t hashval;

  newhash = (hash_element_t *)(malloc(sizeof(hash_element_t)));
  memset(newhash, 0, sizeof(hash_element_t));
  newhash->key  = key;
  newhash->data = data;
  newhash->free = free;

  hashval = ht->makehash(key) % ht->nelems;

  if (ht->elems[hashval] == NULL) {
    ht->elems[hashval] = newhash;
    ht->elems[hashval]->child = NULL;
    ++ht->noccupied;
  } else {
    curhash = ht->elems[hashval];
    while (curhash->child != NULL)
      curhash = curhash->child;

    curhash->child = newhash;
    newhash->child = NULL;
  }
  ht->nentries++;
}

void *ht_find(hash_table_t *ht, void *key)
{
  size_t hashval;
  hash_element_t *curhash;

  // Save us jumping into the hash table at all.
  if (ht->lastelem != NULL) {
    if (ht->compare(ht->lastelem, key))
      return ht->lastelem->data;
  }

  hashval = ht->makehash(key) % ht->nelems;

  if ((curhash = ht->elems[hashval]) == NULL)
    return NULL;

  do {
    if (ht->compare(curhash->key, key)) {
      ht->lastelem = curhash;
      return curhash->data;
    }
    curhash = curhash->child;
  } while (curhash != NULL);

  return NULL;
}


void ht_print_stats(hash_table_t *ht)
{
  emslog(LINFO,"Hash table size         = %d\n", ht->nelems);
  emslog(LINFO,"  Number occupied cells = %d\n", ht->noccupied);
  emslog(LINFO,"  Total number elements = %d\n", ht->nentries);
}


/*UR-ADDED missing destroy function, caused substantial memory leak
 * at the same time make this more self-sufficient and provide the option
 * of freeing more complex objects
 */
static void ht_element_free(hash_table_t *ht, hash_element_t* el)
{
  hash_element_t* prev;
  if(el == NULL)
    return;

  do{
    /* check if there is a individual function for this object */
    if(el->free != NULL)
      el->free(el->key,el->data);
    /* is there a general one */
    else if(ht->free != NULL)
      ht->free(el->key,el->data);

    /*now get rid of it */
    prev = el;
    el = el->child;
    free(prev);
  }while( el != NULL);
}

/* Deletes an entry from the table.  Returns a pointer to the data that
 * was associated with the key so that the calling code can dispose it
 * properly.
 *
 * This should be the exact opposite of ht_add
 *
 * @param table The hash table
 * @param key The key
 * @return The associated data or NULL
 */
void *ht_delete(hash_table_t* ht, void* key)
{
  hash_element_t *curhash;
  void *data = NULL;
  size_t hashval = ht->makehash(key) % ht->nelems;
  
  // Found nothing, so return nothing
  if ((curhash = ht->elems[hashval]) == NULL)
    return NULL;

  /*
   * Found a match
   */
  if (ht->compare(curhash->key, key)) {
    // Grab the data to return
    data = curhash->data;
    /*
     * See if it has any children
     */
    if (curhash->child != NULL) {
      // Re-jig the list
      ht->elems[hashval] = curhash->child;
      ht->lastelem = curhash->child;
    } else {
      /*
       * Remove from the hash
       */ 
      ht->elems[hashval] = NULL;
      ht->lastelem = NULL;
      ht->noccupied--;
    }
  } else {
    /*
     * Loop over the list
     */
    while ( (curhash = curhash->child) != NULL ) {
      if (ht->compare(curhash->key, key)) {
	// found it, so re-jig the list
	curhash->child = curhash->child->child;
	data = curhash->data;
	break;
      }
    }
  }
  
  /*
   * Delete the element and return the data
   */
  if (curhash != NULL) {
    ht_element_free(ht, curhash);
    ht->nentries--;
    return data;
  }
  return(NULL); // we found nothing
}

/*
 * Destroy the whole hash table, frees element and key memory
 */
void ht_destroy(hash_table_t *ht)
{
  size_t i;
  if(is_log_enabled(LTRACE))
    ht_print_stats(ht);

  for(i = 0;i < ht->nelems ;i++)
    ht_element_free(ht,ht->elems[i]);

  free(ht);
}


void ht_void_destroy(void* ht)
{
  ht_destroy((hash_table_t*)ht);
}

/* For each entry, calls a specified function with corresponding data as a
 * parameter.
 *
 * @param table The hash table
 * @param func The action function
 */
void ht_process(hash_table_t* ht, void (*func) (void*))
{
  int i;
  
  for (i = 0; i < ht->nelems; ++i) {
    hash_element_t *ht_elem = ht->elems[i];
    if (ht_elem != NULL)
      do {
	func(ht_elem->data);
      } while ((ht_elem = ht_elem->child) != NULL);
  }
}

int ht_getnentries(hash_table_t* ht)
{
    return ht->nentries;
}

int ht_getnelems(hash_table_t* ht)
{
    return ht->nelems;
}

int ht_getnoccupied(hash_table_t* ht)
{
    return ht->noccupied;
}

/*********************/
/* PRIVATE FUNCTIONS */
/*********************/
/*
 * Note: for the _hash functions, the return type must be size_t but
 * we cannot dereference the key as such because on 64 bit platforms
 * size_t is typedef'd to an unsigned long which is 8 byte wide.
 * The implicit cast of the return value is okay
 */

/* functions for double keys */
static size_t d1_hash(void* key)
{
  unsigned int* v = key;

#if INT_PER_DOUBLE == 2
  return v[0] + v[1];
#else
#error not implemented
#endif
}

static int d1_compare(void* key1, void* key2)
{
  return *(double*) key1 == *(double*) key2;
}

static void* d1_copy(void* key)
{
  double* newkey = malloc(sizeof(double));
  
  *newkey = *(double*) key;
  
  return newkey;
}

/* 
 * functions for for double[2] keys 
 */
static size_t d2_hash(void* key)
{
  unsigned int* v = key;
  
#if INT_PER_DOUBLE == 2
  /*
   * PS: here multiplications suppose to make (a,b) and (b,a) generate
   * different hash values 
   */
  return v[0] + v[1] + v[2] * 3 + v[3] * 7;
#else
#error not implemented
#endif
}

static int d2_compare(void* key1, void* key2)
{
    return (((double*) key1)[0] == ((double*) key2)[0]) && (((double*) key1)[1] == ((double*) key2)[1]);
}

static void* d2_copy(void* key)
{
    double* newkey = malloc(sizeof(double) * 2);

    newkey[0] = ((double*) key)[0];
    newkey[1] = ((double*) key)[1];

    return newkey;
}

/* 
 * functions for for double[3] keys 
 */
static size_t d3_hash(void* key)
{
  unsigned int* v = key;
  
#if INT_PER_DOUBLE == 2
  return v[0] + v[1] + v[2] * 3 + v[3] * 7 + v[4] * 11 + v[5] * 13;
#else
#error not implemented
#endif
}

static int d3_compare(void* key1, void* key2)
{
    return (
	    (((double*) key1)[0] == ((double*) key2)[0]) && 
	    (((double*) key1)[1] == ((double*) key2)[1]) && 
	    (((double*) key1)[2] == ((double*) key2)[2])
	    );
}

static void* d3_copy(void* key)
{
    double* newkey = malloc(sizeof(double) * 3);

    newkey[0] = ((double*) key)[0];
    newkey[1] = ((double*) key)[1];
    newkey[2] = ((double*) key)[2];

    return newkey;
}

/* 
 * functions for for int[1] keys 
 */
static size_t i1_hash(void* key)
{
  unsigned int *v = key;

  return v[0];
}


static int i1_compare(void* key1, void* key2)
{
  return (((int*) key1)[0] == ((int*) key2)[0]);
}

static void* i1_copy(void* key)
{
  int* newkey = malloc(sizeof(int));
  
  newkey[0] = ((int*) key)[0];
  
  return newkey;
}

/* 
 * functions for for int[2] keys 
 */
static size_t i2_hash(void* key)
{
#if BYTE_PER_INT >= 4
  unsigned int* v = key;
  
  return v[0] + (v[1] << 16);
#else
#error not implemented
#endif
}

static int i2_compare(void* key1, void* key2)
{
  return (((int*) key1)[0] == ((int*) key2)[0]) && 
                        (((int*) key1)[1] == ((int*) key2)[1]);
}

static void* i2_copy(void* key)
{
  int* newkey = malloc(sizeof(int) * 2);
  
  newkey[0] = ((int*) key)[0];
  newkey[1] = ((int*) key)[1];

  return newkey;
}

/*
 * Generic function to free keys only
 * The hashtables that use this free function need to free the data
 * using the process method
 */
static void ht_key_free_only(void *key, void *data)
{
  free(key);
}


/*
 * Generic function to free and value pairs
 */
static void ht_free_pair(void *key, void *data)
{
  free(key);
  free(data);
}

/*
 * Function that free's the key and destroys the hash table value
 */
static void ht_free_destroy(void *key, void *data)
{
  free(key);
  ht_void_destroy(data);
}

/****************************/
/* PUBLIC WRAPPER FUNCTIONS */
/****************************/
/** Create hash table for a single double element
 *  @param size size of hash elements
 */
hash_table_t* ht_create_d1(int size)
{
  assert(sizeof(double) == INT_PER_DOUBLE * sizeof(int));
  return ht_create_complex(size, d1_hash, d1_compare, d1_copy, ht_key_free_only);
}

/** Create hash table for a pair of double elements
 *  @param size size of hash elements
 */
hash_table_t* ht_create_d2(int size)
{
  assert(sizeof(double) == INT_PER_DOUBLE * sizeof(int));
  return ht_create_complex(size, d2_hash, d2_compare, d2_copy, ht_key_free_only);
}

/** Create hash table for triple double elements
 *  @param size size of hash elements
 */
hash_table_t* ht_create_d3(int size)
{
  assert(sizeof(double) == INT_PER_DOUBLE * sizeof(int));
  return ht_create_complex(size, d3_hash, d3_compare, d3_copy, ht_key_free_only);
}

/** Create hash table for a single integer element
 *  @param size size of hash elements
 */
hash_table_t* ht_create_i1(int size)
{
  return ht_create_complex(size, i1_hash, i1_compare, i1_copy, ht_key_free_only);
}

/** Create hash table for a pair of integer elements
 *  @param size size of hash elements
 */
hash_table_t* ht_create_i2(int size)
{
  assert(sizeof(int) == BYTE_PER_INT);
  return ht_create_complex(size, i2_hash, i2_compare, i2_copy, ht_key_free_only);
}


// EOF
