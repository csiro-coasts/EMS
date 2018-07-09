/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/misc/hqueue.c
 *
 *  \brief Hash-queue abstract datatype
 *
 *  Create, maintain and index into a queue.  Actually,
 *  this is a special kind of a hybrid hash/queue :
 *   - It has hash-like indexing where indexing on the first
 *     column i.e. key is allowed but IS NOT EFFICIENT
 *   - Hash fixed length, and any push's at the back means you
 *     lose the front entry.
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: hqueue.c 5831 2018-06-26 23:48:06Z riz008 $
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ems.h"

typedef struct {
  int    num_rows;
  int    num_cols;
  void ***table;
  int    num_entries;
} hqueue_t;

// fwd decl of local functions
static void *move_rows_up(hqueue_t *hq);

/*
 * Create function
 *  2 x 2 table of pointers
 *  NB. This is hard-coded for now as this is the only case that is
 *      currently being used, however, this really ought to be generalised
 */
hqueue *hq_create(int num_rows)
{
  int nrows = num_rows;
  int ncols = 2;
  int c,r;

  hqueue_t *hq = (hqueue_t *)malloc(sizeof(hqueue_t));
  if (hq == NULL)
    quit("hq_create: unable to allocate memory for the hqueue struct");
  
  hq->num_rows = nrows;
  hq->num_cols = ncols;
  hq->table    = p_alloc_2d(ncols, nrows);
  hq->num_entries = 0;

  // Loop over and stick in null pointers
  for (r=0; r<nrows; r++) {
    // allocate key
    int *key_ptr = (int *)malloc(sizeof(int));
    *key_ptr     = -1;
    hq->table[r][0] = key_ptr;

    for (c=1; c<ncols; c++) {
      // allocate value
      hq->table[r][c] = NULL;
    }
  }
    
  return((hqueue *)hq);
}

/*
 * Frees all memory
 */ 
void hq_destroy(void *hq_ptr) 
{
  hqueue_t *hq = (hqueue_t *)hq_ptr;
  p_free_2d(hq->table);
  free(hq);
}


/*
 * Discards the first entry in the table, moves all others up one and
 * then tags the new row at the end
 */
void *hq_push(hqueue *hq_ptr, int key, void *data_ptr)
{

  hqueue_t *hq = (hqueue_t *) hq_ptr;
  int r = hq->num_entries;

  void *old_ptr = NULL;

  if (r == hq->num_rows) {
    // The table is full so move everything up one
    old_ptr = move_rows_up(hq);
    r--;
  }

  // Add to the end
  *(int *)hq->table[r][0] = key;
  hq->table[r][1] = data_ptr;
  r++;
  hq->num_entries = r;

  return old_ptr;
}


/*
 * Loops through the table matching the key and returning the
 * corresponding value.
 * Note: This is NOT effecient as it has to loop through each element
 * one by one. However its not an issue for small sizes.
 */
void *hq_get_value_by_key(hqueue *hq_ptr, int key)
{
  int r;
  hqueue_t *hq = (hqueue_t*) hq_ptr;

  // simple linear search
  for (r=0; r<hq->num_rows; r++) {
    if (*(int *)hq->table[r][0] == key)
      return(hq->table[r][1]);
  }
  
  return(NULL);
}


/*
 * Prints out the table
 */
void print_hqueue(hqueue *hq_ptr)
{
  hqueue_t *hq = (hqueue_t*) hq_ptr;
  int nelems   = hq->num_entries;
  int r;
  
  printf("hq (%d):\n", nelems);
  for (r=0; r<nelems; r++) {
    int    key   = *(int *) hq->table[r][0];
    double value = *(double *) hq->table[r][1];

    printf("\t%d : %5d | %.10f\n", r, key, value);
  }
}

/*
 * Local functions
 */
static void *move_rows_up(hqueue_t *hq)
{
  void *old_ptr = hq->table[0][1];
  // Move everything up one and return the data pointer of the
  // deleted row
  int i;
  for (i=1; i<hq->num_rows; i++) {
    *(int *)hq->table[i-1][0] = *(int *)hq->table[i][0];
    hq->table[i-1][1] = hq->table[i][1];
  }

  return(old_ptr);
}



// EOF
