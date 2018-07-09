/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/hqueue.h
 *
 *  \brief Header for hash-queue ADT
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: hqueue.h 5834 2018-06-27 00:55:37Z riz008 $
 */


#ifndef _HQUEUE_H
#define _HQUEUE_H

typedef struct hqueue_t hqueue;

/*
 * Public API's
 */
hqueue *hq_create(int num_rows);

void *hq_push(hqueue *hq_ptr, int key, void *data_ptr);

void hq_destroy(void *hq_ptr);

void *hq_get_value_by_key(hqueue *hq_ptr, int key);

void print_hqueue(hqueue *hq_ptr);

/*
 * These generic functions have not been implemented yet
 */
/*
double hq_get_key_by_index(hqueue *hq, int index);
double hq_get_value_by_index(hqueue *hq, int index);
double hq_get_first_key(hqueue *hq);
double hq_get_first_value(hqueue *hq);
double hq_get_last_key(hqueue *hq);
double hq_get_last_value(hqueue *hq);
*/


#endif

/* EOF */
