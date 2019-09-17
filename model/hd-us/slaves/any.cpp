/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/slaves/any.cpp
 *  
 *  Description:
 *  Optimised ANY function - use STL map as the look up table
 *
 *  This data structure is effectively a sorted array that does a
 *  binary search
 *
 *  To use:
 *    o) Initialise the array
 *    o) Use the find function within the loop, and
 *    o) Destroy the map
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id:$
 *
 */
#include <stdio.h>
#include <map>

/*
 * Uses the values of the array to create the map
 *  - inserts into a sorted binary tree
 */
extern "C"
void *ANY_map_init(int array[], int ns)
{
  std::map<int, int> *mparr = new std::map<int, int>;
  /* First create the map -no prealloc */
  for (int i=0; i<ns; i++) {
    mparr->insert( std::pair<int,int>(array[i], i) );
  }
  return(mparr);
}

/*
 * Need to use te find() method as the square brackets method [] will
 * create the key if it doesn't exist, whereas find does not
 */
extern "C"
int ANY_map_find(void *arr, int var)
{
  std::map <int, int> *mparr = static_cast< std::map <int, int> *> (arr);
  std::map <int, int>::iterator it;
  it = mparr->find(var);
  if (it != mparr->end()) {
    // value
    return((*it).second);
  }
  return(0);
}

/*
 * Clears memory
 */
extern "C"
void ANY_map_destroy(void *arr)
{
  std::map <int, int> *mparr = static_cast< std::map <int, int> *> (arr);
  mparr->clear();
  delete mparr;
}
