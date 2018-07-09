/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/io/nc_buffered.cpp
 *
 *  \brief Buffered access to some nc_get routines
 * 
 *  Wrappers around some nc_get_* routines to provide buffered access
 *  to data. That is, a circular buffer (of size N) is created for
 *  each timeslice of data and read in a threaded environment in
 *  parallel using OpenMP.
 *
 *  Notes:
 *  1. This assumes sequential reads. If we need to take care of
 *  arbitrary time then we'll need to increase the buffer size and
 *  have a loop that finds the correct buffer, discarding any previous
 *  time points - might need this for SP_INTERP
 *  Related issue: swr is directly read into the master, if it exists
 *  in the transport file. Needs to be specified as "LIGHT file"
 *
 *  2. Semaphores are used for flow of control while mutex locking
 *  for when reading from file. ncid seems to be a static within the
 *  nc_ routines
 *
 *  TODO:
 *   * Move data structure into df rather than v. Not an issue at the
 *   moment but will be if we are concurrently trying to read from
 *   multiple files. Right now its only implemented for the transport files
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id:$
 */

#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>
#include <string.h>
#include "netcdf.h"
#include "ems.h"

/* Globals for this file */
namespace nc_buffered {
  int NUM_NC_BUFFERS = 2;
  enum {
    iTIME,
    iDATA
  };
}
using namespace std;
using namespace nc_buffered;

/*
 * Holds a single buffer
 */
template <class T>
class ncBuf
{
public:
  /* No constructor, use default */
  void init(size_t len) {
    data = new T[len];
    time = -1; /* unset */
    sem_init(&done_read, 0, 0);
    sem_init(&done_copy, 0, 1); // force first read
  }
  /* Destructor */
  ~ncBuf(void) {
    delete[] data;
    sem_destroy(&done_read);
    sem_destroy(&done_copy);
  }

  /* Read data from file */
  void read_data(int ncid, int varid, size_t start[], size_t count[],
		 pthread_mutex_t *imtx) {
    int ret;

    /* Wait for copy to be done */
    sem_wait(&done_copy);

    /* 
     * Call netcdf rotuine to get data, need protection for
     * simultaneous variable reads
     */
    pthread_mutex_lock(imtx);
    ret = nc_get_vara_double(ncid, varid, start, count, data);
    pthread_mutex_unlock(imtx);
    if (ret != NC_NOERR)
      quit("nc_buffered netcdf error %s\n", nc_strerror(ret));

    /* Update time for this buffer */
    time = start[iTIME];

    /* We're done, signal we're done */
    sem_post(&done_read);
  }

  /* Return cached data */
  void get_data(size_t tindex, T *var, size_t len) {
    /* Warn when we have to wait, except for the very first time */
    if (time > 0) {
      int val;
      sem_getvalue(&done_read, &val);
      if (!val)
	warn("(lib:io:nc_buffered) Waiting for buffered read, timeIdx = %d\n", 
	     tindex);
    }
    sem_wait(&done_read);
    if (time != tindex)
      quit((char*)"(lib:io:nc_buffered) Incorrect data (%d,%d)\n", 
	   time, tindex);
    memcpy(var, data, len*sizeof(T));
    sem_post(&done_copy);
  }

private:
  /* Pointer to buffered data */
  T *data;

  /* Keep track of time index */
  int time;

  /* Use semaphores for synchronisation */
  sem_t done_read;
  sem_t done_copy;
};

/*
 * Container class for ncBuf's
 */
template <class T>
class ncBuffers
{
public:
  /* Nominal constructor */
  ncBuffers(int num, int nid, int vid, size_t time, size_t len, 
	    pthread_mutex_t *imtx)
    : tStartIdx(time)
  {
    size_t start[2];
    int unlimdim;

    /* Allocate buffers and initialise each */
    bufs = new ncBuf<T>[num];
    for (int i=0; i<num; i++) {
      bufs[i].init(len);
    }
    num_bufs = num;

    /* Set up variables needed for nc_get_vara */
    ncid  = nid;
    varid = vid;
    count[iTIME] = 1L;
    count[iDATA] = len;

    /* Get the number of records (i.e. unlimited dimension) in this file */
    if (nc_inq_unlimdim(ncid, &unlimdim) != NC_NOERR)
      quit((char*)"(lib:io:nc_buffered) Unlimited dimension not found\n");
    nc_inq_dimlen(ncid,	unlimdim, &nrecords);

    /* Initialise counters */
    curr_pos = next_pos = 0;

    /* Cache mutex */
    mtx = imtx;

    /* Create this thread */
    if (pthread_create(&thd, NULL, pth_loop, this))
      quit((char*)"(lib:io:nc_buffered) Unable to create thread\n");
  }
  /* Destructor - cleans up the buffers and locks*/
  ~ncBuffers(void) {
    delete[] bufs;
    pthread_cancel(thd);
  }
  /* Main public get function for the data */
  void get_data(size_t time, T *var);

private:
  /* 
   * This function needs to be static to conform to pthread_create's
   * start_routine
   */
  static void* pth_loop(void *ptr) {
    (static_cast<ncBuffers<T> *>(ptr))->doNcRead();
    return(NULL);
  }

  /* Main loop */
  void doNcRead(void);

  /* The starting time index */
  const size_t tStartIdx;

  /* The number of records in this nc file */
  size_t nrecords;

  /* Keep track of info for nc_get_vara */
  int    ncid;
  int    varid;
  size_t count[2];

  /* Pointer to the buffers and associated locks */
  ncBuf<T> *bufs;
  int   num_bufs;
  /*
   * Indicies to keep track of the current and next buffers
   * NB. It'd be nice if STL lists were circular
   */
  int curr_pos;
  int next_pos;

  /* Thread */
  pthread_t thd;
  /* Cache mutex */
  pthread_mutex_t *mtx;
};

/*
 * Keeps reading until the end of file and then returns, effectively
 * exiting the thread
 */
template <class T>
void ncBuffers<T>::doNcRead(void)
{
  int tIdx = tStartIdx;
  size_t start[2], stacksz;
  start[iDATA] = 0L;
  bool done = false;
  while(!done) {
    for (int i=0; i<num_bufs && !done; i++) {
      start[iTIME] = tIdx++;
      if (start[iTIME] < nrecords)
	bufs[i].read_data(ncid, varid, start, count, mtx);
      else
	done = true;
    }
  }
}

/*
 * Main workhourse method that keeps track of the next buffer
 */
template <class T>
void ncBuffers<T>::get_data(size_t time, T* var)
{
  /* Cache current position, i.e. data should be valid for this */
  int pos = curr_pos;

  /* Set the current position to next and release its lock */
  next_pos = (curr_pos + 1) % num_bufs;
  curr_pos = next_pos;
  
  /* Blocking get */
  bufs[pos].get_data(time, var, count[iDATA]);
}

/*
 * External C API's
 */
// This function is specific to a time indexed sparse array of type double
extern "C"
void nc_buffered_get_vara_double(int ncid, df_variable_t *v, int varid, 
				 size_t start[], size_t count[], double *var)
{
  /* 
   * Grab the object pointer of the buffers object for this variable
   * in this file
   */
  void              *vbuf = v->nc_bufs;
  ncBuffers<double> *bufs = NULL;

  /*
   * xxx Need to fix this
   * At the moment this is global so *any* read will be locked for all
   * threads. This is not the most optimum behaviour, need only lock
   * one per ncid.
   * Fix: Rework ncBuffers to hang off df and use a map to get at the
   *     variable buffers underneath.
   *     Also makes for destroy easier
   * Timing is a problem ... we may not want consecutive time either
   * due to SP_INTERP type issues (in which case we'll need to keep
   * nieghbouring times hanging around) or overlap in mutli files.
   * Fix: Add a dequeue that grows at one end and shrinks at the other
   *      if time is not available, then go get it and make a jump and
   *      re-init all the previous buffers ... do this with semaphores 
   *      Acutally this is not really an issue ... it all sorts itself
   *      out as a new file will have a new object!!! TRANS_VARS was
   *      the problem!!
   */
  static pthread_mutex_t *nc_mtx = NULL;

  /* Get time index and length of data array */
  size_t time = start[iTIME];
  size_t len  = count[iDATA];
  /* See if we've already created the buffers object */
  if (vbuf != NULL) {
    /* Yes */
    bufs = static_cast< ncBuffers<double> *>(vbuf);
  } else {
    /* Initialise mutex - one per ncid */
    if (nc_mtx == NULL) {
      nc_mtx = new pthread_mutex_t;
      pthread_mutex_init(nc_mtx, NULL);
    }
    /* First time, need to create it */
    bufs = new ncBuffers<double>(NUM_NC_BUFFERS, ncid, varid, time, len, nc_mtx);

    /* Set pointer on variable */
    v->nc_bufs = static_cast <void *>(bufs);
  }
  /* Blocking get */
  bufs->get_data(time, var);
}

/*
 * Cleanup functions
 */
extern "C"
void nc_buffered_clean_up_var(df_variable_t *v)
{
  void *vbuf = v->nc_bufs;
  if (vbuf != NULL) {
    delete (static_cast< ncBuffers<double> *>(vbuf));
    v->nc_bufs = NULL;
  }
}

/*
 * Loops over cleans up ALL variables
 */
extern "C"
void nc_buffered_clean_up_all_vars(datafile_t *df)
{
  int i;
  for (i = 0; i < df->nv; ++i)
    nc_buffered_clean_up_var(&df->variables[i]);
}
