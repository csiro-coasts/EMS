/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/slaves/TransferBuffer.cpp
 *  
 *  Description: Classes to hold transfer buffers and MPI comms routines
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: TransferBuffer.cpp 5841 2018-06-28 06:51:55Z riz008 $
 *
 * The class hierarchy is as follows:
 *
 *                           /----------------\
 *			     | TransferBuffer |    Base class that holds lhs=rhs pointers
 *			     \----------------/
 *			              |
 *	    +-------------------------|
 *	    |			      V
 *	    |		    /-------------------\
 *	    |		    | TransferBufferBuf |  Adds contiguous
 *	    | 		    \-------------------/  buffers to facilitate transfers
 *	    |		              |
 *	    |			      |-------------+---------+
 *	    V			      V             V         V
 *	/-------\               /----------\     /-----\  /---------\
 *	| ShMem |  (others..)   | ShMemBuf |     | MPI |  | MPI-RMA |  etc...
 *	\-------/               \----------/     \-----/  \---------/
 *    
 *    Shared memory             Shared memory     MPI implementation using
 *    using direct              using buffers     buffers, two/one sided
 *       access                                       communications
 *
 *  Its also possible to implement the classes as separate hierarhies
 *  and then use multiple inheritance to mix and match whatever we
 *  need to implement for the final product but for now this does the job
 *
 * Farhan Jul. 2015
 */



#include <string.h>
#include <numeric>
#include <vector>
#include <map>
#include <sys/time.h>

#include "ems_conf.h"
#ifdef HAVE_MPI
#include <mpi.h>
#endif

extern "C" {
#include "hd.h"
}

#define WRITE_BUF_SIZES (1)

using namespace std;

/*
 * Base class to hold pointers for transfers
 */
class TransferBuffer
{
public:
  /* Constructor */
  TransferBuffer(int nwin):m_nwin(nwin) {
    /* Allocate memory - windows are 1-based */
    m_lhs = new vector <void *>[nwin+1];
    m_rhs = new vector <void *>[nwin+1];
    m_sz  = new vector <size_t>[nwin+1];
  }
  /* Destructor */
  virtual ~TransferBuffer(void) {
    delete[] m_lhs;
    delete[] m_rhs;
    delete[] m_sz;
  }
  
  /* Initialisation method */
  virtual void init(void) {
    /* nothing to do here */
  }
  
  /* Method to add pointers, keeping track of the datatype size */
  template <class T>
  inline void add(T* l, T* r, int nw) {
    m_lhs[nw].push_back(l);
    m_rhs[nw].push_back(r);
    m_sz[nw].push_back(sizeof(T));
  }
  /* Overloaded for master_fill */
  template <class T>
  inline void add(T* l, T* r) {
    m_lhs[1].push_back(l);
    m_rhs[1].push_back(r);
    m_sz[1].push_back(sizeof(T));
  }

  /* 
   * Various comms functions
   */
  virtual void send2w(int wn) = 0;
  virtual void send2m(void)   = 0;
  virtual void wrecv(void)    = 0;
  virtual void mrecv(int wn)  = 0;
  virtual void transfer2m(int n) = 0;

  void printSizes(FILE *f, const char *desc) { 
     fprintf(f, " %s\n", desc);
     for (int i=1; i<=m_nwin; i++) {
       size_t sz = getSize(i);
       if (sz)
	 fprintf(f, "  <- %2d = %.2e\n", i, (double)sz);
       else
	 fprintf(f, "  <- %2d = -\n", i);
     }
  }

  size_t getSize(int i) {
    return(accumulate(m_sz[i].begin(), m_sz[i].end(), (size_t)0));
  }
  
protected:
  /*
   * The total number of windows
   */
  const int m_nwin;
  /*
   * Vectors to hold pointers for the lhs, rhs and size of the element
   * they are pointing (i.e. datatype size) respectively
   */
  vector <void *> *m_lhs;
  vector <void *> *m_rhs;
  vector <size_t> *m_sz;
};

/*
 * Subclass that uses a contiguous buffer to faciliate transfers
 */
class TransferBufferBuf : public TransferBuffer
{
public:
  /* Constructor */
  TransferBufferBuf(int nwin):TransferBuffer(nwin) {
    // one-based
    m_bufPtrs = new char*[nwin+1];
    m_bufSzs  = new size_t[nwin+1];
  }
  /* Destructor */
  ~TransferBufferBuf(void) {
    delete[] m_bufPtrs;
    delete[] m_bufSzs;
  } 
  /* 
   * Initialise the message passing buffers 
   * Note : This can only be called *after* all the add operations are
   *        finished with
   */
  inline void init(void) {
    for (int i=1; i<=m_nwin; i++) {
      /* Get the total size for each buffer */
      size_t sz = getSize(i);
    
      m_bufPtrs[i] = new char[sz];
      m_bufSzs[i]  = sz;
    }
  }

protected:
  /*
   * Methods to fill-in/take-out from the contiguous buffers
   */
  void bufIn(int n)
  {
    /* RHS values go into the buffer */
    char *ptr = m_bufPtrs[n];
    for (int i=0; i<m_rhs[n].size(); i++) {
      memcpy(ptr, m_rhs[n][i], m_sz[n][i]);
      ptr += m_sz[n][i];
    }
  }
  void bufOut(int n)
  {
    /* LHS values are filled in from the buffer */
    char *ptr = m_bufPtrs[n];
    for (int i=0; i<m_lhs[n].size(); i++) {
      memcpy(m_lhs[n][i], ptr, m_sz[n][i]);
      ptr += m_sz[n][i];
    }
  }

  /*
   * Contiguous buffers for message passing 
   */
  char  **m_bufPtrs;
  size_t *m_bufSzs;

};


/*
 * Class that implements the shared memory protocol
 *   i.e. 
 */
class TransferBufferShMem : public TransferBuffer
{
public:
  TransferBufferShMem(int nwin):TransferBuffer(nwin)
  { 
    /* Nothing specific here */
  }
  
  /* No sends */
  virtual void send2w(int wn) { }
  virtual void send2m(void)   { }

  /*
   * This is a shared memory object, so no send/recv's -
   * straight off memcpy for all buffers
   */
  inline void directXfer(void) {
    for (int n=1; n <= m_nwin; n++)
      for (int i=0; i<m_lhs[n].size(); i++)
	memcpy(m_lhs[n][i], m_rhs[n][i], m_sz[n][i]);
  }
  
  inline virtual void wrecv(void)       { directXfer(); }
  inline virtual void mrecv(int wn)     { directXfer(); }
  inline virtual void transfer2m(int n) { directXfer(); }

};


/*
 * Class that implements the shared memory protocol but uses buffers.
 * 
 * There is no reason to use this class other than to test buffered
 * transfers on their own without the transfer protocol overhead
 */
class TransferBufferShMemBuf : public TransferBufferBuf
{
public:
  TransferBufferShMemBuf(int nwin):TransferBufferBuf(nwin)
  {
    /* Nothing specific here */
  }

  /* No sends */
  virtual void send2w(int wn) { }
  virtual void send2m(void)   { }

  inline void directXfer(void) {
    for (int n=1; n <= m_nwin; n++)
      bufIn(n); // and send
    
    // recv and ...
    for (int n=1; n <= m_nwin; n++)
      bufOut(n);
  }

  inline virtual void wrecv(void)       { directXfer(); }
  inline virtual void mrecv(int wn)     { directXfer(); }
  inline virtual void transfer2m(int n) { directXfer(); }

};


#ifdef HAVE_MPI
/*
 * Class that implements the two-sided communication MPI protocol
 *
 * Implememted using buffers
 */
class TransferBufferMPI : public TransferBufferBuf
{
public:
  TransferBufferMPI(int nwin, int rank):TransferBufferBuf(nwin), m_rank(rank)
  {
    
  }

  /*
   * Send this processes data to the master
   */
  inline virtual void send2m(void) {
    MPI_Request req;
    int tag, dst;
    const int one = 1;

    if (m_rank) {
      /* Fill contiguous buffer from win data */
      bufIn(one);

      dst = 0;
      tag = m_rank;

      /* Non-blocking */
      MPI_Isend(m_bufPtrs[one], m_bufSzs[one], MPI_CHAR, dst, tag, MPI_COMM_WORLD, &req);
      MPI_Request_free(&req);
    }
  }

  /*
   * Master receive from a window
   */
  inline virtual void mrecv(int wn) {
    MPI_Status  status;
    int tag, src;
    const int one = 1;

    if (wn != m_rank + 1) {
      src = wn-1;
      tag = wn-1;
      
      /* Blocking */
      MPI_Recv(m_bufPtrs[one], m_bufSzs[one], MPI_CHAR, src, tag, MPI_COMM_WORLD, &status);
      // {int foo; char buf[MAXSTRLEN]; MPI_Error_string(status.MPI_ERROR, buf, &foo); }

      /* Empty contiguous buffer into master */
      bufOut(one);
    }
  }

  /*
   * Send this processes data to another window
   */
  inline virtual void send2w(int wn) {
    MPI_Request req;
    int tag, dst;
    int self = m_rank+1;

    if (wn != self && m_bufSzs[self] > 0) {
      /* Fill contiguous buffer from win data */
      bufIn(self);

      dst = wn-1;
      tag = m_rank;

      /* Non-blocking */
      MPI_Isend(m_bufPtrs[self], m_bufSzs[self], MPI_CHAR, dst, tag, MPI_COMM_WORLD, &req);
      MPI_Request_free(&req);
    }
  }

  /*
   * Window receive from all other windows
   */
  inline virtual void wrecv(void) {
    int n;
    int self = m_rank+1;

    /* Do all receives */
#pragma omp parallel for
    for (n=1; n <= m_nwin; n++) {
      MPI_Status status;
      int tag, src;
      /* 
       * Skip self - although strictly speaking it would be a no-op
       *             anyway due to its zero size 
       */
      if (n == self || m_bufSzs[n] == 0) continue;
      
      src = n-1;
      tag = n-1;

      /* Blocking */
      MPI_Recv(m_bufPtrs[n], m_bufSzs[n], MPI_CHAR, src, tag, MPI_COMM_WORLD, &status);
      // {int foo; char buf[MAXSTRLEN]; MPI_Error_string(status.MPI_ERROR, buf, &foo); }

      /* Empty contiguous buffer into win data */
      bufOut(n); 
    }
  }

  /*
   * Window transfer to master - send/recv pair
   *
   * Used for master_fill
   */
  inline void transfer2m(int wn) { 
    MPI_Status  status;
    MPI_Request req;
    int tag, src, dst;
    int self = m_rank+1;
    const int one = 1;

    /* All non-master processes send */
    if (m_rank && (wn == self)) {
      /* Fill contiguous buffer from win data */
      bufIn(one);
      
      dst = 0;
      tag = m_rank;
    
      /* Non-blocking send */
      MPI_Isend(m_bufPtrs[one], m_bufSzs[one], MPI_CHAR, dst, tag, MPI_COMM_WORLD, &req);
      MPI_Request_free(&req);
    }

    /* Only the master receives */
    if (m_rank == 0) {
      /* This window */
      if (wn == 1) {
	/* Do direct transfer */
	bufIn(one);
	bufOut(one);
      } else {
	/* All others */
	src = wn - 1;
	tag = wn - 1;
	
	/* Blocking */
	MPI_Recv(m_bufPtrs[one], m_bufSzs[one], MPI_CHAR, src, tag, MPI_COMM_WORLD, &status);

	/* Empty into master */
	bufOut(one);
      }
    }
  }

private:
  /* MPI rank or ID of this process */
  const int m_rank;

};


/*
 * Class that implements the one-sided communication MPI protocol
 *
 * Also known as Remote Memory Access
 *
 * Recommended to use v3.x
 */
class TransferBufferMPI_RMA : public TransferBufferBuf
{
  /* WIP */
  inline void send2w(int wn) { }
  inline void send2m(void)   { }
  inline void wrecv(void)    { }
  inline void mrecv(int wn)  { }
  inline void transfer2m(int n) { }

};
#endif // HAVE_MPI


/*
 * Instantiate the appropriate TransferBuffer object
 */
static TransferBuffer *create_transfer_obj(master_t *master, int single = 0)
{
  /* eg. master_fill only requires one buffer */
  int num = (single ? 1 : master->nwindows);

#ifdef HAVE_MPI

  if (master->dp_mode & DP_MPI)
    return new TransferBufferMPI(num, master->mpi_rank);

#endif

  /* Check level of buffering to use */
  if (master->dp_mode & DP_BUFFEREDv2)
    return new TransferBufferShMemBuf(num);
  
  if (master->dp_mode & DP_BUFFEREDv1)
    return new TransferBufferShMem(num);
      
  /* Should never get here */
  hd_quit("TransferBuffer:create_transfer_obj, unknown DP_MODE of %d supplied\n", master->dp_mode);
  
  return(NULL);
}

/* Helper typedef used in container */
typedef map <int, TransferBuffer *> MultiBuf;

/*
 * Container class to hold the various transfer buffer object pointers   
 */
class TransferBufferObjs
{
public:
  /* Destructor */
  ~TransferBufferObjs(void) {
    delete fill_3d;
    delete fill_2d;
    delete master_fill;
    delete master_fill_ts;
    cleanup(refill_3d);
    cleanup(refill_2d);
    cleanup(empty_3d);
    cleanup(empty_2d);
  }
  void cleanup(MultiBuf &mb) {
    MultiBuf::iterator pos;
    for (pos = mb.begin(); pos != mb.end(); pos++)
      delete pos->second;
  }

  void dumpSizes(int rank, int win) {
    FILE *f = NULL;
    MultiBuf::iterator pos;
    char buf[MAXSTRLEN];

    if (WRITE_BUF_SIZES && rank == 0) {
      f = fopen("buffers.txt", "a");
      
      fprintf(f, "\nFor window :%d\n", win);
      master_fill->printSizes(f, "master_fill");
      master_fill_ts->printSizes(f, "master_fill_ts");
      fill_3d->printSizes(f, "fill_3d");
      fill_2d->printSizes(f, "fill_2d");

      for (pos = refill_3d.begin(); pos != refill_3d.end(); pos++) {
	sprintf(buf, "refill_3d[%d]", pos->first);
	pos->second->printSizes(f, buf);
      }

      for (pos = refill_2d.begin(); pos != refill_2d.end(); pos++) {
	sprintf(buf, "refill_2d[%d]", pos->first);
	pos->second->printSizes(f, buf);
      }

      for (pos = empty_3d.begin(); pos != empty_3d.end(); pos++) {
	sprintf(buf, "empty_3d[%d]", pos->first);
	pos->second->printSizes(f, buf);
      }

      fclose(f);
    }
  }
  
  /* One each */
  TransferBuffer *fill_3d;
  TransferBuffer *fill_2d;
  TransferBuffer *master_fill;
  TransferBuffer *master_fill_ts;
  /* These have multiple objects */
  MultiBuf refill_3d;
  MultiBuf refill_2d;
  MultiBuf empty_3d;
  MultiBuf empty_2d;
};


/*
 * Local functions to initialise transfer buffers
 */
static TransferBuffer *win_data_fill_3d_init(master_t *master,  geometry_t **window,
					     window_t **windat, int n)
{
  int c, cc, lc;       /* Local sparse coordinate / counter */
  int tn;              /* Tracer counter                    */
  int llc, wwn;

  /* Create objects */
  TransferBuffer *obj = create_transfer_obj(master);

  /*-----------------------------------------------------------------*/
  /* Auxiliary cells only. This includes variables that require      */
  /* auxiliary cells only to be transfered to the slave since the    */
  /* wet cells already exist on the slave.                           */
  /* 3D variables:                                                   */
  for (cc = 1; cc <= window[n]->nm2s; cc++) {
    lc = window[n]->m2s[cc];
    c = window[n]->wsa[lc];

    wwn = master->geom->fm[c].wn;
    llc = master->geom->fm[c].sc;
    
    obj->add(&windat[n]->u1[lc], &windat[wwn]->u1[llc], wwn);
    obj->add(&windat[n]->u2[lc], &windat[wwn]->u2[llc], wwn);
    obj->add(&windat[n]->w[lc],  &windat[wwn]->w[llc], wwn);
    obj->add(&windat[n]->Kz[lc], &windat[wwn]->Kz[llc], wwn);
    obj->add(&windat[n]->Vz[lc], &windat[wwn]->Vz[llc], wwn);
    obj->add(&windat[n]->dens[lc], &windat[wwn]->dens[llc], wwn);
    for (tn = 0; tn < windat[n]->ntr; tn++) {
      obj->add(&windat[n]->tr_wc[tn][lc], &windat[wwn]->tr_wc[tn][llc], wwn);
    }
    obj->add(&windat[n]->u1b[lc], &windat[wwn]->u1b[llc], wwn);
    obj->add(&windat[n]->u2b[lc], &windat[wwn]->u2b[llc], wwn);
    /* Note : xxx zoom */
    obj->add(&windat[n]->u1flux3d[lc], &windat[wwn]->u1flux3d[llc], wwn);
    obj->add(&windat[n]->u2flux3d[lc], &windat[wwn]->u2flux3d[llc], wwn);
  }

  /*-----------------------------------------------------------------*/
  /* 2D variables:                                                   */
  for (cc = 1; cc <= window[n]->nm2sS; cc++) {
    lc = window[n]->m2s[cc];
    c = window[n]->wsa[lc];

    wwn = master->geom->fm[c].wn;
    llc = master->geom->fm[c].sc;

    obj->add(&windat[n]->wtop[lc], &windat[wwn]->wtop[llc], wwn);
    obj->add(&windat[n]->wbot[lc], &windat[wwn]->wbot[llc], wwn);
    obj->add(&windat[n]->u1bot[lc], &windat[wwn]->u1bot[llc], wwn);
    obj->add(&windat[n]->u2bot[lc], &windat[wwn]->u2bot[llc], wwn);
    /* The updated value of eta is required at the auxiliary cells */
    /* to set dzu1 and dzu2, used in the calculation of the */
    /* advection terms.  */
    obj->add(&windat[n]->eta[lc], &windat[wwn]->eta[llc],wwn);
  }

  /* Initialise buffer */
  obj->init();

  return(obj);

} /* end fill_3d_init */

static TransferBuffer *win_data_fill_2d_init(master_t *master,  geometry_t **window,
					     window_t **windat, int n)
{
  int c, cc, lc;       /* Local sparse coordinate / counter */
  int llc, wwn;

  /* Create objects */
  TransferBuffer *obj = create_transfer_obj(master);

  /*-----------------------------------------------------------------*/
  /* Auxiliary cells only. This includes variables that require */
  /* auxiliary cells only to be transfered to the slave since the */
  /* wet cells already exist on the slave.  */
  for (cc = 1; cc <= window[n]->nm2sS; cc++) {
    lc = window[n]->m2s[cc];
    c = window[n]->wsa[lc];
    
    wwn = master->geom->fm[c].wn;
    llc = master->geom->fm[c].sc;

    obj->add(&windat[n]->eta[lc], &windat[wwn]->eta[llc],wwn);
    obj->add(&windat[n]->detadt[lc], &windat[wwn]->detadt[llc],wwn);
    obj->add(&windat[n]->u1av[lc], &windat[wwn]->u1av[llc],wwn);
    obj->add(&windat[n]->u2av[lc], &windat[wwn]->u2av[llc],wwn);
    obj->add(&windat[n]->depth_e1[lc], &windat[wwn]->depth_e1[llc],wwn);
    obj->add(&windat[n]->depth_e2[lc], &windat[wwn]->depth_e2[llc],wwn);
    obj->add(&windat[n]->etab[lc], &windat[wwn]->etab[llc],wwn);
  }

  /* Initialise buffer */
  obj->init();

  return(obj);
}

static TransferBuffer *win_data_master_fill_ts_init(master_t *master, geometry_t *window,
						    window_t *windat, win_priv_t *wincon)
{
  int c, cc, lc;
  int tt, tn;
  
  /* Create objects */
  TransferBuffer *obj = create_transfer_obj(master, 1);

  for (cc = 1; cc <= window->ns2m_ts; cc++) {
    lc = window->s2m_ts[cc];
    c = window->wsa[lc];
    if (lc <= window->ewetS) {
      obj->add(&master->eta[c], &windat->eta[lc]);
      obj->add(&master->u1av[c], &windat->u1av[lc]);
      obj->add(&master->u2av[c], &windat->u2av[lc]);
      for (tt = 0; tt < master->ntrS; tt++) {
	obj->add(&master->tr_wcS[tt][c], &windat->tr_wcS[tt][lc]);
      }
    }
    obj->add(&master->u1[c], &windat->u1[lc]);
    obj->add(&master->u2[c], &windat->u2[lc]);
    for (tt = 0; tt < master->ntrmap_s2m_3d; tt++) {
      tn = master->trmap_s2m_3d[tt];
      obj->add(&master->tr_wc[tn][c], &windat->tr_wc[tn][lc]);
    }
  }

  /* Initialise buffer */
  obj->init();

  return(obj);
}

static TransferBuffer *win_data_master_fill_init(master_t *master, geometry_t *window,
						 window_t *windat, win_priv_t *wincon)
{
  int c, cb, cc, lc;
  int tt, tn, n, k = 0;
  int zse1 = (int)(window->zmfe1 / 2);
  int zse2 = (int)(window->zmfe2 / 2);

  /* Create objects */
  TransferBuffer *obj = create_transfer_obj(master, 1);

  /**********/
  /* S2M_2D */
  /**********/
  for (cc = 1; cc <= window->b2_t; cc++) {
    lc = window->w2_t[cc];
    c = window->wsa[lc];
    obj->add(&master->eta[c], &windat->eta[lc]);
    obj->add(&master->topz[c], &windat->topz[lc]);
    obj->add(&master->wtop[c], &windat->wtop[lc]);
    for (tn = 0; tn < master->ntrS; tn++)
      obj->add(&master->tr_wcS[tn][c], &windat->tr_wcS[tn][lc]);
    for (tn = 0; tn < windat->nsed; tn++) {
      for (k = 0; k < window->sednz; k++)
        obj->add(&master->tr_sed[tn][k][c], &windat->tr_sed[tn][k][lc]);
    }
  }
  
  /* e1 face centered variables */
  for (cc = 1; cc <= window->b2_e1; cc++) {
    lc = window->w2_e1[cc];
    cb = window->wsa[window->bot_e1[cc]];
    c = window->wsa[lc];
    if (zse1) c = geom->zse1[c];
    obj->add(&master->u1av[c], &windat->u1av[lc]);
  }
  
  /* e2 face centered variables */
  for (cc = 1; cc <= window->b2_e2; cc++) {
    lc = window->w2_e2[cc];
    cb = window->wsa[window->bot_e1[cc]];
    c = window->wsa[lc];
    if (zse2) c = geom->zse2[c];
    obj->add(&master->u2av[c], &windat->u2av[lc]);
  }

  /**********/
  /* S2M_3D */
  /**********/

  /* Grid centered variables */
  for (cc = 1; cc <= window->b3_t; cc++) {
    lc = window->w3_t[cc];
    c = window->wsa[lc];
    
    obj->add(&master->w[c], &windat->w[lc]);
    obj->add(&master->dens[c], &windat->dens[lc]);
    obj->add(&master->dens_0[c], &windat->dens_0[lc]);
    obj->add(&master->Vz[c], &windat->Vz[lc]);
    obj->add(&master->Kz[c], &windat->Kz[lc]);
    obj->add(&master->dz[c], &wincon->dz[lc]);
    obj->add(&master->u1vh[c], &wincon->u1vh[lc]);
    obj->add(&master->u2vh[c], &wincon->u2vh[lc]);
    for (tt = 0; tt < master->ntrmap_s2m_3d; tt++) {
      tn = master->trmap_s2m_3d[tt];
      obj->add(&master->tr_wc[tn][c], &windat->tr_wc[tn][lc]);
    }
  }
  
  /* e1 face centered variables */
  for (cc = 1; cc <= window->b3_e1; cc++) {
    lc = window->w3_e1[cc];
    c = window->wsa[lc];
    if (zse1) c = geom->zse1[c];
    obj->add(&master->u1[c], &windat->u1[lc]);
    obj->add(&master->dzu1[c], &windat->dzu1[lc]);
  }
  
  /* e2 face centered variables */
  for (cc = 1; cc <= window->b3_e2; cc++) {
    lc = window->w3_e2[cc];
    c = window->wsa[lc];
    if (zse2) c = geom->zse2[c];
    obj->add(&master->u2[c], &windat->u2[lc]);
    obj->add(&master->dzu2[c], &windat->dzu2[lc]);
  }

  /* Fill the R_EDGE and F_EDGE cell centres if required */
  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    if (open->ocodex & R_EDGE) {
      for (tt = 0; tt < master->ntrmap_s2m_3d; tt++) {
	tn = master->trmap_s2m_3d[tt];
	if (master->trinfo_3d[tn].type & E1VAR) {
	  for (cc = 1; cc <= open->no3_e1; cc++) {
	    lc = open->obc_e1[cc];
	    c = window->wsa[lc];
	    obj->add(&master->tr_wc[tn][c], &windat->tr_wc[tn][lc]);
	  }
	}
      }
    }
    if (open->ocodey & F_EDGE) {
      for (tt = 0; tt < master->ntrmap_s2m_3d; tt++) {
	tn = master->trmap_s2m_3d[tt];
	if (master->trinfo_3d[tn].type & E2VAR) {
	  for (cc = 1; cc <= open->no3_e2; cc++) {
	    lc = open->obc_e2[cc];
	    c = window->wsa[lc];
	    obj->add(&master->tr_wc[tn][c], &windat->tr_wc[tn][lc]);
	  }
	}
      }
    }
  }

  /* Initialise buffer */
  obj->init();

  return(obj);
}

static TransferBuffer *win_data_refill_3d_init(master_t *master,  geometry_t **window,
					       window_t **windat, int n, int mode)
{
  int c, cc, lc;       /* Local sparse coordinate / counter */
  int llc, wwn;

  /* Create objects */
  TransferBuffer *obj = create_transfer_obj(master);

  if (mode & MIXING) {
    /* 3D cells */
    for (cc = 1; cc <= window[n]->nm2s; cc++) {
      lc = window[n]->m2s[cc];
      c = window[n]->wsa[lc];

      wwn = master->geom->fm[c].wn;
      llc = master->geom->fm[c].sc;

      obj->add(&windat[n]->Kz[lc], &windat[wwn]->Kz[llc], wwn);
      obj->add(&windat[n]->Vz[lc], &windat[wwn]->Vz[llc], wwn);
      obj->add(&windat[n]->dzu1[lc], &windat[wwn]->dzu1[llc], wwn);
      obj->add(&windat[n]->dzu2[lc], &windat[wwn]->dzu2[llc], wwn);

      if (windat[n]->tke && windat[n]->diss) {
	obj->add(&windat[n]->tke[lc], &windat[wwn]->tke[llc], wwn);
	obj->add(&windat[n]->diss[lc], &windat[wwn]->diss[llc], wwn);
      }
      if (windat[n]->tke && windat[n]->omega) {
	obj->add(&windat[n]->tke[lc], &windat[wwn]->tke[llc], wwn);
	obj->add(&windat[n]->omega[lc], &windat[wwn]->omega[llc], wwn);
      }
      if (windat[n]->Q2 && windat[n]->Q2L) {
	obj->add(&windat[n]->Q2[lc], &windat[wwn]->Q2[llc], wwn);
	obj->add(&windat[n]->Q2L[lc], &windat[wwn]->Q2L[llc], wwn);
      }
      if (windat[n]->sdc) {
	obj->add(&windat[n]->sdc[lc], &windat[wwn]->sdc[llc], wwn);
      }
    }
    /* Surface cells */
    for (cc = 1; cc <= window[n]->nm2sS; cc++) {
      lc = window[n]->m2s[cc];
      c = window[n]->wsa[lc];
      
      wwn = master->geom->fm[c].wn;
      llc = master->geom->fm[c].sc;
      
      obj->add(&windat[n]->sur_e1[lc], &windat[wwn]->sur_e1[llc], wwn);
      obj->add(&windat[n]->sur_e2[lc], &windat[wwn]->sur_e2[llc], wwn);
    }
  }

  if (mode & VELOCITY) {
    for (cc = 1; cc <= window[n]->nm2s; cc++) {
      lc = window[n]->m2s[cc];
      c = window[n]->wsa[lc];
      
      wwn = master->geom->fm[c].wn;
      llc = master->geom->fm[c].sc;
      
      obj->add(&windat[n]->u1[lc], &windat[wwn]->u1[llc], wwn);
      obj->add(&windat[n]->u2[lc], &windat[wwn]->u2[llc], wwn);
      /* Cell thickness is required for tracer horizontal diffusion */
      obj->add(&windat[n]->dzu1[lc], &windat[wwn]->dzu1[llc], wwn);
      obj->add(&windat[n]->dzu2[lc], &windat[wwn]->dzu2[llc], wwn);
      // xxx zoom
      obj->add(&windat[n]->u1flux3d[lc], &windat[wwn]->u1flux3d[llc], wwn);
      obj->add(&windat[n]->u2flux3d[lc], &windat[wwn]->u2flux3d[llc], wwn);
    }
    /* 2D */
    for (cc = 1; cc <= window[n]->nm2sS; cc++) {
      lc = window[n]->m2s[cc];
      c = window[n]->wsa[lc];
      
      wwn = master->geom->fm[c].wn;
      llc = master->geom->fm[c].sc;
      
      obj->add(&windat[n]->eta[lc], &windat[wwn]->eta[llc],wwn);
    }
  }

  if (mode & TRACERS) {
    for (cc = 1; cc <= window[n]->nm2s; cc++) {
      lc = window[n]->m2s[cc];
      c = window[n]->wsa[lc];
      
      wwn = master->geom->fm[c].wn;
      llc = master->geom->fm[c].sc;
      
      obj->add(&windat[n]->dens[lc], &windat[wwn]->dens[llc], wwn);
    }
  }

  /* Initialise buffer */
  obj->init();
  
  return(obj);
}

static TransferBuffer *win_data_refill_2d_init(master_t *master,  geometry_t **window,
					       window_t **windat, int n, int mode)
{
  int c, cc, lc;       /* Local sparse coordinate / counter */
  int llc, wwn;

  /* Create objects */
  TransferBuffer *obj = create_transfer_obj(master);

  if (mode & NVELOCITY) {
    /* 2D cells */
    for (cc = 1; cc <= window[n]->nm2sS; cc++) {
      lc = window[n]->m2s[cc];
      c = window[n]->wsa[lc];
      
      wwn = master->geom->fm[c].wn;
      llc = master->geom->fm[c].sc;
      
      obj->add(&windat[n]->nu1av[lc], &windat[wwn]->nu1av[llc], wwn);
      obj->add(&windat[n]->nu2av[lc], &windat[wwn]->nu2av[llc], wwn);
      // xxx zoom
    }
  }
  if (mode & DEPTH) {
    /* 2D cells */
    for (cc = 1; cc <= window[n]->nm2sS; cc++) {
      lc = window[n]->m2s[cc];
      c = window[n]->wsa[lc];
      
      wwn = master->geom->fm[c].wn;
      llc = master->geom->fm[c].sc;
      
      obj->add(&windat[n]->depth_e1[lc], &windat[wwn]->depth_e1[llc], wwn);
      obj->add(&windat[n]->depth_e2[lc], &windat[wwn]->depth_e2[llc], wwn);
    }
  }
  if (mode & ETA_A) {
    /* 2D cells */
    for (cc = 1; cc <= window[n]->nm2sS; cc++) {
      lc = window[n]->m2s[cc];
      c = window[n]->wsa[lc];
      
      wwn = master->geom->fm[c].wn;
      llc = master->geom->fm[c].sc;
      
      obj->add(&windat[n]->eta[lc], &windat[wwn]->eta[llc], wwn);
    }
  }

  /* Initialise buffer */
  obj->init();

  return(obj);
}

/*
 * This function collects all remote windows fields into the master's
 */
static TransferBuffer *win_data_empty_3d_init(master_t *master,  geometry_t **window,
					      window_t **windat, int n, int mode)
{
  /* Create objects */
  TransferBuffer *obj = create_transfer_obj(master, 1);

  if (mode & CFL) {
    if (!(master->cfl & NONE)) {
      obj->add(&windat[n]->mcfl2d, &windat[n]->mcfl2d);
      obj->add(&windat[n]->mcfl3d, &windat[n]->mcfl3d);
      obj->add(&windat[n]->cflc, &windat[n]->cflc);
    }
  }
  if (mode & TMASS) {
    /* Set the totals */
    // xxx similar to above
  }
  // xxx flags etc?
  
  if (mode & TRACERS) {
    obj->add(&windat[n]->wclk, &windat[n]->wclk);
  }
  
  /* Initialise buffer */
  obj->init();

  return(obj);
}


/****************/
/* PUBLIC API's */
/****************/
extern "C"
void win_data_buffered_recv_fill_3d(master_t *master, geometry_t *window, 
				    window_t *windat, int nwindows)
{
  TransferBufferObjs *obj = (TransferBufferObjs *)windat->win_data_buffered_obj;

  /* Update slave from master */
  //TIMING_SET;
  win_data_slave_update_fill_3d(master, window, windat, nwindows);
  //TIMING_DUMP(2, "   slave_update_3d");

#ifdef HAVE_MPI
  //TIMING_SET;
  //  MPI_Barrier(MPI_COMM_WORLD);
  //TIMING_DUMP(2, "   Barrier fill_3d");
#endif  
  //  TIMING_SET;
  /* Do the slave-slave transfer */
  obj->fill_3d->wrecv();
  //  TIMING_DUMP(2, "   wrecv_3d");
  
}

extern "C"
void win_data_buffered_send_fill_3d(window_t *windat, int wn)
{
  TransferBufferObjs *obj = (TransferBufferObjs *)windat->win_data_buffered_obj;
  
  /* Do the slave-slave transfer */
  obj->fill_3d->send2w(wn);
  
}

extern "C"
void win_data_buffered_recv_fill_2d(master_t *master, geometry_t *window,
				    window_t *windat, int nwindows)
{
  TransferBufferObjs *obj = (TransferBufferObjs *)windat->win_data_buffered_obj;

  /* Update from master */
  win_data_slave_update_fill_2d(master, window, windat, nwindows);

  /* Do the slave-slave transfer */
  obj->fill_2d->wrecv();
  
}

extern "C"
void win_data_buffered_send_fill_2d(window_t *windat, int wn)
{
  TransferBufferObjs *obj = (TransferBufferObjs *)windat->win_data_buffered_obj;

  /* Do the slave-slave transfer */
  obj->fill_2d->send2w(wn);
  
}

extern "C"
void win_data_buffered_recv_refill_3d(master_t *master, geometry_t *window,
				      window_t *windat, int nwindows, int mode)
{
  TransferBufferObjs *obj = (TransferBufferObjs *)windat->win_data_buffered_obj;

  /* Update from master */
  win_data_slave_update_refill_3d(master, window, windat, nwindows, mode);

  /* Do the slave-slave transfer */
  obj->refill_3d[mode]->wrecv();
}

extern "C"
void win_data_buffered_send_refill_3d(window_t *windat, int wn, int mode)
{
  TransferBufferObjs *obj = (TransferBufferObjs *)windat->win_data_buffered_obj;

  /* Do the slave-slave transfer */
  obj->refill_3d[mode]->send2w(wn);
}

extern "C"
void win_data_buffered_recv_refill_2d(master_t *master, geometry_t *window,
				      window_t *windat, int nwindows, int mode)
{
  TransferBufferObjs *obj = (TransferBufferObjs *)windat->win_data_buffered_obj;

  /* Do the slave-slave transfer */
  obj->refill_2d[mode]->wrecv();

}

extern "C"
void win_data_buffered_send_refill_2d(window_t *windat, int wn, int mode)
{
  TransferBufferObjs *obj = (TransferBufferObjs *)windat->win_data_buffered_obj;

  /* Do the slave-slave transfer */
  obj->refill_2d[mode]->send2w(wn);

}

extern "C"
void win_data_buffered_none(master_t *master, geometry_t *window,
			    window_t *windat, int mode)
{
  /* null function */
}

extern "C"
void win_data_buffered_empty_2d(master_t *master, geometry_t *window,
				window_t *windat, int mode)
{
  if (mode & VELOCITY) {
    memcpy(master->trincS, windat->trincS, master->ntrS * sizeof(int));
  }
}


extern "C"
void win_data_buffered_empty_3d(master_t *master, geometry_t *window,
				window_t *windat, int mode)
{
  int lc, cc, c;
  int tn, tt;
  TransferBufferObjs *obj = (TransferBufferObjs *)windat->win_data_buffered_obj;
  /* Find avoids creating an implicit entry with operator[] */
  MultiBuf::iterator pos = obj->empty_3d.find(mode);

  /* Do the slave-master transfer */
  if (pos != obj->empty_3d.end()) {
    pos->second->send2m();
  }

  if (mode & TRACERS) {
    /*
     * The following fills the local master with this windows data, the
     * whole domain goes through relaxation and resetting and then data
     * is copied back when win_data_slave_update_fill_3d runs
     */
    
    /*---------------------------------------------------------------*/
    /* All wet cells (cell centered). Include relaxation and reset   */
    /* tracers since tracers may be altered by updates. Include      */
    /* OBC cells.  */
    for (cc = 1; cc <= window->b3_t; cc++) {
      lc = window->w3_t[cc];
      c = window->wsa[lc];
      for (tt = 0; tt < master->nrlx; tt++) {
	tn = master->relax[tt];
	master->tr_wc[tn][c] = windat->tr_wc[tn][lc];
      }
      for (tt = 0; tt < master->nres; tt++) {
	tn = master->reset[tt];
	master->tr_wc[tn][c] = windat->tr_wc[tn][lc];
      }
    }
    
    memcpy(master->trinc, windat->trinc, master->ntr * sizeof(int));
    /* Mean iteration counter */
    if (master->means & RESET) master->means &= ~RESET;
    for (cc = 1; cc <= window->ns2m; cc++) {
      lc = window->s2m[cc];
      c = window->wsa[lc];
      master->dens[c] = windat->dens[lc];
      for (tt = 0; tt < master->ntrmap_s2m_3d; tt++) {
	tn = master->trmap_s2m_3d[tt];
	master->tr_wc[tn][c] = windat->tr_wc[tn][lc];
      }
    }
    if (!(master->means & NONE)) {
      for (cc = 1; cc <= window->enonS; cc++) {
	c = window->wsa[cc];
	master->meanc[c] = windat->meanc[cc];
      }
    }
  }

}

extern "C"
void win_data_buffered_recv_empty_3d(window_t *windat, int wn, int mode)
{
  TransferBufferObjs *obj = (TransferBufferObjs *)windat->win_data_buffered_obj;
  /* Find avoids creating an implicit entry with operator[] */
  MultiBuf::iterator pos = obj->empty_3d.find(mode);

  /* Do the slave-master transfer */
  if (pos != obj->empty_3d.end()) {
    pos->second->mrecv(wn);
  }
}

extern "C"
void win_data_buffered_master_fill(master_t *master,  geometry_t **window, 
				   window_t **windat, win_priv_t **wincon)
{

  /* Need to loop over all windows and collect */
#pragma omp parallel for
  for (int n=1; n<=master->nwindows; n++) {
    TransferBufferObjs *obj = (TransferBufferObjs *)windat[n]->win_data_buffered_obj;
    /* Do the slave-master transfer */
    obj->master_fill->transfer2m(n);
  }
}

extern "C"
void win_data_buffered_master_fill_ts(master_t *master,  geometry_t **window, 
				      window_t **windat, win_priv_t **wincon)
{
  TransferBufferObjs *obj;

  /* Need to loop over all windows and collect */
  for (int n=1; n<=master->nwindows; n++) {
    obj = (TransferBufferObjs *)windat[n]->win_data_buffered_obj;
    /* Do the slave-master transfer */
    obj->master_fill_ts->transfer2m(n);
  }
}

extern "C"
void win_data_buffered_update_master(master_t *master, window_t **windat, int mode)
{
  int n;

  /* The master process is filled with all the windows' data so we can do this */
  if (mode & TRACERS) {
    for (n=1; n<=master->nwindows; n++) {
      master->wclk[n] += windat[n]->wclk;
    }
  }
  
  if (mode & CFL) {
    if (!(master->cfl & NONE)) {
      for (n=1; n<=master->nwindows; n++) {
	if (windat[n]->mcfl2d < master->mcfl2d)
	  master->mcfl2d = windat[n]->mcfl2d;
	if (windat[n]->mcfl3d < master->mcfl3d) {
	  master->mcfl3d = windat[n]->mcfl3d;
	  master->cflc = windat[n]->cflc;
	}
      }
    }
  }
}

/*
 * Initialise objects for buffered transfers and setup the transfer
 * function pointers
 */
extern "C"
void win_data_init_transfer_buf_funcs(master_t *master,  geometry_t **window,
				      window_t **windat, win_priv_t **wincon,
				      int nwindows)
{
  int n;
  /*
   * Buffered mode : OPENMP or MPI
   *
   * Currently only supports slave-slaves transfers
   *
   */
  if (master->dp_mode & DP_BUFFERED) {
    /* Clear buffers file */
    if (WRITE_BUF_SIZES && master->mpi_rank == 0)
      unlink("buffers.txt");

    /* Each process will contain information for all windows */
    for (n=1; n<=nwindows; n++) {
      /* container object */
      TransferBufferObjs *obj = new TransferBufferObjs;
      /*
       * Create and cache buffer objects for each transfer mode
       */
      obj->fill_3d = win_data_fill_3d_init(master, window, windat, n);
      obj->fill_2d = win_data_fill_2d_init(master, window, windat, n);
      /* 3D refill modes */
      obj->refill_3d[MIXING]   = win_data_refill_3d_init(master, window, windat, n, MIXING);
      obj->refill_3d[VELOCITY] = win_data_refill_3d_init(master, window, windat, n, VELOCITY);
      /* 2D refill modes */
      obj->refill_2d[NVELOCITY] = win_data_refill_2d_init(master, window, windat, n, NVELOCITY);
      obj->refill_2d[DEPTH]     = win_data_refill_2d_init(master, window, windat, n, DEPTH);
      obj->refill_2d[ETA_A]     = win_data_refill_2d_init(master, window, windat, n, ETA_A);
      /* 3D empty to master */
      obj->empty_3d[TRACERS] = win_data_empty_3d_init(master, window, windat, n, TRACERS);
      obj->empty_3d[CFL]     = win_data_empty_3d_init(master, window, windat, n, CFL);
      /* Empty all to master */
      obj->master_fill = win_data_master_fill_init(master, window[n], windat[n], wincon[n]);
      /* Empty ts subset to master */
      obj->master_fill_ts = win_data_master_fill_ts_init(master, window[n], windat[n], wincon[n]);

      obj->dumpSizes(master->mpi_rank, n);

      windat[n]->win_data_buffered_obj = (void *)obj;
    }

    /* Setup function pointers */
    master->win_data_fill_3d   = win_data_buffered_recv_fill_3d;
    master->win_data_fill_2d   = win_data_buffered_recv_fill_2d;
    master->win_data_refill_3d = win_data_buffered_recv_refill_3d;
    master->win_data_refill_2d = win_data_buffered_recv_refill_2d;
    master->win_data_empty_3d  = win_data_buffered_empty_3d;
    master->win_data_empty_2d  = win_data_buffered_empty_2d;
    master->master_fill_ts = win_data_buffered_master_fill_ts;
    master->master_fill    = win_data_buffered_master_fill;
    // 
    master->update_master = win_data_buffered_update_master;

    // if (master) print sizes

  } else {
    /*
     * Classic mode: NONE PTHREADS OPENMP
     *
     * Wire-up the original master-slave transfer routines
     */
    master->win_data_fill_3d   = win_data_fill_3d;
    master->win_data_fill_2d   = win_data_fill_2d;
    master->win_data_refill_3d = win_data_refill_3d;
    master->win_data_refill_2d = win_data_refill_2d;
    master->win_data_empty_3d  = win_data_empty_3d;
    master->win_data_empty_2d  = win_data_empty_2d;
    // No need to update master - the empty functions handle this
    master->update_master  = NULL;
    master->master_fill_ts = master_fill_ts;
    master->master_fill    = master_fill;
  }
}  

/*
 * Free memory
 */
extern "C"
void win_data_init_transfer_buf_cleanup(window_t *windat)
{
  if (master->dp_mode & DP_BUFFERED)
    delete ((TransferBufferObjs *)windat->win_data_buffered_obj);
}
