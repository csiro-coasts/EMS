/*
 * DA interface to Intels MKL
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: da_mkl.cpp 5850 2018-06-29 05:18:59Z riz008 $
 */

#ifdef USE_MKL
#include <mkl_cblas.h>
void call_mkl(int A1, int A2, double *A,
	      int B1, int B2, double *B,
	      int C1, int C2, double *C)
{
  cblas_dgemm(CblasRowMajor,  /* const  CBLAS_ORDER Order */
	      CblasNoTrans,   /* const  CBLAS_TRANSPOSE TransA */
	      CblasNoTrans,   /* const  CBLAS_TRANSPOSE TransB */
	      A1,      /* const MKL_INT M */
	      B2,       /* const MKL_INT N */
	      A2,     /* const MKL_INT K */
	      1.0,            /* const double alpha */
	      A,        /* const double *A */
	      A1, /* const MKL_INT lda */
	      B,          /* const double *B */
	      B1,   /* const MKL_INT ldb */
	      0.0,            /* const double beta */
	      C,          /* double *C */
	      A1);    /* const MKL_INT ldc */
}
#endif

