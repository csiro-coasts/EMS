/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/da/da_obj.cpp
 *  
 *  Description:
 *  Data assimilation classes
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: da_obj.cpp 5850 2018-06-29 05:18:59Z riz008 $
 *
 */

/*
 * TODO:
 *  o) Error function
 *  o) Make this lib optional
 *
 *  o) if num for depth else variable as well as lat/longs
 *  o) cleanup memory - error conditions
 *  o) Introduce children for da_obs; need this for the glider
 */
#include <list>
#include <vector>
#include <iterator>
#include <string>
using namespace std;

#ifdef USE_CUBLAS
/* Includes, cuda */
#include <cuda_runtime.h>
#include <cublas_v2.h>
#endif

#include <sys/time.h>
#include <iostream>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include "ems.h"
#include "da_utils.hpp"
#include "da_obj.hpp"

/**
 * Constructor
 */
da_obj::da_obj(double ls)
{
  f_A    = NULL;
  f_R    = NULL;
  f_HA   = NULL;
  f_K    = NULL;
  f_wo   = NULL;
  f_wb   = NULL;
  f_wa   = NULL;
  f_Hwb  = NULL;
  f_maps = NULL;
  f_x    = NULL;
  f_y    = NULL;

  f_full_wo_val = NULL;
  f_full_wo_err = NULL;
  
  f_ls = ls;
}


/**
 * Destructor
 */
da_obj::~da_obj(void)
{
  int i;

  if (f_A != NULL)
    gsl_matrix_free(f_A);

  if (f_R != NULL)
    gsl_matrix_free(f_R);

  if (f_HA != NULL)
    gsl_matrix_free(f_HA);

  if (f_K != NULL)
    gsl_matrix_free(f_K);

  if (f_wo != NULL)
    gsl_vector_free(f_wo);

  if (f_full_wo_val != NULL)
    gsl_vector_free(f_full_wo_val);

  if (f_full_wo_err != NULL)
    gsl_vector_free(f_full_wo_err);

  if (f_wb != NULL)
    gsl_vector_free(f_wb);

  if (f_x != NULL)
    gsl_vector_free(f_x);

  if (f_y != NULL)
    gsl_vector_free(f_y);

  if (f_Hwb != NULL)
    gsl_vector_free(f_Hwb);

  if (f_wa != NULL)
    gsl_vector_free(f_wa);

  if (f_maps != NULL)
    delete(f_maps);

  for (i=0; i<f_obs_arr.size(); i++)
    delete(f_obs_arr[i]);

  /*
   * f_obs_exec and its contents will be automatically deleted
   */
}


/**
 * Allocates the A matrix and sets all values to NaN
 * @param nrows number of rows
 * @param ncols number of columns
 * @return 0 on success, 1 failure
 */
int da_obj::allocA(int nrows, int ncols)
{
  f_A = gsl_matrix_alloc((size_t)nrows, (size_t)ncols);
  
  if (f_A == NULL)
    return(1);

  gsl_matrix_set_all(f_A, DA_INVALID_VAL);

  return(0);
}


/**
 * Add the state name and offset
 */
void da_obj::add_state(const char *str, int offset)
{
  f_states[str] = offset;
}


/**
 * Gets the offset value from the map
 */
int da_obj::get_state_offset(char *ref_var)
{
  map<string, int>::iterator pos;
  if ( (pos = f_states.find(ref_var)) == f_states.end() )
    quit("DA state '%s' does not exist\n", ref_var);

  return(pos->second);
}


/**
 * Fill one element of A
 */
void da_obj::fillA(int row, int col, double val)
{
  gsl_matrix_set(f_A, (size_t)row, (size_t)col, val);
}

/**
 * Fill one column of the A matrix
 */
void da_obj::fillA(int col, double *vals)
{
  size_t r;
  for (r=0; r<f_A->size1; r++)
    gsl_matrix_set(f_A, r, (size_t)col, vals[r]);
}

/**
 * Check to make sure all values are finite
 * We could just grab the first pointer and loop over the whole array
 * but not sure if that'll work with the tda concept in gsl
 */
int da_obj::checkA(void)
{
  int r,c;
  for (r=0; r<f_A->size1; r++) {
    for (c=0; c<f_A->size2; c++) {
      double val = gsl_matrix_get(f_A, (size_t)r, (size_t)c);
      if (isnan(val) || val == DA_INVALID_VAL)
	return(0);
    }
  }
  
  return(1);
}


/**
 * Writes out the given matrix to a plain text file. This is useful for
 * debugging purposes
 */
void da_obj::dumpMatrix(const char *fname, const gsl_matrix *M)
{
  size_t r,c;
  FILE *F = fopen(fname, "w");
  
  for (r=0; r<M->size1; r++) {
    for (c=0; c<M->size2; c++) {
      double val = gsl_matrix_get(M, r, c);
      fprintf(F, "%f ", val);
    }
    fprintf(F, "\n");
  }

  fclose(F);
  sync();
}

/**
 * Wrapper for A
 */
void da_obj::writeA(const char *fname)
{
  dumpMatrix(fname, f_A);
}

/**
 * Adds the observation object to the array
 */
void da_obj::add_obs(da_obs *obs)
{
  f_obs_arr.push_back(obs);
}


/**
 * Loops through to read all of the observations and sticks them in an
 * observations execution list
 */
int da_obj::read_all_obs(double t)
{
  int i;
  
  /* Clear out any previous observations */
  f_obs_exec.clear();

  for (i=0; i<f_obs_arr.size(); i++) {
    /* Splice the list for this observation into the master */
    f_obs_exec.splice(f_obs_exec.end(), f_obs_arr[i]->read_obs(t, f_maps));
  }

  return(f_obs_exec.size());
}


/**
 * Constructs the R matrix
 */
void da_obj::constructR(void)
{
  size_t nobs = f_obs_exec.size();
  list <da_obs_exec>::iterator obs;
  int i = 0;

  /* Clear old memory */
  if (f_R != NULL)
    gsl_matrix_free(f_R);

  /* Allocate memory for the matrix */
  f_R = gsl_matrix_calloc(nobs, nobs);
  if (f_R == NULL)
    quit("Unable to allocate memory for the R matrix\n");

  /* Loop over the observations */
  for (obs = f_obs_exec.begin(); obs != f_obs_exec.end(); obs++) {
    gsl_matrix_set(f_R, i, i, (*obs).f_err);
    i++;
  }
}


/**
 * Constructs the HA matrix
 */
void da_obj::constructHA_Hwb(void)
{
  size_t nobs     = f_obs_exec.size();
  size_t nmembers = f_A->size2;
  list <da_obs_exec>::iterator obs;
  int r = 0;

  /*
   * Make sure wb is set
   */
  if (f_wb == NULL)
    quit("wb has not been populated\n");

  /* 
   * Clear old memory
   */
  if (f_HA != NULL)
    gsl_matrix_free(f_HA);
  
  /* Allocate memory for the matrix */
  f_HA = gsl_matrix_calloc(nobs, nmembers);
  if (f_HA == NULL)
    quit("Unable to allocate memory for the HA matrix\n");


  /* Clear old and allocate new */
  if (f_Hwb != NULL)
    gsl_vector_free(f_Hwb);
  f_Hwb = gsl_vector_calloc(nobs);
  if (f_Hwb == NULL)
    quit("Unable to allocate Hwb vector\n");

  /* Loop over the observations */
  for (obs = f_obs_exec.begin(); obs != f_obs_exec.end(); obs++) {
    int row = (*obs).f_gs;

    /* Quick check */
    if (row >= f_A->size1)
      quit("Invalid row observation calculation\n");
    
    /* Go ahead and assign */
    {
      gsl_vector_const_view rowH = gsl_matrix_const_row(f_A, row);
      gsl_matrix_set_row(f_HA, r, &rowH.vector);
      gsl_vector_set(f_Hwb, r, gsl_vector_get(f_wb, row));
    }

    /* Increment counter */
    r++;
  }
}

/**
 * Applies localisation, overwriting the 
 */
void da_obj::applyLocalisaton(gsl_matrix *AHA)
{
  size_t len  = AHA->size1;
  gsl_vector *lf  = gsl_vector_calloc(len);
  gsl_vector *col = gsl_vector_calloc(len);
  list <da_obs_exec>::iterator obs;
  int i,j = 0;
  /*
   * Matlab code from EJ
   * distance = sqrt((x_centre-lon_obs).^2+(y_centre-lat_obs).^2);
   * lf_2d = 1./(1+0.5.*(distance./LS).^2); %localisation factor
   */
  //FILE *f = fopen("local.txt", "w");
  //int done = 0;
  for (obs = f_obs_exec.begin(); obs != f_obs_exec.end(); obs++,j++) {
    for (i=0; i<len; i++) {
      /* Get the differences */
      double dx = gsl_vector_get(f_x, i) - obs->f_x;
      double dy = gsl_vector_get(f_y, i) - obs->f_y;
      /* Calculate distance */
      double dist = gsl_hypot(dx, dy);
      /* Calculate localisation factor */
      double fact = 1.0 / (1 + 0.5* gsl_pow_2(dist / f_ls));
      /* Set factor */
      gsl_vector_set(lf, i, fact);
      // xxx
      /*
	{
	int I,J,K;
	f_maps->getIJK(i,  &I, &J, &K);
	// write surface only
	if (K == 46 && !done) {
	fprintf(f, "dist(%d,%d,%d)=%f;\n", I+1, J+1, K+1, dist);
	fprintf(f, "fact(%d,%d,%d)=%f;\n", I+1, J+1, K+1, fact);
	}
	}
      */
    }
    /* Apply and replace column */
    gsl_matrix_get_col(col, AHA, j);
    gsl_vector_mul(lf, col);
    gsl_matrix_set_col (AHA, j, lf);
    /*
      for (i=0; i<len; i++) {
      int I,J,K;
      f_maps->getIJK(i,  &I, &J, &K);
      // write surface only
      if (K == 46 && !done) {
      fprintf(f, "P0(%d,%d,%d)=%f;\n", I+1, J+1, K+1, gsl_vector_get(col, i));
      fprintf(f, "P1(%d,%d,%d)=%f;\n", I+1, J+1, K+1, gsl_vector_get(lf,  i));
      }
      }
      done = 1;
    */
  }
  //fclose(f);
  /* Clean up */
  gsl_vector_free(lf);
  gsl_vector_free(col);
}


#ifdef USE_CUBLAS
static void doCUBLAS(gsl_matrix *AHA, gsl_matrix *P, gsl_matrix *f_K)
{
  double *d_A, *d_B, *d_C;
  cublasStatus_t status;
  cublasHandle_t handle;
  double alpha = 1.0;
  double beta  = 0.0;
  struct timeval tm1, tm2;

  printf("Start GPU\n");
  /* Initialise CUBLAS */
  gettimeofday(&tm1, NULL);
  status = cublasCreate(&handle);
  if (status != CUBLAS_STATUS_SUCCESS)
    fprintf (stderr, "!!!! CUBLAS initialization error\n");
  gettimeofday(&tm2, NULL);
  printf(" start cublas %.5f\n", 
	 ((tm2.tv_sec  - tm1.tv_sec)+(tm2.tv_usec - tm1.tv_usec)*1e-6));
  
  /* Allocate device memory for the matrices */
  gettimeofday(&tm1, NULL);
  if (cudaMalloc((void**)&d_A, AHA->size1*AHA->size2 * sizeof(double)) != cudaSuccess)
    fprintf (stderr, "!!!! device memory allocation error (allocate A)\n");
  if (cudaMalloc((void**)&d_B, P->size1*P->size2 * sizeof(double)) != cudaSuccess)
    fprintf (stderr, "!!!! device memory allocation error (allocate B)\n");
  if (cudaMalloc((void**)&d_C, f_K->size1*f_K->size2 * sizeof(double)) != cudaSuccess)
    fprintf (stderr, "!!!! device memory allocation error (allocate C)\n");
  gettimeofday(&tm2, NULL);
  printf(" allocate d_ arrays%.5f\n", 
	 ((tm2.tv_sec  - tm1.tv_sec)+(tm2.tv_usec - tm1.tv_usec)*1e-6));
  
  /* Initialize the device matrices with the host matrices */
  gettimeofday(&tm1, NULL);
  status = cublasSetMatrix(AHA->size2, AHA->size1,sizeof(double),AHA->data, AHA->tda, d_A, AHA->size2);
  if (status != CUBLAS_STATUS_SUCCESS)
    fprintf (stderr, "!!!! device access error (write A)\n");
  status = cublasSetMatrix(P->size2, P->size1,sizeof(double),P->data, P->tda, d_B, P->size2);
  if (status != CUBLAS_STATUS_SUCCESS)
    fprintf (stderr, "!!!! device access error (write B)\n");
  gettimeofday(&tm2, NULL);
  printf(" cudaSetMatrices %.5f\n", 
	 ((tm2.tv_sec  - tm1.tv_sec)+(tm2.tv_usec - tm1.tv_usec)*1e-6));

  /* Perform the calcuation */
  gettimeofday(&tm1, NULL);
  status = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, P->size2, AHA->size1, P->size1, 
		       &alpha, d_B, P->size2, d_A, AHA->size2, &beta, d_C, f_K->size2);
  if (status != CUBLAS_STATUS_SUCCESS)
    fprintf (stderr, "!!!! kernel execution error.\n");
  cudaThreadSynchronize();
  gettimeofday(&tm2, NULL);
  printf(" dgemm %.5f\n", 
	 ((tm2.tv_sec  - tm1.tv_sec)+(tm2.tv_usec - tm1.tv_usec)*1e-6));
  
  /* Now copy the matrix back */
  gettimeofday(&tm1, NULL);
  status = cublasGetMatrix(f_K->size2, f_K->size1, sizeof(double), d_C,
			   f_K->size2, f_K->data, f_K->tda);
  if (status != CUBLAS_STATUS_SUCCESS)
    fprintf (stderr, "!!!! device access error (read C)\n");
  gettimeofday(&tm2, NULL);
  printf(" Get matrix %.5f\n", 
	 ((tm2.tv_sec  - tm1.tv_sec)+(tm2.tv_usec - tm1.tv_usec)*1e-6));

  /* Cleanup */
  gettimeofday(&tm1, NULL);
  if (cudaFree(d_A) != cudaSuccess)
    fprintf (stderr, "!!!! memory free error (A)\n");
  if (cudaFree(d_B) != cudaSuccess)
    fprintf (stderr, "!!!! memory free error (B)\n");
  if (cudaFree(d_C) != cudaSuccess)
    fprintf (stderr, "!!!! memory free error (C)\n");
  gettimeofday(&tm2, NULL);
  printf(" Cleanup %.5f\n", 
	 ((tm2.tv_sec  - tm1.tv_sec)+(tm2.tv_usec - tm1.tv_usec)*1e-6));

  /* Shutdown */
  gettimeofday(&tm1, NULL);
  status = cublasDestroy(handle);
  if (status != CUBLAS_STATUS_SUCCESS)
    fprintf (stderr, "!!!! shutdown error (A)\n");
  gettimeofday(&tm2, NULL);
  printf("Shutdown %.5f\n", 
	 ((tm2.tv_sec  - tm1.tv_sec)+(tm2.tv_usec - tm1.tv_usec)*1e-6));
}
#endif

/*
 * Calculate the K matrix using blas functions
 *
 *  K = (1/n-1)*A(HA)' * inv((1/n-1)(HA)(HA)' + R)
 *              ^^^^^         ^^^^^^^^^^^^^^
 *               AHA                P 
 */
void da_obj::constructK(void)
{
  /* Pointer for K and intermediate matrices */
  gsl_matrix *AHA, *P;

  /* various size elements */
  size_t rowsA;
  size_t nmem;
  size_t nobs;
  double i_nmem;
  int ret;
  size_t i,j;
  struct timeval tm1, tm2;

  /* Need R and HA to proceed */
  constructR();
  constructHA_Hwb();
  
  /* Get some sizes */
  rowsA  = f_A->size1;
  nmem   = f_A->size2;
  nobs   = f_HA->size1;
  i_nmem = 1.0 / (nmem-1);

  /* Allocate matrices */
  if (f_K != NULL)
    gsl_matrix_free(f_K);

  f_K = gsl_matrix_calloc(rowsA, nobs);
  if (f_K == NULL)
    quit("Error creating K matrix\n");

  P = gsl_matrix_calloc(nobs, nobs);
  if (P == NULL)
    quit("Error creating intermediate P matrix\n");

  AHA = gsl_matrix_calloc(rowsA, nobs);
  if (AHA == NULL)
    quit("Error creating intermediate AHA matrix\n");

  /* AHA */
  ret = gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, f_A, f_HA, 0.0, AHA);

  /* Modify AHA for localisation, if using */
  //  dumpMatrix("aha_before.txt", AHA);
  if (f_ls)
    applyLocalisaton(AHA);
  // dumpMatrix("aha_after.txt", AHA);

  /* P */
  ret = gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, f_HA, f_HA, 0.0, P);
  ret = gsl_matrix_scale(P, i_nmem);

  /* P = P + R */
  ret = gsl_matrix_add(P, f_R);

  /* P = inv(P) */
  ret = gsl_linalg_cholesky_decomp(P);
  ret = gsl_linalg_cholesky_invert(P);

  /* Finally construct K and then scale */
  // gettimeofday(&tm1, NULL);
#ifdef USE_CUBLAS
  doCUBLAS(AHA, P, f_K);
#else
  ret = gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, AHA, P, 0.0, f_K);
#endif
  //gettimeofday(&tm2, NULL);
  //printf("Took %.5f\n", 
  //((tm2.tv_sec  - tm1.tv_sec)+(tm2.tv_usec - tm1.tv_usec)*1e-6));
  //fflush(NULL);

  ret = gsl_matrix_scale(f_K, i_nmem);

  /* free up intermediate variables */
  gsl_matrix_free(AHA);
  gsl_matrix_free(P);
}


/**
 * Constructs the observation vector
 */
void da_obj::constructWo(void)
{
  size_t i = 0;
  size_t nobs = f_obs_exec.size();
  list <da_obs_exec>::iterator sta = f_obs_exec.begin();
  list <da_obs_exec>::iterator fin = f_obs_exec.end();

  /* Clear old data and allocate a new one */
  if (f_wo != NULL) {
    gsl_vector_free(f_wo);
    f_wo = NULL;
  }
  f_wo = gsl_vector_calloc(nobs);
  if (f_wo == NULL)
    quit("Unable to allocate observations vector\n");

  /* Loop through and fill values */
  while(sta != fin)
    gsl_vector_set(f_wo, i++, (*sta++).f_val);
}


/**
 * Constructs the full observations vector
 */
void da_obj::construct_full_wo_val(void)
{
  size_t i = 0;
  list <da_obs_exec>::iterator obs;

  /* Clear old data and allocate a new one */
  if (f_full_wo_val != NULL) {
    gsl_vector_free(f_full_wo_val);
    f_full_wo_val = NULL;
  }
  if (f_wa == NULL)
    quit("Cannot construct full wo before analysis\n");

  /* Allocate the same size as the analysis vector */
  f_full_wo_val = gsl_vector_calloc(f_wa->size);
  if (f_full_wo_val == NULL)
    quit("Unable to allocate full observations vector\n");

  /* Set them all to nan's */
  gsl_vector_set_all(f_full_wo_val, NaN);

  /* Loop through and fill values */
  for (obs = f_obs_exec.begin(); obs != f_obs_exec.end(); obs++)
    gsl_vector_set(f_full_wo_val, (*obs).f_gs, (*obs).f_val);

}


/**
 * Constructs the full observations vector of errors
 */
void da_obj::construct_full_wo_err(void)
{
  size_t i = 0;
  list <da_obs_exec>::iterator obs;

  /* Clear old data and allocate a new one */
  if (f_full_wo_err != NULL) {
    gsl_vector_free(f_full_wo_err);
    f_full_wo_err = NULL;
  }
  if (f_wa == NULL)
    quit("Cannot construct full wo before analysis\n");

  /* Allocate the same size as the analysis vector */
  f_full_wo_err = gsl_vector_calloc(f_wa->size);
  if (f_full_wo_err == NULL)
    quit("Unable to allocate full observations vector\n");

  /* Set them all to nan's */
  gsl_vector_set_all(f_full_wo_err, NaN);

  /* Loop through and fill values */
  for (obs = f_obs_exec.begin(); obs != f_obs_exec.end(); obs++)
    gsl_vector_set(f_full_wo_err, (*obs).f_gs, (*obs).f_err);

}


/**
 * Sets up the full state vector
 */
void da_obj::set_wb(int index, double data)
{
  size_t len = f_A->size1;
  
  /* Only need to allocate once */
  if (f_wb == NULL) {
    f_wb = gsl_vector_calloc(len);
  }
  
  /* Quick check */
  if (index >= len)
    quit("Index (%d) out of range for state vector\n", index);
  
  /* All good, go ahead and set value */
  gsl_vector_set(f_wb, index, data);

}

/**
 * Fills in the location vectors
 */
void da_obj::set_xy(int index, double x, double y)
{
  size_t len = f_A->size1;
  
  /* Only need to allocate once */
  if (f_x == NULL) {
    f_x = gsl_vector_calloc(len);
    f_y = gsl_vector_calloc(len);
  }

  /* Quick check */
  if (index >= len)
    quit("Index (%d) out of range for state vector\n", index);
  
  /* All good, go ahead and set value */
  gsl_vector_set(f_x, index, x);
  gsl_vector_set(f_y, index, y);
}

/**
 * Calculate analysis field
 *
 *    wa = wb + K(wo - Hwb)
 *              ^^^^^^^^^^
 *                  y : wa - intermediate
 */
void da_obj::do_analysis(void)
{
  size_t len;
  size_t nobs;
  size_t n;
  list <da_obs_exec>::iterator sta = f_obs_exec.begin();
  list <da_obs_exec>::iterator fin = f_obs_exec.end();
  gsl_vector *wo_minus_Hwb;

  /* Need K and the observations vector to proceed */
  constructK();
  constructWo();

  len  = f_A->size1;
  nobs = f_K->size2;
  
  /* Allocate vectors, if needed */
  if (f_wa == NULL) {
    f_wa = gsl_vector_calloc(len);
    if (f_wa == NULL)
      quit("Unable to allocate analysis vector\n");
  }

  wo_minus_Hwb = gsl_vector_calloc(nobs);
  if (wo_minus_Hwb == NULL)
    quit("Unable to allocate intermediate analysis vector\n");

  /* do the difference */
  n = 0;
  while (sta != fin) {
    double a = (*sta).f_val;
    double b = gsl_vector_get(f_Hwb, n);
    gsl_vector_set(wo_minus_Hwb, n, (a-b));
    n++;
    sta++;
  }

  /* Construct y but store in wa*/
  gsl_blas_dgemv(CblasNoTrans, 1.0, f_K, wo_minus_Hwb, 0.0, f_wa);
     
  /* Finally sum up wb and y(wa) to form wa */
  gsl_vector_add(f_wa, f_wb);

  /* Free intermediate mem */
  gsl_vector_free(wo_minus_Hwb);
}


/** Copies background into analysis
 *
 */
void da_obj::copy_wb_wa(void)
{
  int len = f_A->size1;
  /* Allocate vectors, if needed */
  if (f_wa == NULL) {
    f_wa = gsl_vector_calloc(len);
    if (f_wa == NULL)
      quit("Unable to allocate analysis vector\n");
  }
  /* This should've been allocated and set */
  if (f_wb == NULL) {
    quit("Background vector not set\n");
  }
  /* Copy over */
  gsl_vector_memcpy(f_wa, f_wb);
}


/** Gets a single value of the analysis field. Opposite of set_wb
 *
 */
double da_obj::get_wa(int index)
{
  return(gsl_vector_get(f_wa, index));
}


/** Gets a single value of the full observations vector
 *
 */
double da_obj::get_wo_val(int index)
{
  if (index == 0)
    construct_full_wo_val();

  return(gsl_vector_get(f_full_wo_val, index));
}


/** Gets a single value of the full observations vector
 *
 */
double da_obj::get_wo_err(int index)
{
  if (index == 0)
    construct_full_wo_err();

  return(gsl_vector_get(f_full_wo_err, index));
}


/* end da_obj.cpp */
