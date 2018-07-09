/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/misc/memory.c
 *
 *  \brief Array memory allocation
 *
 *  Routines for allocation of 1, 2 and 3
 *  dimensional arrays of various types
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: memory.c 5831 2018-06-26 23:48:06Z riz008 $
 */

#include <stdio.h>
#include <stdlib.h>
#include "ems.h"
#include <string.h>
#include <assert.h>


/** Allocate and clear a 1d array of double values.
  *
  * @param n1 size of array.
  * @return an array n1 of double value.
  */
double *d_alloc_1d(int n1)
{
  return (double *)alloc_1d(n1,sizeof(double));
}

/** Deallocates the memory for a 1d array of doubles.
  * Was allocated using d_alloc1d.
  *
  * @param p pointer to aray of doubles.
  */
void d_free_1d(double *p)
{
  free_1d(p);
}

/** Allocate and clear a 2d array. Note that the storage for the values
  * is one large piece of memory so that it can be read or written
  * easily. The i,j value is accessed by saying array[j][i], because
  * array[j] is a pointer to the jth row of the matrix.
  * n1 is the number of columns, n2 is the number of rows.
  *
  * @param n1 size of first dimension of the array.
  * @param n2 size of second dimension of the array.
  * @return a 2d array n1*n2 double values.
  */
double **d_alloc_2d(int n1, int n2)
{
  return (double **)alloc_2d(n1,n2,sizeof(double));
}

/** Deallocates the memory for a 2d array of doubles.
  * Was allocated using d_alloc2d.
  *
  * @param pp pointer to aray of doubles.
  */
void d_free_2d(double **pp)
{
  free_2d(pp);
}


/** Allocate a 3d array. Note that the storage for the values
  * is one large piece of memory so that it can be read or written
  * easily. The i,j,k value is accessed by saying array[k][j][i].
  * array[k] is a pointer to the kth plane of the array. array[k][j]
  * is a pointer to the jth row in the kth plane of the array.
  *
  * @param n1 size of first dimension of the array.
  * @param n2 size of second dimension of the array.
  * @param n3 size of third dimension of the array.
  * @return a 3d array n1*n2*n3 double values.
  */
double ***d_alloc_3d(int n1, int n2, int n3)
{
  return (double ***)alloc_3d(n1,n2,n3,sizeof(double));
}

/** Deallocates the memory for a 3d array of doubles.
  * Was allocated using d_alloc3d.
  *
  * @param ppp pointer to aray of doubles.
  */
void d_free_3d(double ***ppp)
{	
  free_3d(ppp);
}

/** Allocate a 4d array. Note that the storage for the values
  * is one large piece of memory so that it can be read or written
  * easily. The i,j,k,l value is accessed by saying array[l][k][j][i].
  *
  * @param n1 size of first dimension of the array.
  * @param n2 size of second dimension of the array.
  * @param n3 size of third dimension of the array.
  * @param n4 size of third dimension of the array.
  * @return a 4d array n1*n2*n3*n4 double values.
  */
double ****d_alloc_4d(int n1, int n2, int n3, int n4)
{
  return alloc_4d(n1, n2, n3, n4, sizeof(double));
}

/** Deallocates the memory for a 4d array of doubles.
  * Was allocated using d_alloc4d.
  *
  * @param pppp pointer to aray of doubles.
  */
void d_free_4d(double ****pppp)
{
  free_4d(pppp);
}

/** Allocate and clear a 1d array of float values.
  *
  * @param n1 size of array.
  * @return an array n1 of float value.
  */
float *f_alloc_1d(int n1)
{
  return (float *)alloc_1d(n1,sizeof(float));
}


/** Deallocates the memory for a 1d array of floats.
  * Was allocated using f_alloc1d.
  *
  * @param p pointer to aray of floats.
  */
void f_free_1d(float *p)
{
  free_1d(p);
}

/** Allocate and clear a 2d array. Note that the storage for the values
  * is one large piece of memory so that it can be read or written
  * easily. The i,j value is accessed by saying array[j][i], because
  * array[j] is a pointer to the jth row of the matrix.
  * n1 is the number of columns, n2 is the number of rows.
  *
  * @param n1 size of first dimension of the array.
  * @param n2 size of second dimension of the array.
  * @return a 2d array n1*n2 float values.
  */
float **f_alloc_2d(int n1, int n2)
{
  return alloc_2d(n1,n2,sizeof(float));
}

/** Deallocates the memory for a 2d array of floats.
  * Was allocated using f_alloc2d.
  *
  * @param pp pointer to aray of floats.
  */
void f_free_2d(float **pp)
{
  free_2d(pp);
}


/** Allocate a 3d array. Note that the storage for the values
  * is one large piece of memory so that it can be read or written
  * easily. The i,j,k value is accessed by saying array[k][j][i].
  * array[k] is a pointer to the kth plane of the array. array[k][j]
  * is a pointer to the jth row in the kth plane of the array.
  *
  * @param n1 size of first dimension of the array.
  * @param n2 size of second dimension of the array.
  * @param n3 size of third dimension of the array.
  * @return a 3d array n1*n2*n3 float values.
  */
float ***f_alloc_3d(int n1, int n2, int n3)
{
  return (float ***)alloc_3d(n1, n2, n3,sizeof(float));
}

/** Deallocates the memory for a 3d array of floats.
  * Was allocated using f_alloc3d.
  *
  * @param ppp pointer to aray of floats.
  */
void f_free_3d(float ***ppp)
{
  free_3d(ppp);
}

/** Allocate a 4d array. Note that the storage for the values
  * is one large piece of memory so that it can be read or written
  * easily. The i,j,k,l value is accessed by saying array[l][k][j][i].
  *
  * @param n1 size of first dimension of the array.
  * @param n2 size of second dimension of the array.
  * @param n3 size of third dimension of the array.
  * @param n4 size of third dimension of the array.
  * @return a 4d array n1*n2*n3*n4 float values.
  */
float ****f_alloc_4d(int n1, int n2, int n3, int n4)
{
  return alloc_4d(n1, n2, n3, n4, sizeof(float));
}

/** Deallocates the memory for a 4d array of floats.
  * Was allocated using f_alloc4d.
  *
  * @param pppp pointer to aray of floats.
  */
void f_free_4d(float ****pppp)
{
  free_4d(pppp);
}

/** Allocate and clear a 1d array of long values.
  *
  * @param n1 size of array.
  * @return an array n1 of long values.
  */
long *l_alloc_1d(int n1)
{
  return (long *)alloc_1d(n1,sizeof(long));
}

/** Deallocates the memory for a 1d array of longs.
  * Was allocated using l_alloc1d.
  *
  * @param p pointer to aray of longs.
  */
void l_free_1d(long *p)
{
  free_1d(p);
}


/** Allocate and clear a 2d array. Note that the storage for the values
  * is one large piece of memory so that it can be read or written
  * easily. The i,j value is accessed by saying array[j][i], because
  * array[j] is a pointer to the jth row of the matrix.
  * n1 is the number of columns, n2 is the number of rows.
  *
  * @param n1 size of first dimension of the array.
  * @param n2 size of second dimension of the array.
  * @return a 2d array n1*n2 long values.
  */
long **l_alloc_2d(int n1, int n2)
{
  return (long **)alloc_2d(n1, n2,sizeof(long));
}


/** Deallocates the memory for a 2d array of longs.
  * Was allocated using l_alloc2d.
  *
  * @param pp pointer to aray of longs.
  */
void l_free_2d(long **pp)
{
  free_2d(pp);
}


/** Allocate a 3d array. Note that the storage for the values
  * is one large piece of memory so that it can be read or written
  * easily. The i,j,k value is accessed by saying array[k][j][i].
  * array[k] is a pointer to the kth plane of the array. array[k][j]
  * is a pointer to the jth row in the kth plane of the array.
  *
  * @param n1 size of first dimension of the array.
  * @param n2 size of second dimension of the array.
  * @param n3 size of third dimension of the array.
  * @return a 3d array n1*n2*n3 long values.
  */
long ***l_alloc_3d(int n1, int n2, int n3)
{
  return (long ***)alloc_3d(n1, n2, n3, sizeof(long));
}


/** Deallocates the memory for a 3d array of longs.
  * Was allocated using l_alloc3d.
  *
  * @param ppp pointer to aray of longs.
  */
void l_free_3d(long ***ppp)
{
  free_3d(ppp);
}


/** Allocate and clear a 1d array of integer values.
  *
  * @param n1 size of array.
  * @return an array n1 of integer values.
  */
int *i_alloc_1d(int n1)
{
  return (int *)alloc_1d(n1,sizeof(int));
}


/** Deallocates the memory for a 1d array of integers.
  * Was allocated using i_alloc1d.
  *
  * @param p pointer to aray of integers.
  */
/* free a 1d array */
void i_free_1d(int *p)
{
  free_1d(p);
}


/** Allocate and clear a 2d array. Note that the storage for the values
  * is one large piece of memory so that it can be read or written
  * easily. The i,j value is accessed by saying array[j][i], because
  * array[j] is a pointer to the jth row of the matrix.
  * n1 is the number of columns, n2 is the number of rows.
  *
  * @param n1 size of first dimension of the array.
  * @param n2 size of second dimension of the array.
  * @return a 2d array n1*n2 integer values.
  */
int **i_alloc_2d(int n1, int n2)
{
  return (int **)alloc_2d(n1, n2, sizeof(int));
}


/** Deallocates the memory for a 2d array of integers.
  * Was allocated using i_alloc2d.
  *
  * @param pp pointer to aray of integers.
  */
void i_free_2d(int **pp)
{
  free_2d(pp);
}


/** Allocate a 3d array. Note that the storage for the values
  * is one large piece of memory so that it can be read or written
  * easily. The i,j,k value is accessed by saying array[k][j][i].
  * array[k] is a pointer to the kth plane of the array. array[k][j]
  * is a pointer to the jth row in the kth plane of the array.
  *
  * @param n1 size of first dimension of the array.
  * @param n2 size of second dimension of the array.
  * @param n3 size of third dimension of the array.
  * @return a 3d array n1*n2*n3 integer values.
  */
int ***i_alloc_3d(int n1, int n2, int n3)
{
  return (int ***)alloc_3d(n1, n2, n3,sizeof(int));
}


/** Deallocates the memory for a 3d array of integers.
  * Was allocated using i_alloc3d.
  *
  * @param ppp pointer to aray of integers.
  */
void i_free_3d(int ***ppp)
{
  free_3d(ppp);
}


/** Allocate a 4d array. Note that the storage for the values
  * is one large piece of memory so that it can be read or written
  * easily. The i,j,k,l value is accessed by saying array[l][k][j][i].
  *
  * @param n1 size of first dimension of the array.
  * @param n2 size of second dimension of the array.
  * @param n3 size of third dimension of the array.
  * @param n4 size of third dimension of the array.
  * @return a 4d array n1*n2*n3*n4 integer values.
  */
int ****i_alloc_4d(int n1, int n2, int n3, int n4)
{
  return alloc_4d(n1, n2, n3, n4, sizeof(int));
}


/** Deallocates the memory for a 4d array of integers.
  * Was allocated using i_alloc4d.
  *
  * @param pppp pointer to aray of integers.
  */
void i_free_4d(int ****pppp)
{
  free_4d(pppp);
}


/** Allocate and clear a 1d array of short values.
  *
  * @param n1 size of array.
  * @return an array n1 of short values.
  */
short *s_alloc_1d(int n1)
{
  return (short *)alloc_1d(n1,sizeof(short));
}


/** Deallocates the memory for a 1d array of shorts.
  * Was allocated using s_alloc1d.
  *
  * @param p pointer to aray of shorts.
  */
void s_free_1d(short int *p)
{
  free_1d(p);
}


/** Allocate and clear a 2d array. Note that the storage for the values
  * is one large piece of memory so that it can be read or written
  * easily. The i,j value is accessed by saying array[j][i], because
  * array[j] is a pointer to the jth row of the matrix.
  * n1 is the number of columns, n2 is the number of rows.
  *
  * @param n1 size of first dimension of the array.
  * @param n2 size of second dimension of the array.
  * @return a 2d array n1*n2 short values.
  */
short **s_alloc_2d(int n1, int n2)
{
  return (short **)alloc_2d(n1, n2, sizeof(short));
}

/** Deallocates the memory for a 2d array of shorts.
  * Was allocated using s_alloc2d.
  *
  * @param pp pointer to aray of shorts.
  */
void s_free_2d(short int **pp)
{
  free_2d(pp);
}


/** Allocate a 3d array. Note that the storage for the values
  * is one large piece of memory so that it can be read or written
  * easily. The i,j,k value is accessed by saying array[k][j][i].
  * array[k] is a pointer to the kth plane of the array. array[k][j]
  * is a pointer to the jth row in the kth plane of the array.
  *
  * @param n1 size of first dimension of the array.
  * @param n2 size of second dimension of the array.
  * @param n3 size of third dimension of the array.
  * @return a 3d array n1*n2*n3 short values.
  */
short ***s_alloc_3d(int n1, int n2, int n3)
{
  return (short ***)alloc_3d(n1, n2, n3,sizeof(short));
}

/** Deallocates the memory for a 3d array of shorts.
  * Was allocated using s_alloc3d.
  *
  * @param ppp pointer to aray of shorts.
  */
void s_free_3d(short int ***ppp)
{
  free_3d(ppp);
}


/** Allocate a 4d array. Note that the storage for the values
  * is one large piece of memory so that it can be read or written
  * easily. The i,j,k,l value is accessed by saying array[l][k][j][i].
  *
  * @param n1 size of first dimension of the array.
  * @param n2 size of second dimension of the array.
  * @param n3 size of third dimension of the array.
  * @param n4 size of third dimension of the array.
  * @return a 4d array n1*n2*n3*n4 short values.
  */
short ****s_alloc_4d(int n1, int n2, int n3, int n4)
{
  return alloc_4d(n1, n2, n3, n4, sizeof(short));
}

/** Deallocates the memory for a 4d array of shorts.
  * Was allocated using s_alloc4d.
  *
  * @param pppp pointer to aray of shorts.
  */
void s_free_4d(short int ****pppp)
{
  free_4d(pppp);
}


/** Allocate and clear a 1d array of char values.
  *
  * @param n1 size of array.
  *
  * @return an array n1 of char values.
  */
char *c_alloc_1d(int n1)
{
  return (char *)alloc_1d(n1,sizeof(char));
}

/** Deallocates the memory for a 1d array of chars.
  * Was allocated using c_alloc1d.
  *
  * @param p pointer to aray of chars.
  */
void c_free_1d(char *p)
{
  free_1d(p);
}


/** Allocate and clear a 1d array of unsigned char values.
  *
  * @param n1 size of array.
  *
  * @return an array n1 of char values.
  */
unsigned char *uc_alloc_1d(int n1)
{
  return (unsigned char *)alloc_1d(n1,sizeof(unsigned char));
}


/** Deallocates the memory for a 1d array of unsigned chars.
  * Was allocated using uc_alloc1d.
  *
  * @param p pointer to aray of chars.
  */
void uc_free_1d(char *p)
{
  free_1d(p);
}


/** Allocate and clear a 2d array. Note that the storage for the values
  * is one large piece of memory so that it can be read or written
  * easily. The i,j value is accessed by saying array[j][i], because
  * array[j] is a pointer to the jth row of the matrix.
  * n1 is the number of columns, n2 is the number of rows.
  *
  * @param n1 size of first dimension of the array.
  * @param n2 size of second dimension of the array.
  * @return a 2d array n1*n2 char values.
  */
char **c_alloc_2d(int n1, int n2)
{
  return (char **)alloc_2d(n1, n2, sizeof(char));
}


/** Deallocates the memory for a 2d array of chars.
  * Was allocated using c_alloc2d.
  *
  * @param pp pointer to aray of chars.
  */
void c_free_2d(char **pp)
{
  free_2d(pp);
}


/** Allocate a 3d array. Note that the storage for the values
  * is one large piece of memory so that it can be read or written
  * easily. The i,j,k value is accessed by saying array[k][j][i].
  * array[k] is a pointer to the kth plane of the array. array[k][j]
  * is a pointer to the jth row in the kth plane of the array.
  *
  * @param n1 size of first dimension of the array.
  * @param n2 size of second dimension of the array.
  * @param n3 size of third dimension of the array.
  * @return a 3d array n1*n2*n3 char values.
  */
char ***c_alloc_3d(int n1, int n2, int n3)
{
  return (char ***)alloc_3d(n1, n2, n3, sizeof(char));
}


/** Deallocates the memory for a 3d array of chars.
  * Was allocated using c_alloc3d.
  *
  * @param ppp pointer to aray of chars.
  */
void c_free_3d(char ***ppp)
{
  free_3d(ppp);
}


/** Allocate a 4d array. Note that the storage for the values
  * is one large piece of memory so that it can be read or written
  * easily. The i,j,k,l value is accessed by saying array[l][k][j][i].
  *
  * @param n1 size of first dimension of the array.
  * @param n2 size of second dimension of the array.
  * @param n3 size of third dimension of the array.
  * @param n4 size of third dimension of the array.
  * @return a 4d array n1*n2*n3*n4 char values.
  */
char ****c_alloc_4d(int n1, int n2, int n3, int n4)
{
  return (char ****) alloc_4d(n1, n2, n3, n4, sizeof(char));
}


/** Deallocates the memory for a 4d array of chars.
  * Was allocated using c_alloc4d.
  *
  * @param pppp pointer to aray of chars.
  */
void c_free_4d(char ****pppp)
{
  free_4d(pppp);
}


 /*PS*/
/** Allocate and clear a 1d array of pointers.
  *
  * @param n1 size of array.
  * @return a 1d array of pointers.
  */
void **p_alloc_1d(int n1)
{
  return alloc_1d(n1, sizeof(void *));
}


/** Deallocates the memory for a 1d array of pointers.
  * Was allocated using p_alloc1d.
  *
  * @param p pointer to aray of pointers.
  */
void p_free_1d(void **p)
{
  free_1d(p);
}


/** Allocate and clear a 2d array of pointers. Note that the storage for the 
  * values is one large piece of memory so that it can be read or written
  * easily. The i,j value is accessed by saying array[j][i], because
  * array[j] is a pointer to the jth row of the matrix.
  * n1 is the number of columns, n2 is the number of rows.
  *
  * @param n1 size of first dimension of the array.
  * @param n2 size of second dimension of the array.
  * @return a 2d array n1*n2 pointers.
  */
void ***p_alloc_2d(int n1, int n2)
{
  return alloc_2d(n1, n2, sizeof(void *));
}


/** Deallocates the memory for a 2d array of pointers.
  * Was allocated using p_alloc2d.
  *
  * @param pp pointer to aray of pointers.
  */
void p_free_2d(void ***pp)
{
  free_2d(pp);
}


/** Allocate a 3d array of pointers. Note that the storage for the values
  * is one large piece of memory so that it can be read or written
  * easily. The i,j,k value is accessed by saying array[k][j][i].
  * array[k] is a pointer to the kth plane of the array. array[k][j]
  * is a pointer to the jth row in the kth plane of the array.
  *
  * @param n1 size of first dimension of the array.
  * @param n2 size of second dimension of the array.
  * @param n3 size of third dimension of the array.
  * @return a 3d n1*n2*n3 array of pointers.
  */
void ****p_alloc_3d(int n1, int n2, int n3)
{
  return alloc_3d(n1, n2, n3, sizeof(void *));
}


/** Deallocates the memory for a 3d array of pointers.
  * Was allocated using p_alloc3d.
  *
  * @param ppp pointer to aray of pointers.
  */
void p_free_3d(void ****ppp)
{
  free_3d(ppp);
}


/** Allocate a 4d array of pointers. Note that the storage for the values
  * is one large piece of memory so that it can be read or written
  * easily. The i,j,k,l value is accessed by saying array[l][k][j][i].
  *
  * @param n1 size of first dimension of the array.
  * @param n2 size of second dimension of the array.
  * @param n3 size of third dimension of the array.
  * @param n4 size of third dimension of the array.
  * @return a 4d array n1*n2*n3*n4 pointers.
  */
void *****p_alloc_4d(int n1, int n2, int n3, int n4)
{
  return alloc_4d(n1, n2, n3, n4, sizeof(void *));
}


/** Deallocates the memory for a 4d array of pointers.
  * Was allocated using p_alloc4d.
  *
  * @param pppp pointer to aray of pointers.
  */
void p_free_4d(void *****pppp)
{
  free_4d(pppp);
}




/* Allocates n1xn2xn3xn4 4D array of something. Note that it will be accessed as 
 * [n4][n3][n2][n1].
 * @param n1 Size of the first dimension
 * @param n2 Size of the second dimension
 * @param n3 Size of the third dimension
 * @param n4 Size of the fourth dimension
 * @param unitsize Size of one array element in bytes
 * @return 4D array
 */
void* alloc_4d(int n1, int n2, int n3, int n4, size_t unitsize)
{
  size_t size,ind;
  const size_t MAX_ELEMS = -1; // Works because size_t will always be unsigned
    char* p;
    char** pp;
    char*** ppp;
    char**** pppp;
    int i;

    assert(n1 > 0);
    assert(n2 > 0);
    assert(n3 > 0);
    assert(n4 > 0);
    /*
     * We should probably roll this out everywhere but I want to keep tabs on every case
     * FR: May 2015
     // assert((double) n1 * (double) n2 * (double) n3 * (double) n4 < (double) UINT_MAX);
    */
    assert((double) n1 * (double) n2 * (double) n3 * (double) n4 < (double) MAX_ELEMS);

    size = ((size_t)n1) * ((size_t)n2) * ((size_t)n3) * ((size_t)n4);
    if ((p = calloc(size, unitsize)) == NULL)
        quit("alloc_4d() out of memory:  %s %d\n", strerror(errno),(size * unitsize));

    assert((double) n2 * (double) n3 * (double) n4 * (double) sizeof(void*) < (double) UINT_MAX);

    size = ((size_t)n2) * ((size_t)n3) * ((size_t)n4);
    if ((pp = malloc(size * sizeof(void*))) == NULL)
        quit("alloc_4d() assigning rows failed:  %s\n", strerror(errno));
    for (i = 0; i < size; i++) {
      ind = ((size_t)i) * ((size_t)n1) * ((size_t)unitsize);
      pp[i] = &p[ind];
    }

    assert((double) n3  * (double) n4 * (double) sizeof(void*) < (double) UINT_MAX);

    size = ((size_t)n3) * ((size_t)n4);
    if ((ppp = malloc(size * sizeof(void*))) == NULL)
        quit("alloc_4d() assigning planes failed:  %s %d\n", strerror(errno),(size * sizeof(void*)));
    for (i = 0; i < size; i++)
        ppp[i] = &pp[i * n2];

    assert((double) n4 * (double) sizeof(void*) < (double) UINT_MAX);

    size = ((size_t)n4) * sizeof(void*);
    if ((pppp = malloc(size)) == NULL)
        quit("alloc_4d() assigning cube failed:  %s %d\n", strerror(errno),size);
    for (i = 0; i < n4; i++)
        pppp[i] = &ppp[i * n3];
        
    return pppp;
}


/* Destroys 4D array.
 * @param pppp 4D array
 */
void free_4d(void* pppp)
{
    void*** ppp;
    void** pp;
    void* p;

    assert(pppp != NULL);
    p = ((void****) pppp)[0][0][0];
    pp = ((void****) pppp)[0][0];
    ppp = ((void****) pppp)[0];
    free(pppp);
    assert(ppp != NULL);
    free(ppp);
    assert(pp != NULL);
    free(pp);
    assert(p != NULL);
    free(p);
    pppp = NULL;
}


/* Allocates n1xn2xn3 3D array of something. Note that it will be accessed as 
 * [n3][n2][n1].
 * @param n1 Size of the first dimension
 * @param n2 Size of the second dimension
 * @param n3 Size of the third dimension
 * @param unitsize Size of one array element in bytes
 * @return 3D array
 */
void* alloc_3d(int n1, int n2, int n3, size_t unitsize)
{
    size_t size;
    char* p;
    char** pp;
    char*** ppp;
    int i;
    
    assert(n1 > 0);
    assert(n2 > 0);
    assert(n3 > 0);
    assert((double) n1 * (double) n2 * (double) n3 < (double) UINT_MAX);

   size = ((size_t)n1) * ((size_t)n2) * ((size_t)n3);
    if ((p = calloc(size, unitsize)) == NULL)
        quit("alloc_3d() memory failed:  %s %d\n", strerror(errno),(size * unitsize));

    assert((double) n2 * (double) n3 * (double) sizeof(void*) < (double) UINT_MAX);

    size = ((size_t)n2) * ((size_t)n3);
    if ((pp = malloc(size * sizeof(void*))) == NULL)
        quit("alloc_3d()assigning rows failed:  %s\n", strerror(errno));
    for (i = 0; i < size; i++)
        pp[i] = &p[i * n1 * unitsize];

    assert((double) n3 * (double) sizeof(void*) < (double) UINT_MAX);

    size = ((size_t)n3) * sizeof(void*);
    if ((ppp = malloc(size)) == NULL)
        quit("alloc_3d() assigning planes failed:  %s\n", strerror(errno));
    for (i = 0; i < n3; i++)
        ppp[i] = &pp[i * n2];

    return ppp;
}

/* Destroys 3D array.
 * @param ppp 3D array
 */
void free_3d(void* ppp)
{
    void** pp;
    void* p;
    assert(ppp != NULL);
    p = ((void***) ppp)[0][0];
    pp = ((void***) ppp)[0];
    free(ppp);
    assert(pp != NULL);
    free(pp);
    assert(p != NULL);
    free(p);
    ppp = NULL;
}

/* Allocates n1xn2 matrix of something. Note that it will be accessed as 
 * [n2][n1].
 * @param n1 Number of columns
 * @param n2 Number of rows
 * @param unitsize Size of one matrix element in bytes
 * @return Matrix
 */
void* alloc_2d(int n1, int n2, size_t unitsize)
{
    size_t size;
    char* p;
    char** pp;
    int i;

    assert(n1 > 0);
    assert(n2 > 0);
    assert((double) n1 * (double) n2 <= (double) UINT_MAX);

    size = ((size_t)n1) * ((size_t)n2);
    if ((p = calloc(size, unitsize)) == NULL)
        quit("alloc_2d(): %s\n", strerror(errno));

    assert((double) n2 * (double) sizeof(void*) <= (double) UINT_MAX);

    size =  ((size_t)n2) * sizeof(void*);
    if ((pp = malloc(size)) == NULL)
        quit("alloc_2d(): %s\n", strerror(errno));
    for (i = 0; i < n2; i++)
        pp[i] = &p[i * n1 * unitsize];

    return pp;
}

/* Destroys n1xn2 matrix.
 * @param pp Matrix
 */
void free_2d(void* pp)
{
    void* p;

    assert(pp != NULL);
    p = ((void**) pp)[0];
    free(pp);
    assert(p != NULL);
    free(p);
    pp = NULL;
}


/* Allocates n1 array of something. Note that it will be accessed as 
 * [n1].
 * @param n1 Number of columns
 * @param unitsize Size of one array element in bytes
 * @return Array
 */
void* alloc_1d(int n1, size_t unitsize)
{
    char* p;
		size_t size =  (size_t)n1;
    assert(n1 > 0);
    assert((double) n1 <= (double) UINT_MAX);

    if ((p = calloc(size, unitsize)) == NULL)
        quit("alloc_1d(): %s unit %d\n", strerror(errno),unitsize);

    return p;
}


/* Destroys n1 array.
 * @param p Array
 */
void free_1d(void* p)
{
    assert(p != NULL);
    free(p);
    p = NULL;
}
