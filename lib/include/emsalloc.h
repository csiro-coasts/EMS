/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/emsalloc.h
 *
 *  \brief Header file for memory allocation library
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: emsalloc.h 5834 2018-06-27 00:55:37Z riz008 $
 */


#ifdef __cplusplus
extern "C" {
#endif

#if HAVE_ALLOCA_H
# include <alloca.h>
#elif defined _MSC_VER || defined __BORLANDC__ || defined __MINGW32__
# include <malloc.h>
#elif defined __GNUC__
# define alloca __builtin_alloca
#elif defined _AIX
# define alloca __alloca
#else
# include <stddef.h>
# ifdef  __cplusplus
extern "C"
# endif
void *alloca (size_t);
#endif


  double *d_alloc_1d(int n1);
  void d_free_1d(double *p);
  double **d_alloc_2d(int n1, int n2);
  void d_free_2d(double **pp);
  double ***d_alloc_3d(int n1, int n2, int n3);
  void d_free_3d(double ***ppp);
  double ****d_alloc_4d(int n1, int n2, int n3, int n4);
  void d_free_4d(double ****ppp);
  float *f_alloc_1d(int n1);
  void f_free_1d(float *p);
  float **f_alloc_2d(int n1, int n2);
  void f_free_2d(float **pp);
  float ***f_alloc_3d(int n1, int n2, int n3);
  void f_free_3d(float ***ppp);
  float ****f_alloc_4d(int n1, int n2, int n3, int n4);
  void f_free_4d(float ****ppp);
  long *l_alloc_1d(int n1);
  void l_free_1d(long *p);
  long **l_alloc_2d(int n1, int n2);
  void l_free_2d(long **pp);
  long ***l_alloc_3d(int n1, int n2, int n3);
  void l_free_3d(long ***ppp);
  int *i_alloc_1d(int n1);
  void i_free_1d(int *p);
  int **i_alloc_2d(int n1, int n2);
  void i_free_2d(int **pp);
  int ***i_alloc_3d(int n1, int n2, int n3);
  void i_free_3d(int ***ppp);
  int ****i_alloc_4d(int n1, int n2, int n3, int n4);
  void i_free_4d(int ****ppp);
  short *s_alloc_1d(int n1);
  void s_free_1d(short *p);
  short **s_alloc_2d(int n1, int n2);
  void s_free_2d(short **pp);
  short ***s_alloc_3d(int n1, int n2, int n3);
  void s_free_3d(short ***ppp);
  short ****s_alloc_4d(int n1, int n2, int n3, int n4);
  void s_free_4d(short ****ppp);
  char *c_alloc_1d(int n1);
  void c_free_1d(char *p);
  unsigned char *uc_alloc_1d(int n1);
  void uc_free_1d(char *p);
  char **c_alloc_2d(int n1, int n2);
  void c_free_2d(char **pp);
  char ***c_alloc_3d(int n1, int n2, int n3);
  void c_free_3d(char ***ppp);
  char ****c_alloc_4d(int n1, int n2, int n3, int n4);
  void c_free_4d(char ****ppp);
  void **p_alloc_1d(int n1);
  void p_free_1d(void **p);
  void ***p_alloc_2d(int n1, int n2);
  void p_free_2d(void ***pp);
  void ****p_alloc_3d(int n1, int n2, int n3);
  void p_free_3d(void ****ppp);
  void *****p_alloc_4d(int n1, int n2, int n3, int n4);
  void p_free_4d(void *****ppp);

/* Added UR 3/11/2004 */
  void *alloc_1d(int n1, size_t unitsize);
  void *alloc_2d(int n1, int n2, size_t unitsize);
  void *alloc_3d(int n1, int n2, int n3, size_t unitsize);
  void *alloc_4d(int n1, int n2, int n3, int n4, size_t unitsize);
  void free_1d(void* p);
  void free_2d(void* p);
  void free_3d(void* p);
  void free_4d(void* p);
/* UR */
#ifdef  __cplusplus
}
#endif
