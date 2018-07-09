/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/da/da_utils.cpp
 *  
 *  Description:
 *  Data assimilation utility functions
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: da_utils.cpp 5850 2018-06-29 05:18:59Z riz008 $
 *
 */

#include "da_utils.hpp"

/**
 * Constructor
 */
da_maps::da_maps(double **gx, double **gy, double *gz, int nx, int ny, int nz) 
{
  int i,j,k;

  /* Cache the grid centre sizes */
  f_ni = nx;
  f_nj = ny;
  f_nk = nz;

  /* Create the xytoij tree */
  f_tree = grid_xytoij_init(gx, gy, nx, ny);
  
  /* Copy the vertical array - remember layers is z_grid */
  f_layers = d_alloc_1d(nz+1);
  memcpy(f_layers, gz, (nz+1)*sizeof(double));

  /* Allocate the map arrays */
  f_map_2D = i_alloc_2d(f_ni, f_nj);
  f_map_3D = i_alloc_3d(f_ni, f_nj, f_nk);

  /* Initialise the maps */
  for (k=0; k<f_nk; k++)
    for (j=0; j<f_nj; j++)
      for (i=0; i<f_ni; i++) {
	f_map_3D[k][j][i] = DA_INVALID_VAL;
	if (k==0)
	  f_map_2D[j][i] = DA_INVALID_VAL;
      }

  /* Reserve more than enough space */
  f_map_2D_s2i.reserve(f_ni * f_nj);
  f_map_2D_s2j.reserve(f_ni * f_nj);

  f_map_3D_s2i.reserve(f_ni * f_nj * f_nk);
  f_map_3D_s2j.reserve(f_ni * f_nj * f_nk);
  f_map_3D_s2k.reserve(f_ni * f_nj * f_nk);
}


/**
 * Destructor
 */
da_maps::~da_maps(void)
{
  if (f_tree != NULL)
    tree_destroy(f_tree);

  if (f_layers != NULL)
    d_free_1d(f_layers);

  if (f_map_3D != NULL)
    i_free_3d(f_map_3D);

  if (f_map_2D != NULL)
    i_free_2d(f_map_2D);
}


/******************/
/* PUBLIC METHODS */
/******************/

/** 
 * Sets a local sparse map value for 2D
 */
void da_maps::fill_map(int ls, int i, int j) {

  if (i < f_ni && j < f_nj) {
    /* Store in the 2D array */
    f_map_2D[j][i] = ls;
    /* 
     * Cache backwards lookup 
     * Note: Its assumed that this is embeddeded within a counting loop
     */
    f_map_2D_s2i.push_back(i);
    f_map_2D_s2j.push_back(j);
    
  } else
    quit("DA-lib: 2D Dimension map out of bounds error");
}


/**
 * Sets a local sparse map value for 3D
 */
void da_maps::fill_map(int ls, int i, int j, int k) {

  if (i < f_ni && j < f_nj && k < f_nk) {
    f_map_3D[k][j][i] = ls;
    /* Cache reverse maps */
    f_map_3D_s2i.push_back(i);
    f_map_3D_s2j.push_back(j);
    f_map_3D_s2k.push_back(k);
  }
  else
    quit("DA-lib: 3D Dimension map out of bounds error");
}


/**
 * Get the sparse index for 2D from coordinate space
 */
int da_maps::getS(double x, double y)
{
  int i,j;
  
  if (!getIJ(x, y, &i, &j))
    return(DA_INVALID_VAL);

  return(getS(i, j));
}


/**
 * Get the sparse index for 3D from coordinate space
 */
int da_maps::getS(double x, double y, double z)
{
  int i,j;
  
  if (!getIJ(x, y, &i, &j))
    return(DA_INVALID_VAL);

  return(getS(i, j, getK(z)));
}


/**
 * Get the sparse index, given 2D index plus depth
 */
int da_maps::getS(int i, int j, double z)
{
  return(getS(i, j, getK(z)));
}

/**
 * Get ij from xy
 */
int da_maps::getIJ(double x, double y, int *i, int *j)
{
  return(grid_xytoij(f_tree, x, y, i, j));
}

/**
 * Get ijk from sparse location
 */
int da_maps::getIJK(int s, int *i, int *j, int *k)
{
  if (s >= f_map_3D_s2i.size())
    return(0);

  *i = f_map_3D_s2i[s];
  *j = f_map_3D_s2j[s];
  *k = f_map_3D_s2k[s];

  return(1);
}


/**
 *  Get xy from a sparse location
 *
 *  Note: This should not be called for 3D data
 */
int da_maps::getXY(int s, double *x, double *y)
{
  int i = f_map_2D_s2i[s];
  int j = f_map_2D_s2j[s];
  
  return(getXY(i, j, x, y));
}

/**
 *  Get xy from ij
 */
int da_maps::getXY(int i, int j, double *x, double *y)
{
  /* This should always return true */
  if (grid_ijtoxy(f_tree, i, j, x, y)) {
    /* All good */
    return(1);
  }
  
  return(0);
}

/**
 * Get the surface coordinate
 */
int da_maps::getCS(double x, double y)
{
  int k = f_nk-1;
  int i, j;
  /* Start from the top */
  if (grid_xytoij(f_tree, x, y, &i, &j)) {
    for ( ; k >= 0; k--) {
      int s = f_map_3D[k][j][i];
      if (s > DA_INVALID_VAL)
	return(s);
    }
  }
  /* Fail condition */
  return(DA_INVALID_VAL);
}

/*******************/
/* PRIVATE METHODS */
/*******************/
/**
 * Get the sparse coordinate for 3D from index space
 */
int da_maps::getS(int i, int j, int k)
{
  int s;

  if (i < f_ni && j < f_nj && k < f_nk && (s = f_map_3D[k][j][i]) > DA_INVALID_VAL)
    return(s);
  else
    return(DA_INVALID_VAL);
}


/** 
 * Get the sparse coordinate for 2D from index space
 */
int da_maps::getS(int i, int j)
{
  int s;

  if (i < f_ni && j < f_nj && (s = f_map_2D[j][i]) > DA_INVALID_VAL)
    return(s);
  else
    return(DA_INVALID_VAL);
}


/** Search for the appropriate k level
 *
 */
int da_maps::getK(double z)
{
  int k;

  /* See if this is too deep */
  if (z <= f_layers[0])
    return(-1);

  for (k=0; k<f_nk+1; k++) {
    if (f_layers[k] > z)
      break;
  }
  return(k-1);
}

/* end da_utils.cpp */
