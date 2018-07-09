/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/da/da_utils.hpp
 *  
 *  Description:
 *  Header file for data assimilation utility functions
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: da_utils.hpp 5850 2018-06-29 05:18:59Z riz008 $
 *
 */

extern "C" {
#include <stdlib.h>
#include <math.h>
#include "ems.h"
}

#include <list>
#include <vector>
using namespace std;

#define DA_INVALID_VAL (-99)


/** Class to facilitate cartesian to sparse conversions
 *
 *  This encompasses both the 2D and 3D maps
 *
 */
class da_maps
{
public:
  /** Constructor
   */
  da_maps(double **gx, double **gy, double *gz, int nx, int ny, int nz);
  ~da_maps(void);

  /** Public methods */

  /** These are the main exports. Both from index and coord spaces
   */
  int getS(double x, double y);
  int getS(double x, double y, double z);
  int getS(int i, int j, double z);
  int getCS(double x, double y);
  /* Returns 1 on success, 0 failure */
  int getIJ(double x, double y, int *i, int *j);
  int getIJK(int s, int *i, int *j, int *k);
  /* Returns 1 on success, 0 failure */
  int getXY(int s, double *x, double *y);
  int getXY(int i, int j, double *x, double *y);
  
  /** Fill valid data
   */
  void fill_map(int ls, int i, int j);
  void fill_map(int ls, int i, int j, int k);

  /**
   * Size
   */
  inline int get_s2d_len(void) { 
    return(f_map_2D_s2i.size());
  }

  /**
   * Returns true is the given X,Y point is in the domain
   */
  inline bool inDomain(double X, double Y) {
    if (getS(X,Y) > DA_INVALID_VAL)
      return(true);
    return(false);
  }

private:
  /** Tree pointer from the EMS lib
   */
  xytoij_tree_t *f_tree;

  /**
   * 3D cartesian sizes
   */
  int f_ni;
  int f_nj;
  int f_nk;

  /** Vertical geometry
   */
  double *f_layers;

  /**
   * 2D & 3D Cartesian to local sparse maps
   */
  int  **f_map_2D;
  int ***f_map_3D;

  /**
   * 2D sparse to index maps. We could cache the x,y as well here
   */
  vector<int> f_map_2D_s2i;
  vector<int> f_map_2D_s2j;

  /**
   * 3D sparse to index maps. We could cache the x,y as well here
   */
  vector<int> f_map_3D_s2i;
  vector<int> f_map_3D_s2j;
  vector<int> f_map_3D_s2k;

  /** Methods */
  int getS(int i, int j);
  int getS(int i, int j, int k);
  int getK(double z);
};


/* end da_utils.hpp */
