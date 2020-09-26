/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/master/pp.c
 *  
 *  Description: Sparse preprocessor
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: pp.c 5841 2018-06-28 06:51:55Z riz008 $
 *
 */


/****************************************************************************

This module contains preprocessor sample code which will generate a sparse
(compressed) array, associated mappings and vector subscripts given a
bathymetry array, layer structure and open boundary information as input. 
The user may have to modify nomenclature to render compatable with 
existing structures etc.

The routine build_sparse_grid() generates the sparse array, mappings and 
vectors. Sparse grid information is stored in the data structure 'sgrid', 
defined below. Open boundary information for boundary n is stored in 
sgrid->onc[n].

Arguments to buld_sparse_grid() are:

int nce1;        : Grid size in x direction. The actual array size of
the model should be nce1+1 to account for velocity faces on the eastern
limit of the domain.

int nce2;        : Grid size in y direction. The actual array size of
the model should be nce2+1 to account for velocity faces on the northern
limit of the domain.

int nz;          : Size of the grid in the vertical direction.

double **bathy;  : Cell centered bathymetry array of size bathy[nce2][nce1].
All bathymety values should be negative. The value 0 is assumed to be
mean sea level. A cell is identified as a land cell if its value > 10. 

double *layers;  : Layer structure of size layers[nz];
Layers should be an array of size nz, where layers[nz] = maximum depth and
layers[0] =0.

int nobc;        : Number of open boundaries.
int npts[nobc];  : Number of surface cells in each open boundary.
int iloc[nobc][npts]; : List of i locations for each open boundary.
int jloc[nobc][npts]; : List of j locations for each open boundary.
int type[nobc];       : Boundary orientation for each open boundary. If
the u1 component of velocity is normal to the open boundary, the boundary
is denoted as U1BDRY and the type is given the value 1. If the u2 
component of velocity is normal to the open boundary, the boundary is 
denoted as U2BDRY and the type is given the value 2.

e.g. the following domain has three open boundaries, two cross-shelf
U1BDRY and one along-shore U2BDRY.

 int nce1 = 4;
 int nce2 = 3;
 int nz = 11;
 double bathy[3][4] = {{-60, -70, -80, -90},   <- wet cells
                       {-55, -65, -75, -85},   <- wet cells
                       { 99,  99,  99,  99}};  <- land cells    
 double layers[11] = {-100, -90, -80, -70, -60, -50, -40, -30, -20, -10, 0};
 int nobc = 3;
 int npts[3] = {2, 2, 4};
 int iloc[3][4] = {{0, 0},         <- boundary 0
                   {3, 3},         <- boundary 1
                   {0, 1, 2, 3}};  <- boundary 2
 int jloc[3][4] = {{1, 2},         <- boundary 0
                   {1, 2},         <- boundary 1
                   {2, 2, 2, 2}};  <- boundary 2
 int type[3] = {1, 1, 2};


  Written by M. Herzfeld, 2006.


*****************************************************************************/

#include <stdlib.h>
#include <string.h>
#include "hd.h"
#include "tracer.h"

#define NOTWET (SOLID|OUTSIDE)
#define W_SE 1
#define E_SE 2
#define S_SE 4
#define N_SE 8

#define NE_D 1
#define SE_D 2
#define NW_D 4
#define SW_D 8

#define NE_IC 16
#define SE_IC 32
#define NW_IC 64
#define SW_IC 128

#define WET 256


/*-------------------------------------------------------------------*/
/* Sparse mapping from Cartesian space to sparse space. The sparse   */
/* map is filled in the x direction first, followed by the y and z   */
/* directions. The ghost points are located after the wet points for */
/* each layer. This allows the first layer to act as a map from 2D   */
/* cartesian space to the sparse grid.                               */
/* The grid geometry is defined below, where 'o' corresponds to a    */
/* wet cell and 'x' to a land cell (i.e. a ghost point).             */
/*                                                                   */
/*                     o o o o                                       */
/* Straight edges :    x x x x                                       */
/* (south)                                                           */
/*                                                                   */
/* Outside corner :    o o o o                                       */
/* (south-west)        x x x o                                       */
/*                     x x x o                                       */
/*                     x x x o                                       */
/*                                                                   */
/* Inside corner :     x o o o                                       */
/* (south-west)        x o o o                                       */
/*                     x o o o                                       */
/*                     x x x x                                       */
/*                                                                   */
/* Diagonal :          x x o o                                       */
/*                     x x o o                                       */
/*                     o o x x                                       */
/*                     o o x x                                       */
/*                                                                   */
/* All direction are defined with refererence to the compass points, */
/* north being towards the page top.                                 */
/*                                                                   */
/*                               north                               */
/*                                 |                                 */
/*                             NW  |  NE                             */
/*                                 |                                 */
/*                     west ------ o ------ east                     */
/*                                 |                                 */
/*                             SW  |  SE                             */
/*                                 |                                 */
/*                               south                               */
/*                                                                   */
/*-------------------------------------------------------------------*/
void build_sparse_grid(parameters_t *params,
		       geometry_t *sgrid,
		       int nce1,       /* Grid size, x direction     */
		       int nce2,       /* Grid size, y direction     */
		       int nz,         /* Grid size, z direction     */
		       double **bathy, /* Bathymetry array           */
		       double *layers, /* Layer structure            */
		       int nobc,       /* Number of open boundaries  */
		       int *npts,      /* No. cells in each boundary */
		       int **iloc,     /* Boundary x location list   */
		       int **jloc,     /* Boundary y location list   */
		       int *type       /* Boundary orientation       */
		       )
{
  int c = 1;         /* Locations of wet points in the sparse array  */
  int cc;            /* Locations of ghost points in sparse array    */
  int lc;            /* Loc of boundary points in the boundary array */
  int i, j, k;       /* x,y,z counters                               */
  int ii, jj;        /* x,y counters                                 */
  int n;             /* General counters                             */
  int c1, c2, c3;    /* Counters                                     */
  int b1, b2, b3;    /* Counters                                     */
  int xp1;           /* xp1 = map[k][j][i+1]                         */
  int xm1;           /* xm1 = map[k][j][i-1]                         */
  int yp1;           /* yp1 = map[k][j+1][i]                         */
  int ym1;           /* ym1 = map[k][j-1][i]                         */
  int *end_lc;       /* End loaction of the wet points in each layer */
  int *num_2d;       /* Number of ghost points in each layer         */
  int num_dc = 0;    /* Number of diagonal connection ghost cells    */
  int num_se = 0;    /* # straight edge/outside corner ghost cells   */
  int num_ic = 0;    /* Number of inside corner ghost cells          */
  int num_mc = 0;    /* Number of multiple ghost cells               */
  int num_wc = 0;    /* Number of 3D wet grid cells                  */
  int num_gc = 0;    /* Number of ghost grid cells                   */
  int num_sc = 0;    /* Number of sediment grid cells                */
  int num_scg = 0;   /* Number of sediment ghost grid cells          */
  int num_gc2D = 0;  /* Number of 2D ghost grid cells                */
  int num_obc = 0;   /* Number of OBC ghost cells                    */
  int num_obc2D = 0; /* Number of 2D OBC ghost cells                 */
  short **kbot;      /* k index of the bottom                        */
  long scen;         /* Flag for ghost cell type                     */
  int gchck = 0;     /* Checks locations of ghost points             */
  int nfe1;          /* Grid dimension in the e1 direction           */
  int nfe2;          /* Grid dimension in the e2 direction           */
  unsigned long ***flg; /* Flag for Cartesian grid                   */
  unsigned long ***sw;  /* South-west corner corner array            */
  unsigned long ***se;  /* South-east corner corner array            */
  unsigned long ***nw;  /* North-west corner corner array            */
  unsigned long ***ne;  /* North-east corner corner array            */
  int *ocodex;       /* x location of OBCs                           */
  int *ocodey;       /* y location of OBCs                           */
  double bmax;       /* Maximum depth                                */


  int sigma = 0;     /* Set to 1 for sigma model                     */

  /* Read the grid dimensions                                        */
  nfe1 = nce1 + 1;
  nfe2 = nce2 + 1;

  /* Set OUTSIDE cells                                               */
  for (i = 0; i < nce1; i++)
    for (j = 0; j < nce2; j++)
      if (isnan(bathy[j][i])) bathy[j][i] = NOTVALID;

  /* Reset flags for sigma calculations if required                  */
  if (sigma) {
    /* Find the maximum depth                                        */
    bmax = 1e10;
    for (i = 0; i < nce1; i++)
      for (j = 0; j < nce2; j++)
	if (bathy[j][i] < bmax)
	  bmax = bathy[j][i];
    for (n = 0; n < nz; n++)
      layers[n] *= fabs(bmax);
  }

  /* Allocate memory                                                 */
  flg = (unsigned long ***)l_alloc_3d(nce1 + 1, nce2 + 1, nz);

  /* Make the flags array */
  make_flags(params, flg, bathy, layers, nce1, nce2, nz);
  u1_flags(flg, nce1, nce2, nz);
  u2_flags(flg, nce1, nce2, nz);
  /* For a sigma model                                               */
  if (sigma)
    sigma_flags(nfe1, nfe2, nz - 1, flg);

  /* Initialise                                                      */
  end_lc = i_alloc_1d(nz + 1);
  num_2d = i_alloc_1d(nz + 1);
  sgrid->nce1 = nce1;
  sgrid->nce2 = nce2;
  sgrid->nfe1 = nfe1;
  sgrid->nfe2 = nfe2;
  sgrid->nz = nz;
  sgrid->nobc = nobc;
  for (k = 0; k <= nz; k++)
    end_lc[k] = 0;

  /*-----------------------------------------------------------------*/
  /* Set up the array of the bottom coordinate, including the edges  */
  /* defined by i=nce1 and j=nce2.  */
  kbot = s_alloc_2d(nfe1, nfe2);
  /* Set the bottom coordinate in the interior                       */
  for (j = 0; j < nce2; j++) {
    for (i = 0; i < nce1; i++) {
      if (flg[nz - 1][j][i] & NOTWET) {
        kbot[j][i] = -1;
      } else {
        /* Check bottom is inside model                              */
        if (bathy[j][i] < layers[0]) {
            hd_quit("build_sparse_grid() : Bottom lower than lowest model level\n");
	}
        /* Check bottom and top below MAXGRIDZ */
        if (bathy[j][i] >= MAXGRIDZ) {
          hd_quit("build_sparse_grid() : Bottom above MAXGRIDZ\n");
	}
        /* Loop vertically to find bottom and top */
        kbot[j][i] = -1;
        if (sigma)
          kbot[j][i] = 0;
        else {
          for (k = 0; k < nz; k++) {
            if (bathy[j][i] >= layers[k] && bathy[j][i] < layers[k + 1])
              kbot[j][i] = k;
          }
        }
        if (kbot[j][i] < 0) {
	  hd_quit("build_sparse_grid() : bottom not found at (%d %d) : %5.1f m\n",
             i, j, bathy[j][i]);
	}
      }
    }
  }

  /* Set the bottom coordinate on the edges                          */
  k = nz - 1;
  i = nce1;
  for (j = 0; j < nfe2; j++)
    if (!(flg[k][j][i] & NOTWET))
      kbot[j][i] = kbot[j][i - 1];
  j = nce2;
  for (i = 0; i < nfe1; i++)
    if (!(flg[k][j][i] & NOTWET))
      kbot[j][i] = kbot[j - 1][i];

  /* Set the flags for the open boundaries                           */
  for (n = 0; n < nobc; n++) {
    for (cc = 0; cc < npts[n]; cc++) {
      i = iloc[n][cc];
      j = jloc[n][cc];

      kbot[j][i] = set_kbot(i, j, flg[nz - 1], kbot, type[n], nce1, nce2);
      for (k = kbot[j][i]; k < nz; k++) {
        if (k == -1) {
          hd_quit("Land cell found in boundary %d list at %d %d\n", n, i, j);
	}
        c2 = b_flags(nce1, nce2, nz, flg, i, j, k, &ii, &jj, 1, type[n]);
      }

    }
  }

  /* Reset the flag on U1BDRY|R_EDGE and U2BDRY|F_EDGE to wet.       */ 
  for (n = 0; n < nobc; n++) {
    for (cc = 0; cc < npts[n]; cc++) {
      i = iloc[n][cc];
      j = jloc[n][cc];
      for (k = kbot[j][i]; k < nz; k++) {
	if (type[n] == U1BDRY || type[n] == U2BDRY) {
	  if (flg[k][j][i] & (U1BDRY | R_EDGE) || 
	      flg[k][j][i] & (U2BDRY | F_EDGE)) {
	    if (flg[k][j][i] & OUTSIDE) {
	      flg[k][j][i] &= ~OUTSIDE;
	    }
	  }
	}
      }
    }
  }
  corner_flags(flg, nce1, nce2, nz);

  /*-----------------------------------------------------------------*/
  /* Generate the map from Cartesian space to sparse space,          */
  /* sgrid->map. First map the cells in the entire grid that may      */
  /* potentially be wet. The sparse grid at this stage consists of   */
  /* wet cells stored in consecutive order, i.e. a 'wet' sparse map  */
  /* is created. The number of wet cells in each layer is stored in  */
  /* end_lc[]. This is used later to insert the ghost cells for each */
  /* layer into locations following the wet cells in each layer. The */
  /* 'wet' map only exists so that the ghost cells can be identified */
  /* at this stage and is reset when the ghost points are inserted.  */
  /* The total number of wet cells is stored in num_wc here.         */
  sgrid->ewetS = 0;
  sgrid->map = (unsigned long ***)l_alloc_3d(nfe1, nfe2, nz);
  for (k = nz - 1; k >= 0; k--) {
    for (j = 0; j < nfe2; j++) {
      for (i = 0; i < nfe1; i++) {
        sgrid->map[k][j][i] = 0;
        if (!(flg[k][j][i] & NOTWET)) {
          sgrid->map[k][j][i] = c;
          c++;
          num_wc++;
          if (k == nz - 1)
            sgrid->ewetS++;
          end_lc[k]++;
        }
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* If no wet points are found then issue a warning and quit.       */
  if (num_wc == 0) {
    hd_quit("No wet points located in the domain.\n");
  }

  /*-----------------------------------------------------------------*/
  /* Find the locations of the ghost cells. These are the first land */
  /* (dry) cells adjacent to a wet cell in the x, y or z direction.  */
  /* The ghost cells in the z direction constitute the sediment      */
  /* layer. Note : in certain instances (e.g. outside corners) there */
  /* may exist more than one wet cell associated with the same ghost */
  /* cell -> multiple ghost cells). First count the total number of  */
  /* ghost cells in the grid and store the number of ghost cells in  */
  /* each layer in num_2d[].                                         */
  c1 = 0;
  for (k = nz - 1; k >= 0; k--) {
    num_2d[k] = 0;
    for (j = 0; j < nfe2; j++) {
      for (i = 0; i < nfe1; i++) {

        /*-----------------------------------------------------------*/
        /* Count the number of diagonal & inside corner ghost cells  */
        /* Note : it is possible for a ghost cell to constitute more */
        /* than one inside corner or diagonal, hence each corner or  */
        /* diagonal must be checked individually.                    */
        scen = di_ic(sgrid, i, j, k, &xp1, &xm1, &yp1, &ym1);
        if (scen & NW_D) {
          num_dc++;
          num_2d[k]++;
        }
        if (scen & NE_D) {
          num_dc++;
          num_2d[k]++;
        }
        if (scen & SW_D) {
          num_dc++;
          num_2d[k]++;
        }
        if (scen & SE_D) {
          num_dc++;
          num_2d[k]++;
        }
        if (scen & NW_IC) {
          num_ic++;
          num_2d[k]++;
        }
        if (scen & NE_IC) {
          num_ic++;
          num_2d[k]++;
        }
        if (scen & SW_IC) {
          num_ic++;
          num_2d[k]++;
        }
        if (scen & SE_IC) {
          num_ic++;
          num_2d[k]++;
        }

        /*-----------------------------------------------------------*/
        /* Count the number of straight edge and outside corner      */
        /* ghost cells.                                              */
        scen = se_oc(sgrid, i, j, k, &xp1, &xm1, &yp1, &ym1, 1);
        if (scen & W_SE)
          num_2d[k]++;
        if (scen & E_SE)
          num_2d[k]++;
        if (scen & S_SE)
          num_2d[k]++;
        if (scen & N_SE)
          num_2d[k]++;
      }
    }
    num_gc += num_2d[k];
  }
  num_se = num_gc - (num_dc + num_ic);
  sgrid->nbpt = sgrid->nbpte1 = sgrid->nbpte2 = num_gc;
  sgrid->nbptS = sgrid->nbpte1S = sgrid->nbpte2S = num_gc2D = num_2d[nz - 1];

  /*-----------------------------------------------------------------*/
  /* Calculate the number of sediment cells and add to the total     */
  /* number of ghost cells. Note: the number of sediment cells is    */
  /* equivalent to the number of potentially wet cells in the        */
  /* surface layer.                                                  */
  /* First add the number of ghost cells to sedimet cells (this is   */
  /* the same as the number of ghost cells for the surface layer).   */
  num_scg = 0;
  if (!(params->compatible & V1598) && (params->trasc & LAGRANGE))
    num_scg = sgrid->nbptS;    
  num_sc += num_scg;
  for (j = 0; j < nfe2; j++)
    for (i = 0; i < nfe1; i++)
      if (!(flg[nz - 1][j][i] & NOTWET))
        num_sc += 1;
  num_gc += num_sc;

  /*-----------------------------------------------------------------*/
  /* Get the number of OBC ghost cells. There are laux cells for     */
  /* open boundary point so that the advection scheme may be used on */
  /* the boundary.                                                   */

  ocodex = i_alloc_1d(nobc);
  ocodey = i_alloc_1d(nobc);
  for (n = 0; n < nobc; n++) {
    if (params->compatible & V1957) continue;
    for (cc = 0; cc < npts[n]; cc++) {
      i = iloc[n][cc];
      j = jloc[n][cc];
      if (type[n] & U1BDRY) {
	if (i == nce1 - 1 || 
	    (i < nce1 && (bathy[j][i] < layers[0] || bathy[j][i] > params->etamax))) {
	  ocodex[n] = R_EDGE;
	  for (k = kbot[j][i]; k < nz; k++) {
	    num_2d[k] += laux - 1;
	    num_gc += laux - 1;
	    num_obc += laux - 1;
	  }
	  num_gc2D += laux - 1;
	  num_obc2D += laux - 1;
	} else {
	  ocodex[n] = L_EDGE;
	  for (k = kbot[j][i]; k < nz; k++) {
	    num_2d[k] += laux;
	    num_gc += laux;
	    num_obc += laux;
	  }
	  num_gc2D += laux;
	  num_obc2D += laux;
	}
      }
      if (type[n] & U2BDRY) {
	if (j == nce2 - 1 || 
	    (j < nce2 && (bathy[j][i] < layers[0] || bathy[j][i] > params->etamax))) {
	  ocodey[n] = F_EDGE;
	  for (k = kbot[j][i]; k < nz; k++) {
	    num_2d[k] += laux - 1;
	    num_gc += laux - 1;
	    num_obc += laux - 1;
	  }
	  num_gc2D += laux - 1;
	  num_obc2D += laux - 1;
	} else {
	  ocodey[n] = B_EDGE;
	  for (k = kbot[j][i]; k < nz; k++) {
	    num_2d[k] += laux;
	    num_gc += laux;
	    num_obc += laux;
	  }
	  num_gc2D += laux;
	  num_obc2D += laux;
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Re-organize the sparse grid so that the ghost cells in each     */
  /* layer are inserted at the end of the wet cells for that layer.  */
  /* This results in the wet cells corresponding to new sparse       */
  /* locations in the sparse grid, so the sgrid->map must be reset.   */
  /* Only maps from Cartesian to sparse space for wet cells are set; */
  /* ghost cells in the sparse grid are set to zero at this stage to */
  /* be assigned values later. The sparse locations of the sediment  */
  /* layer lies at the end of the sparse array.                      */

  /* Find the end location of wet cells in each layer. This is also  */
  /* used later to find the start location of ghost cells in each    */
  /* layer to assign the ghost to wet and wet to ghost spatial maps. */
  for (k = nz - 2; k >= 0; k--) {
    c = end_lc[k + 1] + num_2d[k + 1];
    end_lc[k] += c;
  }

  /* Reinitialise the sparse map                                     */
  for (k = nz - 1; k >= 0; k--)
    for (j = 0; j < nfe2; j++)
      for (i = 0; i < nfe1; i++)
        sgrid->map[k][j][i] = 0;

  /* Reset the Cartestian map to point to new sparse locations       */
  c = 1;
  for (k = nz - 1; k >= 0; k--) {
    for (j = 0; j < nfe2; j++) {
      for (i = 0; i < nfe1; i++) {
        if (!(flg[k][j][i] & NOTWET)) {
          sgrid->map[k][j][i] = c;
          if ((gchck == c))
            hd_warn("%d = wet cell at (%d %d %d)\n", gchck, i, j, k);
          c = c + 1;
        }
      }
    }
    c = end_lc[k] + num_2d[k] + 1;
  }

  /*-----------------------------------------------------------------*/
  /* Save the total number of cells in the sparse grid.              */
  /* Note : these cell locations are for the re-organized array (2D  */
  /* wet + ghost cells followed by 3D wet + ghost) and hence do not  */
  /* reflect the number of cells in the sparse array; e.g. ewet is   */
  /* the number of wet cells in the grid but the number of wet + 2D  */
  /* ghost cells.                                                    */
  sgrid->sgnum = sgrid->enon = num_wc + num_gc;
  sgrid->sgnumS = sgrid->enonS = sgrid->ewetS + num_gc2D;
  sgrid->ewet = num_wc + num_gc2D;
  sgrid->snon = sgrid->ewet + 1;
  sgrid->snonS = sgrid->ewetS + 1;
  /* The sparse grid holds an empty location at 0 (unused) hence the */
  /* size of memory for arrays is sgrid->sgnum+1.                     */
  sgrid->sgsiz = sgrid->sgnum + 1;
  sgrid->sgsizS = sgrid->sgnumS + 1;

  /*-----------------------------------------------------------------*/
  /* Allocate memory for the spatial maps and initialize             */
  alloc_sgrid(sgrid);

  /* All the maps are made self-mapping at this stage. Those maps    */
  /* which are not re-assigned below remain self pointing.           */
  for (c = 1; c <= sgrid->sgnum; c++) {
    sgrid->xp1[c] = c;
    sgrid->yp1[c] = c;
    sgrid->zp1[c] = c;
    sgrid->xm1[c] = c;
    sgrid->ym1[c] = c;
    sgrid->zm1[c] = c;
  }

  /* Allocate memory for the corner arrays (used to define multiple  */
  /* ghost cells) and initialise. Since the possibility exists for   */
  /* wet cells to have multiple ghost points, by splitting the ghost */
  /* cells into those that occupy west, east, north and south edges, */
  /* then a unique ghost cell can be identified and placed in the    */
  /* corner arrays, bl, br, tl and tr. Therefore, using the four     */
  /* corner arrays a unique map for ghost points from Cartesian      */
  /* space to sparse space is achieved.                              */
  sw = (unsigned long ***)l_alloc_3d(nfe1, nfe2, nz);
  se = (unsigned long ***)l_alloc_3d(nfe1, nfe2, nz);
  nw = (unsigned long ***)l_alloc_3d(nfe1, nfe2, nz);
  ne = (unsigned long ***)l_alloc_3d(nfe1, nfe2, nz);
  for (k = nz - 1; k >= 0; k--)
    for (j = 0; j < nfe2; j++)
      for (i = 0; i < nfe1; i++) {
        sw[k][j][i] = 0;
        se[k][j][i] = 0;
        nw[k][j][i] = 0;
        ne[k][j][i] = 0;
      }

  /*-----------------------------------------------------------------*/
  /* Generate the spatial maps for the sparse grid. This includes    */
  /* wet cells that map to wet cells (unique map), wet cells that    */
  /* map to ghost cells (unique map), ghost cells that map to wet    */
  /* cells (multiple maps), and ghost cells that map to ghost cells  */
  /* (unique map).                                                   */
  /* First assign the maps of wet cells that map to other wet cells. */
  /* These are defined by cells with non-zero neighbours in the map  */
  /* sgrid->map. The bottom layer maps to the sediment layer and      */
  /* vice versa. The sediment layer maps downward to itself (i.e.    */
  /* self pointing).                                                 */
  for (k = nz - 1; k >= 0; k--) {
    for (j = 0; j < nfe2; j++) {
      for (i = 0; i < nfe1; i++) {

        if (k < kbot[j][i])
          continue;
        scen = se_oc(sgrid, i, j, k, &xp1, &xm1, &yp1, &ym1, 0);

        if (scen != 0) {

          /* c = the sparse location of wet cells at (i,j,k)         */
          c = sgrid->map[k][j][i];

          if (xp1 != 0)
            sgrid->xp1[c] = sgrid->map[k][j][i + 1];
          if (xm1 != 0)
            sgrid->xm1[c] = sgrid->map[k][j][i - 1];
          if (yp1 != 0)
            sgrid->yp1[c] = sgrid->map[k][j + 1][i];
          if (ym1 != 0)
            sgrid->ym1[c] = sgrid->map[k][j - 1][i];
          if (k > kbot[j][i] && sgrid->map[k - 1][j][i] != 0) {
            sgrid->zm1[c] = sgrid->map[k - 1][j][i];
          }
          if (k < nz - 1 && sgrid->map[k + 1][j][i] != 0) {
            sgrid->zp1[c] = sgrid->map[k + 1][j][i];
          }
        }
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Allocate memory for the boundary vectors. The vector bpt        */
  /* contains the sparse locations of all ghost cells (land boundary */
  /* cells) laterally adjacent to wet cells. The vector bin contains */
  /* the sparse locations of the first wet cell adjacent to the land */
  /* boundary.                                                       */
  sgrid->bpt = i_alloc_1d(sgrid->nbpt + 1);
  sgrid->bin = i_alloc_1d(sgrid->nbpt + 1);
  sgrid->dbin = i_alloc_1d(sgrid->nbpt + 1);
  sgrid->wgst = i_alloc_1d(sgrid->sgsiz);

  for (c = 1; c <= sgrid->nbpt; c++) {
    sgrid->bpt[c] = 0;
    sgrid->bin[c] = 0;
  }

  /*-----------------------------------------------------------------*/
  /* Assign the spatial maps of ghost cells that map to wet cells    */
  /* and wet cells which map to ghost cells. Also assign the sparse  */
  /* locations to bpt and bin. In each layer the ghost cell counter, */
  /* cc, is set to end_lc[] and incremented. Wet cells corresponding */
  /* to diagonals, straight edges, inside/outside corners are        */
  /* identified in sgrid->map on the basis of their neighbour's       */
  /* values (zero = ghost; nonzero = wet) and then assigned to cc    */
  /* for the map.                                                    */
  lc = 1;
  for (k = nz - 1; k >= 0; k--) {

    /* cc = the sparse location of ghost cells for layer k           */
    cc = end_lc[k] + 1;

    for (j = 0; j < nfe2; j++) {
      for (i = 0; i < nfe1; i++) {

        /*-----------------------------------------------------------*/
        /* Map the straight edges and outside corners.               */
        /* Note that ghost points currently are set to zero in       */
        /* sgrid->map, so to get the ghost-wet map a land cell is     */
        /* first checked to see if its neighbours are wet. If so     */
        /* then set the relevant wet-ghost spatial map equal to the  */
        /* ghost cell counter, cc, and the ghost-wet spatial map to  */
        /* the sparse coordinate of the wet neighbour cell.          */
        scen = se_oc(sgrid, i, j, k, &xp1, &xm1, &yp1, &ym1, 1);

        /* Ghost cells on west edges and west outside corners        */
        if (scen & W_SE) {
          /* Save the ghost cell location in the corner array        */
          sw[k][j][i] = se[k][j][i] = nw[k][j][i] = ne[k][j][i] = cc;

          /* Check for maps on east straight edges                   */
          sgrid->xm1[xp1] = cc;  /* Wet-ghost spatial map             */
          sgrid->xp1[cc] = xp1;  /* Ghost-wet spatial map             */

          /* Check for ghost-wet maps on east outside corners, ie.   */
          /* if the east edge comprises an outside corner there      */
          /* will also exist ghost-wet maps in the north/south       */
          /* direction depending on the orientation of the corner.   */
          /* The wet-ghost map for these corners will be defined     */
          /* with the relevant straight edge map below.              */
          if (yp1 != 0)
            sgrid->yp1[cc] = yp1;
          if (ym1 != 0)
            sgrid->ym1[cc] = ym1;

          /* Add the ghost cell to the boundary vectors              */
          sgrid->bpt[lc] = cc;
          sgrid->bin[lc] = xp1;
          sgrid->dbin[lc] = 2;
          if ((gchck == cc))
            hd_warn("%d = west edge/outside corner at (%d %d %d)\n",
                 gchck, i, j, k);
          lc++;
          cc++;
        }

        /* Ghost cells on east edges and east outside corners        */
        if (scen & E_SE) {
          /* Save the ghost cell location in the corner array        */
          sw[k][j][i] = se[k][j][i] = nw[k][j][i] = ne[k][j][i] = cc;

          /* Check for maps on west straight edges                   */
          sgrid->xp1[xm1] = cc;  /* Wet-ghost spatial map             */
          sgrid->xm1[cc] = xm1;  /* Ghost-wet spatial map             */
          /* Check for ghost-wet maps on west outside corners        */
          if (yp1 != 0)
            sgrid->yp1[cc] = yp1;
          if (ym1 != 0)
            sgrid->ym1[cc] = ym1;

          /* Add the ghost cell to the boundary vectors              */
          sgrid->bpt[lc] = cc;
          sgrid->bin[lc] = xm1;
          sgrid->dbin[lc] = 6;
          if ((gchck == cc))
            hd_warn("%d = east edge/outside corner at (%d %d %d)\n",
                 gchck, i, j, k);
          lc++;
          cc++;
        }

        /* Ghost cells on south edges and south outside corners      */
        if (scen & S_SE) {
          /* Save the ghost cell location in the corner array        */
          sw[k][j][i] = se[k][j][i] = nw[k][j][i] = ne[k][j][i] = cc;

          /* Check for maps on south straight edges                  */
          sgrid->ym1[yp1] = cc;  /* Wet-ghost spatial map             */
          sgrid->yp1[cc] = yp1;  /* Ghost-wet spatial map             */
          /* Check for ghost-wet maps on south outside corners       */
          if (xp1 != 0)
            sgrid->xp1[cc] = xp1;
          if (xm1 != 0)
            sgrid->xm1[cc] = xm1;

          /* Add the ghost cell to the boundary vectors              */
          sgrid->bpt[lc] = cc;
          sgrid->bin[lc] = yp1;
          sgrid->dbin[lc] = 0;
          if ((gchck == cc))
            hd_warn("%d = south edge/outside corner at (%d %d %d)\n", gchck,
                 i, j, k);
          lc++;
          cc++;
        }

        /* Ghost cells on north edges and north outside corners      */
        if (scen & N_SE) {
          /* Save the ghost cell location in the corner array        */
          sw[k][j][i] = se[k][j][i] = nw[k][j][i] = ne[k][j][i] = cc;

          /* Check for maps on north straight edges                  */
          sgrid->yp1[ym1] = cc;  /* Wet-ghost spatial map             */
          sgrid->ym1[cc] = ym1;  /* Ghost-wet spatial map             */
          /* Check for ghost-wet maps on north outside corners       */
          if (xp1 != 0)
            sgrid->xp1[cc] = xp1;
          if (xm1 != 0)
            sgrid->xm1[cc] = xm1;

          /* Add the ghost cell to the boundary vectors              */
          sgrid->bpt[lc] = cc;
          sgrid->bin[lc] = ym1;
          sgrid->dbin[lc] = 4;
          if ((gchck == cc))
            hd_warn("%d = north edge/outside corner at (%d %d %d)\n", gchck,
                 i, j, k);
          lc++;
          cc++;
        }

        /*-----------------------------------------------------------*/
        /* Map diagonals and inside corners. This must be performed  */
        /* after the straight edge part so that the corner arrays    */
        /* are written with unique diagonal/inside corner sparse     */
        /* locations. Note: do not use an else if() construct here   */
        /* so that any multiple inside corners/diagonals are set.    */
        /* Map the diagonal ghost cells. Here a wet cell is checked  */
        /* to see if its neighbours are not wet and if its diagonal  */
        /* is wet. If so the ghost counter, cc, is incrmented, bpt   */
        /* and the corner array are set to cc, bin is set to the wet */
        /* cell checked. Note: no direct map exist from the wet cell */
        /* checked to its diagonal ghost; these are achieved through */
        /* a wet-ghost and ghost-ghost map.                          */

        /* Get the ghost locations of the inside corners. These are  */
        /* located by checking if the land cell's neighbours are     */
        /* also land cells (non-wet) but the diagonal is wet. As     */
        /* with the diagonal ghost cells, a direct map doesn't exist */
        /* from land cell in the corner to its diagonal wet cell;    */
        /* these are achieved through the ghost-ghost (defined       */
        /* below) and ghost wet (defined above) maps. However, the   */
        /* boundary vectors and corner arrays need to be assigned    */
        /* here.                                                     */

        scen = di_ic(sgrid, i, j, k, &xp1, &xm1, &yp1, &ym1);
        if (scen & (SW_D | SW_IC)) {
          sw[k][j][i] = cc;
          sgrid->bpt[lc] = cc;
          sgrid->bin[lc] = sgrid->map[k][j + 1][i + 1];
          sgrid->dbin[lc] = 1;
          if (scen & SW_D && gchck == cc)
            hd_warn("%d = south-west diagonal at (%d %d %d)\n",
                 gchck, i, j, k);
          if (scen & SW_IC && gchck == cc)
            hd_warn("%d = south-west inside corner at (%d %d %d)\n",
                 gchck, i, j, k);
          lc++;
          cc++;
        }
        if (scen & (NW_D | NW_IC)) {
          nw[k][j][i] = cc;
          sgrid->bpt[lc] = cc;
          sgrid->bin[lc] = sgrid->map[k][j - 1][i + 1];
          sgrid->dbin[lc] = 3;
          if (scen & NW_D && gchck == cc)
            hd_warn("%d = north-west diagonal at (%d %d %d)\n",
                 gchck, i, j, k);
          if (scen & NW_IC && gchck == cc)
            hd_warn("%d = north-west inside corner at (%d %d %d)\n",
                 gchck, i, j, k);
          lc++;
          cc++;
        }
        if (scen & (SE_D | SE_IC)) {
          se[k][j][i] = cc;
          sgrid->bpt[lc] = cc;
          sgrid->bin[lc] = sgrid->map[k][j + 1][i - 1];
          sgrid->dbin[lc] = 7;
          if (scen & SE_D && gchck == cc)
            hd_warn("%d = south-east diagonal at (%d %d %d)\n",
                 gchck, i, j, k);
          if (scen & SE_IC && gchck == cc)
            hd_warn("%d = south-east inside corner at (%d %d %d)\n",
                 gchck, i, j, k);
          lc++;
          cc++;
        }
        if (scen & (NE_D | NE_IC)) {
          ne[k][j][i] = cc;
          sgrid->bpt[lc] = cc;
          sgrid->bin[lc] = sgrid->map[k][j - 1][i - 1];
          sgrid->dbin[lc] = 5;
          if (scen & NE_D && gchck == cc)
            hd_warn("%d = north-east diagonal at (%d %d %d)\n",
                 gchck, i, j, k);
          if (scen & NE_IC && gchck == cc)
            hd_warn("%d = north-east inside corner at (%d %d %d)\n",
                 gchck, i, j, k);
          lc++;
          cc++;
        }
      }
    }
    /* Save the last sparse location of ghost cells for layer k */
    end_lc[k] = cc;
  }

  /* Set the locations of the corner arrays for the case where the   */
  /* cell above the corner cell is wet to the wet cell. This allows  */
  /* inside corners that lie below wet cells to map upwards to the   */
  /* wet cell. If the cell below is a inside corner also, then the   */
  /* corner array will be defined and the downward map will be set   */
  /* via the ghost-ghost maps. If the cell below is not a ghost cell */
  /* then the inside corner remains self pointing. Note that         */
  /* although the inside corner maps upward to the wet cell, the wet */
  /* cell will map downward to the sediment cell.                    */
  for (k = nz - 2; k >= 0; k--) {
    for (j = 0; j < nfe2; j++) {
      for (i = 0; i < nfe1; i++) {
        if (!(flg[nz - 1][j][i] & NOTWET)) {
          if (sw[k][j][i] != 0 && sw[k + 1][j][i] == 0)
            sw[k + 1][j][i] = sgrid->map[k + 1][j][i];
          if (se[k][j][i] != 0 && se[k + 1][j][i] == 0)
            se[k + 1][j][i] = sgrid->map[k + 1][j][i];
          if (nw[k][j][i] != 0 && nw[k + 1][j][i] == 0)
            nw[k + 1][j][i] = sgrid->map[k + 1][j][i];
          if (ne[k][j][i] != 0 && ne[k + 1][j][i] == 0)
            ne[k + 1][j][i] = sgrid->map[k + 1][j][i];
        }
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Assign the spatial maps of ghost cells that map to ghost cells. */
  /* Ghost cells are still set to zero (and not unique) in the the   */
  /* Cartesian-sparse map, sgrid->map, but have been assigned unique  */
  /* values in the sparse grid and corner arrays in the above        */
  /* mapping procedures.                                             */

  /* Map the ghost-ghost straight edges, outside corners and maps    */
  /* from these locations to inside corners and diagonals.           */
  /* To obtain the ghost-ghost map for straight edge/outside corners */
  /* a wet cell is first checked to see if its neighbours are non-   */
  /* wet ghost cells using the maps sgrid->map and flg. If            */
  /* so, the map from this ghost cell to other neighbouring ghost    */
  /* cells (if any) is made by approaching the ghost cells via the   */
  /* wet cells. Maps from straight edge ghost points to inside       */
  /* corner ghost points are assigned also. The corner arrays are    */
  /* used to retrieve the ghost cell in the corner and again this    */
  /* cell is approached via the wet cell.                            */
  for (k = nz - 1; k >= 0; k--) {
    for (j = 0; j < nfe2; j++) {
      for (i = 0; i < nfe1; i++) {
        if (k < kbot[j][i])
          continue;
        if (!(flg[k][j][i] & NOTWET)) {

          /* c = the sparse location of wet cells at (i,j,k)         */
          c = sgrid->map[k][j][i];
          scen = se_oc(sgrid, i, j, k, &xp1, &xm1, &yp1, &ym1, 0);

          /* Ghost cells on west edges and west outside corners      */
          if (scen & W_SE) {
            /* Lateral maps for straight coast/outside corners       */
            /* Here (i,j,k) is wet, (i-1,j,k) is dry.  */
            if (yp1 != 0 && flg[k][j + 1][i - 1] & NOTWET)
              sgrid->yp1[sgrid->xm1[c]] = sgrid->xm1[yp1];
            if (ym1 != 0 && flg[k][j - 1][i - 1] & NOTWET)
              sgrid->ym1[sgrid->xm1[c]] = sgrid->xm1[ym1];

            /* Lateral maps from west edge ghost cells to inside     */
            /* corners.                                              */
            /* Here (i,j,k) is wet, (i-1,j,k) is dry and (i,j+1,k)   */
            /* & (i-1,j+1,k) are dry for north-west inside corners   */
            /* and (i,j-1,k) & (i-1,j-1,k) are dry for south-west    */
            /* inside corners.                                       */
            if (j < nce2 && flg[k][j + 1][i] & NOTWET &&
                flg[k][j + 1][i - 1] & NOTWET)
              sgrid->yp1[sgrid->xm1[c]] = nw[k][j + 1][i - 1];
            if (j > 0 && flg[k][j - 1][i] & NOTWET &&
                flg[k][j - 1][i - 1] & NOTWET)
              sgrid->ym1[sgrid->xm1[c]] = sw[k][j - 1][i - 1];

            /* Lateral maps from west edge ghost cells to diagonals. */
            /* Here (i,j,k) is wet, (i-1,j,k) is dry and (i,j+1,k)   */
            /* & (i-1,j+1,k) are dry/wet for north-west diagonals    */
            /* (i,j-1,k) & (i-1,j-1,k) are dry/wet for south-west    */
            /* and diagonals.                                        */
            if (j < nce2 && flg[k][j + 1][i] & NOTWET &&
                !(flg[k][j + 1][i - 1] & NOTWET))
              sgrid->yp1[sgrid->xm1[c]] = nw[k][j + 1][i - 1];
            if (j > 0 && flg[k][j - 1][i] & NOTWET &&
                !(flg[k][j - 1][i - 1] & NOTWET))
              sgrid->ym1[sgrid->xm1[c]] = sw[k][j - 1][i - 1];

            /* Ghost-ghost vertical maps                             */
            if (k > kbot[j][i]) {
              sgrid->zm1[sgrid->xm1[c]] = sgrid->xm1[sgrid->map[k - 1][j][i]];
            }
            if (k < nz - 1)
              sgrid->zp1[sgrid->xm1[c]] = sgrid->xm1[sgrid->map[k + 1][j][i]];
          }

          /* Ghost cells on east edges                               */
          if (scen & E_SE) {
            /* Lateral maps for straight coast/outside corners       */
            /* Here (i,j,k) is wet, (i+1,j,k) is dry.                */
            if (yp1 != 0 && flg[k][j + 1][i + 1] & NOTWET)
              sgrid->yp1[sgrid->xp1[c]] = sgrid->xp1[yp1];
            if (ym1 != 0 && flg[k][j - 1][i + 1] & NOTWET)
              sgrid->ym1[sgrid->xp1[c]] = sgrid->xp1[ym1];

            /* Lateral maps from east edge ghost cells to inside     */
            /* corners.                                              */
            /* Here (i,j,k) is wet, (i+1,j,k) is dry and (i,j+1,k)   */
            /* & (i+1,j+1,k) are dry for north-east inside corners   */
            /* and (i,j-1,k) & (i+1,j-1,k) are dry for south         */
            /* east inside corners.                                  */
            if (j < nce2 && flg[k][j + 1][i] & NOTWET &&
                flg[k][j + 1][i + 1] & NOTWET)
              sgrid->yp1[sgrid->xp1[c]] = ne[k][j + 1][i + 1];
            if (j > 0 && flg[k][j - 1][i] & NOTWET &&
                flg[k][j - 1][i + 1] & NOTWET)
              sgrid->ym1[sgrid->xp1[c]] = se[k][j - 1][i + 1];

            /* Lateral maps from east edge ghost cells to diagonals. */
            /* Here (i,j,k) is wet, (i+1,j,k) is dry and (i,j+1,k)   */
            /* & (i+1,j+1,k) are dry/wet for north-east diagonals    */
            /* and (i,j-1,k) & (i+1,j-1,k) are dry/wet for south-east*/
            /* diagonals.                                            */
            if (j < nce2 && flg[k][j + 1][i] & NOTWET &&
                !(flg[k][j + 1][i + 1] & NOTWET))
              sgrid->yp1[sgrid->xp1[c]] = ne[k][j + 1][i + 1];
            if (j > 0 && flg[k][j - 1][i] & NOTWET &&
                !(flg[k][j - 1][i + 1] & NOTWET))
              sgrid->ym1[sgrid->xp1[c]] = se[k][j - 1][i + 1];

            /* Ghost-ghost vertical maps                             */
            if (k > kbot[j][i])
              sgrid->zm1[sgrid->xp1[c]] = sgrid->xp1[sgrid->map[k - 1][j][i]];
            if (k < nz - 1)
              sgrid->zp1[sgrid->xp1[c]] = sgrid->xp1[sgrid->map[k + 1][j][i]];
          }

          /* Ghost cells on north edges                              */
          if (scen & N_SE) {

            /* Lateral maps for straight coast/outside corners       */
            /* Here (i,j,k) is wet, (i,j+1,k) is dry.                */
            if (xp1 != 0 && flg[k][j + 1][i + 1] & NOTWET)
              sgrid->xp1[sgrid->yp1[c]] = sgrid->yp1[xp1];
            if (xm1 != 0 && flg[k][j + 1][i - 1] & NOTWET)
              sgrid->xm1[sgrid->yp1[c]] = sgrid->yp1[xm1];

            /* Lateral maps from north edge ghost cells to inside    */
            /* corners.  */
            /* Here (i,j,k) is wet, (i,j+1,k) is dry and (i+1,j,k)   */
            /* & (i+1,j+1,k) are dry for north-east inside corners   */
            /* and (i-1,j,k) & (i-1,j+1,k) are dry for north-west    */
            /* inside corners.                                       */
            if (i < nce1 && flg[k][j][i + 1] & NOTWET &&
                flg[k][j + 1][i + 1] & NOTWET)
              sgrid->xp1[sgrid->yp1[c]] = ne[k][j + 1][i + 1];
            if (i > 0 && flg[k][j][i - 1] & NOTWET &&
                flg[k][j + 1][i - 1] & NOTWET)
              sgrid->xm1[sgrid->yp1[c]] = nw[k][j + 1][i - 1];

            /* Lateral maps from north edge ghost cells to diagonals.*/
            /* Here (i,j,k) is wet, (i,j+1,k) is dry and (i+1,j,k)   */
            /* & (i+1,j+1,k) are dry/wet for north-east diagonals    */
            /* and (i-1,j,k) & (i-1,j+1,k) are dry/wet for north-west*/
            /* diagonals.                                            */
            if (i < nce1 && flg[k][j][i + 1] & NOTWET &&
                !(flg[k][j + 1][i + 1] & NOTWET))
              sgrid->xp1[sgrid->yp1[c]] = ne[k][j + 1][i + 1];
            if (i > 0 && flg[k][j][i - 1] & NOTWET &&
                !(flg[k][j + 1][i - 1] & NOTWET))
              sgrid->xm1[sgrid->yp1[c]] = nw[k][j + 1][i - 1];

            /* Ghost-ghost vertical maps                             */
            if (k > kbot[j][i])
              sgrid->zm1[sgrid->yp1[c]] = sgrid->yp1[sgrid->map[k - 1][j][i]];
            if (k < nz - 1)
              sgrid->zp1[sgrid->yp1[c]] = sgrid->yp1[sgrid->map[k + 1][j][i]];
          }

          /* Ghost cells on south edges                              */
          if (scen & S_SE) {

            /* Lateral maps for straight coast/outside corners       */
            /* Here (i,j,k) is wet, (i,j-1,k) is dry.                */
            if (xp1 != 0 && !(flg[k][j][i + 1] & NOTWET) &&
                flg[k][j - 1][i + 1] & NOTWET)
              sgrid->xp1[sgrid->ym1[c]] = sgrid->ym1[xp1];
            if (xm1 != 0 && !(flg[k][j][i - 1] & NOTWET) &&
                flg[k][j - 1][i - 1] & NOTWET)
              sgrid->xm1[sgrid->ym1[c]] = sgrid->ym1[xm1];

            /* Lateral maps from south edge ghost cells to inside    */
            /* corners.                                              */
            /* Here (i,j,k) is wet, (i,j-1,k) is dry and (i+1,j,k)   */
            /* & (i+1,j-1,k) are dry for south-east inside           */
            /* corners and (i-1,j,k) & (i-1,j-1,k) are dry for       */
            /* south-west inside corners.                            */
            if (i < nce1 && flg[k][j][i + 1] & NOTWET &&
                flg[k][j - 1][i + 1] & NOTWET)
              sgrid->xp1[sgrid->ym1[c]] = se[k][j - 1][i + 1];
            if (i > 0 && flg[k][j][i - 1] & NOTWET &&
                flg[k][j - 1][i - 1] & NOTWET)
              sgrid->xm1[sgrid->ym1[c]] = sw[k][j - 1][i - 1];

            /* Lateral maps from south edge ghost cells to           */
            /* diagonals.                                            */
            /* Here (i,j,k) is wet, (i,j-1,k) is dry and (i+1,j,k)   */
            /* & (i+1,j-1,k) are dry/wet for south-east diagonals    */
            /* corners and (i-1,j,k) & (i-1,j-1,k) are dry & wet for */
            /* south-west diagonals.                                 */
            if (i < nce1 && flg[k][j][i + 1] & NOTWET &&
                !(flg[k][j - 1][i + 1] & NOTWET))
              sgrid->xp1[sgrid->ym1[c]] = se[k][j - 1][i + 1];
            if (i > 0 && flg[k][j][i - 1] & NOTWET &&
                !(flg[k][j - 1][i - 1] & NOTWET))
              sgrid->xm1[sgrid->ym1[c]] = sw[k][j - 1][i - 1];

            /* Ghost-ghost vertical maps */
            if (k > kbot[j][i]) {
              sgrid->zm1[sgrid->ym1[c]] = sgrid->ym1[sgrid->map[k - 1][j][i]];
            }

            if (k < nz - 1) {
              sgrid->zp1[sgrid->ym1[c]] = sgrid->ym1[sgrid->map[k + 1][j][i]];
            }
          }
        }
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Assign the ghost-ghost maps from inside corners and diagonals   */
  /* to straight edges.                                              */
  /* To obtain the map for diagonal ghost cells to adjacent ghost    */
  /* cells a wet cell is first checked to see if its neighbours are  */
  /* dry using flg. Then if the diagonal is wet the diagonal         */
  /* ghost cell corresponding to the initial wet cell is retrieved   */
  /* via the corner array, and the maps from this ghost cell to      */
  /* other neighbouring ghost cells is made by approaching the ghost */
  /* cells via the wet cells.                                        */
  /* To obtain the map for inside corner ghost cells to adjacent     */
  /* ghost cells a wet cell's neighbours are first checked to see if */
  /* they are dry using flg. Then if the diagonal is also dry        */
  /* the maps are assigned in the same manner as the diagonal ghost  */
  /* cells.                                                          */
  /* Note : diagonal and inside corner ghost cells are generally     */
  /* only required to set the boundary values for Semi-Lagrangain    */
  /* type advection schemes and do not enter in calculations for     */
  /* most finite difference schemes.                                 */
  for (k = nz - 1; k >= 0; k--) {
    for (j = 0; j < nfe2; j++) {
      for (i = 0; i < nfe1; i++) {

        scen = di_ic(sgrid, i, j, k, &xp1, &xm1, &yp1, &ym1);

        /* Maps from diagonal ghost cells to other ghost cells.      */
        /* Lateral maps from south-west diagonals.                   */
        /* Here (i+1,j+1,k) is wet, (i+1,j,k) & (i,j+1,k) are dry    */
        /* and (i,j,k), the diagonal, is wet.                        */
        /* Lateral maps from south-west inside corners.              */
        /* Here (i+1,j+1,k) is wet, (i+1,j,k) & (i,j+1,k) are dry    */
        /* and (i,j,k), the inside corner, is dry.                   */
        if (i < nce1 && j < nce2 && scen & (SW_D | SW_IC)) {
          c = sgrid->map[k][j + 1][i + 1];
          sgrid->xp1[sw[k][j][i]] = sgrid->ym1[c];
          sgrid->yp1[sw[k][j][i]] = sgrid->xm1[c];
          /* Ghost-ghost vertical maps                               */
          if (k < nz - 1)
            sgrid->zp1[sw[k][j][i]] = sw[k + 1][j][i];
          if (k > kbot[j + 1][i + 1])
            sgrid->zm1[sw[k][j][i]] = sw[k - 1][j][i];
        }

        /* Lateral maps from north-west diagonals.                   */
        /* Here (i+1,j-1,k) is wet, (i+1,j,k) & (i,j-1,k) are dry    */
        /* and (i,j,k), the diagonal, is wet.                        */
        /* Lateral maps from south-east inside corners.              */
        /* Here (i+1,j-1,k) is wet, (i+1,j,k) & (i,j-1,k) are dry    */
        /* and (i,j,k), the inside corner, is dry.                   */
        if (i < nce1 && j > 0 && scen & (NW_D | NW_IC)) {
          c = sgrid->map[k][j - 1][i + 1];
          sgrid->xp1[nw[k][j][i]] = sgrid->yp1[c];
          sgrid->ym1[nw[k][j][i]] = sgrid->xm1[c];
          /* Ghost-ghost vertical maps                               */
          if (k < nz - 1)
            sgrid->zp1[nw[k][j][i]] = nw[k + 1][j][i];
          if (k > kbot[j - 1][i + 1])
            sgrid->zm1[nw[k][j][i]] = nw[k - 1][j][i];
        }

        /* Lateral maps from south-east diagonals.                   */
        /* Here (i-1,j+1,k) is wet, (i-1,j,k) & (i,j+1,k) are dry    */
        /* and (i,j,k), the diagonal, is wet.                        */
        /* Lateral maps from north-west inside corners.              */
        /* Here (i-1,j+1,k) is wet, (i-1,j,k) & (i,j+1,k) are dry    */
        /* and (i,j,k), the inside corner, is dry.                   */
        if (i > 0 && j < nce2 && scen & (SE_D | SE_IC)) {
          c = sgrid->map[k][j + 1][i - 1];
          sgrid->xm1[se[k][j][i]] = sgrid->ym1[c];
          sgrid->yp1[se[k][j][i]] = sgrid->xp1[c];
          /* Ghost-ghost vertical maps                               */
          if (k < nz - 1)
            sgrid->zp1[se[k][j][i]] = se[k + 1][j][i];
          if (k > kbot[j + 1][i - 1])
            sgrid->zm1[se[k][j][i]] = se[k - 1][j][i];
        }

        /* Lateral maps from north-east diagonals.                   */
        /* Here (i-1,j-1,k) is wet, (i-1,j,k) & (i,j-1,k) are dry    */
        /* and (i,j,k), the diagonal, is wet.                        */
        /* Lateral maps from north-east inside corners.              */
        /* Here (i-1,j-1,k) is wet, (i-1,j,k) & (i,j-1,k) are dry    */
        /* and (i,j,k), the inside corner, is dry.                   */
        if (i > 0 && j > 0 && scen & (NE_D | NE_IC)) {
          c = sgrid->map[k][j - 1][i - 1];
          sgrid->xm1[ne[k][j][i]] = sgrid->yp1[c];
          sgrid->ym1[ne[k][j][i]] = sgrid->xp1[c];
          /* Ghost-ghost vertical maps                               */
          if (k < nz - 1)
            sgrid->zp1[ne[k][j][i]] = ne[k + 1][j][i];
          if (k > kbot[j - 1][i - 1])
            sgrid->zm1[ne[k][j][i]] = ne[k - 1][j][i];
        }
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the sediment - wet maps                                     */
  cc = sgrid->sgnum - num_sc + 1;
  for (j = 0; j < nfe2; j++)
    for (i = 0; i < nfe1; i++) {
      if (!(flg[nz - 1][j][i] & NOTWET)) {
        k = kbot[j][i];
        c = sgrid->map[k][j][i];
        sgrid->zm1[c] = cc;
        sgrid->zp1[cc] = c;
	/* Horizontal maps for the sediment */
	if (!(params->compatible & V1598) && (params->trasc & LAGRANGE)) { /* Used for Lagrange scheme */
	  sgrid->xp1[cc] = sgrid->zm1[sgrid->xp1[c]];
	  sgrid->xm1[cc] = sgrid->zm1[sgrid->xm1[c]];
	  sgrid->yp1[cc] = sgrid->zm1[sgrid->yp1[c]];
	  sgrid->ym1[cc] = sgrid->zm1[sgrid->ym1[c]];
	} else {  /* Self mapping */
	  sgrid->xp1[cc] = cc;
	  sgrid->xm1[cc] = cc;
	  sgrid->yp1[cc] = cc;
	  sgrid->ym1[cc] = cc;
	}
        if (gchck == cc && k > 0 )
          hd_warn("%d = sediment location at (%d %d %d)\n", gchck,
               i, j, k - 1);
        cc += 1;
      }
    }

  /*-----------------------------------------------------------------*/
  /* gsed_t are cells which lie below lateral ghosts, adjacent to a  */
  /* sediment cell.                                                  */
  /* ised_t are wet cells diagonally upwards (1 cell up, 1 cell      */
  /* accross) from gsed_t.                                           */
  /* Get the maps for sediment ghost cells (for Lagrange scheme) and */
  /* reset the corresponding vertical maps.                          */
  sgrid->ngsed = num_scg;
  if (sgrid->ngsed) {
    sgrid->gsed_t = i_alloc_1d(sgrid->ngsed + 1);
    sgrid->ised_t = i_alloc_1d(sgrid->ngsed + 1);
    for (i = 1; i <= sgrid->ngsed; i++) {
      sgrid->gsed_t[i] = cc;
      cc += 1;
    }
    cc = sgrid->gsed_t[1];
    /* Vertical maps first */
    for (i = 1; i <= sgrid->nbptS; i++) {
      c = sgrid->bpt[i];
      while (c != sgrid->zm1[c])
	c = sgrid->zm1[c];
      /* Get the corresponding interior cells to gsed_t */
      j = ANY(c, sgrid->bpt, sgrid->nbpt);
      if (j) sgrid->ised_t[i] = sgrid->bin[j];
      sgrid->zm1[c] = cc;
      sgrid->zp1[cc] = c;
      cc += 1;
    }

    /* Horizontal maps */
    for (i = 1; i <= sgrid->ngsed; i++) {
      cc = sgrid->gsed_t[i];    /* Sediment ghost */
      c = sgrid->zp1[cc];       /* Ghost */

      lc = sgrid->zm1[sgrid->xp1[c]];
      sgrid->xp1[cc] = lc;
      if (sgrid->xm1[lc] == lc) sgrid->xm1[lc] = cc;

      lc = sgrid->zm1[sgrid->xm1[c]];
      sgrid->xm1[cc] = lc;
      if (sgrid->xp1[lc] == lc) sgrid->xp1[lc] = cc;

      lc = sgrid->zm1[sgrid->yp1[c]];
      sgrid->yp1[cc] = lc;
      if (sgrid->ym1[lc] == lc) sgrid->ym1[lc] = cc;

      lc = sgrid->zm1[sgrid->ym1[c]];
      sgrid->ym1[cc] = lc;
      if (sgrid->yp1[lc] == lc) sgrid->yp1[lc] = cc;
    }
    num_scg = cc;
  }

  /*-----------------------------------------------------------------*/
  /* Assign a value to ghost cells in the Cartesian-sparse map.      */
  /* Since there exists the possibility that multiple ghost cells    */
  /* exist, a unique sparse location cannot be assigned to           */
  /* sgrid->map, and these locations in the map are assigned values   */
  /* to indicate that the cell exists in sparse space. The corner    */
  /* arrays are used to assign these values.                         */
  for (j = 1; j < nce2; j++) {
    for (i = 1; i < nce1; i++) {
      for (k = nz - 1; k >= 0; k--) {
        if (k < kbot[j][i])
          continue;
        if (!(flg[k][j][i] & NOTWET)) {
          /* Ghost cells on west straight edges                      */
          if (flg[k][j][i - 1] & NOTWET)
            sgrid->map[k][j][i - 1] = sw[k][j][i - 1];
          /* Ghost cells on east straight edges                      */
          if (flg[k][j][i + 1] & NOTWET)
            sgrid->map[k][j][i + 1] = se[k][j][i + 1];
          /* Ghost cells on south straight edges                     */
          if (flg[k][j - 1][i] & NOTWET)
            sgrid->map[k][j - 1][i] = sw[k][j - 1][i];
          /* Ghost cells on north straight edges                     */
          if (flg[k][j + 1][i] & NOTWET)
            sgrid->map[k][j + 1][i] = nw[k][j + 1][i];
          /* Ghost cells on south-west corners                       */
          if (flg[k][j][i - 1] & NOTWET &&
              flg[k][j - 1][i] & NOTWET && flg[k][j - 1][i - 1] & NOTWET)
            sgrid->map[k][j - 1][i - 1] = sw[k][j - 1][i - 1];
          /* Ghost cells on north-west corners                       */
          if (flg[k][j][i - 1] & NOTWET &&
              flg[k][j + 1][i] & NOTWET && flg[k][j + 1][i - 1] & NOTWET)
            sgrid->map[k][j + 1][i - 1] = nw[k][j + 1][i - 1];
          /* Ghost cells on south-east corners                       */
          if (flg[k][j][i + 1] & NOTWET &&
              flg[k][j - 1][i] & NOTWET && flg[k][j - 1][i + 1] & NOTWET)
            sgrid->map[k][j - 1][i + 1] = se[k][j - 1][i + 1];
          /* Ghost cells on north-east corners                       */
          if (flg[k][j][i + 1] & NOTWET &&
              flg[k][j + 1][i] & NOTWET && flg[k][j + 1][i + 1] & NOTWET)
            sgrid->map[k][j + 1][i + 1] = ne[k][j + 1][i + 1];
        }
      }
    }
  }
  /* Count the number of multiple ghost cells                        */
  for (k = nz - 1; k >= 0; k--) {
    for (j = 1; j < nce2; j++) {
      for (i = 1; i < nce1; i++) {
        int a[5];
        a[1] = sw[k][j][i];
        a[2] = se[k][j][i];
        a[3] = nw[k][j][i];
        a[4] = ne[k][j][i];
        for (c = 1; c < 4; c++) {
          for (cc = c + 1; cc <= 4; cc++)
            if (a[c] == a[cc])
              a[cc] = 0;
        }
        cc = 0;
        for (c = 1; c <= 4; c++)
          if (a[c])
            cc++;
        if (cc > 0)
          num_mc += (cc - 1);
      }
    }
  }

  l_free_3d((long ***)sw);
  l_free_3d((long ***)se);
  l_free_3d((long ***)nw);
  l_free_3d((long ***)ne);
  i_free_1d(num_2d);

  /*-----------------------------------------------------------------*/
  /* Save the locations of sparse cells on open boundaries to the    */
  /* boundary arrays.                                                */
  /* The open boundary location vectors are:                         */
  /* obc_t = cell centered sparse location of the open boundary      */
  /* oi1_t = one cell into the interior from obc_t                   */
  /* oi2_t = two cells into the interior from obc_t                  */
  /* obc_e1 = u1 face centered sparse location of the open boundary  */
  /* oi1_e1 = one cell into the interior from obc_e1                 */
  /* oi2_e1 = two cells into the interior from obc_e1                */
  /* obc_e2 = u2 face centered sparse location of the open boundary  */
  /* oi1_e2 = one cell into the interior from obc_e2                 */
  /* oi2_e2 = two cells into the interior from obc_e2                */

  /* Allocate memory for the open boundary data structure            */
  sgrid->open = (open_bdrys_t **)malloc(sizeof(open_bdrys_t *) * sgrid->nobc);

  /* Allocate memory for the boundary structures                     */
  for (n = 0; n < sgrid->nobc; n++) {

    sgrid->open[n] = (open_bdrys_t *)malloc(sizeof(open_bdrys_t));
    memset(sgrid->open[n], 0, sizeof(open_bdrys_t));
    sgrid->open[n]->type = type[n];
    sgrid->open[n]->npts = npts[n];
    sgrid->open[n]->ocodex = ocodex[n];
    sgrid->open[n]->ocodey = ocodey[n];
    sgrid->open[n]->iloc = i_alloc_1d(npts[n]);
    sgrid->open[n]->jloc = i_alloc_1d(npts[n]);
    sgrid->open[n]->bgz = 0;
    for (cc = 0; cc < sgrid->open[n]->npts; cc++) {
      sgrid->open[n]->iloc[cc] = iloc[n][cc];
      sgrid->open[n]->jloc[cc] = jloc[n][cc];
    }
  }
  i_free_1d(ocodex);
  i_free_1d(ocodey);

  /* Set up and save sparse locations for the OBC vectors */
  for (n = 0; n < sgrid->nobc; n++) {

    sgrid->open[n]->ntr = 0;
    set_OBC_cells(sgrid, sgrid->open[n], sgrid->open[n], n, nce1, nce2, nz, flg);

    /* Point the interior cell maps on the boundaries to the correct */
    /* spatial map.                                                  */
    if (sgrid->open[n]->ocodex & L_EDGE) {
      sgrid->open[n]->nmap = sgrid->xp1;
      sgrid->open[n]->omap = sgrid->xm1;
      sgrid->open[n]->tmpp = sgrid->yp1;
      sgrid->open[n]->tmpm = sgrid->ym1;
    }
    if (sgrid->open[n]->ocodex & R_EDGE) {
      sgrid->open[n]->nmap = sgrid->xm1;
      sgrid->open[n]->omap = sgrid->xp1;
      sgrid->open[n]->tmpp = sgrid->yp1;
      sgrid->open[n]->tmpm = sgrid->ym1;
    }
    if (sgrid->open[n]->ocodey & B_EDGE) {
      sgrid->open[n]->nmap = sgrid->yp1;
      sgrid->open[n]->omap = sgrid->ym1;
      sgrid->open[n]->tmpp = sgrid->xp1;
      sgrid->open[n]->tmpm = sgrid->xm1;
    }
    if (sgrid->open[n]->ocodey & F_EDGE) {
      sgrid->open[n]->nmap = sgrid->ym1;
      sgrid->open[n]->omap = sgrid->yp1;
      sgrid->open[n]->tmpp = sgrid->xp1;
      sgrid->open[n]->tmpm = sgrid->xm1;
    }
  }

  /* Set maps to be self-pointing over R_EDGE and F_EDGE OBCs        */
  set_map_bdry(sgrid);

  /* Get the cyclic boundary locations.                              */
  for (n = 0; n < sgrid->nobc; n++) {
    /* Tracers                                                       */
    for (cc = 1; cc <= sgrid->open[n]->no3_t; cc++) {
      c = isvalidc(sgrid->open[n]->obc_t[cc], sgrid->sgnum, "tracer OBCs");
      sgrid->open[n]->cyc_t[cc] = cyclic_m2(sgrid, sgrid->open[n]->ocodex,
					   sgrid->open[n]->ocodey,
					   sgrid->open[n]->nmap, c);
    }

    /* u1 velocity                                                   */
    /* Normal component                                              */
    for (cc = 1; cc <= sgrid->open[n]->no3_e1; cc++) {
      c = isvalidc(sgrid->open[n]->obc_e1[cc], sgrid->sgnum, 
		   "e1 normal velocity OBCs");
      sgrid->open[n]->cyc_e1[cc] = cyclic_m1(sgrid, sgrid->open[n]->ocodex,
					    sgrid->open[n]->ocodey,
					    sgrid->open[n]->nmap, c);
    }

    /* Tangential component                                          */
    for (cc = sgrid->open[n]->no3_e1 + 1; cc <= sgrid->open[n]->to3_e1; cc++) {
      c = isvalidc(sgrid->open[n]->obc_e1[cc], sgrid->sgnum, 
		   "e1 tangential velocity OBCs");
      sgrid->open[n]->cyc_e1[cc] = cyclic_m2(sgrid, sgrid->open[n]->ocodex,
					    sgrid->open[n]->ocodey,
					    sgrid->open[n]->nmap, c);
    }

    /* u2 velocity                                                   */
    /* Normal component                                              */
    for (cc = 1; cc <= sgrid->open[n]->no3_e2; cc++) {
      c = isvalidc(sgrid->open[n]->obc_e2[cc], sgrid->sgnum, 
		   "e2 normal velocity OBCs");
      sgrid->open[n]->cyc_e2[cc] = cyclic_m1(sgrid, sgrid->open[n]->ocodex,
					    sgrid->open[n]->ocodey,
					    sgrid->open[n]->nmap, c);
    }

    /* Tangential component                                          */
    for (cc = sgrid->open[n]->no3_e2 + 1; cc <= sgrid->open[n]->to3_e2; cc++) {
      c = isvalidc(sgrid->open[n]->obc_e2[cc], sgrid->sgnum, 
		   "e2 tangential velocity OBCs");
      sgrid->open[n]->cyc_e2[cc] = cyclic_m2(sgrid, sgrid->open[n]->ocodex,
					    sgrid->open[n]->ocodey,
					    sgrid->open[n]->nmap, c);
    }

    /* Get the bottom coordinate vector                              */
    sgrid->open[n]->bot_t = i_alloc_1d(sgrid->open[n]->no2_t + 1);
    for (cc = 1; cc <= sgrid->open[n]->no2_t; cc++) {

      c = c2 = sgrid->open[n]->obc_t[cc];
      while(c != sgrid->zm1[c])
	c = sgrid->zm1[c];
      sgrid->open[n]->bot_t[cc] = sgrid->zp1[c];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the OBC ghost cell maps.                                    */
  for (n = 0; n < sgrid->nobc; n++) {
    open_bdrys_t *open = sgrid->open[n];
    int nlaux = laux;

    if ((open->type & U1BDRY && open->ocodex & R_EDGE) ||
	(open->type & U2BDRY && open->ocodey & F_EDGE))
      nlaux--;

    if (!(params->compatible & V1957)) {
      if (open->bgz) {
	for (cc = 0; cc < open->npts; cc++) {
	  i = open->iloc[cc];
	  j = open->jloc[cc];
	  for (k = kbot[j][i]; k < nz; k++) {
	    c = c2 = sgrid->map[k][j][i];
	    c1 = end_lc[k];
	    for (ii = 0; ii < nlaux; ii++) {
	      sgrid->open[n]->omap[c] = c1;
	      sgrid->open[n]->nmap[c1] = c;
	      sgrid->wgst[c1] = c2;
	      if (gchck == c1 && DEBUG("init_m")) {
		if (open->ocodex & L_EDGE)
		  dlog("init_m", "%d = OBC ghost location at (%d %d %d)\n", gchck, i-ii-1, j, k);
		if (open->ocodex & R_EDGE)
		  dlog("init_m", "%d = OBC ghost location at (%d %d %d)\n", gchck, i+ii+1, j, k);
		if (open->ocodey & B_EDGE)
		  dlog("init_m", "%d = OBC ghost location at (%d %d %d)\n", gchck, i, j-ii-1, k);
		if (open->ocodey & F_EDGE)
		  dlog("init_m", "%d = OBC ghost location at (%d %d %d)\n", gchck, i, j+ii+1, k);
	      }
	      c = c1;
	      c1++;
	    }
	    end_lc[k] += nlaux;
	  }
	}
      }

      /* Get the boundary ghost cells if required */
      if (open->ocodex & R_EDGE)
	open->ogc_t = open->obc_e1;
      else if (open->ocodey & F_EDGE)
	open->ogc_t = open->obc_e2;
      else {
	open->ogc_t = i_alloc_1d(open->no3_t + 1);
	for (cc = 1; cc <= open->no3_t; cc++) {
	  c = open->obc_t[cc];
	  open->ogc_t[cc] = open->omap[c];
	}
      }
    }
  }
  i_free_1d(end_lc);

  /*-----------------------------------------------------------------*/
  /* Get the wet cell mask : note do not include the edges at nfe1   */
  /* and nfe2 in this mask.                                          */
  /* Separate wet point vectors are required for tracers and each    */
  /* velocity component, since a cell may contain both land and      */
  /* water locations depending on whether the cell face or center is */
  /* considered.                                                     */
  /* v3_t = number of non-boundary cell centered wet cells           */
  /* b3_t = number of cell centered wet cells                        */
  /* a3_t = number of cell centered wet cells including cell centers */
  /* corresponding to e1 and e2 OBC cells. Note: this may not        */
  /* be num_wc since OUTSIDE cells whose faces are open              */
  /* boundaries are also included.                                   */
  /* n2_t = number of cell centered wet + ghost cells                */
  /* v3_e1 = number of non-boundary wet cells on u1 faces            */
  /* v3_e2 = number of non-boundary wet cells on u2 faces            */
  /* v2_e1 = number of non-boundary wet cells on u1av faces          */
  /* v2_e2 = number of non-boundary wet cells on u2av faces          */
  /* w3_t = global cells to process vector for tracers               */
  /* w3_e1 = global cells to process vector for u1 velocity          */
  /* w3_e2 = global cells to process vector for u2 velocity          */
  /* wsa=global cell centered sparse locations listed consecutively, */
  /* 1 to a3_t with 2D cells listed first a2_t cells being the       */
  /* surface.                                                        */

  /* Initialise vector positions not used in the global vectors to   */
  /* zero.                                                           */
  sgrid->a3_t = 0;
  sgrid->a2_t = 0;

  /* Count the number of cell centered locations that are wet and    */
  /* not adjacent to a u1 or u2 boundary. Note: this is not simply   */
  /* the quantity num_wc-nu1bdry-nu2bdry since some cells are both   */
  /* u1 and u2 boundaries (e.g. corner cells).                       */
  sgrid->v3_t = sgrid->v2_t = 0;
  sgrid->b3_t = sgrid->b2_t = 0;
  b2 = b3 = 0;
  for (k = nz - 1; k >= 0; k--) {
    for (j = 0; j < nfe2; j++) {
      for (i = 0; i < nfe1; i++) {
        int fi = i + 1;
        int fj = j + 1;
        /* Count the e1 or e2 boundary cells whose cell centers are  */
        /* OUTSIDE cells.                                            */
        if (flg[k][j][i] & OUTSIDE &&
            (flg[k][j][i] & U1BDRY || flg[k][j][i] & U2BDRY)) {
          sgrid->a3_t++;
          if (k == nz - 1)
            sgrid->a2_t++;
        }
        /* Count the cells centers which are wet                     */
        if (!(flg[k][j][i] & (NOTWET))) {
          /* Count wet cell centers on the whole nfe1 x nfe2 grid    */
          /* (i.e. sparse cells which may have a valid tracer, e1 or */
          /* e2 velocity value).                                     */
          sgrid->a3_t++;
          if (k == nz - 1)
            sgrid->a2_t++;
          if (i < nce1 && j < nce2) {
            /* Count wet cell centers on the whole nce1 x nce2 grid  */
            /* (i.e. valid tracer / elevation cells).                */
	    sgrid->b3_t++;
	    if (k == nz - 1)
	      sgrid->b2_t++;
            /* Count non-boundary wet cell centers                   */
            if (!(flg[k][j][i] & (U1BDRY | U2BDRY)) &&
                !(flg[k][fj][i] & U2BDRY) && !(flg[k][j][fi] & U1BDRY)) {
              sgrid->v3_t++;
              if (k == nz - 1)
                sgrid->v2_t++;
            }
          }
        }
      }
    }
  }

  /* Count the number of wet cells for e1 and e2 velocity            */
  sgrid->v3_e1 = sgrid->v2_e1 = sgrid->b3_e1 = sgrid->b2_e1 = 0;
  sgrid->v3_e2 = sgrid->v2_e2 = sgrid->b3_e2 = sgrid->b2_e2 = 0;
  for (k = nz - 1; k >= 0; k--) {
    for (j = 0; j < nfe2; j++) {
      for (i = 0; i < nfe1; i++) {
        if (j < nce2 && !(flg[k][j][i] & (U1SOLID | U1OUTSIDE))) {
	  sgrid->b3_e1 += 1;
          if (k == nz - 1)
            sgrid->b2_e1++;
          if (!(flg[k][j][i] & (U1BDRY))) {
	    sgrid->v3_e1 += 1;
            if (k == nz - 1)
              sgrid->v2_e1++;
          }
        }
        if (i < nce1 && !(flg[k][j][i] & (U2SOLID | U2OUTSIDE))) {
	  sgrid->b3_e2 += 1;
          if (k == nz - 1)
            sgrid->b2_e2++;
          if (!(flg[k][j][i] & (U2BDRY))) {
	    sgrid->v3_e2 += 1;
            if (k == nz - 1)
              sgrid->v2_e2++;
          }
        }
      }
    }
  }

  /* Allocate the cells to process vectors                           */
  sgrid->w3_t = i_alloc_1d(sgrid->a3_t + sgrid->nbpt + 1);
  sgrid->w3_e1 = i_alloc_1d(sgrid->b3_e1 + sgrid->nbpt + 1);
  sgrid->w3_e2 = i_alloc_1d(sgrid->b3_e2 + sgrid->nbpt + 1);
  sgrid->wsa = i_alloc_1d(sgrid->a3_t + 1);

  sgrid->w2_t = i_alloc_1d(sgrid->a2_t + sgrid->nbptS + 1);
  sgrid->w2_e1 = i_alloc_1d(sgrid->b2_e1 + sgrid->nbptS + 1);
  sgrid->w2_e2 = i_alloc_1d(sgrid->b2_e2 + sgrid->nbptS + 1);

  sgrid->bot_t = i_alloc_1d(sgrid->b2_t + 1);
  sgrid->bot_e1 = i_alloc_1d(sgrid->b2_e1 + 1);
  sgrid->bot_e2 = i_alloc_1d(sgrid->b2_e2 + 1);
  sgrid->sur_t = i_alloc_1d(sgrid->b2_t + 1);
  c = 0;
  for (j = 0; j < nfe2; j++)
    for (i = 0; i < nfe1; i++)
      if (sgrid->map[nz - 1][j][i] > (unsigned long)c)
        c = sgrid->map[nz - 1][j][i];
  sgrid->c2cc = i_alloc_1d(c + 1);

  /* Fill the cell centered vector with non-boundary cell locations  */
  sgrid->a3_t = sgrid->b3_t + 1;
  sgrid->a2_t = sgrid->b2_t + 1;
  sgrid->b3_t = sgrid->v3_t + 1;
  sgrid->b2_t = sgrid->v2_t + 1;
  sgrid->v3_t = 1;
  sgrid->v2_t = 1;
  c1 = 1;
  for (k = nz - 1; k >= 0; k--) {
    for (j = 0; j < nfe2; j++) {
      for (i = 0; i < nfe1; i++) {
        int fi = i + 1;
        int fj = j + 1;
        c = sgrid->map[k][j][i];
        /* Map the sparse grid, 2D cells first followed by 3D cells  */
        if ((!(flg[k][j][i] & (NOTWET))) ||
            (flg[k][j][i] & OUTSIDE && (flg[k][j][i] & U1BDRY ||
                                        flg[k][j][i] & U2BDRY))) {
          sgrid->wsa[c1] = c;
          c1++;
        }
        /* Map the e1 or e2 boundary cells whose cell centers lie on */
        /* the nce1 or nce2 limits of the grid.                      */
        if (!(flg[k][j][i] & (NOTWET)) && (i == nce1 || j == nce2)) {
          sgrid->w3_t[sgrid->a3_t] = c;
          sgrid->a3_t += 1;
          if (k == nz - 1) {
            sgrid->w2_t[sgrid->a2_t] = c;
            sgrid->a2_t++;
          }
        }
        /* Map the cell centered wet cells                           */
        if (!(flg[k][j][i] & (NOTWET)) && i < nce1 && j < nce2) {
          /* Cell centers that do not lie on an open boundary        */
          if (!(flg[k][j][i] & (U1BDRY | U2BDRY)) &&
              !(flg[k][fj][i] & U2BDRY) && !(flg[k][j][fi] & U1BDRY)) {
            sgrid->w3_t[sgrid->v3_t] = c;
            sgrid->v3_t += 1;
            if (k == nz - 1) {
              sgrid->w2_t[sgrid->v2_t] = c;
              c2 = sgrid->map[kbot[j][i]][j][i];
              sgrid->bot_t[sgrid->v2_t] = c2;
              sgrid->v2_t++;
            }
          }
          /* Cell centers on open boundaries                         */
	  else {
            sgrid->w3_t[sgrid->b3_t] = c;
            sgrid->b3_t += 1;
            if (k == nz - 1) {
              sgrid->w2_t[sgrid->b2_t] = c;
              c2 = sgrid->map[kbot[j][i]][j][i];
              sgrid->bot_t[sgrid->b2_t] = c2;
              sgrid->b2_t++;
            }
          }
        }
      }
    }
  }
  sgrid->v3_t -= 1;
  sgrid->v2_t -= 1;
  sgrid->b3_t -= 1;
  sgrid->b2_t -= 1;
  sgrid->a3_t -= 1;
  sgrid->a2_t -= 1;
  sgrid->ns2 = sgrid->a2_t;
  sgrid->ns3 = sgrid->a3_t;

  for (cc = 1; cc <= sgrid->b2_t; cc++) {
    c = sgrid->w2_t[cc];
    sgrid->sur_t[cc] = c;
    sgrid->c2cc[c] = cc;
  }

  /* e1 face centered wet cells. This must be done in the order of x */
  /* direction first, followed by y then z.                          */
  sgrid->b3_e1 = sgrid->v3_e1 + 1;
  sgrid->b2_e1 = sgrid->v2_e1 + 1;
  sgrid->v3_e1 = 1;
  sgrid->v2_e1 = 1;
  for (k = nz - 1; k >= 0; k--) {
    for (j = 0; j < nce2; j++) {
      for (i = 0; i < nfe1; i++) {
        c = sgrid->map[k][j][i];
        if (!(flg[k][j][i] & (U1SOLID | U1OUTSIDE))) {
          if (!(flg[k][j][i] & U1BDRY)) {
	    sgrid->w3_e1[sgrid->v3_e1] = c;
	    sgrid->v3_e1++;
            if (k == nz - 1) {
              sgrid->w2_e1[sgrid->v2_e1] = c;
              sgrid->v2_e1++;
            }
          } else {
	    sgrid->w3_e1[sgrid->b3_e1] = c;
	    sgrid->b3_e1++;
            if (k == nz - 1) {
              sgrid->w2_e1[sgrid->b2_e1] = c;
              sgrid->b2_e1++;
            }
          }
        }
      }
    }
  }
  sgrid->v3_e1 -= 1;
  sgrid->v2_e1 -= 1;
  sgrid->b3_e1 -= 1;
  sgrid->b2_e1 -= 1;

  /* e2 face centered wet cells.                                     */
  sgrid->b3_e2 = sgrid->v3_e2 + 1;
  sgrid->b2_e2 = sgrid->v2_e2 + 1;
  sgrid->v3_e2 = 1;
  sgrid->v2_e2 = 1;
  for (k = nz - 1; k >= 0; k--) {
    for (j = 0; j < nfe2; j++) {
      for (i = 0; i < nce1; i++) {
        c = sgrid->map[k][j][i];
        if (!(flg[k][j][i] & (U2SOLID | U2OUTSIDE))) {
          if (!(flg[k][j][i] & U2BDRY)) {
	    sgrid->w3_e2[sgrid->v3_e2] = c;
	    sgrid->v3_e2++;
            if (k == nz - 1) {
              sgrid->w2_e2[sgrid->v2_e2] = c;
              sgrid->v2_e2++;
            }
          } else {
	    sgrid->w3_e2[sgrid->b3_e2] = c;
	    sgrid->b3_e2++;
            if (k == nz - 1) {
              sgrid->w2_e2[sgrid->b2_e2] = c;
              sgrid->b2_e2++;
            }
          }
        }
      }
    }
  }
  sgrid->v3_e2 -= 1;
  sgrid->v2_e2 -= 1;
  sgrid->b3_e2 -= 1;
  sgrid->b2_e2 -= 1;

  /* Fill the bottom vectors with the sparse coordinate of the       */
  /* bottom. This must be done in the order of z direction first,    */
  /* followed by x then y (and hence cannot be incorporated into the */
  /* wet cell mask loop above).                                      */
  b1 = b2 = 1;
  c1 = sgrid->v2_e1 + 1;
  c2 = sgrid->v2_e2 + 1;
  for (j = 0; j < nfe2; j++) {
    for (i = 0; i < nfe1; i++) {
      for (k = nz - 1; k >= 0; k--) {
        c = sgrid->map[k][j][i];
        if (!(flg[k][j][i] & NOTWET)) {
          if (j < nce2 && !(flg[k][j][i] & (U1SOLID | U1OUTSIDE))) {
            if (k == 0) {
              if (!(flg[k][j][i] & U1BDRY)) {
                sgrid->bot_e1[b1] = c;
                b1++;
              } else {
                sgrid->bot_e1[c1] = c;
                c1++;
              }
            } else {
              if (flg[k - 1][j][i] & (U1SOLID | U1OUTSIDE)) {
                if (!(flg[k][j][i] & U1BDRY)) {
                  sgrid->bot_e1[b1] = c;
                  b1++;
                } else {
                  sgrid->bot_e1[c1] = c;
                  c1++;
                }
              }
            }
          }
          if (i < nce1 && !(flg[k][j][i] & (U2SOLID | U2OUTSIDE))) {
            if (k == 0) {
              if (!(flg[k][j][i] & U2BDRY)) {
                sgrid->bot_e2[b2] = c;
                b2++;
              } else {
                sgrid->bot_e2[c2] = c;
                c2++;
              }
            } else {
              if (flg[k - 1][j][i] & (U2SOLID | U2OUTSIDE)) {
                if (!(flg[k][j][i] & U2BDRY)) {
                  sgrid->bot_e2[b2] = c;
                  b2++;
                } else {
                  sgrid->bot_e2[c2] = c;
                  c2++;
                }
              }
            }
          }
        }
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the sparse to Cartesian maps                                */
  for (c = 1; c <= sgrid->sgnum; c++) {
    sgrid->s2i[c] = NOTVALID;
    sgrid->s2j[c] = NOTVALID;
    sgrid->s2k[c] = NOTVALID;
  }
  for (k = nz - 1; k >= 0; k--)
    for (j = 0; j < nfe2; j++)
      for (i = 0; i < nfe1; i++) {
        if (!(flg[k][j][i] & NOTWET)) {
          c = sgrid->map[k][j][i];
          sgrid->s2i[c] = i;
          sgrid->s2j[c] = j;
          sgrid->s2k[c] = k;
        }
      }

  /* Sparse - Cartesian maps for the sediment                        */
  c = sgrid->sgnum - num_sc + 1;
  for (j = 0; j < nfe2; j++)
    for (i = 0; i < nfe1; i++) {
      if (!(flg[nz - 1][j][i] & NOTWET)) {
        k = kbot[j][i] - 1;
        sgrid->s2i[c] = i;
        sgrid->s2j[c] = j;
        if (k >= 0) {
          sgrid->map[k][j][i] = c;
          sgrid->s2k[c] = k;
        }
	if (k == -1) {
	  sgrid->s2i[c] = NOTVALID;
	  sgrid->s2j[c] = NOTVALID;
	}
        c++;
      }
    }

  /*-----------------------------------------------------------------*/
  /* Get the lateral boundary masks.                                 */
  /* nbpt = number of lateral boundary cells (ghost cells)           */
  /* nbptS = number of ghost cells in the 2D grid                    */
  /* bpt = cell centered lateral boundary mask                       */
  /* bin = one point into the interior from bpt                      */
  /* bpte1 = u1 face centered lateral boundary mask                  */
  /* bine1 = one point into the interior from bpte1                  */
  /* bpte2 = u2 face centered lateral boundary mask                  */
  /* bine2 = one point into the interior from bpte2                  */
  sgrid->bin2 = i_alloc_1d(sgrid->nbpt + 1);
  sgrid->bpte1 = i_alloc_1d(sgrid->nbpt + 1);
  sgrid->bpte2 = i_alloc_1d(sgrid->nbpt + 1);
  sgrid->bine1 = i_alloc_1d(sgrid->nbpt + 1);
  sgrid->bine2 = i_alloc_1d(sgrid->nbpt + 1);
  sgrid->bpte1S = i_alloc_1d(sgrid->nbptS + 1);
  sgrid->bpte2S = i_alloc_1d(sgrid->nbptS + 1);
  sgrid->bine1S = i_alloc_1d(sgrid->nbptS + 1);
  sgrid->bine2S = i_alloc_1d(sgrid->nbptS + 1);
  memset(sgrid->wgst, 0, sgrid->sgsiz * sizeof(int));

#if defined(HAVE_SEDIMENT_MODULE)
  sgrid->sed_t = i_alloc_1d(sgrid->ewetS + 1);
#endif

  /* Count the lateral ghost cells which are associated with normal  */
  /* velocities. These are cells on eastern and western edges for u1 */
  /* velocities and northern and southern edges for u2 velocities.   */
  /* These land boundary locations are organised so that the first   */
  /* sgrid->be#n cells contain the normal velocities (which have a    */
  /* zero flux boundary condition) followed by the tangential        */
  /* velocities (which have a free/no slip boundary condition.       */
  sgrid->nbe1 = sgrid->nbe1S = 0;
  sgrid->nbe2 = sgrid->nbe2S = 0;
  for (cc = 1; cc <= sgrid->nbpt; cc++) {
    c = sgrid->bpt[cc];
    lc = sgrid->bin[cc];
    sgrid->wgst[c] = lc;
    if (sgrid->xp1[lc] == c || sgrid->xm1[lc] == c) {
      sgrid->nbe1 += 1;
      if (cc <= sgrid->nbptS)
        sgrid->nbe1S += 1;
    }
    if (sgrid->yp1[lc] == c || sgrid->ym1[lc] == c) {
      sgrid->nbe2 += 1;
      if (cc <= sgrid->nbptS)
        sgrid->nbe2S += 1;
    }
  }

  /* Get the locations of the lateral ghost cells                    */
  for (cc = 1; cc <= sgrid->ngsed; cc++) {
    c = sgrid->gsed_t[cc];    /* Sediment ghost */
    sgrid->wgst[c] = sgrid->ised_t[cc];
  }
  b1 = sgrid->nbe1 + 1;
  b2 = sgrid->nbe2 + 1;
  b3 = sgrid->nbe1S + 1;
  c3 = sgrid->nbe2S + 1;
  c1 = c2 = 1;
  for (cc = 1; cc <= sgrid->nbpt; cc++) {
    c = sgrid->bpt[cc];
    lc = sgrid->bin[cc];

    if (sgrid->xm1[lc] == c) {
      sgrid->bin2[cc] = sgrid->xp1[lc];
      /* Do not use the staggered e1 cell for locations where an     */
      /* open boundary is adjacent to an OUTSIDE cell in the L_EDGE. */
      i = sgrid->s2i[lc];
      j = sgrid->s2j[lc];
      k = sgrid->s2k[lc];
      if (flg[k][j][i] & U1BDRY)
        sgrid->bpte1[c1] = c;
      else
        sgrid->bpte1[c1] = sgrid->xp1[c];
      sgrid->bine1[c1] = sgrid->xp1[sgrid->bpte1[c1]];
      if (cc <= sgrid->nbptS) {
        sgrid->bpte1S[c1] = sgrid->bpte1[c1];
        sgrid->bine1S[c1] = sgrid->bine1[c1];
      }
      c1++;
    } else if (sgrid->xp1[lc] == c) {
      sgrid->bin2[cc] = sgrid->xm1[lc];
      sgrid->bpte1[c1] = c;
      sgrid->bine1[c1] = lc;
      if (cc <= sgrid->nbptS) {
        sgrid->bpte1S[c1] = sgrid->bpte1[c1];
        sgrid->bine1S[c1] = sgrid->bine1[c1];
      }
      c1++;
    } else {
      sgrid->bpte1[b1] = c;
      sgrid->bine1[b1] = lc;
      if (cc <= sgrid->nbptS) {
        sgrid->bpte1S[b3] = sgrid->bpte1[b1];
        sgrid->bine1S[b3] = sgrid->bine1[b1];
        b3++;
      }
      b1++;
    }

    if (sgrid->ym1[lc] == c) {
      sgrid->bin2[cc] = sgrid->yp1[lc];
      /* Do not use the staggered e2 cell for locations where an     */
      /* open boundary is adjacent to an OUTSIDE cell in the B_EDGE. */
      i = sgrid->s2i[lc];
      j = sgrid->s2j[lc];
      k = sgrid->s2k[lc];
      if (flg[k][j][i] & U2BDRY)
        sgrid->bpte2[c2] = c;
      else
        sgrid->bpte2[c2] = sgrid->yp1[c];
      sgrid->bine2[c2] = sgrid->yp1[sgrid->bpte2[c2]];
      if (cc <= sgrid->nbptS) {
        sgrid->bpte2S[c2] = sgrid->bpte2[c2];
        sgrid->bine2S[c2] = sgrid->bine2[c2];
      }
      c2++;
    } else if (sgrid->yp1[lc] == c) {
      sgrid->bin2[cc] = sgrid->ym1[lc];
      sgrid->bpte2[c2] = c;
      sgrid->bine2[c2] = lc;
      if (cc <= sgrid->nbptS) {
        sgrid->bpte2S[c2] = sgrid->bpte2[c2];
        sgrid->bine2S[c2] = sgrid->bine2[c2];
      }
      c2++;
    } else {
      sgrid->bpte2[b2] = c;
      sgrid->bine2[b2] = lc;
      if (cc <= sgrid->nbptS) {
        sgrid->bpte2S[c3] = sgrid->bpte2[b2];
        sgrid->bine2S[c3] = sgrid->bine2[b2];
        c3++;
      }
      b2++;
    }
  }

  /* Fill the cells to process vectors with ghost cells.             */
  /* Cell centered vectors.                                          */
  sgrid->n3_t = sgrid->a3_t + 1;
  sgrid->n2_t = sgrid->a2_t + 1;
  for (cc = 1; cc <= sgrid->nbpt; cc++) {
    c = sgrid->bpt[cc];
    sgrid->w3_t[sgrid->n3_t] = c;
    sgrid->n3_t++;
    if (cc <= sgrid->nbptS) {
      sgrid->w2_t[sgrid->n2_t] = c;
      sgrid->n2_t++;
    }
  }
  /* e1 face centered vectors                                        */
  sgrid->n3_e1 = sgrid->b3_e1 + 1;
  sgrid->n2_e1 = sgrid->b2_e1 + 1;
  for (cc = 1; cc <= sgrid->nbpt; cc++) {
    c = sgrid->bpte1[cc];
    sgrid->w3_e1[sgrid->n3_e1] = c;
    sgrid->n3_e1++;
    if (cc <= sgrid->nbptS) {
      sgrid->w2_e1[sgrid->n2_e1] = c;
      sgrid->n2_e1++;
    }
  }
  /* e2 face centered vectors                                        */
  sgrid->n3_e2 = sgrid->b3_e2 + 1;
  sgrid->n2_e2 = sgrid->b2_e2 + 1;
  for (cc = 1; cc <= sgrid->nbpt; cc++) {
    c = sgrid->bpte2[cc];
    sgrid->w3_e2[sgrid->n3_e2] = c;
    sgrid->n3_e2++;
    if (cc <= sgrid->nbptS) {
      sgrid->w2_e2[sgrid->n2_e2] = c;
      sgrid->n2_e2++;
    }
  }
  sgrid->n3_t -= 1;
  sgrid->n2_t -= 1;
  sgrid->n3_e1 -= 1;
  sgrid->n2_e1 -= 1;
  sgrid->n3_e2 -= 1;
  sgrid->n2_e2 -= 1;

  /*-----------------------------------------------------------------*/
  /* Get the 3D - 2D map (including the sediment)                    */
  sgrid->m2d = i_alloc_1d(sgrid->sgsiz);
  for (cc = 1; cc <= sgrid->enon; cc++)
    sgrid->m2d[cc] = 0;
  for (cc = 1; cc <= sgrid->enonS; cc++) {
    c = cc;
    while (sgrid->zm1[c] && c != sgrid->zm1[c]) {
      sgrid->m2d[c] = cc;
      c = sgrid->zm1[c];
    }
    sgrid->m2d[c] = cc;
  }

  /* Get the 3D - 2D map for sub-surface ghost cells (i.e. ghost     */
  /* cells where the cell above is wet).                             */
  for (cc = 1; cc <= sgrid->enon; cc++) {
    c = cc;
    c1 = sgrid->zp1[cc];
    while (sgrid->m2d[cc] == 0 && c != c1) {
      sgrid->m2d[cc] = sgrid->m2d[c1];
      c = c1;
      c1 = sgrid->zp1[c1];
    }
  }

  /* Map boundary ghosts if required */
  for (n = 0; n < sgrid->nobc; n++) {
    open_bdrys_t *open = sgrid->open[n];
    if (params->compatible & V1957) continue;
    for (cc = 1; cc <= open->no3_t; cc++) {
      c = open->obc_t[cc];
      c1 = sgrid->m2d[c];
      for(ii = 0; ii < open->bgz; ii++) {
	c = open->omap[c];
	c1 = open->omap[c1];
	sgrid->m2d[c] = c1;
      }
    }
  }

  check_sparse(sgrid, sgrid->a3_t, num_gc, num_sc, num_se, num_ic, num_dc,
               num_mc, flg);

  /*-----------------------------------------------------------------*/
  /* Free memory                                                     */
  l_free_3d((long ***)flg);
}

/* END build_sparse_grid()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to allocate memory from the geometry data structure       */
/* arrays.                                                           */
/*-------------------------------------------------------------------*/
void alloc_sgrid(geometry_t *sgrid /* Model geometry structure     */
  )
{
  int size;                        /* Size of the 3D sparse array    */
  int sizeS;                       /* Size of the 2D sparse array    */

  size = sgrid->enon + 1;
  sizeS = sgrid->enonS + 1;

  sgrid->xp1 = i_alloc_1d(size);
  sgrid->yp1 = i_alloc_1d(size);
  sgrid->zp1 = i_alloc_1d(size);
  sgrid->xm1 = i_alloc_1d(size);
  sgrid->ym1 = i_alloc_1d(size);
  sgrid->zm1 = i_alloc_1d(size);
  sgrid->s2i = i_alloc_1d(size);
  sgrid->s2j = i_alloc_1d(size);
  sgrid->s2k = i_alloc_1d(size);
  sgrid->xmyp1 = i_alloc_1d(size);
  sgrid->xpym1 = i_alloc_1d(size);
}

/* END alloc_sgrid()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to allocate memory from the geometry data structure       */
/* arrays.                                                           */
/*-------------------------------------------------------------------*/
void clear_sgrid(geometry_t *sgrid   /* Model geometry structure     */
  )
{
  int n;

  i_free_1d(sgrid->xp1);
  i_free_1d(sgrid->yp1);
  i_free_1d(sgrid->zp1);
  i_free_1d(sgrid->xm1);
  i_free_1d(sgrid->ym1);
  i_free_1d(sgrid->zm1);
  i_free_1d(sgrid->xmyp1);
  i_free_1d(sgrid->xpym1);
  i_free_1d(sgrid->s2i);
  i_free_1d(sgrid->s2j);
  i_free_1d(sgrid->s2k);
  i_free_1d(sgrid->wsa);
  for (n = 0; n < sgrid->nobc; n++) {
    i_free_1d(sgrid->open[n]->iloc);
    i_free_1d(sgrid->open[n]->jloc);
    free((open_bdrys_t **)sgrid->open[n]);
  }
  free((geometry_t *)sgrid);
}

/* END clear_sgrid()                                                 */
/*-------------------------------------------------------------------*/


