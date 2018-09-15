/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/master/preprocess.c
 *  
 *  Description: SHOC preprocessor
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: preprocess.c 5901 2018-08-28 02:10:22Z riz008 $
 *
 */

#include <stdlib.h>
#include <string.h>
#include "hd.h"
#include "tracer.h"

/*-------------------------------------------------------------------*/
/* Routines contained in this module                                 */
/* build_sparse_map(): Generates the global sparse coordinate system */
/* check_sparse() : Checks the sparse system is OK                   */
/* se_oc()        : Locates straight edges and outside corners       */
/* di_ic()        : Locates diagonals and inside corners             */
/* get_sloc()     : Finds (i,j,k) from sparse coordinates            */
/*-------------------------------------------------------------------*/

#define SIZETYPE	size_t

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

tide_details_t *sreadtideforce(FILE * fp, int bnum, int n);
int *set_explicit_maps(geometry_t *geom, int nem, int *emis, int *emjs,
		       int *emid, int *emjd, int *kti, int *kbi,
		       unsigned long ***flag, short **kbot, int *p1,
		       int *m1, char *tag, int *emapf);
void reset_points(parameters_t *params, unsigned long **flag, 
		  int *nlist, int *listi, int *listj);
int iswetijk(int c, geometry_t *geom, unsigned long ***flag);

typedef struct sparse_grid sparse_grid_t;
typedef struct obc_bdry obc_bdry_t;


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
void build_sparse_map(parameters_t *params  /* Input parameters data
                                               structure */
  )
{
  int c = 1;                    /* Locations of wet points in the sparse
                                   array */
  int cc;                       /* Locations of ghost points in sparse
                                   array */
  int lc;                       /* Loc of boundary points in the boundary
                                   array */
  int i, j, k, ii, jj;          /* x,y,z counters */
  int n;                        /* General counters */
  int xp1;                      /* xp1 = geom->map[k][j][i+1] */
  int xm1;                      /* xm1 = geom->map[k][j][i-1] */
  int yp1;                      /* yp1 = geom->map[k][j+1][i] */
  int ym1;                      /* ym1 = geom->map[k][j-1][i] */
  int *end_lc;                  /* End loaction of the wet points in each
                                   layer */
  int *num_2d;                  /* Number of ghost points in each layer */
  int num_dc = 0;               /* Number of diagonal connection ghost
                                   points */
  int num_se = 0;               /* # straight edge/outside corner ghost
                                   points */
  int num_ic = 0;               /* Number of inside corner ghost points */
  int num_mc = 0;               /* Number of multiple ghost cells */
  int num_wc = 0;               /* Number of 3D wet grid cells */
  int num_gc = 0;               /* Number of ghost grid cells */
  int num_sc = 0;               /* Number of sediment grid cells */
  int num_scg = 0;              /* Number of sediment ghost grid cells */
  int num_gc2D = 0;             /* Number of 2D ghost grid cells */
  int num_obc = 0;              /* Number of OBC ghost cells */
  int num_obc2D = 0;            /* Number of 2D OBC ghost cells */
  short **kbot;                 /* k index of the bottom */
  int c1, c2, c3;               /* Counters */
  int b1, b2, b3;               /* Counters */
  long scen;                    /* Flag for ghost cell type */
  unsigned long ***sw;          /* South-west corner corner array */
  unsigned long ***se;          /* South-east corner corner array */
  unsigned long ***nw;          /* North-west corner corner array */
  unsigned long ***ne;          /* North-east corner corner array */
  int ntr;                      /* Number of tracers */
  int gchck = 0;                /* Checks locations of ghost points */
  int nce1, nfe1;               /* Grid dimension in the e1 direction */
  int nce2, nfe2;               /* Grid dimension in the e2 direction */
  int nz;                       /* Number of vertical layers */
  unsigned long ***flg;         /* Flag for Cartesian grid */
  double *layers;               /* Vertical layers array */
  double *olayers;              /* Original vertical layers array */
  double **bathy;               /* Bathymetry array */
  double d1;                    /* Dummy */

  if (DEBUG("init_m"))
    dlog("init_m", "\nStart preprocessor\n\n");

  /*-----------------------------------------------------------------*/
  /* Read the grid dimensions                                        */
  nce1 = params->nce1;
  nfe1 = nce1 + 1;
  nce2 = params->nce2;
  nfe2 = nce2 + 1;
  olayers = d_alloc_1d(params->nz + 1);
  memcpy(olayers, params->layers, (params->nz + 1) * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Reset flags for sigma calculations if required                  */
  if (params->sigma) {
    if (!(params->stab & (NONE | SUB_STEP_TRACER)))
      params->stab = SUB_STEP;
    params->hmin = 0.0;
    params->thin_merge = 0;
    set_sigma_distrib(params, 0, 0);
    for (n = 0; n < params->nz; n++)
      params->layers[n] *= params->bmax;
  }

  /*-----------------------------------------------------------------*/
  /* Allocate memory                                                 */
  nz = params->nz;
  layers = params->layers;
  flg = (unsigned long ***)l_alloc_3d(nce1 + 1, nce2 + 1, nz);
  bathy = params->topo;
  /* bathy=d_alloc_2d(nce1,nce2); */
  if (params->runmode & DUMP) {
    for (n = 1; n <= params->noutside; n++) {
      i = params->oute1[n];
      j = params->oute2[n];
      bathy[j][i] = NOTVALID;
    }
    for (n = 1; n <= params->nland; n++) {
      i = params->lande1[n];
      j = params->lande2[n];
      bathy[j][i] = LANDCELL;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Make the flags array                                            */
  make_flags(params, flg, bathy, olayers, nce1, nce2, nz);
  u1_flags(flg, nce1, nce2, nz);
  u2_flags(flg, nce1, nce2, nz);
  if (params->sigma)
    sigma_flags(params->nfe1, params->nfe2, params->nz-1, flg);

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */
  geom = (geometry_t *)malloc(sizeof(geometry_t));
  memset(geom, 0, sizeof(geometry_t));
  end_lc = i_alloc_1d(nz + 1);
  num_2d = i_alloc_1d(nz + 1);
  geom->nwindows = params->nwindows;
  if (params->win_size) {
    geom->win_size = d_alloc_1d(geom->nwindows + 1);
    for (n = 1; n <= geom->nwindows; n++) {
      geom->win_size[n] = params->win_size[n - 1];
    }
  }
  if (params->nwn) {
    geom->nwn = params->nwn;
    geom->wnx = params->wnx;  
    geom->wny = params->wny;
  }

  geom->nce1 = nce1;
  geom->nfe1 = nfe1;
  geom->nce2 = nce2;
  geom->nfe2 = nfe2;
  geom->nz = nz;
  geom->sednz = params->sednz;
  geom->layers = d_alloc_1d(nz + 1);
  memcpy(geom->layers, olayers, (params->nz + 1) * sizeof(double));
  d_free_1d(olayers);

  geom->nobc = params->nobc;
  ntr = params->ntr;
  for (k = 0; k <= nz; k++)
    end_lc[k] = 0;

  geom->compatible = params->compatible;

  /*-----------------------------------------------------------------*/
  /* Set up the array of the bottom coordinate, including the edges  */
  /* defined by i=nce1 and j=nce2.                                   */
  kbot = s_alloc_2d(nfe1, nfe2);
  /* Set the bottom coordinate in the interior                       */
  for (j = 0; j < nce2; j++) {
    for (i = 0; i < nce1; i++) {
      if (flg[nz - 1][j][i] & NOTWET) {
        kbot[j][i] = -1;
      } else {
        /* Check bottom is inside model                              */
        if (bathy[j][i] < layers[0])
          hd_quit_and_dump
            ("build_sparse_map() : Bottom lower than lowest model level\n");
        /* Check bottom and top below MAXGRIDZ                       */
        if (bathy[j][i] >= MAXGRIDZ)
          hd_quit_and_dump("build_sparse_map() : Bottom above MAXGRIDZ\n");
        /* Loop vertically to find bottom and top                    */
        kbot[j][i] = -1;
        if (params->sigma)
          kbot[j][i] = 0;
        else {
          for (k = 0; k < nz; k++) {
            if (bathy[j][i] >= layers[k] && bathy[j][i] < layers[k + 1])
              kbot[j][i] = (short)k;
          }
        }
        if (kbot[j][i] < 0)
          hd_quit_and_dump
            ("build_sparse_map() : bottom not found at (%d %d) : %5.1f m\n",
             i, j, bathy[j][i]);
      }
    }
  }

  /*-----------------------------------------------------------------*/
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

  /*-----------------------------------------------------------------*/
  /* Set the flags for the open boundaries                           */
  for (n = 0; n < geom->nobc; n++) {
    for (cc = 0; cc < params->open[n]->npts; cc++) {
      i = params->open[n]->iloc[cc];
      j = params->open[n]->jloc[cc];

      kbot[j][i] =
        set_kbot(i, j, flg[nz - 1], kbot, params->open[n]->type, nce1,
                 nce2);

      for (k = kbot[j][i]; k < nz; k++) {
        if (k == -1)
          hd_quit("Land cell found in boundary%d list at (%d,%d) : bathy = %6.2f\n", n, i, j, bathy[j][i]);
        c2 =
          b_flags(nce1, nce2, nz, flg, i, j, k, &ii, &jj, 1,
                  params->open[n]->type);
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Reset the flag on U1BDRY|R_EDGE and U2BDRY|F_EDGE to wet.       */ 
  /* This is set to wet at this stage since it must be included in   */
  /* the wet domain. Reset to OUTSIDE prior to writing an input file */
  /* in AUTO mode.                                                   */
  for (n = 0; n < geom->nobc; n++) {
    for (cc = 0; cc < params->open[n]->npts; cc++) {
      i = params->open[n]->iloc[cc];
      j = params->open[n]->jloc[cc];
      for (k = kbot[j][i]; k < nz; k++) {
	if (params->open[n]->type == U1BDRY || 
	    params->open[n]->type == U2BDRY) {
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
  /* Reset any (i,j) lists to exclude dry cells                      */
  reset_points(params, flg[nz-1], &params->nzoom, params->zci, params->zcj);

  /*-----------------------------------------------------------------*/
  /* Generate the map from Cartesian space to sparse space,          */
  /* geom->map. First map the cells in the entire grid that may      */
  /* potentially be wet. The sparse grid at this stage consists of   */
  /* wet cells stored in consecutive order, i.e. a 'wet' sparse map  */
  /* is created. The number of wet cells in each layer is stored in  */
  /* end_lc[]. This is used later to insert the ghost cells for each */
  /* layer into locations following the wet cells in each layer. The */
  /* 'wet' map only exists so that the ghost cells can be identified */
  /* at this stage and is reset when the ghost points are inserted.  */
  /* The total number of wet cells is stored in num_wc here.         */
  geom->ewetS = 0;
  geom->map = (unsigned long ***)l_alloc_3d(nfe1, nfe2, nz);
  for (k = nz - 1; k >= 0; k--) {
    for (j = 0; j < nfe2; j++) {
      for (i = 0; i < nfe1; i++) {
        geom->map[k][j][i] = 0;
        if (!(flg[k][j][i] & NOTWET)) {
          geom->map[k][j][i] = c;
          c++;
          num_wc++;
          if (k == nz - 1)
            geom->ewetS++;
          end_lc[k]++;
        }
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* If no wet points are found then issue a warning and quit.  */
  if (num_wc == 0) {
    hd_quit("No wet points located in the domain.\n");
  }
  if (DEBUG("init_m"))
    dlog("init_m", "Wet point counting OK\n");

  /*-----------------------------------------------------------------*/
  /* Find the locations of the ghost cells. These are the first land */
  /* (dry) cells adjacent to a wet cell in the x, y or z direction.  */
  /* The ghost cells in the z direction constitute the sediment */
  /* layer. Note : in certain instances (e.g. outside corners) there */
  /* may exist more than one wet cell associated with the same ghost */
  /* cell -> multiple ghost cells). First count the total number of */
  /* ghost cells in the grid and store the number of ghost cells in */
  /* each layer in num_2d[].  */
  c1 = 0;
  for (k = nz - 1; k >= 0; k--) {
    num_2d[k] = 0;
    for (j = 0; j < nfe2; j++) {
      for (i = 0; i < nfe1; i++) {

        /*-----------------------------------------------------------*/
        /* Count the number of diagonal & inside corner ghost cells */
        /* Note : it is possible for a ghost cell to constitute more */
        /* than one inside corner or diagonal, hence each corner or */
        /* diagonal must be checked individually.  */
        scen = di_ic(geom, i, j, k, &xp1, &xm1, &yp1, &ym1);
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
        /* Count the number of straight edge and outside corner */
        /* ghost cells.  */
        scen = se_oc(geom, i, j, k, &xp1, &xm1, &yp1, &ym1, 1);
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
  geom->nbpt = geom->nbpte1 = geom->nbpte2 = num_gc;
  geom->nbptS = geom->nbpte1S = geom->nbpte2S = num_gc2D = num_2d[nz - 1];

  /*-----------------------------------------------------------------*/
  /* Calculate the number of sediment cells and add to the total */
  /* number of ghost cells. Note: the number of sediment cells is */
  /* equivalent to the number of potentially wet cells in the */
  /* surface layer.  */
  /* First add the number of ghost cells to sedimet cells (this is   */
  /* the same as the number of ghost cells for the surface layer).   */
  num_scg = 0;
  if (!(params->compatible & V1598) && (params->trasc & LAGRANGE))
    num_scg = geom->nbptS;    
  num_sc += num_scg;
  /* Get the sediment cells corresponding to wet cells               */
  for (j = 0; j < nfe2; j++)
    for (i = 0; i < nfe1; i++)
      if (!(flg[nz - 1][j][i] & NOTWET))
        num_sc += 1;
  num_gc += num_sc;

  /*-----------------------------------------------------------------*/
  /* Get the number of OBC ghost cells. There are laux cells for     */
  /* open boundary point so that the advection scheme may be used on */
  /* the boundary.                                                   */
  for (n = 0; n < params->nobc; n++) {
    open_bdrys_t *open = params->open[n];
    if (params->compatible & V1957) continue;
    for (cc = 0; cc < open->npts; cc++) {
      i = open->iloc[cc];
      j = open->jloc[cc];
      if (open->type & U1BDRY) {
	if (i == nce1 - 1 || 
	    (i < nce1 && (bathy[j][i] < layers[0] || bathy[j][i] > params->etamax))) {
	  open->ocodex = R_EDGE;
	  for (k = kbot[j][i]; k < nz; k++) {
	    num_2d[k] += laux - 1;
	    num_gc += laux - 1;
	    num_obc += laux - 1;
	  }
	  num_gc2D += laux - 1;
	  num_obc2D += laux - 1;
	} else {
	  open->ocodex = L_EDGE;
	  for (k = kbot[j][i]; k < nz; k++) {
	    num_2d[k] += laux;
	    num_gc += laux;
	    num_obc += laux;
	  }
	  num_gc2D += laux;
	  num_obc2D += laux;
	}
      }
      if (open->type & U2BDRY) {
	if (j == nce2 - 1 || 
	    (j < nce2 && (bathy[j][i] < layers[0] || bathy[j][i] > params->etamax))) {
	  open->ocodey = F_EDGE;
	  for (k = kbot[j][i]; k < nz; k++) {
	    num_2d[k] += laux - 1;
	    num_gc += laux - 1;
	    num_obc += laux - 1;
	  }
	  num_gc2D += laux - 1;
	  num_obc2D += laux - 1;
	} else {
	  open->ocodey = B_EDGE;
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
  if (DEBUG("init_m"))
    dlog("init_m", "Ghost point counting OK\n");

  /*-----------------------------------------------------------------*/
  /* Re-organize the sparse grid so that the ghost cells in each */
  /* layer are inserted at the end of the wet cells for that layer.  */
  /* This results in the wet cells corresponding to new sparse */
  /* locations in the sparse grid, so the geom->map must be reset.  */
  /* Only maps from Cartesian to sparse space for wet cells are set; */
  /* ghost cells in the sparse grid are set to zero at this stage to */
  /* be assigned values later. The sparse locations of the sediment */
  /* layer lies at the end of the sparse array.  */

  /* Find the end location of wet cells in each layer. This is also */
  /* used later to find the start location of ghost cells in each */
  /* layer to assign the ghost to wet and wet to ghost spatial maps. */
  for (k = nz - 2; k >= 0; k--) {
    c = end_lc[k + 1] + num_2d[k + 1];
    end_lc[k] += c;
  }

  /* Reinitialise the sparse map */
  for (k = nz - 1; k >= 0; k--)
    for (j = 0; j < nfe2; j++)
      for (i = 0; i < nfe1; i++)
        geom->map[k][j][i] = 0;

  /* Reset the Cartestian map to point to new sparse locations */
  c = 1;
  for (k = nz - 1; k >= 0; k--) {
    for (j = 0; j < nfe2; j++) {
      for (i = 0; i < nfe1; i++) {
        if (!(flg[k][j][i] & NOTWET)) {
          geom->map[k][j][i] = c;
          if ((gchck == c) && DEBUG("init_m"))
            dlog("init_m", "%d = wet cell at (%d %d %d)\n", gchck, i, j, k);
          c = c + 1;
        }
      }
    }
    c = end_lc[k] + num_2d[k] + 1;
  }

  if (DEBUG("init_m"))
    dlog("init_m", "Wet point map OK\n");

  /*-----------------------------------------------------------------*/
  /* Save the total number of cells in the sparse grid.              */
  /* Note : these cell locations are for the re-organized array (2D  */
  /* wet + ghost cells followed by 3D wet + ghost) and hence do not  */
  /* reflect the number of cells in the sparse array; e.g. ewet is   */
  /* not just the number of wet cells in the grid but the number of  */
  /* wet + 2D ghost cells. An array subscripted (1:ewet) will not    */
  /* contain all the valid wet cells, but a mixture of wet and       */
  /* ghost. Use geom->a3_t to extract wet cells only.                */
  geom->sgnum = geom->enon = num_wc + num_gc;
  geom->sgnumS = geom->enonS = geom->ewetS + num_gc2D;
  geom->ewet = num_wc + num_gc2D;
  geom->snon = geom->ewet + 1;
  geom->snonS = geom->ewetS + 1;

  /* The sparse grid holds an empty location at 0 (unused) hence the */
  /* size of memory for arrays is geom->sgnum+1.  */
  geom->sgsiz = geom->sgnum + 1;
  geom->sgsizS = geom->sgnumS + 1;

  /*-----------------------------------------------------------------*/
  /* Set the layer maps */
  /* 
     vlr=i_alloc_2d(geom->sgsizS+1,nz);

     for(k=nz-1; k>=0; k--) { for(c=1; c<=geom->sgnumS; c++)vlr[k][c]=0;
     c=1; for(j=0; j<nfe2; j++) { for(i=0; i<nfe1; i++) {
     if(!(flg[k][j][i]&NOTWET)) { vlr[k][c]=geom->map[k][j][i]; c++; } } }
     } */

  /*-----------------------------------------------------------------*/
  /* Allocate memory for the spatial maps and initialize */
  alloc_geom(geom, (MAP_A | GRID_A | MASTER_A));

  /* All the maps are made self-mapping at this stage. Those maps */
  /* which are not re-assigned below remain self pointing.  */
  for (c = 1; c <= geom->sgnum; c++) {
    geom->xp1[c] = c;
    geom->yp1[c] = c;
    geom->zp1[c] = c;
    geom->xm1[c] = c;
    geom->ym1[c] = c;
    geom->zm1[c] = c;
  }

  /* Allocate memory for the corner arrays (used to define multiple */
  /* ghost cells) and initialise. Since the possibility exists for */
  /* wet cells to have multiple ghost points, by splitting the ghost */
  /* cells into those that occupy west, east, north and south edges, */
  /* then a unique ghost cell can be identified and placed in the */
  /* corner arrays, bl, br, tl and tr. Therefore, using the four */
  /* corner arrays a unique map for ghost points from Cartesian */
  /* space to sparse space is achieved.  */
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
  /* Generate the spatial maps for the sparse grid. This includes */
  /* wet cells that map to wet cells (unique map), wet cells that */
  /* map to ghost cells (unique map), ghost cells that map to wet */
  /* cells (multiple maps), and ghost cells that map to ghost cells */
  /* (unique map).  */
  /* First assign the maps of wet cells that map to other wet cells. */
  /* These are defined by cells with non-zero neighbours in the map */
  /* geom->map. The bottom layer maps to the sediment layer and */
  /* vice versa. The sediment layer maps downward to itself (i.e.  */
  /* self pointing).  */
  for (k = nz - 1; k >= 0; k--) {
    for (j = 0; j < nfe2; j++) {
      for (i = 0; i < nfe1; i++) {

        if (k < kbot[j][i])
          continue;
        scen = se_oc(geom, i, j, k, &xp1, &xm1, &yp1, &ym1, 0);

        if (scen != 0) {

          /* c = the sparse location of wet cells at (i,j,k) */
          c = geom->map[k][j][i];

          if (xp1 != 0)
            geom->xp1[c] = geom->map[k][j][i + 1];
          if (xm1 != 0)
            geom->xm1[c] = geom->map[k][j][i - 1];
          if (yp1 != 0)
            geom->yp1[c] = geom->map[k][j + 1][i];
          if (ym1 != 0)
            geom->ym1[c] = geom->map[k][j - 1][i];
          if (k > kbot[j][i] && geom->map[k - 1][j][i] != 0) {
            geom->zm1[c] = geom->map[k - 1][j][i];
          }
          if (k < nz - 1 && geom->map[k + 1][j][i] != 0) {
            geom->zp1[c] = geom->map[k + 1][j][i];
          }
        }
      }
    }
  }
  if (DEBUG("init_m"))
    dlog("init_m", "Wet point spatial maps OK\n");

  /*-----------------------------------------------------------------*/
  /* Allocate memory for the boundary vectors. The vector bpt */
  /* contains the sparse locations of all ghost cells (land boundary */
  /* cells) laterally adjacent to wet cells. The vector bin contains */
  /* the sparse locations of the first wet cell adjacent to the land */
  /* boundary.  */
  geom->bpt = i_alloc_1d(geom->nbpt + 1);
  geom->bin = i_alloc_1d(geom->nbpt + 1);
  geom->dbin = i_alloc_1d(geom->nbpt + 1);
  geom->wgst = i_alloc_1d(geom->sgsiz);

  for (c = 1; c <= geom->nbpt; c++) {
    geom->bpt[c] = 0;
    geom->bin[c] = 0;
  }

  /*-----------------------------------------------------------------*/
  /* Assign the spatial maps of ghost cells that map to wet cells */
  /* and wet cells which map to ghost cells. Also assign the sparse */
  /* locations to bpt and bin. In each layer the ghost cell counter, */
  /* cc, is set to end_lc[] and incremented. Wet cells corresponding */
  /* to diagonals, straight edges, inside/outside corners are */
  /* identified in geom->map on the basis of their neighbour's */
  /* values (zero = ghost; nonzero = wet) and then assigned to cc */
  /* for the map.  */
  lc = 1;
  for (k = nz - 1; k >= 0; k--) {

    /* cc = the sparse location of ghost cells for layer k */
    cc = end_lc[k] + 1;

    for (j = 0; j < nfe2; j++) {
      for (i = 0; i < nfe1; i++) {

        /*-----------------------------------------------------------*/
        /* Map the straight edges and outside corners.  */
        /* Note that ghost points currently are set to zero in */
        /* geom->map, so to get the ghost-wet map a land cell is */
        /* first checked to see if its neighbours are wet. If so */
        /* then set the relevant wet-ghost spatial map equal to the */
        /* ghost cell counter, cc, and the ghost-wet spatial map to */
        /* the sparse coordinate of the wet neighbour cell.  */
        scen = se_oc(geom, i, j, k, &xp1, &xm1, &yp1, &ym1, 1);

        /* Ghost cells on west edges and west outside corners */
        if (scen & W_SE) {
          /* Save the ghost cell location in the corner array */
	  /*sw[k][j][i] = se[k][j][i] = nw[k][j][i] = ne[k][j][i] = cc;*/
	  sw[k][j][i] = nw[k][j][i] = cc;

	  /* Check for maps on east straight edges */
          geom->xm1[xp1] = cc;  /* Wet-ghost spatial map */
          geom->xp1[cc] = xp1;  /* Ghost-wet spatial map */

          /* Check for ghost-wet maps on east outside corners, ie.  */
          /* if the east edge comprises an outside corner there */
          /* will also exist ghost-wet maps in the north/south */
          /* direction depending on the orientation of the corner.  */
          /* The wet-ghost map for these corners will be defined */
          /* with the relevant straight edge map below.  */
          if (yp1 != 0)
            geom->yp1[cc] = yp1;
          if (ym1 != 0)
            geom->ym1[cc] = ym1;

          /* Add the ghost cell to the boundary vectors */
          geom->bpt[lc] = cc;
          geom->bin[lc] = xp1;
          geom->dbin[lc] = 2;
          if ((gchck == cc) && DEBUG("init_m"))
            dlog("init_m", "%d = west edge/outside corner at (%d %d %d)\n",
                 gchck, i, j, k);
          lc++;
          cc++;
        }

        /* Ghost cells on east edges and east outside corners */
        if (scen & E_SE) {
          /* Save the ghost cell location in the corner array */
          /*sw[k][j][i] = se[k][j][i] = nw[k][j][i] = ne[k][j][i] = cc;*/
          se[k][j][i] = ne[k][j][i] = cc;

          /* Check for maps on west straight edges */
          geom->xp1[xm1] = cc;  /* Wet-ghost spatial map */
          geom->xm1[cc] = xm1;  /* Ghost-wet spatial map */
          /* Check for ghost-wet maps on west outside corners */
          if (yp1 != 0)
            geom->yp1[cc] = yp1;
          if (ym1 != 0)
            geom->ym1[cc] = ym1;

          /* Add the ghost cell to the boundary vectors */
          geom->bpt[lc] = cc;
          geom->bin[lc] = xm1;
          geom->dbin[lc] = 6;
          if ((gchck == cc) && DEBUG("init_m"))
            dlog("init_m", "%d = east edge/outside corner at (%d %d %d)\n",
                 gchck, i, j, k);
          lc++;
          cc++;
        }

        /* Ghost cells on south edges and south outside corners */
        if (scen & S_SE) {
          /* Save the ghost cell location in the corner array */
          /*sw[k][j][i] = se[k][j][i] = nw[k][j][i] = ne[k][j][i] = cc;*/
          se[k][j][i] = sw[k][j][i] = cc;

          /* Check for maps on south straight edges */
          geom->ym1[yp1] = cc;  /* Wet-ghost spatial map */
          geom->yp1[cc] = yp1;  /* Ghost-wet spatial map */
          /* Check for ghost-wet maps on south outside corners */
          if (xp1 != 0)
            geom->xp1[cc] = xp1;
          if (xm1 != 0)
            geom->xm1[cc] = xm1;

          /* Add the ghost cell to the boundary vectors */
          geom->bpt[lc] = cc;
          geom->bin[lc] = yp1;
          geom->dbin[lc] = 0;
          if ((gchck == cc) && DEBUG("init_m"))
            dlog("init_m",
                 "%d = south edge/outside corner at (%d %d %d)\n", gchck,
                 i, j, k);
          lc++;
          cc++;
        }

        /* Ghost cells on north edges and north outside corners */
        if (scen & N_SE) {
          /* Save the ghost cell location in the corner array */
          /*sw[k][j][i] = se[k][j][i] = nw[k][j][i] = ne[k][j][i] = cc;*/
          ne[k][j][i] = nw[k][j][i] = cc;

          /* Check for maps on north straight edges */
          geom->yp1[ym1] = cc;  /* Wet-ghost spatial map */
          geom->ym1[cc] = ym1;  /* Ghost-wet spatial map */
          /* Check for ghost-wet maps on north outside corners */
          if (xp1 != 0)
            geom->xp1[cc] = xp1;
          if (xm1 != 0)
            geom->xm1[cc] = xm1;

          /* Add the ghost cell to the boundary vectors */
          geom->bpt[lc] = cc;
          geom->bin[lc] = ym1;
          geom->dbin[lc] = 4;
          if ((gchck == cc) && DEBUG("init_m"))
            dlog("init_m",
                 "%d = north edge/outside corner at (%d %d %d)\n", gchck,
                 i, j, k);
          lc++;
          cc++;
        }

        /*-----------------------------------------------------------*/
        /* Map diagonals and inside corners. This must be performed */
        /* after the straight edge part so that the corner arrays */
        /* are written with unique diagonal/inside corner sparse */
        /* locations. Note: do not use an else if() construct here */
        /* so that any multiple inside corners/diagonals are set.  */
        /* Map the diagonal ghost cells. Here a wet cell is checked */
        /* to see if its neighbours are not wet and if its diagonal */
        /* is wet. If so the ghost counter, cc, is incrmented, bpt */
        /* and the corner array are set to cc, bin is set to the wet */
        /* cell checked. Note: no direct map exist from the wet cell */
        /* checked to its diagonal ghost; these are achieved through */
        /* a wet-ghost and ghost-ghost map.  */

        /* Get the ghost locations of the inside corners. These are */
        /* located by checking if the land cell's neighbours are */
        /* also land cells (non-wet) but the diagonal is wet. As */
        /* with the diagonal ghost cells, a direct map doesn't exist */
        /* from land cell in the corner to its diagonal wet cell; */
        /* these are achieved through the ghost-ghost (defined */
        /* below) and ghost wet (defined above) maps. However, the */
        /* boundary vectors and corner arrays need to be assigned */
        /* here.  */

        scen = di_ic(geom, i, j, k, &xp1, &xm1, &yp1, &ym1);
        if (scen & (SW_D | SW_IC)) {
          sw[k][j][i] = cc;
          geom->bpt[lc] = cc;
          geom->bin[lc] = geom->map[k][j + 1][i + 1];
          geom->dbin[lc] = 1;
          if (scen & SW_D && gchck == cc && DEBUG("init_m"))
            dlog("init_m", "%d = south-west diagonal at (%d %d %d)\n",
                 gchck, i, j, k);
          if (scen & SW_IC && gchck == cc && DEBUG("init_m"))
            dlog("init_m", "%d = south-west inside corner at (%d %d %d)\n",
                 gchck, i, j, k);
          lc++;
          cc++;
        }
        if (scen & (NW_D | NW_IC)) {
          nw[k][j][i] = cc;
          geom->bpt[lc] = cc;
          geom->bin[lc] = geom->map[k][j - 1][i + 1];
          geom->dbin[lc] = 3;
          if (scen & NW_D && gchck == cc && DEBUG("init_m"))
            dlog("init_m", "%d = north-west diagonal at (%d %d %d)\n",
                 gchck, i, j, k);
          if (scen & NW_IC && gchck == cc && DEBUG("init_m"))
            dlog("init_m", "%d = north-west inside corner at (%d %d %d)\n",
                 gchck, i, j, k);
          lc++;
          cc++;
        }
        if (scen & (SE_D | SE_IC)) {
          se[k][j][i] = cc;
          geom->bpt[lc] = cc;
          geom->bin[lc] = geom->map[k][j + 1][i - 1];
          geom->dbin[lc] = 7;
          if (scen & SE_D && gchck == cc && DEBUG("init_m"))
            dlog("init_m", "%d = south-east diagonal at (%d %d %d)\n",
                 gchck, i, j, k);
          if (scen & SE_IC && gchck == cc && DEBUG("init_m"))
            dlog("init_m", "%d = south-east inside corner at (%d %d %d)\n",
                 gchck, i, j, k);
          lc++;
          cc++;
        }
        if (scen & (NE_D | NE_IC)) {
          ne[k][j][i] = cc;
          geom->bpt[lc] = cc;
          geom->bin[lc] = geom->map[k][j - 1][i - 1];
          geom->dbin[lc] = 5;
          if (scen & NE_D && gchck == cc && DEBUG("init_m"))
            dlog("init_m", "%d = north-east diagonal at (%d %d %d)\n",
                 gchck, i, j, k);
          if (scen & NE_IC && gchck == cc && DEBUG("init_m"))
            dlog("init_m", "%d = north-east inside corner at (%d %d %d)\n",
                 gchck, i, j, k);
          lc++;
          cc++;
        }
      }
    }
    /* Save the last sparse location of ghost cells for layer k */
    end_lc[k] = cc;
  }

  /* Set the locations of the corner arrays for the case where the */
  /* cell above the corner cell is wet to the wet cell. This allows */
  /* inside corners that lie below wet cells to map upwards to the */
  /* wet cell. If the cell below is a inside corner also, then the */
  /* corner array will be defined and the downward map will be set */
  /* via the ghost-ghost maps. If the cell below is not a ghost cell */
  /* then the inside corner remains self pointing. Note that */
  /* although the inside corner maps upward to the wet cell, the wet */
  /* cell will map downward to the sediment cell.  */
  for (k = nz - 2; k >= 0; k--) {
    for (j = 0; j < nfe2; j++) {
      for (i = 0; i < nfe1; i++) {
        if (!(flg[nz - 1][j][i] & NOTWET)) {
          if (sw[k][j][i] != 0 && sw[k + 1][j][i] == 0)
            sw[k + 1][j][i] = geom->map[k + 1][j][i];
          if (se[k][j][i] != 0 && se[k + 1][j][i] == 0)
            se[k + 1][j][i] = geom->map[k + 1][j][i];
          if (nw[k][j][i] != 0 && nw[k + 1][j][i] == 0)
            nw[k + 1][j][i] = geom->map[k + 1][j][i];
          if (ne[k][j][i] != 0 && ne[k + 1][j][i] == 0)
            ne[k + 1][j][i] = geom->map[k + 1][j][i];
        }
      }
    }
  }
  if (DEBUG("init_m"))
    dlog("init_m", "Wet point - ghost point spatial maps OK\n");

  /*-----------------------------------------------------------------*/
  /* Assign the spatial maps of ghost cells that map to ghost cells. */
  /* Ghost cells are still set to zero (and not unique) in the the */
  /* Cartesian-sparse map, geom->map, but have been assigned unique */
  /* values in the sparse grid and corner arrays in the above */
  /* mapping procedures.  */

  /* Map the ghost-ghost straight edges, outside corners and maps */
  /* from these locations to inside corners and diagonals.  */
  /* To obtain the ghost-ghost map for straight edge/outside corners */
  /* a wet cell is first checked to see if its neighbours are non- */
  /* wet ghost cells using the maps geom->map and flg. If */
  /* so, the map from this ghost cell to other neighbouring ghost */
  /* cells (if any) is made by approaching the ghost cells via the */
  /* wet cells. Maps from straight edge ghost points to inside */
  /* corner ghost points are assigned also. The corner arrays are */
  /* used to retrieve the ghost cell in the corner and again this */
  /* cell is approached via the wet cell.  */
  for (k = nz - 1; k >= 0; k--) {
    for (j = 0; j < nfe2; j++) {
      for (i = 0; i < nfe1; i++) {
        if (k < kbot[j][i])
          continue;
        if (!(flg[k][j][i] & NOTWET)) {

          /* c = the sparse location of wet cells at (i,j,k) */
          c = geom->map[k][j][i];
          scen = se_oc(geom, i, j, k, &xp1, &xm1, &yp1, &ym1, 0);

          /* Ghost cells on west edges and west outside corners */
          if (scen & W_SE) {
            /* Lateral maps for straight coast/outside corners */
            /* Here (i,j,k) is wet, (i-1,j,k) is dry.  */
            if (yp1 != 0 && flg[k][j + 1][i - 1] & NOTWET)
              geom->yp1[geom->xm1[c]] = geom->xm1[yp1];
            if (ym1 != 0 && flg[k][j - 1][i - 1] & NOTWET)
              geom->ym1[geom->xm1[c]] = geom->xm1[ym1];

            /* Lateral maps from west edge ghost cells to inside */
            /* corners.  */
            /* Here (i,j,k) is wet, (i-1,j,k) is dry and (i,j+1,k) */
            /* & (i-1,j+1,k) are dry for north-west inside corners */
            /* and (i,j-1,k) & (i-1,j-1,k) are dry for south-west */
            /* inside corners.  */
            if (j < nce2 && flg[k][j + 1][i] & NOTWET &&
                flg[k][j + 1][i - 1] & NOTWET)
              geom->yp1[geom->xm1[c]] = nw[k][j + 1][i - 1];
            if (j > 0 && flg[k][j - 1][i] & NOTWET &&
                flg[k][j - 1][i - 1] & NOTWET)
              geom->ym1[geom->xm1[c]] = sw[k][j - 1][i - 1];

            /* Lateral maps from west edge ghost cells to diagonals. */
            /* Here (i,j,k) is wet, (i-1,j,k) is dry and (i,j+1,k) */
            /* & (i-1,j+1,k) are dry/wet for north-west diagonals and */
            /* (i,j-1,k) & (i-1,j-1,k) are dry/wet for south-west */
            /* diagonals.  */
            if (j < nce2 && flg[k][j + 1][i] & NOTWET &&
                !(flg[k][j + 1][i - 1] & NOTWET))
              geom->yp1[geom->xm1[c]] = nw[k][j + 1][i - 1];
            if (j > 0 && flg[k][j - 1][i] & NOTWET &&
                !(flg[k][j - 1][i - 1] & NOTWET))
              geom->ym1[geom->xm1[c]] = sw[k][j - 1][i - 1];

            /* Ghost-ghost vertical maps */
            if (k > kbot[j][i]) {
              geom->zm1[geom->xm1[c]] = geom->xm1[geom->map[k - 1][j][i]];
            }
            if (k < nz - 1)
              geom->zp1[geom->xm1[c]] = geom->xm1[geom->map[k + 1][j][i]];
          }

          /* Ghost cells on east edges */
          if (scen & E_SE) {
            /* Lateral maps for straight coast/outside corners */
            /* Here (i,j,k) is wet, (i+1,j,k) is dry.  */
            if (yp1 != 0 && flg[k][j + 1][i + 1] & NOTWET)
              geom->yp1[geom->xp1[c]] = geom->xp1[yp1];
            if (ym1 != 0 && flg[k][j - 1][i + 1] & NOTWET)
              geom->ym1[geom->xp1[c]] = geom->xp1[ym1];

            /* Lateral maps from east edge ghost cells to inside */
            /* corners.  */
            /* Here (i,j,k) is wet, (i+1,j,k) is dry and (i,j+1,k) */
            /* & (i+1,j+1,k) are dry for north-east inside corners */
            /* and (i,j-1,k) & (i+1,j-1,k) are dry for south */
            /* east inside corners.  */
            if (j < nce2 && flg[k][j + 1][i] & NOTWET &&
                flg[k][j + 1][i + 1] & NOTWET)
              geom->yp1[geom->xp1[c]] = ne[k][j + 1][i + 1];
            if (j > 0 && flg[k][j - 1][i] & NOTWET &&
                flg[k][j - 1][i + 1] & NOTWET)
              geom->ym1[geom->xp1[c]] = se[k][j - 1][i + 1];

            /* Lateral maps from east edge ghost cells to diagonals. */
            /* Here (i,j,k) is wet, (i+1,j,k) is dry and (i,j+1,k) */
            /* & (i+1,j+1,k) are dry/wet for north-east diagonals and */
            /* (i,j-1,k) & (i+1,j-1,k) are dry/wet for south-east */
            /* diagonals.  */
            if (j < nce2 && flg[k][j + 1][i] & NOTWET &&
                !(flg[k][j + 1][i + 1] & NOTWET))
              geom->yp1[geom->xp1[c]] = ne[k][j + 1][i + 1];
            if (j > 0 && flg[k][j - 1][i] & NOTWET &&
                !(flg[k][j - 1][i + 1] & NOTWET))
              geom->ym1[geom->xp1[c]] = se[k][j - 1][i + 1];

            /* Ghost-ghost vertical maps */
            if (k > kbot[j][i])
              geom->zm1[geom->xp1[c]] = geom->xp1[geom->map[k - 1][j][i]];
            if (k < nz - 1)
              geom->zp1[geom->xp1[c]] = geom->xp1[geom->map[k + 1][j][i]];
          }

          /* Ghost cells on north edges */
          if (scen & N_SE) {

            /* Lateral maps for straight coast/outside corners */
            /* Here (i,j,k) is wet, (i,j+1,k) is dry.  */
            if (xp1 != 0 && flg[k][j + 1][i + 1] & NOTWET)
              geom->xp1[geom->yp1[c]] = geom->yp1[xp1];
            if (xm1 != 0 && flg[k][j + 1][i - 1] & NOTWET)
              geom->xm1[geom->yp1[c]] = geom->yp1[xm1];

            /* Lateral maps from north edge ghost cells to inside */
            /* corners.  */
            /* Here (i,j,k) is wet, (i,j+1,k) is dry and (i+1,j,k) */
            /* & (i+1,j+1,k) are dry for north-east inside corners */
            /* and (i-1,j,k) & (i-1,j+1,k) are dry for north-west */
            /* inside corners.  */
            if (i < nce1 && flg[k][j][i + 1] & NOTWET &&
                flg[k][j + 1][i + 1] & NOTWET)
              geom->xp1[geom->yp1[c]] = ne[k][j + 1][i + 1];
            if (i > 0 && flg[k][j][i - 1] & NOTWET &&
                flg[k][j + 1][i - 1] & NOTWET)
              geom->xm1[geom->yp1[c]] = nw[k][j + 1][i - 1];

            /* Lateral maps from north edge ghost cells to diagonals. */
            /* Here (i,j,k) is wet, (i,j+1,k) is dry and (i+1,j,k) */
            /* & (i+1,j+1,k) are dry/wet for north-east diagonals and */
            /* (i-1,j,k) & (i-1,j+1,k) are dry/wet for north-west */
            /* diagonals.  */
            if (i < nce1 && flg[k][j][i + 1] & NOTWET &&
                !(flg[k][j + 1][i + 1] & NOTWET))
              geom->xp1[geom->yp1[c]] = ne[k][j + 1][i + 1];
            if (i > 0 && flg[k][j][i - 1] & NOTWET &&
                !(flg[k][j + 1][i - 1] & NOTWET))
              geom->xm1[geom->yp1[c]] = nw[k][j + 1][i - 1];

            /* Ghost-ghost vertical maps */
            if (k > kbot[j][i])
              geom->zm1[geom->yp1[c]] = geom->yp1[geom->map[k - 1][j][i]];
            if (k < nz - 1)
              geom->zp1[geom->yp1[c]] = geom->yp1[geom->map[k + 1][j][i]];
          }

          /* Ghost cells on south edges */
          if (scen & S_SE) {

            /* Lateral maps for straight coast/outside corners */
            /* Here (i,j,k) is wet, (i,j-1,k) is dry.  */
            if (xp1 != 0 && !(flg[k][j][i + 1] & NOTWET) &&
                flg[k][j - 1][i + 1] & NOTWET)
              geom->xp1[geom->ym1[c]] = geom->ym1[xp1];
            if (xm1 != 0 && !(flg[k][j][i - 1] & NOTWET) &&
                flg[k][j - 1][i - 1] & NOTWET)
              geom->xm1[geom->ym1[c]] = geom->ym1[xm1];

            /* Lateral maps from south edge ghost cells to inside */
            /* corners.  */
            /* Here (i,j,k) is wet, (i,j-1,k) is dry and (i+1,j,k) */
            /* & (i+1,j-1,k) are dry for south-east inside */
            /* corners and (i-1,j,k) & (i-1,j-1,k) are dry for */
            /* south-west inside corners.  */
            if (i < nce1 && flg[k][j][i + 1] & NOTWET &&
                flg[k][j - 1][i + 1] & NOTWET)
              geom->xp1[geom->ym1[c]] = se[k][j - 1][i + 1];
            if (i > 0 && flg[k][j][i - 1] & NOTWET &&
                flg[k][j - 1][i - 1] & NOTWET)
              geom->xm1[geom->ym1[c]] = sw[k][j - 1][i - 1];

            /* Lateral maps from south edge ghost cells to */
            /* diagonals.  */
            /* Here (i,j,k) is wet, (i,j-1,k) is dry and (i+1,j,k) */
            /* & (i+1,j-1,k) are dry/wet for south-east diagonals */
            /* corners and (i-1,j,k) & (i-1,j-1,k) are dry & wet for */
            /* south-west diagonals.  */
            if (i < nce1 && flg[k][j][i + 1] & NOTWET &&
                !(flg[k][j - 1][i + 1] & NOTWET))
              geom->xp1[geom->ym1[c]] = se[k][j - 1][i + 1];
            if (i > 0 && flg[k][j][i - 1] & NOTWET &&
                !(flg[k][j - 1][i - 1] & NOTWET))
              geom->xm1[geom->ym1[c]] = sw[k][j - 1][i - 1];

            /* Ghost-ghost vertical maps */
            if (k > kbot[j][i]) {
              geom->zm1[geom->ym1[c]] = geom->ym1[geom->map[k - 1][j][i]];
            }

            if (k < nz - 1) {
              geom->zp1[geom->ym1[c]] = geom->ym1[geom->map[k + 1][j][i]];
            }
          }
        }
      }
    }
  }
  if (DEBUG("init_m"))
    dlog("init_m", "Ghost point - ghost point spatial maps OK\n");

  /*-----------------------------------------------------------------*/
  /* Assign the ghost-ghost maps from inside corners and diagonals */
  /* to straight edges.  */
  /* To obtain the map for diagonal ghost cells to adjacent ghost */
  /* cells a wet cell is first checked to see if its neighbours are */
  /* dry using flg. Then if the diagonal is wet the diagonal */
  /* ghost cell corresponding to the initial wet cell is retrieved */
  /* via the corner array, and the maps from this ghost cell to */
  /* other neighbouring ghost cells is made by approaching the ghost */
  /* cells via the wet cells.  */
  /* To obtain the map for inside corner ghost cells to adjacent */
  /* ghost cells a wet cell's neighbours are first checked to see if */
  /* they are dry using flg. Then if the diagonal is also dry */
  /* the maps are assigned in the same manner as the diagonal ghost */
  /* cells.  */
  /* Note : diagonal and inside corner ghost cells are generally only */
  /* required to set the boundary values for Semi-Lagrangain type */
  /* advection schemes and do not enter in calculations for most */
  /* finite difference schemes.  */
  for (k = nz - 1; k >= 0; k--) {
    for (j = 0; j < nfe2; j++) {
      for (i = 0; i < nfe1; i++) {

        /* if(k<kbot[j][i])continue; */

        scen = di_ic(geom, i, j, k, &xp1, &xm1, &yp1, &ym1);

	/*-----------------------------------------------------------*/
        /* Maps from diagonal ghost cells to other ghost cells.  */
        /* Lateral maps from south-west diagonals.  */
        /* Here (i+1,j+1,k) is wet, (i+1,j,k) & (i,j+1,k) are dry */
        /* and (i,j,k), the diagonal, is wet.  */
        /* Lateral maps from south-west inside corners.  */
        /* Here (i+1,j+1,k) is wet, (i+1,j,k) & (i,j+1,k) are dry */
        /* and (i,j,k), the inside corner, is dry.  */
        if (i < nce1 && j < nce2 && scen & (SW_D | SW_IC)) {
          c = geom->map[k][j + 1][i + 1];
          geom->xp1[sw[k][j][i]] = geom->ym1[c];
          geom->yp1[sw[k][j][i]] = geom->xm1[c];
          /* Ghost-ghost vertical maps */
          if (k < nz - 1)
            geom->zp1[sw[k][j][i]] = sw[k + 1][j][i];
          if (k > kbot[j + 1][i + 1])
            geom->zm1[sw[k][j][i]] = sw[k - 1][j][i];
        }

        /* Lateral maps from north-west diagonals.  */
        /* Here (i+1,j-1,k) is wet, (i+1,j,k) & (i,j-1,k) are dry */
        /* and (i,j,k), the diagonal, is wet.  */
        /* Lateral maps from south-east inside corners.  */
        /* Here (i+1,j-1,k) is wet, (i+1,j,k) & (i,j-1,k) are dry */
        /* and (i,j,k), the inside corner, is dry.  */
        if (i < nce1 && j > 0 && scen & (NW_D | NW_IC)) {
          c = geom->map[k][j - 1][i + 1];
          geom->xp1[nw[k][j][i]] = geom->yp1[c];
          geom->ym1[nw[k][j][i]] = geom->xm1[c];
          /* Ghost-ghost vertical maps */
          if (k < nz - 1)
            geom->zp1[nw[k][j][i]] = nw[k + 1][j][i];
          if (k > kbot[j - 1][i + 1])
            geom->zm1[nw[k][j][i]] = nw[k - 1][j][i];
        }

        /* Lateral maps from south-east diagonals.  */
        /* Here (i-1,j+1,k) is wet, (i-1,j,k) & (i,j+1,k) are dry */
        /* and (i,j,k), the diagonal, is wet.  */
        /* Lateral maps from north-west inside corners.  */
        /* Here (i-1,j+1,k) is wet, (i-1,j,k) & (i,j+1,k) are dry */
        /* and (i,j,k), the inside corner, is dry.  */
        if (i > 0 && j < nce2 && scen & (SE_D | SE_IC)) {
          c = geom->map[k][j + 1][i - 1];
          geom->xm1[se[k][j][i]] = geom->ym1[c];
          geom->yp1[se[k][j][i]] = geom->xp1[c];
          /* Ghost-ghost vertical maps */
          if (k < nz - 1)
            geom->zp1[se[k][j][i]] = se[k + 1][j][i];
          if (k > kbot[j + 1][i - 1])
            geom->zm1[se[k][j][i]] = se[k - 1][j][i];
        }

        /* Lateral maps from north-east diagonals.  */
        /* Here (i-1,j-1,k) is wet, (i-1,j,k) & (i,j-1,k) are dry */
        /* and (i,j,k), the diagonal, is wet.  */
        /* Lateral maps from north-east inside corners.  */
        /* Here (i-1,j-1,k) is wet, (i-1,j,k) & (i,j-1,k) are dry */
        /* and (i,j,k), the inside corner, is dry.  */
        if (i > 0 && j > 0 && scen & (NE_D | NE_IC)) {
          c = geom->map[k][j - 1][i - 1];
          geom->xm1[ne[k][j][i]] = geom->yp1[c];
          geom->ym1[ne[k][j][i]] = geom->xp1[c];
          /* Ghost-ghost vertical maps */
          if (k < nz - 1)
            geom->zp1[ne[k][j][i]] = ne[k + 1][j][i];
          if (k > kbot[j - 1][i - 1])
            geom->zm1[ne[k][j][i]] = ne[k - 1][j][i];
        }
      }
    }
  }
  if (DEBUG("init_m"))
    dlog("init_m",
         "Inside corner / diagonal connection spatial maps OK\n");

  /* Set the higher order spatial maps */
  /* 
     for(c=1; c<=geom->sgnum; c++) {
     geom->xm12[c]=geom->xm1[geom->xm1[c]];
     geom->ym12[c]=geom->ym1[geom->ym1[c]];
     geom->zm12[c]=geom->zm1[geom->zm1[c]]; } */

  /*-----------------------------------------------------------------*/
  /* Set the sediment - wet maps */
  cc = geom->sgnum - num_sc + 1;
  for (j = 0; j < nfe2; j++)
    for (i = 0; i < nfe1; i++) {
      if (!(flg[nz - 1][j][i] & NOTWET)) {
        k = kbot[j][i];
        c = geom->map[k][j][i];
        geom->zm1[c] = cc;
        geom->zp1[cc] = c;
        /* Horizontal maps for the sediment */
	if (!(params->compatible & V1598) && (params->trasc & LAGRANGE)) { /* Used for Lagrange scheme */
	  geom->xp1[cc] = geom->zm1[geom->xp1[c]];
	  geom->xm1[cc] = geom->zm1[geom->xm1[c]];
	  geom->yp1[cc] = geom->zm1[geom->yp1[c]];
	  geom->ym1[cc] = geom->zm1[geom->ym1[c]];
	} else {  /* Self mapping */
	  geom->xp1[cc] = cc;
	  geom->xm1[cc] = cc;
	  geom->yp1[cc] = cc;
	  geom->ym1[cc] = cc;
	}
        if (gchck == cc && k > 0 && DEBUG("init_m"))
          dlog("init_m", "%d = sediment location at (%d %d %d)\n", gchck,
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
  geom->ngsed = num_scg;
  if (geom->ngsed) {
    geom->gsed_t = i_alloc_1d(geom->ngsed + 1);
    geom->ised_t = i_alloc_1d(geom->ngsed + 1);
    for (i = 1; i <= geom->ngsed; i++) {
      geom->gsed_t[i] = cc;
      cc += 1;
    }
    cc = geom->gsed_t[1];
    /* Vertical maps first */
    for (i = 1; i <= geom->nbptS; i++) {
      c = geom->bpt[i];
      while (c != geom->zm1[c])
	c = geom->zm1[c];
      /* Get the corresponding interior cells to gsed_t */
      j = ANY(c, geom->bpt, geom->nbpt);
      if (j) geom->ised_t[i] = geom->bin[j];
      geom->zm1[c] = cc;
      geom->zp1[cc] = c;
      cc += 1;
    }
    /* Horizontal maps */
    for (i = 1; i <= geom->ngsed; i++) {
      cc = geom->gsed_t[i];    /* Sediment ghost */
      c = geom->zp1[cc];       /* Ghost */

      lc = geom->zm1[geom->xp1[c]];
      geom->xp1[cc] = lc;
      if (geom->xm1[lc] == lc) geom->xm1[lc] = cc;

      lc = geom->zm1[geom->xm1[c]];
      geom->xm1[cc] = lc;
      if (geom->xp1[lc] == lc) geom->xp1[lc] = cc;

      lc = geom->zm1[geom->yp1[c]];
      geom->yp1[cc] = lc;
      if (geom->ym1[lc] == lc) geom->ym1[lc] = cc;

      lc = geom->zm1[geom->ym1[c]];
      geom->ym1[cc] = lc;
      if (geom->yp1[lc] == lc) geom->yp1[lc] = cc;
    }
    num_scg = cc;
  }
  if (DEBUG("init_m"))
    dlog("init_m", "Sediment - wet cell map OK\n");

  /*-----------------------------------------------------------------*/
  /* Assign a value to ghost cells in the Cartesian-sparse map.  */
  /* Since there exists the possibility that multiple ghost cells */
  /* exist, a unique sparse location cannot be assigned to */
  /* geom->map, and these locations in the map are assigned values */
  /* to indicate that the cell exists in sparse space. The corner */
  /* arrays are used to assign these values.  */
  for (j = 1; j < nce2; j++) {
    for (i = 1; i < nce1; i++) {
      for (k = nz - 1; k >= 0; k--) {
        if (k < kbot[j][i])
          continue;
        if (!(flg[k][j][i] & NOTWET)) {
          /* Ghost cells on west straight edges */
          if (flg[k][j][i - 1] & NOTWET)
            geom->map[k][j][i - 1] = sw[k][j][i - 1];
          /* Ghost cells on east straight edges */
          if (flg[k][j][i + 1] & NOTWET)
            geom->map[k][j][i + 1] = se[k][j][i + 1];
          /* Ghost cells on south straight edges */
          if (flg[k][j - 1][i] & NOTWET)
            geom->map[k][j - 1][i] = sw[k][j - 1][i];
          /* Ghost cells on north straight edges */
          if (flg[k][j + 1][i] & NOTWET)
            geom->map[k][j + 1][i] = nw[k][j + 1][i];
          /* Ghost cells on south-west corners */
          if (flg[k][j][i - 1] & NOTWET &&
              flg[k][j - 1][i] & NOTWET && flg[k][j - 1][i - 1] & NOTWET)
            geom->map[k][j - 1][i - 1] = sw[k][j - 1][i - 1];
          /* Ghost cells on north-west corners */
          if (flg[k][j][i - 1] & NOTWET &&
              flg[k][j + 1][i] & NOTWET && flg[k][j + 1][i - 1] & NOTWET)
            geom->map[k][j + 1][i - 1] = nw[k][j + 1][i - 1];
          /* Ghost cells on south-east corners */
          if (flg[k][j][i + 1] & NOTWET &&
              flg[k][j - 1][i] & NOTWET && flg[k][j - 1][i + 1] & NOTWET)
            geom->map[k][j - 1][i + 1] = se[k][j - 1][i + 1];
          /* Ghost cells on north-east corners */
          if (flg[k][j][i + 1] & NOTWET &&
              flg[k][j + 1][i] & NOTWET && flg[k][j + 1][i + 1] & NOTWET)
            geom->map[k][j + 1][i + 1] = ne[k][j + 1][i + 1];
        }
      }
    }
  }
  /* Count the number of multiple ghost cells */
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

  /* Get the multiple ghost cell mappings */
  geom->mgc = i_alloc_1d(geom->sgsiz);
  memset(geom->mgc, 0, geom->sgsiz * sizeof(int));
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
	c = lc = geom->map[k][j][i];
        for (cc = 1; cc <= 4; cc++) {
	  if (a[cc] && a[cc] != lc && !geom->mgc[c]) {
	    geom->mgc[c] = a[cc];
	    c = a[cc];
	  }
	}
      }
    }
  }

  l_free_3d((long ***)sw);
  l_free_3d((long ***)se);
  l_free_3d((long ***)nw);
  l_free_3d((long ***)ne);
  i_free_1d(num_2d);
  if (DEBUG("init_m"))
    dlog("init_m", "Multiple ghost points set OK\n");

  /*-----------------------------------------------------------------*/
  /* Set the flags and maps for explicitly defined maps.             */
  /* e1 direction                                                    */
  geom->sm_e1 = set_explicit_maps(geom, params->nemx, params->emisx, 
				  params->emjsx, params->emidx, 
				  params->emjdx, params->ktx, params->kbx, 
				  flg, kbot, geom->xp1, geom->xm1, "E1",
				  &params->exmapf);
  
  /* e2 direction                                                    */
  geom->sm_e2 = set_explicit_maps(geom, params->nemy, params->emisy, 
				  params->emjsy, params->emidy, 
				  params->emjdy, params->kty, params->kby,
				  flg, kbot, geom->yp1, geom->ym1, "E2",
				  &params->exmapf);

  /*-----------------------------------------------------------------*/
  /* Save the locations of sparse cells on open boundaries to the */
  /* boundary arrays.  */
  /* The open boundary location vectors are: */
  /* obc_t = cell centered sparse location of the open boundary */
  /* oi1_t = one cell into the interior from obc_t */
  /* oi2_t = two cells into the interior from obc_t */
  /* obc_e1 = u1 face centered sparse location of the open boundary */
  /* oi1_e1 = one cell into the interior from obc_e1 */
  /* oi2_e1 = two cells into the interior from obc_e1 */
  /* obc_e2 = u2 face centered sparse location of the open boundary */
  /* oi1_e2 = one cell into the interior from obc_e2 */
  /* oi2_e2 = two cells into the interior from obc_e2 */

  /* Allocate memory for the open boundary data structure */
  geom->open =
    (open_bdrys_t **)malloc(sizeof(open_bdrys_t *) * geom->nobc);

  /* Allocate memory for the boundary structures */
  for (n = 0; n < geom->nobc; n++) {

    geom->open[n] = OBC_alloc();

    /* 
       geom->open[n]->bdata=(bdry_details_t
       *)malloc(sizeof(bdry_details_t));
       geom->open[n]->bdata_t=(bdry_details_t
       *)malloc(sizeof(bdry_details_t)*ntr); */
    geom->open[n]->ntr = ntr;

    copy_OBC_conds(params->open[n], geom->open[n], n, params->trinfo_3d);

    /* Set the flags for this open boundary */
    /* 
       for(cc=0; cc<params->open[n]->npts; cc++) {
       i=params->open[n]->iloc[cc]; j=params->open[n]->jloc[cc];
       kb=set_kbot(i,j,flg[nz-1],kbot,geom->open[n]->type,nce1,nce2);
       for(k=kb; k<nz; k++) { if(k==-1)hd_quit("Land cell found in
       boundary %d list at %d %d\n",n,i,j);
       c2=b_flags(nce1,nce2,nz,flg,i,j,k,&ii,&jj,1,geom->open[n]->type); }
       } */
  }

  /* Set up and save sparse locations for the OBC vectors */
  for (n = 0; n < geom->nobc; n++) {

    set_OBC_cells(geom, geom->open[n], params->open[n], n, nce1, nce2, nz, flg);

    /* Set the interior face stagger for normal component */
    /* Inner staggerd cells adjacent to land are omitted */
    if (geom->open[n]->stagger & INFACE) {
      c1 = c2 = c3 = 1;
      /* u1 velocity */
      for (cc = 1; cc <= geom->open[n]->to3_e1; cc++) {
	b1 = geom->open[n]->obc_e1[cc];
	b2 = geom->open[n]->oi1_e1[cc];
	b3 = geom->open[n]->oi2_e1[cc];
	if (cc <= geom->open[n]->no3_e1) {
	  if (b2 != b3 && b3 != geom->xm1[b3]) {
	  /*if (b2 != b3) {*/
	    geom->open[n]->obc_e1[c1] = b2;
	    geom->open[n]->oi1_e1[c1] = b3;
	    if (geom->open[n]->ocodex & L_EDGE)
	      geom->open[n]->oi2_e1[c1] = geom->xp1[b3];
	    if (geom->open[n]->ocodex & R_EDGE)
	      geom->open[n]->oi2_e1[c1] = geom->xm1[b3];
	    c1 += 1;
	    if (geom->zp1[b1] == b1) c2 += 1;
	    c3 += 1;
	  }
	} else {
	  geom->open[n]->obc_e1[c3] = b1;
	  geom->open[n]->oi1_e1[c3] = b2;
	  geom->open[n]->oi2_e1[c3] = b3;
	  c3 += 1;
	} 
      }
      c1 = geom->open[n]->no3_e1 - c1 + 1;
      c2 = geom->open[n]->no2_e1 - c2 + 1;
      if (c1) {
	geom->open[n]->no3_e1 -= c1;
	geom->open[n]->no2_e1 -= c2;
	geom->open[n]->to2_e1 -= c1;
	geom->open[n]->to3_e1 -= c1;
      }
      /* u2 velocity */
      c1 = c2 = c3 = 1;
      for (cc = 1; cc <= geom->open[n]->to3_e2; cc++) {
	b1 = geom->open[n]->obc_e2[cc];
	b2 = geom->open[n]->oi1_e2[cc];
	b3 = geom->open[n]->oi2_e2[cc];
	if (cc <= geom->open[n]->no3_e2) {
	  if (b2 != b3&& b3 != geom->ym1[b3]) {
	    /*if (b2 != b3) {*/
	    geom->open[n]->obc_e2[c1] = b2;
	    geom->open[n]->oi1_e2[c1] = b3;
	    if (geom->open[n]->ocodey & B_EDGE)
	      geom->open[n]->oi2_e2[c1] = geom->yp1[b3];
	    if (geom->open[n]->ocodey & F_EDGE)
	      geom->open[n]->oi2_e2[c1] = geom->ym1[b3];
	    c1 += 1;
	    if (geom->zp1[b1] == b1) c2 += 1;
	    c3 += 1;
	  }
	} else {
	  geom->open[n]->obc_e2[c3] = b1;
	  geom->open[n]->oi1_e2[c3] = b2;
	  geom->open[n]->oi2_e2[c3] = b3;
	  c3 += 1;
	}
      }
      c1 = geom->open[n]->no3_e2 - c1 + 1;
      c2 = geom->open[n]->no2_e2 - c2 + 1;
      if (c1) {
	geom->open[n]->no3_e2 -= c1;
	geom->open[n]->no2_e2 -= c2;
	geom->open[n]->to2_e2 -= c1;
	geom->open[n]->to3_e2 -= c1;
      }
    }

    /* Point the interior cell maps on the boundaries to the correct */
    /* spatial map.  */
    if (geom->open[n]->ocodex & (L_EDGE|LO_EDGE)) {
      geom->open[n]->nmap = geom->xp1;
      geom->open[n]->omap = geom->xm1;
      geom->open[n]->tmpp = geom->yp1;
      geom->open[n]->tmpm = geom->ym1;
    }
    if (geom->open[n]->ocodex & R_EDGE) {
      geom->open[n]->nmap = geom->xm1;
      geom->open[n]->omap = geom->xp1;
      geom->open[n]->tmpp = geom->yp1;
      geom->open[n]->tmpm = geom->ym1;
    }
    if (geom->open[n]->ocodey & (B_EDGE|BO_EDGE)) {
      geom->open[n]->nmap = geom->yp1;
      geom->open[n]->omap = geom->ym1;
      geom->open[n]->tmpp = geom->xp1;
      geom->open[n]->tmpm = geom->xm1;
    }
    if (geom->open[n]->ocodey & F_EDGE) {
      geom->open[n]->nmap = geom->ym1;
      geom->open[n]->omap = geom->yp1;
      geom->open[n]->tmpp = geom->xp1;
      geom->open[n]->tmpm = geom->xm1;
    }
  }

  /* Set maps to be self-pointing over R_EDGE and F_EDGE OBCs */
  set_map_bdry(geom);

  /* Get the cyclic boundary locations. */
  for (n = 0; n < geom->nobc; n++) {
    char buf[MAXSTRLEN];

    /* Set any CYCLED maps. */
    make_cycled(geom, geom->open[n], flg);

    /* Tracers */
    for (cc = 1; cc <= geom->open[n]->no3_t; cc++) {
      /* i[j]loc are defined over no2_t */
      k = (cc-1) % geom->open[n]->no2_t;
      sprintf(buf, "tracer OBC%d, cell %d (%d %d)", n, cc, 
	      params->open[n]->iloc[k], params->open[n]->jloc[k]);
      c = isvalidc(geom->open[n]->obc_t[cc], geom->sgnum, buf);
      if (!(geom->open[n]->bcond_ele & CYCLED) &&
	  !(ANY0(CYCLED, geom->open[n]->bcond_tra, geom->open[n]->ntr)))
	geom->open[n]->cyc_t[cc] = cyclic_m2(geom, geom->open[n]->ocodex,
					     geom->open[n]->ocodey,
					     geom->open[n]->nmap, c);
    }

    /* u1 velocity */
    /* Normal component */
    for (cc = 1; cc <= geom->open[n]->no3_e1; cc++) {
      c = isvalidc(geom->open[n]->obc_e1[cc], geom->sgnum, 
		   "e1 normal velocity OBCs");
      if (!(geom->open[n]->bcond_nor & CYCLED))
	geom->open[n]->cyc_e1[cc] = cyclic_m1(geom, geom->open[n]->ocodex,
					      geom->open[n]->ocodey,
					      geom->open[n]->nmap, c);
    }

    /* Tangential component */
    for (cc = geom->open[n]->no3_e1 + 1; cc <= geom->open[n]->to3_e1; cc++) {
      c = isvalidc(geom->open[n]->obc_e1[cc], geom->sgnum, 
		   "e1 tangential velocity OBCs");
      if (!(geom->open[n]->bcond_tan & CYCLED))
	geom->open[n]->cyc_e1[cc] = cyclic_m2(geom, geom->open[n]->ocodex,
					      geom->open[n]->ocodey,
					      geom->open[n]->nmap, c);
    }

    /* u2 velocity */
    /* Normal component */
    for (cc = 1; cc <= geom->open[n]->no3_e2; cc++) {
      c = isvalidc(geom->open[n]->obc_e2[cc], geom->sgnum, 
		   "e2 normal velocity OBCs");
      if (!(geom->open[n]->bcond_nor & CYCLED))
	geom->open[n]->cyc_e2[cc] = cyclic_m1(geom, geom->open[n]->ocodex,
					      geom->open[n]->ocodey,
					      geom->open[n]->nmap, c);
    }

    /* Tangential component */
    for (cc = geom->open[n]->no3_e2 + 1; cc <= geom->open[n]->to3_e2; cc++) {
      c = isvalidc(geom->open[n]->obc_e2[cc], geom->sgnum, 
		   "e2 tangential velocity OBCs");
      if (!(geom->open[n]->bcond_tan & CYCLED))
	geom->open[n]->cyc_e2[cc] = cyclic_m2(geom, geom->open[n]->ocodex,
					      geom->open[n]->ocodey,
					      geom->open[n]->nmap, c);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the OBC ghost cell maps.                                    */
  for (n = 0; n < params->nobc; n++) {
    open_bdrys_t *open = params->open[n];
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
	    c = c2 = geom->map[k][j][i];
	    c1 = end_lc[k];
	    for (ii = 0; ii < nlaux; ii++) {
	      geom->open[n]->omap[c] = c1;
	      geom->open[n]->nmap[c1] = c;
	      geom->wgst[c1] = c2;
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
      if (geom->open[n]->ocodex & R_EDGE)
	geom->open[n]->ogc_t = geom->open[n]->obc_e1;
      else if (geom->open[n]->ocodey & F_EDGE)
	geom->open[n]->ogc_t = geom->open[n]->obc_e2;
      else {
	geom->open[n]->ogc_t = i_alloc_1d(geom->open[n]->no3_t + 1);
	for (cc = 1; cc <= geom->open[n]->no3_t; cc++) {
	  c = geom->open[n]->obc_t[cc];
	  geom->open[n]->ogc_t[cc] = geom->open[n]->omap[c];
	}
      }
    }
  }
  i_free_1d(end_lc);

  if (DEBUG("init_m"))
    dlog("init_m", "Open boundary maps OK\n");

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
  /* n3_t = number of cell centered wet + ghost cells                */
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
  geom->a3_t = 0;
  geom->a3_e1 = 0;
  geom->x3_e1 = 0;
  geom->a3_e2 = 0;
  geom->x3_e2 = 0;
  geom->a2_t = 0;
  geom->a2_e1 = 0;
  geom->x2_e1 = 0;
  geom->a2_e2 = 0;
  geom->x2_e2 = 0;

  /* Count the number of cell centered locations that are wet and    */
  /* not adjacent to a u1 or u2 boundary. Note: this is not simply   */
  /* the quantity num_wc-nu1bdry-nu2bdry since some cells are both   */
  /* u1 and u2 boundaries (e.g. corner cells).                       */

  /* Note that OUTSIDE status has been removed on U1BDRY|R_EDGE and  */
  /* U2BDRY|F_EDGE at this stage, and is re-established later.       */

  /* pre-V1246 included a2_/a3_ cells in the b2_/b3_ vector.         */
  /* post-V1246 split a2_/a3_ cells and b2_/b3_ cells.               */

  geom->v3_t = geom->v2_t = 0;
  geom->b3_t = geom->b2_t = 0;
  for (k = nz - 1; k >= 0; k--) {
    for (j = 0; j < nfe2; j++) {
      for (i = 0; i < nfe1; i++) {
        int fi = i + 1;
        int fj = j + 1;
        int bi = (i == 0) ? i : i - 1;
        int bj = (j == 0) ? j : j - 1;

        /* Count the cells that are wet or boundary cells, including */
	/* OBC cells at R_EDGE or F_EDGEs.                           */
	/* Note: this should never be invoked if OUTSIDE status is   */
	/* removed - retain for compatibility.                       */
        if (flg[k][j][i] & OUTSIDE &&
            (flg[k][j][i] & U1BDRY || flg[k][j][i] & U2BDRY)) {
          geom->a3_t++;
          if (k == nz - 1)
            geom->a2_t++;
        }

        /* Count the cells centers which are wet                     */
        if (!(flg[k][j][i] & (NOTWET))) {

          /* Count wet cell centers on the whole nfe1 x nfe2 grid    */
          /* (i.e. sparse cells which may have a valid tracer, e1 or */
          /* e2 velocity value).                                     */
          geom->a3_t++;
          if (k == nz - 1)
            geom->a2_t++;

          if (i < nce1 && j < nce2) {
            /* Count wet cell centers on the whole nce1 x nce2 grid  */
            /* (i.e. valid tracer / elevation cells).                */
	    if (params->compatible & V1246) {
		geom->b3_t++;
		if (k == nz - 1)
		  geom->b2_t++;
	    } else {
	      if (!(flg[k][j][i]&(R_EDGE|F_EDGE))) {
		geom->b3_t++;
		if (k == nz - 1) {
		  geom->b2_t++;		  
		}
	      }
	    }
            /* Count non-boundary wet cell centers                   */
            if (!(flg[k][j][i] & (U1BDRY | U2BDRY)) &&
                !(flg[k][fj][i] & U2BDRY) && !(flg[k][j][fi] & U1BDRY)) {
              geom->v3_t++;
              if (k == nz - 1)
                geom->v2_t++;
            }
          }
        }
      }
    }
  }

  /* Count the number of wet cells for e1 and e2 velocity            */
  geom->v3_e1 = geom->v2_e1 = geom->b3_e1 = geom->b2_e1 = 0;
  geom->v3_e2 = geom->v2_e2 = geom->b3_e2 = geom->b2_e2 = 0;
  for (k = nz - 1; k >= 0; k--) {
    for (j = 0; j < nfe2; j++) {
      for (i = 0; i < nfe1; i++) {
        if (j < nce2 && !(flg[k][j][i] & (U1SOLID | U1OUTSIDE))) {
	  if (!(flg[k][j][i] & U1GEN))
	    geom->b3_e1 += 1;
          if (k == nz - 1)
            geom->b2_e1++;
          if (!(flg[k][j][i] & (U1BDRY))) {
	    if (!(flg[k][j][i] & U1GEN))
	      geom->v3_e1 += 1;
            if (k == nz - 1)
              geom->v2_e1++;
          }
        }
        if (i < nce1 && !(flg[k][j][i] & (U2SOLID | U2OUTSIDE))) {
	  if (!(flg[k][j][i] & U2GEN))
	    geom->b3_e2 += 1;
          if (k == nz - 1)
            geom->b2_e2++;
          if (!(flg[k][j][i] & (U2BDRY))) {
	    if (!(flg[k][j][i] & U2GEN))
	      geom->v3_e2 += 1;
            if (k == nz - 1)
              geom->v2_e2++;
          }
        }
      }
    }
  }

  /* Allocate the cells to process vectors                           */
  geom->w3_t = i_alloc_1d(geom->a3_t + geom->nbpt + 1);
  geom->w3_e1 = i_alloc_1d(geom->b3_e1 + geom->nbpt + 1);
  geom->w3_e2 = i_alloc_1d(geom->b3_e2 + geom->nbpt + 1);
  geom->wsa = i_alloc_1d(geom->a3_t + 1);

  geom->w2_t = i_alloc_1d(geom->a2_t + geom->nbptS + 1);
  geom->w2_e1 = i_alloc_1d(geom->b2_e1 + geom->nbptS + 1);
  geom->w2_e2 = i_alloc_1d(geom->b2_e2 + geom->nbptS + 1);

  geom->bot_t = i_alloc_1d(geom->b2_t + 1);
  geom->bot_e1 = i_alloc_1d(geom->b2_e1 + 1);
  geom->bot_e2 = i_alloc_1d(geom->b2_e2 + 1);
  geom->sur_t = i_alloc_1d(geom->b2_t + 1);
  geom->sur_e1 = i_alloc_1d(geom->sgsizS);
  geom->sur_e2 = i_alloc_1d(geom->sgsizS);
  params->ns2 = geom->a2_t;
  params->ns3 = geom->a3_t;

  c = 0;
  for (j = 0; j < nfe2; j++)
    for (i = 0; i < nfe1; i++)
      if (geom->map[nz - 1][j][i] > (unsigned long)c)
        c = geom->map[nz - 1][j][i];
  geom->c2cc = i_alloc_1d(c + 1);

  /* Fill the cell centered vector with non-boundary cell locations  */
  geom->a3_t = geom->b3_t + 1;
  geom->a2_t = geom->b2_t + 1;
  geom->b3_t = geom->v3_t + 1;
  geom->b2_t = geom->v2_t + 1;
  geom->v3_t = 1;
  geom->v2_t = 1;
  c1 = 1;
  for (k = nz - 1; k >= 0; k--) {
    for (j = 0; j < nfe2; j++) {
      for (i = 0; i < nfe1; i++) {
        int fi = i + 1;
        int fj = j + 1;
        c = geom->map[k][j][i];

        /* Map the sparse grid, 2D cells first followed by 3D cells  */
        if ((!(flg[k][j][i] & (NOTWET))) ||
            (flg[k][j][i] & OUTSIDE && (flg[k][j][i] & U1BDRY ||
                                        flg[k][j][i] & U2BDRY))) {
          geom->wsa[c1] = c;
          c1++;
        }

        /* Map the e1 or e2 boundary cells whose cell centers are    */
        /* OUTSIDE cells.                                            */
	/* Old code
	if(flg[k][j][i] & OUTSIDE && 
	   (flg[k][j][i]&U1BDRY || flg[k][j][i]&U2BDRY)) {
				      
	  geom->w3_t[geom->a3_t] = c;
	  geom->a3_t += 1; 
	  if(k==nz-1) { 
	    geom->w2_t[geom->a2_t] = c;
	    geom->a2_t++; 
	  } 
	}
	*/

        /* Map the e1 or e2 boundary cells whose cell centers lie on */
        /* the nce1 or nce2 limits of the grid and those in the      */
	/* interior that correspond to open boundaries on e1 or e2   */
	/* faces.                                                    */
        if (!(flg[k][j][i] & (NOTWET))) {
	  if (i == nce1 || j == nce2) {
	    geom->w3_t[geom->a3_t] = c;
	    geom->a3_t += 1;
	    if (k == nz - 1) {
	      geom->w2_t[geom->a2_t] = c;
	      geom->a2_t++;
	    }
	  } else if (!(params->compatible & V1246) &&
		     ((flg[k][j][i] & (U1BDRY|U2BDRY)) &&
		      (flg[k][j][i] & (R_EDGE|F_EDGE)))) {
	    geom->w3_t[geom->a3_t] = c;
	    geom->a3_t += 1; 
	    if(k==nz-1) { 
	      geom->w2_t[geom->a2_t] = c;
	      geom->a2_t++;
	    }
	  }
	}

        /* Map the cell centered wet cells                           */
        if (!(flg[k][j][i] & (NOTWET)) && i < nce1 && j < nce2) {
          /* Cell centers that do not lie on an open boundary        */
          if (!(flg[k][j][i] & (U1BDRY | U2BDRY)) &&
              !(flg[k][fj][i] & U2BDRY) && !(flg[k][j][fi] & U1BDRY)) {
            geom->w3_t[geom->v3_t] = c;
            geom->v3_t += 1;
            if (k == nz - 1) {
              geom->w2_t[geom->v2_t] = c;
              c2 = geom->map[kbot[j][i]][j][i];
              geom->bot_t[geom->v2_t] = c2;
              geom->v2_t++;
            }
          }
          /* Cell centers on eta/tracer open boundaries              */
          else if (params->compatible & V1246 || 
		   (flg[k][j][i]&(U1BDRY|U2BDRY) &&
		    !(flg[k][j][i]&(R_EDGE|F_EDGE))) ||
		   (flg[k][j][fi]&U1BDRY && flg[k][j][fi]&R_EDGE) ||
		   (flg[k][fj][i]&U2BDRY && flg[k][fj][i]&F_EDGE)) {
            geom->w3_t[geom->b3_t] = c;
            geom->b3_t += 1;
            if (k == nz - 1) {
              geom->w2_t[geom->b2_t] = c;
              c2 = geom->map[kbot[j][i]][j][i];
              geom->bot_t[geom->b2_t] = c2;
              geom->b2_t++;
            }
          }
        }
      }
    }
  }
  geom->v3_t -= 1;
  geom->v2_t -= 1;
  geom->b3_t -= 1;
  geom->b2_t -= 1;
  geom->a3_t -= 1;
  geom->a2_t -= 1;
  geom->ns2 = geom->a2_t;
  geom->ns3 = geom->a3_t;

  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    geom->sur_t[cc] = c;
    geom->c2cc[c] = cc;
  }

  /* e1 face centered wet cells. This must be done in the order of x */
  /* direction first, followed by y then z.  */
  geom->b3_e1 = geom->v3_e1 + 1;
  geom->b2_e1 = geom->v2_e1 + 1;
  geom->v3_e1 = 1;
  geom->v2_e1 = 1;
  for (k = nz - 1; k >= 0; k--) {
    for (j = 0; j < nce2; j++) {
      for (i = 0; i < nfe1; i++) {
        c = geom->map[k][j][i];
        if (!(flg[k][j][i] & (U1SOLID | U1OUTSIDE))) {
          if (!(flg[k][j][i] & U1BDRY)) {
	    if (!(flg[k][j][i] & U1GEN)) {
	      geom->w3_e1[geom->v3_e1] = c;
	      geom->v3_e1++;
	    }
            if (k == nz - 1) {
              geom->w2_e1[geom->v2_e1] = c;
              geom->v2_e1++;
            }
          } else {
	    if (!(flg[k][j][i] & U1GEN)) {
	      geom->w3_e1[geom->b3_e1] = c;
	      geom->b3_e1++;
	    }
            if (k == nz - 1) {
              geom->w2_e1[geom->b2_e1] = c;
              geom->b2_e1++;
            }
          }
        }
      }
    }
  }
  geom->v3_e1 -= 1;
  geom->v2_e1 -= 1;
  geom->b3_e1 -= 1;
  geom->b2_e1 -= 1;

  /* e2 face centered wet cells.  */
  geom->b3_e2 = geom->v3_e2 + 1;
  geom->b2_e2 = geom->v2_e2 + 1;
  geom->v3_e2 = 1;
  geom->v2_e2 = 1;
  for (k = nz - 1; k >= 0; k--) {
    for (j = 0; j < nfe2; j++) {
      for (i = 0; i < nce1; i++) {
        c = geom->map[k][j][i];
        if (!(flg[k][j][i] & (U2SOLID | U2OUTSIDE))) {
          if (!(flg[k][j][i] & U2BDRY)) {
	    if (!(flg[k][j][i] & U2GEN)) {
	      geom->w3_e2[geom->v3_e2] = c;
	      geom->v3_e2++;
	    }
            if (k == nz - 1) {
              geom->w2_e2[geom->v2_e2] = c;
              geom->v2_e2++;
            }
          } else {
	    if (!(flg[k][j][i] & U2GEN)) {
	      geom->w3_e2[geom->b3_e2] = c;
	      geom->b3_e2++;
	    }
            if (k == nz - 1) {
              geom->w2_e2[geom->b2_e2] = c;
              geom->b2_e2++;
            }
          }
        }
      }
    }
  }
  geom->v3_e2 -= 1;
  geom->v2_e2 -= 1;
  geom->b3_e2 -= 1;
  geom->b2_e2 -= 1;

  /* Fill the bottom vectors with the sparse coordinate of the */
  /* bottom. This must be done in the order of z direction first, */
  /* followed by x then y (and hence cannot be incorporated into the */
  /* wet cell mask loop above).  */
  b1 = b2 = 1;
  c1 = geom->v2_e1 + 1;
  c2 = geom->v2_e2 + 1;
  for (j = 0; j < nfe2; j++) {
    for (i = 0; i < nfe1; i++) {
      for (k = nz - 1; k >= 0; k--) {
        c = geom->map[k][j][i];
        if (!(flg[k][j][i] & NOTWET)) {
          if (j < nce2 && !(flg[k][j][i] & (U1SOLID | U1OUTSIDE))) {
            if (k == 0) {
              if (!(flg[k][j][i] & U1BDRY)) {
                geom->bot_e1[b1] = c;
                b1++;
              } else {
                geom->bot_e1[c1] = c;
                c1++;
              }
            } else {
              if (flg[k - 1][j][i] & (U1SOLID | U1OUTSIDE)) {
                if (!(flg[k][j][i] & U1BDRY)) {
                  geom->bot_e1[b1] = c;
                  b1++;
                } else {
                  geom->bot_e1[c1] = c;
                  c1++;
                }
              }
            }
          }
          if (i < nce1 && !(flg[k][j][i] & (U2SOLID | U2OUTSIDE))) {
            if (k == 0) {
              if (!(flg[k][j][i] & U2BDRY)) {
                geom->bot_e2[b2] = c;
                b2++;
              } else {
                geom->bot_e2[c2] = c;
                c2++;
              }
            } else {
              if (flg[k - 1][j][i] & (U2SOLID | U2OUTSIDE)) {
                if (!(flg[k][j][i] & U2BDRY)) {
                  geom->bot_e2[b2] = c;
                  b2++;
                } else {
                  geom->bot_e2[c2] = c;
                  c2++;
                }
              }
            }
          }
        }
      }
    }
  }
  if (DEBUG("init_m"))
    dlog("init_m", "Wet cell masks OK\n");

  /*-----------------------------------------------------------------*/
  /* Set the sparse to Cartesian maps */
  for (c = 1; c <= geom->sgnum; c++) {
    geom->s2i[c] = NOTVALID;
    geom->s2j[c] = NOTVALID;
    geom->s2k[c] = NOTVALID;
  }
  for (k = nz - 1; k >= 0; k--)
    for (j = 0; j < nfe2; j++)
      for (i = 0; i < nfe1; i++) {
        if (!(flg[k][j][i] & NOTWET)) {
          c = geom->map[k][j][i];
          geom->s2i[c] = (short)i;
          geom->s2j[c] = (short)j;
          geom->s2k[c] = (short)k;
        }
      }

  /* Sparse - Cartesian maps for the sediment */
  c = geom->sgnum - num_sc + 1;
  for (j = 0; j < nfe2; j++)
    for (i = 0; i < nfe1; i++) {
      if (!(flg[nz - 1][j][i] & NOTWET)) {
        k = kbot[j][i] - 1;
        geom->s2i[c] = (short)i;
        geom->s2j[c] = (short)j;
        if (k >= 0) {
          geom->map[k][j][i] = c;
          geom->s2k[c] = (short)k;
        }
	if (k == -1) {
	  geom->s2i[c] = NOTVALID;
	  geom->s2j[c] = NOTVALID;
	}
        c++;
      }
    }
  if (DEBUG("init_m"))
    dlog("init_m", "Sparse - Cartesian maps OK\n");

  /*-----------------------------------------------------------------*/
  /* Get the lateral boundary masks.  */
  /* nbpt = number of lateral boundary cells (ghost cells) */
  /* nbptS = number of ghost cells in the 2D grid */
  /* bpt = cell centered lateral boundary mask */
  /* bin = one point into the interior from bpt */
  /* bpte1 = u1 face centered lateral boundary mask */
  /* bine1 = one point into the interior from bpte1 */
  /* bpte2 = u2 face centered lateral boundary mask */
  /* bine2 = one point into the interior from bpte2 */
  geom->bin2 = i_alloc_1d(geom->nbpt + 1);
  geom->bpte1 = i_alloc_1d(geom->nbpt + 1);
  geom->bpte2 = i_alloc_1d(geom->nbpt + 1);
  geom->bine1 = i_alloc_1d(geom->nbpt + 1);
  geom->bine2 = i_alloc_1d(geom->nbpt + 1);
  geom->bpte1S = i_alloc_1d(geom->nbptS + 1);
  geom->bpte2S = i_alloc_1d(geom->nbptS + 1);
  geom->bine1S = i_alloc_1d(geom->nbptS + 1);
  geom->bine2S = i_alloc_1d(geom->nbptS + 1);
  geom->brsm = i_alloc_1d(geom->sgsizS);

#if defined(HAVE_SEDIMENT_MODULE)
  geom->sed_t = i_alloc_1d(geom->ewetS + 1);
#endif

  /* Count the lateral ghost cells which are associated with normal */
  /* velocities. These are cells on eastern and western edges for u1 */
  /* velocities and northern and southern edges for u2 velocities.  */
  /* These land boundary locations are organised so that the first */
  /* geom->be#n cells contain the normal velocities (which have a */
  /* zero flux boundary condition) followed by the tangential */
  /* velocities (which have a free/no slip boundary condition.  */
  geom->nbe1 = geom->nbe1S = 0;
  geom->nbe2 = geom->nbe2S = 0;
  for (cc = 1; cc <= geom->nbpt; cc++) {
    c = geom->bpt[cc];
    lc = geom->bin[cc];
    if (geom->xp1[lc] == c || geom->xm1[lc] == c) {
      geom->nbe1 += 1;
      if (cc <= geom->nbptS)
        geom->nbe1S += 1;
    }
    if (geom->yp1[lc] == c || geom->ym1[lc] == c) {
      geom->nbe2 += 1;
      if (cc <= geom->nbptS)
        geom->nbe2S += 1;
    }
  }

#if defined(HAVE_SEDIMENT_MODULE)
  /* Get the locations of the seiment layer vector */
  cc = geom->sgnum - num_sc + 1;
  for (c = 1; c <= geom->ewetS; c++) {
    geom->sed_t[c] = cc;
    cc++;
  }
#endif

  /* Get the locations of the lateral ghost cells */
  for (cc = 1; cc <= geom->ngsed; cc++) {
    c = geom->gsed_t[cc];    /* Sediment ghost */
    geom->wgst[c] = geom->ised_t[cc];
  }
  b1 = geom->nbe1 + 1;
  b2 = geom->nbe2 + 1;
  b3 = geom->nbe1S + 1;
  c3 = geom->nbe2S + 1;
  c1 = c2 = 1;
  for (cc = 1; cc <= geom->nbpt; cc++) {
    c = geom->bpt[cc];
    lc = geom->bin[cc];
    geom->wgst[c] = lc;
    if (geom->xm1[lc] == c) {
      geom->bin2[cc] = geom->xp1[lc];
      /* Do not use the staggered e1 cell for locations where an */
      /* open boundary is adjacent to an OUTSIDE cell in the L_EDGE. */
      i = geom->s2i[lc];
      j = geom->s2j[lc];
      k = geom->s2k[lc];
      if (iswetijk(lc, geom, flg) && flg[k][j][i] & U1BDRY)
	geom->bpte1[c1] = c;
      else
	geom->bpte1[c1] = geom->xp1[c];
      geom->bine1[c1] = geom->xp1[geom->bpte1[c1]];
      if (cc <= geom->nbptS) {
        geom->bpte1S[c1] = geom->bpte1[c1];
        geom->bine1S[c1] = geom->bine1[c1];
      }
      c1++;
    } else if (geom->xp1[lc] == c) {
      geom->bin2[cc] = geom->xm1[lc];
      geom->bpte1[c1] = c;
      geom->bine1[c1] = lc;
      if (cc <= geom->nbptS) {
        geom->bpte1S[c1] = geom->bpte1[c1];
        geom->bine1S[c1] = geom->bine1[c1];
      }
      c1++;
    } else {
      geom->bpte1[b1] = c;
      geom->bine1[b1] = lc;
      if (cc <= geom->nbptS) {
        geom->bpte1S[b3] = geom->bpte1[b1];
        geom->bine1S[b3] = geom->bine1[b1];
        b3++;
      }
      b1++;
    }

    if (geom->ym1[lc] == c) {
      geom->bin2[cc] = geom->yp1[lc];
      /* Do not use the staggered e2 cell for locations where an */
      /* open boundary is adjacent to an OUTSIDE cell in the B_EDGE. */
      i = geom->s2i[lc];
      j = geom->s2j[lc];
      k = geom->s2k[lc];
      if (iswetijk(lc, geom, flg) && flg[k][j][i] & U2BDRY)
	geom->bpte2[c2] = c;
      else
	geom->bpte2[c2] = geom->yp1[c];
      geom->bine2[c2] = geom->yp1[geom->bpte2[c2]];
      if (cc <= geom->nbptS) {
        geom->bpte2S[c2] = geom->bpte2[c2];
        geom->bine2S[c2] = geom->bine2[c2];
      }
      c2++;
    } else if (geom->yp1[lc] == c) {
      geom->bin2[cc] = geom->ym1[lc];
      geom->bpte2[c2] = c;
      geom->bine2[c2] = lc;
      if (cc <= geom->nbptS) {
        geom->bpte2S[c2] = geom->bpte2[c2];
        geom->bine2S[c2] = geom->bine2[c2];
      }
      c2++;
    } else {
      geom->bpte2[b2] = c;
      geom->bine2[b2] = lc;
      if (cc <= geom->nbptS) {
        geom->bpte2S[c3] = geom->bpte2[b2];
        geom->bine2S[c3] = geom->bine2[b2];
        c3++;
      }
      b2++;
    }
  }

  /* Get the radiation stress mask */
  memset(geom->brsm, 0, geom->sgsizS * sizeof(int));
  for (cc = geom->nbe1S + 1; cc <= geom->nbpte1S; cc++) {
    c = geom->bpte1S[cc];
    lc = geom->bine1S[cc];
    if (lc == geom->ym1[c]) geom->brsm[lc] |= (U2BDRY|F_EDGE);
    if (lc == geom->yp1[c]) geom->brsm[lc] |= (U2BDRY|B_EDGE);
  }
  for (cc = geom->nbe2S + 1; cc <= geom->nbpte2S; cc++) {
    c = geom->bpte2S[cc];
    lc = geom->bine2S[cc];
    if (lc == geom->xm1[c]) geom->brsm[lc] |= (U1BDRY|R_EDGE);
    if (lc == geom->xp1[c]) geom->brsm[lc] |= (U1BDRY|L_EDGE);
  }

  /* Fill the cells to process vectors with ghost cells.  */
  /* Cell centered vectors.  */
  geom->n3_t = geom->a3_t + 1;
  geom->n2_t = geom->a2_t + 1;
  for (cc = 1; cc <= geom->nbpt; cc++) {
    c = geom->bpt[cc];
    geom->w3_t[geom->n3_t] = c;
    geom->n3_t++;
    if (cc <= geom->nbptS) {
      geom->w2_t[geom->n2_t] = c;
      geom->n2_t++;
    }
  }
  /* e1 face centered vectors */
  geom->n3_e1 = geom->b3_e1 + 1;
  geom->n2_e1 = geom->b2_e1 + 1;
  for (cc = 1; cc <= geom->nbpt; cc++) {
    c = geom->bpte1[cc];
    geom->w3_e1[geom->n3_e1] = c;
    geom->n3_e1++;
    if (cc <= geom->nbptS) {
      geom->w2_e1[geom->n2_e1] = c;
      geom->n2_e1++;
    }
  }
  /* e2 face centered vectors */
  geom->n3_e2 = geom->b3_e2 + 1;
  geom->n2_e2 = geom->b2_e2 + 1;
  for (cc = 1; cc <= geom->nbpt; cc++) {
    c = geom->bpte2[cc];
    geom->w3_e2[geom->n3_e2] = c;
    geom->n3_e2++;
    if (cc <= geom->nbptS) {
      geom->w2_e2[geom->n2_e2] = c;
      geom->n2_e2++;
    }
  }
  geom->n3_t -= 1;
  geom->n2_t -= 1;
  geom->n3_e1 -= 1;
  geom->n2_e1 -= 1;
  geom->n3_e2 -= 1;
  geom->n2_e2 -= 1;
  if (DEBUG("init_m"))
    dlog("init_m", "Lateral boundary masks OK\n");

  /*-----------------------------------------------------------------*/
  /* Get the 3D - 2D map (including the sediment) */
  geom->m2d = i_alloc_1d(geom->sgsiz);
  for (cc = 1; cc <= geom->enon; cc++)
    geom->m2d[cc] = 0;
  for (cc = 1; cc <= geom->enonS; cc++) {
    c = cc;
    while (geom->zm1[c] && c != geom->zm1[c]) {
      geom->m2d[c] = cc;
      c = geom->zm1[c];
    }
    geom->m2d[c] = cc;
  }

  /* Get the 3D - 2D map for sub-surface ghost cells (i.e. ghost */
  /* cells where the cell above is wet).  */
  for (cc = 1; cc <= geom->enon; cc++) {
    c = cc;
    c1 = geom->zp1[cc];
    while (geom->m2d[cc] == 0 && c != c1) {
      geom->m2d[cc] = geom->m2d[c1];
      c = c1;
      c1 = geom->zp1[c1];
    }
  }

  /* Map boundary ghosts if required */
  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    if (params->compatible & V1957) continue;
    for (cc = 1; cc <= open->no3_t; cc++) {
      c = open->obc_t[cc];
      c1 = geom->m2d[c];
      for(ii = 0; ii < open->bgz; ii++) {
	c = open->omap[c];
	c1 = open->omap[c1];
	geom->m2d[c] = c1;
      }
    }
  }

  if (DEBUG("init_m"))
    dlog("init_m", "3D - 2D maps OK\n");

  /*-----------------------------------------------------------------*/
  /* Set the bathymetry value. This is overwritten from the input */
  /* file for MANUAL operation.  */
  memset(geom->botz, 0, geom->sgsizS * sizeof(double));
  for (c = 1; c <= geom->sgnumS; c++) {
    i = geom->s2i[c];
    j = geom->s2j[c];
    k = geom->s2k[c];
    if (i == NOTVALID && j == NOTVALID && k == NOTVALID)
      continue;                 /* Continue on ghosts */
    if (i < params->nce1 && j < params->nce2) {
      geom->botz[c] = bathy[j][i];
    }
  }
  if (params->sigma) {
    for (n = 0; n < params->nz; n++)
      params->layers[n] /= params->bmax;
  }

  /* Set a no-gradient condition over the sediment for gridz */
  geom->topgrid = layers[nz];
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->bot_t[cc];
    geom->gridz[geom->zm1[c]] = geom->botz[geom->m2d[c]];
  }

  /* Find the deepest coordinte in the grid */
  d1 = 0.0;
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    if (geom->botz[c] < d1) {
      d1 = geom->botz[c];
      geom->cdeep = geom->bot_t[cc];
    }
  }

  /* Count the number of dry cells beneath the surface */
  geom->bdry = 0;
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    i = geom->s2i[c];
    j = geom->s2j[c];
    geom->bdry += kbot[j][i];
  }

  check_sparse(geom, geom->a3_t, num_gc, num_sc, num_se, num_ic, num_dc,
               num_mc, flg);

  if (DEBUG("init_m"))
    dlog("init_m", "Preprocessor done\n\n");

  /*-----------------------------------------------------------------*/
  /* Set up the windows */
  window_build(geom, params);
  write_windows(geom, flg);
  if (DEBUG("init_m"))
    dlog("init_m", "\nWindows created OK\n");

  /*-----------------------------------------------------------------*/
  /* Make and initialise the master data structure */
  master = master_build(params, geom);

  /*-----------------------------------------------------------------*/
  /* Set up the dump data structure */
  dumpdata = dumpdata_build(params, geom, master, flg, layers, bathy);

  /* Set geographical flag */
  geom->is_geog = 0;
  if (!(params->compatible & V5895))
    geom->is_geog = (strlen(params->projection) > 0) &&
      (strcasecmp(params->projection, GEOGRAPHIC_TAG) == 0);
  
  /*-----------------------------------------------------------------*/
  /* Free memory */
  l_free_3d((long ***)flg);
  /*d_free_2d(bathy);*/
}

/* END build_sparse_map()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set up the open boundary vectors                       */
/*-------------------------------------------------------------------*/
void set_OBC_cells(geometry_t *geom,   /* Window geometry            */
		   open_bdrys_t *open, /* Open boundary structure    */
		   open_bdrys_t *io,   /* Input parameter OB structure */
		   int n,       /* Boundary number                   */
		   int nce1,    /* Grid size in e1 direction         */
		   int nce2,    /* Grid size in e2 direction         */
		   int nz,      /* Grid size in vertical direction   */
		   unsigned long ***flag  /* Cell flag array         */
  )
{
  int i, j, k, tn;              /* Open boundary cell locations */
  int ii, jj, kk;               /* Cell centered open boundary locations */
  int ie1, je1, ie2, je2;       /* Cell centered open boundary locations */
  int c, cc, c1, c2;            /* Sparse counters / coordinates */
  int c1n, c1t;                 /* Sparse e1 velocity coordinates */
  int c2n, c2t;                 /* Sparse e2 velocity coordinates */
  short **tag;                  /* Non-zero when boundary cell counted */
  unsigned long bflag;          /* Boundary cell flag */

  /*-----------------------------------------------------------------*/
  /* Initialise the tag array to check for duplicate boundary */
  /* cells.  */
  open->ocodex = open->ocodey = 0;
  open->no2_t = 0;
  open->no2_e1 = 0;
  open->no2_e2 = 0;
  open->to2_e1 = 0;
  open->to2_e2 = 0;
  open->no3_t = 0;
  open->no3_e1 = 0;
  open->no3_e2 = 0;
  open->to3_e1 = 0;
  open->to3_e2 = 0;
  tag = s_alloc_2d(nce1 + 1, nce2 + 1);
  for (j = 0; j < nce2 + 1; j++)
    for (i = 0; i < nce1 + 1; i++)
      tag[j][i] = 0;

  /*-----------------------------------------------------------------*/
  /* Read the number of cells in this boundary */
  open->no2_t = io->npts;

  /*-----------------------------------------------------------------*/
  /* Read i and j indices for this boundary and count the cell */
  /* center and cell face sparse boundary locations.  */
  for (cc = 0; cc < io->npts; cc++) {
    i = io->iloc[cc];
    j = io->jloc[cc];
    k = nz - 1;
    c = c1n = c1t = c2n = c2t = 0;

    /*---------------------------------------------------------------*/
    /* Velocity boundaries. Here b_flag() returns the location of */
    /* the normal open boundary face in ii and jj.  */
    if (open->type == (U1BDRY | U2BDRY)) {
      bflag =
        b_flags(nce1, nce2, nz, flag, i, j, k, &ii, &jj, 1, open->type);
      /* U1BDRY on the right edge */
      if (bflag & R_EDGE && !(tag[j][i] & R_EDGE)) {
        c = geom->map[nz - 1][j][i];
        c1n = geom->map[nz - 1][jj][ii];
        tag[j][i] |= R_EDGE;
        open->no2_e1++;
        open->ocodex = R_EDGE;
      }
      /* U1BDRY on the left edge */
      if (bflag & L_EDGE && !(tag[j][i] & L_EDGE)) {
        c = geom->map[nz - 1][j][i];
        c1n = geom->map[nz - 1][jj][ii];
        tag[j][i] |= L_EDGE;
        open->no2_e1++;
        open->ocodex = L_EDGE;
      }
      /* Set the e2 tangential sparse location if this exists */
      if (!(flag[k][j][i] & (U2BDRY | U2SOLID | U2OUTSIDE))) {
        open->to2_e2++;
        c2t = c;
      }

      bflag =
        b_flags(nce1, nce2, nz, flag, i, j, k, &ii, &jj, 0, open->type);
      /* U2BDRY on the front edge */
      if (bflag & F_EDGE && !(tag[j][i] & F_EDGE)) {
        c = geom->map[nz - 1][j][i];
        c2 = geom->map[nz - 1][jj][ii];
        tag[j][i] |= F_EDGE;
        open->no2_e2++;
        open->ocodey = F_EDGE;
      }
      /* U2BDRY on the back edge */
      if (bflag & B_EDGE && !(tag[j][i] & B_EDGE)) {
        c = geom->map[nz - 1][j][i];
        c2 = geom->map[nz - 1][jj][ii];
        tag[j][i] |= B_EDGE;
        open->no2_e2++;
        open->ocodey = B_EDGE;
      }
      /* Set the e1 tangential sparse location if this exists */
      if (!(flag[k][j][i] & (U1BDRY | U1SOLID | U1OUTSIDE))) {
        open->to2_e1++;
        c1t = c;
      }
    }
    /*---------------------------------------------------------------*/
    /* U1 or U2 boundaries. Here b_flag() returns the location of */
    /* the corresponding cell center in ii and jj.  */
    else {
      /* U1 open boundaries */
      bflag =
        b_flags(nce1, nce2, nz, flag, i, j, k, &ie1, &je1, 1, open->type);
      /* U1BDRY on the right edge */
      if (bflag & R_EDGE && !(tag[j][i] & R_EDGE)) {
        c1n = geom->map[nz - 1][j][i];
        c = geom->map[nz - 1][je1][ie1];
        tag[j][i] |= R_EDGE;
        open->no2_e1++;
        open->ocodex = R_EDGE;
      }
      /* U1BDRY on the left edge */
      if (bflag & L_EDGE && !(tag[j][i] & L_EDGE)) {
        c1n = geom->map[nz - 1][j][i];
        c = geom->map[nz - 1][je1][ie1];
        tag[j][i] |= L_EDGE;
        open->no2_e1++;
        open->ocodex = L_EDGE;
      }
      /* Set the e2 tangential sparse location if this exists */
      if (!(flag[k][je1][ie1] & (U2BDRY | U2SOLID | U2OUTSIDE))) {
        open->to2_e2++;
        c2t = c;
      }
      /* U2 open boundaries */
      bflag =
        b_flags(nce1, nce2, nz, flag, i, j, k, &ie2, &je2, 0, open->type);
      /* U2BDRY on the front edge */
      if (bflag & F_EDGE && !(tag[j][i] & F_EDGE)) {
        c2n = geom->map[nz - 1][j][i];
        c = geom->map[nz - 1][je2][ie2];
        tag[j][i] |= F_EDGE;
        open->no2_e2++;
        open->ocodey = F_EDGE;
      }
      /* U2BDRY on the back edge */
      if (bflag & B_EDGE && !(tag[j][i] & B_EDGE)) {
        c2n = geom->map[nz - 1][j][i];
        c = geom->map[nz - 1][je2][ie2];
        tag[j][i] |= B_EDGE;
        open->no2_e2++;
        open->ocodey = B_EDGE;
      }
      /* Set the e1 tangential sparse location if this exists */
      if (!(flag[k][je2][ie2] & (U1BDRY | U1SOLID | U1OUTSIDE))) {
        open->to2_e1++;
        c1t = c;
      }
    }

    if (!c) {
      /*
      printf("WARNING : Can't find sparse coordinate for boundary cell %d %d\n",
           i, j);
      */
      if (io->type & U1BDRY) {
	ie1 = i;
	je1 = j;
	c1n = c = geom->map[nz - 1][j][i];
        tag[j][i] |= L_EDGE;
        io->ocodex = open->ocodex = LO_EDGE;
	io->bgz = open->bgz = 0;
        open->no2_e1++;
	if (!(flag[k][je1][ie1] & (U2BDRY | U2SOLID | U2OUTSIDE)))
	  c2t = c;
	open->to2_e1--;
      }
      if (io->type & U2BDRY) {
	ie2 = i;
	je2 = j;
	c2n = c = geom->map[nz - 1][j][i];
        tag[j][i] |= B_EDGE;
        io->ocodey = open->ocodey = BO_EDGE;
	io->bgz = open->bgz = 0;
        open->no2_e2++;
	if (!(flag[k][je2][ie2] & (U1BDRY | U1SOLID | U1OUTSIDE)))
	  c1t = c;
	open->to2_e2--;
      }
    }

    /* Loop through the water column and count the 3D tracer cells */
    while (c && c != geom->zm1[c]) {
      open->no3_t++;
      c = geom->zm1[c];
    }
    /* Loop down and count 3D e1 normal velocity cells */
    while (c1n && c1n != geom->zm1[c1n]) {
      open->no3_e1++;
      c1n = geom->zm1[c1n];
    }
    /* Loop down and count 3D e1 tangential velocity cells. Note */
    /* that tangential velocities may occupy a solid face when the */
    /* cell center is wet : in this case the cell should not be */
    /* included in the tangential boundary array. This configuration */
    /* is identified by flag[kk][je1][ie1] = U1SOLID, where kk */
    /* decrements at the same rate as c1t.  */
    kk = k;
    while (c1t && c1t != geom->zm1[c1t] && !(flag[kk][je2][ie2] & U1SOLID)) {
      open->to3_e1++;
      c1t = geom->zm1[c1t];
      kk--;
    }
    /* Loop down and count 3D e2 normal velocity cells */
    while (c2n && c2n != geom->zm1[c2n]) {
      open->no3_e2++;
      c2n = geom->zm1[c2n];
    }
    /* Loop down and count 3D e2 tangential velocity cells */
    kk = k;
    while (c2t && c2t != geom->zm1[c2t] && !(flag[kk][je1][ie1] & U2SOLID)) {
      open->to3_e2++;
      c2t = geom->zm1[c2t];
      kk--;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Allocate menory for the boundary vectors */
  /* Open boundary vectors */
  c1 = open->no3_e1 + open->to3_e1 + 1;
  c2 = open->no3_e2 + open->to3_e2 + 1;
  open->obc_t = i_alloc_1d(open->no3_t + 1);
  open->obc_e1 = i_alloc_1d(c1);
  open->obc_e2 = i_alloc_1d(c2);

  /* Sparse locations one cell into the interior */
  open->oi1_t = i_alloc_1d(open->no3_t + 1);
  open->oi1_e1 = i_alloc_1d(c1);
  open->oi1_e2 = i_alloc_1d(c2);

  /* Sparse locations two cells into the interior */
  open->oi2_t = i_alloc_1d(open->no3_t + 1);
  open->oi2_e1 = i_alloc_1d(c1);
  open->oi2_e2 = i_alloc_1d(c2);

  /* Sparse locations for cyclic boundary conditions */
  open->cyc_t = i_alloc_1d(open->no3_t + 1);
  open->cyc_e1 = i_alloc_1d(c1);
  open->cyc_e2 = i_alloc_1d(c2);

  /* Reset the vector sizes for filling below */
  /* Tangential velocity positions */
  open->to3_e1 = open->no3_e1 + open->to2_e1 + 1;
  open->to3_e2 = open->no3_e2 + open->to2_e2 + 1;
  open->to2_e1 = open->no3_e1 + 1;
  open->to2_e2 = open->no3_e2 + 1;
  /* Normal velocity and cell centered positions */
  open->no3_t = open->no2_t + 1;
  open->no3_e1 = open->no2_e1 + 1;
  open->no3_e2 = open->no2_e2 + 1;
  open->no2_t = 1;
  open->no2_e1 = 1;
  open->no2_e2 = 1;

  /* Reset the tag array to check for duplicate OBC cells */
  for (j = 0; j < nce2 + 1; j++)
    for (i = 0; i < nce1 + 1; i++)
      tag[j][i] = 0;

  /*-----------------------------------------------------------------*/
  /* Assign the sparse coordinates to the boundary vectors. Fill the */
  /* first open->no2_ cells with surface layer coordinates.  */
  /* u1 boundaries.  */
  for (cc = 0; cc < io->npts; cc++) {
    i = io->iloc[cc];
    j = io->jloc[cc];
    k = nz - 1;
    c = c1n = c1t = c2n = c2t = 0;

    /*---------------------------------------------------------------*/
    /* Velocity boundaries. Here b_flag() returns the location of */
    /* the normal open boundary face in ii and jj.  */
    if (open->type == (U1BDRY | U2BDRY)) {
      bflag =
        b_flags(nce1, nce2, nz, flag, i, j, k, &ii, &jj, 1, open->type);
      /* U1BDRY on the right edge */
      if (bflag & R_EDGE && !(tag[j][i] & R_EDGE)) {
        c = geom->map[nz - 1][j][i];
        c1n = geom->map[nz - 1][jj][ii];
        obc_loc(geom, c, geom->xm1, open->obc_t, open->oi1_t, open->oi2_t,
                &open->no2_t, &open->no3_t);
        obc_loc(geom, c1n, geom->xm1, open->obc_e1, open->oi1_e1, open->oi2_e1,
                &open->no2_e1, &open->no3_e1);
        /* Set the e2 tangential sparse location if this exists */
        if (!(flag[k][j][i] & (U2BDRY | U2SOLID | U2OUTSIDE)))
          obc_loc_tan(geom, c, geom->xm1, open->obc_e2, open->oi1_e2,
                      open->oi2_e2, &open->to2_e2, &open->to3_e2, flag,
                      U2SOLID, i, j, k);
        tag[j][i] |= R_EDGE;
      }
      /* U1BDRY on the left edge */
      if (bflag & L_EDGE && !(tag[j][i] & L_EDGE)) {
        c = geom->map[nz - 1][j][i];
        c1n = geom->map[nz - 1][jj][ii];
        obc_loc(geom, c, geom->xp1, open->obc_t, open->oi1_t, open->oi2_t,
                &open->no2_t, &open->no3_t);
        obc_loc(geom, c1n, geom->xp1, open->obc_e1, open->oi1_e1, open->oi2_e1,
                &open->no2_e1, &open->no3_e1);
        /* Set the e2 tangential sparse location if this exists */
        if (!(flag[k][j][i] & (U2BDRY | U2SOLID | U2OUTSIDE)))
          obc_loc_tan(geom, c, geom->xp1, open->obc_e2, open->oi1_e2,
                      open->oi2_e2, &open->to2_e2, &open->to3_e2, flag,
                      U2SOLID, i, j, k);
        tag[j][i] |= L_EDGE;
      }

      bflag =
        b_flags(nce1, nce2, nz, flag, i, j, k, &ii, &jj, 0, open->type);
      /* U2BDRY on the front edge */
      if (bflag & F_EDGE && !(tag[j][i] & F_EDGE)) {
        c = geom->map[nz - 1][j][i];
        c2 = geom->map[nz - 1][jj][ii];
        obc_loc(geom, c, geom->ym1, open->obc_t, open->oi1_t, open->oi2_t,
                &open->no2_t, &open->no3_t);
        obc_loc(geom, c1n, geom->ym1, open->obc_e2, open->oi1_e2, open->oi2_e2,
                &open->no2_e2, &open->no3_e2);
        /* Set the e1 tangential sparse location if this exists */
        if (!(flag[k][j][i] & (U1BDRY | U1SOLID | U1OUTSIDE)))
          obc_loc_tan(geom, c, geom->ym1, open->obc_e1, open->oi1_e1,
                      open->oi2_e1, &open->to2_e1, &open->to3_e1, flag,
                      U1SOLID, i, j, k);
        tag[j][i] |= F_EDGE;
      }
      /* U2BDRY on the back edge */
      if (bflag & B_EDGE && !(tag[j][i] & B_EDGE)) {
        c = geom->map[nz - 1][j][i];
        c2 = geom->map[nz - 1][jj][ii];
        obc_loc(geom, c, geom->yp1, open->obc_t, open->oi1_t, open->oi2_t,
                &open->no2_t, &open->no3_t);
        obc_loc(geom, c1n, geom->yp1, open->obc_e2, open->oi1_e2, open->oi2_e2,
                &open->no2_e2, &open->no3_e2);
        /* Set the e1 tangential sparse location if this exists */
        if (!(flag[k][j][i] & (U1BDRY | U1SOLID | U1OUTSIDE)))
          obc_loc_tan(geom, c, geom->yp1, open->obc_e1, open->oi1_e1,
                      open->oi2_e1, &open->to2_e1, &open->to3_e1, flag,
                      U1SOLID, i, j, k);
        tag[j][i] |= B_EDGE;
      }
    }
    /*---------------------------------------------------------------*/
    /* U1 or U2 boundaries. Here b_flag() returns the location of */
    /* the corresponding cell center in ii and jj.  */
    else {
      /* U1 open boundaries */
      bflag =
        b_flags(nce1, nce2, nz, flag, i, j, k, &ii, &jj, 1, open->type);
      /* U1BDRY on the right edge */
      if (bflag & R_EDGE && !(tag[j][i] & R_EDGE)) {
        c1n = geom->map[nz - 1][j][i];
        c = geom->map[nz - 1][jj][ii];
        obc_loc(geom, c, geom->xm1, open->obc_t, open->oi1_t, open->oi2_t,
                &open->no2_t, &open->no3_t);
        obc_loc(geom, c1n, geom->xm1, open->obc_e1, open->oi1_e1, open->oi2_e1,
                &open->no2_e1, &open->no3_e1);
        /* Set the e2 tangential sparse location if this exists */
        if (!(flag[k][jj][ii] & (U2BDRY | U2SOLID | U2OUTSIDE))) {
          obc_loc_tan(geom, c, geom->xm1, open->obc_e2, open->oi1_e2,
                      open->oi2_e2, &open->to2_e2, &open->to3_e2, flag,
                      U2SOLID, ii, jj, k);
        }
        tag[j][i] |= R_EDGE;
      }

      /* U1BDRY on the left edge */
      if (bflag & L_EDGE && !(tag[j][i] & L_EDGE)) {
        c1n = geom->map[nz - 1][j][i];
        c = geom->map[nz - 1][jj][ii];
        obc_loc(geom, c, geom->xp1, open->obc_t, open->oi1_t, open->oi2_t,
                &open->no2_t, &open->no3_t);
        obc_loc(geom, c1n, geom->xp1, open->obc_e1, open->oi1_e1, open->oi2_e1,
                &open->no2_e1, &open->no3_e1);
        /* Set the e2 tangential sparse location if this exists */
        if (!(flag[k][jj][ii] & (U2BDRY | U2SOLID | U2OUTSIDE)))
          obc_loc_tan(geom, c, geom->xp1, open->obc_e2, open->oi1_e2,
                      open->oi2_e2, &open->to2_e2, &open->to3_e2, flag,
                      U2SOLID, ii, jj, k);
        tag[j][i] |= L_EDGE;
      }
      /* U2 open boundaries */
      bflag =
        b_flags(nce1, nce2, nz, flag, i, j, k, &ii, &jj, 0, open->type);
      /* U2BDRY on the front edge */
      if (bflag & F_EDGE && !(tag[j][i] & F_EDGE)) {
        c2n = geom->map[nz - 1][j][i];
        c = geom->map[nz - 1][jj][ii];
        obc_loc(geom, c, geom->ym1, open->obc_t, open->oi1_t, open->oi2_t,
                &open->no2_t, &open->no3_t);
        obc_loc(geom, c2n, geom->ym1, open->obc_e2, open->oi1_e2, open->oi2_e2,
                &open->no2_e2, &open->no3_e2);
        /* Set the e1 tangential sparse location if this exists */
        if (!(flag[k][jj][ii] & (U1BDRY | U1SOLID | U1OUTSIDE)))
          obc_loc_tan(geom, c, geom->ym1, open->obc_e1, open->oi1_e1,
                      open->oi2_e1, &open->to2_e1, &open->to3_e1, flag,
                      U1SOLID, ii, jj, k);
        tag[j][i] |= F_EDGE;
      }
      /* U2BDRY on the back edge */
      if (bflag & B_EDGE && !(tag[j][i] & B_EDGE)) {
        c2n = geom->map[nz - 1][j][i];
        c = geom->map[nz - 1][jj][ii];
        obc_loc(geom, c, geom->yp1, open->obc_t, open->oi1_t, open->oi2_t,
                &open->no2_t, &open->no3_t);
        obc_loc(geom, c2n, geom->yp1, open->obc_e2, open->oi1_e2, open->oi2_e2,
                &open->no2_e2, &open->no3_e2);
        /* Set the e1 tangential sparse location if this exists */
        if (!(flag[k][jj][ii] & (U1BDRY | U1SOLID | U1OUTSIDE)))
          obc_loc_tan(geom, c, geom->yp1, open->obc_e1, open->oi1_e1,
                      open->oi2_e1, &open->to2_e1, &open->to3_e1, flag,
                      U1SOLID, ii, jj, k);
        tag[j][i] |= B_EDGE;
      }

      if (!c) {
	if (io->type & U1BDRY) {
	  ii = i;
	  jj = j;
	  c1n = c = geom->map[nz - 1][j][i];
	  obc_loc(geom, c, geom->xm1, open->obc_t, open->oi1_t, open->oi2_t,
		  &open->no2_t, &open->no3_t);
	  obc_loc(geom, c1n, geom->xm1, open->obc_e1, open->oi1_e1, open->oi2_e1,
		  &open->no2_e1, &open->no3_e1);
	  if (!(flag[k][jj][ii] & (U2BDRY | U2SOLID | U2OUTSIDE)))
	    obc_loc_tan(geom, c, geom->xp1, open->obc_e2, open->oi1_e2,
			open->oi2_e2, &open->to2_e2, &open->to3_e2, flag,
			U2SOLID, ii, jj, k);
	  tag[j][i] |= L_EDGE;
	}
	if (io->type & U2BDRY) {
	  ii = i;
	  jj = j;
	  c2n = c = geom->map[nz - 1][j][i];
	  obc_loc(geom, c, geom->yp1, open->obc_t, open->oi1_t, open->oi2_t,
		  &open->no2_t, &open->no3_t);
	  obc_loc(geom, c2n, geom->yp1, open->obc_e2, open->oi1_e2, open->oi2_e2,
		  &open->no2_e2, &open->no3_e2);
	  if (!(flag[k][jj][ii] & (U1BDRY | U1SOLID | U1OUTSIDE)))
	    obc_loc_tan(geom, c, geom->yp1, open->obc_e1, open->oi1_e1,
			open->oi2_e1, &open->to2_e1, &open->to3_e1, flag,
			U1SOLID, ii, jj, k);
	  tag[j][i] |= B_EDGE;
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Allocate memory for master - slave transfer vectors */
  if (open->bcond_ele & (FILEIN | CUSTOM) || open->bcond_nor2d & FLATHR) {
    open->transfer_eta = d_alloc_1d(open->no2_t);
  }
  if (open->type & U1BDRY) {
    i = (open->relax_zone_nor) ? open->relax_zone_nor : 1;
    if (open->bcond_nor & (FILEIN | CUSTOM))
      open->transfer_u1 = d_alloc_1d(i * open->no3_e1 + 1);
    i = (open->relax_zone_tan) ? open->relax_zone_tan : 1;
    if (open->bcond_tan & (FILEIN | CUSTOM))
      open->transfer_u2 = d_alloc_1d(i * open->to3_e2 + 1);
    if (open->bcond_nor2d & FILEIN)
      open->transfer_u1av = d_alloc_1d(open->no2_e1 + 1);
    if (open->bcond_nor2d & CUSTOM)
      open->transfer_u1av = d_alloc_1d(open->no3_e1 + 1);
    if (open->bcond_tan2d & (FILEIN | CUSTOM))
      open->transfer_u2av = d_alloc_1d(open->to3_e2 + 1);
  }
  else if (open->type & U2BDRY) {
    i = (open->relax_zone_nor) ? open->relax_zone_nor : 1;
    if (open->bcond_nor & (FILEIN | CUSTOM))
      open->transfer_u2 = d_alloc_1d(i * open->no3_e2 + 1);
    i = (open->relax_zone_tan) ? open->relax_zone_tan : 1;
    if (open->bcond_tan & (FILEIN | CUSTOM))
      open->transfer_u1 = d_alloc_1d(i * open->to3_e1 + 1);
    if (open->bcond_nor2d & FILEIN)
      open->transfer_u2av = d_alloc_1d(open->no2_e2 + 1);
    if (open->bcond_nor2d & CUSTOM)
      open->transfer_u2av = d_alloc_1d(open->no3_e2 + 1);
    if (open->bcond_tan2d & (FILEIN | CUSTOM))
      open->transfer_u1av = d_alloc_1d(open->to3_e1 + 1);
  }

  i = 1;
  open->ntt = 0;
  if (open->ntr)
    open->trm = i_alloc_1d(open->ntr);
  for (tn = 0; tn < open->ntr; tn++) {
    if (open->bcond_tra[tn] & (FILEIN | CUSTOM)) {
      open->trm[tn] = open->ntt;
      if (open->bcond_tra[tn] & (TRCONC|TRCONF)) i = open->bgz;
      open->ntt++;
    }
  }

  if (open->ntt)
    open->t_transfer = d_alloc_2d((i * open->no3_t) + 1, open->ntt);

  /*-----------------------------------------------------------------*/
  /* Set the boundary vectors to the correct size */
  open->no3_t -= 1;
  open->no3_e1 -= 1;
  open->no3_e2 -= 1;
  open->to3_e1 -= 1;
  open->to3_e2 -= 1;
  open->no2_t -= 1;
  open->no2_e1 -= 1;
  open->no2_e2 -= 1;
  open->to2_e1 -= 1;
  open->to2_e2 -= 1;

  s_free_2d(tag);
}

/* END set_OBC_cells()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to fill the open boundary vectors for cell center and     */
/* normal velocity sparse locations.                                 */
/*-------------------------------------------------------------------*/
void obc_loc(geometry_t *geom,  /* Window geometry                   */
	     int c,             /* Cell centered tracer OBC location */
             int *map,          /* Map in the interior direction     */
             int *obc,          /* Open boundary vector              */
             int *oi1,          /* Interior boundary vector : 1 cell */
             int *oi2,          /* Interior boundary vector : 2 cells */
             int *nb2,          /* 2D open boundary counter          */
             int *nb3           /* 3D open boundary counter          */
  )
{
  int c1, c2;                   /* Interior cells */

  /* Get the locations of the interior cells */
  c1 = map[c];
  c2 = map[c1];

  /* Assign a surface boundary sparse coordinate to the first */
  /* locations in the boundary vectors.  */
  obc[*nb2] = c;
  oi1[*nb2] = c1;
  oi2[*nb2] = c2;
  *nb2 = *nb2 + 1;

  /* Assign a all non-surface boundary sparse coordinate to the */
  /* the locations in the boundary vectors subsequent to the 2D end */
  /* position, open->n2d_ */
  /* Loop down and maps 3D tracer cells */
  c = geom->zm1[c];
  while (c != geom->zm1[c]) {
    c1 = map[c];
    c2 = map[c1];
    obc[*nb3] = c;
    oi1[*nb3] = c1;
    oi2[*nb3] = c2;
    *nb3 = *nb3 + 1;

    c = geom->zm1[c];
  }
}

/* END obc_loc()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to fill the open boundary vectors for tangential velocity */
/* sparse locations. Note that tangential velocities may occupy a    */
/* solid face when the cell center is wet : in this case the cell    */
/* should not be included in the tangential boundary array. This     */
/* configuration is identified by flag[k][jj][ii] = U1SOLID, where   */
/* kk decrements at the same rate as c, for e1 tangential velocity.  */
/*-------------------------------------------------------------------*/
void obc_loc_tan(geometry_t *geom, /* Window geometry                */
		 int c,         /* Cell centered tracer OBC location */
                 int *map,      /* Map in the interior direction     */
                 int *obc,      /* Open boundary vector              */
                 int *oi1,      /* Interior boundary vector : 1 cell */
                 int *oi2,      /* Interior boundary vector : 2 cells */
                 int *nb2,      /* 2D open boundary counter          */
                 int *nb3,      /* 3D open boundary counter          */
                 unsigned long ***flag,
                 unsigned long mask, int i, int j, int k)
{
  int c1, c2;                   /* Interior cells */

  /* Get the locations of the interior cells */
  c1 = map[c];
  c2 = map[c1];

  /* Assign a surface boundary sparse coordinate to the first */
  /* locations in the boundary vectors.  */
  obc[*nb2] = c;
  oi1[*nb2] = c1;
  oi2[*nb2] = c2;
  *nb2 = *nb2 + 1;

  /* Assign a all non-surface boundary sparse coordinate to the */
  /* the locations in the boundary vectors subsequent to the 2D end */
  /* position, open->n2d_ */
  /* Loop down and maps 3D tracer cells */
  c = geom->zm1[c];
  k--;
  while (c != geom->zm1[c] && !(flag[k][j][i] & mask)) {
    c1 = map[c];
    c2 = map[c1];
    obc[*nb3] = c;
    oi1[*nb3] = c1;
    oi2[*nb3] = c2;
    *nb3 = *nb3 + 1;
    c = geom->zm1[c];
    k--;
  }
}

/* END obc_loc_tan()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to make the cyclic maps for CYCLED boundary conditions    */
/*-------------------------------------------------------------------*/
void make_cycled(geometry_t *geom,      /* Global geometry           */
		 open_bdrys_t *open,    /* Global open boundary      */
		 unsigned long ***flag  /* Cell flag array           */
		 )
{
  int i, j, k, nz = geom->nz;
  int cc, bb, c, b, c3;

  if (open->ncyc) {

    /* Cell centered cycled locations */
    c3 = open->no2_t + 1;
    for (cc = 1; cc <= open->no2_t; cc++) {
      i = open->ilocc[cc - 1];
      j = open->jlocc[cc - 1];
      c = open->obc_t[cc];
      b = geom->map[nz - 1][j][i];

      if (flag[nz - 1][j][i] & NOTWET)
	hd_quit("make_cycled: CYCLED location (%d %d) is not a wet cell\n", 
		i, j);

      open->cyc_t[cc] = b;
      c = geom->zm1[c];
      b = geom->zm1[b];
      while (c != geom->zm1[c]) {
	open->cyc_t[c3] = b;
	c3++;
	c = geom->zm1[c];
	b = geom->zm1[b];
	if (b == geom->zm1[b])
	  b = geom->zp1[b];
      }
    }

    /* Tangential cycled locations in the e1 direction               */
    c3 = open->to2_e1 + 1;
    for (cc = open->no3_e1 + 1; cc <= open->to2_e1; cc++) {
      c = open->obc_e1[cc];
      i = geom->s2i[c];
      j = geom->s2j[c];
      /* Find the index in obc_t corresponding to the sparse         */
      /* location c.                                                 */
      for (bb = 1; bb <= open->no2_t; bb++)
	if (c == open->obc_t[bb])
	  break;
      /* Get the corresponding cycled sparse location                */
      b = open->cyc_t[bb];
      open->cyc_e1[cc] = b;
      c = geom->zm1[c];
      b = geom->zm1[b];
      k = nz - 2;
      while (c != geom->zm1[c] && !(flag[k][j][i] & U1SOLID)) {
	open->cyc_e1[c3] = b;
	c3++;
	c = geom->zm1[c];
	b = geom->zm1[b];
	if (b == geom->zm1[b])
	  b = geom->zp1[b];
	k--;
      }
    }

    /* Tangential cycled locations in the e2 direction */
    c3 = open->to2_e2 + 1;
    for (cc = open->no3_e2 + 1; cc <= open->to2_e2; cc++) {
      c = open->obc_e2[cc];
      i = geom->s2i[c];
      j = geom->s2j[c];
      /* Find the index in obc_t corresponding to the sparse         */
      /* location c.                                                 */
      for (bb = 1; bb <= open->no2_t; bb++)
	if (c == open->obc_t[bb])
	  break;
      /* Get the corresponding cycled sparse location                */
      b = open->cyc_t[bb];
      open->cyc_e2[cc] = b;
      c = geom->zm1[c];
      b = geom->zm1[b];
      k = nz - 2;
      while (c != geom->zm1[c] && !(flag[k][j][i] & U2SOLID)) {
	open->cyc_e2[c3] = b;
	c3++;
	c = geom->zm1[c];
	b = geom->zm1[b];
	if (b == geom->zm1[b])
	  b = geom->zp1[b];
	k--;
      }
    }

    /* Normal cycled locations in the e1 direction                   */
    c3 = open->no2_e1 + 1;
    for (cc = 1; cc <= open->no2_e1; cc++) {
      c = open->obc_t[cc];
      b = open->cyc_t[cc];
      if (open->ocodex & R_EDGE )
	b = geom->xp1[b];
      open->cyc_e1[cc] = b;
      c = geom->zm1[c];
      b = geom->zm1[b];
      while (c != geom->zm1[c]) {
	open->cyc_e1[c3] = b;
	c3++;
	c = geom->zm1[c];
	b = geom->zm1[b];
	if (b == geom->zm1[b])
	  b = geom->zp1[b];
      }
    }

    /* Normal cycled locations in the e2 direction                   */
    c3 = open->no2_e2 + 1;
    for (cc = 1; cc <= open->no2_e2; cc++) {
      c = open->obc_t[cc];
      b = open->cyc_t[cc];
      if (open->ocodey & F_EDGE )
	b = geom->yp1[b];
      open->cyc_e2[cc] = b;
      c = geom->zm1[c];
      b = geom->zm1[b];
      while (c != geom->zm1[c]) {
	open->cyc_e2[c3] = b;
	c3++;
	c = geom->zm1[c];
	b = geom->zm1[b];
	if (b == geom->zm1[b])
	  b = geom->zp1[b];
      }
    }
  }
}

/* END make_cycled()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets maps to be self-pointing across open boundaries              */
/*-------------------------------------------------------------------*/
void set_map_bdry(geometry_t *geom      /* Global geometry           */
  )
{
  int n;
  int cc, c;
  open_bdrys_t **open = geom->open;

  /* Set any self-mapping maps interior to eastern open boundaries   */
  /* to be self-mapping on the open boundary.                        */
  for (n = 0; n < geom->nobc; n++) {
    if (!(ANY0(TRCONC|FILEIN, open[n]->bcond_tra, open[n]->ntr))) {
      if (open[n]->type & U1BDRY && open[n]->ocodex & R_EDGE) {
	for (cc = 1; cc <= open[n]->no3_e1; cc++) {
	  c = open[n]->obc_e1[cc];
	  geom->xp1[c] = c;
	}
      }
      if (open[n]->type & U2BDRY && open[n]->ocodey & F_EDGE) {
	for (cc = 1; cc <= open[n]->no3_e2; cc++) {
	  c = open[n]->obc_e2[cc];
	  geom->yp1[c] = c;
	}
      }
    }
  }
}

/* END set_map_bdry()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to check the sparse map                                   */
/*-------------------------------------------------------------------*/
void check_sparse(geometry_t *geom, /* Model geometry structure */
                  int nwc,      /* Number of wet cells in the grid */
                  int ngc,      /* Number of ghost cells in the grid */
                  int nsc,      /* Number of sediment cells in the grid */
                  int num_se,   /* Number of straight edge ghost cells */
                  int num_ic,   /* Number of inside corner ghost cells */
                  int num_dc,   /* Number of diagonal connection ghost
                                   cells */
                  int num_mc,   /* Number of multiple ghost cells */
                  unsigned long ***flg  /* Cell status flag */
  )
{
  int i, j, k;                  /* Cartesian indices */
  int ii, jj, kk;               /* Cartesian indices */
  int n, c, cc, cs, ci;         /* Sparse indices */
  char buf1[8];                 /* String buffer */
  char buf2[8];                 /* String buffer */
  open_bdrys_t **open = geom->open; /* OBC data structure */

  if (!DEBUG("init_m"))
    return;

  /* Check that all sparse locations are within the correct range */
  /* and that the Cartesian to sparse maps will map back to the */
  /* correct Cartesian locations.  */

  for (j = 0; j <= geom->nce2; j++)
    for (i = 0; i <= geom->nce1; i++)
      for (k = 0; k < geom->nz; k++) {
        if (flg[k][j][i] & (SOLID | OUTSIDE))
          continue;
        c = geom->map[k][j][i];
        ii = geom->s2i[c];
        jj = geom->s2j[c];
        kk = geom->s2k[c];
        if (c > geom->sgnum || c < 1)
          dlog("init_m", "Invalid sparse coordinate %d : %d %d %d\n", c, i,
               j, k);
        if (ii != i)
          dlog("init_m", "Sparse i mapping error : %d -> %d\n", c, i);
        if (jj != j)
          dlog("init_m", "Sparse j mapping error : %d -> %d (%d %d %d)\n",
               c, j, i, j, k);
        if (kk != k)
          dlog("init_m", "Sparse k mapping error : %d -> %d\n", c, k);
      }

  /* Check that the directional maps are within the correct range */
  for (c = 1; c <= geom->sgnum; c++) {
    ii = geom->s2i[c];
    jj = geom->s2j[c];
    kk = geom->s2k[c];
    ci = geom->xm1[c];
    if (ci > geom->sgnum || ci < 1)
      dlog("init_m", "Invalid west map at %d : %d\n", c, ci);
    ci = geom->xp1[c];
    if (ci > geom->sgnum || ci < 1)
      dlog("init_m", "Invalid east map at %d : %d\n", c, ci);
    ci = geom->yp1[c];
    if (ci > geom->sgnum || ci < 1)
      dlog("init_m", "Invalid north map at %d : %d\n", c, ci);
    ci = geom->ym1[c];
    if (ci > geom->sgnum || ci < 1)
      dlog("init_m", "Invalid south map at %d : %d\n", c, ci);
    ci = geom->zp1[c];
    if (ci > geom->sgnum || ci < 1)
      dlog("init_m", "Invalid up map at %d (%d %d %d) : %d \n", c, ii, jj,
           kk, ci);
    ci = geom->zm1[c];
    if (ci > geom->sgnum || ci < 1)
      dlog("init_m", "Invalid down map at %d (%d %d %d) : %d\n", c, ii, jj,
           kk, ci);
  }

  /* Check that the wet cell map is within the correct range */
  for (ci = 1; ci <= geom->b3_t; ci++) {
    c = geom->w3_t[ci];
    if (c > geom->sgnum || c < 1) {
      dlog("init_m", "Invalid wet cell coordinate %d : %d ", c, ci);
      get_sloc(geom, ci);
      dlog("init_m", "\n");
    }
    i = geom->s2i[c];
    j = geom->s2j[c];
    k = geom->s2k[c];
    if (flg[k][j][i] & NOTWET)
      dlog("init_m", "Dry cell in wet cell mask at %d : %d %d %d\n", c, i,
           j, k);
  }

  /* Check that the wet cell map is within the correct range */
  for (ci = 1; ci <= geom->v3_t; ci++) {
    c = geom->w3_t[ci];
    if (c > geom->sgnum || c < 1) {
      dlog("init_m", "Invalid tracer wet cell coordinate %d : %d ", c, ci);
      get_sloc(geom, c);
      dlog("init_m", "\n");
    }
    i = geom->s2i[c];
    j = geom->s2j[c];
    k = geom->s2k[c];
    ii = i + 1;
    jj = j + 1;

    if (flg[k][j][i] & NOTWET)
      dlog("init_m", "Dry cell in wet cell mask at %d : %d %d %d\n", c, i,
           j, k);

    if (flg[k][j][i] & U1BDRY)
      dlog("init_m", "West boundary in wet cell mask at %d : %d %d %d\n",
           c, i, j, k);
    if (flg[k][j][i] & U2BDRY)
      dlog("init_m", "South boundary in wet cell mask at %d : %d %d %d\n",
           c, i, j, k);
    if (flg[k][j][ii] & U1BDRY)
      dlog("init_m", "East boundary in wet cell mask at %d : %d %d %d\n",
           c, i, j, k);
    if (flg[k][jj][i] & U2BDRY)
      dlog("init_m", "North boundary in wet cell mask at %d : %d %d %d\n",
           c, i, j, k);
  }

  /* Check that any spatial map will always map back to the original */
  /* sparse coordinate. This is only performed for wet cells since */
  /* it is possible for multiple ghost cells to point at each other */
  /* if a directional map is reversed.  */
  for (ci = 1; ci <= geom->a3_t; ci++) {
    c = geom->w3_t[ci];
    cc = geom->yp1[c];
    cs = geom->ym1[cc];
    if (cs != c && cc != c) {
      dlog("init_m", "Invalid north-south map at %d (%d %d %d)\n",
           c, geom->s2i[c], geom->s2j[c], geom->s2k[c]);
      get_sloc(geom, c);
      dlog("init_m", "              north map : %d (%d %d %d)\n",
           cc, geom->s2i[cc], geom->s2j[cc], geom->s2k[cc]);
      dlog("init_m", "               south map : %d (%d %d %d)\n",
           cs, geom->s2i[cs], geom->s2j[cs], geom->s2k[cs]);
    }
    cc = geom->ym1[c];
    cs = geom->yp1[cc];
    if (cs != c && cc != c) {
      dlog("init_m", "Invalid south-north map at %d (%d %d %d)\n",
           c, geom->s2i[c], geom->s2j[c], geom->s2k[c]);
      get_sloc(geom, c);
      dlog("init_m", "              south map : %d (%d %d %d)\n",
           cc, geom->s2i[cc], geom->s2j[cc], geom->s2k[cc]);
      dlog("init_m", "               north map : %d (%d %d %d)\n",
           cs, geom->s2i[cs], geom->s2j[cs], geom->s2k[cs]);
    }
    cc = geom->xm1[c];
    cs = geom->xp1[cc];
    if (cs != c && cc != c) {
      dlog("init_m", "Invalid west-east map at %d (%d %d %d)\n",
           c, geom->s2i[c], geom->s2j[c], geom->s2k[c]);
      dlog("init_m", "              west map : %d (%d %d %d)\n",
           cc, geom->s2i[cc], geom->s2j[cc], geom->s2k[cc]);
      dlog("init_m", "               east map : %d (%d %d %d)\n",
           cs, geom->s2i[cs], geom->s2j[cs], geom->s2k[cs]);
    }
    cc = geom->xp1[c];
    cs = geom->xm1[cc];
    if (cs != c && cc != c) {
      dlog("init_m", "Invalid east-west map at %d (%d %d %d)\n",
           c, geom->s2i[c], geom->s2j[c], geom->s2k[c]);
      dlog("init_m", "              east map : %d (%d %d %d)\n",
           cc, geom->s2i[cc], geom->s2j[cc], geom->s2k[cc]);
      dlog("init_m", "               west map : %d (%d %d %d)\n",
           cs, geom->s2i[cs], geom->s2j[cs], geom->s2k[cs]);
    }
    cc = geom->zp1[c];
    cs = geom->zm1[cc];
    /* 
       if(cs!=c && cc!=c) { dlog("init_m", "Invalid up-down map at %d (%d
       %d %d)\n", c,geom->s2i[c],geom->s2j[c],geom->s2k[c]);
       dlog("init_m", " up map : %d (%d %d %d)\n",
       cc,geom->s2i[cc],geom->s2j[cc],geom->s2k[cc]); dlog("init_m", "
       down map : %d (%d %d %d)\n",
       cs,geom->s2i[cs],geom->s2j[cs],geom->s2k[cs]); } */
    cc = geom->zm1[c];
    cs = geom->zp1[cc];
    if (cs != c && cc != c) {
      dlog("init_m", "Invalid down-up map at %d (%d %d %d)\n",
           c, geom->s2i[c], geom->s2j[c], geom->s2k[c]);
      dlog("init_m", "              down map : %d (%d %d %d)\n",
           cc, geom->s2i[cc], geom->s2j[cc], geom->s2k[cc]);
      dlog("init_m", "              up map : %d (%d %d %d)\n",
           cs, geom->s2i[cs], geom->s2j[cs], geom->s2k[cs]);
    }
  }

  /* Print information relating to the map */
  dlog("init_m", "\nGrid diagnostics\n");
  cc = geom->nce1 * geom->nce2 * geom->nz;
  dlog("init_m",
       "Number of 3D sparse cells = %d (%d %% of 3D Cartesian grid)\n",
       geom->sgnum, 100 * geom->sgnum / cc);
  cc = geom->nce1 * geom->nce2;
  dlog("init_m",
       "Number of 2D sparse cells = %d (%d %% of 2D Cartesian grid)\n",
       geom->sgnumS, 100 * geom->sgnumS / cc);
  dlog("init_m", "Number of wet cells = %d (%5.2f %% sparse cells)\n", nwc,
       100 * (double)nwc / (double)geom->sgnum);
  dlog("init_m", "Number of non-boundary wet cells = %d (%5.2f %%)\n",
       geom->v3_t, 100 * (double)geom->v3_t / (double)geom->sgnum);
  dlog("init_m", "Number of non-boundary e1 cells = %d (%5.2f %%)\n",
       geom->v3_e1, 100 * (double)geom->v3_e1 / (double)geom->sgnum);
  dlog("init_m", "Number of non-boundary e2 cells = %d (%5.2f %%)\n",
       geom->v3_e2, 100 * (double)geom->v3_e2 / (double)geom->sgnum);

  dlog("init_m", "Number of 2D non-boundary wet cells = %d (%5.2f %%)\n",
       geom->v2_t, 100 * (double)geom->v2_t / (double)geom->sgnumS);
  dlog("init_m", "Number of 2D non-boundary e1 cells = %d (%5.2f %%)\n",
       geom->v2_e1, 100 * (double)geom->v2_e1 / (double)geom->sgnumS);
  dlog("init_m", "Number of 2D non-boundary e2 cells = %d (%5.2f %%)\n",
       geom->v2_e2, 100 * (double)geom->v2_e2 / (double)geom->sgnumS);

  dlog("init_m", "Number of ghost cells = %d (%5.2f %%)\n", ngc,
       100 * (double)ngc / (double)geom->sgnum);
  dlog("init_m", "Number of sediment cells = %d (%5.2f %% ghost cells)\n",
       nsc, 100 * (double)nsc / (double)ngc);
  dlog("init_m", "Number of straight edge ghost cells = %d (%5.2f %%)\n",
       num_se, 100 * (double)num_se / (double)ngc);
  dlog("init_m", "Number of inside corner ghost cells = %d (%5.2f %%)\n",
       num_ic, 100 * (double)num_ic / (double)ngc);
  dlog("init_m", "Number of diagonal ghost cells = %d (%5.2f %%)\n",
       num_dc, 100 * (double)num_dc / (double)ngc);
  dlog("init_m", "Number of multiple ghost cells = %d (%5.2f %%)\n",
       num_mc, 100 * (double)num_mc / (double)ngc);
  if (ngc != nsc + num_se + num_ic + num_dc)
    dlog("init_m",
         "WARNING : Sum of ghost cell constituents != number of ghost cells\n");
  dlog("init_m", "\n");
  for (n = 0; n < geom->nobc; n++) {
    if (open[n]->type & U1BDRY)
      dlog("init_m", "Boundary %d (%s) : U1\n", n, open[n]->name);
    else
      dlog("init_m", "Boundary %d (%s) : U2\n", n, open[n]->name);
    strcpy(buf1, "");
    strcpy(buf2, "");
    if (open[n]->ocodex & L_EDGE)
      strcpy(buf1, "west");
    if (open[n]->ocodex & R_EDGE)
      strcpy(buf1, "east");
    if (open[n]->ocodey & F_EDGE)
      strcpy(buf2, "north");
    if (open[n]->ocodey & B_EDGE)
      strcpy(buf2, "south");
    dlog("init_m", "Orientation = %s %s edge\n", buf2, buf1);
    dlog("init_m", "Number of 3D u1 cells = %d\n", open[n]->no3_e1);
    dlog("init_m", "Number of 3D u2 cells = %d\n", open[n]->no3_e2);
    dlog("init_m", "Number of 3D tracer cells = %d\n", open[n]->no3_t);
    dlog("init_m", "Number of 2D u1 cells = %d\n", open[n]->no2_e1);
    dlog("init_m", "Number of 2D u2 cells = %d\n", open[n]->no2_e2);
    dlog("init_m", "Number of 2D elevation cells = %d\n", open[n]->no2_t);
    dlog("init_m", "\n");
  }
}

/* END check_sparse()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to check for straight edge / outside corner ghost cells.  */
/* This can be done assuming cell (i,j,k) is a wet cell (mode=0) or  */
/* a ghost cell (mode=1).                                            */
/*-------------------------------------------------------------------*/
unsigned long se_oc(geometry_t *geom, /* Model geometry structure */
                    int i, int j, int k,  /* Indicies */
                    int *xp1, int *xm1, /* x direction spatial indicies */
                    int *yp1, int *ym1, /* y direction spatial indicies */
                    int mode    /* Calculation mode */
  )
{
  unsigned long ret = 0;

  *xp1 = 0;
  *yp1 = 0;
  *xm1 = 0;
  *ym1 = 0;

  if (i < geom->nce1)
    *xp1 = geom->map[k][j][i + 1];
  if (j < geom->nce2)
    *yp1 = geom->map[k][j + 1][i];
  if (i > 0)
    *xm1 = geom->map[k][j][i - 1];
  if (j > 0)
    *ym1 = geom->map[k][j - 1][i];

  if (mode) {
    if (geom->map[k][j][i] == 0) {
      /* Check for ghost cells on west edges and west outside corners */
      if (i < geom->nce1 && *xp1 != 0)
        ret |= (W_SE);
      /* Check for ghost cells on east edges and east outside corners */
      if (i > 0 && *xm1 != 0)
        ret |= (E_SE);
      /* Check for ghost cells on south edges and south outside corners */
      if (j < geom->nce2 && *yp1 != 0)
        ret |= S_SE;
      /* Check for ghost cells on north edges and north outside corners */
      if (j > 0 && *ym1 != 0)
        ret |= (N_SE);
    } else
      ret |= (WET);
  } else {
    if (geom->map[k][j][i] != 0) {
      ret |= (WET);
      /* Check for ghost cells on west edges and west outside corners */
      if (i < geom->nce1 && *xp1 == 0)
        ret |= (E_SE);
      /* Check for ghost cells on east edges and east outside corners */
      if (i > 0 && *xm1 == 0)
        ret |= (W_SE);
      /* Check for ghost cells on south edges and south outside corners */
      if (j < geom->nce2 && *yp1 == 0)
        ret |= (N_SE);
      /* Check for ghost cells on north edges and north outside corners */
      if (j > 0 && *ym1 == 0)
        ret |= (S_SE);
    }
  }
  return (ret);
}

/* END se_oc()                                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to check for diagonal / inside corner ghost cells         */
/*-------------------------------------------------------------------*/
unsigned long di_ic(geometry_t *geom, /* Model geometry structure */
                    int i, int j, int k,  /* Indicies */
                    int *xp1, int *xm1, /* x direction spatial indicies */
                    int *yp1, int *ym1  /* y direction spatial indicies */
  )
{
  unsigned long ret = 0;
  *xp1 = 0;
  *yp1 = 0;
  *xm1 = 0;
  *ym1 = 0;

  if (i < geom->nce1)
    *xp1 = geom->map[k][j][i + 1];
  if (j < geom->nce2)
    *yp1 = geom->map[k][j + 1][i];
  if (i > 0)
    *xm1 = geom->map[k][j][i - 1];
  if (j > 0)
    *ym1 = geom->map[k][j - 1][i];

  if (geom->map[k][j][i] != 0) {
    if (i < geom->nce1) {
      if (j < geom->nce2)
        /* Check for north-east diagonals */
        if (geom->map[k][j + 1][i + 1] != 0 && *xp1 == 0 && *yp1 == 0)
          ret |= (SW_D);
      if (j > 0)
        /* Check for south-east diagonals */
        if (geom->map[k][j - 1][i + 1] != 0 && *xp1 == 0 && *ym1 == 0)
          ret |= (NW_D);
    }
    if (i > 0) {
      if (j < geom->nce2)
        /* Check for north-west diagonals */
        if (geom->map[k][j + 1][i - 1] != 0 && *xm1 == 0 && *yp1 == 0)
          ret |= (SE_D);
      if (j > 0)
        /* Check for south-west diagonals */
        if (geom->map[k][j - 1][i - 1] != 0 && *xm1 == 0 && *ym1 == 0)
          ret |= (NE_D);
    }
  }

  if (geom->map[k][j][i] == 0) {
    if (i < geom->nce1) {
      if (j < geom->nce2)
        /* Check for south-west inside corners */
        if (geom->map[k][j + 1][i + 1] != 0 && *xp1 == 0 && *yp1 == 0)
          ret |= (SW_IC);
      if (j > 0)
        /* Check for north-west inside corners */
        if (geom->map[k][j - 1][i + 1] != 0 && *xp1 == 0 && *ym1 == 0)
          ret |= (NW_IC);
    }
    if (i > 0) {
      if (j < geom->nce2)
        /* Check for south-east inside corners */
        if (geom->map[k][j + 1][i - 1] != 0 && *xm1 == 0 && *yp1 == 0)
          ret |= (SE_IC);
      if (j > 0)
        /* Check for north-east inside corners */
        if (geom->map[k][j - 1][i - 1] != 0 && *xm1 == 0 && *ym1 == 0)
          ret |= (NE_IC);
    }
  }
  return (ret);
}

/* END di_ic()                                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Prints the Cartesian location given a sparse location             */
/*-------------------------------------------------------------------*/
void get_sloc(geometry_t *geom, /* Model geometry structure */
              int c             /* Sparse location to map */
  )
{
  int i, j, k;

  for (j = 0; j < geom->nfe2; j++)
    for (i = 0; i < geom->nfe1; i++)
      for (k = 0; k < geom->nz; k++) {
        if (geom->map[k][j][i] == (unsigned long)c)
          dlog("init_m", "(%d %d %d)", i, j, k);
      }
}

/* END get_sloc()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to map into the interior until an open boundary is found  */
/* (the map becomes self mapping : the map[map[c]] is used to also   */
/* account for inner explicit maps) then set up the interior sparse  */
/* location for normal cyclic open boundaries.                       */
/*-------------------------------------------------------------------*/
int cyclic_m1(geometry_t *geom, /* Window geometry                   */
	      int codex,        /* Boundary code in e1 direction     */
              int codey,        /* Boundary code in e2 direction     */
              int *map,         /* Map in direction normal to boundary */
              int c             /* Self mapping sparse location      */
  )
{
  int cyc;                      /* Cyclic sparse coordinate */
  cyc = c;

  if (codex) {
    while (c != map[map[c]])
      c = map[c];
    if (codex & L_EDGE)
      cyc = geom->xm1[c];
    if (codex & R_EDGE)
      cyc = geom->xp1[c];
  } else if (codey) {
    while (c != map[map[c]])
      c = map[c];
    if (codey & B_EDGE)
      cyc = geom->ym1[c];
    if (codey & F_EDGE)
      cyc = geom->yp1[c];
  }
  return (cyc);
}

/* END cyclic_m1()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to map into the interior until an open boundary is found  */
/* (the map becomes self mapping) then set up the interior sparse    */
/* location for tracer and tangential cyclic open boundaries.        */
/*-------------------------------------------------------------------*/
int cyclic_m2(geometry_t *geom, /* Window geometry                   */
	      int codex,        /* Boundary code in e1 direction */
              int codey,        /* Boundary code in e2 direction */
              int *map,         /* Map in direction normal to boundary */
              int c             /* Self mapping sparse location */
  )
{
  int cyc;                      /* Cyclic sparse coordinate */
  cyc = c;

  if (codex) {
    while (c != map[map[c]])
      c = map[c];
    if (codex & L_EDGE)
      cyc = geom->xm1[geom->xm1[c]];
    if (codex & R_EDGE)
      cyc = geom->xp1[c];
  } else if (codey) {
    while (c != map[map[c]])
      c = map[c];
    if (codey & B_EDGE)
      cyc = geom->ym1[geom->ym1[c]];
    if (codey & F_EDGE)
      cyc = geom->yp1[c];
  }
  return (cyc);
}

/* END cyclic_m2()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to establish which type of boundaries surround the cell   */
/* centered at location (i,j).                                       */
/*-------------------------------------------------------------------*/
unsigned long b_flags(int nce1, /* Size of the grid in the e1 direction */
                      int nce2, /* Size of the grid in the e2 direction */
                      int nz,   /* Size of the grid in the z direction */
                      unsigned long ***flag,  /* Flags array */
                      int i,    /* Cell centered x location */
                      int j,    /* Cell centered y location */
                      int k,    /* Cell centered z location */
                      int *fi,  /* x location of boundary cell face */
                      int *fj,  /* y location of boundary cell face */
                      int mode, /* mode=1; set U1 location, mode=0; set U2 
                                   location */
                      int btype /* Type of boundary edge */
  ) {
  unsigned long bt = 0;
  unsigned long mask, omask = (OUTSIDE | SOLID);
  int mi, mj, pi, pj, insf;

  *fi = i;
  *fj = j;
  if (flag[k][j][i] & SOLID)
    return (SOLID);

  /*-----------------------------------------------------------------*/
  /* Set the boundary flag for VELOCITY boundaries */
  if (btype == (U1BDRY | U2BDRY)) {
    mask = (OUTSIDE | U1OUTSIDE);
    /* Check for U1BDRY on the R_EDGE */
    if (i >= 0 && i < nce1 - 1) {
      pi = i + 1;
      if (flag[k][j][pi] & mask && !(flag[k][j][i] & mask)) {
        bt |= R_EDGE;
        flag[k][j][i] |= (U1BDRY | R_EDGE);
        if (flag[k][j][i] & U1SOLID)
          flag[k][j][i] &= ~U1SOLID;
        if (mode)
          *fi = pi;
      }
    } else if (i == nce1 - 1) {
      pi = i + 1;
      if (!(flag[k][j][i] & mask)) {
        bt |= R_EDGE;
        flag[k][j][i] |= (U1BDRY | R_EDGE);
        if (flag[k][j][i] & U1SOLID)
          flag[k][j][i] &= ~U1SOLID;
        if (mode)
          *fi = pi;
      }
    }

    /* Check for U1BDRY on the L_EDGE */
    if (i > 0 && i <= nce1 - 1) {
      mi = i;
      if (flag[k][j][mi] & mask && !(flag[k][j][i] & mask)) {
        bt |= L_EDGE;
        flag[k][j][i] |= (U1BDRY | L_EDGE);
        if (flag[k][j][i] & U1SOLID)
          flag[k][j][i] &= ~U1SOLID;
        if (mode)
          *fi = i;
      }
    } else if (i == 0) {
      if (!(flag[k][j][i] & mask)) {
        bt |= L_EDGE;
        flag[k][j][i] |= (U1BDRY | L_EDGE);
        if (flag[k][j][i] & U1SOLID)
          flag[k][j][i] &= ~U1SOLID;
        if (mode)
          *fi = i;
      }
    }

    /* Check for U2BDRY on the F_EDGE */
    mask = (OUTSIDE | U2OUTSIDE);
    if (j >= 0 && j < nce2 - 1) {
      pj = j + 1;
      if (flag[k][pj][i] & mask && !(flag[k][j][i] & mask)) {
        bt |= F_EDGE;
        flag[k][j][i] |= (U2BDRY | F_EDGE);
        if (flag[k][j][i] & U2SOLID)
          flag[k][j][i] &= ~U2SOLID;
        if (!mode)
          *fj = pj;
      }
    } else if (j == nce2 - 1) {
      pj = j + 1;
      if (!(flag[k][j][i] & mask)) {
        bt |= F_EDGE;
        flag[k][j][i] |= (U2BDRY | F_EDGE);
        if (flag[k][j][i] & U2SOLID)
          flag[k][j][i] &= ~U2SOLID;
        if (!mode)
          *fj = pj;
      }
    }

    /* Check for U2BDRY on the B_EDGE */
    if (j > 0 && j <= nce2 - 1) {
      mj = j;
      if (flag[k][mj][i] & mask && !(flag[k][j][i] & mask)) {
        bt |= B_EDGE;
        flag[k][j][i] |= (U2BDRY | B_EDGE);
        if (flag[k][j][i] & U2SOLID)
          flag[k][j][i] &= ~U2SOLID;
        if (!mode)
          *fj = j;
      }
    } else if (j == 0) {
      if (!(flag[k][j][i] & mask)) {
        bt |= B_EDGE;
        flag[k][j][i] |= (U2BDRY | B_EDGE);
        if (flag[k][j][i] & U2SOLID)
          flag[k][j][i] &= ~U2SOLID;
        if (!mode)
          *fj = j;
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the boundary flag for U1BDRY boundaries */
  else if (btype == U1BDRY) {
    mask = (OUTSIDE | U1OUTSIDE | U1SOLID | U1BDRY);
    if (!(flag[k][j][i] & mask)) {
      if (flag[k][j][i] & R_EDGE)
        bt |= R_EDGE;
      if (flag[k][j][i] & L_EDGE)
        bt |= L_EDGE;
    }
    /* Check for U1BDRY on the R_EDGE. Note: include the flag insf   */
    /* for the case where U2BDRY|F_EDGE OBC's exist on top-left      */
    /* inside corners.                                               */
    insf = 0;
    if (i > 0 && i <= nce1) {
      if ((flag[k][j][i-1] & U2BDRY) && (flag[k][j][i-1] & F_EDGE)) insf = 1;
      mi = i - 1;
      if (flag[k][j][i] & mask && !(flag[k][j][i - 1] & omask) && !insf) {
        bt |= R_EDGE;
        flag[k][j][i] |= (U1BDRY | R_EDGE);
        if (flag[k][j][i] & U1SOLID)
          flag[k][j][i] &= ~U1SOLID;
        if (flag[k][j][i] & U1OUTSIDE)
          flag[k][j][i] &= ~U1OUTSIDE;
        /* This is set to wet at this stage since it must be */
        /* included in the wet domain. Reset to OUTSIDE prior to */
        /* writing an input file in AUTO mode.  */
	/*
        if (flag[k][j][i] & OUTSIDE)
          flag[k][j][i] &= ~OUTSIDE;
	*/
        if (mode)
          *fi = mi;
      }
    }
    /* Check for U1BDRY on the L_EDGE. Note: include the flag insf   */
    /* for the case where U2BDRY|F_EDGE OBC's exist on top-right     */
    /* inside corners.                                               */
    insf = 0;
    if (i > 0 && i < nce1) {
      if ((flag[k][j][i+1] & U2BDRY) && (flag[k][j][i+1] & F_EDGE)) insf = 1;
      pi = i;
      if (flag[k][j][i] & mask && !(flag[k][j][i + 1] & omask) && !insf) {
        bt |= L_EDGE;
        flag[k][j][i] |= (U1BDRY | L_EDGE);
        if (flag[k][j][i] & U1SOLID)
          flag[k][j][i] &= ~U1SOLID;
        if (flag[k][j][i] & U1OUTSIDE)
          flag[k][j][i] &= ~U1OUTSIDE;
        if (mode)
          *fi = pi;
      }
    }
    if (i == 0) {
      bt |= L_EDGE;
      flag[k][j][i] |= (U1BDRY | L_EDGE);
      if (flag[k][j][i] & U1SOLID)
        flag[k][j][i] &= ~U1SOLID;
      if (flag[k][j][i] & U1OUTSIDE)
        flag[k][j][i] &= ~U1OUTSIDE;
      if (mode)
        *fi = 0;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the boundary flag for U2BDRY boundaries */
  else if (btype == U2BDRY) {
    mask = (OUTSIDE | U2OUTSIDE | U2SOLID | U2BDRY);
    if (!(flag[k][j][i] & mask)) {
      if (flag[k][j][i] & F_EDGE)
        bt |= F_EDGE;
      if (flag[k][j][i] & B_EDGE)
        bt |= B_EDGE;
    }
    /* Check for U2BDRY on the F_EDGE. Note: include the flag insf   */
    /* for the case where U1BDRY|R_EDGE OBC's exist on bottom-right  */
    /* inside corners.                                               */
    insf = 0;
    if (j > 0 && j <= nce2) {
      if ((flag[k][j-1][i] & U1BDRY) && (flag[k][j-1][i] & R_EDGE)) insf = 1;
      mj = j - 1;
      if (flag[k][j][i] & mask && !(flag[k][j - 1][i] & omask) && !insf) {
        bt |= F_EDGE;
        flag[k][j][i] |= (U2BDRY | F_EDGE);
        if (flag[k][j][i] & U2SOLID)
          flag[k][j][i] &= ~U2SOLID;
        if (flag[k][j][i] & U2OUTSIDE)
          flag[k][j][i] &= ~U2OUTSIDE;
	/*
        if (flag[k][j][i] & OUTSIDE)
          flag[k][j][i] &= ~OUTSIDE;
	*/
        if (!mode)
          *fj = mj;
      }
    }
    /* Check for U2BDRY on the B_EDGE. Note: include the flag insf   */
    /* for the case where U1BDRY|R_EDGE OBC's exist on top-right     */
    /* inside corners.                                               */
    insf = 0;
    if (j > 0 && j < nce2) {
      if ((flag[k][j+1][i] & U1BDRY) && (flag[k][j+1][i] & R_EDGE)) insf = 1;
      pj = j;
      if (flag[k][j][i] & mask && !(flag[k][j + 1][i] & omask) && !insf) {
        bt |= B_EDGE;
        flag[k][j][i] |= (U2BDRY | B_EDGE);
        if (flag[k][j][i] & U2SOLID)
          flag[k][j][i] &= ~U2SOLID;
        if (flag[k][j][i] & U2OUTSIDE)
          flag[k][j][i] &= ~U2OUTSIDE;
        if (!mode)
          *fj = pj;
      }
    }
    if (j == 0) {
      bt |= B_EDGE;
      flag[k][j][i] |= (U2BDRY | B_EDGE);
      if (flag[k][j][i] & U2SOLID)
        flag[k][j][i] &= ~U2SOLID;
      if (flag[k][j][i] & U2OUTSIDE)
        flag[k][j][i] &= ~U2OUTSIDE;
      if (!mode)
        *fj = 0;
    }
  }

  return (bt);
}

/* END b_flags()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set a no-gradient condition on kbot over OUTSIDE cells */
/* on R_EDGE and F_EDGE open boundaries.                             */
/*-------------------------------------------------------------------*/
int
set_kbot(int i, int j, unsigned long **flag, short **kbot, int btype,
         int nce1, int nce2)
{
  unsigned long mask, omask = (OUTSIDE | SOLID);
  int mi, mj;

  if (btype & U1BDRY) {
    mask = (OUTSIDE | U1OUTSIDE | U1SOLID | U1BDRY);
    if (i > 0 && i <= nce1) {
      mi = i - 1;
      if (flag[j][i] & mask && !(flag[j][mi] & omask)) {
        return (kbot[j][mi]);
      }
    }
  }
  if (btype & U2BDRY) {
    mask = (OUTSIDE | U2OUTSIDE | U2SOLID | U2BDRY);
    if (j > 0 && j <= nce2) {
      mj = j - 1;
      if (flag[j][i] & mask && !(flag[mj][i] & omask))
        return (kbot[mj][i]);
    }
  }
  return (kbot[j][i]);
}

/* END set_kbot()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read the batymetry and set SOLID and OUTSIDE cells     */
/*-------------------------------------------------------------------*/
void make_flags(parameters_t *params, /* Input parameters data structure */
                unsigned long ***flag,  /* Flags array */
                double **bathy, /* Bathymetry array */
                double *layers, /* Layer spacing array */
                int nce1,       /* Size of the grid in the e1 direction */
                int nce2,       /* Size of the grid in the e2 direction */
                int nz          /* Size of the grid in the z direction */
  )
{
  int i, j, k;                  /* Counters */
  int is, ie, js, je;           /* Limits of grid */
  int xsize, ysize;             /* Size of original grid */
  int bathylimit = 0;           /* Flag to set min and max bathymetry */
  int percent = 0;              /* Flag to set maximum elevation */
  double min_cell_thickness = 0.0;  /* Minimum cell thickness */
  unsigned long topflag = 0;    /* Surface flag value */
  double bmin, bmax;            /* Minimum and maximum bathymetry */
  double etamax;                /* Maximum surface elevation */
  double *mincelldz;            /* Minimum layer thickness */

  /*-----------------------------------------------------------------*/
  /* Read parameters required for adjusting the batymetry */
  is = js = 0;
  ie = xsize = nce1;
  je = ysize = nce2;
  if (params->edgef) {
    is = js = params->edgef;
    ie = nce1 - params->edgef;
    je = nce2 - params->edgef;
    xsize = nce1 - 2 * params->edgef;
    ysize = nce2 - 2 * params->edgef;
  }
  etamax = params->etamax;
  min_cell_thickness = atof(params->mct);
  if (params->mct[strlen(params->mct) - 1] == '%')
    percent = 1;
  mincelldz = d_alloc_1d(nz);
  for (k = 0; k < nz; k++) {
    if (percent)
      mincelldz[k] =
        (layers[k + 1] - layers[k]) * min_cell_thickness / 100.0;
    else
      mincelldz[k] = min_cell_thickness;
  }

  /*-----------------------------------------------------------------*/
  /* Read the batymetry array */
  bmin = bmax = 0.0;
  if (params->bmin)
    bmin = params->bmin;
  if (params->bmax)
    bmax = params->bmax;
  if (bmin && bmax) {
    if (bmin > 0 && bmax > 0 && bmin > bmax)
      hd_quit("make_flags: BATHYMIN > BATHYMAX\n");
    else
      bathylimit = 1;
  }

  /*-----------------------------------------------------------------*/
  /* Initialise */
  for (j = 0; j <= nce2; j++)
    for (i = 0; i <= nce1; i++)
      for (k = 0; k < nz; k++) {
        flag[k][j][i] = 0;
      }

  /*-----------------------------------------------------------------*/
  /* Set the flags to solid on the nce1 and nce2 edges */
  /* flag = l_alloc_3d(nce1+1,nce2+1,nz); */
  for (k = 0; k < nz; k++) {
    i = nce1;
    for (j = 0; j < nce2; j++) {
      if (bathy[j][i - 1] >= layers[k + 1] || bathy[j][i - 1] >= etamax)
        flag[k][j][i] |= SOLID;
      else
        flag[k][j][i] |= U1OUTSIDE;
    }
    j = nce2;
    for (i = 0; i < nce1; i++) {
      if (bathy[j - 1][i] >= layers[k + 1] || bathy[j - 1][i] >= etamax)
        flag[k][j][i] |= SOLID;
      else
        flag[k][j][i] |= U2OUTSIDE;
    }
  }

  for (j = 0; j < nce2; j++) {
    for (i = 0; i < nce1; i++) {
      /* Adjust botz to min and max bathymetry, checking for outside */
      /* cells (botz < gridz[0]) and land cells (botz > etamax).  */
      if (bathylimit && bathy[j][i] < etamax && bathy[j][i] > -bmin)
        bathy[j][i] = -bmin;
      if (bathylimit && bathy[j][i] >= layers[0] && bathy[j][i] < -bmax)
        bathy[j][i] = -bmax;

      /* Mark outside cells */
      if (bathy[j][i] < layers[0]) {
        flag[nz - 1][j][i] |= OUTSIDE;
      }

      /* Load the cell values */
      topflag = flag[nz - 1][j][i];
      for (k = 0; k < nz; k++) {
        if (topflag & OUTSIDE) {
          flag[k][j][i] |= OUTSIDE;
        } else {
          /* Don't allow cells to be less than MIN_CELL_THICKNESS */
          if (bathy[j][i] < layers[k + 1] && 
	      bathy[j][i] > layers[k + 1] - mincelldz[k]) { 
            double newbotz;
            if ((layers[k + 1] - bathy[j][i]) > 
		(layers[k + 1] - layers[k]) / 2)                
              newbotz = layers[k];
            else
              newbotz = layers[k + 1];
            bathy[j][i] = newbotz;
          }

          if (bathy[j][i] >= layers[k + 1] || bathy[j][i] >= etamax) {
            /* A solid cell */
            flag[k][j][i] |= SOLID;
          }
        }
      }
    }
  }
  /* Reset explicitly defined OUTSIDE cells                          */
  for (is = 1; is <= params->noutside; is++) {
    i = params->oute1[is];
    j = params->oute2[is];
    for (k = 0; k < nz; k++)
      flag[k][j][i] |= OUTSIDE;
  }
}

/* END make_flags()                                                  */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Routine to set the U1SOLID and U1OUTSIDE flags                    */
/*-------------------------------------------------------------------*/
void u1_flags(unsigned long ***flag,  /* Flags array */
              int nce1,         /* Size of the grid in the e1 direction */
              int nce2,         /* Size of the grid in the e2 direction */
              int nz            /* Size of the grid in the z direction */
  )
{
  int i, cl, cr;
  int j;
  int k;
  unsigned long topflag;

  /* Loop over u1 points */
  for (j = 0; j < nce2; j++) {
    for (i = 0; i <= nce1; i++) {
      cr = (i < nce1) ? i : NOTVALID;
      cl = (i > 0) ? i - 1 : NOTVALID;

      /* Check that u1 point has at least 1 neighbouring cell */
      if (cl == NOTVALID && cr == NOTVALID)
        hd_quit("do_u1: u1 point with no cell to left or right : %d %d\n",
                i, j);
      /* Flag boundary points or points outside model */
      if (cl == NOTVALID) {
        /* u1 point only has cell on right */
        if (flag[nz - 1][j][cr] & OUTSIDE)
          /* u1 point not in model */
          topflag = U1OUTSIDE;
        else
          /* Left boundary - solid for now */
          topflag = U1SOLID | L_EDGE;
      } else if (cr == NOTVALID) {
        /* u1 point only has cell on left */
        if (flag[nz - 1][j][cl] & OUTSIDE)
          /* u1 point not in model */
          topflag = U1OUTSIDE;
        else
          /* right boundary - solid for now */
          topflag = U1SOLID | R_EDGE;
      } else {
        /* u1 point has cells on both left and right */
        if ((flag[nz - 1][j][cl] & OUTSIDE) &&
            (flag[nz - 1][j][cr] & OUTSIDE))
          /* u1 point not in model */
          topflag = U1OUTSIDE;
        else if (flag[nz - 1][j][cl] & OUTSIDE)
          /* Left boundary - solid for now */
          topflag = U1SOLID | L_EDGE;
        else if (flag[nz - 1][j][cr] & OUTSIDE)
          /* Right boundary - solid for now */
          topflag = U1SOLID | R_EDGE;
        else
          /* u1 interior point */
          topflag = 0;
      }
      for (k = 0; k < nz; k++) {
        if (topflag != 0)
          /* Exterior boundary or outside */
          flag[k][j][i] |= topflag;
        else if (flag[k][j][cl] & SOLID)
          /* Interior solid wall */
          flag[k][j][i] |= U1SOLID | L_EDGE;
        else if (flag[k][j][cr] & SOLID)
          /* Interior solid wall */
          flag[k][j][i] |= U1SOLID | R_EDGE;
      }
    }
  }

  /* Redefine the masks on the edges */
  for (k = nz - 1; k >= 0; k--) {
    for (j = 0; j <= nce2; j++) {
      i = nce1;
      if (flag[k][j][i - 1] & SOLID)
        flag[k][j][i] |= SOLID;
      /* F1 : Set U1OUTSIDE as OUTSIDE also */
      if (flag[k][j][i - 1] & OUTSIDE)
        flag[k][j][i] |= OUTSIDE;
      if (flag[k][j][i - 1] & U1SOLID)
        flag[k][j][i] |= U1SOLID;
      if (flag[k][j][i - 1] & U2SOLID)
        flag[k][j][i] |= U2SOLID;
    }
  }
}

/* END u1_flags()                                                    */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Routine to set the U2SOLID and U2OUTSIDE flags                    */
/*-------------------------------------------------------------------*/
void u2_flags(unsigned long ***flag,  /* Flags array */
              int nce1,         /* Size of the grid in the e1 direction */
              int nce2,         /* Size of the grid in the e2 direction */
              int nz            /* Size of the grid in the z direction */
  )
{
  int i, cb, cf;
  int j;
  int k;
  unsigned long topflag;

  /* Loop over u2 points */
  for (j = 0; j <= nce2; j++) {
    for (i = 0; i < nce1; i++) {
      cf = (j < nce2) ? j : NOTVALID;
      cb = (j > 0) ? j - 1 : NOTVALID;

      /* Check that u2 point has at least 1 neighbouring cell */
      if (cb == NOTVALID && cf == NOTVALID)
        hd_quit
          ("buildgrid: u2 point with no cell to back or forward : %d %d\n",
           i, j);
      /* Flag boundary points or points outside model */
      if (cb == NOTVALID) {
        /* u2 point only has cell on forward */
        if (flag[nz - 1][cf][i] & OUTSIDE)
          /* u2 point not in model */
          topflag = U2OUTSIDE;
        else
          /* Back boundary - solid for now */
          topflag = U2SOLID | B_EDGE;
      } else if (cf == NOTVALID) {
        /* u2 point only has cell on back */
        if (flag[nz - 1][cb][i] & OUTSIDE)
          /* u2 point not in model */
          topflag = U2OUTSIDE;
        else
          /* forward boundary - solid for now */
          topflag = U2SOLID | F_EDGE;
      } else {
        /* u2 point has cells on both back and forward */
        if ((flag[nz - 1][cb][i] & OUTSIDE) &&
            (flag[nz - 1][cf][i] & OUTSIDE))
          /* u2 point not in model */
          topflag = U2OUTSIDE;
        else if (flag[nz - 1][cb][i] & OUTSIDE)
          /* Back boundary - solid for now */
          topflag = U2SOLID | B_EDGE;
        else if (flag[nz - 1][cf][i] & OUTSIDE)
          /* Forward boundary - solid for now */
          topflag = U2SOLID | F_EDGE;
        else
          /* u2 interior point */
          topflag = 0;
      }
      for (k = 0; k < nz; k++) {
        if (topflag != 0)
          /* Exterior boundary or outside */
          flag[k][j][i] |= topflag;
        else if (flag[k][cb][i] & SOLID)
          /* Interior solid wall */
          flag[k][j][i] |= U2SOLID | B_EDGE;
        else if (flag[k][cf][i] & SOLID)
          /* Interior solid wall */
          flag[k][j][i] |= U2SOLID | F_EDGE;
      }
    }
  }
  /* Redefine the masks on the edges */
  for (k = nz - 1; k >= 0; k--) {
    for (i = 0; i <= nce1; i++) {
      j = nce2;
      if (flag[k][j - 1][i] & SOLID)
        flag[k][j][i] |= SOLID;
      /* F1 : Set U2OUTSIDE as OUTSIDE also */
      if (flag[k][j - 1][i] & OUTSIDE)
        flag[k][j][i] |= OUTSIDE;
      if (flag[k][j - 1][i] & U1SOLID)
        flag[k][j][i] |= U1SOLID;
      if (flag[k][j - 1][i] & U2SOLID)
        flag[k][j][i] |= U2SOLID;
    }
  }
}

/* END u2_flags()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set the corner flags                                   */
/*-------------------------------------------------------------------*/
void corner_flags(unsigned long ***flag,  /* Flags array */
                  int nce1,     /* Size of the grid in the e1 direction */
                  int nce2,     /* Size of the grid in the e2 direction */
                  int nz        /* Size of the grid in the z direction */
  )
{
  int i, j, k;
  int cr, cl, cb, cf;
  int nfe1 = nce1 + 1;
  int nfe2 = nce2 + 1;

  /* Loop over grid corners */
  for (k = 0; k < nz; k++) {
    for (j = 0; j < nfe2; j++) {
      for (i = 0; i < nfe1; i++) {
	cr = (i < nce1) ? i : NOTVALID;
	cl = (i > 0) ? i - 1 : NOTVALID;
        cf = (j < nce2) ? j : NOTVALID;
        cb = (j > 0) ? j - 1 : NOTVALID;

        /* Is this in the middle of the water? */
        if (validcelli(cl) && validcelli(cr) &&
            validcellj(cb) && validcellj(cf) &&
            !(flag[k][cb][cl] & (SOLID | OUTSIDE)) &&
            !(flag[k][cb][cr] & (SOLID | OUTSIDE)) &&
            !(flag[k][cf][cl] & (SOLID | OUTSIDE)) &&
            !(flag[k][cf][cr] & (SOLID | OUTSIDE)))
          flag[k][j][i] |= ALLWATER;

        /* Is this in the middle of a back edge? */
        if (validcelli(cl) && validcelli(cr) &&
            flag[k][j][cl] & U2BDRY && flag[k][j][cl] & B_EDGE &&
            flag[k][j][cr] & U2BDRY && flag[k][j][cr] & B_EDGE)
          flag[k][j][i] |= MID_B_EDGE;

        /* Is this in the middle of a front edge? */
        if (validcelli(cl) && validcelli(cr) &&
            flag[k][j][cl] & U2BDRY && flag[k][j][cl] & F_EDGE &&
            flag[k][j][cr] & U2BDRY && flag[k][j][cr] & F_EDGE)
          flag[k][j][i] |= MID_F_EDGE;

        /* Is this in the middle of a left edge? */
        if (validcellj(cb) && validcellj(cf) &&
            flag[k][cb][i] & U1BDRY && flag[k][cb][i] & L_EDGE &&
            flag[k][cf][i] & U1BDRY && flag[k][cf][i] & L_EDGE)
          flag[k][j][i] |= MID_L_EDGE;

        /* Is this in the middle of a right edge? */
        if (validcellj(cb) && validcellj(cf) &&
            flag[k][cb][i] & U1BDRY && flag[k][cb][i] & R_EDGE &&
            flag[k][cf][i] & U1BDRY && flag[k][cf][i] & R_EDGE)
          flag[k][j][i] |= MID_R_EDGE;
      }
    }
  }

  /* Set cells to solid if required on the end faces */
  j = nfe2 - 1;
  for (k = 0; k < nz; k++) {
    int jj = nce2 - 1;
    for (i = 0; i < nce1; i++) {
      if (flag[k][jj][i] & SOLID || flag[k][j][i] & U2SOLID)
        flag[k][j][i] |= (SOLID | U2SOLID);
      if (flag[k][jj][i] & OUTSIDE || flag[k][j][i] & U2OUTSIDE)
        flag[k][j][i] |= (OUTSIDE | U2OUTSIDE);
    }
  }
  i = nfe1 - 1;
  for (k = 0; k < nz; k++) {
    int ii = nce1 - 1;
    for (j = 0; j < nce2; j++) {
      if (flag[k][j][ii] & SOLID || flag[k][j][i] & U1SOLID)
        flag[k][j][i] |= (SOLID | U1SOLID);
      if (flag[k][j][ii] & OUTSIDE || flag[k][j][i] & U1OUTSIDE)
        flag[k][j][i] |= (OUTSIDE | U1OUTSIDE);
    }
  }
}

/* END corner_flags()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set the maps and flags for explicit mappings           */
/*-------------------------------------------------------------------*/
int *set_explicit_maps(geometry_t *geom, /* Global geometry          */
		       int nem,     /* Number of explicit maps       */
		       int *emis,   /* i location of source map      */
		       int *emjs,   /* j location of source map      */
		       int *emid,   /* i location of destination map */
		       int *emjd,   /* j location of destination map */
		       int *kti,    /* Surface layer for map         */
                       int *kbi,    /* Bottom layer for map          */
		       unsigned long ***flag, /* Flag array          */
		       short **kbot, /* k index of the bottom        */
		       int *ip1,     /* Forward map                  */
		       int *im1,     /* Backward map                 */
		       char *tag,    /* Tag for E1 or E2             */
		       int *emapf    /* Explicit map flag            */
		       )
{
  int i, j, k;
  int cc, cs, cd, cb;
  int is, js, id, jd;
  int kt, kb;
  int nz = geom->nz - 1;
  unsigned long M_SOLID = U1SOLID;     /* Solid flag                 */
  unsigned long M_OUTSIDE = U1OUTSIDE; /* Outside flag               */
  unsigned long M_BEDGE = L_EDGE;      /* Near edge flag             */
  unsigned long M_FEDGE = R_EDGE;      /* Far edge flag              */
  unsigned long D_EDGE = U1AVZERO;     /* Sub-surface map flag       */
  int *sm = NULL;                      /* Sub-surface map array      */
  int *p1, *m1;
  int mode = 0;

  if(strcmp(tag, "E2") == 0) {
    M_SOLID = U2SOLID;
    M_OUTSIDE = U2OUTSIDE;
    M_BEDGE = B_EDGE;
    M_FEDGE = F_EDGE;
    D_EDGE = U2AVZERO;
    mode = 1;
  }

  /* Allocate memory for sub-surface explicit mappings if required   */
  if (nem) {
    sm = i_alloc_1d(geom->sgsizS);
    memset(sm, 0, geom->sgsizS * sizeof(int));
  }

  /* Allocate memory for the zoomed grid inner explicit map mask     */
  if (nem && geom->emmask == NULL) {
    geom->emmask = i_alloc_1d(geom->sgsiz);
    memset(geom->emmask, 0, geom->sgsiz * sizeof(int));
  }

  for (cc = 0; cc < nem; cc++) {

    /* Get the source and destination cells and make sure they are   */
    /* wet.                                                          */
    k = nz;
    is = emis[cc];
    js = emjs[cc];
    id = emid[cc];
    jd = emjd[cc];
    p1 = ip1;
    m1 = im1;

    if (flag[k][js][is] & NOTWET)
      hd_quit("preprocess: line %d MAP_POINTS_%s (%d %d) is not a wet cell\n", 
	      cc, tag, is, js);
    if (flag[k][jd][id] & NOTWET)
      hd_quit("preprocess: line %d MAP_POINTS_%s (%d %d) is not a wet cell\n", 
	      cc, tag, id, jd);
    
    /* Set the maps for cyclic style OBCs                            */
    if (flag[k][js][is] & M_BEDGE && flag[k][jd][id] & M_BEDGE) {
      p1 = im1;
      *emapf |= E1_INNER;
    }
    if (flag[k][js][is] & M_FEDGE && flag[k][jd][id] & M_FEDGE) {
      m1 = ip1;
      *emapf |= E2_INNER;
    }

    /* Read the source and destination cells and make source         */
    /* correspond to a R_EDGE or F_EDGE (i.e. swap indicies if       */
    /* neccessary). Exit with a warning if the cells do not have     */
    /* the correctly oriented solid edges.                           */
    if (flag[k][js][is] & M_BEDGE) {
      i = id; j = jd;
      id = is; jd = js;
      is = i; js = j;
    } else if (!(flag[k][jd][id] & M_BEDGE)) {
      if(strcmp(tag, "E1") == 0)
	hd_quit("preprocess: MAP_POINTS_%s line %d (%d %d) must contain a cell with solid left edge.\n", tag, cc, id, jd);
      else
	hd_quit("preprocess: MAP_POINTS_%s line %d (%d %d) must contain a cell with solid back edge.\n", tag, cc, id, jd);
    }

    /* Reset the flags                                               */
    kb = max(kbot[jd][id], kbot[js][is]);
    kb = max(kbi[cc], kb);
    kt = max(kti[cc], kb);
    for (k = nz; k >= kb; k--) {
      if (flag[k][js][is] & M_BEDGE && flag[k][jd][id] & M_BEDGE &&
	  *emapf & (E1_INNER|E2_INNER)) {
	if (flag[k][js][is] & M_SOLID)
	  flag[k][js][is] &= ~M_SOLID;
	if (flag[k][js][is] & M_OUTSIDE)
	  flag[k][js][is] &= ~M_OUTSIDE;
      }
      if (flag[k][jd][id] & M_SOLID)
	flag[k][jd][id] &= ~M_SOLID;
      if (flag[k][jd][id] & M_OUTSIDE)
	flag[k][jd][id] &= ~M_OUTSIDE;
    }

    /* Map the cell centers. Note; although the explicit maps are    */
    /* set from the surface (nz) to kb, when the code executes the   */
    /* cells to process for velocity are only set from sm_e1 / sm_e2 */
    /* to kb, hence cells above sm_e1/sm_e2 remain at zero velocity. */
    /* This also means no tracer transport occurs above sm_e1/sm_e2  */
    /* (since velocity is zero).                                     */
    cs = geom->map[nz][js][is];
    cd = geom->map[nz][jd][id];
    cb = geom->map[kb][js][is];
    p1[cs] = cd;
    m1[cd] = cs;

    if (flag[nz][js][is] & M_BEDGE && flag[nz][jd][id] & M_BEDGE) {
      geom->emmask[cs] = cd;
      geom->emmask[cd] = cs;
    }
    if (flag[nz][js][is] & M_FEDGE && flag[nz][jd][id] & M_FEDGE) {
      geom->emmask[cs] = cd;
      geom->emmask[cd] = cs;
    }

    /* Map sub-surface cells                                         */
    cs = geom->zm1[cs];
    cd = geom->zm1[cd];
    k = nz;
    while (cs <= cb) {
      if (flag[k][js][is] & M_BEDGE && flag[k][jd][id] & M_BEDGE) {
	geom->emmask[cs] = cd;
	geom->emmask[cd] = cs;
      }
      if (flag[k][js][is] & M_FEDGE && flag[k][jd][id] & M_FEDGE) {
	geom->emmask[cs] = cd;
	geom->emmask[cd] = cs;
      }
      p1[cs] = cd;
      m1[cd] = cs;
      cs = geom->zm1[cs];
      cd = geom->zm1[cd];
      k--;	
    }

    /* Find the first sub-surface cell for explicit mapping          */
    if (kti[cc] < nz) {
      cd = geom->map[nz][jd][id];
      sm[cd] = geom->map[kt][jd][id];
    }
    for (k = kt; k >= kb; k--) {
      flag[k][jd][id] |= D_EDGE;
    }
  }
  return(sm);
}

/* END set_explicit_maps()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to allocate memory from the geometry data structure       */
/* arrays.                                                           */
/*-------------------------------------------------------------------*/
void alloc_geom(geometry_t *geom, /* Model geometry structure */
                unsigned long mask  /* Mask for which arrays to allocate */
  )
{
  int size;                     /* Size of the 3D sparse array */
  int sizeS;                    /* Size of the 2D sparse array */

  size = geom->enon + 1;
  sizeS = geom->enonS + 1;

  /* Mapping arrays */
  if (mask & MAP_A) {
    geom->xp1 = i_alloc_1d(size);
    geom->yp1 = i_alloc_1d(size);
    geom->zp1 = i_alloc_1d(size);
    geom->xm1 = i_alloc_1d(size);
    geom->ym1 = i_alloc_1d(size);
    geom->zm1 = i_alloc_1d(size);
    /* 
       geom->xm12 = i_alloc_1d(size); geom->ym12 = i_alloc_1d(size);
       geom->zm12 = i_alloc_1d(size); */
  }
  if (mask & GRID_A) {
    /* 3D arrays */
    geom->gridz = d_alloc_1d(size);
    geom->cellz = d_alloc_1d(size);
    /* 2D arrays */
    geom->h1acell = d_alloc_1d(sizeS);
    geom->h2acell = d_alloc_1d(sizeS);
    geom->h1au1 = d_alloc_1d(sizeS);
    geom->h2au2 = d_alloc_1d(sizeS);
    geom->h2au1 = d_alloc_1d(sizeS);
    geom->h1au2 = d_alloc_1d(sizeS);
    geom->cellarea = d_alloc_1d(sizeS);
    geom->dHde1 = d_alloc_1d(sizeS);
    geom->dHde2 = d_alloc_1d(sizeS);
    geom->botz = d_alloc_1d(sizeS);
    geom->botzu1 = d_alloc_1d(sizeS);
    geom->botzu2 = d_alloc_1d(sizeS);
    geom->botzgrid = d_alloc_1d(sizeS);
    geom->cellx = d_alloc_1d(sizeS);
    geom->celly = d_alloc_1d(sizeS);
    geom->gridx = d_alloc_1d(sizeS);
    geom->gridy = d_alloc_1d(sizeS);
    geom->u1x = d_alloc_1d(sizeS);
    geom->u1y = d_alloc_1d(sizeS);
    geom->u2x = d_alloc_1d(sizeS);
    geom->u2y = d_alloc_1d(sizeS);

    if (geom->sednz) {
      geom->gridz_sed = d_alloc_2d(sizeS, geom->sednz + 1);
      geom->cellz_sed = d_alloc_2d(sizeS, geom->sednz);
    }
  }
  if (mask & MASTER_A) {
    geom->s2i = i_alloc_1d(size);
    geom->s2j = i_alloc_1d(size);
    geom->s2k = i_alloc_1d(size);
    geom->thetau1 = d_alloc_1d(sizeS);
    geom->sinthcell = d_alloc_1d(sizeS);
    geom->costhcell = d_alloc_1d(sizeS);
    geom->sinthu1 = d_alloc_1d(sizeS);
    geom->costhu1 = d_alloc_1d(sizeS);
    geom->thetau2 = d_alloc_1d(sizeS);
    geom->sinthu2 = d_alloc_1d(sizeS);
    geom->costhu2 = d_alloc_1d(sizeS);
  }
  if (mask & WINDOW_A) {
    geom->wsa = i_alloc_1d(size);
    geom->xmyp1 = i_alloc_1d(size);
    geom->xpym1 = i_alloc_1d(size);
  }
}

/* END alloc_geom()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to check that a sparse coordinate is in the valid range,  */
/* and exit if not.                                                  */
/*-------------------------------------------------------------------*/
int isvalidc(int c, int sz, char *warn) {

  if (c <= 0 || c > sz)
    hd_quit("invalid sparse coordinate located in: %s\n", warn);

  return(c);
}

int iswetijk(int c, geometry_t *geom, unsigned long ***flag) {
  int i = geom->s2i[c];
  int j = geom->s2j[c];
  int k = geom->s2k[c];
  if (geom->compatible & V1957) return 1;
  if (i == NOTVALID || j == NOTVALID || k == NOTVALID)
    return 0;
  else
    return 1;
}

/* END isvalidc()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to point the common arrays in a window geometry structure */
/* to the corresponding array in the master geometery.               */
/*-------------------------------------------------------------------*/
void point_geom(geometry_t *window, /* Window geometry structure */
                geometry_t *geom  /* Model geometry structure */
  )
{
  if (geom->enon != window->enon)
    hd_quit
      ("point_geom: Master and window 3D sparse array sizes incompatible\n");

  if (geom->enonS != window->enonS)
    hd_quit
      ("point_geom: Master and window 2D sparse array sizes incompatible\n");

  /* Point the window maps at the master */
  window->xp1 = geom->xp1;
  window->xm1 = geom->xm1;
  window->yp1 = geom->yp1;
  window->ym1 = geom->ym1;
  window->zp1 = geom->zp1;
  window->zm1 = geom->zm1;

  /* Point the window geometric arrays at the master */
  /* 3D arrays */
  window->gridz = geom->gridz;
  window->cellz = geom->cellz;

  if (window->sednz) {
    window->gridz_sed = geom->gridz_sed;
    window->cellz_sed = geom->cellz_sed;
  }

  /* 2D arrays */
  window->h1acell = geom->h1acell;
  window->h2acell = geom->h2acell;
  window->h1au1 = geom->h1au1;
  window->h2au2 = geom->h2au2;
  window->h2au1 = geom->h2au1;
  window->h1au2 = geom->h1au2;
  window->cellarea = geom->cellarea;
  window->dHde1 = geom->dHde1;
  window->dHde2 = geom->dHde2;
  window->botz = geom->botz;
  window->botzu1 = geom->botzu1;
  window->botzu2 = geom->botzu2;
  window->botzgrid = geom->botzgrid;
  window->cellx = geom->cellx;
  window->celly = geom->celly;
  window->gridx = geom->gridx;
  window->gridy = geom->gridy;
  window->u1x = geom->u1x;
  window->u1y = geom->u1y;
  window->u2x = geom->u2x;
  window->u2y = geom->u2y;
  window->sinthcell = geom->sinthcell;
  window->costhcell = geom->costhcell;
}

/* END point_geom()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to reset the flags array for sigma                        */
/*-------------------------------------------------------------------*/
void sigma_flags(int nfe1, int nfe2, int nz, unsigned long ***flg)
{
  int i, j, k;

  /* Set the flag as wet through the entire water column if the */
  /* surface layer is wet.  */
  for (j = 0; j < nfe2; j++)
    for (i = 0; i < nfe1; i++) {
      if (!(flg[nz][j][i] & SOLID)) {
        for (k = nz; k >= 0; k--)
          if (flg[k][j][i] & SOLID)
            flg[k][j][i] &= ~(SOLID);
      }
      if (!(flg[nz][j][i] & U1SOLID)) {
        for (k = nz; k >= 0; k--)
          if (flg[k][j][i] & U1SOLID)
            flg[k][j][i] &= ~(U1SOLID);
      }
      if (!(flg[nz][j][i] & U2SOLID)) {
        for (k = nz; k >= 0; k--)
          if (flg[k][j][i] & U2SOLID)
            flg[k][j][i] &= ~(U2SOLID);
      }

      if (flg[nz][j][i] & U2BDRY)
        for (k = nz; k >= 0; k--)
          if (!(flg[k][j][i] & U2BDRY))
            flg[k][j][i] |= U2BDRY;
      if (flg[nz][j][i] & F_EDGE)
        for (k = nz; k >= 0; k--)
          if (!(flg[k][j][i] & F_EDGE))
            flg[k][j][i] |= F_EDGE;
      if (flg[nz][j][i] & B_EDGE)
        for (k = nz; k >= 0; k--)
          if (!(flg[k][j][i] & B_EDGE))
            flg[k][j][i] |= B_EDGE;
      if (!(flg[nz][j][i] & F_EDGE))
        for (k = nz; k >= 0; k--)
          if (flg[k][j][i] & F_EDGE)
            flg[k][j][i] &= ~(F_EDGE);
      if (!(flg[nz][j][i] & B_EDGE))
        for (k = nz; k >= 0; k--)
          if (flg[k][j][i] & B_EDGE)
            flg[k][j][i] &= ~(B_EDGE);
      if (!(flg[nz][j][i] & MID_B_EDGE))
        for (k = nz; k >= 0; k--)
          if (flg[k][j][i] & MID_B_EDGE)
            flg[k][j][i] &= ~(MID_B_EDGE);
      if (!(flg[nz][j][i] & MID_F_EDGE))
        for (k = nz; k >= 0; k--)
          if (flg[k][j][i] & MID_F_EDGE)
            flg[k][j][i] &= ~(MID_F_EDGE);

      if (flg[nz][j][i] & U1BDRY)
        for (k = nz; k >= 0; k--)
          if (!(flg[k][j][i] & U1BDRY))
            flg[k][j][i] |= U1BDRY;
      if (flg[nz][j][i] & L_EDGE)
        for (k = nz; k >= 0; k--)
          if (!(flg[k][j][i] & L_EDGE))
            flg[k][j][i] |= L_EDGE;
      if (flg[nz][j][i] & R_EDGE)
        for (k = nz; k >= 0; k--)
          if (!(flg[k][j][i] & R_EDGE))
            flg[k][j][i] |= R_EDGE;
      if (!(flg[nz][j][i] & L_EDGE))
        for (k = nz; k >= 0; k--)
          if (flg[k][j][i] & L_EDGE)
            flg[k][j][i] &= ~(L_EDGE);
      if (!(flg[nz][j][i] & R_EDGE))
        for (k = nz; k >= 0; k--)
          if (flg[k][j][i] & R_EDGE)
            flg[k][j][i] &= ~(R_EDGE);
      if (!(flg[nz][j][i] & MID_R_EDGE))
        for (k = nz; k >= 0; k--)
          if (flg[k][j][i] & MID_R_EDGE)
            flg[k][j][i] &= ~(MID_R_EDGE);
      if (!(flg[nz][j][i] & MID_L_EDGE))
        for (k = nz; k >= 0; k--)
          if (flg[k][j][i] & MID_L_EDGE)
            flg[k][j][i] &= ~(MID_L_EDGE);

      if (flg[nz][j][i] & ALLWATER) {
        for (k = nz; k >= 0; k--)
          if (!(flg[k][j][i] & ALLWATER))
            flg[k][j][i] |= ALLWATER;
      }
    }
}

/* END sigma_flags()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Removes dry cells from a points list                              */
/*-------------------------------------------------------------------*/
void reset_points(parameters_t *params, unsigned long **flag, 
		  int *nlist, int *listi, int *listj) {
  int n, nn, i, j;
  int *ci, *cj;

  nn = 0;

  for (n = 1; n <= *nlist; n++) {
    i = listi[n];
    j = listj[n];
    if (!(flag[j][i] & (SOLID|OUTSIDE))) nn++;
  }
  if (nn) {
    ci = i_alloc_1d(nn +1);
    cj = i_alloc_1d(nn +1);
    nn = 1;
    for (n = 1; n <= *nlist; n++) {
      i = listi[n];
      j = listj[n];
      if (!(flag[j][i] & (SOLID|OUTSIDE))) {
	ci[nn] = i;
	cj[nn] = j;
	nn++;
      }
    }
    *nlist = nn - 1;
    for (n = 1; n <= *nlist; n++) {
      listi[n] = ci[n];
      listj[n] = cj[n];
    }
    i_free_1d(ci);
    i_free_1d(cj);
  }
}

/* END reset_points()                                                */
/*-------------------------------------------------------------------*/

