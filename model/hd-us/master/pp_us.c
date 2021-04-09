/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/master/pp_us.c
 *  
 *  Description: Sparse preprocessor
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: pp_us.c 6736 2021-03-30 00:39:49Z her127 $
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

int nz;          : Size of the grid in the vertical direction.

double *bathy;  : Cell centered bathymetry array of size bathy[ns2].
All bathymety values should be negative. The value 0 is assumed to be
mean sea level. A cell is identified as a land cell if its value > 10. 

double *layers;  : Layer structure of size layers[nz];
Layers should be an array of size nz, where layers[nz] = maximum depth and
layers[0] =0.

A mesh_t structure *m that contains the coordinate information (xloc, yloc),
a list of integers pointing to these coordinate index locations to define the 
  polygon edges and centres, and open boundary information:

int nobc;        : Number of open boundaries.
int npts[nobc];  : Number of surface cells in each open boundary.
int loc[nobc][npts]; : List of locations for each open boundary.
int obc[0:1][npts][nobc]; : Mapping function containing the coordinate
                            indes locations that define the edge of the OBC.

e.g. the following domain has three open boundaries, two cross-shelf
U1BDRY and one along-shore U2BDRY.

"b3 int ns2 = 12;
 int nz = 11;
 double bathy[12] =   {-60, -70, -80, -90,   <- wet cells
                       -55, -65, -75, -85,   <- wet cells
                        99,  99,  99,  99};  <- land cells    
 double layers[11] = {-100, -90, -80, -70, -60, -50, -40, -30, -20, -10, 0};
 int nobc = 3;
 int npts[3] = {2, 2, 4};
 int iloc[3][4] = {{0, 0},         <- boundary 0
                   {3, 3},         <- boundary 1
                   {0, 1, 2, 3}};  <- boundary 2
 int jloc[3][4] = {{1, 2},         <- boundary 0
                   {1, 2},         <- boundary 1
                   {2, 2, 2, 2}};  <- boundary 2


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

/* Codes for vertices */
#define V_W  1
#define V_B  2
#define V_G  4
#define V_A  V_W|V_B|V_G

void make_flags_us(parameters_t *params, unsigned long **flag, double *bathy,
		   double *layers, int ns2, int nz);
int is_obc(int c, int nobc, int *npts, int **loc);
int obc_num(int c, int nobc, int *npts, int **loc);
void obc_num_all(int c, int nobc, int *npts, int **loc, int *mask);
int is_obce(int cc, int j, mesh_t *m);
int is_obceo(int npe, int cc, int j, double **x, double **y, double ***posx, double ***posy,
             int nobc, int *npts);
void write_us_map(geometry_t *sgrid, parameters_t *params);
void reorder_edges(geometry_t *sgrid);
int cyc_m2(geometry_t *sgrid, int *nmap, int *omap, int c);
int oedge(int npe, int n);
void swap_edge(geometry_t *sgrid, int e);
int check_vert2d(geometry_t *sgrid, int e, int e1, int i1, int i2,
	       double **locx, double **locy, int **emap, int *maske, int *b1, int *b2, int *b3);
int check_vert3d(geometry_t *sgrid, int vs, int vv, int ed, int eb, int *c1, int *c2, int *c3);
void find_vertices(geometry_t *sgrid, double **locx, double **locy, int *maske, int **emap);
int find_vertex(int c, double x, double y, double **xloc, double **yloc, int *mask, int ns2, int npe);
int find_edges_tan(geometry_t *sgrid, mesh_t *m, int k, int *e, int *ee, int *c4, int *b4,
		   int *maske, int **vle, point *tegl);
int get_nve(int npe);
void vertex_map_4(geometry_t *sgrid);
void create_delaunay_cell(geometry_t *geom, parameters_t *params);
void create_delaunay_cent(geometry_t *geom, parameters_t *params);
void build_delaunay_cell(geometry_t *geom, parameters_t *params, point *tegl);
int get_limit_obc(parameters_t *params, int ns2, int **neic);
int e2eo(geometry_t *sgrid, int e, int j);
int get_bind(geometry_t *sgrid, int c, int *mask);

int usejo = 0; /* Method of finding opposite edge of a polygon */

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
void build_sparse_grid_us(parameters_t *params,
			  geometry_t *sgrid,
			  int nz,         /* Grid size, z direction  */
			  double *bathy,  /* Input bathymetry array  */
			  double *layers  /* Layer structure         */
			  )
{
  int ns2;           /* Grid size                                    */
  int e, v, c = 1;   /* Locations of wet points in the sparse array  */
  int cc, ee, vv;    /* Counter                                      */
  int cs, es, vs;    /* Surface cell centre / edge                   */
  int cb, eb;        /* Boundary/bottom cell centre/edge             */
  int gc;            /* Locations of ghost points in sparse array    */
  int cn, cns;       /* Sparse neighbour location                    */
  int lc;            /* Loc of boundary points in the boundary array */
  int i, j, k;       /* x,y,z counters                               */
  int ii, jj;        /* x,y counters                                 */
  int n, tn;         /* General counters                             */
  int c1, c2, c3, c4;/* Counters                                     */
  int b1, b2, b3, b4;/* Counters                                     */
  int ei, e1, e2;    /* Counters                                     */
  int v1, v2;        /* Counters                                     */
  int cp1;           /* Spatial map in positive direction            */
  int cm1;           /* Spatial map in negative direction            */
  int *end_wet;      /* End loaction of the wet points in each layer */
  int *num_gst;      /* Number of ghost points in each layer         */
  int *num_bdy;      /* Number of OBC points in each layer           */
  int num_mc = 0;    /* Number of multiple ghost cells               */
  int num_wc = 0;    /* Number of 3D wet grid cells                  */
  int num_gc = 0;    /* Number of ghost grid cells                   */
  int num_sc = 0;    /* Number of sediment grid cells                */
  int num_scg = 0;   /* Number of sediment ghost grid cells          */
  int num_gc2D = 0;  /* Number of 2D ghost grid cells                */
  int num_obc = 0;   /* Number of OBC ghost cells                    */
  int num_obc2D = 0; /* Number of 2D OBC ghost cells                 */
  int num_vert;      /* Number of vertices                           */
  short *kbot;       /* k index of the bottom                        */
  long scen;         /* Flag for ghost cell type                     */
  int gchck = 0;     /* Checks locations of ghost points             */
  int npe;           /* Number of nodes                              */
  int npem;          /* Maximum number of nodes                      */
  int nve;           /* Number of vertices                           */
  int *mask;         /* Centre mask                                  */
  int *maske;        /* Edge mask                                    */
  int *maskv;        /* Vertex mask                                  */
  int *maskb;        /* Boundary mask                                */
  int **emap;        /* Edge ghost mask                              */
  int *imapc;        /* Interior centre map                          */
  int *omapc;        /* Outward centre map                           */
  int *imape;        /* Interior edge map                            */
  point *tegl;       /* Tangential ghost edge location               */
  int **vlm;         /* Vertical location centre map                 */
  int **vle;         /* Vertical location edge map                   */
  unsigned long **flg; /* Flag for unstructured grid                 */
  unsigned long ***flag;/* Flag for Cartesian grid                   */
  double **topo;     /* Cartesian bathymetry                         */
  double bmax;       /* Maximum depth                                */
  int ewet;          /* Number of 3D wet cells                       */
  int ewetS;         /* Number of 2D wet cells                       */
  int sigma = 0;     /* Set to 1 for sigma model                     */
  int laus = 2;      /* Ghost cells adjacent to OBCs                 */
  int rtype = 2;     /* Type of Thuburn (2009) weights               */
  int *cc2s;         /* Input list to sparse array map               */
  int *e2ee;         /* Edge to index map                            */
  int *c2cc;         /* Cell to index map                            */
  mesh_t *m;         /* Input mesh                                   */
  int **neic;
  int **neij;
  int *dume;
  double *len, **rw;
  double **locx, **locy;
  double *olayers;
  int **v2e, **v2c;
  int dof = 0;
  double *xloc, *yloc;
  int ***eloc;
  int nvert;
  int tage = 1;      /* Include tangential ghost edges               */
  int pc = 0;        /* Print cell centre information                */
  int pe = 0;        /* Print cell edge information                  */
  int printcell = 0; /* Print cell number coordinates                */
  double printlocx = NOTVALID;  /* Cell x location for debugging     */
  double printlocy = NOTVALID;  /* Cell y location for debugging     */
  int newcode = 1;
  int newvert = 1;
  int *wgsts;

  if (DEBUG("init_m"))
    dlog("init_m", "\nStart unstructured preprocessor\n");

  m = params->mesh;
  ns2 = m->ns2;
  topo = NULL;
  flag = NULL;
  olayers = d_alloc_1d(params->nz + 1);
  memcpy(olayers, params->layers, (params->nz + 1) * sizeof(double));

  /* Set OUTSIDE cells                                               */
  for (c = 1; c <= ns2; c++) {
    if (isnan(bathy[c])) bathy[c] = NOTVALID;
  }

  /* Reset flags for sigma calculations if required                  */
  if (sigma) {
    /* Find the maximum depth                                        */
    bmax = 1e10;
    for (c = 1; n <= ns2; c++)
	if (bathy[c] < bmax)
	  bmax = bathy[c];
    for (n = 0; n < nz; n++)
      layers[n] *= fabs(bmax);
  }

  /* Allocate memory                                                 */
  flg = (unsigned long **)l_alloc_2d(ns2+1, nz);

  /* Make the flags array */
  make_flags_us(params, flg, bathy, layers, ns2, nz);

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */
  sgrid->us_type = params->us_type;
  end_wet = i_alloc_1d(nz + 1);
  num_gst = i_alloc_1d(nz + 1);
  num_bdy = i_alloc_1d(nz + 1);
  sgrid->nz = nz;
  sgrid->sednz = params->sednz;
  npe = m->mnpe;
  npem = sgrid->npem = m->mnpe;
  sgrid->nvem = get_nve(npem);
  sgrid->neem = 2 * sgrid->npem - 1;
  sgrid->nobc = m->nobc;

  for (k = 0; k <= nz; k++)
    end_wet[k] = num_gst[k] = num_bdy[k] = 0;
  sgrid->nwindows = params->nwindows;
  if (params->win_size) {
    sgrid->win_size = d_alloc_1d(sgrid->nwindows + 1);
    for (n = 1; n <= sgrid->nwindows; n++) {
      sgrid->win_size[n] = params->win_size[n - 1];
    }
  }
  if (params->nwn) {
    sgrid->nwn = params->nwn;
    sgrid->wnx = params->wnx;  
    sgrid->wny = params->wny;
  }
  sgrid->layers = d_alloc_1d(nz + 1);
  memcpy(sgrid->layers, olayers, (params->nz + 1) * sizeof(double));
  d_free_1d(olayers);
  sgrid->compatible = params->compatible;
 
  /*-----------------------------------------------------------------*/
  /* Set up the array of the bottom coordinate                       */
  kbot = s_alloc_1d(ns2 + 1);
  /* Set the bottom coordinate in the interior                       */
  for (c = 1; c <= ns2; c++) {
    if (flg[nz - 1][c] & NOTWET) {
      kbot[c] = -1;
    } else {
      /* Check bottom is inside model                                */
      if (bathy[c] < layers[0]) {
	hd_quit("build_sparse_grid() : Bottom lower than lowest model level\n");
      }
      /* Check bottom and top below MAXGRIDZ                         */
      if (bathy[c] >= MAXGRIDZ) {
	hd_quit("build_sparse_grid() : Bottom above MAXGRIDZ\n");
      }
      
      /* Loop vertically to find bottom and top                      */
      kbot[c] = -1;
      for (k = 0; k < nz; k++) {
	if (bathy[c] >= layers[k] && bathy[c] < layers[k + 1])
	  kbot[c] = k;
      }
      if (kbot[c] < 0) {
	hd_quit("build_sparse_grid() : bottom not found at (%d) : %5.1f m\n",
		c, bathy[c]);
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Map the cells in the entire grid that may potentially be wet.   */
  /* The sparse grid at this stage consists of wet cells stored in   */
  /* consecutive order, i.e. a 'wet' sparse map is created. The      */
  /* number of wet cells in each layer is stored in end_wet[]. This  */
  /* is used later to insert the ghost cells for each layer into     */
  /* locations following the wet cells in each layer.                */
  /* The total number of wet cells is stored in num_wc here.         */
  /*
  sgrid->ewetS = 0;
  for (k = nz - 1; k >= 0; k--) {
    for (cc = 1; cc <= ns2; cc++) {
      if (k < kbot[cc]) continue;
      num_wc++;
      if (k == nz - 1)
	sgrid->ewetS++;
      if (is_obc(cc, m->nobc, m->npts, m->loc))
	num_bdy[k]++;
      else
	end_wet[k]++;
    }
  }
  */
  /*-----------------------------------------------------------------*/
  /* If no wet points are found then issue a warning and quit.       */
  /*
  if (num_wc == 0) {
    hd_quit("No wet points located in the domain.\n");
  }
  */

  /*-----------------------------------------------------------------*/
  /* Find the neighbour cells to each cell centre. For cell cc, this */
  /* is defined the cell centre that shares a common edge with cc.   */
  /* This is done using find_neighbour_l() which assumes a map       */
  /* eloc[c][j][] exits that relates vertices in a list of locations */
  /* (xloc, yloc) for each cell c and edge. If for and edges i & j:  */ 
  /* eloc[0][c1][j]==eloc[0][c2][i]&&eloc[1][c1][j]==eloc[1][c2][i]||*/
  /* eloc[0][c1][j]==eloc[1][c2][i]&&eloc[1][c1][j]==eloc[0][c2][i]  */
  /* then c1 and c2 are neighbours. This method allows edges to be   */
  /* listed with different sense and start vertex.                   */
  /* Currently the eloc map is created from params->x and params->y  */
  /* in convert_mesh_input() and stored in params->mesh.             */
  /* Alternatively, find_neighbour() can be used with params->x and  */
  /* params->y. Here the (x,y) location are compared and             */
  /* x[c1][j]==x[c2][i] && x[c1][j+1]==x[c2][i-1] then c1 and c2 are */
  /* neighbours. This assumes that all cells have their vertices     */
  /* listed in the same sense (clockwise or anti-clockwise) and      */
  /* start from the same vertex.                                     */
  xloc = m->xloc;
  yloc = m->yloc;
  eloc = m->eloc;
  if (newcode) {
    neic = m->neic;
    neij = m->neij;

  } else {
  neic = i_alloc_2d(ns2+1, npem+1);
  neij = i_alloc_2d(ns2+1, npem+1);
  for (cc = 1; cc <= ns2; cc++) {
    npe = m->npe[cc];
    for (j = 1; j <= npe; j++) {
      jj = j;
      /*cn = find_neighbour(cc, params->x, params->y, npe, ns2, &jj);*/
      /*cn = find_neighbour_l(cc, eloc, npe, ns2, &jj);*/
      cn = find_neighbour_l(cc, eloc, m->npe, ns2, &jj);
      neic[j][cc] = cn;
      neij[j][cc] = jj;
    }
  }
  }
  get_mesh_obc(params, neic);

  /* Old code using (double) boundary locations in params->posx
     and params->posy.
  if (get_limit_obc(params, ns2, neic)) {
    nobc = params->nbu;
    npts = params->nptsu;
    loc = params->locu;
  }
  */
  if (DEBUG("init_m"))
    dlog("init_m", "\nMappings generated OK\n");

  /*-----------------------------------------------------------------*/
  /* Map the cells in the entire grid that may potentially be wet.   */
  /* The sparse grid at this stage consists of wet cells stored in   */
  /* consecutive order, i.e. a 'wet' sparse map is created. The      */
  /* number of wet cells in each layer is stored in end_wet[]. This   */
  /* is used later to insert the ghost cells for each layer into     */
  /* locations following the wet cells in each layer.                */
  /* The total number of wet cells is stored in num_wc here.         */
  sgrid->ewetS = 0;
  for (k = nz - 1; k >= 0; k--) {
    for (cc = 1; cc <= ns2; cc++) {
      if (k < kbot[cc]) continue;
      num_wc++;
      if (k == nz - 1)
	sgrid->ewetS++;
      if (is_obc(cc, m->nobc, m->npts, m->loc)) {
	num_bdy[k]++;
      }
      else
	end_wet[k]++;
    }
  }

  /*-----------------------------------------------------------------*/
  /* If no wet points are found then issue a warning and quit.       */
  if (num_wc == 0) {
    hd_quit("No wet points located in the domain.\n");
  }
  if (DEBUG("init_m"))
    dlog("init_m", "\nWet cell centres generated OK\n");

  /*-----------------------------------------------------------------*/
  /* Find the locations of the ghost cells. These are the first land */
  /* (dry) cells adjacent to a wet cell in the x, y or z direction.  */
  /* The ghost cells in the z direction constitute the sediment      */
  /* layer. Note : in certain instances (e.g. outside corners) there */
  /* may exist more than one wet cell associated with the same ghost */
  /* cell -> multiple ghost cells). First count the total number of  */
  /* ghost cells in the grid and store the number of ghost cells in  */
  /* each layer in num_gst[]. Open boundary cells are included so    */
  /* the land cell adjacent to the start/end of the OBC is included  */
  /* as a ghost cell. The OBC ghosts are removed from the ghost      */
  /* vector later.                                                   */
  sgrid->nbe2 = 0;
  for (k = nz - 1; k >= 0; k--) {
    num_gst[k] = 0;
    for (cc = 1; cc <= ns2; cc++) {
      /* ERROR : should be kbot[cc] below, but this causes seg fault */
      /* in some instances, indicating not all ghost cells are       */
      /* accounted for. - fixed 2/3/2017.                            */
      if (k < kbot[cc]) continue;
      npe = m->npe[cc];
      for (j = 1; j <= npe; j++) {
	cn = neic[j][cc];
	jj = neij[j][cc];
	/*if (!cn && !is_obc(cc, m->nobc, m->npts, m->loc)) {*/
	if (!cn) {
	  num_gst[k]++;
	} else {
	  /* Neighbour cn is shallower than cell cc                  */
	  if (k < kbot[cn] && kbot[cc] < kbot[cn]) {
	    num_gst[k]++;
	  }
	}
      }
    }
    num_gc += num_gst[k];
  }
  sgrid->nbpt = num_gc;
  sgrid->nbptS = num_gc2D = num_gst[nz - 1];
  if (DEBUG("init_m"))
    dlog("init_m", "\nGhost cell centres generated OK\n");

  /*-----------------------------------------------------------------*/
  /* Calculate the number of sediment cells and add to the total     */
  /* number of ghost cells. Note: the number of sediment cells is    */
  /* equivalent to the number of potentially wet cells in the        */
  /* surface layer.                                                  */
  /* First add the number of ghost cells to sedimet cells (this is   */
  /* the same as the number of ghost cells for the surface layer).   */
  num_sc = ns2;
  num_gc += num_sc;

  /*-----------------------------------------------------------------*/
  /* Get the number of OBC ghost cells. There are laus cells for     */
  /* open boundary point so that the advection scheme may be used on */
  /* the boundary. The first of these has been included in the loop  */
  /* above.                                                          */
  for (k = nz - 1; k >= 0; k--) {
    for (cc = 1; cc <= ns2; cc++) {
      if (k < kbot[cc]) continue;
      npe = m->npe[cc];
      for (j = 1; j <= npe; j++) {
	cn = neic[j][cc];
	jj = neij[j][cc];
	if (!cn && is_obc(cc, m->nobc, m->npts, m->loc)) {
	  num_gst[k] += (laus - 1);
	  num_gc += (laus - 1);
	  num_obc += 1;
	  if (k == nz - 1) {
	    num_gc2D += (laus - 1);
	    num_obc2D += 1;
	  }
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
    c = end_wet[k + 1] + num_bdy[k + 1] + num_gst[k + 1];
    end_wet[k] += c;
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
  sgrid->szc = sgrid->sgsiz;
  sgrid->szcS = sgrid->sgsizS;
  if (DEBUG("init_m"))
    dlog("init_m", "\nCell centres counted OK\n");

  /*-----------------------------------------------------------------*/
  /* Allocate memory for the boundary vectors. The vector bpt        */
  /* contains the sparse locations of all ghost cells (land boundary */
  /* cells) laterally adjacent to wet cells. The vector bin contains */
  /* the sparse locations of the first wet cell adjacent to the land */
  /* boundary.                                                       */
  sgrid->bpt = i_alloc_1d(sgrid->nbpt + 1);
  sgrid->bin = i_alloc_1d(sgrid->nbpt + 1);
  sgrid->dbin = i_alloc_1d(sgrid->nbpt + 1);
  sgrid->dbpt = i_alloc_1d(sgrid->nbpt + 1);
  sgrid->wgst = i_alloc_1d(sgrid->szc);
  sgrid->m2d = i_alloc_1d(sgrid->szc);
  sgrid->c2cc = i_alloc_1d(sgrid->szc);
  sgrid->zp1 = i_alloc_1d(sgrid->szc);
  sgrid->zm1 = i_alloc_1d(sgrid->szc);
  memset(sgrid->wgst, 0, sgrid->szc * sizeof(int));
  memset(sgrid->bpt, 0, (sgrid->nbpt+1) * sizeof(int));
  memset(sgrid->bin, 0, (sgrid->nbpt+1) * sizeof(int));

  /*
  n = end_wet[0] + num_bdy[0] + num_gst[0];
  sgrid->w3_t = i_alloc_1d(n + 1);
  n = end_wet[nz - 1] + num_bdy[nz - 1] + num_gst[nz - 1];
  sgrid->w2_t = i_alloc_1d(n + 1);
  */

  sgrid->v2_t = 1;
  sgrid->b2_t = end_wet[nz - 1];
  sgrid->n2_t = end_wet[nz - 1] + num_bdy[nz - 1];

  /*-----------------------------------------------------------------*/
  /* Get the mapping from the input list to sparse array.            */
  /* Cells are categorized according to:                             */
  /* Wet cell: surrounded by npe cells.                              */
  /* OBC cell: surrounded by < npe cells and listed in vector loc.   */
  /* Land boundary cell: surrounded by < npe cells and not  listed   */
  /* in loc.                                                         */
  /* Get the index to sparse map, cc2s, where cc ranges from 1:ns2.  */
  /* Also get the vertical layer map, vlm[c][k] where c ranges from  */
  /* 1:ns2.                                                          */
  /* Ghost cell: cell adjacent to a land boundary cell.              */
  /* Here cells are listed in order of wet, land boundary, OBC.      */
  cc2s = i_alloc_1d(ns2 + 1);
  cc2s[0] = 0;
  c = 1;
  vlm = i_alloc_2d(nz+1, sgrid->szcS);
  /* Wet cells                                                       */
  for (k = nz - 1; k >= 0; k--) {
    vlm[0][k] = 0;
    /* Wet cells surrounded by wet cells                             */
    for (cc = 1; cc <= ns2; cc++) {
      if (k < kbot[cc]) continue;
      npe = m->npe[cc];
      n = 0;
      for (j = 1; j <= npe; j++) {
	cn = neic[j][cc];
	jj = neij[j][cc];
	if (cn) n++;
      }
      if (n == npe && !(is_obc(cc, m->nobc, m->npts, m->loc))) {
	sgrid->c2cc[c] = cc;
	if (k == nz - 1) cc2s[cc] = c;
	vlm[cc2s[cc]][k] = c;
	c++;
      }
    }

    /* Land boundary cells                                           */
    for (cc = 1; cc <= ns2; cc++) {
      if (k < kbot[cc]) continue;
      npe = m->npe[cc];
      n = 0;
      for (j = 1; j <= npe; j++) {
	cn = neic[j][cc];
	jj = neij[j][cc];
	if (cn) n++;
      }
      if (n < npe && !(is_obc(cc, m->nobc, m->npts, m->loc))) {
	sgrid->c2cc[c] = cc;
	if (k == nz - 1) cc2s[cc] = c;
	vlm[cc2s[cc]][k] = c;
	c++;
      }
    }

    /* Open boundary cells                                           */
    for (cc = 1; cc <= ns2; cc++) {
      if (k < kbot[cc]) continue;
      npe = m->npe[cc];
      n = 0;
      for (j = 1; j <= npe; j++) {
	cn = neic[j][cc];
	jj = neij[j][cc];
	if (cn) n++;
      }
      if (n < npe && is_obc(cc, m->nobc, m->npts, m->loc)) {
	sgrid->c2cc[c] = cc;
	if (k == nz - 1) cc2s[cc] = c;
	vlm[cc2s[cc]][k] = c;
	c++;
      }
    }
    c = end_wet[k] + num_bdy[k] + num_gst[k] + 1;
  }
  sgrid->cc2s = cc2s;

  /*-----------------------------------------------------------------*/
  /* Get the spatial maps. This defines where cells lie in the       */
  /* sparse array relative to each other.                            */
  sgrid->c2c = i_alloc_2d(sgrid->szc, npem+1);
  sgrid->npe = i_alloc_1d(sgrid->szcS);
  memset(sgrid->npe, 0, sgrid->szcS * sizeof(int));
  sgrid->s2k = i_alloc_1d(sgrid->szc);
  memset(sgrid->s2k, 0, sgrid->szc * sizeof(int));
  sgrid->v2_t = sgrid->b2_t = sgrid->n2_t = 0;
  sgrid->v3_t = sgrid->b3_t = sgrid->n3_t = 0;

  /* Get the number of nodes for each cell                           */
  for (cc = 1; cc <= ns2; cc++) {
    for (j = 1; j <= m->npe[cc]; j++) {
      c = cc2s[cc];
      sgrid->npe[c] = m->npe[cc];
    }
  }

  /* Initialise to self-pointing                                     */
  for (c = 1; c <= sgrid->sgnum; c++) {
    for (j = 1; j <= npem; j++)
      sgrid->c2c[j][c] = c;
    sgrid->zp1[c] = c;
    sgrid->zm1[c] = c;
  }

  /* Get the cell to cell mapping.                                   */
  c1 = 1;	         /* Ghost vector counter                     */
  wgsts = i_alloc_1d(sgrid->sgsizS + 1); /* work array */
  memset(wgsts, 0, (sgrid->sgsizS + 1) * sizeof(int));
  for (k = nz - 1; k >= 0; k--) {
    gc = end_wet[k] + num_bdy[k] + 1;  /* 1st ghost cell in layer k  */
    c2 = end_wet[nz - 1] + num_bdy[nz - 1] + 1; /* Surface ghost     */
    for (cc = 1; cc <= ns2; cc++) {
      if (k < kbot[cc]) continue;
      c = vlm[cc2s[cc]][k];     /* Wet cell in layer k               */
      cs = vlm[cc2s[cc]][nz-1]; /* Surface layer wet cell            */
      npe = m->npe[cc];
      sgrid->s2k[c] = k;
      for (j = 1; j <= npe; j++) {

	cn = neic[j][cc];
	jj = neij[j][cc];
	cns = vlm[cc2s[cn]][nz-1]; /* Surface neighbour to c         */
	cn = vlm[cc2s[cn]][k];     /* Neighbour to c                 */
	if (!cn) {
	  /* Cell c is adjacent to a ghost cell. This includes the   */
	  /* first OBC cell (of the laus OBC ghost cells).           */
	  sgrid->c2c[j][c] = gc;
	  /* Note: npe and neic are not defined for ghost cells. We  */
	  /* are required to choose a direction for the ghost that   */
	  /* maps to the interior; the convention we choose is jo(). */
	  sgrid->c2c[jo(j, npe)][gc] = c;
	  if (cns)
	    sgrid->m2d[gc] = cns;
	  else {
	    if (wgsts[cs])
	      sgrid->m2d[gc] = wgsts[cs];
	    else
	      wgsts[cs] = sgrid->m2d[gc] = c2;
	  }
	  sgrid->m2d[c] = cs;
	  sgrid->wgst[gc] = c;
	  if (k == nz - 1) sgrid->npe[gc] = sgrid->npe[c];
	  /* Check the ghost mapping doesn't exceed 2D size          */
	  if(sgrid->m2d[gc]>sgrid->szcS)
	    hd_quit("pp_us: Ghost mapping > 2D size(%d): k=%d cc=%d c=%d cn=%d cns=%d gc=%d cs=%d\n",
		    sgrid->szcS,k,cc,c,cn,cns,gc,c2);
	  /*if (!(is_obc(cc, m->nobc, m->npts, m->loc))) {*/
	  sgrid->bpt[c1] = gc;
	  sgrid->bin[c1] = c;
	  /* The same convention used for c2c for ghost cells must   */
	  /* be used for dbin.                                       */
	  sgrid->dbin[c1] = jo(j, npe);
	  sgrid->dbpt[c1] = j;
	  c1++;
	  /*
	  }
	  */
	  gc++;
	  if (!cns) c2++;
	  /* Count the ghost cells (including the remaining laus-1   */
	  /* ghosts adjacent to OBCs).                               */
	  if (!cns && is_obc(cc, m->nobc, m->npts, m->loc)) {  /* OBC ghosts */
	    if (k == nz - 1) sgrid->n2_t += laus;
	    sgrid->n3_t += laus;
	    for (n = 1; n < laus; n++) {
	      cn = sgrid->c2c[j][c];
	      sgrid->c2c[j][cn] = gc;
	      /* Ghost direction to interior uses jo() convention    */
	      sgrid->c2c[jo(j, npe)][gc] = cn;
	      sgrid->m2d[gc] = (cns) ? cns : c2;
	      sgrid->wgst[gc] = c;
	      if (k == nz - 1) sgrid->npe[gc] = sgrid->npe[c];
	      gc++;
	      c2++;
	    }
	  } else {                             /* Land ghosts        */
	    if (k == nz - 1) sgrid->n2_t++;
	    sgrid->n3_t++;
	  }
	} else {
	  /* Cell c is a wet cell                                    */
	  sgrid->c2c[j][c] = cn;
	  sgrid->c2c[jj][cn] = c;
	  sgrid->wgst[c] = 0;
	  sgrid->m2d[c] = cs;
	}
      }
      /* Count the wet and OBC cells to process                      */
      if (!is_obc(cc, m->nobc, m->npts, m->loc)) {
	if (k == nz - 1) sgrid->v2_t++;
	sgrid->v3_t++;
      } else {
	if (k == nz - 1) sgrid->b2_t++;
	sgrid->b3_t++;
      }
    }
  }
  /* Free temp array */
  i_free_1d(wgsts);
  
  for (cc = 1; cc <= sgrid->nbpt; cc++) {
    c = sgrid->bpt[cc];
    c2 = sgrid->bin[cc];
    sgrid->s2k[c] = sgrid->s2k[c2];
  }

  if(sgrid->sgnum - num_sc != sgrid->v3_t + sgrid->b3_t + sgrid->n3_t)
    hd_quit("Inconsistent unstructured cell array and processing vector sizes : %d != %d\n",
	    sgrid->sgnum - num_sc, sgrid->v3_t + sgrid->b3_t + sgrid->n3_t);
  if (DEBUG("init_m"))
    dlog("init_m", "\nCell centre mappings created OK\n");

  /*-----------------------------------------------------------------*/
  /* Fill the cells to process vectors                               */
  sgrid->b2_t += sgrid->v2_t;
  sgrid->n2_t += sgrid->b2_t;
  sgrid->b3_t += sgrid->v3_t;
  sgrid->n3_t += sgrid->b3_t;

  sgrid->w3_t = i_alloc_1d(sgrid->n3_t + 1);
  sgrid->w2_t = i_alloc_1d(sgrid->n2_t + 1);
  maskb = i_alloc_1d(sgrid->szcS);

  c1 = 1;
  c2 = sgrid->v2_t + 1;
  c3 = sgrid->b2_t + 1;
  b1 = 1;
  b2 = sgrid->v3_t + 1;
  b3 = sgrid->b3_t + 1;
  memset(maskb, 0, sgrid->szcS * sizeof(int));
  for (k = nz - 1; k >= 0; k--) {
    gc = end_wet[k] + num_bdy[k] + 1;
    for (cc = 1; cc <= ns2; cc++) {
      if (k < kbot[cc]) continue;
      c = vlm[cc2s[cc]][k];     /* Wet cell in layer k               */
      cs = vlm[cc2s[cc]][nz-1]; /* Surface layer wet cell            */
      npe = m->npe[cc];
      for (j = 1; j <= npe; j++) {
	cn = neic[j][cc];
	jj = neij[j][cc];
	cns = vlm[cc2s[cn]][nz-1];
	cn = vlm[cc2s[cn]][k];
	if (!cn) {
	  /* Open boundary ghost cells */
	  if (!cns && is_obc(cc, m->nobc, m->npts, m->loc)) {
	    for (jj = 0; jj < laus; jj++) {
	      sgrid->w3_t[b3++] = gc;
	      if (k == nz - 1) {
		if (pc) printf("OBC ghost c=%d, laus=%d, cc=%d\n", gc, jj, c3);
		sgrid->w2_t[c3++] = gc;
	      }
	      gc++;
	    }
	  } else {
	    /* Land ghost cells */
	    sgrid->w3_t[b3++] = gc;
	    /*printf("Land ghost c=%d, cc=%d, c=%d, dir=%d %d %d %d %d %d\n", gc, b3, c, jj, cc, cs, k, m->iloc[cc], m->jloc[cc]);*/
	    if (k == nz - 1) {
	      if (pc) printf("Land ghost gc=%d, cc=%d, c=%d, dir=%d\n", gc, c3, c, jj);
	      sgrid->w2_t[c3++] = gc;
	    }
	    gc++;
	  }
	}
      }
      if (!is_obc(cc, m->nobc, m->npts, m->loc)) {
	sgrid->w3_t[b1++] = c;
	/*printf("centre wet k=%d c=%d, cs=%d cc=%d\n", k, c, cs, c1);*/
	if (k == nz - 1) {
	  if (pc) printf("Wet c=%d, cc=%d\n", c, c1);
	  sgrid->w2_t[c1++] = c;
	}
      } else {
	sgrid->w3_t[b2++] = c;
	if (k == nz - 1) {
	  if (pc) printf("OBC c=%d, cc=%d\n", c, c2);
	  sgrid->w2_t[c2++] = c;
	  maskb[c] = 1;
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the vertical maps                                           */
  gc = sgrid->sgnum - num_sc + 1;              /* Sediment cells     */
  for (k = nz - 1; k >= 0; k--) {
    for (cc = 1; cc <= ns2; cc++) {
      if (k < kbot[cc]) continue;
      c = vlm[cc2s[cc]][k];
      cs = vlm[cc2s[cc]][nz-1];
      npe = m->npe[cc];
      /* Wet cells                                                   */
      cp1 = (k == nz - 1) ? k : k + 1;
      cm1 = (k == 0) ? k : k - 1;
      cp1 = vlm[cc2s[cc]][cp1];
      cm1 = vlm[cc2s[cc]][cm1];
      sgrid->zp1[c] = cp1;
      sgrid->zm1[c] = cm1;
      /* Sediment cells                                              */
      if (k == kbot[cc]) {
	sgrid->zm1[c] = gc;
	sgrid->zp1[gc] = c;
	sgrid->wgst[gc] = c;
	sgrid->m2d[gc] = cs;
	gc++;
      }
      /* Ghost cells. Note: ghost sediments map up to the interior   */
      /* cell.                                                       */ 
      for (j = 1; j <= npe; j++) {
	cn = sgrid->c2c[j][c];
	if (sgrid->wgst[cn]) {
	  sgrid->zp1[cn] = sgrid->c2c[j][sgrid->zp1[c]];
	  sgrid->zm1[cn] = sgrid->c2c[j][sgrid->zm1[c]];
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the surface and bottom maps for centres                     */
  sgrid->bot_t = i_alloc_1d(sgrid->n2_t + 1);
  sgrid->sur_t = i_alloc_1d(sgrid->n2_t + 1);
  for (cc = 1; cc <= sgrid->n2_t; cc++) {
    c = cs = sgrid->w2_t[cc];
    sgrid->sur_t[cc] = c;
    while(c != sgrid->zm1[c]) {
      c =  sgrid->zm1[c];
    }
    sgrid->bot_t[cc] = sgrid->zp1[c];
  }
  if (DEBUG("init_m"))
    dlog("init_m", "\nCells to process vectors created OK\n");

  /*-----------------------------------------------------------------*/
  /* Find the number of cell edges                                   */
  locx = d_alloc_2d(npem+1, ns2+1);
  locy = d_alloc_2d(npem+1, ns2+1);
  sgrid->v3_e1 = sgrid->b3_e1 = sgrid->n3_e1 = sgrid->x3_e1 = 0;
  sgrid->v2_e1 = sgrid->b2_e1 = sgrid->n2_e1 = sgrid->x2_e1 = 0;
  mask = i_alloc_1d(sgrid->szc);
  for (k = 0; k <= nz; k++) {
    end_wet[k] = 0;
    num_gst[k] = 0;
  }
  /* Non OBC faces                                                   */
  memset(mask, 0, sgrid->szc * sizeof(int));
  mask[0] = 1;
  for (k = nz - 1; k >= 0; k--) {
    for (cc = 1; cc <= ns2; cc++) {
      npe = m->npe[cc];
      for (j = 1; j <= npe; j++) {
	locx[cc][j] = NOTVALID;
	locy[cc][j] = NOTVALID;
      }
    }
    for (cc = 1; cc <= ns2; cc++) {
      if (k < kbot[cc]) continue;
      c = vlm[cc2s[cc]][k];
      npe = m->npe[cc];
      for (j = 1; j <= npe; j++) {
	cn = cb = neic[j][cc];
	jj = neij[j][cc];
	cn = vlm[cc2s[cn]][k];
	if (cn) {
	  if (!mask[cn]) {
	    if (is_obc(cc, m->nobc, m->npts, m->loc) && is_obc(cb, m->nobc, m->npts, m->loc)) {
	      /* Tangential OBC edges                                */
	      if (k == nz - 1) sgrid->b2_e1++;
	      sgrid->b3_e1++;
	    } else {
	      /* Wet edges                                           */
	      if (k == nz - 1) sgrid->v2_e1++;
	      sgrid->v3_e1++;
	    }
	    end_wet[k]++;
	  }
	} else {
	  if (is_obce(cc, j, m) >= 0) {
	  /*if (is_obc(cc, m->nobc, m->npts, m->loc)) {*/
	    /* Normal OBC edges                                      */
	    if (k == nz - 1) sgrid->b2_e1++;
	    sgrid->b3_e1++;
	    end_wet[k]++;
	  } else {
	    double x = m->xloc[m->eloc[0][cc][j]];
	    double y = m->yloc[m->eloc[0][cc][j]];
	    int found = 0;
	    /* Land edges (ghost cells)                              */
	    if (k == nz - 1)
	      sgrid->n2_e1++;
	    sgrid->n3_e1++;
	    num_gst[k]++;
	  }
	}
      }
      mask[c] = 1;
    }
    /* Tangential ghost cells                                        */
    if (tage)
      num_gst[k] += find_edges_tan(sgrid, m, k, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
  }

  /* Find the end location of wet cells in each layer.               */
  for (k = nz - 2; k >= 0; k--) {
    c = end_wet[k + 1] + num_gst[k + 1];
    end_wet[k] += c;
  }

  /*-----------------------------------------------------------------*/
  /* Get the edge vertical layer map. In this process the edge cells */
  /* are re-ordered. Here we order so that wet cells precede ghost   */
  /* (land boundary) cells for every layer.                          */
  /* Wet cells                                                       */
  sgrid->nbe1 = sgrid->n3_e1 + 1;
  sgrid->nbe1S = sgrid->n2_e1 + 1;
  sgrid->nbpte1 = sgrid->n3_e1 + sgrid->x3_e1;
  sgrid->nbpte1S = sgrid->n2_e1 + sgrid->x2_e1;
  sgrid->b2_e1 += sgrid->v2_e1;
  sgrid->n2_e1 += sgrid->b2_e1;
  sgrid->x2_e1 += sgrid->n2_e1;
  num_sc = sgrid->x2_e1;
  sgrid->b3_e1 += sgrid->v3_e1;
  sgrid->n3_e1 += sgrid->b3_e1;
  sgrid->x3_e1 += sgrid->n3_e1;

  sgrid->sze = sgrid->x3_e1 + num_sc + 1;
  sgrid->szeS = sgrid->x2_e1 + 1;

  vle = i_alloc_2d(nz+1, sgrid->szeS);
  sgrid->e2c = i_alloc_2d(2, sgrid->sze);
  sgrid->e2e = i_alloc_2d(2, sgrid->sze);
  sgrid->c2e = i_alloc_2d(sgrid->szc, npem+1);
  maske = i_alloc_1d(sgrid->sze);
  emap = i_alloc_2d(npem+1, sgrid->szcS);
  tegl = malloc(sgrid->szeS * sizeof(point));

  c1 = 1;                   /* Wet cell counter                      */
  c2 = sgrid->v2_e1 + 1;    /* Open boundary counter                 */
  c3 = sgrid->b2_e1 + 1;    /* Closed boundary counter               */
  c4 = sgrid->n2_e1 + 1;    /* Tangential boundary counter           */
  b1 = 1;                   /* Wet cell counter                      */
  b2 = sgrid->v3_e1 + 1;    /* Open boundary counter                 */
  b3 = sgrid->b3_e1 + 1;    /* Closed boundary counter               */
  b4 = sgrid->n3_e1 + 1;    /* Tangential boundary counter           */
  sgrid->w2_e1 = i_alloc_1d(sgrid->szeS);
  sgrid->w3_e1 = i_alloc_1d(sgrid->sze);
  sgrid->zp1e = i_alloc_1d(sgrid->sze);
  sgrid->zm1e = i_alloc_1d(sgrid->sze);
  sgrid->bpte1 = i_alloc_1d(sgrid->nbpte1 + 1);  /* Edge boundary    */
  sgrid->bine1 = i_alloc_1d(sgrid->nbpte1 + 1);  /* Interior edge    */
  sgrid->bine2 = i_alloc_1d(sgrid->nbpte1 + 1);  /* Interior centre  */
  sgrid->bpte1S = i_alloc_1d(sgrid->nbpte1S + 1); /* Edge boundary   */
  sgrid->bine1S = i_alloc_1d(sgrid->nbpte1S + 1); /* Interior centre */
  sgrid->nbpte1 = sgrid->nbpte1S = 1;
  e = 1;
  memset(mask, 0, sgrid->szc * sizeof(int));
  sgrid->e2k = i_alloc_1d(sgrid->sze);

  /* Get the edge vectors                                            */
  for (k = nz - 1; k >= 0; k--) {
    ee = 1;
    for (cc = 1; cc <= sgrid->n2_t; cc++) {
      c = sgrid->w2_t[cc];
      for (j = 1; j <= npem; j++)
	emap[c][j] = 0;
    }
    /* Edges surrounded by wet cells (including OBC cells)           */
    for (cc = 1; cc <= ns2; cc++) {
      if (k < kbot[cc]) continue;
      c = vlm[cc2s[cc]][k];
      npe = m->npe[cc];
      for (j = 1; j <= npe; j++) {
	cn = cb = neic[j][cc];
	jj = neij[j][cc];
	cn = vlm[cc2s[cn]][k];
	if (cn) {
	  if (!mask[cn]) {
	    vle[ee][k] = e;
	    sgrid->e2c[e][0] = c;
	    sgrid->e2c[e][1] = cn;
	    sgrid->c2e[j][c] = e;
	    sgrid->c2e[jj][cn] = e;
	    sgrid->e2e[e][0] = j;
	    sgrid->e2e[e][1] = jj;
	    if (is_obc(cc, m->nobc, m->npts, m->loc) && is_obc(cb, m->nobc, m->npts, m->loc)) {
	      /* Edges between OBC cells (tangential OBC cells)    */
	      if (k == nz - 1) {
		sgrid->w2_e1[c2] = e;
		if (pe) printf("tan obc c1=%d c2=%d counter=%d e=%d\n",c,cn,c2,e);
		c2++;
	      }
	      maske[e] = max(obc_num(cc, m->nobc, m->npts, m->loc), obc_num(cb, m->nobc, m->npts, m->loc));
	      sgrid->w3_e1[b2] = e;
	      b2++;
	      sgrid->e2k[e] = k;
	    } else {

	      /* Edges between wet cells                           */
	      if (k == nz - 1) {
		sgrid->w2_e1[c1] = e;
		if (pe) printf("wet c1=%d c2=%d counter=%d edge=%d\n",c,cn,c1,e);
		c1++;
	      }
	      maske[e] = -1;
	      sgrid->w3_e1[b1] = e;
	      b1++;
	      sgrid->e2k[e] = k;
	    }
	    e++;
	    ee++;
	  }
	}
      }
      mask[c] = 1;
    }

    /* Edges adjacent to a land boundary (including OBC cells)       */
    for (cc = 1; cc <= ns2; cc++) {
      if (k < kbot[cc]) continue;
      c = vlm[cc2s[cc]][k];
      npe = m->npe[cc];
      for (j = 1; j <= npe; j++) {
	cn = neic[j][cc];
	jj = neij[j][cc];
	cn = vlm[cc2s[cn]][k];
	if (!cn) {
	  gc = sgrid->c2c[j][c];
	  vle[ee][k] = e;
	  sgrid->e2c[e][0] = c;
	  sgrid->e2c[e][1] = gc;
	  sgrid->c2e[j][c] = e;
	  sgrid->e2e[e][0] = j;
	  /* Use the jo() convention as was used to define the       */
	  /* direction of the interior cell for c2c.                 */
	  sgrid->c2e[jo(j, npe)][gc] = e;
	  sgrid->e2e[e][1] = jo(j, npe);
	  /* Edges between OBC cells and land (normal OBC cells)     */
	  if ((ei = is_obce(cc, j, m)) >= 0) {
	    /*if (is_obc(cc, m->nobc, m->npts, m->loc)) {*/
	    if (k == nz - 1) {
	      sgrid->w2_e1[c2] = e;
	      if (pe) printf("nor obc c=%d gc=%d counter=%d edge=%d j=%d\n",c,gc,c2,e,j);
	      c2++;
	    }
	    maske[e] = ei;
	    sgrid->w3_e1[b2] = e;
	    b2++;
	    sgrid->e2k[e] = k;
	  } else {
	    /* Normal ghost cells                                    */
	    if (k == nz - 1) {
	      sgrid->bpte1S[sgrid->nbpte1S] = e;
	      sgrid->bine1S[sgrid->nbpte1S++] = c;
	      sgrid->w2_e1[c3] = e;
	      if (pe) printf("ghost c=%d gc=%d counter=%d ie=%d\n",c,gc,c3,e);
	      c3++;
	    }
	    maske[e] = -2;
	    sgrid->bpte1[sgrid->nbpte1] = e;
	    sgrid->bine1[sgrid->nbpte1++] = c;
	    sgrid->w3_e1[b3] = e;
	    b3++;
	    sgrid->e2k[e] = k;
	  }
	  e++;
	  ee++;
	}
      }
    }

    /* Tangential ghost cells.                                       */
    if (tage)
      find_edges_tan(sgrid, m, k, &e, &ee, &c4, &b4, maske, vle, tegl);
    e = end_wet[k] + num_gst[k] + 1;
  }
  d_free_2d(locx);
  d_free_2d(locy);
  i_free_2d(emap);
  reorder_edges(sgrid);
  sgrid->nbpte1S--;
  sgrid->nbpte1--;
  sgrid->nbe1S--;
  sgrid->nbe1--;

  /*-----------------------------------------------------------------*/
  /* Get the surface and bottom maps for edges                       */
  /* At this stage the bottom edge is set to be the edge adjacent to */
  /* the deepest of the cell centres sharing the edge. This is done  */
  /* so that ghost vertices are included in the vertex vector. The   */
  /* bottom edge is reset to the edge adjacent to shallowest cell    */
  /* after the vertex vector is established.                         */
  c2cc = i_alloc_1d(sgrid->szc); 
  for (cc = 1; cc <= sgrid->n2_t; cc++) {
    c = sgrid->w2_t[cc];
    c2cc[c] = cc;
  }
  sgrid->bot_e1 = i_alloc_1d(sgrid->szeS);
  sgrid->sur_e1 = i_alloc_1d(sgrid->szeS);
  for (ee = 1; ee <= sgrid->n2_e1; ee++) {
    int cp;
    e = sgrid->w2_e1[ee];
    sgrid->sur_e1[ee] = e;
    c1 = sgrid->e2c[e][0];
    j = sgrid->e2e[e][0];
    c2 = sgrid->e2c[e][1];
    /* Set any ghost cells to wet cells */
    if (sgrid->wgst[c1]) {
      c1 = c2;
      j = sgrid->e2e[e][1];
    }
    if (sgrid->wgst[c2]) c2 = c1;
    /* Find the deepest cell */
    if (kbot[sgrid->c2cc[c2]] < kbot[sgrid->c2cc[c1]]) {
      c1 = c2;
      j = sgrid->e2e[e][1];
    }
    cp = c1;
    /*while (!sgrid->wgst[sgrid->zm1[c1]] && !sgrid->wgst[sgrid->zm1[c2]]) {*/
    while (c1 != sgrid->zm1[c1]) {
      cp = c1;
      c1 = sgrid->zm1[c1];
      c2 = sgrid->zm1[c2];
      if(c1==0||c2==0)hd_quit("Error in surface and bottom edge maps");
    }
    c1 = sgrid->zp1[c1];
    sgrid->bot_e1[ee] = sgrid->c2e[j][c1];
  }

  /* Check the edge to centre maps                                   */
  for (ee = 1; ee <= sgrid->n3_e1; ee++) {
    int found = 0;
    e = sgrid->w3_e1[ee];
    c = sgrid->e2c[e][0];
    if (sgrid->wgst[c]) c = sgrid->e2c[e][1];
    cc = sgrid->c2cc[c];
    npe = m->npe[cc];
    for (n = 1; n <= npe; n++) {
      e1 = sgrid->c2e[n][c];
      if (e == e1) {
	found = 1;
	break;
      }
    }
    if (!found) hd_quit("Edge mapping error: e = %d, cs = %d, cc = %d\n",e, sgrid->m2d[c], cc);
  }

  /*-----------------------------------------------------------------*/
  /* Get the maps in the e2e[1] (ep) and e2e[0] (em) directions.     */
  sgrid->ep = i_alloc_1d(sgrid->sze);
  sgrid->em = i_alloc_1d(sgrid->sze);
  for (ee = 1; ee <= sgrid->n3_e1; ee++) {
    e = sgrid->w3_e1[ee];
    j = sgrid->e2e[e][1];
    sgrid->ep[e] = e2e(sgrid, e, j);
    j = sgrid->e2e[e][0];
    sgrid->em[e] = e2e(sgrid, e, j);
  }
  /* Get the intereior edges to ghosts from the cell centres         */
  for (ee = 1; ee <= sgrid->nbpte1; ee++) {
    e = sgrid->bpte1[ee];
    c = sgrid->bine1[ee];
    cs = sgrid->m2d[c];
    sgrid->bine2[ee] = sgrid->bine1[ee];
    for (j = 1; j <= sgrid->npe[cs]; j++) {
      if (sgrid->c2e[j][c] == e) {
	/*sgrid->bine1[ee] = e2e(sgrid, e, jo(j, sgrid->npe[cs]));*/
	if (usejo)
	  sgrid->bine1[ee] = sgrid->c2e[jo(j, sgrid->npe[cs])][c];
	else
	  sgrid->bine1[ee] = sgrid->c2e[jocw(sgrid, c, j)][c];
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Vertical maps                                                   */
  gc = sgrid->x3_e1 + 1;
  /* Initialize                                                      */
  for (ee = 1; ee <= sgrid->n3_e1; ee++) {
    e = sgrid->w3_e1[ee];
    sgrid->zp1e[e] =  e;
    sgrid->zm1e[e] =  e;
  }
  /* Wet maps                                                        */
  for (ee = 1; ee <= sgrid->n3_e1; ee++) {
    e = sgrid->w3_e1[ee];
    sgrid->zm1e[e] =  eb = zm1e(sgrid, e);
    /*sgrid->zp1e[e] = zp1e(sgrid, e);*/
    if (e != eb) sgrid->zp1e[eb] = e;
  }
  /* Sediment maps                                                   */
  for (ee = 1; ee <= sgrid->n2_e1; ee++) {
    e = sgrid->bot_e1[ee];
    sgrid->zm1e[e] = gc;
    sgrid->zm1e[gc] = gc;
    sgrid->zp1e[gc] = e;
    gc++;
  }
  sgrid->m2de = i_alloc_1d(sgrid->sze);
  e2ee = i_alloc_1d(sgrid->sze);
  for (ee = 1; ee <= sgrid->n2_e1; ee++) {
    int ed = 0;
    e = es = sgrid->w2_e1[ee];
    e2ee[e] = ee;
    while (e != sgrid->zm1e[e]) {
      sgrid->m2de[e] = es;
      e = sgrid->zm1e[e];
    }
    sgrid->m2de[e] = es; /* Sediment                                 */
  }

  /* Check edges map to the bottom coordinate                        */
  for (ee = 1; ee <= sgrid->n2_e1; ee++) {
    e = es = sgrid->w2_e1[ee];
    eb = sgrid->bot_e1[ee];
    n = 1;
    while (e != eb) {
      e = sgrid->zm1e[e];
      if (n > sgrid->nz) {
	c = sgrid->e2c[e][0];
	cs = sgrid->m2d[c];
	cc = sgrid->c2cc[cs];
	hd_quit("Edge mapping error: es=%d eb=%d e=%d zm1e[e]=%d cc=%d\n",
				 es, eb, e, sgrid->zm1e[e], cc);
      }
      n++;
    }
  }
  /* Check the surface mappings                                        */
  for (ee = 1; ee <= sgrid->n3_e1; ee++) {
    e = sgrid->w3_e1[ee];
    es = sgrid->m2de[e];
    c = sgrid->e2c[e][0];
    cs = sgrid->m2d[c];
    c1 = sgrid->c2cc[cs];
    c = sgrid->e2c[e][1];
    cs = sgrid->m2d[c];
    c2 = sgrid->c2cc[cs];    
    if (es == 0) {
      hd_warn("Can't map edge between cells %d[%f %f] bathy = %f and\n", 
	      c1, m->xloc[m->eloc[0][c1][0]], m->yloc[m->eloc[0][c1][0]], bathy[c1]);
      hd_warn("                             %d[%f %f] bathy = %f.\n", 
	      c2, m->xloc[m->eloc[0][c2][0]], m->yloc[m->eloc[0][c2][0]], bathy[c2]);
      hd_warn("Check that bathy or try different BATHY_INTERP_RULE!\n");
      hd_quit("Can't map edge %d to the surface. Check the runlog.\n", e);
    }
  }

  /* Get the edge sign vector                                        */
  sgrid->eSc = i_alloc_2d(sgrid->szcS, npem + 1);
  for (cc = 1; cc <= sgrid->b2_t; cc++) {
    c = sgrid->w2_t[cc];
    npe = sgrid->npe[c];
    for (n = 1; n <= npe; n++) {
      e = sgrid->c2e[n][c];
      if (c == sgrid->e2c[e][1])
	sgrid->eSc[n][c] = 1;    /* e points out of cell c */
      else
	sgrid->eSc[n][c] = -1;   /* e points into cell c   */
    }
  }

  /* The edge vectors do not contain the sediment edges; these cells */
  /* are accessed using zme1[bot_e1]. Truncate the edge vector size  */
  /* so that sediment edges are omitted.                             */
  /*sgrid->n3_e1 -= num_sc;*/
  if (DEBUG("init_m"))
    dlog("init_m", "\nEdge mappings created OK\n");

  /*-----------------------------------------------------------------*/
  /* Get the centre and edge locations                               */
  sgrid->cellx = d_alloc_1d(sgrid->szcS);
  sgrid->celly = d_alloc_1d(sgrid->szcS);

  for (cc =sgrid->b2_t+1; cc <= sgrid->n2_t; cc++) {
    c = sgrid->w2_t[cc];
    sgrid->cellx[c] = NaN;
    sgrid->celly[c] = NaN;
  }
  locx = d_alloc_2d(2, sgrid->szeS);
  locy = d_alloc_2d(2, sgrid->szeS);
  for (cc = 1; cc <= sgrid->b2_t; cc++) {
    c = sgrid->w2_t[cc];
    c1 = sgrid->c2cc[c];
    npe = sgrid->npe[c];
    /*
    sgrid->cellx[c] = params->x[c1][0];
    sgrid->celly[c] = params->y[c1][0];
    */
    sgrid->cellx[c] = m->xloc[m->eloc[0][c1][0]];
    sgrid->celly[c] = m->yloc[m->eloc[0][c1][0]];

    for (j = 1; j <= npe; j++) {
      e = sgrid->c2e[j][c];
      locx[e][0] = m->xloc[m->eloc[0][c1][j]];
      locy[e][0] = m->yloc[m->eloc[0][c1][j]];
      locx[e][1] = m->xloc[m->eloc[1][c1][j]];
      locy[e][1] = m->yloc[m->eloc[1][c1][j]];
      /*
      locx[e][0] = params->x[c1][j];
      locy[e][0] = params->y[c1][j];
      locx[e][1] = params->x[c1][jp(j, npe)];
      locy[e][1] = params->y[c1][jp(j, npe)];
      */
    }
  }

  /* Reset the perimeter cells to be the centre of mass of the       */
  /* polygon.                                                        */
  /*
  for (cc = 1; cc <= sgrid->b2_t; cc++) {
    c = sgrid->w2_t[cc];
    if (maskb[c]) {
      double x = 0.0, y = 0.0;
      c1 = sgrid->c2cc[c];
      npe = sgrid->npe[c];
      for (j = 1; j <= npe; j++) {
	x += m->xloc[m->eloc[0][c1][j]];
	y += m->yloc[m->eloc[0][c1][j]];
      }
      sgrid->cellx[c] = x / (double)npe;
      sgrid->celly[c] = y / (double)npe;
    }
  }
  */

  /* Get the lengths of the edges. params->h1 stores the length of   */
  /* the edges as inferred from the vertex locations in the input    */
  /* file. params->h2 stores the distance from the cell centre to    */
  /* the middle of each edge. 
  /* h1au1 : edge length (length between vertices)                   */
  /* h2au1 : length between centres for each edge.                   */
  /* h1acell : length between edges, in the upstream (index 0) dir.  */
  /*           i.e. distance across a cell                           */
  /* Length between edges: cell centred                              */
  /* thetau1 : angle between u1 vector and x axis                    */
  /* thetau2 : angle between u1 vector and y axis                    */
  sgrid->h2au1 = d_alloc_1d(sgrid->szeS);
  sgrid->h1au1 = d_alloc_1d(sgrid->szeS);
  sgrid->cellarea = d_alloc_1d(sgrid->n2_t + 1);
  sgrid->edgearea = d_alloc_1d(sgrid->szeS);
  sgrid->h1acell = d_alloc_1d(sgrid->szeS);
  sgrid->hacell = d_alloc_2d(sgrid->n2_t + 1, sgrid->npem + 1);
  sgrid->thetau1 = d_alloc_1d(sgrid->szeS);
  sgrid->thetau2 = d_alloc_1d(sgrid->szeS);
  /* Edge lengths                                                    */
  for (ee = 1; ee <= sgrid->n2_e1; ee++) {
    e = sgrid->w2_e1[ee];
    c = sgrid->e2c[e][0];
    j = sgrid->e2e[e][0];
    sgrid->h2au1[e] = 0.0;
    if (sgrid->wgst[c]) {
      c = sgrid->e2c[e][1];
      j = sgrid->e2e[e][1];
    }
    cc = sgrid->c2cc[c];
    /* Note: h1[][0:npe-1] and h2[][0:npe-1]                         */
    sgrid->h1au1[e] = params->h1[cc][j-1];
  }

  /* Lengths between centres                                         */
  for (cc = 1; cc <= sgrid->b2_t; cc++) {
    c = sgrid->w2_t[cc];
    c1 = sgrid->c2cc[c];
    npe = sgrid->npe[c];
    for (j = 1; j <= npe; j++) {
      e = sgrid->c2e[j][c];
      if (e) sgrid->h2au1[e] += params->h2[c1][j-1];
    }
  }

  /* Lengths between edges in the e2c[0] edge direction. (only valid */
  /* for symmetrical polygons, e.g. square, hex, oct).               */
  for (ee = 1; ee <= sgrid->n2_e1; ee++) {
    e = sgrid->w2_e1[ee];
    c = sgrid->e2c[e][0];
    npe = sgrid->npe[c];
    if (!sgrid->wgst[c]) {
      j = sgrid->e2e[e][0];
      cc = sgrid->c2cc[c];
      n = j + npe / 2;
      if (n > npe) n -= npe;
      sgrid->h1acell[e] = params->h2[cc][j-1] + params->h2[cc][n-1];
    }
  }
  /* Lengths between edges and the centre, cell centred. The length  */
  /* between two edges is twice this value, but is strictly only     */
  /* valid for symmetrical polygons.                                 */
  for (cc = 1; cc <= sgrid->b2_t; cc++) {
    c = sgrid->w2_t[cc];
    npe = sgrid->npe[c];
    for (n = 1; n <= npe; n++) {
      if (usejo)
	j = jo(n, npe);
      else
	j = jocw(sgrid, c, n);
      e = sgrid->c2e[n][c];
      e1 = sgrid->c2e[j][c];
      sgrid->hacell[n][c] = 0.5 * sgrid->h2au1[e];
      /*sgrid->hacell[n][c] = 0.5 * (sgrid->h2au1[e] + sgrid->h2au1[e1]);*/
      /*
      if (c == sgrid->e2c[e][0])
	sgrid->hacell[n][c] = sgrid->hacell[j][c] = sgrid->h1acell[e];
      */
    }
  }

  /* Cell area                                                       */
  memset(sgrid->edgearea, 0, sgrid->szeS * sizeof(double));
  for (cc = 1; cc <= sgrid->b2_t; cc++) {
    c = sgrid->w2_t[cc];
    c1 = sgrid->c2cc[c];
    npe = sgrid->npe[c];
    sgrid->cellarea[c] = 0.0;
    for (j = 1; j <= npe; j++) {
      e = sgrid->c2e[j][c];
      /* Edge area : Ringler et al, (2010) Eq. 48                    */
      sgrid->edgearea[e] = 0.5 * (sgrid->h1au1[e] * sgrid->h2au1[e]);
      sgrid->cellarea[c] += 0.5 * (sgrid->h1au1[e] * params->h2[c1][j-1]);
    }
  }

  /* Angle to the horizontal. Here the direction of the normal       */
  /* vector is defined; this must be consistent with the conventions */
  /* used in reoder_edges().                                         */
  /*
  for (ee = 1; ee <= sgrid->n2_e1; ee++) {
    double sgn, eps = 1e-12;
    double pid = acos(-1.0);
    double pi2 = pid / 2.0;
    double pi3 = 3.0 * pid / 2.0;
    e = sgrid->w2_e1[ee];
    j = sgrid->e2e[e][0];
    c = sgrid->e2c[e][0];
    c1 = sgrid->c2cc[c];
    sgn = 0.0;
    if(sgrid->wgst[c]) {
      c = sgrid->e2c[e][1];
      j = sgrid->e2e[e][1];
      c1 = sgrid->c2cc[c];
      sgn = 1.0;
    }
    sgrid->thetau1[e] = params->a2[c1][j-1];

    if (sgn) sgrid->thetau1[e] += pid;
    if (sgrid->thetau1[e] >= 2.0 * pid) sgrid->thetau1[e] -= 2.0 * pid;
    sgn = (sgrid->thetau1[e] >= 0.0) ? 1.0 : -1.0;
    if (fabs(sgrid->thetau1[e] - pi2) < eps)
      sgrid->thetau1[e] = sgn * pi2;
    if (fabs(sgrid->thetau1[e] - pid) < eps)
      sgrid->thetau1[e] = sgn * pid;
    if (fabs(sgrid->thetau1[e] - pi3) < eps)
      sgrid->thetau1[e] = sgn * pi3;
  }
  */

  for (cc = 1; cc <= sgrid->b2_t; cc++) {
    double sgn, eps = 1e-12;
    double pid = acos(-1.0);
    double pi2 = pid / 2.0;
    double pi3 = 3.0 * pid / 2.0;
    c = sgrid->w2_t[cc];
    c1 = sgrid->c2cc[c];
    npe = sgrid->npe[c];
    for (j = 1; j <= npe; j++) {
      e = sgrid->c2e[j][c];
      if (e) {
	sgrid->thetau1[e] = params->a2[c1][j-1];
	if(sgrid->eSc[j][c] == 1) sgrid->thetau1[e] += pid;

	if (sgrid->thetau1[e] >= 2.0 * pid) sgrid->thetau1[e] -= 2.0 * pid;
	sgn = (sgrid->thetau1[e] >= 0.0) ? 1.0 : -1.0;
	if (fabs(sgrid->thetau1[e] - pi2) < eps)
	  sgrid->thetau1[e] = sgn * pi2;
	if (fabs(sgrid->thetau1[e] - pid) < eps)
	  sgrid->thetau1[e] = sgn * pid;
	if (fabs(sgrid->thetau1[e] - pi3) < eps)
	  sgrid->thetau1[e] = sgn * pi3;

	sgrid->thetau2[e] = sgrid->thetau1[e] + PI / 2.0;
	if (sgrid->thetau2[e] >= 2.0 * pid) sgrid->thetau2[e] -= 2.0 * pid;
	sgn = (sgrid->thetau2[e] >= 0.0) ? 1.0 : -1.0;
	if (fabs(sgrid->thetau2[e] - pi2) < eps)
	  sgrid->thetau2[e] = sgn * pi2;
	if (fabs(sgrid->thetau2[e] - pid) < eps)
	  sgrid->thetau2[e] = sgn * pid;
	if (fabs(sgrid->thetau2[e] - pi3) < eps)
	  sgrid->thetau2[e] = sgn * pi3;
      }
    }
  }
  if (DEBUG("init_m"))
    dlog("init_m", "\nMesh metrics created OK\n");

  /*-----------------------------------------------------------------*/
  /* Find the number of cell vertices                                */
  sgrid->v3_e2 = sgrid->b3_e2 = sgrid->n3_e2 = 0;
  sgrid->v2_e2 = sgrid->b2_e2 = sgrid->n2_e2 = 0;
  emap = i_alloc_2d(2, sgrid->sze);
  for (ee = 1; ee <= sgrid->n3_e1; ee++) {
    e = sgrid->w3_e1[ee];
    emap[e][0] = 0;
    emap[e][1] = 0;
  }
  sgrid->szvS = 0;

  if (newvert) {
    find_vertices(sgrid, locx, locy, maske, emap);
  } else {
  /* Count the number of 2D vertices                                */
  for (cc = 1; cc <= sgrid->b2_t; cc++) {
    c = sgrid->w2_t[cc];
    npe = sgrid->npe[c];
    for (j = 1; j <= npe; j++) {
      e = sgrid->c2e[j][c];
      if (emap[e][0] && emap[e][1]) continue;
      for (ee = 1; ee <= sgrid->n2_e1; ee++) {
	e1 = sgrid->w2_e1[ee];
	if (e1 == e) continue;
	check_vert2d(sgrid, e, e1, 0, 0, locx, locy, emap, maske, NULL, NULL, NULL);
	check_vert2d(sgrid, e, e1, 0, 1, locx, locy, emap, maske, NULL, NULL, NULL);
	check_vert2d(sgrid, e, e1, 1, 0, locx, locy, emap, maske, NULL, NULL, NULL);
	check_vert2d(sgrid, e, e1, 1, 1, locx, locy, emap, maske, NULL, NULL, NULL);
      }
    }
  }

  /* Populate the 2D vertices. These are ordered as:                 */
  /* 1 to v2_e2 : vertices with both edges that are wet              */
  /* v2_e2 to b2_e2 : vertices with both edges that are OBCs         */
  /* b2_e2 to n2_e2 : vertices with one edge that is a ghost         */
  sgrid->b2_e2 += sgrid->v2_e2;
  sgrid->n2_e2 += sgrid->b2_e2;
  sgrid->szvS = sgrid->n2_e2 + 1;
  sgrid->w2_e2 = i_alloc_1d(sgrid->szvS);
  sgrid->nve = i_alloc_1d(sgrid->szvS);
  memset(sgrid->nve, 0, sgrid->szvS * sizeof(int));
  for (ee = 1; ee <= sgrid->n3_e1; ee++) {
    e = sgrid->w3_e1[ee];
    emap[e][0] = 0;
    emap[e][1] = 0;
  }
  b1 = 1;
  b2 = sgrid->v2_e2 + 1;
  b3 = sgrid->b2_e2 + 1;

  sgrid->n2_e2 = sgrid->b2_e2;
  sgrid->b2_e2 = sgrid->v2_e2;
  sgrid->v2_e2 = 0;

  /* Set the vertices to process                                     */
  for (cc = 1; cc <= sgrid->b2_t; cc++) {
    c = sgrid->w2_t[cc];
    npe = sgrid->npe[c];
    for (j = 1; j <= npe; j++) {
      e = sgrid->c2e[j][c];
      if (emap[e][0] && emap[e][1]) continue;
      for (ee = 1; ee <= sgrid->n2_e1; ee++) {
	e1 = sgrid->w2_e1[ee];
	if (e1 == e) continue;
	check_vert2d(sgrid, e, e1, 0, 0, locx, locy, emap, maske, &b1, &b2, &b3);
	check_vert2d(sgrid, e, e1, 0, 1, locx, locy, emap, maske, &b1, &b2, &b3);
	check_vert2d(sgrid, e, e1, 1, 0, locx, locy, emap, maske, &b1, &b2, &b3);
	check_vert2d(sgrid, e, e1, 1, 1, locx, locy, emap, maske, &b1, &b2, &b3);
      }
    }
  }
  }

  /* Get the maximum number of edges from vertices                   */
  sgrid->nvem = 0;
  for (vv = 1; vv <= sgrid->n2_e2; vv++) {
    v = sgrid->w2_e2[vv];
    if (sgrid->nve[v] > sgrid->nvem) sgrid->nvem = sgrid->nve[v];
  }

  /* Get the 2D map from vertices to edges                           */
  v2e = i_alloc_2d(sgrid->nvem+1, sgrid->szvS);
  maskv = i_alloc_1d(sgrid->szvS);
  for (vv = 1; vv <= sgrid->n2_e2; vv++) {
    v = sgrid->w2_e2[vv];
    maskv[v] = 1;
  }
  for (ee = 1; ee <= sgrid->n2_e1; ee++) {
    e = sgrid->w2_e1[ee];
    for (j = 0; j <= 1; j++) {
      v = emap[e][j];
      v2e[v][maskv[v]] = e;
      maskv[v]++;
    }
  }

  /* Count the number of 3D vertices                                 */
  for (vv = 1; vv <= sgrid->n2_e2; vv++) {
    int found;
    v = sgrid->w2_e2[vv];
    /* Get the deepest edge surrounding v                            */
    found = 0;
    for (n = 1; n <= sgrid->nve[v]; n++) {
      e = v2e[v][n];
      if (e) {
	e2 = sgrid->bot_e1[e2ee[sgrid->m2de[e]]];
	if (!found) {
	  ei = e;
	  eb = e2;
	  found = 1;;
	} else {
	  if (sgrid->e2k[e2] < sgrid->e2k[eb]) {
	    ei = e;
	    eb = e2;
	  }
	}
      }
    }
    check_vert3d(sgrid, v, vv, ei, eb, NULL, NULL, NULL);
  }

  sgrid->b3_e2 += sgrid->v3_e2;
  sgrid->n3_e2 += sgrid->b3_e2;
  sgrid->szv = sgrid->n3_e2 + 1;

  sgrid->w3_e2 = i_alloc_1d(sgrid->szv);
  sgrid->v2e = i_alloc_2d(sgrid->nvem+1, sgrid->szv);
  sgrid->e2v = i_alloc_2d(2, sgrid->sze);
  sgrid->m2dv = i_alloc_1d(sgrid->szv);
  sgrid->zp1v = i_alloc_1d(sgrid->szv);
  sgrid->zm1v = i_alloc_1d(sgrid->szv);
  i_free_1d(maskv);
  maskv = i_alloc_1d(sgrid->szv);

  for (v = 1; v <= sgrid->n3_e2; v++) {
    sgrid->zp1v[v] = v;
    sgrid->zm1v[v] = v;
    sgrid->m2dv[v] = v;
    maskv[v] = 1;
  }
  for (v = 1; v < sgrid->szv; v++) maskv[v] = 1;

  c1 = sgrid->n2_e2 + 1;
  c2 = c1 + sgrid->v3_e2 - sgrid->v2_e2;
  c3 = c1 + sgrid->b3_e2 - sgrid->b2_e2;

  sgrid->n3_e2 = sgrid->b3_e2;
  sgrid->b3_e2 = sgrid->v3_e2;
  sgrid->v3_e2 = 0;

  /* Set the 3D vertices                                             */
  for (vv = 1; vv <= sgrid->n2_e2; vv++) {
    int found;
    v = sgrid->w2_e2[vv];
    /* Get the deepest edge surrounding v                            */
    found = 0;
    for (n = 1; n <= sgrid->nve[v]; n++) {
      e = v2e[v][n];
      if (e) {
	e2 = sgrid->bot_e1[e2ee[sgrid->m2de[e]]];
	if (!found) {
	  ei = e;
	  eb = e2;
	  found = 1;;
	} else {
	  if (sgrid->e2k[e2] < sgrid->e2k[eb]) {
	    ei = e;
	    eb = e2;
	  }
	}
      }
    }
    check_vert3d(sgrid, v, vv, ei, eb, &c1, &c2, &c3);
  }

  /* Get the maps from edges to vertices and vertices to edges.      */
  /* Note: vertices below the bottom edge (sediment vertices) are    */
  /* not mapped (i.e. they map to 0), hence if a loop extends over   */
  /* n3_e2 and v2e is used then the edge should be checked using     */
  /* if(e) before it is used.                                        */
  for (ee = 1; ee <= sgrid->n2_e1; ee++) {
    es = sgrid->w2_e1[ee];
    eb = sgrid->bot_e1[ee];
    for (j = 0; j <= 1; j++) {
      e = es;
      c = sgrid->e2c[es][0];
      jj = sgrid->e2e[e][0];
      v = vs = emap[e][j];
      sgrid->e2v[e][j] = v;
      sgrid->v2e[v][maskv[v]] = e;
      maskv[v]++;
      while (e < eb) {
	e = sgrid->zm1e[e];
	v = sgrid->zm1v[v];
	sgrid->e2v[e][j] = v;
	sgrid->v2e[v][maskv[v]] = e;
	maskv[v]++;
      }
    }
  }
  for (vv = 1; vv <= sgrid->n2_e2; vv++) {
    v = sgrid->w2_e2[vv];
    if(maskv[v]-1 > sgrid->nve[v])
      hd_quit("Vertex mapping error: v = %d(%d edges); v2e has %d maps\n",v, sgrid->nve[v], maskv[v]-1);
  }
  i_free_1d(maskv);
  i_free_2d(v2e);

  /* Re-order the e2v maps so that index 0 is upstream. The upstream */
  /* vertex is defined as that which to the left of the normal       */
  /* velocity vector, i.e. in the direction of k x u1, (i.e. the     */
  /* normal vector rotated 90 degrees anti-clockwise).               */
  /* See Thuburn et al, (2009), p 8324.                              */
  for (ee = 1; ee <= sgrid->n3_e1; ee++) {
    int nn, n1 , n2;
    e = sgrid->w3_e1[ee];
    v1 = sgrid->e2v[e][0];
    v2 = sgrid->e2v[e][1];
    c = sgrid->e2c[e][0];
    n = sgrid->e2e[e][0];
    if (sgrid->wgst[c]) {
      c = sgrid->e2c[e][1];
      n = sgrid->e2e[e][1];
    }
    cs = sgrid->m2d[c];
    npe = sgrid->npe[cs];
    /* Get the index n that points to the edge e from the centre c   */
    for (nn = 1; nn <= npe; nn++) {
      e1 = sgrid->c2e[nn][c];
      if (e1 == e) {
	n = nn;
	break;
      }
    }
    /* Get the index n1 for the edge that has v1 as a vertex         */
    for (nn = 1; nn <= npe; nn++) {
      e1 = sgrid->c2e[nn][c];
      if (e == e1) continue;
      if (sgrid->e2v[e1][0] == v1 || sgrid->e2v[e1][1] == v1) {
	n1 = nn;
	break;
      }
    }
    /* Get the index n2 for the edge that has v2 as a vertex         */
    for (nn = 1; nn <= npe; nn++) {
      e2 = sgrid->c2e[nn][c];
      if (e == e2) continue;
      if (sgrid->e2v[e2][0] == v2 || sgrid->e2v[e2][1] == v2) {
	n2 = nn;
	break;
      }
    }

    /* Default: set vertices using inward point vectors              */
    if (sgrid->eSc[n][cs] == -1) {
      /* Special case for a vertex on the indexing transition        */
      if (n == 1 || n == npe) {
	if (n1 > n2) {
	  sgrid->e2v[e][0] = v2;
	  sgrid->e2v[e][1] = v1;
	} else {
	  sgrid->e2v[e][0] = v1;
	  sgrid->e2v[e][1] = v2;
	}
      } else {
	/* All other vertex indices                                  */
	if (n1 < n2) {
	  sgrid->e2v[e][0] = v2;
	  sgrid->e2v[e][1] = v1;
	} else {
	  sgrid->e2v[e][0] = v1;
	  sgrid->e2v[e][1] = v2;
	}
      }
    } else {  
      /* Outward pointing vectors if inward vector points to         */
      /* a ghost cell.                                               */
      if (n == 1 || n == npe) {
	if (n1 < n2) {
	  sgrid->e2v[e][0] = v2;
	  sgrid->e2v[e][1] = v1;
	} else {
	  sgrid->e2v[e][0] = v1;
	  sgrid->e2v[e][1] = v2;
	}
      } else {
	/* All other vertex indices                                  */
	if (n1 > n2) {
	  sgrid->e2v[e][0] = v2;
	  sgrid->e2v[e][1] = v1;
	} else {
	  sgrid->e2v[e][0] = v1;
	  sgrid->e2v[e][1] = v2;
	}
      }
    }
  }

  /* Check the vertex to edge mappings                               */
  for (ee = 1; ee <= sgrid->n3_e1; ee++) {
    int verbose = 0;
    e = sgrid->w3_e1[ee];
    es = sgrid->m2de[e];
    if (verbose) {
      if (e == verbose) {
	c = sgrid->e2c[es][0];
	n = sgrid->e2e[es][0];
	if (sgrid->wgst[c]) {
	  c = sgrid->e2c[es][1];
	  n = sgrid->e2e[es][1];
	}
	c1 = sgrid->c2cc[c];
	printf("e=%d es=%d c=[%f %f] e=[%f %f]\n",e, es, sgrid->cellx[c], sgrid->cellx[c],
	       m->xloc[m->eloc[0][c1][n]], m->yloc[m->eloc[0][c1][n]]);
      }
      printf("e=%d es=%d\n", e, es);
    }
    for (j = 0; j <= 1; j++) {
      int found = 0;
      v = sgrid->e2v[e][j];
      vs = sgrid->m2dv[v];
      if (verbose) printf("  j=%d v=%d\n",j, v);
      for (vv = 1; vv <= sgrid->nve[vs]; vv++) {
	e1 = sgrid->v2e[v][vv];
	if (verbose) printf("    vs=%d vv=%d e1=%d\n",vs, vv,e1);
	if (e1 == e) {
	  found = 1;
	  break;
	}
      }
      if (!found) hd_quit("vertex map error: vertex %d on edge %d (direction %d) does not remap to edge.\n", v, e, j);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the gridded mappings if required                            */
  /*
  sgrid->s2k = i_alloc_1d(sgrid->szc);
  memset(sgrid->s2k, 0, sgrid->szc * sizeof(int));
  for (k = nz - 1; k >= 0; k--) {
    for (cc = 1; cc <= ns2; cc++) {
      if (k < kbot[cc]) continue;
      c = vlm[cc2s[cc]][k];
      sgrid->s2k[c] = k;
    }
  }
  for (cc = 1; cc <= sgrid->nbpt; cc++) {
    c = sgrid->bpt[cc];
    c2 = sgrid->bin[cc];
    sgrid->s2k[c] = sgrid->s2k[c2];
  }
  */
  /*-----------------------------------------------------------------*/
  /* Get the centre to vertex map                                    */
  sgrid->c2v = i_alloc_2d(sgrid->szc, sgrid->npem+1);
  for (cc = 1; cc <= sgrid->n3_t; cc++) {
    int nn, jn, found;
    c = sgrid->w3_t[cc];
    cs = sgrid->m2d[c];
    npe = sgrid->npe[cs];
    jj = 1;
    for (n = 1; n <= npe; n++) {
      e = sgrid->c2e[n][c];
      if (sgrid->wgst[c] && e == 0) continue;
      jn = 0;
      if (sgrid->eSc[n][cs] == -1) jn = 1;
      for (j = 0; j <= 1; j++) {
        v = sgrid->e2v[e][jn];
	if (!sgrid->wgst[c] && v <= 0)
	  if (params->us_type & US_IJ)
	    hd_quit("c2v: Zero vertex: cell %d (cs=%d, cc=%d [%d %d %d]), edge %d(%d:%d) direction %d\n", c, sgrid->m2d[c], sgrid->c2cc[sgrid->m2d[c]], 
		    m->iloc[sgrid->c2cc[sgrid->m2d[c]]],
		    m->jloc[sgrid->c2cc[sgrid->m2d[c]]], sgrid->s2k[c], e, 
		    sgrid->m2de[e], n, j);
	  else
	    hd_quit("c2v: Zero vertex: cell %d (cs=%d, cc=%d [layer %d]), edge %d(%d:%d) direction %d\n", c, sgrid->m2d[c], sgrid->c2cc[sgrid->m2d[c]], 
		    sgrid->s2k[c], e, sgrid->m2de[e], n, j);
	found = 0;
	for (nn = 1; nn <= npe; nn++) {
	  if (sgrid->c2v[nn][c] == v) {
	    found = 1;
	    break;
	  }
	}
	if (!found) {
	  if(jj > npe)
	    if (params->us_type & US_IJ)
	      hd_quit("c2v: jj > %d: cell %d (cs=%d, cc=%d [%d %d %d]), edge %d(%d:%d) direction %d\n", 
		      npe, c, sgrid->m2d[c], sgrid->c2cc[sgrid->m2d[c]], 
		      m->iloc[sgrid->c2cc[sgrid->m2d[c]]],
		      m->jloc[sgrid->c2cc[sgrid->m2d[c]]], sgrid->s2k[c], e, 
		      sgrid->m2de[e], n, j);
	    else
	      hd_quit("c2v: jj > %d: cell %d (cs=%d, cc=%d [layer %d]), edge %d(%d:%d) direction %d\n", 
		      npe, c, sgrid->m2d[c], sgrid->c2cc[sgrid->m2d[c]], 
		      sgrid->s2k[c], e, sgrid->m2de[e], n, j);
	  sgrid->c2v[jj][c] = v;
	  jj++;
	}
	jn = (jn == 0) ? 1 : 0;
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the vertex to centre map. This is only mapped for wet       */
  /* centres; i.e. vertices do not map to ghost centres. This is     */
  /* consistent with the partial dual areas (dualareap) for each     */
  /* vertex summing to the dual area.                                */
  sgrid->nvc = i_alloc_1d(sgrid->szvS);
  memset(sgrid->nvc, 0, sgrid->szvS * sizeof(int));
  sgrid->nvcm = 0;

  /* Count the centres adjacent to a vertex                          */
  /* Note: 7 is the maximum number of centres expected to be         */
  /* adjacent to a vertex (6 centres for trianges).                  */
  v2c = i_alloc_2d(7, sgrid->szv);  
  for (vv = 1; vv <= sgrid->n2_e2; vv++) {
    int found, nn;
    v = sgrid->w2_e2[vv];
    sgrid->nvc[v] = 0;
    for (n = 1; n <= sgrid->nve[v]; n++) {
      e = sgrid->v2e[v][n];
      for (j = 0; j <= 1; j++) {
	c = sgrid->e2c[e][j];
	/* If ghost cells are to be included in the v2c mapping,     */
	/* then comment the next line (also in next loop).           */
	if (sgrid->wgst[c]) continue;
	found = 0;
	for (nn = 1; nn <= sgrid->nvc[v]; nn++) {
	  if (v2c[v][nn] == c) {
	    found = 1;
	    break;
	  }
	}
	if (!found && e) {
	  sgrid->nvc[v]++;
	  v2c[v][sgrid->nvc[v]] = c;
	  if (sgrid->nvc[v] > sgrid->nvcm) sgrid->nvcm = sgrid->nvc[v];
	}
      }
    }
  }
  i_free_2d(v2c);

  /* Set the mappings                                                */
  sgrid->v2c = i_alloc_2d(sgrid->nvcm+1, sgrid->szv);
  for (vv = 1; vv <= sgrid->n3_e2; vv++) {
    v = sgrid->w3_e2[vv];
    vs = sgrid->m2dv[v];
    for (n = 1; n <= sgrid->nvc[vs]; n++)
      sgrid->v2c[vv][n] = 0;
  }
  for (vv = 1; vv <= sgrid->n3_e2; vv++) {
    int found, nn;
    v = sgrid->w3_e2[vv];
    vs = sgrid->m2dv[v];
    jj = 1;
    for (n = 1; n <= sgrid->nve[vs]; n++) {
      e = sgrid->v2e[v][n];
      for (j = 0; j <= 1; j++) {
	c = sgrid->e2c[e][j];
	if (sgrid->wgst[c]) continue;
	found = 0;
	for (nn = 1; nn <= sgrid->nvc[vs]; nn++) {
	  if (sgrid->v2c[v][nn] == c) {
	    found = 1;
	    break;
	  }
	}
	if (!found && e) {
	  sgrid->v2c[v][jj] = c;
	  jj++;
	}
      }
    }
  }

  /* Get the vertical maps for vertices                              */
  /*
  for (ee = 1; ee <= sgrid->n3_e1; ee++) {
    int zm1;
    e = sgrid->w3_e1[ee];
    zm1 = sgrid->zm1e[e];
    if (e != sgrid->bot_e1[ee]) {
      v1 = sgrid->e2v[e][0];
      v2 = sgrid->e2v[zm1][0];
      sgrid->zm1v[v1] = v2;
      sgrid->zp1v[v2] = v1;
      v1 = sgrid->e2v[e][1];
      v2 = sgrid->e2v[zm1][1];
      sgrid->zm1v[v1] = v2;
      sgrid->zp1v[v2] = v1;
    }
  }

  for (vv = 1; vv <= sgrid->n2_e2; vv++) {
    v = vs = sgrid->w2_e2[vv];
    while (v != sgrid->zm1v[v]) {
      sgrid->m2dv[v] = vs;
      v = sgrid->zm1v[v];
    }
  }
  */

  /* Set the tangential ghost edges to map to vertices               */
  /*
  for (ee = sgrid->n2_e1 + 1; ee <= sgrid->x2_e1; ee++) {
    int c1, g1, c2, g2, v1, v2, i1, i2;
    e = sgrid->w2_e1[ee];
    g1 = sgrid->e2c[e][0];
    c1 = sgrid->wgst[g1];
    g2 = sgrid->e2c[e][1];
    c2 = sgrid->wgst[g2];
    if (c1 != c2) {
      int found = 0;
      for (i1 = 1; i1 <= sgrid->npe[c1]; i1++) {
	v1 = sgrid->c2v[i1][c1];
	if (found) break;
	for (i2 = 1; i2 <= sgrid->npe[c2]; i2++) {
	  v2 = sgrid->c2v[i2][c2];
	  if (found) break;
	  if (v1 == v2) {
	    for (n = 1; n <= sgrid->nve[v1]; n++) {
	      if (sgrid->nvc[v1] < sgrid->nve[v1]) {
		sgrid->e2v[e][0] = v1;
		found = 1;
		break;
	      }
	    }
	  }
	}
      }
    }
  }
  */

  if (DEBUG("init_m"))
    dlog("init_m", "\nVertex mappings created OK\n");

  /* Reset the bottom edge to be the edge adjacent to the shallowest */
  /* of the cell centres sharing the edge.                           */
  for (ee = 1; ee <= sgrid->n2_e1; ee++) {
    e = sgrid->w2_e1[ee];
    c1 = sgrid->e2c[e][0];
    j = sgrid->e2e[e][0];
    c2 = sgrid->e2c[e][1];

    /* Set any ghost cells to wet cells */
    if (sgrid->wgst[c1]) {
      c1 = c2;
      j = sgrid->e2e[e][1];
    }
    if (sgrid->wgst[c2]) c2 = c1;
    /* Find the shallowest cell */
    if (kbot[sgrid->c2cc[c2]] > kbot[sgrid->c2cc[c1]]) {
      c1 = c2;
      j = sgrid->e2e[e][1];
    }
    while (c1 != sgrid->zm1[c1]) {
      c1 = sgrid->zm1[c1];
    }
    c1 = sgrid->zp1[c1];
    sgrid->bot_e1[ee] = sgrid->c2e[j][c1];
  }

  /* Reset the vertical maps                                         */
  /* Re-included on 17.03 : OK check required.                       */

  for (ee = 1; ee <= sgrid->n2_e1; ee++) {
    e = es = sgrid->w2_e1[ee];
    eb = sgrid->bot_e1[ee];
    while (e != sgrid->zm1e[e])
      e = sgrid->zm1e[e];
    sgrid->zm1e[eb] = e;
    sgrid->zp1e[e] = eb;
  }

  /* Get the vertex geographic locations by finding the edges        */
  /* that emanate from a vertex and have the same locx and locy for  */
  /* each edge and index 0 or 1.                                     */
  sgrid->gridx = d_alloc_1d(sgrid->szvS);
  sgrid->gridy = d_alloc_1d(sgrid->szvS);
  for (vv = 1; vv <= sgrid->n2_e2; vv++) {
    int nn;
    double lat1, lon1, lat2, lon2;
    int found = 0;
    v = sgrid->w2_e2[vv];
    for (j = 1; j <= sgrid->nve[v]; j++) {
      e1 = sgrid->v2e[v][j];
      if (!e1) continue;
      for (n = 0; n <= 1; n++) {
        lon1 = locx[e1][n]; lat1 = locy[e1][n];
        for (jj = 1; jj <= sgrid->nve[v]; jj++) {
	  if (jj == j) continue;
	  e2 = sgrid->v2e[v][jj];
	  if (!e2) continue;
	  for (nn = 0; nn <= 1; nn++) {
	    lon2 = locx[e2][nn]; lat2 = locy[e2][nn];
	    if (lat1 == lat2 && lon1 == lon2) {
	      sgrid->gridx[v] = lon1;
	      sgrid->gridy[v] = lat1;
	      found = 1;
	      break;
	    }
	  }
	  if (found) break;
	}
	if (found) break;
      }
      if (found) break;
    }
  }

  /* Get the edge geographic locations. These are the means of the   */
  /* edges's verticies.                                              */
  sgrid->u1x = d_alloc_1d(sgrid->szeS);
  sgrid->u1y = d_alloc_1d(sgrid->szeS);
  for (ee = 1; ee <= sgrid->n2_e1; ee++) {
    e = sgrid->w2_e1[ee];
    v1 = sgrid->e2v[e][0];
    v2 = sgrid->e2v[e][1];
    sgrid->u1x[e] = 0.5 * (sgrid->gridx[v1] + sgrid->gridx[v2]);
    sgrid->u1y[e] = 0.5 * (sgrid->gridy[v1] + sgrid->gridy[v2]);
  }

  /* Get the vertex sign and index vectors                           */
  sgrid->vIc = i_alloc_2d(sgrid->szcS, npem + 1);
  sgrid->eSv = i_alloc_2d(sgrid->szv, sgrid->nvem + 1);
  for (cc = 1; cc <= sgrid->b2_t; cc++) {
    c = sgrid->w2_t[cc];
    npe = sgrid->npe[c];
    for (n = 1; n <= npe; n++) {
      v = sgrid->c2v[n][c];
      for (vv = 1; vv <= sgrid->nvc[v]; vv++) {
	if (c == sgrid->v2c[v][vv])
	  sgrid->vIc[n][c] = vv;
      }
    }
  }
  for (vv = 1; vv <= sgrid->n3_e2; vv++) {
    v = sgrid->w3_e2[vv];
    vs = sgrid->m2dv[v];
    for (n = 1; n <= sgrid->nve[vs]; n++) {
      e = sgrid->v2e[v][n];
      if (e) {
	if (v == sgrid->e2v[e][0])
	  sgrid->eSv[n][v] = 1;    /* v is on the left of u1  */
	else
	  sgrid->eSv[n][v] = -1;   /* v is on the right of u1 */
      }
    }
  }

  /* Get the area of the dual cell                                   */
  sgrid->dualarea = d_alloc_1d(sgrid->szvS);
  sgrid->dualareap = d_alloc_2d(sgrid->nvcm + 1, sgrid->szvS);
  memset(sgrid->dualarea, 0, (sgrid->szvS) * sizeof(double));
  /* Initialise                                                      */
  for (vv = 1; vv <= sgrid->n2_e2; vv++) {
    v = sgrid->w2_e2[vv];
    sgrid->dualarea[v] = 0.0;
    for (j = 1; j <= sgrid->nvc[v]; j++)
      sgrid->dualareap[v][j] = 0.0;
  }
  /* Get the partial dual areas                                    */
  for (cc = 1; cc <= sgrid->b2_t; cc++) {
    int nn;
    c = sgrid->w2_t[cc];
    c1 = sgrid->c2cc[c];
    npe = sgrid->npe[c];
    for (n = 1; n <= npe; n++) {
      v = sgrid->c2v[n][c];
      j = sgrid->vIc[n][c];
      for (nn = 1; nn <= npe; nn++) {
	e = sgrid->c2e[nn][c];
	if (sgrid->e2v[e][0] == v || sgrid->e2v[e][1] == v) {
	  sgrid->dualareap[v][j] +=  0.25 * sgrid->h1au1[e] * params->h2[c1][nn-1];
	}
      }
    }
  }

  /* Get the dual areas                                            */
  for (cc = 1; cc <= sgrid->b2_t; cc++) {
    c = sgrid->w2_t[cc];
    npe = sgrid->npe[c];
    for (n = 1; n <= npe; n++) {
      v = sgrid->c2v[n][c];
      j = sgrid->vIc[n][c];
      sgrid->dualarea[v] += sgrid->dualareap[v][j];
    }
  }

  /* Get the edges surrounding an edge map and weights.              */
  /* Uses Thuburn et al, (2009) J. Comp. Phys., 228 Eq. 33.          */
  /* Note that k x u rotates the velocity vector by 90 degrees in    */
  /* clockwise direction (Ringler: Momentum, vorticity and transport */
  /* Considerations in the design of a finite-volume dynamical core) */
  /* p31. Here we need the tangential velocity component in the      */
  /* direction -k x u (Thuburn (2009), p8324).                       */
  sgrid->nee = i_alloc_1d(sgrid->sze);
  sgrid->eSe = i_alloc_2d(sgrid->sze, sgrid->neem + 1);
  sgrid->wAe = d_alloc_2d(sgrid->sze, sgrid->neem + 1);
  sgrid->wSe = d_alloc_2d(sgrid->szeS, sgrid->neem + 1);
  dume = i_alloc_1d(sgrid->sze);
  len = d_alloc_1d(sgrid->szcS);
  rw = d_alloc_2d(sgrid->szcS, npem + 1);
  memset(len, 0, sgrid->szcS * sizeof(double));

  for (ee = 1; ee <= sgrid->n2_e1; ee++) {
    e = sgrid->w2_e1[ee];
    c1 = sgrid->e2c[e][0];
    c2 = sgrid->e2c[e][1];
    sgrid->nee[e] = sgrid->npe[c1] + sgrid->npe[c2] - 1;
  }
  for (ee = 1; ee <= sgrid->n3_e1; ee++) {
    e = sgrid->w3_e1[ee];
    dume[e] = 0;
    for (j = 1; j <= sgrid->neem; j++) 
      sgrid->wAe[j][e] = 0.0;
  }
  for (j = 1; j <= sgrid->neem; j++) 
    sgrid->wAe[j][0] = 0.0;
  /* Perimeter of each cell                                          */
  for (cc = 1; cc <= sgrid->b2_t; cc++) {
    c = sgrid->w2_t[cc];
    for (n = 1; n <= sgrid->npe[c]; n++) {
      e = sgrid->c2e[n][c];
      len[c] += sgrid->h1au1[e];
      rw[n][c] = 0.0;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Compute the cell weights used in the Thuburn et al. (2009)      */
  /* algorithm.                                                      */
  /* Thuburn (2009) Section 4 equally distributes rw across vertices */
  /* for quad and hex examples (i.e. rtype=0).                       */
  /* Ringler et al. (2010) p 3073 states rw should be the area of    */
  /* intersection between primal and dual mesh normalized by the     */
  /* primal mesh area (i.e. rtype=2).                                */
  /* rtype = 0 : Equally distributed over vertices                   */
  /* rtype = 1 : Scaled by edge length                               */
  /* rtype = 2 : Scaled by partial cell area                         */
  /* rtype = 3 : Includes contributions of all cells common to a     */
  /*             vertex.                                             */
  for (cc = 1; cc <= sgrid->b2_t; cc++) {
    double d1 = 0.0, d2 = 0.0;
    c = sgrid->w2_t[cc];
    npe = sgrid->npe[c];
    for (j = 1; j <= npe; j++) {
      e = sgrid->c2e[j][c];
      if (rtype == 0)
	rw[j][c] = 1.0 / (double)npe;
      else if (rtype == 1) {
	rw[j][c] = 0.5 * sgrid->h1au1[e] / len[c];
	n = (j == 1) ? npe : j - 1;
	e = sgrid->c2e[n][c];
	rw[j][c] += 0.5 * sgrid->h1au1[e] / len[c];
      } else if (rtype == 2) {
	c1 = sgrid->c2cc[c];
	rw[j][c] += (0.25 * sgrid->h1au1[e] * params->h2[c1][j-1]) / sgrid->cellarea[c];
	n = (j == 1) ? npe : j - 1;
	e = sgrid->c2e[n][c];
	rw[j][c] += (0.25 * sgrid->h1au1[e] * params->h2[c1][j-1]) / sgrid->cellarea[c];
      } else if (rtype == 3) {
	v = sgrid->c2v[j][c];
	for (vv = 1; vv <= sgrid->nvc[v]; vv++) {
	  c1 = sgrid->v2c[v][vv];
	  rw[j][c] += 1.0 / ((double)sgrid->npe[c1] * (double)sgrid->nvc[v]);
	}
      }
      d1 += rw[j][c];
    }
    /* Normalize sum of rw to 1                                      */
    for (j = 1; j <= npe; j++) rw[j][c] /= d1; 
  }

  /* Loop over all cell centres and consider edges surrounding that  */
  /* centre individually.                                            */
  for (cc = 1; cc <= sgrid->b3_t; cc++) {
    double u, vt, d1;
    int de = 0;             /* Debugging edge index                  */
    double rsum;            /* Sum of rw around the cell             */

    c = sgrid->w3_t[cc];
    cs = sgrid->m2d[c];
    npe = sgrid->npe[cs];   /* Number of edges for this centre       */

    /* Loop over edges                                               */
    for (n = 1; n <= npe; n++) {
      e = sgrid->c2e[n][c];    /* Cell edge                          */

      /* The first weight for this edge is always zero; see Thuburn  */
      /* (2009) p8324.                                               */
      if (dume[e] == 0) {      
	dume[e]++;
	sgrid->eSe[dume[e]][e] = e;
	sgrid->wAe[dume[e]][e] = 0.0;
	vt = 0.0;
      }

      /* Loop over edges in a clockwise direction until the original */
      /* edge e is encountered.                                      */
      ee = (n + 1 > npe) ? 1 : n + 1;

      /* Find the last vertex traversing from e1 to e. This is the   */
      /* first vertex from e to e1, or the vertex at the current ee. */
      v = sgrid->c2v[ee][c];
      vs = sgrid->m2dv[v];

      /* Find the index of vertex v corresponding to edge e          */
      for (vv = 1; vv <= sgrid->nve[vs]; vv++) {
	if (e == sgrid->v2e[v][vv]) break;
      }

      /* Initialize the sum of vertex weights to rw at vertex v      */
      rsum = rw[ee][cs];

      while (ee != n) {
	e1 = sgrid->c2e[ee][c];	

	/* Include the edge e1 in the array and compute the weight.  */
	/* Note: the weights of Thuburn (2008) deliver tangential    */
	/* velocities in the -k x ne direction (p8324). i.e. rotated */
	/* clockwise from the normal component. We take the negative */
	/* of this to deliver the tangential velocity in the k x ne  */
	/* direction (rotated anti-clockwise). This means we must    */
	/* add the nonlinear Coriolis tendency (Ringler (2010)       */
	/* Eq. 24) in nonlin_coriolis_3d() to restore the correct    */
	/* direction of tangential velocity in this term (-k x u).   */
	dume[e]++;
	sgrid->eSe[dume[e]][e] = e1;
	if (sgrid->eSv[vv][v])
	  sgrid->wAe[dume[e]][e] = -((rsum - 0.5) * (double)sgrid->eSc[ee][cs]) / (double)sgrid->eSv[vv][v];

	if (de == e) {
	  printf("ee=%d eoe=%d rsum=%f %f : eSc=%d eSv=%d sgn=%f : w=%f\n",dume[e],e1,rsum,rsum-0.5,sgrid->eSc[ee][cs],sgrid->eSv[vv][v],-(double)sgrid->eSc[ee][cs] / (double)sgrid->eSv[vv][v],sgrid->wAe[dume[e]][e]);
	}

	/* Move on to the next edge; increment the edge counter, and */
	/* add rw to sum of vertex weights.                          */
	ee = (ee + 1 > npe) ? 1 : ee + 1;
	rsum += rw[ee][cs];
      }
    }
  }


  /* Oldcode
  for (ee = 1; ee <= sgrid->n3_e1; ee++) {
    e = sgrid->w3_e1[ee];
    dume[e] = 0;
    for (j = 1; j <= sgrid->neem; j++) 
      sgrid->wAe[j][e] = 0.0;
  }
  for (cc = 1; cc <= sgrid->b3_t; cc++) {
    int nn;
    double r;
    double rsum;
    c = sgrid->w3_t[cc];
    cs = sgrid->m2d[c];
    npe = sgrid->npe[cs];
    r = 1.0 / (double)npe;
    for (n = 1; n <= npe; n++) {
      e = sgrid->c2e[n][c];
      if (dume[e] == 0) {
	dume[e]++;
	sgrid->eSe[dume[e]][e] = e;
	sgrid->wAe[dume[e]][e] = 0.0;
      }
      rsum = r;
      e2 = e;
      nn = (n + 1 > npe) ? 1 : n + 1;
      while (nn != n) {
	e1 = sgrid->c2e[nn][c];	
	for (j = 0; j <= 1; j++) {
	  v1 = sgrid->e2v[e2][j];
	  for (jj = 0; jj <= 1; jj++) {
	    v2 = sgrid->e2v[e1][jj];
	    if (v1 == v2) v = v1;
	  }
	}
	vs = sgrid->m2dv[v];
	for (j = 1; j <= sgrid->nve[vs]; j++) {
	  if (e1 == sgrid->v2e[v][j]) break;
	}
	dume[e]++;
	sgrid->eSe[dume[e]][e] = e1;
	if (sgrid->eSv[j][v])
	  sgrid->wAe[dume[e]][e] = (rsum - 0.5) * sgrid->eSc[n][cs] / sgrid->eSv[j][v];

	rsum += r;
	e2 = e1;
	nn = (nn + 1 > npe) ? 1 : nn + 1;
      }
    }
  }
  */

  i_free_1d(dume);
  d_free_1d(len);
  d_free_2d(rw);

  /*
  for (cc = 1; cc <= sgrid->b2_t; cc++) {
  double d1, u, vt;
    c = sgrid->w2_t[cc];
    for (n = 1; n <= sgrid->npe[c]; n++) {
      if (n==3&&c == sgrid->cc2s[17422]) {
	FILE *op = fopen("e.site", "w");
	e = sgrid->c2e[n][c];
	for (ee = 1; ee <= sgrid->nee[e]; ee++) {
	  es = sgrid->eSe[ee][e];	  
	  d1 = sgrid->thetau1[sgrid->m2de[es]] - sgrid->thetau1[sgrid->m2de[e]];
	  u = 1.0 * cos(d1) + 0.0 * sin(d1);
	  vt += u * sgrid->wAe[ee][e];

	  printf("6s c=%d: n=%d e=%d eoe%d=%d : %f %f %f\n",c,n,e,ee,es,sgrid->wAe[ee][e],u,vt);
	  fprintf(op, "%f %f e%d\n",sgrid->u1x[es],sgrid->u1y[es],ee);
	}
	fclose(op);
      }
    }
  }

  for (ee = 1; ee <= sgrid->n2_e1; ee++) {
    e = sgrid->w2_e1[ee];
    for (n = 1; n <= sgrid->nee[e]; n++)

      if(e==7)
	printf("%d %d %d %f\n",e,n,sgrid->eSe[n][e],sgrid->wAe[n][e]);
  }
  */

  if (DEBUG("init_m"))
    dlog("init_m", "\nDual mappings created OK\n");

  /*-----------------------------------------------------------------*/
  /* Get the gridded mappings if required                            */
  sgrid->e2ijk = i_alloc_1d(sgrid->sze);
  for (ee = 1; ee <= sgrid->n3_e1; ee++) {
    e = sgrid->w3_e1[ee];
    sgrid->e2ijk[e] = NOTVALID;
  }

  /* Edge maps to cell centre to retrieve (i,j,k)                    */
  for (ee = 1; ee <= sgrid->n3_e1; ee++) {
    e = sgrid->w3_e1[ee];
    c = sgrid->e2c[e][0];
    if (sgrid->wgst[c])
      c = sgrid->e2c[e][1];
    sgrid->e2ijk[e] = c;
  }

  /* If grid information is not supplied, then use the cell index    */
  /* and surface coordinate for s2i and s2j respectively.            */
  sgrid->s2i = i_alloc_1d(sgrid->szc);
  sgrid->s2j = i_alloc_1d(sgrid->szc);
  for (cc = 1; cc <= sgrid->b2_t; cc++) {
    c = sgrid->w2_t[cc];
    sgrid->s2i[c] = sgrid->c2cc[c];
    sgrid->s2j[c] = sgrid->m2d[c];
  }
  if (params->us_type & US_IJ) {
    int *e2i, *e2j;
    sgrid->nce1 = m->nce1;
    sgrid->nce2 = m->nce2;
    sgrid->nfe1 = m->nce1 + 1;
    sgrid->nfe2 = m->nce2 + 1;
    sgrid->map = (unsigned long ***)l_alloc_3d(m->nce1+1, m->nce2+1, nz);
    flag = (unsigned long ***)l_alloc_3d(sgrid->nce1 + 1, sgrid->nce2 + 1, nz);
    topo = d_alloc_2d(sgrid->nce1, sgrid->nce2);
    /* Initialize                                                      */
    for (j = 0; j < sgrid->nce2; j++)
      for (i = 0; i < sgrid->nce1; i++) {
	if (params->topo)
	  topo[j][i] = params->topo[j][i];
	else
	  topo[j][i] = NaN;
	for (k = nz - 1; k >= 0; k--) {
	  if (params->flag)
	    flag[k][j][i] = params->flag[j][i];
	  else
	    flag[k][j][i] = SOLID;
	}
      }
    for (k = nz - 1; k >= 0; k--) {
      for (cc = 1; cc <= ns2; cc++) {
	if (k < kbot[cc]) continue;
	c = vlm[cc2s[cc]][k];     /* Wet cell in layer k               */
	i = m->iloc[cc];
	j = m->jloc[cc];
	sgrid->s2i[c] = i;
	sgrid->s2j[c] = j;
	sgrid->map[k][j][i] = c;
	topo[j][i] = bathy[cc];
	flag[k][j][i] = flg[k][cc];
      }
    }
    /* Ghost cells */

    for (cc = 1; cc <= sgrid->n2_t; cc++) {
      c = sgrid->w2_t[cc];
      c1 = sgrid->wgst[c];
      if (c1) {
	if (sgrid->c2c[1][c1] == c) {
	  i = sgrid->s2i[c1] - 1;
	  j = sgrid->s2j[c1];
	} else if (sgrid->c2c[3][c1] == c) {
	  i = sgrid->s2i[c1] + 1;
	  j = sgrid->s2j[c1];
	} else if (sgrid->c2c[4][c1] == c) {
	  i = sgrid->s2i[c1];
	  j = sgrid->s2j[c1] - 1;
	} else if (sgrid->c2c[2][c1] == c) {
	  i = sgrid->s2i[c1];
	  j = sgrid->s2j[c1] + 1;
	}
	sgrid->s2i[c] = i;
	sgrid->s2j[c] = j;
	sgrid->s2k[c] = sgrid->s2k[c1];
	/*
	sgrid->s2i[c] = NOTVALID;
	sgrid->s2j[c] = NOTVALID;
	sgrid->s2k[c] = NOTVALID;
	*/
	sgrid->map[sgrid->s2k[c]][sgrid->s2j[c]][sgrid->s2i[c]] = c;
      }
    }

    /* Edge maps to (i,j)                                             */
    e2i = i_alloc_1d(sgrid->szeS);
    e2j = i_alloc_1d(sgrid->szeS);
    for (ee = 1; ee <= sgrid->n2_e1; ee++) {
      e = sgrid->w2_e1[ee];
      c = sgrid->e2c[e][0];
      j = sgrid->e2e[e][0];
      e2i[e] = sgrid->s2i[c];
      e2j[e] = sgrid->s2j[c];
      if (sgrid->wgst[c] && j == 1) {
	e2i[e] = sgrid->s2i[sgrid->e2c[e][1]] + 1;
	e2j[e] = sgrid->s2j[sgrid->e2c[e][1]];
      }
      if (sgrid->wgst[c] && j == 4) {
	e2i[e] = sgrid->s2i[sgrid->e2c[e][1]];
	e2j[e] = sgrid->s2j[sgrid->e2c[e][1]] + 1;
      }
    }

    /* Vertex maps to (i,j)                                           */
    if (sgrid->npem == 4)
      vertex_map_4(sgrid);

    i_free_1d(e2i);
    i_free_1d(e2j);
  }

  /*-----------------------------------------------------------------*/
  /* Check the mappings                                              */
  for (cc = 1; cc <= sgrid->b2_t; cc++) {
    c = sgrid->w2_t[cc];
    npe = sgrid->npe[c];
    for (j = 1; j <= npe; j++) {
      e = sgrid->c2e[j][c];
      if (sgrid->e2c[e][0] != c && sgrid->e2c[e][1] != c) {
	printf("Edge map error @ c=%d e2c(%d,0:1) = (%d,%d) j = %d\n", 
	       c, e, sgrid->e2c[e][0], sgrid->e2c[e][1], j);
      }
      if (sgrid->e2e[e][0] != j && sgrid->e2e[e][1] != j) 
	printf("Edge map error @ c=%d: e2e(%d,0:1) = (%d,%d), j = %d\n", 
	       c, e, sgrid->e2e[e][0], sgrid->e2e[e][1], j);
    }
  }
  if (DEBUG("init_m"))
    dlog("init_m", "\nMappings checked OK\n");

  /*-----------------------------------------------------------------*/
  /* Save the locations of sparse cells on open boundaries to the    */
  /* boundary arrays.                                                */
  /* The open boundary location vectors are:                         */
  /* obc_t = cell centered sparse location of the open boundary      */
  /* oi1_t = one cell into the interior from obc_t                   */
  /* oi2_t = two cells into the interior from obc_t                  */
  /* obc_e1 = edge centered sparse location of the open boundary     */
  /* obc_e2 = cell center corresponding to the normal OBC edge       */
  /* oi1_e1 = one cell into the interior from obc_e1                 */
  /* oi2_e1 = two cells into the interior from obc_e1                */
  /* e2c_e1 = map from edge boundary index ee (for obc_e1) to cell   */
  /*          index cc (for obc_t).                                  */

  /* Allocate memory for the open boundary data structure            */
  if (sgrid->nobc)
    sgrid->open = (open_bdrys_t **)malloc(sizeof(open_bdrys_t *) * sgrid->nobc);

  /* Allocate memory for the boundary structures                     */
  for (n = 0; n < sgrid->nobc; n++) {

    sgrid->open[n] = OBC_alloc();

    sgrid->open[n]->ntr = params->ntr;
    sgrid->open[n]->npts = m->npts[n];
    sgrid->open[n]->iloc = i_alloc_1d(m->npts[n]+1);
    sgrid->open[n]->bgz = 0;
    for (cc = 1; cc <= sgrid->open[n]->npts; cc++)
      sgrid->open[n]->iloc[cc] = m->loc[n][cc];
    sgrid->open[n]->no2_t = 0;
    sgrid->open[n]->no3_t = 0;
    sgrid->open[n]->no2_e1 = 0;
    sgrid->open[n]->no3_e1 = 0;
    sgrid->open[n]->to2_e1 = 0;
    sgrid->open[n]->to3_e1 = 0;
    sgrid->open[n]->bgz = laus;
    copy_OBC_conds(params->open[n], sgrid->open[n], n, params->trinfo_3d);
    /* All unstructured OBCs are of type U1BDRY                      */
    sgrid->open[n]->type = U1BDRY;
  }

  /* Set the OBC mask */
  /*
  memset(mask, 0, sgrid->szc * sizeof(int));
  for (cc = 1; cc <= ns2; cc++) {
    c = cc2s[cc];
    mask[c] = -1;
    for (n = 0; n < m->nobc; n++) {
      for (c1 = 1; c1 <= m->npts[n]; c1++) {
	if (cc == m->loc[n][c1]) {
	  mask[c] = n;
	}
      }
    }
  }
  */

  /* Count the cell centre locations */
  if (m->nobc) {
    i_free_1d(mask);
    mask = i_alloc_1d(m->nobc);
  }
  for (cc = sgrid->v2_t + 1; cc <= sgrid->b2_t; cc++) {
    c = cs = sgrid->w2_t[cc];
    /*
    n = mask[c];
    n = obc_num(sgrid->c2cc[c], m->nobc, m->npts, m->loc);
    */
    /* Several OBCs may share the same cell centre; capture all      */
    /* these instances in the loop below.                            */
    obc_num_all(sgrid->c2cc[c], m->nobc, m->npts, m->loc, mask);
    for(tn = 0; tn < m->nobc; tn++) {
      if ((n = mask[tn]) >= 0) {
	/*if (n >= 0) {*/
	sgrid->open[n]->no2_t++;
	sgrid->open[n]->no3_t++;
	c = sgrid->zm1[cs];
	cb = sgrid->bot_t[cc];
	while (c <= cb) {
	  sgrid->open[n]->no3_t++;
	  c = sgrid->zm1[c];
	}
      }
    }
  }

  /* Count the edge locations                                        */
  /* ocodec = 1 if e2c[0] is a (boundary) ghost cell (e2c[1] is wet) */
  /* ocodec = 0 if e2c[1] is a (boundary) ghost cell (e2c[0] is wet) */
  /* ocodex = index pointing into the domain                         */
  /* ocodey = index pointing out of the domain                       */
  /* Note: these are different for each edge; see ceni[], ini[] and  */
  /* outi[] below.                                                   */
  imape = i_alloc_1d(sgrid->szeS);
  memset(imape, 0, sgrid->szeS * sizeof(int));
  for (cc = sgrid->v2_e1 + 1; cc <= sgrid->b2_e1; cc++) {
    open_bdrys_t *open;
    e = sgrid->w2_e1[cc];
    n = maske[e];
    open = sgrid->open[n];
    c1 = sgrid->e2c[e][0];
    c2 = sgrid->e2c[e][1];
    if (sgrid->wgst[c1]) {      
      open->ocodex = sgrid->e2e[e][0];
      open->ocodey = sgrid->e2e[e][1];
      open->ocodec = 1;
      imape[e] = open->ocodex;
    }
    if (sgrid->wgst[c2]) {
      open->ocodex = sgrid->e2e[e][1];
      open->ocodey = sgrid->e2e[e][0];
      open->ocodec = 0;
      imape[e] = open->ocodex;
    }
    if (sgrid->wgst[c1] || sgrid->wgst[c2]) {
      /* Normal velocity */
      sgrid->open[n]->no2_e1++;
      sgrid->open[n]->no3_e1++;
      e = sgrid->zm1e[e];
      cb = sgrid->bot_e1[cc];
      while (e <= cb) {
	sgrid->open[n]->no3_e1++;
	e = sgrid->zm1e[e];
      }
    } else {
      /* Tangential velocity */
      sgrid->open[n]->to2_e1++;
      sgrid->open[n]->to3_e1++;
      e = sgrid->zm1e[e];
      cb = sgrid->bot_e1[cc];
      while (e <= cb) {
	sgrid->open[n]->to3_e1++;
        e = sgrid->zm1e[e];
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Allocate memory for the boundary vectors                        */
  /* Open boundary vectors                                           */
  for (n = 0; n < sgrid->nobc; n++) {
    open_bdrys_t *open = sgrid->open[n];
    c1 = open->no3_e1 + open->to3_e1 + 1;
    open->obc_t = i_alloc_1d(open->no3_t + 1);
    open->ogc_t = i_alloc_1d(open->no3_e1 + 1);
    open->obc_e1 = i_alloc_1d(c1);
    open->obc_e2 = i_alloc_1d(c1);
    open->bot_t = i_alloc_1d(open->no2_t + 1);
    open->nepc = d_alloc_1d(open->no3_t + 1);
    open->ceni = i_alloc_1d(open->no3_e1 + 1);
    open->dir = i_alloc_1d(open->no3_e1 + 1);
    open->ini = i_alloc_1d(open->no3_e1 + 1);
    open->outi = i_alloc_1d(open->no3_e1 + 1);
    open->inc = i_alloc_1d(open->no3_t + 1);
    open->bec = i_alloc_2d(open->no3_t + 1, npem + 1);
    open->bcc = i_alloc_2d(open->no3_t + 1, npem + 1);
    open->olap = i_alloc_1d(open->no3_t + 1);
    open->nmape = (int **)p_alloc_1d(open->no3_e1+1);
    open->omape = (int **)p_alloc_1d(open->no3_e1+1);
    memset(open->nepc, 0, (open->no3_t + 1) * sizeof(double));

    /* Sparse locations one cell into the interior */
    open->oi1_t = i_alloc_1d(open->no3_t + 1);
    open->oi1_e1 = i_alloc_1d(c1);

    /* Sparse locations two cells into the interior */
    open->oi2_t = i_alloc_1d(open->no3_t + 1);
    open->oi2_e1 = i_alloc_1d(c1);
    
    /* Sparse locations for cyclic boundary conditions */
    open->cyc_t = i_alloc_1d(open->no3_t + 1);
    open->cyc_e1 = i_alloc_1d(c1);
    open->e2c_e1 = i_alloc_1d(c1);

    /* Dummies                                  */
    /*
    open->dum = d_alloc_1d(open->no3_t + 1);
    open->i1 = i_alloc_1d(open->no3_t + 1);
    */
    /* Reset the vector sizes for filling below */
    /* Tangential velocity positions */
    open->to3_e1 = open->no3_e1 + open->to2_e1 + 1;
    open->to2_e1 = open->no3_e1 + 1;
    /* Normal velocity and cell centered positions */
    open->no3_t = open->no2_t + 1;
    open->no3_e1 = open->no2_e1 + 1;
    open->no2_t = 1;
    open->no2_e1 = 1;
  }

  /* Point the interior cell maps on the boundaries to the correct */
  /* spatial map.                                                  */
  for (n = 0; n < sgrid->nobc; n++) {
    sgrid->open[n]->nmap = sgrid->c2c[sgrid->open[n]->ocodex];
    sgrid->open[n]->omap = sgrid->c2c[sgrid->open[n]->ocodey];
    if(sgrid->npem == 4) {
      if (sgrid->open[n]->ocodex == 1 || sgrid->open[n]->ocodex == 3) {
	sgrid->open[n]->tmpp = sgrid->c2c[2];
	sgrid->open[n]->tmpm = sgrid->c2c[4];
      }
      if (sgrid->open[n]->ocodex == 2 || sgrid->open[n]->ocodex == 4) {
	sgrid->open[n]->tmpp = sgrid->c2c[3];
	sgrid->open[n]->tmpm = sgrid->c2c[1];
      }
    }
  }

  /* Set the cell centre locations                                   */
  imapc = i_alloc_1d(sgrid->szcS);
  omapc = i_alloc_1d(sgrid->szcS);
  memset(imapc, 0, sgrid->szcS * sizeof(int));
  memset(omapc, 0, sgrid->szcS * sizeof(int));
  for (cc = sgrid->v2_t + 1; cc <= sgrid->b2_t; cc++) {
    open_bdrys_t *open;
    c = cs = sgrid->w2_t[cc];
    /*
    n = mask[c];
    n = obc_num(sgrid->c2cc[c], m->nobc, m->npts, m->loc);
    */
    j = imapc[cs] = get_bind(sgrid, c, maskb);
    if (usejo)
      omapc[cs] = oedge(sgrid->npe[c], j);
    else
      omapc[cs] = jocw(sgrid, c, j);

    /* Regular quad grids U1 boundary                                */
    if (params->us_type & US_IJ) {
      if(sgrid->wgst[sgrid->c2c[1][cs]]) {
	imapc[cs] = 3;
	omapc[cs] = 1;
      }
      if(sgrid->wgst[sgrid->c2c[3][cs]]) {
	imapc[cs] = 1;
	omapc[cs] = 3;
      }
    }
    /* Regular hex grids U1 boundary                                 */
    if (params->us_type & US_HEX) {
      if(sgrid->wgst[sgrid->c2c[1][cs]]) {
	imapc[cs] = 4;
	omapc[cs] = 1;
      }
      if(sgrid->wgst[sgrid->c2c[4][cs]]) {
	imapc[cs] = 1;
	omapc[cs] = 4;
      }
    }

    /* Find the direction that maps to the interior for this centre  */
    /*
    jj = 0; i = 0;
    for (j = 1; j <= sgrid->npe[c]; j++) {
      cn = sgrid->c2c[j][c];
      if (maskb[cn]) {
	jj += j;
	i += 1;
      }
    }

    jj = omapc[cs] = (i) ? jj / i : 1;
    j = imapc[cs] = oedge(sgrid->npe[c], jj);
    */

    obc_num_all(sgrid->c2cc[c], m->nobc, m->npts, m->loc, mask);
    for(tn = 0; tn < m->nobc; tn++) {
      if ((n = mask[tn]) >= 0) {
	open = sgrid->open[n];
	/*if (n >= 0){*/
	open->obc_t[open->no2_t] = cs;
	open->oi1_t[open->no2_t] = sgrid->c2c[j][cs];
	open->oi2_t[open->no2_t] = sgrid->c2c[j][sgrid->c2c[j][cs]];
	open->no2_t++;
	c = sgrid->zm1[cs];
	cb = sgrid->bot_t[cc];
	while (c <= cb) {
	  open->obc_t[open->no3_t] = c;
	  open->oi1_t[open->no3_t] = sgrid->c2c[j][c];
	  open->oi2_t[open->no3_t] = sgrid->c2c[j][sgrid->c2c[j][c]];
	  open->no3_t++;
	  c = sgrid->zm1[c];
	}
      }
    }
  }

  for (n = 0; n < sgrid->nobc; n++) {
    open_bdrys_t *open = sgrid->open[n];
    open->no2_t--;
    open->no3_t--;
  }

  /* Set the edge locations */
  for (ee = sgrid->v2_e1 + 1; ee <= sgrid->b2_e1; ee++) {
    open_bdrys_t *open;
    e = sgrid->w2_e1[ee];
    es = sgrid->m2de[e];
    n = maske[e];
    open = sgrid->open[n];
    c1 = sgrid->e2c[e][0];
    c2 = sgrid->e2c[e][1];
    if (sgrid->wgst[c1] || sgrid->wgst[c2]) {
      /* Normal velocity */
      open->obc_e1[open->no2_e1] = e;
      if (sgrid->wgst[c1]) {
	open->obc_e2[open->no2_e1] = c2;
	open->ogc_t[open->no2_e1] = c1;
      }
      if (sgrid->wgst[c2]) {
	open->obc_e2[open->no2_e1] = c1;
	open->ogc_t[open->no2_e1] = c2;
      }
      /*
      for(cc = 1; cc <= open->no2_t; cc++)
	if(open->obc_t[cc] == open->obc_e2[open->no2_e1]) {
	  open->nepc[cc]++;
	}
      */
      ei = e2e(sgrid, e, imape[e]);
      open->oi1_e1[open->no2_e1] = ei;
      ei = e2e(sgrid, ei, imape[e]);
      open->oi2_e1[open->no2_e1] = ei;
      open->no2_e1++;
      e = sgrid->zm1e[e];
      cb = sgrid->bot_e1[ee];
      c1 = sgrid->e2c[e][0];
      c2 = sgrid->e2c[e][1];
      while (e <= cb) {
        open->obc_e1[open->no3_e1] = e;
	if (sgrid->wgst[c1]) {
	  open->obc_e2[open->no3_e1] = c2;
	  open->ogc_t[open->no3_e1] = c1;
	}
	if (sgrid->wgst[c2]) {
	  open->obc_e2[open->no3_e1] = c1;
	  open->ogc_t[open->no3_e1] = c2;
	}
	ei = e2e(sgrid, e, imape[es]);
	open->oi1_e1[open->no3_e1] = ei;
	ei = e2e(sgrid, ei, imape[es]);
	open->oi2_e1[open->no3_e1] = ei;
	open->no3_e1++;
        e = sgrid->zm1e[e];
	c1 = sgrid->e2c[e][0];
	c2 = sgrid->e2c[e][1];
      }
    } else {
      /* Tangential velocity */
      open->obc_e1[open->to2_e1] = e;
      open->obc_e2[open->to2_e1] = c1;
      ei = e2e(sgrid, e, imape[e]);
      open->oi1_e1[open->to2_e1] = ei;
      ei = e2e(sgrid, ei, imape[e]);
      open->oi2_e1[open->to2_e1] = ei;
      open->to2_e1++;
      e = sgrid->zm1e[e];
      c1 = sgrid->e2c[e][0];
      cb = sgrid->bot_e1[ee];
      while (e <= cb) {
        open->obc_e1[open->to3_e1] = e;
        open->obc_e2[open->to3_e1] = c1;
	ei = e2e(sgrid, e, imape[es]);
	open->oi1_e1[open->to3_e1] = ei;
	ei = e2e(sgrid, ei, imape[es]);
	open->oi2_e1[open->to3_e1] = ei;
	open->to3_e1++;
        e = sgrid->zm1e[e];
	c1 = sgrid->e2c[e][0];
      }
    }
  }
  i_free_1d(imape);
  for (n = 0; n < sgrid->nobc; n++) {
    open_bdrys_t *open = sgrid->open[n];
    open->no2_e1--;
    open->no3_e1--;
    open->to2_e1--;
    open->to3_e1--;
  }

  /*-----------------------------------------------------------------*/
  /* Get the edge orientation maps                                   */
  /* ceni[ee]=1 if e2c[0] is a (boundary) ghost cell (e2c[1] is wet) */
  /* ceni[ee]=0 if e2c[1] is a (boundary) ghost cell (e2c[0] is wet) */
  /* dir[ee] = 1 if the velocity vector is directed into the cell    */
  /* dir[ee] = -1 if the velocity vector is directed out of the cell */
  /* ini[ee] = index pointing into the domain                        */
  /* outi[ee] = index pointing out of the domain                     */
  /* inc[cc] = average index pointing into the domain                */
  /* bec[j][cc] = e for boundary edges, = 0 for interior edges       */
  /* bcc[j][cc] = c for normal boundary edges, = -c for tangential   */
  /*              boundary edges, = 0 for interior edges.            */
  /* olap[cc] = number of wet (non-boundary, non-ghost) cells        */
  /*            adjacent to the OBC cell at cc.                      */
  /* omape[ee][c] and nmape[ee][c] are arrays of pointers that map   */
  /*           to neighboring cell centres of c across the edge ee.  */
  for (n = 0; n < sgrid->nobc; n++) {
    open_bdrys_t *open = sgrid->open[n];
    for (ee = 1; ee <= open->no3_e1; ee++) {
      e = open->obc_e1[ee];
      c1 = sgrid->e2c[e][0];
      c2 = sgrid->e2c[e][1];
      if (sgrid->wgst[c1]) {
	open->ini[ee] = sgrid->e2e[e][0];
	open->outi[ee] = sgrid->e2e[e][1];
	open->ceni[ee] = 1;
	c = c2;
      }
      if (sgrid->wgst[c2]) {
	open->ini[ee] = sgrid->e2e[e][1];
	open->outi[ee] = sgrid->e2e[e][0];
	open->ceni[ee] = 0;
	c = c1;
      }
      cs = sgrid->m2d[c];
      for (j = 1; j <= sgrid->npe[cs]; j++)
	if (e == sgrid->c2e[j][c])
	  break;
      if (sgrid->eSc[j][cs] == 1)
	open->dir[ee] = -1;
      else
	open->dir[ee] = 1;

      for(cc = 1; cc <= open->no3_t; cc++) {
	if(open->obc_t[cc] == open->obc_e2[ee]) {
	  open->nepc[cc] += 1.0;
	}
      }
      j = open->ini[ee];
      c = open->obc_e2[ee];
      cs = sgrid->m2d[c];
      open->nmape[ee] = sgrid->c2c[j];
      j = open->outi[ee];
      open->omape[ee] = sgrid->c2c[j];
    }

    /* Find the direction that maps to the interior for this centre  */
    for (cc = 1; cc <= open->no3_t; cc++) {
      c = open->obc_t[cc];
      cs = sgrid->m2d[c];

      open->inc[cc] = get_bind(sgrid, cs, maskb);
      open->olap[cc] = 0;
      for (j = 1; j <= sgrid->npe[cs]; j++) {
	open->bec[j][cc] = 0;
	open->bcc[j][cc] = 0;
	cn = sgrid->c2c[j][c];

	if (sgrid->wgst[cn]) {
	  e = sgrid->c2e[j][c];
	  if (maske[e] >= 0) {
	    open->bec[j][cc] = e;
	    open->bcc[j][cc] = cn;

	  }
	} else if (maskb[sgrid->m2d[cn]]) {
	  open->bcc[j][cc] = -cn;
	} else {
	  open->olap[cc]++;
	}
      }
    }
    for(cc = 1; cc <= open->no3_t; cc++)
      open->nepc[cc] = 1.0 / open->nepc[cc];
  }

  /*-----------------------------------------------------------------*/
  /* Get the transfer vectors                                        */
  for (n = 0; n < sgrid->nobc; n++) {
    open_bdrys_t *open = sgrid->open[n];
    /* Allocate memory for master - slave transfer vectors */
    if (open->bcond_ele & (FILEIN | CUSTOM) || open->bcond_nor2d & FLATHR) {
      open->transfer_eta = d_alloc_1d(open->no2_t + 1);
    }
    i = (open->relax_zone_nor) ? open->relax_zone_nor : 1;
    if (open->bcond_nor & (FILEIN | CUSTOM))
      open->transfer_u1 = d_alloc_1d(i * open->no3_e1 + 1);
    i = (open->relax_zone_tan) ? open->relax_zone_tan : 1;
    if (open->bcond_tan & (FILEIN | CUSTOM))
      open->transfer_u2 = d_alloc_1d(i * open->to3_e1 + 1);
    if (open->bcond_nor2d & FILEIN)
      open->transfer_u1av = d_alloc_1d(open->no2_e1 + 1);
    if (open->bcond_nor2d & CUSTOM)
      open->transfer_u1av = d_alloc_1d(open->no3_e1 + 1);
    if (open->bcond_tan2d & (FILEIN | CUSTOM))
      open->transfer_u2av = d_alloc_1d(open->to3_e1 + 1);

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
      if (open->bcond_tra[tn] & NOTHIN) {
	if (open->bcond_tra[tn] & (TRCONC|TRCONF)) i = open->bgz;
    }

    }
    open->ttsz = 0;
    open->t_imap = NULL;
    open->t_transfer = NULL;
    if (open->ntt) {
      open->ttsz = open->no3_t + i * open->no3_e1 + 1;
      open->t_transfer = d_alloc_2d(open->ttsz, open->ntt);
      open->t_imap = i_alloc_2d(i + 1, open->ttsz);
    }
  }

  /* Get the cyclic boundary locations.                              */
  for (n = 0; n < sgrid->nobc; n++) {
    /* Tracers                                                       */

    for (cc = 1; cc <= sgrid->open[n]->no3_t; cc++) {
      c = sgrid->open[n]->obc_t[cc];
      cs = sgrid->m2d[c];
      if (sgrid->open[n]->bcond_ele & CYCLIC || 
	  ANY0(CYCLIC, sgrid->open[n]->bcond_tra, sgrid->open[n]->ntr)) {
	sgrid->open[n]->cyc_t[cc] = cyc_m2(sgrid, sgrid->c2c[imapc[cs]], sgrid->c2c[omapc[cs]], c);

	/* Move outwards an extra cell                               */
	if (params->us_type & US_HEX) 
	  sgrid->open[n]->cyc_t[cc] = sgrid->c2c[omapc[cs]][sgrid->open[n]->cyc_t[cc]];
      }
    }

    /* Normal velocity component                                              */
    for (ee = 1; ee <= sgrid->open[n]->no3_e1; ee++) {
      e = sgrid->open[n]->obc_e1[ee];
      c = sgrid->open[n]->obc_e2[ee];
      es = sgrid->m2de[e];
      cs = sgrid->m2d[c];
      j = (sgrid->open[n]->ceni[ee]) ? sgrid->e2e[e][1] : sgrid->e2e[e][0];
      if (sgrid->open[n]->bcond_nor & CYCLIC) {
	c1 = cyc_m2(sgrid, sgrid->c2c[imapc[cs]], sgrid->c2c[omapc[cs]], c);
	if (params->us_type & US_HEX)
	  c1 = sgrid->c2c[omapc[cs]][c1];
	if (maske[sgrid->c2e[j][c1]] >= 0) c1 = sgrid->c2c[omapc[cs]][c1];
	sgrid->open[n]->cyc_e1[ee] = sgrid->c2e[j][c1];
      }
    }

    /* Tangential velocity component                                          */
    for (ee = sgrid->open[n]->no3_e1 + 1; ee <= sgrid->open[n]->to3_e1; ee++) {
      e = sgrid->open[n]->obc_e1[ee];
      es = sgrid->m2de[e];
      c = sgrid->open[n]->obc_e2[ee];
      cs = sgrid->m2d[c];
      j = sgrid->e2e[e][0];
      if (sgrid->open[n]->bcond_tan & CYCLIC) {
	c1 = cyc_m2(sgrid, sgrid->c2c[imapc[cs]], sgrid->c2c[omapc[cs]], c);
	if (params->us_type & US_HEX)
	  c1 = sgrid->c2c[omapc[cs]][c1];
	if (maske[sgrid->c2e[j][c1]] >= 0) c1 = sgrid->c2c[omapc[cs]][c1];
	sgrid->open[n]->cyc_e1[ee] = sgrid->c2e[j][c1];
      }
    }

    /* Get the bottom coordinate vector                              */
    for (cc = 1; cc <= sgrid->open[n]->no2_t; cc++) {
      int cl;
      c = c2 = sgrid->open[n]->obc_t[cc];
      while(c != sgrid->zm1[c])
	c = sgrid->zm1[c];
      sgrid->open[n]->bot_t[cc] = sgrid->zp1[c];
    }
  }

  /* Set the interior face stagger for normal component */
  /* Inner staggerd cells adjacent to land are omitted */

  for (n = 0; n < sgrid->nobc; n++) {
    open_bdrys_t *open = sgrid->open[n];
    if (open->stagger & INFACE) {
      c1 = c2 = c3 = 1;
      /* u1 velocity */
      for (ee = 1; ee <= open->to3_e1; ee++) {
	b1 = open->obc_e1[ee];
	b2 = open->oi1_e1[ee];
	b3 = open->oi2_e1[ee];
	if (ee <= open->no3_e1) {
	  if (b2 != b3 && b3 != sgrid->em[b3]) {
	    open->obc_e1[c1] = b2;
	    open->oi1_e1[c1] = b3;
	    ei = e2e(sgrid, b3, sgrid->m2de[open->ini[ee]]);
	    open->oi2_e1[c1] = ei;
	    c1 += 1;
	    if (sgrid->zp1e[b1] == b1) c2 += 1;
	    c3 += 1;
	  }
	} else {
	  open->obc_e1[c3] = b1;
	  open->oi1_e1[c3] = b2;
	  open->oi2_e1[c3] = b3;
	  c3 += 1;
	} 
      }
      c1 = open->no3_e1 - c1 + 1;
      c2 = open->no2_e1 - c2 + 1;
      if (c1) {
	open->no3_e1 -= c1;
	open->no2_e1 -= c2;
	open->to2_e1 -= c1;
	open->to3_e1 -= c1;
      }
    }
  }

  /* Edge to centre boundary index map */
  for (n = 0; n < sgrid->nobc; n++) {
    open_bdrys_t *open = sgrid->open[n];
    for (ee = 1; ee <= open->no3_e1; ee++) {
      e = open->obc_e1[ee];
      es = sgrid->m2de[e];
      c = sgrid->e2c[e][open->ini[ee]];
      for (cc = 1; cc <= open->no3_t; cc++) {
	if (c == open->obc_t[cc]) {
	  open->e2c_e1[ee] = cc;
	  break;
	}
      }
    }
  }

  /* Set the flag                                                    */
  for (n = 0; n < sgrid->nobc; n++) {
    open_bdrys_t *open = sgrid->open[n];
    for (cc = 1; cc <= open->no3_t; cc++) {
      c = open->obc_t[cc];
      c1 = sgrid->c2cc[sgrid->m2d[c]];
      k = sgrid->s2k[c];
      flg[k][c1] |= ETASPEC;
      /*
      flg[k][c1] |= (U1BDRY|U2BDRY);
      if (open->ocodec) 
	flg[k][c1] |= (R_EDGE|F_EDGE);
      else
	flg[k][c1] |= (L_EDGE|B_EDGE);
      */
      if (sgrid->us_type & US_IJ)
	flag[k][sgrid->s2j[c]][sgrid->s2i[c]] = flg[k][c1];
    }
  }

  /* Get the ghost tracer transfer vector index mappings             */
  for (n = 0; n < sgrid->nobc; n++) {
    open_bdrys_t *open = sgrid->open[n];
    i = 1;
    if (open->ttsz) {
      for (cc = 1; cc <= open->no3_t; cc++) {
	open->t_imap[cc][0] = i;
	i++;
      }
      for (j = 1; j <= open->bgz; j++) {
	for (ee = 1; ee <= open->no3_e1; ee++) {
	  open->t_imap[ee][j] = i;
	  i++;
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Remove OBC ghosts from the ghost vector. These previously had   */
  /* to be included so that land cells tangential to the OBCs could  */
  /* be included in the ghost vector.                                */
  ii = sgrid->nbpt;
  jj = sgrid->nbptS;
  sgrid->nbpt = sgrid->nbptS = 1;
  for (cc = 1; cc <= ii; cc++) {
    int found = 0;
    c = sgrid->bpt[cc];
    for (n = 0; n < sgrid->nobc; n++) {
      open_bdrys_t *open = sgrid->open[n];
      for (i = 1; i <= open->no3_e1; i++) {
	c1 = open->ogc_t[i];
	if (c1 == c) found = 1;
      }
    }
    if (!found) {
      sgrid->bpt[sgrid->nbpt] = sgrid->bpt[cc];
      sgrid->dbin[sgrid->nbpt] = sgrid->dbin[cc];
      sgrid->dbpt[sgrid->nbpt] = sgrid->dbpt[cc];
      sgrid->bin[sgrid->nbpt++] = sgrid->bin[cc];
      if (cc <= jj) sgrid->nbptS++;
    }
  }
  sgrid->nbpt--;
  sgrid->nbptS--;
  if (DEBUG("init_m"))
    dlog("init_m", "\nOpen boundaries created OK\n");

  /* Compatibility */
  sgrid->wse = sgrid->w3_e1;
  sgrid->wsv = sgrid->w3_e2;

  /*sgrid->wsa = sgrid->w3_t;*/
  sgrid->wsa = i_alloc_1d(sgrid->szc);
  sgrid->a3_t = sgrid->b3_t;
  sgrid->a2_t = sgrid->b2_t;
  c1 = 1;
  for (cc = 1; cc <= sgrid->b2_t; cc++) {
    c = sgrid->w2_t[cc];
    sgrid->wsa[c1++] = c;
  }
  for (cc = 1; cc <= sgrid->v3_t; cc++) {
    c = sgrid->w3_t[cc];
    if (c != sgrid->m2d[c])
      sgrid->wsa[c1++] = c;
  }
  for (cc = sgrid->v3_t+1; cc <= sgrid->b3_t; cc++) {
    c = sgrid->w3_t[cc];
    if (c != sgrid->m2d[c])
      sgrid->wsa[c1++] = c;
  }

  /* Sizes */
  sgrid->szm = max(max(sgrid->szc, sgrid->sze), sgrid->szv);
  sgrid->szmS = max(max(sgrid->szcS, sgrid->szeS), sgrid->szvS);

  /* Set the sponge zones                                            */
  set_sponge_cells(sgrid);

  /*-----------------------------------------------------------------*/
  /* Allocate memory for the geometry vectors                        */
  alloc_geom_us(sgrid, (MAP_A | GRID_A | MASTER_A));

  /*-----------------------------------------------------------------*/
  /* Set the bathymetry value. This is overwritten from the input   */
  /* file for MANUAL operation.  */
  memset(sgrid->botz, 0, sgrid->szcS * sizeof(double));
  for (cc = 1; cc <= sgrid->b2_t; cc++) {
    c = sgrid->w2_t[cc];
    c1 = sgrid->c2cc[c];
    sgrid->botz[c] = bathy[c1];
  }
  if (params->sigma) {
    for (n = 0; n < params->nz; n++)
      params->layers[n] /= bmax;
  }

  /* Set a no-gradient condition over the sediment for gridz */
  sgrid->topgrid = layers[nz];
  for (cc = 1; cc <= sgrid->b2_t; cc++) {
    c = sgrid->bot_t[cc];
    sgrid->gridz[sgrid->zm1[c]] = sgrid->botz[sgrid->m2d[c]];
  }

  /* Find the deepest coordinte in the grid */
  bmax = 0.0;
  for (cc = 1; cc <= sgrid->b2_t; cc++) {
    c = sgrid->w2_t[cc];
    if (sgrid->botz[c] < bmax) {
      bmax = sgrid->botz[c];
      sgrid->cdeep = sgrid->bot_t[cc];
    }
  }

  /* Count the number of dry cells beneath the surface */
  sgrid->bdry = 0;
  for (cc = 1; cc <= ns2; cc++) {
    sgrid->bdry += kbot[cc];
  }

  /* Print cell location info if required                            */
  if (printlocx != NOTVALID && printlocy != NOTVALID) {
    double dist = 0, mdist = HUGE;
    double x1, y1;
    for (cc = 1; cc <= sgrid->b2_t; cc++) {
      c = sgrid->w2_t[cc];
      x1 = sgrid->cellx[c] - printlocx;
      y1 = sgrid->celly[c] - printlocy;
      dist = sqrt(x1 * x1 + y1 * y1);
      if (dist < mdist) {
	mdist = dist;
	c1 = c;
      }
    }
    printcell = c1;
  }

  if (printcell) {
    char buf[MAXSTRLEN];
    c = printcell;;
    cc = sgrid->c2cc[c];
    printf("cell cc=%d c=%d [%f %f]\n",cc,c,m->xloc[m->eloc[0][cc][0]],m->yloc[m->eloc[0][cc][0]]);
    for(n=1; n <= sgrid->npe[c]; n++) {
      v = sgrid->c2v[n][c];
      e = sgrid->c2e[n][c];
      c1 = sgrid->e2c[e][0];
      c2 = sgrid->e2c[e][1];
      strcpy(buf, "out");
      if (sgrid->eSc[n][c] == -1) strcpy(buf, "in");
      printf("c=%d cc=%d n=%d v=%d[%f %f] e=%d[%f %f]:%s co=%d[%f %f] c1=%d[%f %f] theta=%f\n",c,sgrid->c2cc[c],n,v,sgrid->gridx[v],sgrid->gridy[v],e,sgrid->u1x[e],sgrid->u1y[e], buf, c1, sgrid->cellx[c1],sgrid->celly[c1], c2, sgrid->cellx[c2],sgrid->celly[c2],sgrid->thetau1[e]*180.0/PI);
    }
  }
  if (DEBUG("init_m"))
    dlog("init_m", "\nMaster geometry created OK\n");

  /*-----------------------------------------------------------------*/
  /* Set up the windows                                              */
  TIMING_SET;
  window_build(sgrid, params);
  TIMING_DUMP(1, "  window_build");
  
  if (params->us_type & US_IJ)
    write_windows(sgrid, flag);

  /*-----------------------------------------------------------------*/
  /* Get the centre and edge angles                                  */
  write_us_map(sgrid, params);

  /*-----------------------------------------------------------------*/
  /* Make and initialise the master data structure                   */
  master = master_build(params, sgrid);
  geom = sgrid;
  master->geom = sgrid;
  master->sgrid = sgrid;
  get_filter(geom);
  if (DEBUG("init_m"))
    dlog("init_m", "\nMaster data created OK\n");

  /*-----------------------------------------------------------------*/
  /* Set up the dump data structure                                  */
  dumpdata = dumpdata_build(params, sgrid, master, flag, layers, topo);
  if (DEBUG("init_m"))
    dlog("init_m", "\nDumpdata created OK\n");

  /*-----------------------------------------------------------------*/
  /* Make a map from Delaunay triangulation to unstructured          */
  /* coordinate. Used for mapping (x,y) to i.                        */
  /*build_delaunay_cell(geom, params, tegl);*/
  create_delaunay_cell(geom, params);

  /*
  if (params->d) {
    delaunay *d;

    geom->d = d = params->d; 
    geom->tri2c = i_alloc_1d(d->ntriangles);
    memset(geom->tri2c, 0, d->ntriangles * sizeof(int));
    for (cc = 1; cc <= geom->b2_t; cc++) {
      double x, y;
      c = geom->w2_t[cc];
      x = geom->cellx[c];
      y = geom->celly[c];

      for (n = 0; n < d->ntriangles; n++) {
	triangle* t = &d->triangles[n];
	if ((x == d->points[t->vids[0]].x || x == d->points[t->vids[1]].x ||
	     x == d->points[t->vids[2]].x) && 
	    (y == d->points[t->vids[0]].y || y == d->points[t->vids[1]].y ||
	     y == d->points[t->vids[2]].y))
	  geom->tri2c[n] = c;
      }
    }
  }
  */

  /* Set geographical flag */
  geom->is_geog = (strlen(params->projection) > 0) &&
    (strcasecmp(params->projection, GEOGRAPHIC_TAG) == 0);

  /*-----------------------------------------------------------------*/
  /* Free memory                                                     */
  l_free_2d((long **)flg);
  i_free_1d(mask);
  i_free_1d(maskb);
  i_free_2d(neic);
  i_free_2d(neij);
  if (flag) l_free_3d((long ***)flag);
  if (topo) d_free_2d(topo);
}

/* END build_sparse_grid_us()                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Checks if two edges share the same vertex. If the first time this */
/* is encountered then increment the vertex vector (and set a map to */
/* the vertex cell if b1,c1... != NULL).                             */
/*-------------------------------------------------------------------*/
int check_vert2d(geometry_t *sgrid,
	       int e,          /* Edge #1                            */
	       int e1,         /* Edge #2                            */
	       int i1,         /* Index of e to check                */
	       int i2,         /* Index of e1 to check               */
	       double **locx,  /* x position                         */
	       double **locy,  /* y position                         */
	       int **emap,     /* Set to 1 for found vertices        */
	       int *maske,     /* Edge mask                          */
	       int *b1,
	       int *b2,
	       int *b3)
{
  int ret = 0, vs = 0, ef = 0;
  int pf = 0;

  if (locx[e][i1] == locx[e1][i2] && locy[e][i1] == locy[e1][i2]) {
    if (!emap[e][i1] && !emap[e1][i2]) {
      ef = 1;
      /* Wet vector - both vertices of the edge are wet              */
      if (maske[e] == -1 && maske[e1] == -1) {
	sgrid->v2_e2++;
	if (b1 != NULL) {
	  vs = *b1;
	  if (pf) printf("wet v=%d e1=%d e2=%d\n",*b1,e,e1);
	  sgrid->w2_e2[sgrid->v2_e2] = vs;
	  *b1 += 1;
	} else
	  vs = 1;
	ret = 1;
      }
      /* OBC vector - both vertices of the edge are OBCs             */
      if (maske[e] >= 0 && maske[e1] >= 0) {
	sgrid->b2_e2++;
	if (b2 != NULL) {
	  if (pf) printf("obc v=%d e1=%d e2=%d\n",*b2,e,e1);
	  sgrid->w2_e2[sgrid->b2_e2] = *b2;
	  vs = *b2;
	  *b2 += 1;
	} else
	  vs = 1;
	ret = 2;
      }
      /* Ghost vector - at least one vertex of the edge is a ghost   */
      if (maske[e] == -2 || maske[e1] == -2) {
	sgrid->n2_e2++;
	if (b3 != NULL) {
	  if (pf) printf("ghost v=%d e1=%d e2=%d\n",*b3,e,e1);
	  sgrid->w2_e2[sgrid->n2_e2] = *b3;
	  vs = *b3;
	  *b3 += 1;
	} else
	  vs = 1;
	ret = 3;
      }

      /*emap[e][i1] = (b1 == NULL) ? 1 : vs;*/
      emap[e][i1] = vs;
      if (vs == 0) {
	/* Ghost vector - OBC limits at wet-OBC intersections        */
	if ((maske[e] == -1 || maske[e1] >= 0) ||
	    (maske[e] >= 0 || maske[e1] == -1)) {
	  /* Place these in ghost vector. Sould be in the wet vector */
	  /* so that vorticity is computed at wet edges that share a */
	  /* vertex with a tangential OBC edge.
	  sgrid->n2_e2++;
	  if (b3 != NULL) {
	    if (pf) printf("ghost v=%d e1=%d e2=%d\n",*b3,e,e1);
	    sgrid->w2_e2[sgrid->n2_e2] = *b3;
	    vs = *b3;
	    *b3 += 1;
	  } else
	    vs = 1;
	  ret = 3;
	  */

	  /* Wet vector - OBC limits at wet-OBC intersections        */
	  sgrid->v2_e2++;
	  if (b1 != NULL) {
	    vs = *b1;
	    if (pf) printf("wet v=%d e1=%d e2=%d\n",*b1,e,e1);
	    sgrid->w2_e2[sgrid->v2_e2] = vs;
	    *b1 += 1;
	  } else
	    vs = 1;
	  ret = 1;

	}
	emap[e][i1] = vs;
      }
      if (vs == 0) {
	char et1[MAXSTRLEN], et2[MAXSTRLEN];
	if (maske[e] == -1)
	  strcpy(et1, "wet");
	else if (maske[e] == -2)
	  strcpy(et1, "ghost");
	else
	  sprintf(et1, "OBC%d\n", maske[e]);
	if (maske[e1] == -1)
	  strcpy(et2, "wet");
	else if (maske[e1] == -2)
	  strcpy(et2, "ghost");
	else
	  sprintf(et2, "OBC%d\n", maske[e1]);
	hd_quit("check_vert2d: Can't assign vertex at edges %d(%f %f):%s and %d(%f %f):%s. Check OBCs terminate in land cells. %d %d\n",
		e, 0.5*(locx[e][0]+locx[e][1]), 0.5*(locy[e][0]+locy[e][1]), et1,
		e1, 0.5*(locx[e1][0]+locx[e1][1]), 0.5*(locy[e1][0]+locy[e1][1]), et2,maske[e],maske[e1]);
      }
    }
    /* Increment the number of edges emanating from the vertex       */
    if (b1 != NULL && !emap[e1][i2]) {
      if (ef) 
	sgrid->nve[emap[e][i1]] += 2;
      else
	sgrid->nve[emap[e][i1]]++;
    }

    emap[e1][i2] = emap[e][i1];

    return(ret);
  }
  return(0);
}

/* END check_vert2d()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Checks if two edges share the same vertex. If the first time this */
/* is encountered then increment the vertex vector (and set a map to */
/* the vertex cell if b1,c1... != NULL).                             */
/*-------------------------------------------------------------------*/
int check_vert3d(geometry_t *sgrid,
		 int vs,         /* Surface vertex                   */
		 int vv,         /* Vertex counter                   */
		 int ed,         /* Deepest edge                     */
		 int eb,         /* Bottom cell for ei               */
		 int *c1, 
		 int *c2,
		 int *c3)
{
  int ret = 0, zp1;
  int k = sgrid->nz-1;
  int pf = 0;

  /* Count the surface vertex                                        */
  if (vv <= sgrid->v2_e2) {
    sgrid->v3_e2++;
    if (c1 != NULL) {
      sgrid->w3_e2[sgrid->v3_e2] = vs;
      sgrid->m2dv[vs] = vs;
    }
  }
  if (vv > sgrid->v2_e2 && vv <= sgrid->b2_e2) {
    sgrid->b3_e2++;
    if (c2 != NULL) {
      sgrid->w3_e2[sgrid->b3_e2] = vs;
      sgrid->m2dv[vs] = vs;
    }
  }
  if (vv > sgrid->b2_e2 && vv <= sgrid->n2_e2) {
    sgrid->n3_e2++;
    if (c3 != NULL) {
      sgrid->w3_e2[sgrid->n3_e2] = vs;
      sgrid->m2dv[vs] = vs;
    }
  }

  zp1 = vs;
  ed = sgrid->zm1e[ed];
  k--;

  while (ed <= eb) {
    /* Wet vector - both vertices of the edge are wet (or two    */
    /* wet edges intersect a tangential OBC edge - see above in  */
    /* check_vert2d().                                           */
    if (vv <= sgrid->v2_e2) {
      sgrid->v3_e2++;
      if (c1 != NULL) {
	sgrid->w3_e2[sgrid->v3_e2] = *c1;
	sgrid->zp1v[*c1] = zp1;
	sgrid->zm1v[zp1] = *c1;
	sgrid->m2dv[*c1] = vs;
	if(vs==pf)printf("wet edge=%d k=%d v=%d vv=%d\n",ed,k,*c1,sgrid->v3_e2);
	zp1 = *c1;
	*c1 += 1;
      }
    }
    /* OBC vector - both vertices of the edge are OBCs           */
    if (vv > sgrid->v2_e2 && vv <= sgrid->b2_e2) {
      sgrid->b3_e2++;
      if (c2 != NULL) {
	sgrid->w3_e2[sgrid->b3_e2] = *c2;
	sgrid->zp1v[*c2] = zp1;
	sgrid->zm1v[zp1] = *c2;
	sgrid->m2dv[*c2] = vs;
	if(vs==pf)printf("OBC edge=%d k=%d v=%d vv=%d\n",ed,k,*c2,sgrid->b3_e2);
	zp1 = *c2;
	*c2 += 1;
      }
    }
    /* Ghost vector - at least one vertex of the edge is a ghost */
    if (vv > sgrid->b2_e2 && vv <= sgrid->n2_e2) {
      sgrid->n3_e2++;
      if (c3 != NULL) {
	sgrid->w3_e2[sgrid->n3_e2] = *c3;
	sgrid->zp1v[*c3] = zp1;
	sgrid->zm1v[zp1] = *c3;
	sgrid->m2dv[*c3] = vs;
	if(vs==pf)printf("ghost edge=%d k=%d v=%d vv=%d\n",ed,k,*c3,sgrid->n3_e2);
	zp1 = *c3;
	*c3 += 1;
      }
    }
    ed = sgrid->zm1e[ed];
    k--;
  }
}

/* END check_vert3d()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Removes duplicate locations in a list. O(nlog(n)) operations.     */
/*-------------------------------------------------------------------*/
void find_vertices(geometry_t *sgrid, 
		   double **locx, 
		   double **locy,
		   int *maske,
		   int **emap
		   )
{
  int ee, e, e1, e2;
  int i1, i2, nn, n;
  int next, nedge;
  int b1, b2, b3;
  double **edge;
  int checkf = 0;

  /* Map the locations into a continuous vector                      */
  /* Get the vector size                                             */
  /* Allocate                                                        */
  nedge = 2 * sgrid->n2_e1;
  edge = d_alloc_2d(4, nedge);

  /* Make the vector                                                 */
  next = 0;
  for (ee = 1; ee <= sgrid->n2_e1; ee++) {
    e = sgrid->w2_e1[ee];
    edge[next][0] = locx[e][0];
    edge[next][1] = locy[e][0];
    edge[next][2] = (double)e;
    edge[next][3] = 0.0;
    next++;
  }
  for (ee = 1; ee <= sgrid->n2_e1; ee++) {
    e = sgrid->w2_e1[ee];
    edge[next][0] = locx[e][1];
    edge[next][1] = locy[e][1];
    edge[next][2] = (double)e;
    edge[next][3] = 1.0;
    next++;
  }

  /* Sort edges, duplicates are consecutive!                         */
  qsort(edge[0], nedge, sizeof(double)*4, edge_sort_double);

  /* Count the vertices                                              */
  for (n = 0; n < nedge - 1; n++) {
    nn = n + 1;
    if ((edge[n][0] == edge[nn][0]) &&
	(edge[n][1] == edge[nn][1])) {
      e = (int)edge[n][2];
      i1 = (int)edge[n][3];
      e1 = (int)edge[nn][2];
      i2 = (int)edge[nn][3];
      check_vert2d(sgrid, e, e1, i1, i2, locx, locy, emap, maske, NULL, NULL, NULL);
    }
  }

  /* Populate the 2D vertices. These are ordered as:                 */
  /* 1 to v2_e2 : vertices with both edges that are wet              */
  /* v2_e2 to b2_e2 : vertices with both edges that are OBCs         */
  /* b2_e2 to n2_e2 : vertices with one edge that is a ghost         */
  sgrid->b2_e2 += sgrid->v2_e2;
  sgrid->n2_e2 += sgrid->b2_e2;
  sgrid->szvS = sgrid->n2_e2 + 1;
  sgrid->w2_e2 = i_alloc_1d(sgrid->szvS);
  sgrid->nve = i_alloc_1d(sgrid->szvS);
  memset(sgrid->nve, 0, sgrid->szvS * sizeof(int));
  for (ee = 1; ee <= sgrid->n3_e1; ee++) {
    e = sgrid->w3_e1[ee];
    emap[e][0] = 0;
    emap[e][1] = 0;
  }
  b1 = 1;
  b2 = sgrid->v2_e2 + 1;
  b3 = sgrid->b2_e2 + 1;

  sgrid->n2_e2 = sgrid->b2_e2;
  sgrid->b2_e2 = sgrid->v2_e2;
  sgrid->v2_e2 = 0;

  /* Set the vertices to process                                     */
  for (n = 0; n < nedge - 1; n++) {
    nn = n + 1;
    if ((edge[n][0] == edge[nn][0]) &&
	(edge[n][1] == edge[nn][1])) {
      e = (int)edge[n][2];
      i1 = (int)edge[n][3];
      e1 = (int)edge[nn][2];
      i2 = (int)edge[nn][3];
      check_vert2d(sgrid, e, e1, i1, i2, locx, locy, emap, maske, &b1, &b2, &b3);
    }
  }
  d_free_2d(edge);
}

/* END find_vertices()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Removes duplicate locations in a list. O(nlog(n)) operations.     */
/*-------------------------------------------------------------------*/
int find_edges_tan(geometry_t *sgrid,
		   mesh_t *m, 
		   int k,
		   int *e, 
		   int *ee,
		   int *c4,
		   int *b4,
		   int *maske,
		   int **vle,
		   point *tegl
		   )
{
  int cc, c, cg, cco;
  int c1, g1, c2, g2;
  int j1, j2, nn, n, j;
  int next, nedge, npe;
  double **edge;
  int num_gst = 0;
  int pe = 0;

  /* Map the locations into a continuous vector                      */
  /* Get the vector size                                             */
  /* Allocate                                                        */
  nedge = 2 * sgrid->nbpt + 1;
  edge = d_alloc_2d(4, nedge);

  /* Make the vector                                                 */
  next = 0;
  for (cc = 1; cc <= sgrid->nbpt; cc++) {
    cg = sgrid->bpt[cc];
    c = sgrid->bin[cc];
    if (k != sgrid->s2k[c]) continue;
    j = sgrid->dbpt[cc];
    cco = sgrid->c2cc[c];
    if (!is_obc(cco, m->nobc, m->npts, m->loc)) {
      edge[next][0] = m->xloc[m->eloc[0][cco][j]];
      edge[next][1] = m->yloc[m->eloc[0][cco][j]];
      edge[next][2] = (double)c;
      edge[next][3] = (double)j;
      next++;
      edge[next][0] = m->xloc[m->eloc[1][cco][j]];
      edge[next][1] = m->yloc[m->eloc[1][cco][j]];
      edge[next][2] = (double)c;
      edge[next][3] = (double)j;
      next++;
    }
  }
  nedge = next;

  /* Sort edges, duplicates are consecutive!                         */
  qsort(edge[0], nedge, sizeof(double)*4, edge_sort_double);

  /* Count the tangential ghost edges                                */
  if (maske == NULL) {
    for (n = 0; n < nedge - 1; n++) {
      nn = n + 1;
      if ((edge[n][0] == edge[nn][0]) &&
	  (edge[n][1] == edge[nn][1])) {
	c = (int)edge[n][2];
	if (k == sgrid->nz - 1)
	  sgrid->x2_e1++;
	sgrid->x3_e1++;
	num_gst++;
      }
    }
    d_free_2d(edge);
    return (num_gst);
  } else {
    for (n = 0; n < nedge - 1; n++) {
      nn = n + 1;
      if ((edge[n][0] == edge[nn][0]) &&
	  (edge[n][1] == edge[nn][1])) {
	c1 = (int)edge[n][2];
	j1 = (int)edge[n][3];
	g1 = sgrid->c2c[j1][c1];
	c2 = (int)edge[nn][2];
	j2 = (int)edge[nn][3];
	g2 = sgrid->c2c[j2][c2];

	npe = sgrid->npe[sgrid->m2d[c1]];
	for (j = 1; j <= npe; j++)
	  if (c2 == sgrid->c2c[j][c1]) break;
	sgrid->e2c[*e][0] = g1;
	if (j <= npe) {
	  sgrid->c2e[j][g1] = *e;
	  sgrid->c2c[j][g1] = g2;
	  sgrid->e2e[*e][0] = j;
	}

	npe = sgrid->npe[sgrid->m2d[c2]];
	for (j = 1; j <= npe; j++)
	  if (c1 == sgrid->c2c[j][c2]) break;
	sgrid->e2c[*e][1] = g2;
	if (j <= npe) {
	  sgrid->c2e[j][g2] = *e;
	  sgrid->c2c[j][g2] = g1;
	  sgrid->e2e[*e][1] = j;
	}

	if (k == sgrid->nz - 1) {
	  tegl[*e].x = edge[n][0];
	  tegl[*e].y = edge[n][1];
	  sgrid->nbe1S++;
	  sgrid->w2_e1[*c4] = *e;
	  if (pe) printf("tangential ghost edge %d between (%d %d) interior (%d %d)\n",
			 *e, g1, g2, c1, c2);
	  *c4 += 1;
	}
	maske[*e] = -2;
	sgrid->bpte1[sgrid->nbe1] = *e;
	sgrid->bine1[sgrid->nbe1++] = c1;
	sgrid->w3_e1[*b4] = *e;

	if (*ee >= sgrid->szeS) hd_quit("Preprocessor error : check bathymetry specification.\n");
	vle[*ee][k] = *e;
	*b4 += 1;
	*e += 1;
	*ee += 1;
      }
    }
  }
  d_free_2d(edge);
}

/* END find_edges_tan()                                              */
/*-------------------------------------------------------------------*/

int is_obc(int c, int nobc, int *npts, int**loc)
{
  int n, cc;
  for (n = 0; n < nobc; n++) {
    for (cc = 1; cc <= npts[n]; cc++) {
      if (c == loc[n][cc])
	return(c);
    }
  }
  return(0);
}

int is_obce(int cc,            /* Cell centre                        */
	    int j,             /* j index of cell edge               */
	    mesh_t *m
	    )
{
  int n, nc;
  int b1 = m->eloc[0][cc][j];
  int b2 = m->eloc[1][cc][j];
  int o1, o2;
  double eps = 1e-5;

  for (n = 0; n < m->nobc; n++) {
    for (nc = 1; nc <= m->npts[n]; nc++) {
      o1 = m->obc[n][nc][0];
      o2 = m->obc[n][nc][1];
      if ((o1 == b1 && o2 == b2) || (o1 == b1 && o2 == b2)) return(n);
    }
  }
  return(-1);
}

int is_obceo(int npe,           /* Number of edges                    */
	    int cc,            /* Cell centre                        */
	    int j,             /* j index of cell edge               */
	    double **x,        /* Cell vertex x locations            */
	    double **y,        /* Cell vertex y locations            */
	    double ***posx,    /* x locations of OBC edges           */
	    double ***posy,    /* y locations of OBC edges           */
	    int nobc,          /* Number of OBCs                     */
	    int *npts          /* OBC list                           */
	    )
{
  int n, nc;
  double x1 = x[cc][j];
  double x2 = x[cc][jp(j, npe)];
  double y1 = y[cc][j];
  double y2 = y[cc][jp(j, npe)];
  double eps = 1e-5;

  for (n = 0; n < nobc; n++) {
    for (nc = 0; nc < npts[n]; nc++) {
      if (fabs(posx[n][nc][0] - x1) < eps && fabs(posy[n][nc][0] - y1) < eps && 
	  fabs(posx[n][nc][1] - x2) < eps && fabs(posy[n][nc][1] - y2) < eps)
	return(n);
      if (fabs(posx[n][nc][0] - x2) < eps && fabs(posy[n][nc][0] - y2) < eps && 
	  fabs(posx[n][nc][1] - x1) < eps && fabs(posy[n][nc][1] - y1) < eps)
	return(n);
    }
  }
  return(-1);
}

int obc_num(int c, int nobc, int *npts, int **loc)
{
  int n, cc;
  for (n = 0; n < nobc; n++) {
    for (cc = 1; cc <= npts[n]; cc++) {
      if (c == loc[n][cc])
	return(n);
    }
  }
  return(-1);
}

void obc_num_all(int c, int nobc, int *npts, int **loc, int *mask)
{
  int n, cc;

  if (nobc == 0) return;
  for (n = 0; n < nobc; n++) {
    mask[n] = -1;
    for (cc = 1; cc <= npts[n]; cc++) {
      if (c == loc[n][cc])
	mask[n] = n;
    }
  }
}

int zm1e(geometry_t *sgrid, int e)
{
  int ret;
  int c = sgrid->e2c[e][0];
  int j = sgrid->e2e[e][0];
  int cm = sgrid->zm1[c];
  if (sgrid->wgst[cm]) {
    c = sgrid->e2c[e][1];
    j = sgrid->e2e[e][1];
    cm = sgrid->zm1[c];
  }
  ret = sgrid->c2e[j][cm];
  ret = (ret) ? ret : e;
  return(ret);
}

int zm1es(geometry_t *sgrid, int e)
{
  int c1 = sgrid->e2c[e][0];
  int c2 = sgrid->e2c[e][1];
  int j = sgrid->e2e[e][0];
  int ret = sgrid->c2e[j][sgrid->zm1[c1]];
  ret = (ret) ? ret : e;
  c1 = sgrid->zm1[c1];
  c2 = sgrid->zm1[c2];
  if (c1 == sgrid->zm1[c1] || c2 == sgrid->zm1[c2]) ret = e;
  return(ret);
}

int zp1e(geometry_t *sgrid, int e)
{
  int c = sgrid->e2c[e][0];
  int j = sgrid->e2e[e][0];
  int ret = sgrid->c2e[j][sgrid->zp1[c]];
  ret = (ret) ? ret : e;
  return(ret);
}

int e2e(geometry_t *sgrid, int e, int j)
{
  int en, cn, c = sgrid->e2c[e][0];
  int n;
  int cs = sgrid->m2d[c];
  int npe = sgrid->npe[cs];

  for (n = 1; n <= npe; n++) {
    if (e == sgrid->c2e[n][c])
      break;
  }
  if (sgrid->wgst[c]) {
    c = sgrid->e2c[e][1];
    cs = sgrid->m2d[c];
    npe = sgrid->npe[cs];
    if (usejo)
      n = jo(n, npe);
    else
      n = jocw(sgrid, c, n);
  }
  cn = sgrid->c2c[j][c];
  cs = sgrid->m2d[cn];
  n = min(n, sgrid->npe[cs]);
  en = sgrid->c2e[n][cn];
  /*
  int i = (e == sgrid->c2e[j][c]) ? 1 : 0;
  ret = sgrid->c2e[j][sgrid->e2c[e][i]];
  ret = (ret) ? ret : e;
  */
  en = (en) ? en : e;
  return(en);
}

int e2eo(geometry_t *sgrid, int e, int j)
{
  int en, cn, c = sgrid->e2c[e][0];
  int n;
  int cs = sgrid->m2d[c];
  int npe = sgrid->npe[cs];
  printf("a1 c=%d npe=%d j=%d wgst=%d\n",c,npe,j,sgrid->wgst[c]);
  for (n = 1; n <= npe; n++) {
    if (e == sgrid->c2e[n][c])
      break;
  }
  if (sgrid->wgst[c]) {
    c = sgrid->e2c[e][1];
    cs = sgrid->m2d[c];
    npe = sgrid->npe[cs];
    n = jo(n, npe);
  }
  cn = sgrid->c2c[j][c];
  cs = sgrid->m2d[cn];
  n = min(n, sgrid->npe[cs]);
  en = sgrid->c2e[n][cn];
  /*
  int i = (e == sgrid->c2e[j][c]) ? 1 : 0;
  ret = sgrid->c2e[j][sgrid->e2c[e][i]];
  ret = (ret) ? ret : e;
  */
  en = (en) ? en : e;
  return(en);
}


void ete(geometry_t *sgrid, int e, int *e1, int *e2)
{
  int c = sgrid->e2c[e][0];
  int cs = sgrid->m2d[c];
  int npe = sgrid->npe[cs];
  int j = sgrid->e2e[e][0];
  j = (j <= npe/2) ? jm(j, npe) : jp(j, npe);
  *e1 = sgrid->c2e[j][c];
  c = sgrid->e2c[e][1];
  cs = sgrid->m2d[c];
  npe = sgrid->npe[cs];
  j = sgrid->e2e[e][1];
  j = (j <= npe) ? jm(j, npe) : jp(j, npe);
  *e2 = sgrid->c2e[j][c];
}

/*-------------------------------------------------------------------*/
/* Routine to read the batymetry and set SOLID and OUTSIDE cells     */
/*-------------------------------------------------------------------*/
void make_flags_us(parameters_t *params, /* Input parameters data structure */
		   unsigned long **flag,  /* Flags array */
		   double *bathy,  /* Bathymetry array */
		   double *layers, /* Layer spacing array */
		   int ns2,        /* Size of the grid */
		   int nz          /* Size of the grid in the z direction */
  )
{
  int i, j, k, c;               /* Counters */
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
      hd_quit("make_flags_us: BATHYMIN > BATHYMAX (%f > %f)\n", bmin, bmax);
    else
      bathylimit = 1;
  }

  /*-----------------------------------------------------------------*/
  /* Initialise */
  for (c = 1; c <= ns2; c++)
    for (k = 0; k < nz; k++) {
      flag[k][c] = 0;
    }

  for (c = 1; c <= ns2; c++) {
    /* Adjust botz to min and max bathymetry, checking for outside */
    /* cells (botz < gridz[0]) and land cells (botz > etamax).  */
    if (bathylimit && bathy[c] < etamax && bathy[c] > -bmin)
      bathy[c] = -bmin;
    if (bathylimit && bathy[c] >= layers[0] && bathy[c] < -bmax)
      bathy[c] = -bmax;

    /* Mark outside cells */
    if (bathy[c] < layers[0]) {
      flag[nz - 1][c] |= OUTSIDE;
    }

    /* Load the cell values */
    topflag = flag[nz - 1][c];
    for (k = 0; k < nz; k++) {
      if (topflag & OUTSIDE) {
	flag[k][c] |= OUTSIDE;
      } else {
	/* Don't allow cells to be less than MIN_CELL_THICKNESS */
	if (bathy[c] < layers[k + 1] && 
	    bathy[c] > layers[k + 1] - mincelldz[k]) { 
	  double newbotz;
	  if ((layers[k + 1] - bathy[c]) > 
	      (layers[k + 1] - layers[k]) / 2)                
	    newbotz = layers[k];
	  else
	    newbotz = layers[k + 1];
	  bathy[c] = newbotz;
	}
	if (bathy[c] >= layers[k + 1] || bathy[c] >= etamax) {
	  /* A solid cell */
	  flag[k][c] |= SOLID;
	}
      }
    }
  }
}

/* END make_flags_us()                                               */
/*-------------------------------------------------------------------*/



int jp(int j, int npe)
{
  int n = (j == npe) ? 1 : j+1;
  return(n);
}

int jm(int j, int npe)
{
  int n = (j == 1) ? npe : j-1;
  return(n);
}

int jo(int j, int npe)
{
  int n = j + npe /2;
  int jo = (n > npe) ? n - npe : n;
  return(jo);
}

/* Returns the opposite edge number to edge j                        */
int jocc(int **neic, int *npe, int c, int j)
{
  int jj;
  int cn = neic[j][c];

  for (jj = 1; jj <= npe[cn]; jj++) {
    if (c == neic[jj][cn]) {
      return(jj);
    }
  }
  return(jo(j, npe[c]));
}

int jocw(geometry_t *geom, int c, int j)
{
  int jj;
  int e = geom->c2e[j][c];

  jj = geom->e2e[e][0]; 
  if (j == geom->e2e[e][0])
    jj = geom->e2e[e][1];
  return(jj);    
}

/*-------------------------------------------------------------------*/
/* Returns the opposite edge to n                                    */
/*-------------------------------------------------------------------*/
int oedge(int npe, int n)
{
  int ret = (n > npe / 2) ? n - npe / 2 : n + npe / 2;
  return( ret);
}

/* END oedge()                                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to find a neighbour of edge j                             */
/*-------------------------------------------------------------------*/
int find_neighbour(int c, double **x, double **y, int *npe, int ns2, int *j)
{
  int jj, cc;
  int jn = jp(*j,npe[c]);

  for (cc = 1; cc <= ns2; cc++) {
    if (cc == c) continue;
    for (jj = 1; jj <= npe[cc]; jj++) {
      if (x[c][*j] == x[cc][jj] && x[c][jn] == x[cc][jm(jj,npe[cc])] && 
	  y[c][*j] == y[cc][jj] && y[c][jn] == y[cc][jm(jj,npe[cc])]) {
	*j = jm(jj, npe[cc]);
	return(cc);
      }
    }
  }
  return(0);
}

int find_neighbour_l(int c, int ***eloc, int *npe, int ns2, int *j)
{
  int jj, cc;

  for (cc = 1; cc <= ns2; cc++) {
    if (cc == c) continue;
    for (jj = 1; jj <= npe[cc]; jj++) {
      if ((eloc[0][c][*j] == eloc[0][cc][jj] && eloc[1][c][*j] == eloc[1][cc][jj]) ||
	  (eloc[0][c][*j] == eloc[1][cc][jj] && eloc[1][c][*j] == eloc[0][cc][jj])) {
	*j = jj;
	return(cc);
      }
    }
  }
  return(0);
}

/* END find_neighbour()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Re-orders the edge numbers. The index [0] contains the upstream   */
/* cell centre adjacent to cell e, and [1] the downstream centre.    */
/* The gradient at e is then e2c[e][0] - e2c[e][1].                  */
/*-------------------------------------------------------------------*/
void reorder_edges(geometry_t *sgrid) {
  int e, j, c, cs, cc;
  int c1 , c2, j1, j2, npe;

  for (cc = 1; cc <= sgrid->b3_t; cc++) {
    c = sgrid->w3_t[cc];
    cs = sgrid->m2d[c];
    npe = sgrid->npe[cs];
    for (j = 1; j <= npe; j++) {
      e = sgrid->c2e[j][c];
      c1 = sgrid->e2c[e][0];
      c2 = sgrid->e2c[e][1];
      j1 = sgrid->e2e[e][0];
      j2 = sgrid->e2e[e][1];
      if (npe == 4) {
	if ((j == 1 || j == 3) && j1 != 1)
          swap_edge(sgrid, e);
	if ((j == 2 || j == 4) && j1 != 4)
          swap_edge(sgrid, e);
      }
      if (npe == 5) {
	if (j == 1 && j1 != 1)
          swap_edge(sgrid, e);
	if (j == 2 && j1 != 2)
          swap_edge(sgrid, e);
	if (j == 3 && j1 == 3)
          swap_edge(sgrid, e);
	if (j == 4 && j1 == 2)
          swap_edge(sgrid, e);
	if (j == 5 && j1 != 5)
          swap_edge(sgrid, e);
      }
      if (npe == 6) {
	if ((j == 1 || j == 4) && j1 != 1)
          swap_edge(sgrid, e);
	if ((j == 2 || j == 5) && j1 != 2)
          swap_edge(sgrid, e);
	if ((j == 3 || j == 6) && j1 != 6)
          swap_edge(sgrid, e);
      }
    }
  }
}

/* END reorder_edges()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Swaps edge maps                                                   */
/*-------------------------------------------------------------------*/
void swap_edge(geometry_t *sgrid, int e)
{
  int cn = sgrid->e2c[e][0];
  int jn = sgrid->e2e[e][0];

  sgrid->e2c[e][0] = sgrid->e2c[e][1];
  sgrid->e2c[e][1] = cn;
  sgrid->e2e[e][0] = sgrid->e2e[e][1];
  sgrid->e2e[e][1] = jn;
}

/* END swap_edge()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Write the unstructured mapping summary                            */
/*-------------------------------------------------------------------*/
void write_us_map(geometry_t *sgrid,    /* Mapping data              */
		  parameters_t *params  /* Input parameters data     */
		  )
{
  FILE *fp = NULL;       /* File handles for printing data           */
  int n, m, j;           /* Counters                                 */
  int c, cc, e, ee, v, vv; /* Sparse location counters               */
  int cn, npe;
  long t;
  mesh_t *mesh = params->mesh;

  if (windows_log) {

    fp = fopen("window_map_us.txt", "w");
    if (fp == NULL)
      hd_quit("window_build: Can't write to file 'window_map_us.txt'\n");

    fprintf(fp, "\nInput file = %s\n\n", params->prmname);
    time(&t);
    fprintf(fp, "Written at :  %s\n\n", ctime(&t));

    fprintf(fp, "Size of the 3D sparse grid = %d\n", sgrid->sgnum);
    fprintf(fp, "Size of the 2D sparse grid = %d\n", sgrid->sgnumS);
    fprintf(fp, "Number of layers = %d\n", params->nz);
    if (mesh->nce1 > 0 && mesh->nce2 > 0)
      fprintf(fp, "Gridded mesh : size = %d x %d\n", mesh->nce1, mesh->nce2);
    fprintf(fp, "\n");

    fprintf(fp, "Number of 3D wet cells = %d\n", sgrid->v3_t);
    fprintf(fp, "Number of 2D wet cells = %d\n", sgrid->v2_t);

    fprintf(fp, "Number of 3D boundary cells = %d\n", sgrid->b3_t - sgrid->v3_t);
    fprintf(fp, "Number of 2D boundary cells = %d\n", sgrid->b2_t - sgrid->v2_t);

    m = sgrid->n3_t - sgrid->b3_t;
    fprintf(fp, "Number of 3D ghost cells = %d\n", m);
    /*
    if (m != sgrid->nbpt) fprintf(fp,"WARNING : 3D ghost vector sizes incompatible; %d != %d\n",
				   m, sgrid->nbpt);
    */
    n = sgrid->n2_t - sgrid->b2_t;
    fprintf(fp, "Number of 2D ghost cells = %d\n", n);
    /*
    if (n != sgrid->nbptS) fprintf(fp,"WARNING : 2D ghost vector sizes incompatible; %d != %d\n",
				   n, sgrid->nbptS);
    */
    fprintf(fp, "\n");
    fprintf(fp, "Number of 3D edge cells = %d\n", sgrid->n3_e1);
    fprintf(fp, "Number of 2D edge cells = %d\n", sgrid->n2_e1);
    fprintf(fp, "Number of 3D edge wet cells = %d\n", sgrid->v3_e1);
    fprintf(fp, "Number of 2D edge wet cells = %d\n", sgrid->v2_e1);
    fprintf(fp, "Number of 3D edge OBC cells = %d\n", sgrid->b3_e1 - sgrid->v3_e1);
    fprintf(fp, "Number of 2D edge OBC cells = %d\n", sgrid->b2_e1 - sgrid->v2_e1);
    fprintf(fp, "Number of 3D edge boundary cells = %d\n", sgrid->n3_e1 - sgrid->b3_e1);
    fprintf(fp, "Number of 2D edge boundary cells = %d\n", sgrid->n2_e1 - sgrid->b2_e1);

    fprintf(fp, "\n");
    fprintf(fp, "Number of 3D vertex cells = %d\n", sgrid->n3_e2);
    fprintf(fp, "Number of 2D vertex cells = %d\n", sgrid->n2_e2);
    fprintf(fp, "Number of 3D vertex cells with one wet edge = %d\n", sgrid->v3_e2);
    fprintf(fp, "Number of 2D vertex cells with one wet edge = %d\n", sgrid->v2_e2);
    fprintf(fp, "Number of 3D vertex cells with OBC edges = %d\n", sgrid->b3_e2 - sgrid->v3_e2);
    fprintf(fp, "Number of 2D vertex cells with OBC edges = %d\n", sgrid->b2_e2 - sgrid->v2_e2);
    fprintf(fp, "Number of 3D vertex cells with ghost edges = %d\n", sgrid->n3_e2 - sgrid->b3_e2);
    fprintf(fp, "Number of 2D vertex cells with ghost edges = %d\n", sgrid->n2_e2 - sgrid->b2_e2);
    
    fprintf(fp, "\nUnstructured index (c) to input file list (cc) map\n");
    fprintf(fp, "c cc\n");
    for (cc = 1; cc <= sgrid->ewetS; cc++)
      fprintf(fp, "%d  %d\n",cc, sgrid->c2cc[cc]);
    fprintf(fp, "\nCells to process : centres\n");
    fprintf(fp, "    cc    c    ");
    for (n = 1; n <= sgrid->npem; n++)
      fprintf(fp, "e%1d   ", n);
    fprintf(fp, "status\n");
    for (cc = 1; cc <= sgrid->n2_t; cc++) {
      c = sgrid->w2_t[cc];
      npe = sgrid->npe[c];
      if (cc >= 1 && cc <= sgrid->v2_t) {
	fprintf(fp, "%5d %5d ",cc, sgrid->w2_t[cc]);
	for (n = 1; n <= npe; n++)
	  fprintf(fp, "%5d", sgrid->c2e[n][c]);
	fprintf(fp, " : wet/land boundary\n");
      }
      if (cc > sgrid->v2_t && cc <= sgrid->b2_t) {
	fprintf(fp, "%5d %5d ",cc, sgrid->w2_t[cc]);
	for (n = 1; n <= npe; n++)
	  fprintf(fp, "%5d", sgrid->c2e[n][c]);
	fprintf(fp, " : OBC\n");
      }

      if (cc > sgrid->b2_t && cc <= sgrid->n2_t) {
	fprintf(fp, "%5d %5d ",cc, sgrid->w2_t[cc]);
	for (n = 1; n <= npe; n++)
	  fprintf(fp, "%5d", sgrid->c2e[n][c]);
	fprintf(fp, " : ghost\n");
      }
    }

    fprintf(fp, "\nCells to process : edges\n");
    fprintf(fp, "    ee    e    c1    c2    v1    v2   status\n");
    for (ee = 1; ee <= sgrid->n2_e1; ee++) {
      e = sgrid->w2_e1[ee];
      if (ee >= 1 && ee <= sgrid->v2_e1)
	fprintf(fp, "%5d %5d %5d %5d %5d %5d : wet\n",ee, sgrid->w2_e1[ee], 
		sgrid->e2c[e][0],sgrid->e2c[e][1], sgrid->e2v[e][0],sgrid->e2v[e][1]);
      if (ee > sgrid->v2_e1 && ee <= sgrid->b2_e1) {
	if (!sgrid->wgst[sgrid->e2c[e][0]] && !sgrid->wgst[sgrid->e2c[e][1]])
	  fprintf(fp, "%5d %5d %5d %5d %5d %5d : OBC - tan\n",ee, sgrid->w2_e1[ee],
		  sgrid->e2c[e][0],sgrid->e2c[e][1], sgrid->e2v[e][0],sgrid->e2v[e][1]);
	else
	  fprintf(fp, "%5d %5d %5d %5d %5d %5d : OBC - nor\n",ee, sgrid->w2_e1[ee],
		  sgrid->e2c[e][0],sgrid->e2c[e][1], sgrid->e2v[e][0],sgrid->e2v[e][1]);
      }
      if (ee > sgrid->b2_e1 && ee <= sgrid->n2_e1)
	fprintf(fp, "%5d %5d %5d %5d %5d %5d : land boundary\n",ee, sgrid->w2_e1[ee],
		sgrid->e2c[e][0],sgrid->e2c[e][1], sgrid->e2v[e][0],sgrid->e2v[e][1]);
    }

    fprintf(fp, "\nCells to process : vertices\n");
    fprintf(fp, "    vv    v   e1   e2   e3   e4   c1    c2    c3    status\n");
    for (vv = 1; vv <= sgrid->n2_e2; vv++) {
      v = sgrid->w2_e2[vv];
      if (vv >= 1 && vv <= sgrid->v2_e2) {
	fprintf(fp, "%5d %5d",vv, v);
	for (n = 1; n <= sgrid->nve[v]; n++)
	  fprintf(fp, "%5d", sgrid->v2e[v][n]);
	for (n = 1; n <= sgrid->nvc[v]; n++)
	  fprintf(fp, "%5d", sgrid->v2c[v][n]);
	fprintf(fp, " : wet\n");
      }
      if (vv > sgrid->v2_e2 && vv <= sgrid->b2_e2) {
	fprintf(fp, "%5d %5d",vv, v);
	for (n = 1; n <= sgrid->nve[v]; n++)
	  fprintf(fp, "%5d", sgrid->v2e[v][n]);
	for (n = 1; n <= sgrid->nvc[v]; n++)
	  fprintf(fp, "%5d", sgrid->v2c[v][n]);
	fprintf(fp, " : wet/OBC\n");
      }
      if (vv > sgrid->b2_e2 && vv <= sgrid->n2_e2) {
	fprintf(fp, "%5d %5d",vv, v);
	for (n = 1; n <= sgrid->nve[v]; n++)
	  fprintf(fp, "%5d", sgrid->v2e[v][n]);
	for (n = 1; n <= sgrid->nvc[v]; n++)
	  fprintf(fp, "%5d", sgrid->v2c[v][n]);
	fprintf(fp, " : wet/land\n");
      }
    }

    fprintf(fp, "\nGhost cells : centres\n");
    fprintf(fp, "    cc    c    ci\n");
    for (cc = 1; cc <= sgrid->nbptS; cc++) {
      c = sgrid->bpt[cc];
      cn = sgrid->bin[cc];
      fprintf(fp, "%5d %5d %5d\n",cc, c, cn);
    }

    fprintf(fp, "\nGhost cells : edges\n");
    fprintf(fp, "    cc    c    ci\n");
    for (ee = 1; ee <= sgrid->nbpte1S; ee++) {
      e = sgrid->bpte1[ee];
      c = sgrid->bine1[ee];
      fprintf(fp, "%5d %5d %5d\n",ee, e, c);
    }

    fprintf(fp, "\nSurface/bottom cells : centres\n");
    fprintf(fp, "    cc    cs   cb\n");
    for (cc = 1; cc <= sgrid->b2_t; cc++) {
      fprintf(fp, "%5d %5d %5d\n",cc, sgrid->sur_t[cc], sgrid->bot_t[cc]);
    }

    fprintf(fp, "\nSurface/bottom cells : edges\n");
    fprintf(fp, "    ee    es   eb\n");
    for (ee = 1; ee <= sgrid->n2_e1; ee++) {
      fprintf(fp, "%5d %5d %5d\n",ee, sgrid->sur_e1[ee], sgrid->bot_e1[ee]);
    }

    fprintf(fp, "\nOpen boundary edges\n");
    for (n = 0; n < sgrid->nobc; n++) {
      open_bdrys_t *open = sgrid->open[n];
      fprintf(fp, "   obc   ee    e    status\n");
      for (ee = 1; ee <= open->no2_e1; ee++)
	fprintf(fp, "%5d %5d %5d     nor\n",n, ee, open->obc_e1[ee]);
      for (ee = open->no3_e1 + 1; ee <= open->to2_e1; ee++)
	fprintf(fp, "%5d %5d %5d     tan\n",n, ee, open->obc_e1[ee]);
    }

    fprintf(fp, "\nCell dimensions\n");
    fprintf(fp, "    c     j     e   cellx    celly     u1x       u1y       h1au1     h2au1    area\n");
    for (cc = 1; cc <= sgrid->b2_t; cc++) {
      c = sgrid->w2_t[cc];
      npe = sgrid->npe[c];
      for (j = 1; j <= npe; j++) {
	e = sgrid->c2e[j][c];
	fprintf(fp, "%5d %5d %5d %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f\n",c, j, e, 
		sgrid->cellx[c],sgrid->celly[c],sgrid->u1x[e],sgrid->u1y[e],
		sgrid->h1au1[e], sgrid->h2au1[e], sgrid->cellarea[c]);
      }
    }
    fprintf(fp, "\n    v     gridx  gridy  area\n");
    for (cc = 1; cc <= sgrid->b2_e2; cc++) {
      v = sgrid->w2_e2[cc];
      fprintf(fp, "%5d %9.2f %9.2f  %9.2f\n",v, sgrid->gridx[v],sgrid->gridy[v],sgrid->dualarea[v]);
    }
    fprintf(fp, "\n    c     j     e h1au1    h2au1\n");
    for (cc = 1; cc <= sgrid->b2_t; cc++) {
      c = sgrid->w2_t[cc];
      npe = sgrid->npe[c];
      for (j = 1; j <= npe; j++) {
	e = sgrid->c2e[j][c];
	fprintf(fp, "%5d %5d %5d %5.2f %5.2f\n",c, j, e, sgrid->h1au1[e], sgrid->h2au1[e]);
      }
    }
    fclose(fp);
  }
}

/* END write_us_map()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to map into the interior until an open boundary is found  */
/* (the map becomes self mapping) then set up the interior sparse    */
/* location for tracer and tangential cyclic open boundaries.        */
/*-------------------------------------------------------------------*/
int cyc_m2(geometry_t *sgrid, /* Window geometry                   */
	   int *nmap,         /* Map in direction normal to boundary */
	   int *omap,         /* Map in direction normal to boundary */
	   int c             /* Self mapping sparse location */
  )
{
  int cyc;                      /* Cyclic sparse coordinate */
  cyc = c;

  while (c != nmap[nmap[c]]) {
    c = nmap[c];
    printf("%d %d\n",cyc,c);
  }
  cyc = omap[omap[c]];
  return (cyc);
}

/* END cyc_m2()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to allocate memory from the geometry data structure       */
/* arrays.                                                           */
/*-------------------------------------------------------------------*/
void alloc_geom_us(geometry_t *geom, /* Model geometry structure */
		   unsigned long mask  /* Mask for which arrays to allocate */
  )
{
  /* Mapping arrays */
  if (mask & MAP_A) {
    geom->sask = i_alloc_1d(geom->szm);
  }
  if (mask & GRID_A) {
    /* 3D arrays */
    geom->gridz = d_alloc_1d(geom->szc);
    geom->cellz = d_alloc_1d(geom->szc);
    /* 2D arrays */
    geom->botz = d_alloc_1d(geom->szcS);
    geom->botzgrid = d_alloc_1d(geom->szvS);
    geom->dHde1 = d_alloc_1d(geom->szeS);
    geom->botzu1 = d_alloc_1d(geom->szeS);

    if (geom->sednz) {
      geom->gridz_sed = d_alloc_2d(geom->szcS, geom->sednz + 1);
      geom->cellz_sed = d_alloc_2d(geom->szcS, geom->sednz);
    }
  }
  if (mask & MASTER_A) {
    geom->sinthcell = d_alloc_1d(geom->szcS);
    geom->costhcell = d_alloc_1d(geom->szcS);
    geom->sinthu1 = d_alloc_1d(geom->szeS);
    geom->costhu1 = d_alloc_1d(geom->szeS);
    geom->sinthu2 = d_alloc_1d(geom->szeS);
    geom->costhu2 = d_alloc_1d(geom->szeS);
  }
  if (mask & WINDOW_A) {
    geom->wsa = i_alloc_1d(geom->szc);
    geom->wse = i_alloc_1d(geom->sze);
    geom->wsv = i_alloc_1d(geom->szv);
  }
  if (mask & CENTRE_A) {
    geom->gridz = d_alloc_1d(geom->szc);
    geom->cellz = d_alloc_1d(geom->szc);
    if (geom->sednz) {
      geom->gridz_sed = d_alloc_2d(geom->szcS, geom->sednz + 1);
      geom->cellz_sed = d_alloc_2d(geom->szcS, geom->sednz);
    }
    geom->sinthcell = d_alloc_1d(geom->szcS);
    geom->costhcell = d_alloc_1d(geom->szcS);
    geom->botz = d_alloc_1d(geom->szcS);
    geom->wsa = i_alloc_1d(geom->szc);
    geom->c2c = i_alloc_2d(geom->szc, geom->npem+1);
    geom->zp1 = i_alloc_1d(geom->szc);
    geom->zm1 = i_alloc_1d(geom->szc);
    geom->cellarea = d_alloc_1d(geom->szcS);
    geom->cellx = d_alloc_1d(geom->szcS);
    geom->celly = d_alloc_1d(geom->szcS);
    geom->hacell = d_alloc_2d(geom->szcS, geom->npem + 1);
    geom->npe = i_alloc_1d(geom->szcS);
  }
  if (mask & EDGE_A) {
    geom->wse = i_alloc_1d(geom->sze);
    geom->sinthu1 = d_alloc_1d(geom->szeS);
    geom->costhu1 = d_alloc_1d(geom->szeS);
    geom->sinthu2 = d_alloc_1d(geom->szeS);
    geom->costhu2 = d_alloc_1d(geom->szeS);
    geom->dHde1 = d_alloc_1d(geom->szeS);
    geom->botzu1 = d_alloc_1d(geom->szeS);
    geom->h2au1 = d_alloc_1d(geom->szeS);
    geom->h1au1 = d_alloc_1d(geom->szeS);
    geom->h1acell = d_alloc_1d(geom->szeS);
    geom->edgearea = d_alloc_1d(geom->szeS);
    geom->thetau1 = d_alloc_1d(geom->szeS);
    geom->thetau2 = d_alloc_1d(geom->szeS);
    geom->u1x = d_alloc_1d(geom->szeS);
    geom->u1y = d_alloc_1d(geom->szeS);
    geom->nee = i_alloc_1d(geom->szeS);
  }
  if (mask & VERTEX_A) {
    geom->wsv = i_alloc_1d(geom->szv);
    geom->gridx = d_alloc_1d(geom->szvS);
    geom->gridy = d_alloc_1d(geom->szvS);
    geom->botzgrid = d_alloc_1d(geom->szvS);
    geom->dualarea = d_alloc_1d(geom->szvS);
    geom->dualareap = d_alloc_2d(geom->nvcm + 1, geom->szvS);
    geom->nve = i_alloc_1d(geom->szvS);
    geom->nvc = i_alloc_1d(geom->szvS);
  }
}

/* END alloc_geom_us()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to point the common arrays in a window geometry structure */
/* to the corresponding array in the master geometery.               */
/*-------------------------------------------------------------------*/
void point_geom_us(geometry_t *window, /* Window geometry structure */
		   geometry_t *geom  /* Model geometry structure */
  )
{
  if (geom->enon != window->enon)
    hd_quit
      ("point_geom: Master and window 3D sparse array sizes incompatible\n");

  if (geom->enonS != window->enonS)
    hd_quit
      ("point_geom: Master and window 2D sparse array sizes incompatible\n");

  /* Sizes */
  window->szc = geom->szc;
  window->sze = geom->sze;
  window->szv = geom->szv;
  window->szcS = geom->szcS;
  window->szeS = geom->szeS;
  window->szvS = geom->szvS;
  window->szm = geom->szm;
  window->szmS = geom->szmS;

  /* Point the window maps at the master */
  window->npe = geom->npe;
  window->nee = geom->nee;
  window->nve = geom->nve;
  window->nvc = geom->nvc;

  window->c2c = geom->c2c;
  window->c2e = geom->c2e;
  window->c2v = geom->c2v;
  window->e2e = geom->e2e;
  window->e2c = geom->e2c;
  window->e2v = geom->e2v;
  window->v2c = geom->v2c;
  window->v2e = geom->v2e;

  window->eSc = geom->eSc;
  window->eSv = geom->eSv;
  window->eSe = geom->eSe;
  window->vIc = geom->vIc;
  window->wAe = geom->wAe;
  window->wSe = geom->wSe;
  window->dualarea = geom->dualarea;
  window->dualareap = geom->dualareap;

  window->ep = geom->ep;
  window->em = geom->em;
  window->zm1 = geom->zm1;
  window->zp1 = geom->zp1;
  window->zm1e = geom->zm1e;
  window->zp1e = geom->zp1e;
  window->e2k = geom->e2k;
  window->zm1v = geom->zm1v;
  window->zp1v = geom->zp1v;

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
  window->hacell = geom->hacell;
  window->h1au1 = geom->h1au1;
  window->h2au1 = geom->h2au1;
  window->cellarea = geom->cellarea;
  window->edgearea = geom->edgearea;
  window->dHde1 = geom->dHde1;
  window->botz = geom->botz;
  window->botzu1 = geom->botzu1;
  window->botzgrid = geom->botzgrid;
  window->cellx = geom->cellx;
  window->celly = geom->celly;
  window->gridx = geom->gridx;
  window->gridy = geom->gridy;
  window->u1x = geom->u1x;
  window->u1y = geom->u1y;
  window->sinthcell = geom->sinthcell;
  window->costhcell = geom->costhcell;
  window->sinthu1 = geom->sinthu1;
  window->costhu1 = geom->costhu1;
  window->sinthu2 = geom->sinthu2;
  window->costhu2 = geom->costhu2;
  window->sask = geom->sask;
}

/* END point_geom_us()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Prints information associated with edge ed                        */
/*-------------------------------------------------------------------*/
void print_e(geometry_t *sgrid, mesh_t * m, int ed)
{
  int ee, e, e1, k;
  int c1, c2, i, j, n;
  int found = 0;

  for (ee = 1; ee <= sgrid->n3_e1; ee++) {
    e = sgrid->w3_e1[ee];
    k = -1;
    c1 = sgrid->e2c[e][0];
    c2 = sgrid->m2d[c1];
    i = m->iloc[sgrid->c2cc[c2]];
    j = m->jloc[sgrid->c2cc[c2]];
    if (sgrid->s2k) k = sgrid->s2k[c1];
    if (e == ed) {
      printf("e2c0 %d cs=%d (cc=%d)(%d %d %d) e2e=%d\n",c1, c2, sgrid->c2cc[c2], i, j, k, sgrid->e2e[e][0]);
      found = 1;
    }
    c1 = sgrid->e2c[e][1];
    c2 = sgrid->m2d[c1];
    i = m->iloc[sgrid->c2cc[c2]];
    j = m->jloc[sgrid->c2cc[c2]];
    if (sgrid->s2k) k = sgrid->s2k[c1];
    if (e == ed)
      printf("e2c0 %d cs=%d (cc=%d)(%d %d %d) e2e=%d\n",c1, c2, sgrid->c2cc[c2], i, j, k, sgrid->e2e[e][1]);
    if (e == ed) {
      e1 = e;
      n = sgrid->nz-1;
      while (e1 != sgrid->zm1e[e1]) {
	printf("k=%d e=%d\n",n, e1);
	e1 = sgrid->zm1e[e1];
	n--;
      }
    }
  }
  if (!found) printf("Can't find edge %d\n", ed);
}

/* END print_e()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Set up a delaunay triangulation for xytoi interpolation using     */
/* delaunay_xytoi(). The triangulation is based on grid centres and  */
/* grid vertices.                                                    */
/*-------------------------------------------------------------------*/
void create_delaunay_cell(geometry_t *geom, parameters_t *params)
{
  FILE *vf;
  delaunay *d = NULL;
  char key[MAXSTRLEN], buf[MAXSTRLEN];
  int n, cc, c, vv, v;
  int np;
  point *pin;
  int filef = (params->meshinfo) ? 1 : 0;
  int ee, e, c1, c2;
  int ntriangles;
  int inc_gst = 0;    /* Include ghost cell coordinates              */

  np = geom->b2_t + geom->n2_e2;
  if (inc_gst) np += (geom->nbpt + (geom->x2_e1 - geom->n2_e1));
  pin = malloc(np * sizeof(point));

  for (cc = 1, n = 0; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    pin[n].x = geom->cellx[c];
    pin[n++].y = geom->celly[c];
  }
  for (vv = 1; vv <= geom->n2_e2; vv++) {
    v = geom->w2_e2[vv];
    pin[n].x = geom->gridx[v];
    pin[n++].y = geom->gridy[v];
  }

  if (inc_gst) {
    /* Set the location of the ghost cell centres. This is also      */
    /* performed in master_setghosts().                              */
    for (cc = 1; cc <= geom->nbpte1S; cc++) {
      double x, y, d, r;
      int j;
      e = geom->bpte1S[cc];
      c = geom->bine1S[cc];
      x = geom->u1x[e] - geom->cellx[c];
      y = geom->u1y[e] - geom->celly[c];
      d = sqrt(x * x + y * y);
      r = atan2(y, x);
      for (j = 1; j <= geom->npe[c]; j++) {
	if (e == geom->c2e[j][c]) {
	  c1 = geom->c2c[j][c];
	  geom->cellx[c1] = geom->u1x[e] + d * cos(r);
	  geom->celly[c1] = geom->u1y[e] + d * sin(r);
	  pin[n].x = geom->cellx[c1];
	  pin[n++].y = geom->celly[c1];
	}
      }    
    }
    
    /* Set the location of the tangential ghost edges. This is the   */
    /* mid point joining the two cell centres adjacent to the edge.  */
    for (ee = geom->n2_e1+1; ee <= geom->x2_e1; ee++) {
      e = geom->w2_e1[ee];
      c1 = geom->e2c[e][0];
      c2 = geom->e2c[e][1];
      
      if (geom->wgst[c1] != geom->wgst[c2]) {
	geom->u1x[e] = 0.5 * (geom->cellx[c1] + geom->cellx[c2]);
	geom->u1y[e] = 0.5 * (geom->celly[c1] + geom->celly[c2]);
	pin[n].x = geom->u1x[e];
	pin[n++].y = geom->u1y[e];
      }
    }
    np = n;
  }

  if (d) free((delaunay *)d);
  /*d = delaunay_voronoi_build(np, pin, 0, NULL, 0, NULL, 0.0);*/
  geom->d = delaunay_build(np, pin, 0, NULL, 0, NULL);

  d = geom->d;
  geom->tri2c = i_alloc_1d(d->ntriangles);
  geom->c2tri = i_alloc_1d(geom->szcS);
  memset(geom->tri2c, 0, d->ntriangles * sizeof(int));
  TIMING_SET;
#if defined(HAVE_OMP)
#pragma omp parallel for private(c,n)
#endif
  for (cc = 1; cc <= geom->b2_t; cc++) {
    double x, y;
    c = geom->w2_t[cc];
    x = geom->cellx[c];
    y = geom->celly[c];

    for (n = 0; n < d->ntriangles; n++) {
      triangle* t = &d->triangles[n];
      if ((x == d->points[t->vids[0]].x || x == d->points[t->vids[1]].x ||
	   x == d->points[t->vids[2]].x) && 
	  (y == d->points[t->vids[0]].y || y == d->points[t->vids[1]].y ||
	   y == d->points[t->vids[2]].y)) {
	geom->tri2c[n] = c;
	geom->c2tri[c] = n;
      }
    }
  }
  TIMING_DUMP(1, "  create_deluanay nested loop");

  if (inc_gst) {
    for (cc = 1; cc <= geom->nbptS; cc++) {
      double x, y;
      c = geom->bpt[cc];
      x = geom->cellx[c];
      y = geom->celly[c];

      for (n = 0; n < d->ntriangles; n++) {
	triangle* t = &d->triangles[n];
	if ((x == d->points[t->vids[0]].x || x == d->points[t->vids[1]].x ||
	     x == d->points[t->vids[2]].x) && 
	    (y == d->points[t->vids[0]].y || y == d->points[t->vids[1]].y ||
	     y == d->points[t->vids[2]].y)) {
	  geom->tri2c[n] = c;
	  geom->c2tri[c] = n;
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Write to file                                                   */
  if (filef) {
    mesh_ofile_name(params, key);
    sprintf(buf,"%s_v.txt", key);
    if ((vf = fopen(buf, "w")) == NULL) filef = 0;
  }
  if (filef) {
    for (n = 0; n < d->ntriangles; n++) {
      triangle* t = &d->triangles[n];
      fprintf(vf, "%f %f\n", d->points[t->vids[0]].x, d->points[t->vids[0]].y);
      fprintf(vf, "%f %f\n", d->points[t->vids[1]].x, d->points[t->vids[1]].y);
      fprintf(vf, "%f %f\n", d->points[t->vids[2]].x, d->points[t->vids[2]].y);
      fprintf(vf, "%f %f\n", d->points[t->vids[0]].x, d->points[t->vids[0]].y);
      fprintf(vf, "NaN NaN\n");
    }
    fclose(vf);
  }
}

/* END create_delaunay_cell()                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Set up a delaunay triangulation for xytoi interpolation using     */
/* an explicitly defined triangulation based on grid centres and     */
/* grid vertices.                                                    */
/*-------------------------------------------------------------------*/
void build_delaunay_cell(geometry_t *geom, parameters_t *params, point *tegl)
{
  FILE *vf;
  delaunay* d;
  char key[MAXSTRLEN], buf[MAXSTRLEN];
  int n, cc, c, cg, vv, v;
  int **nei;
  point *pin;
  int filef = (params->meshinfo) ? 1 : 0;
  int j, jp;
  int ee, e, c1, c2;
  int inc_gst = 1;    /* Include ghost cell coordinates              */
  int dotest = 0;
  int *c2n, *v2n, *e2n, *eg2n;

  geom->d = delaunay_create();
  d = geom->d;
  d->npoints = geom->b2_t + geom->n2_e2;
  if (inc_gst) d->npoints += (geom->nbpt + (geom->x2_e1 - geom->n2_e1));
  d->points = malloc(d->npoints * sizeof(point));
  c2n = i_alloc_1d(geom->szcS);
  v2n = i_alloc_1d(geom->szvS);
  e2n = i_alloc_1d(geom->szeS);
  eg2n = i_alloc_1d(geom->szeS);

  n = 0;
  d->ntriangles = 0;
  for (cc = 1, n = 0; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    c2n[c] = n;
    d->points[n].x = geom->cellx[c];
    d->points[n++].y = geom->celly[c];
    d->ntriangles += geom->npe[c];
  }
  for (vv = 1; vv <= geom->n2_e2; vv++) {
    v = geom->w2_e2[vv];
    v2n[v] = n;
    d->points[n].x = geom->gridx[v];
    d->points[n++].y = geom->gridy[v];
  }

  if (inc_gst) {
    /* Set the location of the ghost cell centres. This is also      */
    /* performed in master_setghosts().                              */
    for (cc = 1; cc <= geom->nbpte1S; cc++) {
      double x, y, dist, r;
      int j;
      e = geom->bpte1S[cc];
      c = geom->bine1S[cc];
      x = geom->u1x[e] - geom->cellx[c];
      y = geom->u1y[e] - geom->celly[c];
      dist = sqrt(x * x + y * y);
      r = atan2(y, x);
      for (j = 1; j <= geom->npe[c]; j++) {
	if (e == geom->c2e[j][c]) {
	  c1 = geom->c2c[j][c];
	  geom->cellx[c1] = geom->u1x[e] + dist * cos(r);
	  geom->celly[c1] = geom->u1y[e] + dist * sin(r);
	  c2n[c1] = n;
	  d->points[n].x = geom->cellx[c1];
	  d->points[n++].y = geom->celly[c1];
	  d->ntriangles += 3;
	}
      }    
    }
    
    /* Set the location of the tangential ghost edges. This is the   */
    /* mid point joining the two cell centres adjacent to the edge.  */
    for (ee = geom->n2_e1+1; ee <= geom->x2_e1; ee++) {
      e = geom->w2_e1[ee];
      c1 = geom->e2c[e][0];
      c2 = geom->e2c[e][1];
      
      if (geom->wgst[c1] != geom->wgst[c2]) {
	geom->u1x[e] = 0.5 * (geom->cellx[c1] + geom->cellx[c2]);
	geom->u1y[e] = 0.5 * (geom->celly[c1] + geom->celly[c2]);
	e2n[e] = n;
	d->points[n].x = geom->u1x[e];
	d->points[n++].y = geom->u1y[e];
	eg2n[e] = n;
	d->points[n].x = tegl[e].x;
	d->points[n++].y = tegl[e].y;
      }
    }
    d->npoints = n;
  }
  for (n = 0; n < d->npoints; n++) {
    d->xmin = min(d->points[n].x, d->xmin);
    d->xmax = max(d->points[n].x, d->xmax);
    d->ymin = min(d->points[n].y, d->ymin);
    d->ymax = max(d->points[n].y, d->ymax);
  }
  /* Build the triangles                                       */
  n = 0;
  d->triangles = malloc(d->ntriangles * sizeof(triangle));

  geom->tri2c = i_alloc_1d(d->ntriangles);
  geom->c2tri = i_alloc_1d(geom->szcS);
  memset(geom->tri2c, 0, d->ntriangles * sizeof(int));
  for (cc = 1, n = 0; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    for (j = 1; j <= geom->npe[c]; j++) {
      triangle* t = &d->triangles[n];
      jp = (j == geom->npe[c]) ? 1 : j+1;
      t->vids[0] = c2n[c];
      t->vids[1] = v2n[geom->c2v[j][c]];
      t->vids[2] = v2n[geom->c2v[jp][c]];
      geom->tri2c[n] = c;
      geom->c2tri[c] = n;
      n++;
    }
  }
  for (cc = 1; cc <= geom->nbptS; cc++) {
    triangle* t = &d->triangles[n];
    cg = geom->bpt[cc];
    c = geom->bin[cc];
    j = geom->dbpt[cc];
    jp = (j == geom->npe[c]) ? 1 : j+1;
    t->vids[0] = c2n[cg];
    t->vids[1] = v2n[geom->c2v[j][c]];
    t->vids[2] = v2n[geom->c2v[jp][c]];
    geom->tri2c[n] = cg;
    geom->c2tri[cg] = n;
    n++;
  }
  for (ee = geom->n2_e1+1; ee <= geom->x2_e1; ee++) {
    e = geom->w2_e1[ee];
    c1 = geom->e2c[e][0];
    c2 = geom->e2c[e][1];  
    if (geom->wgst[c1] != geom->wgst[c2]) {
      for (j = 0; j <= 1; j++) {
	triangle* t = &d->triangles[n];
	c = geom->e2c[e][j];
	v = geom->e2v[e][0];
	t->vids[0] = c2n[c];
	t->vids[1] = eg2n[e];
	t->vids[2] = e2n[e];
	geom->tri2c[n] = c;
	geom->c2tri[c] = n;
	n++;
      }
    }
  }
  d->ntriangles = n;

  /* Get the neighbours                                              */
  neighbour_finder_b(d, &nei);
  d->neighbours = malloc(d->ntriangles * sizeof(triangle_neighbours));
  for (cc = 0; cc < d->ntriangles; ++cc) {
    triangle* t = &d->triangles[cc];
    triangle_neighbours* ne = &d->neighbours[cc];
    ne->tids[0] = nei[0][cc];
    ne->tids[1] = nei[1][cc];
    ne->tids[2] = nei[2][cc];
  }
  i_free_2d(nei);
  i_free_1d(c2n);
  i_free_1d(v2n);
  i_free_1d(e2n);

  /*-----------------------------------------------------------------*/
  /* Write to file                                                   */
  if (filef) {
    mesh_ofile_name(params, key);
    sprintf(buf,"%s_v.txt", key);
    if ((vf = fopen(buf, "w")) == NULL) filef = 0;
  }
  if (filef) {
    for (n = 0; n < d->ntriangles; n++) {
      triangle* t = &d->triangles[n];
      fprintf(vf, "%f %f\n", d->points[t->vids[0]].x, d->points[t->vids[0]].y);
      fprintf(vf, "%f %f\n", d->points[t->vids[1]].x, d->points[t->vids[1]].y);
      fprintf(vf, "%f %f\n", d->points[t->vids[2]].x, d->points[t->vids[2]].y);
      fprintf(vf, "%f %f\n", d->points[t->vids[0]].x, d->points[t->vids[0]].y);
      fprintf(vf, "NaN NaN\n");
    }
    fclose(vf);
  }

  /* Test the xytoi map                                              */
  if (dotest) {
    point p;
    p.x = 29934.634573;
    p.y = 12342.446276;
    d->first_id = geom->c2tri[533];
    n = delaunay_xytoi_lag(geom->d, &p, geom->d->first_id);
    if (n) {
      triangle* t = &d->triangles[n];
      printf("%d %d\n",n,geom->tri2c[n]);
      printf("%f %f\n", d->points[t->vids[0]].x, d->points[t->vids[0]].y);
      printf("%f %f\n", d->points[t->vids[1]].x, d->points[t->vids[1]].y);
      printf("%f %f\n", d->points[t->vids[2]].x, d->points[t->vids[2]].y);
    }
  }
}

/* END build_delaunay_cell()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Same as create_delaunay_cell() but for windows.                   */
/*-------------------------------------------------------------------*/
void create_delaunay_cell_w(geometry_t *window)
{
  FILE *vf;
  char buf[MAXSTRLEN];
  delaunay *d = NULL;
  int n, cc, c, vv, v;
  int np;
  point *pin;
  int filef = 0;
  /*
  window->d = master->geom->d;
  return;
  */
  np = window->b2_t + window->n2_e2;
  pin = malloc(np * sizeof(point));

  for (cc = 1, n = 0; cc <= window->b2_t; cc++) {
    c = window->w2_t[cc];
    pin[n].x = window->cellx[c];
    pin[n++].y = window->celly[c];
  }
  for (vv = 1; vv <= window->n2_e2; vv++) {
    v = window->w2_e2[vv];
    pin[n].x = window->gridx[v];
    pin[n++].y = window->gridy[v];
  }

  window->d = delaunay_build(np, pin, 0, NULL, 0, NULL);
  d = window->d;

  window->tri2c = i_alloc_1d(d->ntriangles);
  memset(window->tri2c, 0, d->ntriangles * sizeof(int));
  for (cc = 1; cc <= window->b2_t; cc++) {
    double x, y;
    c = window->w2_t[cc];
    x = window->cellx[c];
    y = window->celly[c];

    for (n = 0; n < d->ntriangles; n++) {
      triangle* t = &d->triangles[n];
      if ((x == d->points[t->vids[0]].x || x == d->points[t->vids[1]].x ||
	   x == d->points[t->vids[2]].x) && 
	  (y == d->points[t->vids[0]].y || y == d->points[t->vids[1]].y ||
	   y == d->points[t->vids[2]].y))
	window->tri2c[n] = c;
    }
  }
  if (filef) {
    sprintf(buf,"window%d_v.txt", filef);
    if ((vf = fopen(buf, "w")) == NULL) filef = 0;
  }
  if (filef) {
    for (n = 0; n < d->ntriangles; n++) {
      triangle* t = &d->triangles[n];
      fprintf(vf, "%f %f\n", d->points[t->vids[0]].x, d->points[t->vids[0]].y);
      fprintf(vf, "%f %f\n", d->points[t->vids[1]].x, d->points[t->vids[1]].y);
      fprintf(vf, "%f %f\n", d->points[t->vids[2]].x, d->points[t->vids[2]].y);
      fprintf(vf, "%f %f\n", d->points[t->vids[0]].x, d->points[t->vids[0]].y);
      fprintf(vf, "NaN NaN\n");
    }
    fclose(vf);
  }
}

/* END create_delaunay_cell()                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Set up a delaunay triangulation for interpolation of cell         */
/* centered variables. The triangulation is based on grid centres    */
/* only.                                                             */
/*-------------------------------------------------------------------*/
void create_delaunay_cent(geometry_t *geom, parameters_t *params)
{
  delaunay *d;
  int n, cc, c, vv, v;
  int np;
  point *pin;

  np = geom->b2_t;
  pin = malloc(np * sizeof(point));

  for (cc = 1, n = 0; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    pin[n].x = geom->cellx[c];
    pin[n++].y = geom->celly[c];
  }
  if (d) free((delaunay *)d);
  /*d = delaunay_voronoi_build(np, pin, 0, NULL, 0, NULL, 0.0);*/
  d = delaunay_build(np, pin, 0, NULL, 0, NULL);

}

/* END create_delaunay_cent()                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Interpolates unstructued data horizontally onto a set of points   */
/* using the i_rule method.                                          */
/*-------------------------------------------------------------------*/
void interp_us(geometry_t *geom,    /* Model geometry                */
	       double *vals,        /* Cell centre values            */
	       int nvec,            /* Number of values to interp    */
	       int *vec,            /* Locations to store values     */
	       double *locx,        /* x locations to interp         */
	       double *locy,        /* y locations to interp         */
	       double *ret          /* Array of interpolated values  */
	       )
{
  GRID_SPECS *gs = NULL;
  int cc, c;

  /* Initialise                                                      */
  gs = grid_interp_init(geom->cellx, geom->celly, vals, geom->b2_t, geom->i_rule);

  /* Interpolate the given points                                    */
  for (cc = 1; cc <= nvec; cc++) {
    c = vec[cc];
    ret[c] = grid_interp_on_point(gs, locx[c], locy[c]);
  }

  grid_specs_destroy(gs);
}

/* END interp_us()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Interpolates unstructued edge data horizontally onto a set of     */
/* points using the i_rule method.                                   */
/*-------------------------------------------------------------------*/
void interp_edge_us(geometry_t *geom, /* Model geometry              */
		    double *vals,     /* Cell edge values            */
		    int nvec,         /* Number of values to interp  */
		    int *vec,         /* Locations to store values   */
		    double *locx,     /* x locations to interp       */
		    double *locy,     /* y locations to interp       */
		    double *ret,      /* Interpolated values         */
		    char *i_rule
		    )
{
  GRID_SPECS *gs = NULL;
  int cc, c;

  /* Initialise                                                      */
  gs = grid_interp_init(geom->u1x, geom->u1y, vals, geom->b2_e1, i_rule);

  /* Interpolate the given points                                    */
  for (cc = 1; cc <= nvec; cc++) {
    c = vec[cc];
    ret[c] = grid_interp_on_point(gs, locx[c], locy[c]);
  }

  grid_specs_destroy(gs);
}

/* END interp_edge_us()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Writes a list of segments locations at the grid perimeter from    */
/* which OBC locations may be extracted.                             */
/*-------------------------------------------------------------------*/
int get_limit_obc(parameters_t *params, 
		  int ns2,
		  int **neic
		  )
{
  FILE *fp;
  int n, cc, j, jj, ee, e, m = 0;
  int npe;
  double x1, x2, y1, y2;
  int filef = 1;
  int verbose = 0;
  int ret = 0;

  if (filef)
    if ((fp = fopen("boundary.txt", "w")) == NULL)
      filef = 0;
  
  /*-----------------------------------------------------------------*/
  /* Count the edges within each boundary limits                     */
  for (n = 0; n < params->nbu; n++) {
    open_bdrys_t *open = params->open[n];
    if (open->minlon == NOTVALID || open->minlat == NOTVALID ||
	open->maxlon == NOTVALID || open->maxlat == NOTVALID)
      continue;
    params->nptsu[n] = 0;
    for (cc = 1; cc <= ns2; cc++) {
      npe = params->npe2[cc];
      for (j = 1; j <= npe; j++) {
	if (!neic[j][cc]) {
	  for (ee = 0; ee < open->nedges; ee++) {
	    e = open->edges[ee];
	    if (j != e) continue;
	    x1 = params->x[cc][e];
	    y1 = params->y[cc][e];
	    x2 = params->x[cc][jp(e, npe)];
	    y2 = params->y[cc][jp(e, npe)];
	    if (x1 >= open->minlon && x1 <= open->maxlon && 
		y1 >= open->minlat && y1 <= open->maxlat &&
		x2 >= open->minlon && x2 <= open->maxlon && 
		y2 >= open->minlat && y2 <= open->maxlat) {
	      params->nptsu[n]++;
	      ret = 1;
	    }
	  }
	}
      }
    }
    if (!params->nptsu[n]) 
      hd_quit("get_limit_obc: No boundary edges found for boundary %d, limits (%f,%f)-(%f-%f)\n", n,
	      open->minlon, open->minlat, open->maxlon,open->maxlat);
    if (params->nptsu[n] > m) m = params->nptsu[n];
  }

  /*-----------------------------------------------------------------*/
  /* Allocate memory                                                 */
  if (m) {
    params->locu = i_alloc_2d(m, params->nbu);
    params->posx = d_alloc_3d(2 ,m, params->nbu);
    params->posy = d_alloc_3d(2 ,m, params->nbu);
  }

  /*-----------------------------------------------------------------*/
  /* Assign the edge locations to the boundary vectors               */
  for (n = 0; n < params->nbu; n++) {
    open_bdrys_t *open = params->open[n];
    if (open->minlon == NOTVALID || open->minlat == NOTVALID ||
	open->maxlon == NOTVALID || open->maxlat == NOTVALID)
      continue;
    params->nptsu[n] = 0;
    for (cc = 1; cc <= ns2; cc++) {
      for (j = 1; j <= params->npe2[cc]; j++) {
	if (!neic[j][cc]) {
	  for (ee = 0; ee < open->nedges; ee++) {
	    e = open->edges[ee];
	    if (j != e) continue;
	    x1 = params->x[cc][e];
	    y1 = params->y[cc][e];
	    x2 = params->x[cc][jp(e, npe)];
	    y2 = params->y[cc][jp(e, npe)];
	    if (x1 >= open->minlon && x1 <= open->maxlon && 
		y1 >= open->minlat && y1 <= open->maxlat &&
		x2 >= open->minlon && x2 <= open->maxlon && 
		y2 >= open->minlat && y2 <= open->maxlat) {
	      params->locu[n][params->nptsu[n]] = cc;
	      params->posx[n][params->nptsu[n]][0] = x1;
	      params->posx[n][params->nptsu[n]][1] = x2;
	      params->posy[n][params->nptsu[n]][0] = y1;
	      params->posy[n][params->nptsu[n]][1] = y2;
	      params->nptsu[n]++;
	    }
	  }
	}
      }
    }

    if (verbose) {
      printf("\nBOUNDARY%d     %d\n", n, params->nptsu[n]);
      for (cc = 0; cc < params->nptsu[n]; cc++) {
	printf("%d (%f,%f)-(%f,%f)\n", params->locu[n][cc], params->posx[n][cc][0], params->
	       posy[n][cc][0], params->posx[n][cc][1], params->posy[n][cc][1]);
      }
    }
    if (filef) {
      for (cc = 0; cc < params->nptsu[n]; cc++) {
	fprintf(fp, "%f %f\n", params->posx[n][cc][0], params->posy[n][cc][0]);
	fprintf(fp, "%f %f\n", params->posx[n][cc][1], params->posy[n][cc][1]);
	fprintf(fp, "NaN NaN\n");
      }
    }
  }
  if (filef ) fclose(fp);
  return(ret);
}

/* END get_limit_obc()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Gets the direction that maps to the interior. This is defined as  */
/* a wet (non-boundary, non-ghost) cell that has its opposite edge   */
/* corresponding to a ghost cell.                                    */
/*-------------------------------------------------------------------*/
int get_bind(geometry_t *sgrid, int c, int *mask)
{
  int ret, retn, retg, j, jj;
  int cn, cns, co;
  int cs = sgrid->m2d[c];
  int npe = sgrid->npe[cs];

  ret = retn = retg = 0;
  for (j = 1; j <= npe; j++) {
    cn = sgrid->c2c[j][c];
    cns = sgrid->m2d[cn];
    if (!sgrid->wgst[cn]) {
      if (usejo)
	jj = oedge(npe, j);
      else
	jj = jocw(sgrid, c, j);
      jj = min(jj, npe);
      co = sgrid->c2c[jj][c];
      if (!mask[cns] && sgrid->wgst[co]) ret = j;

      /* If an interior map can't be found, then find the direction  */
      /* of a non-ghost cell whose opposite is a ghost.              */
      if (sgrid->wgst[co]) retn = j;
    } else {
      /* If an interior map still can't be found, use the opposite   */
      /* direction to any ghost cells.                               */
      if (usejo)
	retg = oedge(npe, j);
      else
	retg = jocw(sgrid, c, j);
    }
  }

  if (!ret) ret = retn;
  /* If the interior direction is still a ghost cell, use the        */
  /* average of non OBC, non ghost directions.                       */
  if (mask[sgrid->c2c[ret][c]]) {
    retn = co = 0;
    for (j = 1; j <= npe; j++) {
      cn = sgrid->c2c[j][c];
      cns = sgrid->m2d[cn];
      if (!sgrid->wgst[cn] && !mask[cns]) {
	retn += j;
	co++;
      }
    }
    ret = (co) ? retn / co : 0;
  }
  if (!ret && DEBUG("init_m")) {
    dlog("init_m", "\nget_bind: Can't find interior boundary direction for cell %d\n", c);
  }
  if (ret <=0 || ret > npe) ret = retg;
  return(ret);
}

/* END get_bind()                                                    */
/*-------------------------------------------------------------------*/



int get_nve(int npe)
{
  switch (npe) {
  case 3:  /* Delaunay trianges */
    return(6);
    break;

  case 4:  /* Squares */
    return(4);
    break;

  case 6:  /* Voronoi hexagons */
    return(3);
    break;
  }
}


void vertex_map_4(geometry_t *sgrid)
{
  int cc, c, co;
  int vv, v;
  int i, j;

  if (!(sgrid->us_type & US_IJ)) {
    sgrid->v2i = NULL;
    sgrid->v2j = NULL;
    sgrid->v2ijk = NULL;
    return;
  }

  sgrid->v2i = i_alloc_1d(sgrid->szvS);
  sgrid->v2j = i_alloc_1d(sgrid->szvS);

  for (vv = 1; vv <= sgrid->n2_e2; vv++) {
    v = sgrid->w2_e2[vv];
    sgrid->v2i[v] = NOTVALID;
    sgrid->v2j[v] = NOTVALID;
  }

  for (cc = 1; cc <= sgrid->b2_t; cc++) {
    c = sgrid->w2_t[cc];
    i = sgrid->s2i[c];
    j = sgrid->s2j[c];
    v = sgrid->c2v[1][c];
    sgrid->v2i[v] = i;
    sgrid->v2j[v] = j;
    v = sgrid->c2v[2][c];
    if (sgrid->v2i[v] == NOTVALID && sgrid->v2j[v] == NOTVALID) {
      sgrid->v2i[v] = i;
      sgrid->v2j[v] = j + 1;
    }

    v = sgrid->c2v[3][c];
    if (sgrid->v2i[v] == NOTVALID && sgrid->v2j[v] == NOTVALID) {
      sgrid->v2i[v] = i + 1;
      sgrid->v2j[v] = j + 1;
    }

    v = sgrid->c2v[4][c];
    if (sgrid->v2i[v] == NOTVALID && sgrid->v2j[v] == NOTVALID) {
      sgrid->v2i[v] = i + 1;
      sgrid->v2j[v] = j;
    }
  }

  sgrid->v2ijk = i_alloc_1d(sgrid->szv);
  for (vv = 1; vv <= sgrid->n3_e2; vv++) {
    v = sgrid->w3_e2[vv];
    sgrid->v2ijk[v] = NOTVALID;
  }
  for (cc = 1; cc <= sgrid->b3_t; cc++) {
    c = co = sgrid->w3_t[cc];
    v = sgrid->c2v[1][co];
    sgrid->v2ijk[v] = c;

    c = sgrid->c2c[2][c];
    v = sgrid->c2v[2][co];
    if (sgrid->v2ijk[v] == NOTVALID)
      sgrid->v2ijk[v] = c;

    c = sgrid->c2c[3][c];
    v = sgrid->c2v[3][co];
    if (sgrid->v2ijk[v] == NOTVALID)
      sgrid->v2ijk[v] = c;

    c = sgrid->c2c[3][co];
    v = sgrid->c2v[4][co];
    if (sgrid->v2ijk[v] == NOTVALID)
      sgrid->v2ijk[v] = c;
  }
}


int find_vertex(int c, double x, double y, double **xloc, double **yloc, int *mask, int ns2, int npe)
{
  int cc, j;
  int found = 0;

  for (cc = 1; cc <= ns2; cc++) {
    if(c == cc) continue;
    if (mask[cc]) {
      for (j = 1; j <= npe; j++) {
	if (x == xloc[cc][j] && y == yloc[cc][j]) {
	  found = 1;
	}
      }
    }
  }
  return(found);
}
