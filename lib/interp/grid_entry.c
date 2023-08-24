/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  File: model/lib/grid/grid_entry.c
 *  
 *  Description:
 *  Public API to the Grid interpolation library
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: grid_entry.c 7148 2022-07-07 00:26:44Z her127 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "poly.h"
#include "emsalloc.h"
#include "gridlib.h"

int grid_entry_verbose = 0;

static int nnf = 0;

/*
 * Writes value to the output buffer, if supplied else prints to stdout
 */
static void write_output(FILE *outfile, double val, int *index) {
  if (outfile == NULL) {
    outfile = stdout;
  }

  if (isnan(val)) 
    fprintf(outfile, "NaN\n");
  else
    fprintf(outfile, "%.2f\n", val);
  
  /* 
   * This is for error bounds checking.
   * The idea would be to allow for filling pre-allocated buffers
   * as well as writing to a file handle.
   */
  (*index)++;
}

/*
 * Switch statement for the interpolator
 */
static void set_interpolation_function(GRID_SPECS *gs, gridmap *gm, delaunay *di)
{
  interp_func  interpolate_point = NULL;
  interp_func2 interpolate_point2 = NULL;
  rebuild_func rebuild           = NULL;
  rebuild_func2 rebuild2         = NULL;
  void        *interpolator      = NULL;
  delaunay    *d                 = NULL;
  INTERP_RULE rule               = gs->type;

  if (rule == GRID_CSA) {
    interpolator = csa_create();
    csa_addpoints(interpolator, gs->nbathy, gs->pbathy);
    csa_calculatespline(interpolator);
    interpolate_point = (void (*)(void*, point *)) csa_approximatepoint;
    rebuild = (void (*)(void*, point *)) csa_refindprimarycoeffs;
  } else if (rule == GRID_AVERAGE) {
    if (gm == NULL)
      grid_quit("Error: gridmap is needed for AVERAGING\n");
    interpolator = ga_create(gm);
    ga_addpoints(interpolator, gs->nbathy, gs->pbathy);
    interpolate_point = (void (*)(void*, point *)) ga_getvalue;
    gs->ppe = 1;
  } else {
    /*
     * triangulate 
     */
    if (grid_entry_verbose) {
      fprintf(stderr, "## triangulating...");
      fflush(stdout);
    }
    if (di != NULL)
      d = di;
    else
      d = delaunay_build(gs->nbathy, gs->pbathy, 0, NULL, 0, NULL);
    if (grid_entry_verbose) {
      fprintf(stderr, "done\n");
      fflush(stderr);
    }
    
    if (rule == GRID_NN_SIBSON || rule == GRID_NN_NONSIBSONIAN) {
      if (!nnf) {
	interpolator = nnpi_create(d);
	if (rule == GRID_NN_SIBSON)
	  nnpi_set_rule(interpolator, SIBSON);
	else
	  nnpi_set_rule(interpolator, NON_SIBSONIAN);
	interpolate_point = (void (*)(void*, point *)) nnpi_interpolate_point;
      } else {
	interpolator = nnhpi_create(d, d->npoints);
	if (rule == GRID_NN_SIBSON)
	  nnhpi_set_rule(interpolator, SIBSON);
	else
	  nnhpi_set_rule(interpolator, NON_SIBSONIAN);
	interpolate_point = (void (*)(void*, point *)) nnhpi_interpolate;
	rebuild = (void (*)(void*, point *)) nnhpi_modify_data;
      }
    } else if (rule == GRID_NNA_SIBSON || rule == GRID_NNA_NONSIBSONIAN) {
      interpolator = nnhpi_create(d, d->npoints);
      if (rule == GRID_NNA_SIBSON)
	nnhpi_set_rule(interpolator, SIBSON);
      else
	nnhpi_set_rule(interpolator, NON_SIBSONIAN);
      interpolate_point = (void (*)(void*, point *)) nnhpi_interpolate;
      rebuild = (void (*)(void*, point *)) nnhpi_modify_data;
    } else if (rule == GRID_LINEAR) {
      interpolator = lpi_build(d);
      rebuild = (void (*)(void*, point *)) lpi_rebuild;
      interpolate_point = (void (*)(void*, point *)) lpi_interpolate_point;
    } else if (rule == GRID_BL) {
      interpolator = bl_build(d);
      rebuild = (void (*)(void*, point *)) bl_rebuild;
      rebuild2 = (void (*)(void*, void*, point *)) bl_rebuild2;
      interpolate_point = (void (*)(void*, point *)) bl_interpolate_point;
      interpolate_point2 = (void (*)(void*, void*, point *)) bl_interpolate_point2;
    } else if (rule == GRID_BAL) {
      interpolator = bal_build(d);
      rebuild = (void (*)(void*, point *)) bal_rebuild;
      rebuild2 = (void (*)(void*, void*, point *)) bal_rebuild2;
      interpolate_point = (void (*)(void*, point *)) bal_interpolate_point;
      interpolate_point2 = (void (*)(void*, void*, point *)) bal_interpolate_point2;
    } else if (rule == GRID_LSQQ) {
      interpolator = lsqq_build(d);
      rebuild = (void (*)(void*, point *)) lsqq_rebuild;
      interpolate_point = (void (*)(void*, point *)) lsqq_interpolate_point;
    } else if (rule == GRID_LSQL) {
      interpolator = lsql_build(d);
      rebuild = (void (*)(void*, point *)) lsql_rebuild;
      interpolate_point = (void (*)(void*, point *)) lsql_interpolate_point;
    } else if (rule == GRID_NRST) {
      interpolator = nrst_build(d);
      rebuild = (void (*)(void*, point *)) nrst_rebuild;
      interpolate_point = (void (*)(void*, point *)) nrst_interpolate_point;
    }
  }

  // Assign points
  gs->interpolate_point = interpolate_point;
  gs->interpolate_point2 = interpolate_point2;
  gs->rebuild           = rebuild;
  gs->rebuild2          = rebuild2;
  gs->interpolator      = interpolator;
  gs->d                 = d;
}

/*
 * Char to enum interpolation method
 */
static INTERP_RULE interp_rule_from_char(char *rule)
{
  if (rule == NULL || strlen(rule) == 0)
    return(GRID_LINEAR); // This is the default
  else if (strcasecmp("linear", rule) == 0)
    return(GRID_LINEAR);
  else if (strcasecmp("cubic", rule) == 0)
    return(GRID_CSA);
  else if (strcasecmp("nn_sibson", rule) == 0)
    return(GRID_NN_SIBSON);
  else if (strcasecmp("nn_non_sibson", rule) == 0)
    return(GRID_NN_NONSIBSONIAN);
  else if (strcasecmp("nna_sibson", rule) == 0)
    return(GRID_NNA_SIBSON);
  else if (strcasecmp("nna_non_sibson", rule) == 0)
    return(GRID_NNA_NONSIBSONIAN);
  else if (strcasecmp("average", rule) == 0)
    return(GRID_AVERAGE);
  else if (strcasecmp("quadratic", rule) == 0)
    return(GRID_LSQQ);
  else if (strcasecmp("linearlsq", rule) == 0)
    return(GRID_LSQL);
  else if (strcasecmp("bilinear", rule) == 0)
    return(GRID_BL);
  else if (strcasecmp("baylinear", rule) == 0)
    return(GRID_BAL);
  else if (strcasecmp("nearest", rule) == 0)
    return(GRID_NRST);
  else
    grid_quit("interp_rule_from_char: no such interpolation rule '%s'\n", rule);

  return(GRID_LINEAR);
}


/*
 * Allocates memory for the grid specs struct and populates with some
 * default values
 */
GRID_SPECS *grid_spec_create(void)
{
  GRID_SPECS *gs = (GRID_SPECS *)alloc_1d(1, sizeof(GRID_SPECS));
  
  gs->nbathy = -1;
  gs->pbathy = NULL;
  gs->type   = GRID_CSA;
  gs->gn     = NULL;
  gs->mask   = NULL;
  gs->node_type = NT_DD;
  gs->ppe    = PPE_DEF;
  gs->zmin   = ZMIN_DEF;
  gs->zmax   = ZMAX_DEF;
  gs->cell_i = -1;
  gs->cell_j = -1;
  gs->index_space = 0;
  gs->id = 0;

  gs->output_fileId = NULL;

  gs->interpolate_point = NULL;
  gs->interpolator      = NULL;
  gs->d                 = NULL;

  gs->destroy_pbathy = 0;
  gs->destroy_delaunay = 1;

  return gs;
}

void grid_spec_init(GRID_SPECS *gs)
{
  gs->nbathy = -1;
  gs->pbathy = NULL;
  gs->type   = GRID_CSA;
  gs->gn     = NULL;
  gs->mask   = NULL;
  gs->node_type = NT_DD;
  gs->ppe    = PPE_DEF;
  gs->zmin   = ZMIN_DEF;
  gs->zmax   = ZMAX_DEF;
  gs->cell_i = -1;
  gs->cell_j = -1;
  gs->index_space = 0;
  gs->id = 0;

  gs->output_fileId = NULL;

  gs->interpolate_point = NULL;
  gs->interpolator      = NULL;
  gs->d                 = NULL;

  gs->destroy_pbathy = 0;
  gs->destroy_delaunay = 1;
}

void grid_specs_destroy(GRID_SPECS *gs)
{
  INTERP_RULE rule = gs->type;

  /*
   * clean up, just because 
   */
  if (rule == GRID_CSA)
    csa_destroy(gs->interpolator);
  else if (rule == GRID_AVERAGE)
    ga_destroy(gs->interpolator);
  else {
    if (rule == GRID_NN_SIBSON || rule == GRID_NN_NONSIBSONIAN) {
      if (!nnf)
	nnpi_destroy(gs->interpolator);
      else
	nnhpi_destroy(gs->interpolator);
    }
    else if (rule == GRID_NNA_SIBSON || rule == GRID_NNA_NONSIBSONIAN)
      nnhpi_destroy(gs->interpolator);
    else if (rule == GRID_LINEAR)
      lpi_destroy(gs->interpolator);
    else if (rule == GRID_LINEAR)
      lpi_destroy(gs->interpolator);
    else if (rule == GRID_BL)
      bl_destroy(gs->interpolator);
    else if (rule == GRID_BAL)
      bal_destroy(gs->interpolator);
    else if (rule == GRID_LSQQ)
      lsqq_destroy(gs->interpolator);
    else if (rule == GRID_LSQL)
      lsql_destroy(gs->interpolator);
    else if (rule == GRID_NRST)
      nrst_destroy(gs->interpolator);
  }
  
  // free the points array if flag is set
  if (gs->destroy_pbathy && gs->pbathy != NULL)
    free(gs->pbathy);

  if (gs->destroy_delaunay && gs->d != NULL)
    delaunay_destroy(gs->d);

  free_1d(gs);
}

/*
 * Main entry function for standalone executable, the library should
 * call the _interp_on_point routine
 */
int grid_interp(GRID_SPECS *gs)
{
  int    nbathy = gs->nbathy;
  point *pbathy = gs->pbathy;
  gridnodes *gn = gs->gn;

  int **mask = gs->mask;
  NODETYPE nt = gs->node_type;

  int ppe = gs->ppe;
  
  double zmin = gs->zmin;
  double zmax = gs->zmax;

  int indexspace = gs->index_space;
    
  gridmap *gm = NULL;

  /* Some sanity checking */
  if (nbathy < 3)
    grid_quit("less than 3 input bathymetry values\n");
  if (ppe <= 0 || ppe > PPE_MAX)
    grid_quit("number of points per edge specified = %d greater than %d\n", 
	                                                       ppe, PPE_MAX);
  if (zmin >= zmax)
    grid_quit("min depth = %.3g > max depth = %.3g\n", zmin, zmax);
  if (nt != NT_DD && nt != NT_COR)
    grid_quit("unsupported node type\n");

  /*
   * transform grid nodes to corner type 
   */
  if (gn != NULL && nt != NT_COR) {
    gridnodes* newgn = gridnodes_transform(gn, NT_COR);
    
    gridnodes_destroy(gn);
    gn = newgn;

    // update the specs as this memory is freed outside this function
    gs->gn = gn;
  }

  /*
   * build the grid map for physical <-> index space conversions
   */
  if (gn != NULL)
    gm = gridmap_build(gridnodes_getnce1(gn), gridnodes_getnce2(gn),
		       gridnodes_getx(gn), gridnodes_gety(gn));
  
  /*
   * convert bathymetry to index space if necessary 
   */
  if (indexspace) {
    point* newpbathy = malloc(nbathy * sizeof(point));
    int newnbathy = 0;
    int ii;
    
    for (ii = 0; ii < nbathy; ++ii) {
      point* p = &pbathy[ii];
      point* newp = &newpbathy[newnbathy];
      double ic, jc;
      
      if (gridmap_xy2fij(gm, p->x, p->y, &ic, &jc)) {
	newp->x = ic;
	newp->y = jc;
	newp->z = p->z;
	newnbathy++;
      }
    }
    
    free(pbathy);
    pbathy = newpbathy;
    nbathy = newnbathy;
  }

  /*
   * create interpolator 
   */
  set_interpolation_function(gs, gm, NULL);
  
  /*
   * main cycle -- over grid cells 
   */
  if (gn != NULL) {
    double** gx = gridnodes_getx(gn);
    int jmin, jmax, imin, imax, i=gs->cell_i, j=gs->cell_j;
    int output_index = 0;

    if (i < 0) {
      imin = 0;
      imax = gridnodes_getnce1(gn) - 1;
      jmin = 0;
      jmax = gridnodes_getnce2(gn) - 1;
    } else {
      if (grid_entry_verbose)
	fprintf(stderr, "## calculating depth for cell (%d,%d)\n", i, j);
      imin = i;
      imax = i;
      jmin = j;
      jmax = j;
    }
    
    for (j = jmin; j <= jmax; ++j) {
      for (i = imin; i <= imax; ++i) {
	double sum = 0.0;
	int count = 0;
	int ii, jj;
	double val;

	if ((mask != NULL && mask[j][i] == 0) || 
	    isnan(gx[j][i])         || 
	    isnan(gx[j + 1][i])     || 
	    isnan(gx[j][i + 1])     || 
	    isnan(gx[j + 1][i + 1])) {
	  write_output(gs->output_fileId, NaN, &output_index);
	  continue;
	}
	
	for (ii = 0; ii < ppe; ++ii) {
	  for (jj = 0; jj < ppe; ++jj) {
	    double fi = (double)i + 0.5/(double)ppe * (1.0 + 2.0 * (double)ii);
	    double fj = (double)j + 0.5/(double)ppe * (1.0 + 2.0 * (double)jj);
	    point p;
	    
	    if (!indexspace)
	      gridmap_fij2xy(gm, fi, fj, &p.x, &p.y);
	    else {
	      p.x = fi;
	      p.y = fj;
	    }
	    
	    gs->interpolate_point(gs->interpolator, &p);
	    
	    if (isnan(p.z))
	      continue;
	    else if (p.z < zmin)
	      p.z = zmin;
	    else if (p.z > zmax)
	      p.z = zmax;
	    
	    sum += p.z;
	    count++;
	  }
	}

	if (count == 0) {
	  val = NaN;
	} else {
	  val =  sum / (double)count;
	}

	/* Output bounds check */
	if (output_index < gridnodes_getnce1(gn)*gridnodes_getnce2(gn))
	  write_output(gs->output_fileId, val, &output_index);
	else
	  grid_quit("grid_entry(): output index %d seems to be outside the grid size of %d x %d = %d\n", gridnodes_getnce1(gn), gridnodes_getnce2(gn), gridnodes_getnce1(gn)*gridnodes_getnce2(gn));
      }
    }
  } else {
    if (gs->type == GRID_NN_SIBSON || gs->type == GRID_NN_NONSIBSONIAN ||
	gs->type == GRID_NNA_SIBSON || gs->type == GRID_NNA_NONSIBSONIAN) {
      nnpi_interpolate_points(gs->nbathy, gs->pbathy, WMIN_DEF,
			      gs->npout,  gs->pout);
    } else {
      int n;
      point p;
      // Loop over all points
      for (n=0; n<gs->npout; n++) {
	gs->pout[n].z =  grid_interp_on_point(gs, gs->pout[n].x, gs->pout[n].y);
      }
    }
    
    points_write(gs->output_fileId, gs->npout, gs->pout);
  }

  // cleanup
  if (gm != NULL)
    gridmap_destroy(gm);

  return 0;

} // grid_interp

/*
 * Grid interpolation on the given point
 */
double grid_interp_on_point(GRID_SPECS *gs, double xcoord, double ycoord)
{
  point p;

  p.x = xcoord;
  p.y = ycoord;
  p.z = gs->id;

  gs->interpolate_point(gs->interpolator, &p);

  /*
   * Does this need cleanup here for zmin/zmax and or nan's etc?? -FR
   */
  return p.z;
}

double grid_interp_on_point2(GRID_SPECS **gs, int k1, int k2, double xcoord, double ycoord)
{
  point p;

  p.x = xcoord;
  p.y = ycoord;
  p.z = gs[k1]->id;
  /*printf("k %d %d\n",k1,k2);*/
  gs[k1]->interpolate_point2(gs[k1]->interpolator, gs[k2]->interpolator, &p);

  /*
   * Does this need cleanup here for zmin/zmax and or nan's etc?? -FR
   */
  return p.z;
}

/*
 * Grid interpolation on the given point in 3D
 */
double grid_interp_on_point3d(GRID_SPECS **gs, double xcoord, double ycoord, double depth, double bot)
{
  point p;
  int k1, k2;
  int nz = gs[0]->nz - 1;
  double *z = gs[0]->z;
  double v, v1, v2;

  /* Set up the point */
  p.x = xcoord;
  p.y = ycoord;
  p.z = depth;

  /* Get the depth levels bracketing depth. */
  /* First find the depth level below depth. */
  k2 = nz;
  while (depth < max(z[k2], bot)) {
    k2--;
    if (k2 <= 0)break;
  }
  /* Set the depth level above depth. If k2 is the surface */
  /* layer then set k1 = k2. Note that the surface layer level */
  /* is 0.5(eta+gridz) hence varies in time. */
  k1 = (k2 == nz) ? k2 : k2+1;

  if (bot < z[k1] && bot > z[k2] && depth > bot) k2++;
  /*if (z[k2] < bot && k2 == nz-1) k2++;*/

  /* Interpolate on the layer above depth */
  gs[k1]->interpolate_point(gs[k1]->interpolator, &p);
  v1 = p.z;

  if (k2 == nz || depth < bot) {
    /* Return if depth is above z[nz] or below bot */
    v = v1;
  } else {
    /* Interpolate on the layer below depth, then interpolate */
    /* vertically. */
    gs[k2]->interpolate_point(gs[k2]->interpolator, &p);
    v2 = p.z;
    v = (v2 - v1) * (depth - z[k1]) / (z[k2] - z[k1]) + v1;
  }
  return(v);
}

/*
 * Convenient wrapper function
 *   o) Allocates the GRID_SPECS structure
 *   o) populates the interp-type
 *   o) Allocates the points array and initialises the x,y & z values
 */
GRID_SPECS *grid_interp_init(double *x, double *y, double *z, int npoints,
			     char *rule)
{
  int i;
  GRID_SPECS *gs = grid_spec_create();
  point      *p  = (point *)calloc(npoints, sizeof(point));
  if (p == NULL)
    grid_quit("grid_interp_init: could not allocate points array\n");

  gs->type = interp_rule_from_char(rule);
  
  // fill the points array
  for (i=0; i<npoints; i++) {
    p[i].x = x[i];
    p[i].y = y[i];
    p[i].z = z[i];
  }
  
  gs->nbathy = npoints;
  gs->pbathy = p;
  gs->destroy_pbathy = 1;

  /*
   * Sets-up the internal interpolation function pointer
   * Note: We're passing NULL for the gridmap which won't work for
   *       GRID_AVERAGE
   */
  set_interpolation_function(gs, NULL, NULL);

  return(gs);
}


void grid_interp_init_t(GRID_SPECS *gs, delaunay *d, char *rule, int var)
{
  int i;
  int npoints = d->npoints;
  point *p = (point *)calloc(npoints, sizeof(point));

  if (p == NULL)
    grid_quit("grid_interp_init: could not allocate points array\n");

  gs->type = interp_rule_from_char(rule);

  // fill the points array
  for (i=0; i<npoints; i++) {
    d->points[i].z = d->points[i].v[var];
  }

  gs->nbathy = d->npoints;
  gs->pbathy = d->points;
  gs->destroy_pbathy = 1;
  gs->destroy_delaunay = 0;

  /*
   * Sets-up the internal interpolation function pointer
   * Note: We're passing NULL for the gridmap which won't work for
   *       GRID_AVERAGE
   */
  set_interpolation_function(gs, NULL, d);
  /*gs->d = d;*/
}


/*
 * Re-initializes a grid_spec function with new data.
 * Position data remains the same.
 */
void grid_interp_reinit(GRID_SPECS *gs, double *z, int npoints)
{
  int i;
  point      *p  = gs->pbathy;

  /* Fill the points array. Note this is designed to be used with COMPAS */
  /* arrays, which start at index 1, hence use z[i+1]. */

  for (i=0; i<npoints; i++) {
    p[i].z = z[i+1];
  }  

  gs->rebuild(gs->interpolator, p);
}

// EOF
