/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/outputs/roms_grid.c
 *  
 *  Description:
 *  Routines for conversion of EMS grids to ROMS grids.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: roms_grid.c 5873 2018-07-06 07:23:48Z riz008 $
 *
 */

#include <stdlib.h>
#include <netcdf.h>
#include <string.h>
#include "hd.h"

/*-------------------------------------------------------------------*/
/* Flags                                                             */
/*-------------------------------------------------------------------*/
/* ROMS cell flags                                                   */
#define R_RHO   1
#define R_PSI   2
#define R_U     4
#define R_V     8

/* Bathymetry and layer structure                                    */
#define MAXGRIDZ 1e20
#define NOTVALID -9999
#define LANDCELL 99
#define MAXDEPTH 9000
#define NOTWET (SOLID|OUTSIDE)

/* Defined in mom_grid */
int find_name_marvl(dump_data_t *dumpdata, char *ofile, char *tag, int *p_n);

/* From df_roms */
void df_roms_write_ic(dump_data_t *dumpdata, char *name);
void df_roms_write_bdry(dump_data_t *dumpdata,char *name,double t0,double t1);
void roms_compute_vert(int nz, float *sc, float *Cs, int isW);
void write_roms_bathy(int fid, int nce2, int nce1, double bmin, double **h);

/*------------------------------------------------------------------*/
/* Private routines                                                 */
/*------------------------------------------------------------------*/
romsgrid_t *romsgrid_alloc(void);
void alloc_grid_mem_roms(romsgrid_t *romsgrid, int nce1, int nce2);
void alloc_grid_mem_vert_roms(romsgrid_t *romsgrid);
int dump_create_roms(romsgrid_t *romsgrid, char *name);
void write_roms_attributes(romsgrid_t *romsgrid, int cdfid);
void write_text_att(int cdfid, int varid, const char *name,
                    const char *text);
void write_roms_metrics(romsgrid_t *romsgrid, int mode);
void close_roms_netcdf(romsgrid_t *romsgrid);

/*-------------------------------------------------------------------*/
/* Reads mask from ROMS grid file and masks bathy accordingly        */
/*-------------------------------------------------------------------*/
void mask_bathy_from_roms(char *fname, parameters_t *params, romsgrid_t *romsgrid)
{
  int fid;
  int i, j, n;
  size_t start[4];
  size_t count[4];
  size_t xi_rho, eta_rho;
  size_t xi_psi, eta_psi;
  size_t xi_u, eta_u;
  size_t xi_v, eta_v;
  double **mask_rho;
  double **mask_psi;
  double **mask_u;
  double **mask_v;
  int nce1, nce2;

  /* Note : We're assuming this is always true */
  // int inface = 1;

  memset(start, 0, 4*sizeof(size_t));
  memset(count, 0, 4*sizeof(size_t));

  /*-----------------------------------------------------------------*/
  /* Open the dump file for reading                                  */
  ncw_open(fname, NC_NOWRITE, &fid);

  /*
   * Get ROMS dimensions
   */
  ncw_inq_dimlen(fname, fid, ncw_dim_id(fid, "xi_rho"),  &xi_rho);
  ncw_inq_dimlen(fname, fid, ncw_dim_id(fid, "xi_psi"),  &xi_psi);
  ncw_inq_dimlen(fname, fid, ncw_dim_id(fid, "xi_u"),    &xi_u);
  ncw_inq_dimlen(fname, fid, ncw_dim_id(fid, "xi_v"),    &xi_v);
  ncw_inq_dimlen(fname, fid, ncw_dim_id(fid, "eta_rho"), &eta_rho);
  ncw_inq_dimlen(fname, fid, ncw_dim_id(fid, "eta_psi"), &eta_psi);
  ncw_inq_dimlen(fname, fid, ncw_dim_id(fid, "eta_u"),   &eta_u);
  ncw_inq_dimlen(fname, fid, ncw_dim_id(fid, "eta_v"),   &eta_v);
  nce1 = xi_rho;
  nce2 = eta_rho;

  /* RHO */
  count[0] = eta_rho;
  count[1] = xi_rho;
  mask_rho = d_alloc_2d(count[1], count[0]);
  ncw_get_vara_double(fname, fid, ncw_var_id(fid, "mask_rho"), start, count, mask_rho[0]);

  /* PSI */
  count[0] = eta_psi;
  count[1] = xi_psi;
  mask_psi = d_alloc_2d(count[1], count[0]);
  ncw_get_vara_double(fname, fid, ncw_var_id(fid, "mask_psi"), start, count, mask_psi[0]);

  /* U */
  count[0] = eta_u;
  count[1] = xi_u;
  mask_u = d_alloc_2d(count[1], count[0]);
  ncw_get_vara_double(fname, fid, ncw_var_id(fid, "mask_u"), start, count, mask_u[0]);

  /* V */
  count[0] = eta_v;
  count[1] = xi_v;
  mask_v = d_alloc_2d(count[1], count[0]);
  ncw_get_vara_double(fname, fid, ncw_var_id(fid, "mask_v"), start, count, mask_v[0]);

  /* Fill in the roms struct */
  if (romsgrid != NULL) {
    memcpy(romsgrid->mask_rho[0], mask_rho[0], eta_rho*xi_rho*sizeof(double));
    memcpy(romsgrid->mask_psi[0], mask_psi[0], eta_psi*xi_psi*sizeof(double));
    memcpy(romsgrid->mask_u[0], mask_u[0], eta_u*xi_u*sizeof(double));
    memcpy(romsgrid->mask_v[0], mask_v[0], eta_v*xi_v*sizeof(double));
  } else {
    /* Apply mask */
    for (n=0, j=0; j<nce2; j++)
      for (i=0; i<nce1; i++)
	if (mask_rho[j][i])
	  params->bathy[n++] = -LANDCELL;
  }

  hd_warn("auto_params: Masking applied from ROMS file %s\n", fname);

  /* Cleanup */
  d_free_2d(mask_rho);
  d_free_2d(mask_psi);
  d_free_2d(mask_u);
  d_free_2d(mask_v);
  ncw_close(fname, fid);

}


/*-------------------------------------------------------------------*/
/* Reads bathymetry from a ROMS grid_spec.nc file                    */
/*-------------------------------------------------------------------*/
double *read_bathy_from_roms(char *fname, int nx, int ny, romsgrid_t *romsgrid)
{
  int fid;
  int i, j, n;
  size_t start[4];
  size_t count[4];
  size_t xi_rho, eta_rho;
  size_t bath, one;
  double *bathy = NULL;
  double **h, ***hraw, depthmin, depthmax;
  int nce1, nce2;

  memset(start, 0, 4*sizeof(size_t));
  memset(count, 0, 4*sizeof(size_t));

  /*-----------------------------------------------------------------*/
  /* Open the dump file for reading                                  */
  ncw_open(fname, NC_NOWRITE, &fid);

  /*
   * Get ROMS dimensions
   */
  ncw_inq_dimlen(fname, fid, ncw_dim_id(fid, "bath"), &bath);
  ncw_inq_dimlen(fname, fid, ncw_dim_id(fid, "one"),  &one);
  if (one != 1 && bath !=1)
    hd_quit("Non-scalar one and bath ROMS ids found\n");

  ncw_inq_dimlen(fname, fid, ncw_dim_id(fid, "xi_rho"),   &xi_rho);
  ncw_inq_dimlen(fname, fid, ncw_dim_id(fid, "eta_rho"),  &eta_rho);
  nce1 = xi_rho;
  nce2 = eta_rho;
  if (nce1 != nx || nce2 != ny)
    hd_quit("Inconsistent grid, in bathy, sizes EMS (%dx%d) : ROMS (%dx%d)\n",
	    nx, ny, nce1, nce2);

  /*
   * Read bathy variables
   */
  ncw_get_var_double(fname, fid, ncw_var_id(fid, "depthmin"), &depthmin);
  ncw_get_var_double(fname, fid, ncw_var_id(fid, "depthmax"), &depthmax);
  count[0] = eta_rho;
  count[1] = xi_rho;
  h = d_alloc_2d(count[1], count[0]);
  ncw_get_vara_double(fname, fid, ncw_var_id(fid, "h"), start, count, h[0]);
  
  // xxx need to sort out the order
  count[0] = bath;
  count[1] = eta_rho;
  count[2] = xi_rho;
  hraw = d_alloc_3d(count[2], count[1], count[0]);
  ncw_get_vara_double(fname, fid, ncw_var_id(fid, "hraw"), start, count, hraw[0][0]);

  /* Fill in the roms struct */
  if (romsgrid != NULL) {
    romsgrid->bmin = depthmin;
    romsgrid->bmax = depthmax;
    memcpy(romsgrid->h[0],       h[0],       eta_rho*xi_rho*sizeof(double));
    memcpy(romsgrid->hraw[0][0], hraw[0][0], eta_rho*xi_rho*sizeof(double));
  } else {
    /* Assign bathymetry */
    bathy = d_alloc_1d(nce1 * nce2);
    for (n = 0, j = 0; j < nce2; j++)
      for (i = 0; i < nce1; i++)
	bathy[n++] = h[j][i];
  }

  hd_warn("auto_params: Bathymetry from ROMS file %s\n", fname);

  /* Cleanup */
  d_free_2d(h);
  d_free_3d(hraw);
  ncw_close(fname, fid);

  return(bathy);
}

/*-------------------------------------------------------------------*/
/* Converts ROMS grid_spec.nc file to EMS metrics                    */
/*-------------------------------------------------------------------*/
void read_roms_grid(char *fname, int nx, int ny, double **x, double **y, romsgrid_t *romsgrid)
{
  int fid;
  int i, j;
  size_t start[4];
  size_t count[4];
  size_t xi_rho, eta_rho;
  size_t xi_psi, eta_psi;
  size_t xi_u, eta_u;
  size_t xi_v, eta_v;
  double **lat_rho, **lon_rho;
  double **lat_psi, **lon_psi;
  double **lat_u, **lon_u;
  double **lat_v, **lon_v;
  int nce1, nce2;

  /* Note : We're assuming this is always true */
  int inface = 1;

  memset(start, 0, 4*sizeof(size_t));
  memset(count, 0, 4*sizeof(size_t));

  /*-----------------------------------------------------------------*/
  /* Open the dump file for reading                                  */
  ncw_open(fname, NC_NOWRITE, &fid);

  /*
   * Get ROMS dimensions
   */
  ncw_inq_dimlen(fname, fid, ncw_dim_id(fid, "xi_rho"),  &xi_rho);
  ncw_inq_dimlen(fname, fid, ncw_dim_id(fid, "xi_psi"),  &xi_psi);
  ncw_inq_dimlen(fname, fid, ncw_dim_id(fid, "xi_u"),    &xi_u);
  ncw_inq_dimlen(fname, fid, ncw_dim_id(fid, "xi_v"),    &xi_v);
  ncw_inq_dimlen(fname, fid, ncw_dim_id(fid, "eta_rho"), &eta_rho);
  ncw_inq_dimlen(fname, fid, ncw_dim_id(fid, "eta_psi"), &eta_psi);
  ncw_inq_dimlen(fname, fid, ncw_dim_id(fid, "eta_u"),   &eta_u);
  ncw_inq_dimlen(fname, fid, ncw_dim_id(fid, "eta_v"),   &eta_v);
  nce1 = xi_rho;
  nce2 = eta_rho;
  if (nce1 != nx || nce2 != ny)
    hd_quit("Inconsistent grid sizes EMS (%dx%d) : ROMS (%dx%d)\n",
	    nx, ny, nce1, nce2);

  /* Read cell centre locations */
  count[0] = eta_rho;
  count[1] = xi_rho;
  lat_rho = d_alloc_2d(count[1], count[0]);
  lon_rho = d_alloc_2d(count[1], count[0]);
  ncw_get_vara_double(fname, fid, ncw_var_id(fid, "lat_rho"), start, count, lat_rho[0]);
  ncw_get_vara_double(fname, fid, ncw_var_id(fid, "lon_rho"), start, count, lon_rho[0]);
  for (j=0; j<eta_rho; j++) {
    for (i=0; i<xi_rho; i++) {
      x[2*j+1][2*i+1] = lon_rho[j][i];
      y[2*j+1][2*i+1] = lat_rho[j][i];
    }
  }
  /* Read the cell corners */
  count[0] = eta_psi;
  count[1] = xi_psi;
  lat_psi = d_alloc_2d(count[1], count[0]);
  lon_psi = d_alloc_2d(count[1], count[0]);
  ncw_get_vara_double(fname, fid, ncw_var_id(fid, "lat_psi"), start, count, lat_psi[0]);
  ncw_get_vara_double(fname, fid, ncw_var_id(fid, "lon_psi"), start, count, lon_psi[0]);
  for (j=0; j<eta_psi; j++) {
    for (i=0; i<xi_psi; i++) {
      x[2*j+2][2*i+2] = lon_psi[j][i];
      y[2*j+2][2*i+2] = lat_psi[j][i];
    }
  }
  /* Read the u locations */
  count[0] = eta_u;
  count[1] = xi_u;
  lat_u = d_alloc_2d(count[1], count[0]);
  lon_u = d_alloc_2d(count[1], count[0]);
  ncw_get_vara_double(fname, fid, ncw_var_id(fid, "lat_u"), start, count, lat_u[0]);
  ncw_get_vara_double(fname, fid, ncw_var_id(fid, "lon_u"), start, count, lon_u[0]);
  for (j=0; j<eta_u; j++) {
    for (i=0; i<xi_u; i++) {
      x[2*j+1][2*i+2] = lon_u[j][i];
      y[2*j+1][2*i+2] = lat_u[j][i];
    }
  }
  /* Read the v locations */
  count[0] = eta_v;
  count[1] = xi_v;
  lat_v = d_alloc_2d(count[1], count[0]);
  lon_v = d_alloc_2d(count[1], count[0]);
  ncw_get_vara_double(fname, fid, ncw_var_id(fid, "lat_v"), start, count, lat_v[0]);
  ncw_get_vara_double(fname, fid, ncw_var_id(fid, "lon_v"), start, count, lon_v[0]);
  for (j=0; j<eta_v; j++) {
    for (i=0; i<xi_v; i++) {
      x[2*j+2][2*i+1] = lon_v[j][i];
      y[2*j+2][2*i+1] = lat_v[j][i];
    }
  }
  /*
   * Extrapolate the outer border
   */
  // faces
  // left edge
  if (eta_rho != eta_u)
    hd_quit("Inconsistent ROMS grid dimensions expected eta_rho == eta_u, got eta_rho=%d and eta_u=%d\n", eta_rho, eta_u);
  for (j=0; j<eta_rho; j++) {
    x[2*j+1][0] = lon_rho[j][0] - (lon_u[j][0] - lon_rho[j][0]);
    y[2*j+1][0] = lat_rho[j][0] - (lat_u[j][0] - lat_rho[j][0]);
  }

  // back edge
  if (xi_rho != xi_v)
    hd_quit("Inconsistent ROMS grid dimensions expected xi_rho == xi_v, got xi_rho=%d and xi_v=%d\n", xi_rho, xi_v);
  for (i=0; i<xi_rho; i++) {
    x[0][2*i+1] = lon_rho[0][i] - (lon_v[0][i] - lon_rho[0][i]);
    y[0][2*i+1] = lat_rho[0][i] - (lat_v[0][i] - lat_rho[0][i]);
  }
  // right edge
  for (j=0; j<eta_rho; j++) {
    x[2*j+1][2*nce1] = lon_rho[j][xi_rho-1] + (lon_rho[j][xi_rho-1] - lon_u[j][xi_u-1]);
    y[2*j+1][2*nce1] = lat_rho[j][xi_rho-1] + (lat_rho[j][xi_rho-1] - lat_u[j][xi_u-1]);
  }
  // top edge
  for (i=0; i<xi_rho; i++) {
    x[2*nce2][2*i+1] = lon_rho[eta_rho-1][i] + (lon_rho[eta_rho-1][i] - lon_v[eta_v-1][i]);
    y[2*nce2][2*i+1] = lat_rho[eta_rho-1][i] + (lat_rho[eta_rho-1][i] - lat_v[eta_v-1][i]);
  }

  // corners
  // left edge
  if (eta_psi != eta_v)
    hd_quit("Inconsistent ROMS grid dimensions expected eta_psi == eta_v, got eta_psi=%d and eta_v=%d\n", eta_psi, eta_v);
  for (j=0; j<eta_psi; j++) {
    x[2*j+2][0] = lon_v[j][0] - (lon_psi[j][0] - lon_v[j][0]);
    y[2*j+2][0] = lat_v[j][0] - (lat_psi[j][0] - lat_v[j][0]);
  }
  // back edge
  if (xi_psi != xi_u)
    hd_quit("Inconsistent ROMS grid dimensions expected xi_psi == xi_u, got xi_psi=%d and xi_u=%d\n", xi_psi, xi_u);
  for (i=0; i<xi_psi; i++) {
    x[0][2*i+2] = lon_u[0][i] - (lon_psi[0][i] - lon_u[0][i]);
    y[0][2*i+2] = lat_u[0][i] - (lat_psi[0][i] - lat_u[0][i]);
  }
  // right edge
  for (j=0; j<eta_psi; j++) {
    x[2*j+2][2*nce1] = lon_v[j][xi_v-1] + (lon_v[j][xi_v-1] - lon_psi[j][xi_psi-1]);
    y[2*j+2][2*nce1] = lat_v[j][xi_v-1] + (lat_v[j][xi_v-1] - lat_psi[j][xi_psi-1]);
  }
  // top edge
  for (i=0; i<xi_psi; i++) {
    x[2*nce2][2*i+2] = lon_u[eta_u-1][i] + (lon_u[eta_u-1][i] - lon_psi[eta_psi-1][i]);
    y[2*nce2][2*i+2] = lat_u[eta_u-1][i] + (lat_u[eta_u-1][i] - lat_psi[eta_psi-1][i]);
  }
  // The remaining 4 outer corners
  // bottom left
  x[0][0] = x[0][1] - (x[0][2] - x[0][1]);
  y[0][0] = y[1][0] - (y[2][0] - y[1][0]);
  // bottom right
  x[0][2*nce1] = x[0][2*nce1-1] + (x[0][2*nce1-1] - x[0][2*nce1-2]);
  y[0][2*nce1] = y[1][2*nce1]   - (y[2][2*nce1] - y[1][2*nce1]);
  // top right
  x[2*nce2][2*nce1] = x[2*nce2][2*nce1-1] + (x[2*nce2][2*nce1-1] - x[2*nce2][2*nce1-2]);
  y[2*nce2][2*nce1] = y[2*nce2-1][2*nce1] + (y[2*nce2-1][2*nce1] - y[2*nce2-2][2*nce1]);
  // top left
  x[2*nce2][0] = x[2*nce2][1]   - (x[2*nce2][2] - x[2*nce2][1]);
  y[2*nce2][0] = y[2*nce2-1][0] + (y[2*nce2-1][0] - y[2*nce2-2][0]);

  /* Fill in the romsgrid struct */
  if (romsgrid != NULL) {
    /* Dimensions */
    romsgrid->nfe1 = (inface) ? params->nce1 : params->nce1 + 1;
    romsgrid->nfe2 = (inface) ? params->nce2 : params->nce2 + 1;
    romsgrid->xi_rho  = xi_rho;
    romsgrid->eta_rho = eta_rho;
    romsgrid->xi_psi  = xi_psi;
    romsgrid->eta_psi = eta_psi;
    romsgrid->xi_u  = xi_u;
    romsgrid->eta_u = eta_u;
    romsgrid->xi_v  = xi_v;
    romsgrid->eta_v = eta_v;
    /* Allocate mem */
    alloc_grid_mem_roms(romsgrid, nce1, nce2);
    /* lon and lat variables */
    memcpy(romsgrid->lon_rho[0], lon_rho[0], eta_rho*xi_rho*sizeof(double));
    memcpy(romsgrid->lat_rho[0], lat_rho[0], eta_rho*xi_rho*sizeof(double));
    memcpy(romsgrid->lon_psi[0], lon_psi[0], eta_psi*xi_psi*sizeof(double));
    memcpy(romsgrid->lat_psi[0], lat_psi[0], eta_psi*xi_psi*sizeof(double));
    memcpy(romsgrid->lon_u[0], lon_u[0], eta_u*xi_u*sizeof(double));
    memcpy(romsgrid->lat_u[0], lat_u[0], eta_u*xi_u*sizeof(double));
    memcpy(romsgrid->lon_v[0], lon_v[0], eta_v*xi_v*sizeof(double));
    memcpy(romsgrid->lat_v[0], lat_v[0], eta_v*xi_v*sizeof(double));
    /* metrics */
    count[0] = eta_rho; count[1] = xi_rho;
    ncw_get_vara_double(fname, fid, ncw_var_id(fid, "f"),  start, count, romsgrid->f[0]);
    ncw_get_vara_double(fname, fid, ncw_var_id(fid, "pm"), start, count, romsgrid->pm[0]);
    ncw_get_vara_double(fname, fid, ncw_var_id(fid, "pn"), start, count, romsgrid->pn[0]);
    ncw_get_vara_double(fname, fid, ncw_var_id(fid, "dndx"), start, count, romsgrid->dndx[0]);
    ncw_get_vara_double(fname, fid, ncw_var_id(fid, "dmde"), start, count, romsgrid->dmde[0]);
    ncw_get_vara_double(fname, fid, ncw_var_id(fid, "angle"), start, count, romsgrid->theta[0]);
    ncw_get_var_double(fname, fid, ncw_var_id(fid, "xl"), &romsgrid->xl);
    ncw_get_var_double(fname, fid, ncw_var_id(fid, "el"), &romsgrid->el);
  }

  hd_warn("auto_params: Reading ROMS grid from %s\n", fname);

  /* Cleanup */
  d_free_2d(lon_rho); d_free_2d(lat_rho);
  d_free_2d(lon_psi); d_free_2d(lat_psi);
  d_free_2d(lon_u); d_free_2d(lat_u);
  d_free_2d(lon_v); d_free_2d(lat_v);
  ncw_close(fname, fid);

}

/*------------------------------------------------------------------*/
/* Converts EMS grid to ROMS grid                                   */
/*------------------------------------------------------------------*/
romsgrid_t *convert_roms_grid(parameters_t *params, dump_data_t *dumpdata)
{
  romsgrid_t *romsgrid;         /* Grid data                        */
  double maxb = 0.0;
  int n, i, j, ip, jp, k;       /* Counters                         */
  int nce1, nce2;               /* x and y grid size                */
  int si;
  /* This flag needs to be consistent with df_roms */
  int inface = 1;
  si = (inface) ? 1 : 0;

  if (params->doroms == 0) return(NULL);

  /*-----------------------------------------------------------------*/
  /* Set up the grid structure */
  romsgrid = romsgrid_alloc();
  strcpy(romsgrid->romsfile, params->romsfile);
  romsgrid->doroms = params->doroms;

  nce1 = romsgrid->nce1 = params->nce1;
  nce2 = romsgrid->nce2 = params->nce2;

  /*-----------------------------------------------------------------*/
  /* Read grid attributes */
  if (params->roms_grid_opts & ROMS_GRID_2D) {
    /* Read from ROMS file */
    read_roms_grid(params->roms_grid, nce1, nce2, params->x, params->y, romsgrid);
  } else {
    romsgrid->nfe1 = (inface) ? params->nce1 : params->nce1 + 1;
    romsgrid->nfe2 = (inface) ? params->nce2 : params->nce2 + 1;
    romsgrid->xi_rho = romsgrid->nce1;
    romsgrid->eta_rho = romsgrid->nce2;
    romsgrid->xi_psi = (inface) ? romsgrid->nce1 - 1 : romsgrid->nce1;
    romsgrid->eta_psi = (inface) ? romsgrid->nce2 - 1 : romsgrid->nce2;
    romsgrid->xi_u = (inface) ? romsgrid->nce1 - 1 : romsgrid->nce1;
    romsgrid->eta_u = romsgrid->nce2;
    romsgrid->xi_v = romsgrid->nce1;
    romsgrid->eta_v = (inface) ? romsgrid->nce2 - 1 : romsgrid->nce2;
    romsgrid->vertex = 4;
    romsgrid->nobc = params->nobc;

    alloc_grid_mem_roms(romsgrid, nce1, nce2);
    if (params->doroms & RDGRID) read_grid(params);
    
    if (DEBUG("init_m")) {
      dlog("init_m", "SHOC grid = (%d x %d)\n",romsgrid->nce1, romsgrid->nce2);
      dlog("init_m", "ROMS grid = (%d x %d)\n",romsgrid->xi_rho, romsgrid->eta_rho);
    }
    
    /*-----------------------------------------------------------------*/
    /* Get the grid metrics                                            */
    for (j = 0; j <= 2*romsgrid->nce2; j++)
      for (i = 0; i <= 2*romsgrid->nce1; i++) {
	romsgrid->x[j][i] = params->x[j][i];
	romsgrid->y[j][i] = params->y[j][i];
	romsgrid->h1[j][i] = params->h1[j][i];
	romsgrid->h2[j][i] = params->h2[j][i];
	romsgrid->a1[j][i] = params->a1[j][i];
	romsgrid->a2[j][i] = params->a2[j][i];
      }
    
    /* Insert nearest neighbours for NaNs                              */
    for (j = 0; j <= 2*nce2; j++)
      for (i = 0; i <= 2*nce1; i++) {
	if(isnan(romsgrid->x[j][i])) {
	  find_non_nan(romsgrid->x, 2*nce1, 2*nce2, i, j, &ip, &jp);
	  romsgrid->x[j][i] = romsgrid->x[jp][ip];
	}
	if(isnan(romsgrid->y[j][i])) {
	  find_non_nan(romsgrid->y, 2*nce1, 2*nce2, i, j, &ip, &jp);
	  romsgrid->y[j][i] = romsgrid->y[jp][ip];
	}
	if(isnan(romsgrid->h1[j][i])) {
	  find_non_nan(romsgrid->h1, 2*nce1, 2*nce2, i, j, &ip, &jp);
	  romsgrid->h1[j][i] = romsgrid->h1[jp][ip];
	}
	if(isnan(romsgrid->h2[j][i])) {
	  find_non_nan(romsgrid->h2, 2*nce1, 2*nce2, i, j, &ip, &jp);
	  romsgrid->h2[j][i] = romsgrid->h2[jp][ip];
	}
	if(isnan(romsgrid->a1[j][i])) {
	  find_non_nan(romsgrid->a1, 2*nce1, 2*nce2, i, j, &ip, &jp);
	  romsgrid->a1[j][i] = romsgrid->a1[jp][ip];
	}
	if(isnan(romsgrid->a2[j][i])) {
	  find_non_nan(romsgrid->a2, 2*nce1, 2*nce2, i, j, &ip, &jp);
	  romsgrid->a2[j][i] = romsgrid->a2[jp][ip];
	}
      }
    
    /* Insert dummy values for NaNs                                    */
    /*
      for (j = 0; j <= 2*romsgrid->nce2; j++)
      for (i = 0; i <= 2*romsgrid->nce1; i++) {
      if (isnan(romsgrid->x[j][i]) && i > 1)
      romsgrid->x[j][i] = 2.0*romsgrid->x[j][i-1]-romsgrid->x[j][i-2]; 
      if (isnan(romsgrid->h1[j][i]) && i > 1)
      romsgrid->h1[j][i] = 2.0*romsgrid->h1[j][i-1]-romsgrid->h1[j][i-2]; 
      if (isnan(romsgrid->a1[j][i]) && i > 1)
      romsgrid->a1[j][i] = 2.0*romsgrid->a1[j][i-1]-romsgrid->a1[j][i-2];
      if (isnan(romsgrid->x[j][i]) && j > 1)
      romsgrid->x[j][i] = 2.0*romsgrid->x[j-1][i]-romsgrid->x[j-2][i];
      if (isnan(romsgrid->h1[j][i]) && j > 1)
      romsgrid->h1[j][i] = 2.0*romsgrid->h1[j-1][i]-romsgrid->h1[j-2][i];
      if (isnan(romsgrid->a1[j][i]) && j > 1)
      romsgrid->a1[j][i] = 2.0*romsgrid->a1[j-1][i]-romsgrid->a1[j-2][i];
      
      if (isnan(romsgrid->y[j][i]) && j > 1)
      romsgrid->y[j][i] = 2.0*romsgrid->y[j-1][i]-romsgrid->y[j-2][i];
      
      if (isnan(romsgrid->h2[j][i]) && j > 1)
      romsgrid->h2[j][i] = 2.0*romsgrid->h2[j-1][i]-romsgrid->h2[j-2][i];
      if (isnan(romsgrid->a2[j][i]) && j > 1)
      romsgrid->a2[j][i] = 2.0*romsgrid->a2[j-1][i]-romsgrid->a2[j-2][i];
      
      if (isnan(romsgrid->y[j][i]) && i > 1)
      romsgrid->y[j][i] = 2.0*romsgrid->y[j][i-1]-romsgrid->y[j][i-2];
      
      if (isnan(romsgrid->h2[j][i]) && i > 1)
      romsgrid->h2[j][i] = 2.0*romsgrid->h2[j][i-1]-romsgrid->h2[j][i-2];
      if (isnan(romsgrid->a2[j][i]) && i > 1)
      romsgrid->a2[j][i] = 2.0*romsgrid->a2[j][i-1]-romsgrid->a2[j][i-2];
      }
    */
    
    /*-----------------------------------------------------------------*/
    /* Get the grid positions                                          */
    /* Grid centres                                                    */
    for (j = 0; j < romsgrid->nce2; j++) {
      for (i = 0; i < romsgrid->nce1; i++) {
	romsgrid->lon_rho[j][i] = romsgrid->x[j * 2 + 1][i * 2 + 1];
	romsgrid->lat_rho[j][i] = romsgrid->y[j * 2 + 1][i * 2 + 1];
      }
    }
    /* Grid corners                                                    */
    for (j = si; j < romsgrid->nfe2; j++) {
      for (i = si; i < romsgrid->nfe1; i++) {
	romsgrid->lon_psi[j-si][i-si] = romsgrid->x[j * 2][i * 2];
	romsgrid->lat_psi[j-si][i-si] = romsgrid->y[j * 2][i * 2];
      }
    }
    /* u1 faces                                                        */
    for (j = 0; j < romsgrid->nce2; j++) {
      for (i = si; i < romsgrid->nfe1; i++) {
	romsgrid->lon_u[j][i-si] = romsgrid->x[j * 2 + 1][i * 2];
	romsgrid->lat_u[j][i-si] = romsgrid->y[j * 2 + 1][i * 2];
      }
    }
    /* u2 faces                                                        */
    for (j = si; j < romsgrid->nfe2; j++) {
      for (i = 0; i < romsgrid->nce1; i++) {
	romsgrid->lon_v[j-si][i] = romsgrid->x[j * 2][i * 2 + 1];
	romsgrid->lat_v[j-si][i] = romsgrid->y[j * 2][i * 2 + 1];
      }
    }
    
    /*-----------------------------------------------------------------*/
    /* Get the angles at the cell centre                               */
    for (j = 0; j < romsgrid->nce2; j++) {
      for (i = 0; i < romsgrid->nce1; i++) {
	romsgrid->thetay[j][i] = params->a2[j * 2 + 1][i * 2 + 1] - PI / 2.0;
	romsgrid->thetax[j][i] = params->a1[j * 2 + 1][i * 2 + 1];
	romsgrid->theta[j][i] = 0.5 * (romsgrid->thetax[j][i] + romsgrid->thetay[j][i]);
      }
    }
    
    /*-----------------------------------------------------------------*/
    /* Get the grid metrics                                            */
    for (j = 0; j < romsgrid->nce2; j++)
      for (i = 0; i < romsgrid->nce1; i++) {
	romsgrid->dndx[j][i] = 0.0;
	romsgrid->dmde[j][i] = 0.0;
	romsgrid->pm[j][i] = 1.0 / dumpdata->h1acell[j][i]; 
	romsgrid->pn[j][i] = 1.0 / dumpdata->h2acell[j][i]; 
      }
    for (j = 0; j < romsgrid->nce2; j++)
      for (i = 0; i < romsgrid->nce1; i++) {
	int ip = (i < romsgrid->nce1 - 1) ? i+1 : i;
	int jp = (j < romsgrid->nce2 - 1) ? j+1 : j;
	int im = (i > 0) ? i-1 : i;
	int jm = (j > 0) ? j-1 : j;
	romsgrid->dndx[j][i] = 0.5 * (dumpdata->h2acell[j][ip] - dumpdata->h2acell[j][im]);
	romsgrid->dmde[j][i] = 0.5 * (dumpdata->h1acell[j][ip] - dumpdata->h1acell[j][im]);
      }
    
    for (j = 0; j < romsgrid->nce2; j++)
      for (i = 0; i < romsgrid->nce1; i++)
	romsgrid->f[j][i] = dumpdata->coriolis[j][i];

    romsgrid->xl = romsgrid->el = 0.0;
    for (i = 0; i < romsgrid->nce1; i++)
      romsgrid->xl += dumpdata->h1acell[romsgrid->eta_rho/2][i];
    for (j = 0; j < romsgrid->nce2; j++)
      romsgrid->el += dumpdata->h2acell[j][romsgrid->xi_rho/2];
  } /* end if (ROMS_GRID_2D) */

  /*-----------------------------------------------------------------*/
  /* Get the obc info. Note; ROMS origin is (1,1), SHOC is (0,0), so */
  /* add 1 to the SHOC OBC indices.                                  */
  for (n = 0; n < romsgrid->nobc; n++) {
    i = j = 1;
    if (geom->open[n]->ocodex & R_EDGE) i = 0;
    if (geom->open[n]->ocodey & F_EDGE) j = 0;
    romsgrid->obcis[n] = min(params->open[n]->iloc[0] + i, romsgrid->xi_psi);
    romsgrid->obcjs[n] = min(params->open[n]->jloc[0] + j, romsgrid->eta_psi);
    romsgrid->obcie[n] = min(params->open[n]->iloc[params->open[n]->npts - 1] + i,
			    romsgrid->xi_psi);
    romsgrid->obcje[n] = min(params->open[n]->jloc[params->open[n]->npts - 1] + j,
			    romsgrid->eta_psi);
  }

  /*-----------------------------------------------------------------*/
  /* Read the bathymetry and set the mask at RHO cells               */
  if (params->roms_grid_opts & ROMS_GRID_BATHY) {
    /* Read from ROMS file */
    read_bathy_from_roms(params->roms_grid, nce1, nce2, dumpdata->romsgrid);
  } else {
    romsgrid->bmin = params->bmin;
    romsgrid->bmax = params->bmax;
    for (j = 0; j < romsgrid->nce2; j++)
      for (i = 0; i < romsgrid->nce1; i++) {
	double depth = dumpdata->botz[j][i];
	if (depth == LANDCELL || depth == NOTVALID || isnan(depth)) {
	  romsgrid->hraw[0][j][i] = romsgrid->h[j][i] = NaN;
	  romsgrid->num_levels[j][i] = 0.0;
	}
	else {
	  romsgrid->hraw[0][j][i] = romsgrid->h[j][i] = -depth;
	  romsgrid->num_levels[j][i] = -1.0;
	  for (k = 0; k < params->nz; k++) {
	    if (depth >= params->layers[k] && depth < params->layers[k + 1])
	      romsgrid->num_levels[j][i] = (double)(params->nz - k);
	  }
	  if (depth > maxb) maxb = depth;
	}

	if (romsgrid->num_levels[j][i] < 0)
	  hd_quit("Bottom not found at RHO-cell (%d %d) : %5.1f m\n",i, j, depth);
      }
  } /* end if (ROMS_GRID_BATHY) */
  
  /*-----------------------------------------------------------------*/
  /* Set the mask at PSI cells (grid corners)                        */
  if (params->roms_grid_opts & ROMS_GRID_MASK) {
    /* Read from ROMS file */
    mask_bathy_from_roms(params->roms_grid, params, dumpdata->romsgrid);
  } else {
    /* Mask at cell centres */
    for (j = 0; j < romsgrid->nce2; j++) {
      for (i = 0; i < romsgrid->nce1; i++) {
	double depth = dumpdata->botz[j][i];
	if (depth == LANDCELL || depth == NOTVALID || isnan(depth))
	  romsgrid->mask_rho[j][i] = 0.0;
	else
	  romsgrid->mask_rho[j][i] = 1.0;
      }
    }

    for (j = si; j < romsgrid->nfe2; j++)
      for (i = si; i < romsgrid->nfe1; i++)
	romsgrid->mask_psi[j-si][i-si] = 0;
    /*
     * % Land/Sea mask on PSI-points.
     * pmask(1:L,1:M)=rmask(1:L,1:M ).*rmask(2:Lp,1:M ).* ...
     *                rmask(1:L,2:Mp).*rmask(2:Lp,2:Mp);
     */
    for (j = si; j < romsgrid->nce2; j++)
      for (i = si; i < romsgrid->nce1; i++) {
	int r1 = romsgrid->mask_rho[j-si][i-si];
	int r2 = romsgrid->mask_rho[j-si][i-si+1];
	int r3 = romsgrid->mask_rho[j-si+1][i-si];
	int r4 = romsgrid->mask_rho[j-si+1][i-si+1];
	romsgrid->mask_psi[j-si][i-si] = r1 & r2 & r3 & r4;
      }
    /*
      for (j = si; j < romsgrid->nfe2; j++)
      for (i = si; i < romsgrid->nfe1; i++) {
      int ip = (i < romsgrid->nfe1 - 1) ? i+1 : i;
      int jp = (j < romsgrid->nfe2 - 1) ? j+1 : j;
      if (!isnan(romsgrid->h[j][i])) {
      romsgrid->mask_psi[j-si][i-si] = 1;
      romsgrid->mask_psi[j-si][ip-si] = 1;
      romsgrid->mask_psi[jp-si][i-si] = 1;
      romsgrid->mask_psi[jp-si][ip-si] = 1;
      }
      }
    */
    
    /*-----------------------------------------------------------------*/
    /* Set the mask at U cells                                        */
    for (j = 0; j < romsgrid->nce2; j++)
      for (i = si; i < romsgrid->nfe1; i++)
	romsgrid->mask_u[j][i-si] = 0;
    /*
     * Key-off the rho mask as per ROMS Matlab code
     * [Lp,Mp]=size(rmask);
     * L=Lp-1;
     * M=Mp-1;
     * % Land/Sea mask on U-points.
     * umask(1:L,1:Mp)=rmask(2:Lp,1:Mp).*rmask(1:L,1:Mp);
     */
    for (j = 0; j < romsgrid->nce2; j++)
      for (i = si; i < romsgrid->nce1; i++) {
	int r1 = romsgrid->mask_rho[j][i-si+1];
	int r2 = romsgrid->mask_rho[j][i-si];
	romsgrid->mask_u[j][i-si] = r1 & r2;
      }
    /*
      for (j = 0; j < romsgrid->nce2; j++)
      for (i = si; i < romsgrid->nfe1; i++) {
      int ip = (i < romsgrid->nfe1 - 1) ? i+1 : i;
      if (!isnan(romsgrid->h[j][i])) {
      romsgrid->mask_u[j][i-si] = 1;
      romsgrid->mask_u[j][ip-si] = 1;
      }
      }
    */
    
    /*-----------------------------------------------------------------*/
    /* Set the mask at V cells                                        */
    for (j = si; j < romsgrid->nfe2; j++)
      for (i = 0; i < romsgrid->nce1; i++)
	romsgrid->mask_v[j-si][i] = 0;
    /*
     * % Land/Sea mask on V-points.
     * umask(1:L,1:Mp)=rmask(2:Lp,1:Mp).*rmask(1:L,1:Mp);
     */
    for (j = si; j < romsgrid->nfe2; j++)
      for (i = 0; i < romsgrid->nce1; i++) {
	int r1 = romsgrid->mask_rho[j-si+1][i];
	int r2 = romsgrid->mask_rho[j-si][i];
	romsgrid->mask_v[j-si][i] = r1 & r2;
      }
    /*
      for (j = si; j < romsgrid->nfe2; j++)
      for (i = 0; i < romsgrid->nce1; i++) {
      int jp = (j < romsgrid->nfe2 - 1) ? j+1 : j;
      if (!isnan(romsgrid->h[j][i])) {
      romsgrid->mask_v[j-si][i] = 1;
      romsgrid->mask_v[jp-si][i] = 1;
      }
      }
    */
  } /* end if (ROMS_GRID_MASK) */

  /*-----------------------------------------------------------------*/
  /* Set up the vertical grid                                        */
  romsgrid->nz = params->nz;
  if (params->roms_z2s <= 1.)
    romsgrid->nsigma = romsgrid->nz * params->roms_z2s;
  else if (params->roms_z2s > 1)
    romsgrid->nsigma = (int)params->roms_z2s;
  else
    romsgrid->nsigma = romsgrid->nz;

  /* Allocate memory in the vertical */
  alloc_grid_mem_vert_roms(romsgrid);
  
  for (k = 0; k < romsgrid->nz; k++)
    romsgrid->layers[k] = params->layers[k];
  romsgrid->layers[romsgrid->nz] = 0.0;
  
  for (k = 0; k < romsgrid->nz; k++) {
    i = romsgrid->nz - 1 - k;
    romsgrid->zb[i] = fabs(romsgrid->layers[k]);
    romsgrid->zt[i] = fabs(0.5 * (romsgrid->layers[k] + romsgrid->layers[k + 1]));
  }
  /* Convert to sigma                                               */
  romsgrid->layers[romsgrid->nz-1] = 0.0;
  romsgrid->layers[0] = 1.0;
  for (k = 1; k < romsgrid->nz - 1; k++) {
    romsgrid->layers[k] /= maxb;
  }
  
  /* Create sigma levels and stretching function */
  roms_compute_vert(romsgrid->nsigma,   romsgrid->sc_rho, romsgrid->Cs_rho, 0);
  roms_compute_vert(romsgrid->nsigma+1, romsgrid->sc_w,   romsgrid->Cs_w,   1);
  
  return(romsgrid);
}

/* END convert_roms_grid()                                          */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Writes EMS grid to ROMS grid conversion to netCDF file            */
/*------------------------------------------------------------------*/
void write_roms_grid(dump_data_t *dumpdata)
{
  romsgrid_t *romsgrid = dumpdata->romsgrid;
  FILE *op = NULL;
  int i, j;                     /* Counters                         */
  int mode;                     /* Output mode                      */
  int cfid;
  size_t start[4];
  size_t count[4];
  short v;
  char bv;

  if (romsgrid == NULL) return;
  mode = romsgrid->doroms;

  /* Bail out if not writing to netcdf */
  if (!(mode & CDF)) return;

  /*----------------------------------------------------------------*/
  /* Convert metrics to ROMS format and write to file                */
  if (!find_name_marvl(dumpdata, romsgrid->ofile, "roms_spec", NULL))
    sprintf(romsgrid->ofile, "%s_roms.nc", romsgrid->romsfile);
  cfid = romsgrid->cfid = dump_create_roms(romsgrid, romsgrid->ofile);
  write_roms_metrics(romsgrid, R_RHO);
  write_roms_metrics(romsgrid, R_PSI);
  write_roms_metrics(romsgrid, R_U);
  write_roms_metrics(romsgrid, R_V);

  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  count[0] = 1;
  count[1] = romsgrid->eta_rho;
  count[2] = romsgrid->xi_rho;
  nc_put_vara_double(cfid, ncw_var_id(cfid, "hraw"), start, count, romsgrid->hraw[0][0]);
  /*
    count[0] = romsgrid->eta_rho;
    count[1] = romsgrid->xi_rho;
    nc_put_vara_double(cfid, ncw_var_id(cfid, "h"), start, count, romsgrid->h[0]);
  */
  // Roms cannot handle NaN's in bathymetry
  write_roms_bathy(cfid, romsgrid->eta_rho, 
		   romsgrid->xi_rho, 1.0, romsgrid->h);

  count[0] = romsgrid->eta_rho;
  count[1] = romsgrid->xi_rho;
  nc_put_vara_double(cfid, ncw_var_id(cfid, "f"), start, count, romsgrid->f[0]);
  nc_put_vara_double(cfid, ncw_var_id(cfid, "pm"), start, count, romsgrid->pm[0]);
  nc_put_vara_double(cfid, ncw_var_id(cfid, "pn"), start, count, romsgrid->pn[0]);
  nc_put_vara_double(cfid, ncw_var_id(cfid, "dndx"), start, count, romsgrid->dndx[0]);
  nc_put_vara_double(cfid, ncw_var_id(cfid, "dmde"), start, count, romsgrid->dmde[0]);
  nc_put_vara_double(cfid, ncw_var_id(cfid, "angle"), start, count, romsgrid->theta[0]);

  start[0] = 0;
  start[1] = 0;
  count[0] = 1;
  v = (short)romsgrid->bmin;
  nc_put_vara_short(cfid, ncw_var_id(cfid, "depthmin"), start, count, &v);
  v = (short)romsgrid->bmax;
  nc_put_vara_short(cfid, ncw_var_id(cfid, "depthmax"), start, count, &v);

  /*
   * FR: Changed as advised by John Luick for MARVL 10/10/2014
   */
  // bv = 'F';
  bv = 'T';
  nc_put_var_text(cfid, ncw_var_id(cfid, "spherical"), &bv);

  nc_put_vara_double(cfid, ncw_var_id(cfid, "xl"), start, count, &romsgrid->xl);
  nc_put_vara_double(cfid, ncw_var_id(cfid, "el"), start, count, &romsgrid->el);

  count[0] = romsgrid->nz;
  nc_put_vara_double(cfid, ncw_var_id(cfid, "zt"), start, count, romsgrid->zt);
  nc_put_vara_double(cfid, ncw_var_id(cfid, "zb"), start, count, romsgrid->zb);

  close_roms_netcdf(romsgrid);
  if (DEBUG("init_m"))
    dlog("init_m", "\nROMS grid file %s created OK\n\n", romsgrid->ofile);

  /* Write ROMS IC file */
  if (!find_name_marvl(dumpdata, romsgrid->ofile, "roms_ic", NULL))
    sprintf(romsgrid->ofile, "%s_roms_ic.nc", romsgrid->romsfile);

  df_roms_write_ic(dumpdata, romsgrid->ofile);

  roms_grid_free(romsgrid);
}

/* END write_roms_grid()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to allocate memory for the parameter data structure       */
/*-------------------------------------------------------------------*/
romsgrid_t *romsgrid_alloc(void)
{
  romsgrid_t *romsgrid = (romsgrid_t *)malloc(sizeof(romsgrid_t));
  memset(romsgrid, 0, sizeof(romsgrid_t));
  return romsgrid;
}

/* END romsgrid_alloc()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Allocates memory for grid metric arrays                           */
/*-------------------------------------------------------------------*/
void alloc_grid_mem_roms(romsgrid_t *romsgrid, int nce1, int nce2)
{
  /* Allocate memory for the coordinate and metric arrays which are */
  /* at twice the final resolution.  */
  romsgrid->x = d_alloc_2d(2 * nce1 + 1, 2 * nce2 + 1);
  romsgrid->y = d_alloc_2d(2 * nce1 + 1, 2 * nce2 + 1);
  romsgrid->h1 = d_alloc_2d(2 * nce1 + 1, 2 * nce2 + 1);
  romsgrid->h2 = d_alloc_2d(2 * nce1 + 1, 2 * nce2 + 1);
  romsgrid->a1 = d_alloc_2d(2 * nce1 + 1, 2 * nce2 + 1);
  romsgrid->a2 = d_alloc_2d(2 * nce1 + 1, 2 * nce2 + 1);
  romsgrid->thetax = d_alloc_2d(nce1 + 1, nce2 + 1);
  romsgrid->thetay = d_alloc_2d(nce1 + 1, nce2 + 1);

  /* Allocate the ROMS metric and position arrays */
  romsgrid->lon_psi = d_alloc_2d(romsgrid->xi_psi, romsgrid->eta_psi);
  romsgrid->lat_psi = d_alloc_2d(romsgrid->xi_psi, romsgrid->eta_psi);
  romsgrid->lon_rho = d_alloc_2d(romsgrid->xi_rho, romsgrid->eta_rho);
  romsgrid->lat_rho = d_alloc_2d(romsgrid->xi_rho, romsgrid->eta_rho);
  romsgrid->lon_u = d_alloc_2d(romsgrid->xi_u, romsgrid->eta_u);
  romsgrid->lat_u = d_alloc_2d(romsgrid->xi_u, romsgrid->eta_u);
  romsgrid->lon_v = d_alloc_2d(romsgrid->xi_v, romsgrid->eta_v);
  romsgrid->lat_v = d_alloc_2d(romsgrid->xi_v, romsgrid->eta_v);
  romsgrid->hraw = d_alloc_3d(romsgrid->xi_rho, romsgrid->eta_rho, 1);
  romsgrid->h = d_alloc_2d(romsgrid->xi_rho, romsgrid->eta_rho);
  romsgrid->mask_rho = d_alloc_2d(romsgrid->xi_rho, romsgrid->eta_rho);
  romsgrid->mask_psi = d_alloc_2d(romsgrid->xi_psi, romsgrid->eta_psi);
  romsgrid->mask_u = d_alloc_2d(romsgrid->xi_u, romsgrid->eta_u);
  romsgrid->mask_v = d_alloc_2d(romsgrid->xi_v, romsgrid->eta_v);
  romsgrid->num_levels = d_alloc_2d(romsgrid->xi_rho, romsgrid->eta_rho);
  romsgrid->theta = d_alloc_2d(romsgrid->xi_rho, romsgrid->eta_rho);
  romsgrid->f = d_alloc_2d(romsgrid->xi_rho, romsgrid->eta_rho);
  romsgrid->pm = d_alloc_2d(romsgrid->xi_rho, romsgrid->eta_rho);
  romsgrid->pn = d_alloc_2d(romsgrid->xi_rho, romsgrid->eta_rho);
  romsgrid->dndx = d_alloc_2d(romsgrid->xi_rho, romsgrid->eta_rho);
  romsgrid->dmde = d_alloc_2d(romsgrid->xi_rho, romsgrid->eta_rho);

  /* Allocate open boundary arrays */
  if (romsgrid->nobc) {
    romsgrid->obcis = i_alloc_1d(romsgrid->nobc);
    romsgrid->obcjs = i_alloc_1d(romsgrid->nobc);
    romsgrid->obcie = i_alloc_1d(romsgrid->nobc);
    romsgrid->obcje = i_alloc_1d(romsgrid->nobc);
  }
}

/* END alloc_grid_mem()                                              */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Allocates memory for grid in the vertical                         */
/*-------------------------------------------------------------------*/
void alloc_grid_mem_vert_roms(romsgrid_t *romsgrid)
{
  romsgrid->zt = d_alloc_1d(romsgrid->nz+1);
  romsgrid->zb = d_alloc_1d(romsgrid->nz+1);
  romsgrid->layers = d_alloc_1d(romsgrid->nz + 1);
  
  /* Vertical coordinates */
  romsgrid->sc_rho = f_alloc_1d(romsgrid->nsigma);
  romsgrid->Cs_rho = f_alloc_1d(romsgrid->nsigma);
  romsgrid->sc_w   = f_alloc_1d(romsgrid->nsigma + 1);
  romsgrid->Cs_w   = f_alloc_1d(romsgrid->nsigma + 1);
}
/* END alloc_grid_mem_vert()                                         */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Frees the grid structure                                          */
/*-------------------------------------------------------------------*/
void roms_grid_free(romsgrid_t *romsgrid)
{
  d_free_2d(romsgrid->x);
  d_free_2d(romsgrid->y);
  d_free_2d(romsgrid->h1);
  d_free_2d(romsgrid->h2);
  d_free_2d(romsgrid->a1);
  d_free_2d(romsgrid->a2);
  
  d_free_1d(romsgrid->layers);
  d_free_1d(romsgrid->zb);
  d_free_1d(romsgrid->zt);
  d_free_3d(romsgrid->hraw);
  d_free_2d(romsgrid->h);
  d_free_2d(romsgrid->mask_rho);
  d_free_2d(romsgrid->mask_psi);
  d_free_2d(romsgrid->mask_u);
  d_free_2d(romsgrid->mask_v);
  d_free_2d(romsgrid->num_levels);
  d_free_2d(romsgrid->theta);
  d_free_2d(romsgrid->pm);
  d_free_2d(romsgrid->pn);
  d_free_2d(romsgrid->dndx);
  d_free_2d(romsgrid->dmde);
  if (romsgrid->topo != NULL)
    d_free_2d(romsgrid->topo);
  if (romsgrid->nobc) {
    i_free_1d(romsgrid->obcis);
    i_free_1d(romsgrid->obcjs);
    i_free_1d(romsgrid->obcie);
    i_free_1d(romsgrid->obcje);
  }
  d_free_2d(romsgrid->lon_psi);
  d_free_2d(romsgrid->lat_psi);
  d_free_2d(romsgrid->lon_rho);
  d_free_2d(romsgrid->lat_rho);
  d_free_2d(romsgrid->lon_u);
  d_free_2d(romsgrid->lat_u);
  d_free_2d(romsgrid->lon_v);
  d_free_2d(romsgrid->lat_v);
  f_free_1d(romsgrid->sc_rho);
  f_free_1d(romsgrid->Cs_rho);
  f_free_1d(romsgrid->sc_w);
  f_free_1d(romsgrid->Cs_w);
  
  free((romsgrid_t *)romsgrid);
}

/* END roms_grid_free()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Creates the ROMS output netCDF file                                */
/*-------------------------------------------------------------------*/
int dump_create_roms(romsgrid_t *romsgrid, char *name)
{
  int cdfid;
  int dims[10];
  int vid;
  /* dimension ids */
  int xi_rhoid;                 /* I dimension id at grid centre */
  int eta_rhoid;                /* J dimension id at grid centre */
  int xi_psiid;                 /* I dimension id at grid corner */
  int eta_psiid;                /* J dimension id at grid corner */
  int xi_uid;                   /* I dimension id at u edges     */
  int eta_uid;                  /* J dimension id at u edges     */
  int xi_vid;                   /* I dimension id at v edges     */
  int eta_vid;                  /* J dimension id at v edges     */
  int ztid;                     /* K dimension id at grid center */
  int zbid;                     /* K dimension id at grid face   */
  int oneid;
  int twoid;
  int bathid;

  /*-----------------------------------------------------------------*/
  /* Create the netCDF file                                          */
  if (ncw_create(name,(overwrite_output()?NC_CLOBBER:NC_NOCLOBBER), &cdfid) != NC_NOERR)
    hd_quit("Couldn't create output dump file %s\n", name);

  /*-----------------------------------------------------------------*/
  /* Define dimensions                                               */
  nc_def_dim(cdfid, "xi_rho", romsgrid->xi_rho, &xi_rhoid);
  nc_def_dim(cdfid, "eta_rho", romsgrid->eta_rho, &eta_rhoid);
  nc_def_dim(cdfid, "xi_psi", romsgrid->xi_psi, &xi_psiid);
  nc_def_dim(cdfid, "eta_psi", romsgrid->eta_psi, &eta_psiid);
  nc_def_dim(cdfid, "xi_u", romsgrid->xi_u, &xi_uid);
  nc_def_dim(cdfid, "eta_u", romsgrid->eta_u, &eta_uid);
  nc_def_dim(cdfid, "xi_v", romsgrid->xi_v, &xi_vid);
  nc_def_dim(cdfid, "eta_v", romsgrid->eta_v, &eta_vid);
  nc_def_dim(cdfid, "one", 1, &oneid);
  nc_def_dim(cdfid, "two", 2, &twoid);
  nc_def_dim(cdfid, "bath", 1, &bathid);

  nc_def_dim(cdfid, "zt", romsgrid->nz, &ztid);
  nc_def_dim(cdfid, "zb", romsgrid->nz, &zbid);

  /* Grid cells */
  dims[0] = eta_rhoid;
  dims[1] = xi_rhoid;
  nc_def_var(cdfid, "lon_rho", NC_DOUBLE, 2, dims, &vid);
  nc_def_var(cdfid, "lat_rho", NC_DOUBLE, 2, dims, &vid);
  dims[0] = eta_psiid;
  dims[1] = xi_psiid;
  nc_def_var(cdfid, "lon_psi", NC_DOUBLE, 2, dims, &vid);
  nc_def_var(cdfid, "lat_psi", NC_DOUBLE, 2, dims, &vid);
  dims[0] = eta_uid;
  dims[1] = xi_uid;
  nc_def_var(cdfid, "lon_u", NC_DOUBLE, 2, dims, &vid);
  nc_def_var(cdfid, "lat_u", NC_DOUBLE, 2, dims, &vid);
  dims[0] = eta_vid;
  dims[1] = xi_vid;
  nc_def_var(cdfid, "lon_v", NC_DOUBLE, 2, dims, &vid);
  nc_def_var(cdfid, "lat_v", NC_DOUBLE, 2, dims, &vid);

  /* Vertical grid */
  dims[0] = ztid;
  nc_def_var(cdfid, "zt", NC_DOUBLE, 1, dims, &vid);
  nc_def_var(cdfid, "zb", NC_DOUBLE, 1, dims, &vid);

  dims[0] = bathid;
  dims[1] = eta_rhoid;
  dims[2] = xi_rhoid;
  nc_def_var(cdfid, "hraw", NC_DOUBLE, 3, dims, &vid);
  dims[0] = eta_rhoid;
  dims[1] = xi_rhoid;
  nc_def_var(cdfid, "h", NC_DOUBLE, 2, dims, &vid);
  nc_def_var(cdfid, "mask_rho", NC_DOUBLE, 2, dims, &vid);
  dims[0] = eta_psiid;
  dims[1] = xi_psiid;
  nc_def_var(cdfid, "mask_psi", NC_DOUBLE, 2, dims, &vid);
  dims[0] = eta_uid;
  dims[1] = xi_uid;
  nc_def_var(cdfid, "mask_u", NC_DOUBLE, 2, dims, &vid);
  dims[0] = eta_vid;
  dims[1] = xi_vid;
  nc_def_var(cdfid, "mask_v", NC_DOUBLE, 2, dims, &vid);
  dims[0] = oneid;
  nc_def_var(cdfid, "depthmin", NC_SHORT, 1, dims, &vid);
  nc_def_var(cdfid, "depthmax", NC_SHORT, 1, dims, &vid);
  nc_def_var(cdfid, "spherical", NC_CHAR, 1, dims, &vid);
  dims[0] = eta_rhoid;
  dims[1] = xi_rhoid;
  nc_def_var(cdfid, "f", NC_DOUBLE, 2, dims, &vid);
  nc_def_var(cdfid, "pm", NC_DOUBLE, 2, dims, &vid);
  nc_def_var(cdfid, "pn", NC_DOUBLE, 2, dims, &vid);
  nc_def_var(cdfid, "dndx", NC_DOUBLE, 2, dims, &vid);
  nc_def_var(cdfid, "dmde", NC_DOUBLE, 2, dims, &vid);
  nc_def_var(cdfid, "angle", NC_DOUBLE, 2, dims, &vid);
  dims[0] = oneid;
  nc_def_var(cdfid, "xl", NC_DOUBLE, 1, dims, &vid);
  nc_def_var(cdfid, "el", NC_DOUBLE, 1, dims, &vid);

  write_roms_attributes(romsgrid, cdfid);

  nc_enddef(cdfid);
  return (cdfid);
}

/* END dump_create_roms()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Writes ROMS netCDF attributes                                      */
/*-------------------------------------------------------------------*/
void write_roms_attributes(romsgrid_t *romsgrid, int cdfid)
{
  /* dimension ids */
  int n, vid;
  char buf[MAXSTRLEN];
  char attname[MAXSTRLEN];

  /* time independent variables */
  vid = ncw_var_id(cdfid, "lon_rho");
  write_text_att(cdfid, vid, "long_name", "longitude of RHO-points");
  write_text_att(cdfid, vid, "units", "degree_east");
  vid = ncw_var_id(cdfid, "lat_rho");
  write_text_att(cdfid, vid, "long_name", "latitude of RHO-points");
  write_text_att(cdfid, vid, "units", "degree_east");

  vid = ncw_var_id(cdfid, "lon_psi");
  write_text_att(cdfid, vid, "long_name", "longitude of PSI-points");
  write_text_att(cdfid, vid, "units", "degree_east");
  vid = ncw_var_id(cdfid, "lat_psi");
  write_text_att(cdfid, vid, "long_name", "latitude of PSI-points");
  write_text_att(cdfid, vid, "units", "degree_east");

  vid = ncw_var_id(cdfid, "lon_u");
  write_text_att(cdfid, vid, "long_name", "longitude of U-points");
  write_text_att(cdfid, vid, "units", "degree_east");
  vid = ncw_var_id(cdfid, "lat_u");
  write_text_att(cdfid, vid, "long_name", "latitude of U-points");
  write_text_att(cdfid, vid, "units", "degree_east");

  vid = ncw_var_id(cdfid, "lon_v");
  write_text_att(cdfid, vid, "long_name", "longitude of V-points");
  write_text_att(cdfid, vid, "units", "degree_east");
  vid = ncw_var_id(cdfid, "lat_v");
  write_text_att(cdfid, vid, "long_name", "latitude of V-points");
  write_text_att(cdfid, vid, "units", "degree_east");

  vid = ncw_var_id(cdfid, "hraw");
  write_text_att(cdfid, vid, "long_name", "Working bathymetry at RHO-points");
  write_text_att(cdfid, vid, "units", "meter");
  write_text_att(cdfid, vid, "field", "bath, scalar");
  vid = ncw_var_id(cdfid, "h");
  write_text_att(cdfid, vid, "long_name", "Final bathymetry at RHO-points");
  write_text_att(cdfid, vid, "units", "meter");
  write_text_att(cdfid, vid, "field", "bath, scalar");

  vid = ncw_var_id(cdfid, "mask_rho");
  write_text_att(cdfid, vid, "long_name", "mask on RHO-points");
  write_text_att(cdfid, vid, "option(0)", "land");
  write_text_att(cdfid, vid, "option(1)", "water");
  vid = ncw_var_id(cdfid, "mask_psi");
  write_text_att(cdfid, vid, "long_name", "mask on PSI-points");
  write_text_att(cdfid, vid, "option(0)", "land");
  write_text_att(cdfid, vid, "option(1)", "water");
  vid = ncw_var_id(cdfid, "mask_u");
  write_text_att(cdfid, vid, "long_name", "mask on U-points");
  write_text_att(cdfid, vid, "option(0)", "land");
  write_text_att(cdfid, vid, "option(1)", "water");
  vid = ncw_var_id(cdfid, "mask_v");
  write_text_att(cdfid, vid, "long_name", "mask on V-points");
  write_text_att(cdfid, vid, "option(0)", "land");
  write_text_att(cdfid, vid, "option(1)", "water");

  vid = ncw_var_id(cdfid, "depthmin");
  write_text_att(cdfid, vid, "long_name", "Shallow bathymetry clipping depth");
  write_text_att(cdfid, vid, "units", "meter");
  vid = ncw_var_id(cdfid, "depthmax");
  write_text_att(cdfid, vid, "long_name", "Deep bathymetry clipping depth");
  write_text_att(cdfid, vid, "units", "meter");
  vid = ncw_var_id(cdfid, "spherical");
  write_text_att(cdfid, vid, "long_name", "Grid type logical switch");
  write_text_att(cdfid, vid, "option(T)", "spherical");
  write_text_att(cdfid, vid, "option(F)", "Cartesian");
  vid = ncw_var_id(cdfid, "f");
  write_text_att(cdfid, vid, "long_name", "Coriolis parameter at RHO-points");
  write_text_att(cdfid, vid, "units", "second-1");
  write_text_att(cdfid, vid, "field", "Coriolis, scalar");
  vid = ncw_var_id(cdfid, "pm");
  write_text_att(cdfid, vid, "long_name", "curvilinear coordinate metric in XI");
  write_text_att(cdfid, vid, "units", "meter-1");
  write_text_att(cdfid, vid, "field", "pm, scalar");
  vid = ncw_var_id(cdfid, "pn");
  write_text_att(cdfid, vid, "long_name", "curvilinear coordinate metric in ETA");
  write_text_att(cdfid, vid, "units", "meter-1");
  write_text_att(cdfid, vid, "field", "pn, scalar");
  vid = ncw_var_id(cdfid, "dndx");
  write_text_att(cdfid, vid, "long_name", "xi derivative of inverse metric factor pn");
  write_text_att(cdfid, vid, "units", "meter");
  write_text_att(cdfid, vid, "field", "dndx, scalar");
  vid = ncw_var_id(cdfid, "dmde");
  write_text_att(cdfid, vid, "long_name", "eta derivative of inverse metric factor pm");
  write_text_att(cdfid, vid, "units", "meter");
  write_text_att(cdfid, vid, "field", "dmde, scalar");
  vid = ncw_var_id(cdfid, "angle");
  write_text_att(cdfid, vid, "long_name", "angle between xi axis and east");
  write_text_att(cdfid, vid, "units", "radian");
  vid = ncw_var_id(cdfid, "xl");
  write_text_att(cdfid, vid, "long_name", "domain length in the XI-direction");
  write_text_att(cdfid, vid, "units", "meter");
  vid = ncw_var_id(cdfid, "el");
  write_text_att(cdfid, vid, "long_name", "domain length in the ETA-direction");
  write_text_att(cdfid, vid, "units", "meter");

  /* Vertical grid */
  vid = ncw_var_id(cdfid, "zt");
  write_text_att(cdfid, vid, "long_name", "zt");
  write_text_att(cdfid, vid, "units", "meters");
  write_text_att(cdfid, vid, "cartesian_axis", "z");
  write_text_att(cdfid, vid, "positive", "down");
  vid = ncw_var_id(cdfid, "zb");
  write_text_att(cdfid, vid, "long_name", "zb");
  write_text_att(cdfid, vid, "units", "meters");
  write_text_att(cdfid, vid, "cartesian_axis", "z");
  write_text_att(cdfid, vid, "positive", "down");
    

  /* Global attributes */
  write_text_att(cdfid, NC_GLOBAL, "filename", romsgrid->ofile);
  write_text_att(cdfid, NC_GLOBAL, "type", "ROMS GRID file");
}


/*-------------------------------------------------------------------*/
/* Writes ROMS metric data to netCDF file                             */
/*-------------------------------------------------------------------*/
void write_roms_metrics(romsgrid_t *romsgrid, int mode)
{
  char grd[4];
  size_t start[4];
  size_t count[4];
  char var[MAXSTRLEN];
  int cfid = romsgrid->cfid;
  float vertex[4] = {1, 2, 3, 4};
  int nx = 0, ny = 0;
  double **lat = NULL, **lon = NULL, **mask = NULL;

  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;
  count[0] = 0;
  count[1] = 0;
  count[2] = 0;
  count[3] = 0;
  switch(mode) {
  case R_RHO:
    strcpy(grd, "rho");
    nx = romsgrid->xi_rho;
    ny = romsgrid->eta_rho;
    lat = romsgrid->lat_rho;
    lon = romsgrid->lon_rho;
    mask = romsgrid->mask_rho;
    break;
  case R_PSI:
    strcpy(grd, "psi");
    nx = romsgrid->xi_psi;
    ny = romsgrid->eta_psi;
    lat = romsgrid->lat_psi;
    lon = romsgrid->lon_psi;
    mask = romsgrid->mask_psi;
    break;
  case R_U:
    strcpy(grd, "u");
    nx = romsgrid->xi_u;
    ny = romsgrid->eta_u;
    lat = romsgrid->lat_u;
    lon = romsgrid->lon_u;
    mask = romsgrid->mask_u;
    break;
  case R_V:
    strcpy(grd, "v");
    nx = romsgrid->xi_v;
    ny = romsgrid->eta_v;
    lat = romsgrid->lat_v;
    lon = romsgrid->lon_v;
    mask = romsgrid->mask_v;
    break;
  }

  count[0] = ny;
  count[1] = nx;
  sprintf(var, "lat_%s",grd);
  nc_put_vara_double(cfid, ncw_var_id(cfid, var), start, count, lat[0]);
  sprintf(var, "lon_%s",grd);
  nc_put_vara_double(cfid, ncw_var_id(cfid, var), start, count, lon[0]);
  sprintf(var, "mask_%s",grd);
  nc_put_vara_double(cfid, ncw_var_id(cfid, var), start, count, mask[0]);

}
/* END write_roms_metrics()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Closes a ROMS netCDF file                                          */
/*-------------------------------------------------------------------*/
void close_roms_netcdf(romsgrid_t *romsgrid)
{
  size_t start[4];
  size_t count[4];
  int cfid = romsgrid->cfid;

  ncw_close(romsgrid->ofile, romsgrid->cfid);
}

/* END close_netcdf()                                                */
/*-------------------------------------------------------------------*/
