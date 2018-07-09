/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/outputs/df_mom.c
 *  
 *  Description:
 *  Routines for meco model to allow handling
 *  multiple netCDF output files, with different
 *  sets of time dependent variables and different
 *  output schedules.
 *  Output is written in MOM compatible netCDF.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: mom_grid.c 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

#include <stdlib.h>
#include <netcdf.h>
#include <string.h>
#include "hd.h"

/*-------------------------------------------------------------------*/
/* Flags                                                             */
/*-------------------------------------------------------------------*/
/* MOM cell flags                                                    */
#define TCELL   1
#define NCELL   2
#define ECELL   4
#define CCELL   8

/* Bathymetry and layer structure                                    */
#define MAXGRIDZ 1e20
#define NOTVALID -9999
#define LANDCELL 99
#define MAXDEPTH 9000
#define NOTWET (SOLID|OUTSIDE)

/* Size of EMS grid for MOM - EMS conversion.                       */
/* Limit of EMS grid relative to MOM grid (e.g. include CCELLS is   */
/* the EMS grid: gsize = CCELL or use only MOM TCELLS in the EMS    */
/* grid: gsize = TCELL).                                            */
int gsize = TCELL;

/*------------------------------------------------------------------*/
/* Private routines                                                 */
/*------------------------------------------------------------------*/
momgrid_t *momgrid_alloc(void);
void alloc_grid_mem(momgrid_t *momgrid, int nce1, int nce2);
void alloc_mom_mem(momgrid_t *momgrid, int nce1, int nce2);
void mom_free(momgrid_t *momgrid);
int read2darray(FILE * fp, char *label, double **array, long n1, long n2);
void get_metrics(momgrid_t *momgrid, int mode);
void get_nom_ll(momgrid_t *momgrid);
void print_metrics_cdl(momgrid_t *momgrid, FILE *op, int mode);
void print_array_1d(FILE *op, int nce1, double *a);
void print_farray_1d(FILE *op, int nce1, float *a);
void print_array_2d(FILE *op, int nce1, int nce2, double **a);
void print_array_3d(FILE *op, int nce1, int nce2, int nk, double ***a);
void print_grid_header(momgrid_t *momgrid, FILE *op, int mode);
void print_grid_atts(momgrid_t *momgrid, FILE *op, int mode);
void close_cdl(momgrid_t *momgrid, FILE *op);
int dump_create_mom(momgrid_t *momgrid, char *name, int mode);
void write_attributes(momgrid_t *momgrid, int cdfid, int mode);
void write_text_att(int cdfid, int varid, const char *name,
                    const char *text);
void write_metrics(momgrid_t *momgrid, int mode);
void write_TS(dump_data_t *dumpdata, momgrid_t *momgrid);
void write_eta(dump_data_t *dumpdata, momgrid_t *momgrid);
void close_netcdf(momgrid_t *momgrid);
void print_ems_grid(momgrid_t *momgrid, FILE *op, char *fname);
void get_ij(momgrid_t *momgrid, int nx, int ny, int i, int j, int *ii, int *jj);
int find_name_marvl(dump_data_t *dumpdata, char *ofile, char *tag, int *p_n);

/*------------------------------------------------------------------*/
/* Converts EMS grid to MOM grid                                    */
/*------------------------------------------------------------------*/
momgrid_t *convert_mom_grid(parameters_t *params, dump_data_t *dumpdata)
{
  momgrid_t *momgrid;           /* Grid data                        */
  int n, i, j, ip, jp, k;       /* Counters                         */
  int nce1, nce2;               /* x and y grid size                */

  if (params->domom == 0) return(NULL);

  /*-----------------------------------------------------------------*/
  /* Set up the grid structure */
  momgrid = momgrid_alloc();
  strcpy(momgrid->momfile, params->momfile);
  momgrid->domom = params->domom;

  /*-----------------------------------------------------------------*/
  /* Read grid attributes */
  momgrid->nce1 = params->nce1;
  momgrid->nce2 = params->nce2;
  momgrid->nfe1 = momgrid->nce1 + 1;
  momgrid->nfe2 = momgrid->nce2 + 1;
  nce1 = momgrid->nce1;
  nce2 = momgrid->nce2;
  momgrid->grid_x_T = momgrid->nce1 - 1;
  momgrid->grid_y_T = momgrid->nce2 - 1;
  momgrid->grid_x_C = momgrid->grid_x_T;
  momgrid->grid_y_C = momgrid->grid_y_T;
  momgrid->vertex = 4;
  momgrid->nobc = params->nobc;

  /*-----------------------------------------------------------------*/
  /* Check for cyclic boundaries                                     */
  momgrid->cyoffs = momgrid->cytype = 0;
  i = j = 0;
  if (params->nemx) {
    momgrid->grid_x_T += 1;
    momgrid->grid_x_C += 1;
    momgrid->cyoffs =  params->emjsx[0] - params->emjdx[0];
    momgrid->cytype |= U1BDRY;
    i = 2;
  }
  if (params->nemy) {
    momgrid->grid_y_T += 1;
    momgrid->grid_y_C += 1;
    momgrid->cyoffs = params->emisx[0] - params->emidx[0];
    momgrid->cytype |= U2BDRY;
    j = 2;
  }
  for (n = 0; n < momgrid->nobc; n++) {
    if (params->open[n]->bcond_ele & CYCLIC) {
      if (i == 0 && params->open[n]->type & U1BDRY) {
	momgrid->grid_x_T += 1;
	momgrid->grid_x_C += 1;
	momgrid->cyoffs = params->open[n]->jloc[0];
	ip = params->open[n]->iloc[0];
	momgrid->cytype |= U1BDRY;
	i = 1;
      } else if (i == 1 && params->open[n]->type & U1BDRY) {
	momgrid->cyoffs = (params->open[n]->iloc[0] > ip) ? 
	  momgrid->cyoffs - params->open[n]->jloc[0] : 
	  params->open[n]->jloc[0] - momgrid->cyoffs;
	i = 2;
      }
      if (j == 0 && params->open[n]->type & U2BDRY) {
	momgrid->grid_y_T += 1;
	momgrid->grid_y_C += 1;
	momgrid->cyoffs = params->open[n]->iloc[0];
	jp = params->open[n]->jloc[0];
	momgrid->cytype |= U2BDRY;
	j = 1;
      } else if (j == 1 && params->open[n]->type & U2BDRY) {
	momgrid->cyoffs = (params->open[n]->jloc[0] > jp) ? 
	  momgrid->cyoffs - params->open[n]->iloc[0] : 
	  params->open[n]->iloc[0] - momgrid->cyoffs;
	j = 2;
      }
    }
  }
  if (i != 2)
    hd_warn("mom_grid: can't find cyclic boundaries in the x direction.\n");
  if (j != 2)
    hd_warn("mom_grid: can't find cyclic boundaries in the y direction.\n");

  alloc_grid_mem(momgrid, nce1, nce2);
  if (params->domom & RDGRID) read_grid(params);

  if (DEBUG("init_m")) {
    dlog("init_m", "SHOC grid = (%d x %d)\n",momgrid->nce1, momgrid->nce2);
    dlog("init_m", "MOM grid = (%d x %d)\n",momgrid->grid_x_T, momgrid->grid_y_T);
    if (momgrid->cytype & U1BDRY)
      dlog("init_m", "Cyclic boundary in x direction\n",momgrid->nce1, momgrid->nce2);
    if (momgrid->cytype & U2BDRY)
      dlog("init_m", "Cyclic boundary in y direction\n",momgrid->nce1, momgrid->nce2);
    if (momgrid->cyoffs)
      dlog("init_m", "Cyclic boundary offset = %d\n", momgrid->cyoffs);
  }

  /*-----------------------------------------------------------------*/
  /* Get the obc info. Note; MOM origin is (1,1), SHOC is (0,0), so  */
  /* add 1 to the SHOC OBC indices.                                  */
  for (n = 0; n < momgrid->nobc; n++) {
    i = j = 1;
    if (geom->open[n]->ocodex & R_EDGE) i = 0;
    if (geom->open[n]->ocodey & F_EDGE) j = 0;
    momgrid->obcis[n] = min(params->open[n]->iloc[0] + i, momgrid->grid_x_T);
    momgrid->obcjs[n] = min(params->open[n]->jloc[0] + j, momgrid->grid_y_T);
    momgrid->obcie[n] = min(params->open[n]->iloc[params->open[n]->npts - 1] + i,
			    momgrid->grid_x_T);
    momgrid->obcje[n] = min(params->open[n]->jloc[params->open[n]->npts - 1] + j,
			    momgrid->grid_y_T);
  }

  /*-----------------------------------------------------------------*/
  /* Get the grid metrics                                            */
  for (j = 0; j <= 2*momgrid->nce2; j++)
    for (i = 0; i <= 2*momgrid->nce1; i++) {
      momgrid->x[j][i] = params->x[j][i];
      momgrid->y[j][i] = params->y[j][i];
      momgrid->h1[j][i] = params->h1[j][i];
      momgrid->h2[j][i] = params->h2[j][i];
      momgrid->a1[j][i] = params->a1[j][i];
      momgrid->a2[j][i] = params->a2[j][i];
    }

  /* Insert nearest neighbours for NaNs                              */
  for (j = 0; j <= 2*nce2; j++)
    for (i = 0; i <= 2*nce1; i++) {
      if(isnan(momgrid->x[j][i])) {
	find_non_nan(momgrid->x, 2*nce1, 2*nce2, i, j, &ip, &jp);
	momgrid->x[j][i] = momgrid->x[jp][ip];
      }
      if(isnan(momgrid->y[j][i])) {
	find_non_nan(momgrid->y, 2*nce1, 2*nce2, i, j, &ip, &jp);
	momgrid->y[j][i] = momgrid->y[jp][ip];
      }
      if(isnan(momgrid->h1[j][i])) {
	find_non_nan(momgrid->h1, 2*nce1, 2*nce2, i, j, &ip, &jp);
	momgrid->h1[j][i] = momgrid->h1[jp][ip];
      }
      if(isnan(momgrid->h2[j][i])) {
	find_non_nan(momgrid->h2, 2*nce1, 2*nce2, i, j, &ip, &jp);
	momgrid->h2[j][i] = momgrid->h2[jp][ip];
      }
      if(isnan(momgrid->a1[j][i])) {
	find_non_nan(momgrid->a1, 2*nce1, 2*nce2, i, j, &ip, &jp);
	momgrid->a1[j][i] = momgrid->a1[jp][ip];
      }
      if(isnan(momgrid->a2[j][i])) {
	find_non_nan(momgrid->a2, 2*nce1, 2*nce2, i, j, &ip, &jp);
	momgrid->a2[j][i] = momgrid->a2[jp][ip];
      }
    }

  /* Insert dummy values for NaNs                                    */
  /*
  for (j = 0; j <= 2*momgrid->nce2; j++)
    for (i = 0; i <= 2*momgrid->nce1; i++) {
      if (isnan(momgrid->x[j][i]) && i > 1)
	momgrid->x[j][i] = 2.0*momgrid->x[j][i-1]-momgrid->x[j][i-2]; 
      if (isnan(momgrid->h1[j][i]) && i > 1)
	momgrid->h1[j][i] = 2.0*momgrid->h1[j][i-1]-momgrid->h1[j][i-2]; 
      if (isnan(momgrid->a1[j][i]) && i > 1)
	momgrid->a1[j][i] = 2.0*momgrid->a1[j][i-1]-momgrid->a1[j][i-2];
      if (isnan(momgrid->x[j][i]) && j > 1)
	momgrid->x[j][i] = 2.0*momgrid->x[j-1][i]-momgrid->x[j-2][i];
      if (isnan(momgrid->h1[j][i]) && j > 1)
	momgrid->h1[j][i] = 2.0*momgrid->h1[j-1][i]-momgrid->h1[j-2][i];
      if (isnan(momgrid->a1[j][i]) && j > 1)
	momgrid->a1[j][i] = 2.0*momgrid->a1[j-1][i]-momgrid->a1[j-2][i];

      if (isnan(momgrid->y[j][i]) && j > 1)
	momgrid->y[j][i] = 2.0*momgrid->y[j-1][i]-momgrid->y[j-2][i];

      if (isnan(momgrid->h2[j][i]) && j > 1)
	momgrid->h2[j][i] = 2.0*momgrid->h2[j-1][i]-momgrid->h2[j-2][i];
      if (isnan(momgrid->a2[j][i]) && j > 1)
	momgrid->a2[j][i] = 2.0*momgrid->a2[j-1][i]-momgrid->a2[j-2][i];

      if (isnan(momgrid->y[j][i]) && i > 1)
	momgrid->y[j][i] = 2.0*momgrid->y[j][i-1]-momgrid->y[j][i-2];

      if (isnan(momgrid->h2[j][i]) && i > 1)
	momgrid->h2[j][i] = 2.0*momgrid->h2[j][i-1]-momgrid->h2[j][i-2];
      if (isnan(momgrid->a2[j][i]) && i > 1)
	momgrid->a2[j][i] = 2.0*momgrid->a2[j][i-1]-momgrid->a2[j][i-2];
    }
  */
  /*-----------------------------------------------------------------*/
  /* Read the bathymetry                                             */
  for (j = 0; j < momgrid->grid_y_T; j++)
    for (i = 0; i < momgrid->grid_x_T; i++) {
      double depth = dumpdata->botz[j][i];
      if (depth == LANDCELL || depth == NOTVALID || isnan(depth)) {
	momgrid->depth[j][i] = 0.0;
	momgrid->wet[j][i] = 0.0;
        momgrid->num_levels[j][i] = 0.0;
      }
      else {
	momgrid->depth[j][i] = -depth;
	momgrid->wet[j][i] = 1.0;
        momgrid->num_levels[j][i] = -1.0;
	for (k = 0; k < params->nz; k++) {
	  if (depth >= params->layers[k] && depth < params->layers[k + 1])
	    momgrid->num_levels[j][i] = (double)(params->nz - k);
	}
      }
      if (momgrid->num_levels[j][i] < 0)
	hd_quit("Bottom not found at T cell (%d %d) : %5.1f m\n",i, j, depth);
    }
  for (j = 0; j < momgrid->grid_y_C; j++)
    for (i = 0; i < momgrid->grid_x_C; i++) {
      double depth, d1, d2, d3, d4;
      unsigned long ***flag = dumpdata->flag;
      int nz = params->nz - 1;
      ip = (i < momgrid->nce1 - 1) ? i + 1 : i;
      jp = (j < momgrid->nce2 - 1) ? j + 1 : j;

      /* Reset ip,jp for OBCs adjacent to OUTSIDE cells in interior  */
      if (flag[nz][j][ip] & U1BDRY && flag[nz][j][ip] & R_EDGE) ip = i;
      if (flag[nz][jp][i] & U2BDRY && flag[nz][jp][i] & F_EDGE) jp = j;

      d1 = dumpdata->botz[j][i];
      d2 = dumpdata->botz[j][ip];
      d3 = dumpdata->botz[jp][i]; 
      d4 = dumpdata->botz[jp][ip];
      depth = max(d1,d2);
      depth = max(depth,d3);
      depth = max(depth,d4);
      if (d1 == LANDCELL || d2 == LANDCELL || d3 == LANDCELL || 
	  d4 == LANDCELL || 
	  d1 == NOTVALID || d2 == NOTVALID || d3 == NOTVALID || 
	  d4 == NOTVALID || 
	  isnan(d1) || isnan(d2) || isnan(d3) || isnan(d4)) {
	momgrid->depth_c[j][i] = 0.0;
	momgrid->wet_c[j][i] = 0.0;
        momgrid->num_levels_c[j][i] = 0.0;
      }
      else {
	momgrid->depth_c[j][i] = -depth;
	momgrid->wet_c[j][i] = 1.0;
        momgrid->num_levels_c[j][i] = -1.0;
	for (k = 0; k < params->nz; k++) {
	  if (depth >= params->layers[k] && depth < params->layers[k + 1])
	    momgrid->num_levels_c[j][i] = (double)(params->nz - k);
	}
      }
      if (momgrid->num_levels_c[j][i] < 0)
	hd_quit("Bottom not found at C cell (%d %d) : %5.1f m\n",i, j, depth);
    }

  /* Reset the depth adjacent to outside cells on OBC's              */
  for (n = 0; n < momgrid->nobc; n++) {
    int nm = 2;  /* Number of cells to maintain uniform depth        */
    int nn, sc;
    for (k = 0; k < params->open[n]->npts; k++) {
      i = ip = params->open[n]->iloc[k];
      j = jp = params->open[n]->jloc[k];      
      /* U1 boundaries                                               */
      if (geom->open[n]->type & U1BDRY) {
	sc = 1;
	if (geom->open[n]->ocodex & L_EDGE) {
	  ip = i - 1;
	  sc = -1;
	} else
	  i --;
	for (nn = 0; nn < nm; nn++) {
	  ip += sc * nn;
	  if (ip > 0 && ip < momgrid->grid_x_T && j < momgrid->grid_y_T) {
	    if (momgrid->depth[j][ip] == 0.0) {
	      momgrid->depth[j][ip] = momgrid->depth[j][i];
	      momgrid->wet[j][ip] = 1.0;
	      momgrid->num_levels[j][ip] = momgrid->num_levels[j][i];
	    }
	    if (momgrid->depth_c[j][ip] == 0.0) {
	      momgrid->depth_c[j][ip] = momgrid->depth_c[j][i];
	      momgrid->wet_c[j][ip] = 1.0;
	      momgrid->num_levels_c[j][ip] = momgrid->num_levels_c[j][i];
	    }
	  }
	}
      }
      /* U2 boundaries                                               */
      if (geom->open[n]->type & U2BDRY) {
	sc = 1;
	if (geom->open[n]->ocodey & B_EDGE) {
	  jp = j - 1;
	  sc = -1;
	} else
	  j--;
	for (nn = 0; nn < nm; nn++) {
	  jp += sc * nn;
	  if (jp > 0 && jp < momgrid->grid_y_T && i < momgrid->grid_x_T) {
	    if (momgrid->depth[jp][i] == 0.0) {
	      momgrid->depth[jp][i] = momgrid->depth[j][i];
	      momgrid->wet[jp][i] = 1.0;
	      momgrid->num_levels[jp][i] = momgrid->num_levels[j][i];
	    }
	    if (momgrid->depth_c[jp][i] == 0.0) {
	      momgrid->depth_c[jp][i] = momgrid->depth_c[j][i];
	      momgrid->wet_c[jp][i] = 1.0;
	      momgrid->num_levels_c[jp][i] = momgrid->num_levels_c[j][i];
	    }
	  }
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set up the vertical grid                                        */
  momgrid->nz = params->nz;
  momgrid->layers = d_alloc_1d(momgrid->nz + 1);
  memcpy(momgrid->layers, params->layers, (params->nz + 1) * sizeof(double));
  momgrid->layers[momgrid->nz] = 0.0;

  /* Allocate memory for arrays with vertical dimension              */
  if (momgrid->nz) {
    momgrid->zt = f_alloc_1d(momgrid->nz);
    momgrid->zb = f_alloc_1d(momgrid->nz);
    momgrid->T3 = d_alloc_3d(momgrid->grid_x_T, momgrid->grid_y_T, 
			     momgrid->nz);
    momgrid->C3 = d_alloc_3d(momgrid->grid_x_C, momgrid->grid_y_C, 
			     momgrid->nz);
  }
  for (k = 0; k < momgrid->nz; k++) {
    i = momgrid->nz - 1 - k;
    momgrid->zb[i] = (float)fabs(momgrid->layers[k]);
    momgrid->zt[i] = (float)fabs(0.5 * (momgrid->layers[k] + momgrid->layers[k + 1]));
  }

  return(momgrid);
}

/* END convert_mom_grid()                                           */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Writes EMS grid to MOM grid conversion to netCDF file            */
/*------------------------------------------------------------------*/
void write_mom_grid(dump_data_t *dumpdata)
{
  momgrid_t *momgrid = dumpdata->momgrid;
  FILE *op = NULL;
  int i, j;                     /* Counters                         */
  int mode;                     /* Output mode                      */

  if (momgrid == NULL) return;
  mode = momgrid->domom;

  /*----------------------------------------------------------------*/
  /* Convert metrics to MOM format and write to file                */
  if (mode & CDL) {
    sprintf(momgrid->ofile, "%s_spec.cdl", momgrid->momfile);
    if((op = fopen(momgrid->ofile, "w")) == NULL)
      hd_quit("Can't open file %s\n", momgrid->ofile);
    print_grid_header(momgrid, op, 1);
    print_grid_atts(momgrid, op, TCELL);
    print_grid_atts(momgrid, op, ECELL);
    print_grid_atts(momgrid, op, NCELL);
    print_grid_atts(momgrid, op, CCELL);
  } else if (mode & CDF) {
    if (!find_name_marvl(dumpdata, momgrid->ofile, "mom_spec", NULL))
      sprintf(momgrid->ofile, "%s_spec.nc", momgrid->momfile);
    momgrid->cfid = dump_create_mom(momgrid, momgrid->ofile, METRIC);
  }

  get_nom_ll(momgrid);
  get_metrics(momgrid, TCELL);
  for (j = 0; j < momgrid->grid_y_T; j++)
    for (i = 0; i < momgrid->grid_x_T; i++) {
      momgrid->x_T[j][i] = momgrid->x_g[j][i];
      momgrid->y_T[j][i] = momgrid->y_g[j][i];
    }

  if (mode & CDL) {
    print_grid_header(momgrid, op, 0);
    print_metrics_cdl(momgrid, op, TCELL);
  } else if (mode & CDF) 
    write_metrics(momgrid, TCELL);

  get_metrics(momgrid, ECELL);
  if (mode & CDL)
    print_metrics_cdl(momgrid, op, ECELL);
  else if (mode & CDF) 
    write_metrics(momgrid, ECELL);

  get_metrics(momgrid, NCELL);
  if (mode & CDL)
    print_metrics_cdl(momgrid, op, NCELL);
  else if (mode & CDF) 
    write_metrics(momgrid, NCELL);

  get_metrics(momgrid, CCELL);
  for (j = 0; j < momgrid->grid_y_T; j++)
    for (i = 0; i < momgrid->grid_x_T; i++) {
      momgrid->x_C[j][i] = momgrid->x_g[j][i];
      momgrid->y_C[j][i] = momgrid->y_g[j][i];
    }

  if (mode & CDL)
    print_metrics_cdl(momgrid, op, CCELL);
  else if (mode & CDF) 
    write_metrics(momgrid, CCELL);

  if (mode & CDL)
    close_cdl(momgrid, op);
  else if (mode & CDF)
    close_netcdf(momgrid);
  if (DEBUG("init_m"))
    dlog("init_m", "\nGridspec %s created OK\n\n", momgrid->ofile);

  /* Write the T/S initialisation data to file                      */
  if (mode & CDF) {
    landfillfn_t landfill;
    landfill = locate_landfill_function("cascade_search");
    /*landfill = locate_landfill_function("cascade_search_orig");*/
    /*landfill = locate_landfill_function("default");*/
    landfill(dumpdata);

    if (!find_name_marvl(dumpdata, momgrid->ofile, "mom_ic_temp_salt", NULL))
      sprintf(momgrid->ofile, "ocean_temp_salt.res.nc");
    momgrid->cfid = dump_create_mom(momgrid, momgrid->ofile, TEMSAL);
    write_TS(dumpdata, momgrid);
    if (DEBUG("init_m"))
      dlog("init_m", "\nT/S file %s created OK\n\n", momgrid->ofile);

    momgrid->ofile[0] = '\0';
    if (!find_name_marvl(dumpdata, momgrid->ofile, "mom_ic_eta", NULL) &&
	!(params->runmode & PRE_MARVL))
      sprintf(momgrid->ofile, "eta_t_ic.nc");
    if (strlen(momgrid->ofile)) {
      momgrid->cfid = dump_create_mom(momgrid, momgrid->ofile, ETA_A);
      write_eta(dumpdata, momgrid);
      if (DEBUG("init_m"))
	dlog("init_m", "\neta file %s created OK\n\n", momgrid->ofile);
    }
  }

  mom_free(momgrid);
  mom_grid_free(momgrid, UNUSED);
}

/* END write_mom_grid()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Creates the metrics for a cell and writes to file                 */
/*-------------------------------------------------------------------*/
void get_metrics(momgrid_t *momgrid, int mode)
{
  int i, j, io, jo;
  int im, ic, ip;
  int jm, jc, jp;
  int i0, j0;
  int i1, j1;
  int nx, ny;

  switch(mode) {
  case TCELL:
    io = 0;
    jo = 0;
    nx = momgrid->grid_x_T;
    ny = momgrid->grid_y_T;
    alloc_mom_mem(momgrid, nx, ny);
    break;
  case NCELL:
    io = 0;
    jo = 1;
    nx = momgrid->grid_x_T;
    ny = momgrid->grid_y_C;
    mom_free(momgrid);
    alloc_mom_mem(momgrid, nx, ny);
    break;
  case ECELL:
    io = 1;
    jo = 0;
    nx = momgrid->grid_x_C;
    ny = momgrid->grid_y_T;
    mom_free(momgrid);
    alloc_mom_mem(momgrid, nx, ny);
    break;
  case CCELL:
    io = 1;
    jo = 1;
    nx = momgrid->grid_x_C;
    ny = momgrid->grid_y_C;
    mom_free(momgrid);
    alloc_mom_mem(momgrid, nx, ny);
    break;
  default:
    io = 0;
    jo = 0;
    nx = 0;
    ny = 0;
  }

  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {

      im = io+i*2;
      ic = io+i*2+1;
      ip = io+i*2+2;

      jm = jo+j*2;
      jc = jo+j*2+1;
      jp = jo+j*2+2;

      get_ij(momgrid, nx, ny, ic, jc, &i0, &j0); 
      momgrid->x_g[j][i] = momgrid->x[j0][i0];
      momgrid->y_g[j][i] = momgrid->y[j0][i0];
      /* Verticies begins in SW corner, direction counterclockwise */
      get_ij(momgrid, nx, ny, im, jm, &i0, &j0); 
      momgrid->x_vert[0][j][i] = momgrid->x[j0][i0];
      momgrid->y_vert[0][j][i] = momgrid->y[j0][i0];
      get_ij(momgrid, nx, ny, ip, jm, &i0, &j0); 
      momgrid->x_vert[1][j][i] = momgrid->x[j0][i0];
      momgrid->y_vert[1][j][i] = momgrid->y[j0][i0];
      get_ij(momgrid, nx, ny, ip, jp, &i0, &j0); 
      momgrid->x_vert[2][j][i] = momgrid->x[j0][i0];
      momgrid->y_vert[2][j][i] = momgrid->y[j0][i0];
      get_ij(momgrid, nx, ny, im, jp, &i0, &j0); 
      momgrid->x_vert[3][j][i] = momgrid->x[j0][i0];
      momgrid->y_vert[3][j][i] = momgrid->y[j0][i0];

      get_ij(momgrid, nx, ny, im, jc, &i0, &j0); 
      momgrid->ds_00_02[j][i] = 2.0 * momgrid->h2[j0][i0];
      get_ij(momgrid, nx, ny, ip, jc, &i0, &j0); 
      momgrid->ds_20_22[j][i] = 2.0 * momgrid->h2[j0][i0];
      get_ij(momgrid, nx, ny, ic, jp, &i0, &j0); 
      momgrid->ds_02_22[j][i] = 2.0 * momgrid->h1[j0][i0];
      get_ij(momgrid, nx, ny, ic, jm, &i0, &j0); 
      momgrid->ds_00_20[j][i] = 2.0 * momgrid->h1[j0][i0];
      get_ij(momgrid, nx, ny, im, jc, &i0, &j0);
      get_ij(momgrid, nx, ny, im, jm, &i1, &j1);
      momgrid->ds_00_01[j][i] = 0.5 * (momgrid->h2[j0][i0] +
				       momgrid->h2[j1][i1]);
      get_ij(momgrid, nx, ny, im, jp, &i0, &j0);
      get_ij(momgrid, nx, ny, im, jc, &i1, &j1);
      momgrid->ds_01_02[j][i] = 0.5 * (momgrid->h2[j0][i0] +
				       momgrid->h2[j1][i1]);
      get_ij(momgrid, nx, ny, ic, jp, &i0, &j0);
      get_ij(momgrid, nx, ny, im, jp, &i1, &j1);
      momgrid->ds_02_12[j][i] = 0.5 * (momgrid->h1[j0][i0] +
				       momgrid->h1[j1][i1]);
      get_ij(momgrid, nx, ny, ip, jp, &i0, &j0);
      get_ij(momgrid, nx, ny, ic, jp, &i1, &j1);
      momgrid->ds_12_22[j][i] = 0.5 * (momgrid->h1[j0][i0] +
				       momgrid->h1[j1][i1]);
      get_ij(momgrid, nx, ny, ip, jp, &i0, &j0);
      get_ij(momgrid, nx, ny, ip, jc, &i1, &j1);
      momgrid->ds_21_22[j][i] = 0.5 * (momgrid->h2[j0][i0] +
				       momgrid->h2[j1][i1]);
      get_ij(momgrid, nx, ny, ip, jc, &i0, &j0);
      get_ij(momgrid, nx, ny, ip, jm, &i1, &j1);
      momgrid->ds_20_21[j][i] = 0.5 * (momgrid->h2[j0][i0] +
				       momgrid->h2[j1][i1]);
      get_ij(momgrid, nx, ny, ip, jm, &i0, &j0);
      get_ij(momgrid, nx, ny, ic, jm, &i1, &j1);
      momgrid->ds_10_20[j][i] = 0.5 * (momgrid->h1[j0][i0] +
				       momgrid->h1[j1][i1]);
      get_ij(momgrid, nx, ny, ic, jm, &i0, &j0);
      get_ij(momgrid, nx, ny, im, jm, &i1, &j1);
      momgrid->ds_00_10[j][i] = 0.5 * (momgrid->h1[j0][i0] +
				       momgrid->h1[j1][i1]);
      get_ij(momgrid, nx, ny, ic, jc, &i0, &j0);
      get_ij(momgrid, nx, ny, im, jc, &i1, &j1);
      momgrid->ds_01_11[j][i] = 0.5 * (momgrid->h1[j0][i0] +
				       momgrid->h1[j1][i1]);
      get_ij(momgrid, nx, ny, ic, jp, &i0, &j0);
      get_ij(momgrid, nx, ny, ic, jc, &i1, &j1);
      momgrid->ds_11_12[j][i] = 0.5 * (momgrid->h2[j0][i0] +
				       momgrid->h2[j1][i1]);
      get_ij(momgrid, nx, ny, ip, jc, &i0, &j0);
      get_ij(momgrid, nx, ny, ic, jc, &i1, &j1);
      momgrid->ds_11_21[j][i] = 0.5 * (momgrid->h1[j0][i0] +
				       momgrid->h1[j1][i1]);
      get_ij(momgrid, nx, ny, ic, jc, &i0, &j0);
      get_ij(momgrid, nx, ny, ic, jm, &i1, &j1);
      momgrid->ds_10_11[j][i] = 0.5 * (momgrid->h2[j0][i0] +
				       momgrid->h2[j1][i1]);
      get_ij(momgrid, nx, ny, ic, jc, &i0, &j0);
      momgrid->ds_01_21[j][i] = 2.0 * momgrid->h1[j0][i0];
      momgrid->ds_10_12[j][i] = 2.0 * momgrid->h2[j0][i0];
      
      momgrid->area[j][i] = momgrid->ds_01_21[j][i] * momgrid->ds_10_12[j][i];
      /* Angle clockwise between logical and geographic east */
      /*momgrid->angle[j][i] = momgrid->a1[jc][ic];*/
      /* Note : EMS angles in radians, MOM in degrees */
      momgrid->angle[j][i] = 180.0 * momgrid->a1[j0][i0] / PI;
    }
  }
}

/* END get_metrics()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Returns correct indices accounting for cyclic offsets             */
/*-------------------------------------------------------------------*/
void get_ij(momgrid_t *momgrid, int nx, int ny, int i, int j, int *ii, int *jj) 
{

  *ii = i;
  *jj = j;
  if (i >= 2*nx && momgrid->cytype & U1BDRY) {
    *ii = i - 2*nx;
    *jj = j + 2*momgrid->cyoffs;
  }
  if (j >= 2*ny && momgrid->cytype & U2BDRY) {
    *jj = j - 2*ny;
    *ii = i + 2*momgrid->cyoffs;
  }
  /* The coordinates leading to the offset should be NaN. This is */
  /* a fix to ensure coordinates are >= zero when the offset is   */
  /* added.                                                       */
  *ii = max(*ii, 0);
  *jj = max(*jj, 0);
}

/* END get_ij()                                                      */
/*-------------------------------------------------------------------*/



/*-------------------------------------------------------------------*/
/* Creates the metrics for a cell and writes to file                 */
/*-------------------------------------------------------------------*/
void get_metrics_old(momgrid_t *momgrid, int mode)
{
  int i, j, io, jo;
  int im, ic, ip;
  int jm, jc, jp;

  int ii, jj;

  int nx, ny;

  switch(mode) {
  case TCELL:
    io = 0;
    jo = 0;
    nx = momgrid->grid_x_T;
    ny = momgrid->grid_y_T;
    alloc_mom_mem(momgrid, nx, ny);
    break;
  case NCELL:
    io = 0;
    jo = 1;
    nx = momgrid->grid_x_T;
    ny = momgrid->grid_y_C;
    mom_free(momgrid);
    alloc_mom_mem(momgrid, nx, ny);
    break;
  case ECELL:
    io = 1;
    jo = 0;
    nx = momgrid->grid_x_C;
    ny = momgrid->grid_y_T;
    mom_free(momgrid);
    alloc_mom_mem(momgrid, nx, ny);
    break;
  case CCELL:
    io = 1;
    jo = 1;
    nx = momgrid->grid_x_C;
    ny = momgrid->grid_y_C;
    mom_free(momgrid);
    alloc_mom_mem(momgrid, nx, ny);
    break;
  default:
    io = 0;
    jo = 0;
    nx = 0;
    ny = 0;
  }

  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {

      im = io+i*2;
      ic = io+i*2+1;
      ip = io+i*2+2;

      jm = jo+j*2;
      jc = jo+j*2+1;
      jp = jo+j*2+2;

      /* Set cyclic offsets */
      if (im >= 2*nx && momgrid->cytype & U1BDRY) {
	im -= 2*nx;
	jm += 2*momgrid->cyoffs;
      }
      if (ic >= 2*nx && momgrid->cytype & U1BDRY) {
	ic -= 2*nx;
	jc += 2*momgrid->cyoffs;
      }
      if (ip >= 2*nx && momgrid->cytype & U1BDRY) {
	ip -= 2*nx;
	jp += 2*momgrid->cyoffs;
      }
      if (jm >= 2*ny && momgrid->cytype & U2BDRY) {
	jm -= 2*ny;
	im += 2*momgrid->cyoffs;
      }
      if (jc >= 2*ny && momgrid->cytype & U2BDRY) {
	jc -= 2*ny;
	ic += 2*momgrid->cyoffs;
      }
      if (jp >= 2*ny && momgrid->cytype & U2BDRY) {
	jp -= 2*ny;
	ip += 2*momgrid->cyoffs;
      }

      /* The coordinates leading to the offset should be NaN. This is */
      /* a fix to ensure coordinates are >= zero when the offset is   */
      /* added.                                                       */
      jm = max(jm, 0);
      jc = max(jc, 0);
      jp = max(jp, 0);
      im = max(im, 0);
      ic = max(ic, 0);
      ip = max(ip, 0);

      momgrid->x_g[j][i] = momgrid->x[jc][ic];
      momgrid->y_g[j][i] = momgrid->y[jc][ic];
      /* Verticies begins in SW corner, direction counterclockwise */
      momgrid->x_vert[0][j][i] = momgrid->x[jm][im];
      momgrid->x_vert[1][j][i] = momgrid->x[jm][ip];
      momgrid->x_vert[2][j][i] = momgrid->x[jp][ip];
      momgrid->x_vert[3][j][i] = momgrid->x[jp][im];
      momgrid->y_vert[0][j][i] = momgrid->y[jm][im];
      momgrid->y_vert[1][j][i] = momgrid->y[jm][ip];
      momgrid->y_vert[2][j][i] = momgrid->y[jp][ip];
      momgrid->y_vert[3][j][i] = momgrid->y[jp][im];
      /*
      momgrid->ds_00_02[j][i] = 2.0 * momgrid->h2[jc][im];
      momgrid->ds_20_22[j][i] = 2.0 * momgrid->h2[jc][ip];
      momgrid->ds_02_22[j][i] = 2.0 * momgrid->h1[jp][ic];
      momgrid->ds_00_20[j][i] = 2.0 * momgrid->h1[jm][ic];
      momgrid->ds_00_01[j][i] = momgrid->h2[jc][im];
      momgrid->ds_01_02[j][i] = momgrid->h2[jp][im];
      momgrid->ds_02_12[j][i] = momgrid->h1[jp][ic];
      momgrid->ds_12_22[j][i] = momgrid->h1[jp][ip];
      momgrid->ds_21_22[j][i] = momgrid->h2[jp][ip];
      momgrid->ds_20_21[j][i] = momgrid->h2[jc][ip];
      momgrid->ds_10_20[j][i] = momgrid->h1[jm][ip];
      momgrid->ds_00_10[j][i] = momgrid->h1[jm][ic];
      momgrid->ds_01_11[j][i] = momgrid->h1[jc][ic];
      momgrid->ds_11_12[j][i] = momgrid->h2[jp][ic];
      momgrid->ds_11_21[j][i] = momgrid->h1[jc][ip];
      momgrid->ds_10_11[j][i] = momgrid->h2[jc][ic];
      momgrid->ds_01_21[j][i] = 2.0 * momgrid->h1[jc][ic];
      momgrid->ds_10_12[j][i] = 2.0 * momgrid->h2[jc][ic];
      */
      momgrid->ds_00_02[j][i] = 2.0 * momgrid->h2[jc][im];
      momgrid->ds_20_22[j][i] = 2.0 * momgrid->h2[jc][ip];
      momgrid->ds_02_22[j][i] = 2.0 * momgrid->h1[jp][ic];
      momgrid->ds_00_20[j][i] = 2.0 * momgrid->h1[jm][ic];
      momgrid->ds_00_01[j][i] = 0.5 * (momgrid->h2[jc][im] +
				    momgrid->h2[jm][im]);
      momgrid->ds_01_02[j][i] = 0.5 * (momgrid->h2[jp][im] +
				    momgrid->h2[jc][im]);
      momgrid->ds_02_12[j][i] = 0.5 * (momgrid->h1[jp][ic] +
				    momgrid->h1[jp][im]);
      momgrid->ds_12_22[j][i] = 0.5 * (momgrid->h1[jp][ip] +
				    momgrid->h1[jp][ic]);
      momgrid->ds_21_22[j][i] = 0.5 * (momgrid->h2[jp][ip] +
				    momgrid->h2[jc][ip]);
      momgrid->ds_20_21[j][i] = 0.5 * (momgrid->h2[jc][ip] +
				    momgrid->h2[jm][ip]);
      momgrid->ds_10_20[j][i] = 0.5 * (momgrid->h1[jm][ip] +
				    momgrid->h1[jm][ic]);
      momgrid->ds_00_10[j][i] = 0.5 * (momgrid->h1[jm][ic] +
				    momgrid->h1[jm][im]);
      momgrid->ds_01_11[j][i] = 0.5 * (momgrid->h1[jc][ic] +
				    momgrid->h1[jc][im]);
      momgrid->ds_11_12[j][i] = 0.5 * (momgrid->h2[jp][ic] +
				    momgrid->h2[jc][ic]);
      momgrid->ds_11_21[j][i] = 0.5 * (momgrid->h1[jc][ip] +
				    momgrid->h1[jc][ic]);
      momgrid->ds_10_11[j][i] = 0.5 * (momgrid->h2[jc][ic] +
				    momgrid->h2[jm][ic]);
      momgrid->ds_01_21[j][i] = 2.0 * momgrid->h1[jc][ic];
      momgrid->ds_10_12[j][i] = 2.0 * momgrid->h2[jc][ic];
      momgrid->area[j][i] = momgrid->ds_01_21[j][i] * momgrid->ds_10_12[j][i];
      /* Angle clockwise between logical and geographic east */
      /*momgrid->angle[j][i] = momgrid->a1[jc][ic];*/
      /* Note : EMS angles in radians, MOM in degrees */
      momgrid->angle[j][i] = 180.0 * momgrid->a1[jc][ic] / PI;
    }
  }
}

/* END get_metrics_old()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Gets the nominal latitude and longitude vectors                   */
/*-------------------------------------------------------------------*/
void get_nom_ll(momgrid_t *momgrid)
{
  int i, j, p, m;
  float *a;
  double mnlat = 90.0;
  double mnlon = 180.0;
  double mxlat = -90.0;
  double mxlon = 0.0;
  double inc;
  int mode = 0;       /* mode<0 : consecutive integer labelling      */
                      /* mode=0 : uniform between min/max lat/lon    */
                      /* mode>0 : uses lat long along i=mode,j=mode  */

  if (mode > 0) {
    for (i = 0; i < momgrid->grid_x_T; i++)
      momgrid->grid_x_Ta[i] = momgrid->x[mode][i*2+1];
    for (j = 0; j < momgrid->grid_y_T; j++)
      momgrid->grid_y_Ta[j] = momgrid->y[j*2+1][mode];
    for (i = 0; i < momgrid->grid_x_C; i++)
      momgrid->grid_x_Ca[i] = momgrid->x[mode+1][i*2+2];
    for (j = 0; j < momgrid->grid_y_C; j++)
      momgrid->grid_y_Ca[j] = momgrid->y[j*2+2][mode+1];
  } else if (mode == 0) {
    for (i = 0; i < momgrid->grid_x_T; i++)
      for (j = 0; j < momgrid->grid_y_T; j++) {
	if (momgrid->x[j*2+1][i*2+1] > mxlon) 
	  mxlon = momgrid->x[j*2+1][i*2+1];
	if (momgrid->x[j*2+1][i*2+1] < mnlon) 
	  mnlon = momgrid->x[j*2+1][i*2+1];
      if (momgrid->y[j*2+1][i*2+1] > mxlat) 
	mxlat = momgrid->y[j*2+1][i*2+1];
      if (momgrid->y[j*2+1][i*2+1] < mnlat) 
	mnlat = momgrid->y[j*2+1][i*2+1];
      }
    inc = (mxlon - mnlon) / (double)(momgrid->grid_x_T - 1);
    for (i = 0; i < momgrid->grid_x_T; i++)
      momgrid->grid_x_Ta[i] = mnlon + i * inc;
    inc = (mxlat - mnlat) / (double)(momgrid->grid_y_T - 1);
    for (i = 0; i < momgrid->grid_y_T; i++)
      momgrid->grid_y_Ta[i] = mnlat + i * inc;
    
    mnlat = 90.0;
    mnlon = 180.0;
    mxlat = -90.0;
    mxlon = 0.0;
    for (i = 0; i < momgrid->grid_x_C; i++)
      for (j = 0; j < momgrid->grid_y_C; j++) {
	if (momgrid->x[j*2+2][i*2+2] > mxlon) 
	  mxlon = momgrid->x[j*2+2][i*2+2];
	if (momgrid->x[j*2+2][i*2+2] < mnlon) 
	  mnlon = momgrid->x[j*2+2][i*2+2];
	if (momgrid->y[j*2+2][i*2+2] > mxlat) 
	  mxlat = momgrid->y[j*2+2][i*2+2];
	if (momgrid->y[j*2+2][i*2+2] < mnlat) 
	  mnlat = momgrid->y[j*2+2][i*2+2];
      }
    inc = (mxlon - mnlon) / (double)(momgrid->grid_x_C - 1);
    for (i = 0; i < momgrid->grid_x_C; i++)
      momgrid->grid_x_Ca[i] = mnlon + i * inc;
    inc = (mxlat - mnlat) / (double)(momgrid->grid_y_C - 1);
    for (i = 0; i < momgrid->grid_y_C; i++)
      momgrid->grid_y_Ca[i] = mnlat + i * inc;
  } else if (mode < 0) {
    for (i = 0; i < momgrid->grid_x_T; i++) momgrid->grid_x_Ta[i] = i;
    for (j = 0; j < momgrid->grid_y_T; j++) momgrid->grid_y_Ta[j] = j;
    for (i = 0; i < momgrid->grid_x_C; i++) momgrid->grid_x_Ca[i] = i;
    for (j = 0; j < momgrid->grid_y_C; j++) momgrid->grid_y_Ca[j] = j;
  }

  a = momgrid->grid_x_Ta;
  for (i = 0; i < momgrid->grid_x_T; i++) {
    p = (i < momgrid->grid_x_T - 1) ? i + 1 : i;
    m = (i > 0) ? i - 1 : i;
    if((a[m] < a[i] && a[i] > a[p]) || (a[m] > a[i] && a[i] < a[p]))
      hd_warn("mom_grid: grid_x_T array NOT monotonously ordered at i = %d : %f %f %f\n", i, a[m], a[i], a[p]);
  }
  a = momgrid->grid_y_Ta;
  for (i = 0; i < momgrid->grid_y_T; i++) {
    p = (i < momgrid->grid_y_T - 1) ? i + 1 : i;
    m = (i > 0) ? i - 1 : i;
    if((a[m] < a[i] && a[i] > a[p]) || (a[m] > a[i] && a[i] < a[p]))
      hd_warn("mom_grid: grid_y_T array NOT monotonously ordered at i = %d : %f %f %f\n", i, a[m], a[i], a[p]);
  }
  a = momgrid->grid_x_Ca;
  for (i = 0; i < momgrid->grid_x_C; i++) {
    p = (i < momgrid->grid_x_C - 1) ? i + 1 : i;
    m = (i > 0) ? i - 1 : i;
    if((a[m] < a[i] && a[i] > a[p]) || (a[m] > a[i] && a[i] < a[p]))
      hd_warn("mom_grid: grid_x_C array NOT monotonously ordered at i = %d : %f %f %f\n", i, a[m], a[i], a[p]);
  }
  a = momgrid->grid_y_Ca;
  for (i = 0; i < momgrid->grid_y_C; i++) {
    p = (i < momgrid->grid_y_C - 1) ? i + 1 : i;
    m = (i > 0) ? i - 1 : i;
    if((a[m] < a[i] && a[i] > a[p]) || (a[m] > a[i] && a[i] < a[p]))
      hd_warn("mom_grid: grid_y_C array NOT monotonously ordered at i = %d : %f %f %f\n", i, a[m], a[i], a[p]);
  }
}

/* END get_nom_ll()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to allocate memory for the parameter data structure       */
/*-------------------------------------------------------------------*/
momgrid_t *momgrid_alloc(void)
{
  momgrid_t *momgrid = (momgrid_t *)malloc(sizeof(momgrid_t));
  memset(momgrid, 0, sizeof(momgrid_t));
  return momgrid;
}

/* END momgrid_alloc()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Allocates memory for grid metric arrays                           */
/*-------------------------------------------------------------------*/
void alloc_grid_mem(momgrid_t *momgrid, int nce1, int nce2)
{
  /* Allocate memory for the coordinate and metric arrays which are */
  /* at twice the final resolution.  */
  momgrid->x = d_alloc_2d(2 * nce1 + 1, 2 * nce2 + 1);
  momgrid->y = d_alloc_2d(2 * nce1 + 1, 2 * nce2 + 1);
  momgrid->h1 = d_alloc_2d(2 * nce1 + 1, 2 * nce2 + 1);
  momgrid->h2 = d_alloc_2d(2 * nce1 + 1, 2 * nce2 + 1);
  momgrid->a1 = d_alloc_2d(2 * nce1 + 1, 2 * nce2 + 1);
  momgrid->a2 = d_alloc_2d(2 * nce1 + 1, 2 * nce2 + 1);

  /* Allocate the MOM metric and position arrays */
  momgrid->x_T = d_alloc_2d(momgrid->grid_x_T, momgrid->grid_y_T);
  momgrid->y_T = d_alloc_2d(momgrid->grid_x_T, momgrid->grid_y_T);
  momgrid->x_C = d_alloc_2d(momgrid->grid_x_C, momgrid->grid_y_C);
  momgrid->y_C = d_alloc_2d(momgrid->grid_x_C, momgrid->grid_y_C);
  momgrid->depth = d_alloc_2d(momgrid->grid_x_T, momgrid->grid_y_T);
  momgrid->wet = d_alloc_2d(momgrid->grid_x_T, momgrid->grid_y_T);
  momgrid->num_levels = d_alloc_2d(momgrid->grid_x_T, momgrid->grid_y_T);
  momgrid->depth_c = d_alloc_2d(momgrid->grid_x_C, momgrid->grid_y_C);
  momgrid->wet_c = d_alloc_2d(momgrid->grid_x_C, momgrid->grid_y_C);
  momgrid->num_levels_c = d_alloc_2d(momgrid->grid_x_C, momgrid->grid_y_C);
  momgrid->grid_x_Ta = f_alloc_1d(momgrid->grid_x_T);
  momgrid->grid_y_Ta = f_alloc_1d(momgrid->grid_y_T);
  momgrid->grid_x_Ca = f_alloc_1d(momgrid->grid_x_C);
  momgrid->grid_y_Ca = f_alloc_1d(momgrid->grid_y_C);

  /* Allocate the MOM variable dummy arrays */
  momgrid->T2 = d_alloc_2d(momgrid->grid_x_T, momgrid->grid_y_T);
  momgrid->C2 = d_alloc_2d(momgrid->grid_x_C, momgrid->grid_y_C);

  /* Allocate open boundary arrays */
  if (momgrid->nobc) {
    momgrid->obcis = i_alloc_1d(momgrid->nobc);
    momgrid->obcjs = i_alloc_1d(momgrid->nobc);
    momgrid->obcie = i_alloc_1d(momgrid->nobc);
    momgrid->obcje = i_alloc_1d(momgrid->nobc);
  }
}

/* END alloc_grid_mem()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Allocate memory for mom metrics                                   */
/*-------------------------------------------------------------------*/
void alloc_mom_mem(momgrid_t *momgrid, int nce1, int nce2)
{
  /* Allocate the MOM metric and position arrays */
  momgrid->x_g = d_alloc_2d(nce1, nce2);
  momgrid->y_g = d_alloc_2d(nce1, nce2);
  momgrid->area = d_alloc_2d(nce1, nce2);
  momgrid->angle = d_alloc_2d(nce1, nce2);
  momgrid->ds_00_02 = d_alloc_2d(nce1, nce2);
  momgrid->ds_20_22 = d_alloc_2d(nce1, nce2);
  momgrid->ds_02_22 = d_alloc_2d(nce1, nce2);
  momgrid->ds_00_20 = d_alloc_2d(nce1, nce2);
  momgrid->ds_00_01 = d_alloc_2d(nce1, nce2);
  momgrid->ds_01_02 = d_alloc_2d(nce1, nce2);
  momgrid->ds_02_12 = d_alloc_2d(nce1, nce2);
  momgrid->ds_12_22 = d_alloc_2d(nce1, nce2);
  momgrid->ds_21_22 = d_alloc_2d(nce1, nce2);
  momgrid->ds_20_21 = d_alloc_2d(nce1, nce2);
  momgrid->ds_10_20 = d_alloc_2d(nce1, nce2);
  momgrid->ds_00_10 = d_alloc_2d(nce1, nce2);
  momgrid->ds_01_11 = d_alloc_2d(nce1, nce2);
  momgrid->ds_11_12 = d_alloc_2d(nce1, nce2);
  momgrid->ds_11_21 = d_alloc_2d(nce1, nce2);
  momgrid->ds_10_11 = d_alloc_2d(nce1, nce2);
  momgrid->ds_01_21 = d_alloc_2d(nce1, nce2);
  momgrid->ds_10_12 = d_alloc_2d(nce1, nce2);
  momgrid->x_vert = d_alloc_3d(nce1, nce2, 4);
  momgrid->y_vert = d_alloc_3d(nce1, nce2, 4);
}

/* END alloc_mom_mem()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Frees the grid structure                                          */
/*-------------------------------------------------------------------*/
void mom_grid_free(momgrid_t *momgrid, int mode)
{
  if (mode & UNUSED) {
    d_free_2d(momgrid->x);
    d_free_2d(momgrid->y);
    d_free_2d(momgrid->h1);
    d_free_2d(momgrid->h2);
    d_free_2d(momgrid->a1);
    d_free_2d(momgrid->a2);

    d_free_1d(momgrid->layers);
    d_free_2d(momgrid->depth);
    d_free_2d(momgrid->wet);
    d_free_2d(momgrid->num_levels);
    d_free_2d(momgrid->depth_c);
    d_free_2d(momgrid->wet_c);
    d_free_2d(momgrid->num_levels_c);
    if (momgrid->topo != NULL)
      d_free_2d(momgrid->topo);
    if (momgrid->nobc) {
      i_free_1d(momgrid->obcis);
      i_free_1d(momgrid->obcjs);
      i_free_1d(momgrid->obcie);
      i_free_1d(momgrid->obcje);
    }
  }
  if (mode & ALL) {
    f_free_1d(momgrid->grid_x_Ta);
    f_free_1d(momgrid->grid_y_Ta);
    f_free_1d(momgrid->grid_x_Ca);
    f_free_1d(momgrid->grid_y_Ca);
    d_free_2d(momgrid->x_T);
    d_free_2d(momgrid->y_T);
    d_free_2d(momgrid->x_C);
    d_free_2d(momgrid->y_C);
    f_free_1d(momgrid->zt);
    f_free_1d(momgrid->zb);
    d_free_2d(momgrid->T2);
    d_free_3d(momgrid->T3);
    d_free_2d(momgrid->C2);
    d_free_3d(momgrid->C3);
    free((momgrid_t *)momgrid);
  }
}

/* END mom_grid_free()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Frees the mom metric arrays                                       */
/*-------------------------------------------------------------------*/
void mom_free(momgrid_t *momgrid)
{
  d_free_2d(momgrid->x_g);
  d_free_2d(momgrid->y_g);
  d_free_3d(momgrid->x_vert);
  d_free_3d(momgrid->y_vert);
  d_free_2d(momgrid->area);
  d_free_2d(momgrid->angle);
  d_free_2d(momgrid->ds_00_02);
  d_free_2d(momgrid->ds_20_22);
  d_free_2d(momgrid->ds_02_22);
  d_free_2d(momgrid->ds_00_20);
  d_free_2d(momgrid->ds_00_01);
  d_free_2d(momgrid->ds_01_02);
  d_free_2d(momgrid->ds_02_12);
  d_free_2d(momgrid->ds_12_22);
  d_free_2d(momgrid->ds_21_22);
  d_free_2d(momgrid->ds_20_21);
  d_free_2d(momgrid->ds_10_20);
  d_free_2d(momgrid->ds_00_10);
  d_free_2d(momgrid->ds_01_11);
  d_free_2d(momgrid->ds_11_12);
  d_free_2d(momgrid->ds_11_21);
  d_free_2d(momgrid->ds_10_11);
  d_free_2d(momgrid->ds_01_21);
  d_free_2d(momgrid->ds_10_12);
}

/* END mom_free()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads bathymetry from a MOM grid_spec.nc file                     */
/*-------------------------------------------------------------------*/
double *read_bathy_from_mom(char *fname, int nx, int ny)
{
  int fid, n;
  int i, j;
  size_t start[4];
  size_t count[4];
  size_t grid_x_T, grid_y_T;
  size_t grid_x_C, grid_y_C;
  int nce1, nce2;
  int gs;
  double **d1, **d2;
  double *bathy = NULL;
  int ncerr;

  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;
  /* Time independent variables */
  count[0] = 0;
  count[1] = 0;
  count[2] = 0;
  count[3] = 0;

  /*-----------------------------------------------------------------*/
  /* Open the dump file for reading                                  */
  if ((ncerr = nc_open(fname, NC_NOWRITE, &fid)) != NC_NOERR) {
    printf("Can't find input file %s\n", fname);
    hd_quit((char *)nc_strerror(ncerr));
  }

  /*-----------------------------------------------------------------*/
  /* Get MOM grid dimensions                                         */
  gs = (gsize & TCELL) ? 0 : 1;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "grid_x_T"), &grid_x_T);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "grid_y_T"), &grid_y_T);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "grid_x_C"), &grid_x_C);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "grid_y_C"), &grid_y_C);
  nce1 = grid_x_T + gs;
  nce2 = grid_y_T + gs;
  if (nce1 != nx || nce2 != ny)
    hd_quit("Inconsistent grid sizes EMS (%dx%d) : MOM (%dx%d)\n",
	    nx, ny, nce1, nce2);

  /*-----------------------------------------------------------------*/
  /* Convert the bathymetry                                          */
  bathy = d_alloc_1d(nce1 * nce2);
  count[0] = grid_y_T;
  count[1] = grid_x_T;
  d1 = d_alloc_2d(count[1], count[0]);
  d2 = d_alloc_2d(count[1], count[0]);
  nc_get_vara_double(fid, ncw_var_id(fid, "depth_t"), start, count, d1[0]);
  nc_get_vara_double(fid, ncw_var_id(fid, "wet"), start, count, d2[0]);
  n = 0;
  for (j = 0; j < grid_y_T; j++)
    for (i = 0; i < grid_x_T; i++) {
      bathy[n] = d1[j][i];
      if (d2[j][i] == 0.0)
	bathy[n] = -LANDCELL;
      n++;
      if (gs && i == grid_x_T - 1) {
	bathy[n++] = d1[j][i];
      }
    }
  if (gs) {
    j = grid_y_T - 1;
    for (i = 0; i < grid_x_T; i++) {
      bathy[n] = d1[j][i];
      if (d2[j][i] == 0.0) bathy[n] = -LANDCELL;
      n++;
    }
    bathy[n] = d1[grid_y_T-1][grid_x_T-1];
  }
  d_free_2d(d1);
  d_free_2d(d2);
  ncw_close(fname, fid);
  return(bathy);
}

/* END read_bathy_from_mom()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Creates the MOM output netCDF file                                */
/*-------------------------------------------------------------------*/
int dump_create_mom(momgrid_t *momgrid, char *name, int mode)
{
  int cdfid;
  int dims[10];
  int vid;
  /* dimension ids */
  int grid_x_Tid;               /* I dimension id at grid centre */
  int grid_y_Tid;               /* J dimension id at grid centre */
  int grid_x_Cid;               /* I dimension id at grid corner */
  int grid_y_Cid;               /* J dimension id at grid corner */
  int vertexid;                 /* id for vertices               */
  int ztid;                     /* K dimension id at grid center */
  int zbid;                     /* K dimension id at grid face   */

  /*-----------------------------------------------------------------*/
  /* Create the netCDF file                                          */
  if (ncw_create(name,(overwrite_output()?NC_CLOBBER:NC_NOCLOBBER), &cdfid) != NC_NOERR)
    hd_quit("Couldn't create output dump file %s\n", name);

  /*-----------------------------------------------------------------*/
  /* Define dimensions                                               */
  nc_def_dim(cdfid, "grid_x_T", momgrid->grid_x_T, &grid_x_Tid);
  nc_def_dim(cdfid, "grid_y_T", momgrid->grid_y_T, &grid_y_Tid);
  if (mode & METRIC) {
    nc_def_dim(cdfid, "grid_x_C", momgrid->grid_x_C, &grid_x_Cid);
    nc_def_dim(cdfid, "grid_y_C", momgrid->grid_y_C, &grid_y_Cid);
    nc_def_dim(cdfid, "vertex", 4, &vertexid);
  }
  nc_def_dim(cdfid, "zt", momgrid->nz, &ztid);
  if (mode & METRIC)
    nc_def_dim(cdfid, "zb", momgrid->nz, &zbid);

  /* time independent variables */
  dims[0] = grid_x_Tid;
  nc_def_var(cdfid, "grid_x_T", NC_FLOAT, 1, dims, &vid);
  dims[0] = grid_y_Tid;
  nc_def_var(cdfid, "grid_y_T", NC_FLOAT, 1, dims, &vid);
  if (mode & METRIC) {
    dims[0] = grid_x_Cid;
    nc_def_var(cdfid, "grid_x_C", NC_FLOAT, 1, dims, &vid);
    dims[0] = grid_y_Cid;
    nc_def_var(cdfid, "grid_y_C", NC_FLOAT, 1, dims, &vid);
    dims[0] = vertexid;
    nc_def_var(cdfid, "vertex", NC_FLOAT, 1, dims, &vid);
  }

  /* T grid cells */
  dims[0] = grid_y_Tid;
  dims[1] = grid_x_Tid;
  nc_def_var(cdfid, "x_T", NC_DOUBLE, 2, dims, &vid);
  nc_def_var(cdfid, "y_T", NC_DOUBLE, 2, dims, &vid);
  if (mode & METRIC) {
    dims[0] = vertexid;
    dims[1] = grid_y_Tid;
    dims[2] = grid_x_Tid;
    nc_def_var(cdfid, "x_vert_T", NC_DOUBLE, 3, dims, &vid);
    nc_def_var(cdfid, "y_vert_T", NC_DOUBLE, 3, dims, &vid);
    dims[0] = grid_y_Tid;
    dims[1] = grid_x_Tid;
    nc_def_var(cdfid, "area_T", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "angle_T", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_00_02_T", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_20_22_T", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_02_22_T", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_00_20_T", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_00_01_T", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_01_02_T", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_02_12_T", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_12_22_T", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_21_22_T", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_20_21_T", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_10_20_T", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_00_10_T", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_01_11_T", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_11_12_T", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_11_21_T", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_10_11_T", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_01_21_T", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_10_12_T", NC_DOUBLE, 2, dims, &vid);

    /* E grid cells */
    dims[0] = grid_y_Tid;
    dims[1] = grid_x_Cid;
    nc_def_var(cdfid, "x_E", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "y_E", NC_DOUBLE, 2, dims, &vid);
    dims[0] = vertexid;
    dims[1] = grid_y_Tid;
    dims[2] = grid_x_Cid;
    nc_def_var(cdfid, "x_vert_E", NC_DOUBLE, 3, dims, &vid);
    nc_def_var(cdfid, "y_vert_E", NC_DOUBLE, 3, dims, &vid);
    dims[0] = grid_y_Tid;
    dims[1] = grid_x_Cid;
    nc_def_var(cdfid, "area_E", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "angle_E", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_00_02_E", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_20_22_E", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_02_22_E", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_00_20_E", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_00_01_E", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_01_02_E", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_02_12_E", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_12_22_E", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_21_22_E", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_20_21_E", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_10_20_E", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_00_10_E", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_01_11_E", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_11_12_E", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_11_21_E", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_10_11_E", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_01_21_E", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_10_12_E", NC_DOUBLE, 2, dims, &vid);

    /* N grid cells */
    dims[0] = grid_y_Cid;
    dims[1] = grid_x_Tid;
    nc_def_var(cdfid, "x_N", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "y_N", NC_DOUBLE, 2, dims, &vid);
    dims[0] = vertexid;
    dims[1] = grid_y_Cid;
    dims[2] = grid_x_Tid;
    nc_def_var(cdfid, "x_vert_N", NC_DOUBLE, 3, dims, &vid);
    nc_def_var(cdfid, "y_vert_N", NC_DOUBLE, 3, dims, &vid);
    dims[0] = grid_y_Cid;
    dims[1] = grid_x_Tid;
    nc_def_var(cdfid, "area_N", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "angle_N", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_00_02_N", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_20_22_N", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_02_22_N", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_00_20_N", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_00_01_N", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_01_02_N", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_02_12_N", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_12_22_N", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_21_22_N", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_20_21_N", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_10_20_N", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_00_10_N", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_01_11_N", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_11_12_N", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_11_21_N", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_10_11_N", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_01_21_N", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_10_12_N", NC_DOUBLE, 2, dims, &vid);
    
    /* C grid cells */
    dims[0] = grid_y_Cid;
    dims[1] = grid_x_Cid;
    nc_def_var(cdfid, "x_C", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "y_C", NC_DOUBLE, 2, dims, &vid);
    dims[0] = vertexid;
    dims[1] = grid_y_Cid;
    dims[2] = grid_x_Cid;
    nc_def_var(cdfid, "x_vert_C", NC_DOUBLE, 3, dims, &vid);
    nc_def_var(cdfid, "y_vert_C", NC_DOUBLE, 3, dims, &vid);
    dims[0] = grid_y_Cid;
    dims[1] = grid_x_Cid;
    nc_def_var(cdfid, "area_C", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "angle_C", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_00_02_C", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_20_22_C", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_02_22_C", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_00_20_C", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_00_01_C", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_01_02_C", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_02_12_C", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_12_22_C", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_21_22_C", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_20_21_C", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_10_20_C", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_00_10_C", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_01_11_C", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_11_12_C", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_11_21_C", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_10_11_C", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_01_21_C", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "ds_10_12_C", NC_DOUBLE, 2, dims, &vid);
  }

  /* Vertical grid */
  dims[0] = ztid;
  nc_def_var(cdfid, "zt", NC_FLOAT, 1, dims, &vid);
  if (mode & METRIC) {
    dims[0] = zbid;
    nc_def_var(cdfid, "zb", NC_FLOAT, 1, dims, &vid);

    /* Bathymetry */
    dims[0] = grid_y_Tid;
    dims[1] = grid_x_Tid;
    nc_def_var(cdfid, "depth_t", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "num_levels", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "wet", NC_DOUBLE, 2, dims, &vid);
    // Don't write out c-cell stuff, MOM4p1 will generate this if not
    // present but the MARVL version doesn't even work if present
    // [MOMref1]
    /*
      dims[0] = grid_y_Cid;
      dims[1] = grid_x_Cid;
      nc_def_var(cdfid, "depth_c", NC_DOUBLE, 2, dims, &vid);
      nc_def_var(cdfid, "num_levels_c", NC_DOUBLE, 2, dims, &vid);
      nc_def_var(cdfid, "wet_c", NC_DOUBLE, 2, dims, &vid);
    */
  }

  /* Temperature and salinity                                        */
  if (mode & TEMSAL) {
    dims[0] = ztid;
    dims[1] = grid_y_Tid;
    dims[2] = grid_x_Tid;
    nc_def_var(cdfid, "temp", NC_DOUBLE, 3, dims, &vid);
    nc_def_var(cdfid, "salt", NC_DOUBLE, 3, dims, &vid);
  }

  /* Sea level                                                       */
  if (mode & ETA_A) {
    dims[0] = grid_y_Tid;
    dims[1] = grid_x_Tid;
    nc_def_var(cdfid, "eta_t", NC_DOUBLE, 2, dims, &vid);
  }

  write_attributes(momgrid, cdfid, mode);

  nc_enddef(cdfid);
  return (cdfid);
}

/* END dump_create_mom()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Writes MOM netCDF attributes                                      */
/*-------------------------------------------------------------------*/
void write_attributes(momgrid_t *momgrid, int cdfid, int mode)
{
  /* dimension ids */
  int n, vid;
  char buf[MAXSTRLEN];
  char attname[MAXSTRLEN];

  /* time independent variables */
  vid = ncw_var_id(cdfid, "grid_x_T");
  write_text_att(cdfid, vid, "long_name", "Nominal Longitude of T-cell center");
  write_text_att(cdfid, vid, "units", "degree_east");
  write_text_att(cdfid, vid, "cartesian_axis", "X");
  vid = ncw_var_id(cdfid, "grid_y_T");
  write_text_att(cdfid, vid, "long_name", "Nominal Latitude of T-cell center");
  write_text_att(cdfid, vid, "units", "degree_north");
  write_text_att(cdfid, vid, "cartesian_axis", "Y");
  if (mode & METRIC) {
    vid = ncw_var_id(cdfid, "grid_x_C");
    write_text_att(cdfid, vid, "long_name", "Nominal Longitude of C-cell center");
    write_text_att(cdfid, vid, "units", "degree_east");
    write_text_att(cdfid, vid, "cartesian_axis", "X");
    vid = ncw_var_id(cdfid, "grid_y_C");
    write_text_att(cdfid, vid, "long_name", "Nominal Latitude of C-cell center");
    write_text_att(cdfid, vid, "units", "degree_north");
    write_text_att(cdfid, vid, "cartesian_axis", "Y");
    vid = ncw_var_id(cdfid, "vertex");
    write_text_att(cdfid, vid, "long_name", "Vertex position from southwest couterclockwise");
    write_text_att(cdfid, vid, "units", "none");
  }
  /* T grid cells */
  vid = ncw_var_id(cdfid, "x_T");
  write_text_att(cdfid, vid, "long_name", "Geographic longitude of T_cell centers");
  write_text_att(cdfid, vid, "units", "degree_east");
  vid = ncw_var_id(cdfid, "y_T");
  write_text_att(cdfid, vid, "long_name", "Geographic latitude of T_cell centers");
  write_text_att(cdfid, vid, "units", "degree_north");
  if (mode & METRIC) {
    vid = ncw_var_id(cdfid, "x_vert_T");
    write_text_att(cdfid, vid, "long_name", "Geographic longitude of T_cell vertices begin southwest counterclockwise");
    write_text_att(cdfid, vid, "units", "degree_east");
    vid = ncw_var_id(cdfid, "y_vert_T");
    write_text_att(cdfid, vid, "long_name", "Geographic latitude of T_cell vertices begin southwest counterclockwise");
    write_text_att(cdfid, vid, "units", "degree_north");
    vid = ncw_var_id(cdfid, "area_T");
    write_text_att(cdfid, vid, "long_name", "Area of T_cell");
    write_text_att(cdfid, vid, "units", "m2");
    vid = ncw_var_id(cdfid, "angle_T");
    write_text_att(cdfid, vid, "long_name", "Angle clockwise between logical and geographic east of T_cell");
    write_text_att(cdfid, vid, "units", "degree");
    vid = ncw_var_id(cdfid, "ds_00_02_T");
    write_text_att(cdfid, vid, "long_name", "Length of western face of T_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_20_22_T");
    write_text_att(cdfid, vid, "long_name", "Length of eastern face of T_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_02_22_T");
    write_text_att(cdfid, vid, "long_name", "Length of northern face of T_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_00_20_T");
    write_text_att(cdfid, vid, "long_name", "Length of southern face of T_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_00_01_T");
    write_text_att(cdfid, vid, "long_name", "Distance from southwest corner to western face center of T_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_01_02_T");
    write_text_att(cdfid, vid, "long_name", "Distance from northwest corner to western face center of T_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_02_12_T");
    write_text_att(cdfid, vid, "long_name", "Distance from northwest corner to northern face center of T_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_12_22_T");
    write_text_att(cdfid, vid, "long_name", "Distance from northeast corner to northern face center of T_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_21_22_T");
    write_text_att(cdfid, vid, "long_name", "Distance from northeast corner to eastern face center of T_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_20_21_T");
    write_text_att(cdfid, vid, "long_name", "Distance from southeast corner to eastern face center of T_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_10_20_T");
    write_text_att(cdfid, vid, "long_name", "Distance from southeast corner to southern face center of T_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_00_10_T");
    write_text_att(cdfid, vid, "long_name", "Distance from southwest corner to southern face center of T_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_01_11_T");
    write_text_att(cdfid, vid, "long_name", "Distance from center to western face of T_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_11_12_T");
    write_text_att(cdfid, vid, "long_name", "Distance from center to northern face of T_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_11_21_T");
    write_text_att(cdfid, vid, "long_name", "Distance from center to eastern face of T_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_10_11_T");
    write_text_att(cdfid, vid, "long_name", "Distance from center to southern face of T_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_01_21_T");
    write_text_att(cdfid, vid, "long_name", "width of T_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_10_12_T");
    write_text_att(cdfid, vid, "long_name", "height of T_cell");
    write_text_att(cdfid, vid, "units", "m");

    /* E grid cells */
    vid = ncw_var_id(cdfid, "x_E");
    write_text_att(cdfid, vid, "long_name", "Geographic longitude of E_cell centers");
    write_text_att(cdfid, vid, "units", "degree_east");
    vid = ncw_var_id(cdfid, "y_E");
    write_text_att(cdfid, vid, "long_name", "Geographic latitude of E_cell centers");
    write_text_att(cdfid, vid, "units", "degree_north");
    vid = ncw_var_id(cdfid, "x_vert_E");
    write_text_att(cdfid, vid, "long_name", "Geographic longitude of E_cell vertices begin southwest counterclockwise");
    write_text_att(cdfid, vid, "units", "degree_east");
    vid = ncw_var_id(cdfid, "y_vert_E");
    write_text_att(cdfid, vid, "long_name", "Geographic latitude of E_cell vertices begin southwest counterclockwise");
    write_text_att(cdfid, vid, "units", "degree_north");
    vid = ncw_var_id(cdfid, "area_E");
    write_text_att(cdfid, vid, "long_name", "Area of E_cell");
    write_text_att(cdfid, vid, "units", "m2");
    vid = ncw_var_id(cdfid, "angle_E");
    write_text_att(cdfid, vid, "long_name", "Angle clockwise between logical and geographic east of E_cell");
    write_text_att(cdfid, vid, "units", "degree");
    vid = ncw_var_id(cdfid, "ds_00_02_E");
    write_text_att(cdfid, vid, "long_name", "Length of western face of E_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_20_22_E");
    write_text_att(cdfid, vid, "long_name", "Length of eastern face of E_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_02_22_E");
    write_text_att(cdfid, vid, "long_name", "Length of northern face of E_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_00_20_E");
    write_text_att(cdfid, vid, "long_name", "Length of southern face of E_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_00_01_E");
    write_text_att(cdfid, vid, "long_name", "Distance from southwest corner to western face center of E_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_01_02_E");
    write_text_att(cdfid, vid, "long_name", "Distance from northwest corner to western face center of E_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_02_12_E");
    write_text_att(cdfid, vid, "long_name", "Distance from northwest corner to northern face center of E_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_12_22_E");
    write_text_att(cdfid, vid, "long_name", "Distance from northeast corner to northern face center of E_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_21_22_E");
    write_text_att(cdfid, vid, "long_name", "Distance from northeast corner to eastern face center of E_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_20_21_E");
    write_text_att(cdfid, vid, "long_name", "Distance from southeast corner to eastern face center of E_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_10_20_E");
    write_text_att(cdfid, vid, "long_name", "Distance from southeast corner to southern face center of E_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_00_10_E");
    write_text_att(cdfid, vid, "long_name", "Distance from southwest corner to southern face center of E_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_01_11_E");
    write_text_att(cdfid, vid, "long_name", "Distance from center to western face of E_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_11_12_E");
    write_text_att(cdfid, vid, "long_name", "Distance from center to northern face of E_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_11_21_E");
    write_text_att(cdfid, vid, "long_name", "Distance from center to eastern face of E_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_10_11_E");
    write_text_att(cdfid, vid, "long_name", "Distance from center to southern face of E_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_01_21_E");
    write_text_att(cdfid, vid, "long_name", "width of E_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_10_12_E");
    write_text_att(cdfid, vid, "long_name", "height of E_cell");
    write_text_att(cdfid, vid, "units", "m");
    
    /* N grid cells */
    vid = ncw_var_id(cdfid, "x_N");
    write_text_att(cdfid, vid, "long_name", "Geographic longitude of N_cell centers");
    write_text_att(cdfid, vid, "units", "degree_east");
    vid = ncw_var_id(cdfid, "y_N");
    write_text_att(cdfid, vid, "long_name", "Geographic latitude of N_cell centers");
    write_text_att(cdfid, vid, "units", "degree_north");
    vid = ncw_var_id(cdfid, "x_vert_N");
    write_text_att(cdfid, vid, "long_name", "Geographic longitude of N_cell vertices begin southwest counterclockwise");
    write_text_att(cdfid, vid, "units", "degree_east");
    vid = ncw_var_id(cdfid, "y_vert_N");
    write_text_att(cdfid, vid, "long_name", "Geographic latitude of N_cell vertices begin southwest counterclockwise");
    write_text_att(cdfid, vid, "units", "degree_north");
    vid = ncw_var_id(cdfid, "area_N");
    write_text_att(cdfid, vid, "long_name", "Area of N_cell");
    write_text_att(cdfid, vid, "units", "m2");
    vid = ncw_var_id(cdfid, "angle_N");
    write_text_att(cdfid, vid, "long_name", "Angle clockwise between logical and geographic east of N_cell");
    write_text_att(cdfid, vid, "units", "degree");
    vid = ncw_var_id(cdfid, "ds_00_02_N");
    write_text_att(cdfid, vid, "long_name", "Length of western face of N_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_20_22_N");
    write_text_att(cdfid, vid, "long_name", "Length of eastern face of N_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_02_22_N");
    write_text_att(cdfid, vid, "long_name", "Length of northern face of N_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_00_20_N");
    write_text_att(cdfid, vid, "long_name", "Length of southern face of N_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_00_01_N");
    write_text_att(cdfid, vid, "long_name", "Distance from southwest corner to western face center of N_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_01_02_N");
    write_text_att(cdfid, vid, "long_name", "Distance from northwest corner to western face center of N_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_02_12_N");
    write_text_att(cdfid, vid, "long_name", "Distance from northwest corner to northern face center of N_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_12_22_N");
    write_text_att(cdfid, vid, "long_name", "Distance from northeast corner to northern face center of N_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_21_22_N");
    write_text_att(cdfid, vid, "long_name", "Distance from northeast corner to eastern face center of N_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_20_21_N");
    write_text_att(cdfid, vid, "long_name", "Distance from southeast corner to eastern face center of N_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_10_20_N");
    write_text_att(cdfid, vid, "long_name", "Distance from southeast corner to southern face center of N_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_00_10_N");
    write_text_att(cdfid, vid, "long_name", "Distance from southwest corner to southern face center of N_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_01_11_N");
    write_text_att(cdfid, vid, "long_name", "Distance from center to western face of N_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_11_12_N");
    write_text_att(cdfid, vid, "long_name", "Distance from center to northern face of N_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_11_21_N");
    write_text_att(cdfid, vid, "long_name", "Distance from center to eastern face of N_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_10_11_N");
    write_text_att(cdfid, vid, "long_name", "Distance from center to southern face of N_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_01_21_N");
    write_text_att(cdfid, vid, "long_name", "width of N_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_10_12_N");
    write_text_att(cdfid, vid, "long_name", "height of N_cell");
    write_text_att(cdfid, vid, "units", "m");

    /* C grid cells */
    vid = ncw_var_id(cdfid, "x_C");
    write_text_att(cdfid, vid, "long_name", "Geographic longitude of C_cell centers");
    write_text_att(cdfid, vid, "units", "degree_east");
    vid = ncw_var_id(cdfid, "y_C");
    write_text_att(cdfid, vid, "long_name", "Geographic latitude of C_cell centers");
    write_text_att(cdfid, vid, "units", "degree_north");
    vid = ncw_var_id(cdfid, "x_vert_C");
    write_text_att(cdfid, vid, "long_name", "Geographic longitude of C_cell vertices begin southwest counterclockwise");
    write_text_att(cdfid, vid, "units", "degree_east");
    vid = ncw_var_id(cdfid, "y_vert_C");
    write_text_att(cdfid, vid, "long_name", "Geographic latitude of C_cell vertices begin southwest counterclockwise");
    write_text_att(cdfid, vid, "units", "degree_north");
    vid = ncw_var_id(cdfid, "area_C");
    write_text_att(cdfid, vid, "long_name", "Area of C_cell");
    write_text_att(cdfid, vid, "units", "m2");
    vid = ncw_var_id(cdfid, "angle_C");
    write_text_att(cdfid, vid, "long_name", "Angle clockwise between logical and geographic east of C_cell");
    write_text_att(cdfid, vid, "units", "degree");
    vid = ncw_var_id(cdfid, "ds_00_02_C");
    write_text_att(cdfid, vid, "long_name", "Length of western face of C_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_20_22_C");
    write_text_att(cdfid, vid, "long_name", "Length of eastern face of C_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_02_22_C");
    write_text_att(cdfid, vid, "long_name", "Length of northern face of C_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_00_20_C");
    write_text_att(cdfid, vid, "long_name", "Length of southern face of C_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_00_01_C");
    write_text_att(cdfid, vid, "long_name", "Distance from southwest corner to western face center of C_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_01_02_C");
    write_text_att(cdfid, vid, "long_name", "Distance from northwest corner to western face center of C_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_02_12_C");
    write_text_att(cdfid, vid, "long_name", "Distance from northwest corner to northern face center of C_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_12_22_C");
    write_text_att(cdfid, vid, "long_name", "Distance from northeast corner to northern face center of C_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_21_22_C");
    write_text_att(cdfid, vid, "long_name", "Distance from northeast corner to eastern face center of C_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_20_21_C");
    write_text_att(cdfid, vid, "long_name", "Distance from southeast corner to eastern face center of C_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_10_20_C");
    write_text_att(cdfid, vid, "long_name", "Distance from southeast corner to southern face center of C_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_00_10_C");
    write_text_att(cdfid, vid, "long_name", "Distance from southwest corner to southern face center of C_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_01_11_C");
    write_text_att(cdfid, vid, "long_name", "Distance from center to western face of C_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_11_12_C");
    write_text_att(cdfid, vid, "long_name", "Distance from center to northern face of C_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_11_21_C");
    write_text_att(cdfid, vid, "long_name", "Distance from center to eastern face of C_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_10_11_C");
    write_text_att(cdfid, vid, "long_name", "Distance from center to southern face of C_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_01_21_C");
    write_text_att(cdfid, vid, "long_name", "width of C_cell");
    write_text_att(cdfid, vid, "units", "m");
    vid = ncw_var_id(cdfid, "ds_10_12_C");
    write_text_att(cdfid, vid, "long_name", "height of C_cell");
    write_text_att(cdfid, vid, "units", "m");
  }

  /* Vertical grid */
  vid = ncw_var_id(cdfid, "zt");
  write_text_att(cdfid, vid, "long_name", "zt");
  write_text_att(cdfid, vid, "units", "meters");
  write_text_att(cdfid, vid, "cartesian_axis", "z");
  write_text_att(cdfid, vid, "positive", "down");
  if (mode & METRIC) {
    vid = ncw_var_id(cdfid, "zb");
    write_text_att(cdfid, vid, "long_name", "zb");
    write_text_att(cdfid, vid, "units", "meters");
    write_text_att(cdfid, vid, "cartesian_axis", "z");
    write_text_att(cdfid, vid, "positive", "down");
    
    /* Bathymetry */
    vid = ncw_var_id(cdfid, "depth_t");
    write_text_att(cdfid, vid, "long_name", "topographic depth of T-cell");
    write_text_att(cdfid, vid, "units", "meters");
    vid = ncw_var_id(cdfid, "num_levels");
    write_text_att(cdfid, vid, "long_name", "number of vertical T-cells");
    write_text_att(cdfid, vid, "units", "none");
    vid = ncw_var_id(cdfid, "wet");
    write_text_att(cdfid, vid, "long_name", "land/sea flag (0=land) for T-cell");
    write_text_att(cdfid, vid, "units", "none");

    // See [MOMref1] comments
    /*
      vid = ncw_var_id(cdfid, "depth_c");
      write_text_att(cdfid, vid, "long_name", "topographic depth of C-cell");
      write_text_att(cdfid, vid, "units", "meters");
      vid = ncw_var_id(cdfid, "num_levels_c");
      write_text_att(cdfid, vid, "long_name", "number of vertical C-cells");
      write_text_att(cdfid, vid, "units", "none");
      vid = ncw_var_id(cdfid, "wet_c");
      write_text_att(cdfid, vid, "long_name", "land/sea flag (0=land) for C-cell");
      write_text_att(cdfid, vid, "units", "none");
    */
  }

  /* Temperature and salinity                                        */
  if (mode & TEMSAL) {
    vid = ncw_var_id(cdfid, "temp");
    write_text_att(cdfid, vid, "long_name", "initial potential temp");
    write_text_att(cdfid, vid, "units", "deg_c");
    vid = ncw_var_id(cdfid, "salt");
    write_text_att(cdfid, vid, "long_name", "salinity");
    write_text_att(cdfid, vid, "units", "psu");
  }

  /* Global attributes */
  write_text_att(cdfid, NC_GLOBAL, "filename", momgrid->ofile);
  if (mode & METRIC) {
    write_text_att(cdfid, NC_GLOBAL, "xname", "longitude");
    write_text_att(cdfid, NC_GLOBAL, "yname", "latitude");
    write_text_att(cdfid, NC_GLOBAL, "vertex_convention", "SWCCW");
    if (momgrid->cytype & U2BDRY)
      write_text_att(cdfid, NC_GLOBAL, "y_boundary_type", "cyclic");
    else
      write_text_att(cdfid, NC_GLOBAL, "y_boundary_type", "solid_walls");
    if (momgrid->cytype & U1BDRY)
      write_text_att(cdfid, NC_GLOBAL, "x_boundary_type", "cyclic");
    else
      write_text_att(cdfid, NC_GLOBAL, "x_boundary_type", "solid_walls");
    /*write_text_att(cdfid, NC_GLOBAL, "beta_plane", "y");*/
    write_text_att(cdfid, NC_GLOBAL, "topography", "rectangular_basin");
    write_text_att(cdfid, NC_GLOBAL, "input_file", "/archive/fms/mom4/input_data/OCCAM_p5degree.nc");
    write_text_att(cdfid, NC_GLOBAL, "input_field", "topo");
    write_text_att(cdfid, NC_GLOBAL, "fill_isolated_cells", "y");
    write_text_att(cdfid, NC_GLOBAL, "adjust_topo", "y");
  }

  for (n = 0; n < momgrid->nobc; n++) {
    sprintf(buf,"OBC #%d from (%d,%d) to (%d,%d)",n, momgrid->obcis[n],
	    momgrid->obcjs[n], momgrid->obcie[n], momgrid->obcje[n]);
    sprintf(attname, "obc%d", n);
    write_text_att(cdfid, NC_GLOBAL, attname, buf);
  }
  if (momgrid->cyoffs != 0) {
    if (momgrid->cytype & U1BDRY)
      sprintf(buf,"Cyclic x_boundary has offset of %d", -momgrid->cyoffs);
    else if (momgrid->cytype & U2BDRY)
      sprintf(buf,"Cyclic y_boundary has offset of %d", -momgrid->cyoffs);
    write_text_att(cdfid, NC_GLOBAL, "cyclic_OBC", buf);
  }
}


/*-------------------------------------------------------------------*/
/* Writes MOM metric data to netCDF file                             */
/*-------------------------------------------------------------------*/
void write_metrics(momgrid_t *momgrid, int mode)
{
  char grd = '\0';
  size_t start[4];
  size_t count[4];
  char var[MAXSTRLEN];
  int cfid = momgrid->cfid;
  float vertex[4] = {1, 2, 3, 4};
  int nx = 0, ny = 0;

  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;
  count[0] = 0;
  count[1] = 0;
  count[2] = 0;
  count[3] = 0;
  switch(mode) {
  case TCELL:
    grd = 'T';
    count[0] = momgrid->grid_x_T;
    nc_put_vara_float(cfid, ncw_var_id(cfid, "grid_x_T"), start, count, momgrid->grid_x_Ta);
    count[0] = momgrid->grid_y_T;
    nc_put_vara_float(cfid, ncw_var_id(cfid, "grid_y_T"), start, count, momgrid->grid_y_Ta);
    count[0] = momgrid->grid_x_C;
    nc_put_vara_float(cfid, ncw_var_id(cfid, "grid_x_C"), start, count, momgrid->grid_x_Ca);
    count[0] = momgrid->grid_y_C;
    nc_put_vara_float(cfid, ncw_var_id(cfid, "grid_y_C"), start, count, momgrid->grid_y_Ca);
    count[0] = momgrid->vertex;
    nc_put_vara_float(cfid, ncw_var_id(cfid, "vertex"), start, count, vertex);
    nx = momgrid->grid_x_T;
    ny = momgrid->grid_y_T;
    break;
  case NCELL:
    grd = 'N';
    nx = momgrid->grid_x_T;
    ny = momgrid->grid_y_C;
    break;
  case ECELL:
    grd = 'E';
    nx = momgrid->grid_x_C;
    ny = momgrid->grid_y_T;
    break;
  case CCELL:
    grd = 'C';
    nx = momgrid->grid_x_C;
    ny = momgrid->grid_y_C;
    break;
  }

  count[0] = ny;
  count[1] = nx;
  sprintf(var, "x_%c",grd);
  nc_put_vara_double(cfid, ncw_var_id(cfid, var), start, count, momgrid->x_g[0]);
  sprintf(var, "y_%c",grd);
  nc_put_vara_double(cfid, ncw_var_id(cfid, var), start, count, momgrid->y_g[0]);
  count[0] = momgrid->vertex;
  count[1] = ny;
  count[2] = nx;
  sprintf(var, "x_vert_%c",grd);
  nc_put_vara_double(cfid, ncw_var_id(cfid, var), start, count, momgrid->x_vert[0][0]);
  sprintf(var, "y_vert_%c",grd);
  nc_put_vara_double(cfid, ncw_var_id(cfid, var), start, count, momgrid->y_vert[0][0]);
  count[0] = ny;
  count[1] = nx;
  count[2] = 0;
  sprintf(var, "area_%c",grd);
  nc_put_vara_double(cfid, ncw_var_id(cfid, var), start, count, momgrid->area[0]);
  sprintf(var, "angle_%c",grd);
  nc_put_vara_double(cfid, ncw_var_id(cfid, var), start, count, momgrid->angle[0]);
  sprintf(var, "ds_00_02_%c",grd);
  nc_put_vara_double(cfid, ncw_var_id(cfid, var), start, count, momgrid->ds_00_02[0]);
  sprintf(var, "ds_20_22_%c",grd);
  nc_put_vara_double(cfid, ncw_var_id(cfid, var), start, count, momgrid->ds_20_22[0]);
  sprintf(var, "ds_02_22_%c",grd);
  nc_put_vara_double(cfid, ncw_var_id(cfid, var), start, count, momgrid->ds_02_22[0]);
  sprintf(var, "ds_00_20_%c",grd);
  nc_put_vara_double(cfid, ncw_var_id(cfid, var), start, count, momgrid->ds_00_20[0]);
  sprintf(var, "ds_00_01_%c",grd);
  nc_put_vara_double(cfid, ncw_var_id(cfid, var), start, count, momgrid->ds_00_01[0]);
  sprintf(var, "ds_01_02_%c",grd);
  nc_put_vara_double(cfid, ncw_var_id(cfid, var), start, count, momgrid->ds_01_02[0]);
  sprintf(var, "ds_02_12_%c",grd);
  nc_put_vara_double(cfid, ncw_var_id(cfid, var), start, count, momgrid->ds_02_12[0]);
  sprintf(var, "ds_12_22_%c",grd);
  nc_put_vara_double(cfid, ncw_var_id(cfid, var), start, count, momgrid->ds_12_22[0]);
  sprintf(var, "ds_21_22_%c",grd);
  nc_put_vara_double(cfid, ncw_var_id(cfid, var), start, count, momgrid->ds_21_22[0]);
  sprintf(var, "ds_20_21_%c",grd);
  nc_put_vara_double(cfid, ncw_var_id(cfid, var), start, count, momgrid->ds_20_21[0]);
  sprintf(var, "ds_10_20_%c",grd);
  nc_put_vara_double(cfid, ncw_var_id(cfid, var), start, count, momgrid->ds_10_20[0]);
  sprintf(var, "ds_00_10_%c",grd);
  nc_put_vara_double(cfid, ncw_var_id(cfid, var), start, count, momgrid->ds_00_10[0]);
  sprintf(var, "ds_01_11_%c",grd);
  nc_put_vara_double(cfid, ncw_var_id(cfid, var), start, count, momgrid->ds_01_11[0]);
  sprintf(var, "ds_11_12_%c",grd);
  nc_put_vara_double(cfid, ncw_var_id(cfid, var), start, count, momgrid->ds_11_12[0]);
  sprintf(var, "ds_11_21_%c",grd);
  nc_put_vara_double(cfid, ncw_var_id(cfid, var), start, count, momgrid->ds_11_21[0]);
  sprintf(var, "ds_10_11_%c",grd);
  nc_put_vara_double(cfid, ncw_var_id(cfid, var), start, count, momgrid->ds_10_11[0]);
  sprintf(var, "ds_01_21_%c",grd);
  nc_put_vara_double(cfid, ncw_var_id(cfid, var), start, count, momgrid->ds_01_21[0]);
  sprintf(var, "ds_10_12_%c",grd);
  nc_put_vara_double(cfid, ncw_var_id(cfid, var), start, count, momgrid->ds_10_12[0]);

}
/* END write_metrics()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Closes a MOM netCDF file                                          */
/*-------------------------------------------------------------------*/
void close_netcdf(momgrid_t *momgrid)
{
  size_t start[4];
  size_t count[4];
  int cfid = momgrid->cfid;

  start[0] = 0;
  start[1] = 0;
  count[0] = momgrid->nz;
  nc_put_vara_float(cfid, ncw_var_id(cfid, "zt"), start, count, momgrid->zt);
  nc_put_vara_float(cfid, ncw_var_id(cfid, "zb"), start, count, momgrid->zb);
  count[0] = momgrid->grid_y_T;
  count[1] = momgrid->grid_x_T;
  nc_put_vara_double(cfid, ncw_var_id(cfid, "depth_t"), start, count, momgrid->depth[0]);
  nc_put_vara_double(cfid, ncw_var_id(cfid, "num_levels"), start, count, momgrid->num_levels[0]);
  nc_put_vara_double(cfid, ncw_var_id(cfid, "wet"), start, count, momgrid->wet[0]);
  // See [MOMref1] comments
  /*
    nc_put_vara_double(cfid, ncw_var_id(cfid, "depth_c"), start, count, momgrid->depth_c[0]);
    nc_put_vara_double(cfid, ncw_var_id(cfid, "num_levels_c"), start, count, momgrid->num_levels_c[0]);
    nc_put_vara_double(cfid, ncw_var_id(cfid, "wet_c"), start, count, momgrid->wet_c[0]);
  */
  ncw_close(momgrid->ofile, momgrid->cfid);
}

/* END close_netcdf()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Writes EMS T/S data to MOM netCDF file                            */
/*-------------------------------------------------------------------*/
void write_TS(dump_data_t *dumpdata, momgrid_t *momgrid)
{
  int i, j, k, kk, n, tn, sn;
  size_t start[4] = {0, 0, 0, 0};
  size_t count[4] = {0, 0, 0, 0};
  int cfid = momgrid->cfid;
  double ***d1;

  /* Find the temp and salt tracer indicies                          */
  sn = tn = -1;
  for (n = 0; n < dumpdata->ntr; n++) {
    tracer_info_t *trinfo = dumpdata->trinfo_3d;
    if (strcmp(trinfo[n].name, "salt") == 0) sn = n;
    if (strcmp(trinfo[n].name, "temp") == 0) tn = n;
  }
  if (sn == -1 || tn == -1) {
    hd_warn("Can't find tracers 'temp' or 'salt'\n");
    return;
  }

  /* Write the T/S data to file                                      */
  count[0] = momgrid->grid_x_T;
  nc_put_vara_float(cfid, ncw_var_id(cfid, "grid_x_T"), start, count, momgrid->grid_x_Ta);
  count[0] = momgrid->grid_y_T;
  nc_put_vara_float(cfid, ncw_var_id(cfid, "grid_y_T"), start, count, momgrid->grid_y_Ta);
  count[0] = momgrid->nz;
  nc_put_vara_float(cfid, ncw_var_id(cfid, "zt"), start, count, momgrid->zt);
  count[0] = momgrid->grid_y_T;
  count[1] = momgrid->grid_x_T;
  nc_put_vara_double(cfid, ncw_var_id(cfid, "x_T"), start, count, momgrid->x_T[0]);
  nc_put_vara_double(cfid, ncw_var_id(cfid, "y_T"), start, count, momgrid->y_T[0]);

  count[0] = momgrid->nz;
  count[1] = momgrid->grid_y_T;
  count[2] = momgrid->grid_x_T;
  d1 = d_alloc_3d(momgrid->grid_x_T, momgrid->grid_y_T, momgrid->nz);
  for (k = 0; k < momgrid->nz; k++) {
    kk = momgrid->nz - 1 - k;
    for (j = 0; j < momgrid->grid_y_T; j++)
      for (i = 0; i < momgrid->grid_x_T; i++)
	d1[kk][j][i] = dumpdata->tr_wc[tn][k][j][i];
  }

  nc_put_vara_double(cfid, ncw_var_id(cfid, "temp"), start, count, d1[0][0]);
  for (k = 0; k < momgrid->nz; k++) {
    kk = momgrid->nz - 1 - k;
    for (j = 0; j < momgrid->grid_y_T; j++)
      for (i = 0; i < momgrid->grid_x_T; i++)
	d1[kk][j][i] = dumpdata->tr_wc[sn][k][j][i];
  }
  nc_put_vara_double(cfid, ncw_var_id(cfid, "salt"), start, count, d1[0][0]);

  d_free_3d(d1);
  ncw_close(momgrid->ofile, momgrid->cfid);
}

/* END write_TS()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Writes EMS eta data to MOM netCDF file                            */
/*-------------------------------------------------------------------*/
void write_eta(dump_data_t *dumpdata, momgrid_t *momgrid)
{
  int i, j;
  size_t start[4] = {0, 0, 0, 0};
  size_t count[4] = {0, 0, 0, 0};
  int cfid = momgrid->cfid;
  double **d1;

  /* Write the T/S data to file                                      */
  count[0] = momgrid->grid_x_T;
  nc_put_vara_float(cfid, ncw_var_id(cfid, "grid_x_T"), start, count, momgrid->grid_x_Ta);
  count[0] = momgrid->grid_y_T;
  nc_put_vara_float(cfid, ncw_var_id(cfid, "grid_y_T"), start, count, momgrid->grid_y_Ta);
  count[0] = momgrid->grid_y_T;
  count[1] = momgrid->grid_x_T;
  nc_put_vara_double(cfid, ncw_var_id(cfid, "x_T"), start, count, momgrid->x_T[0]);
  nc_put_vara_double(cfid, ncw_var_id(cfid, "y_T"), start, count, momgrid->y_T[0]);
  count[0] = momgrid->grid_y_T;
  count[1] = momgrid->grid_x_T;
  d1 = d_alloc_2d(momgrid->grid_x_T, momgrid->grid_y_T);
  for (j = 0; j < momgrid->grid_y_T; j++)
    for (i = 0; i < momgrid->grid_x_T; i++) {
      d1[j][i] = dumpdata->eta[j][i];
      if (isnan(d1[j][i])) d1[j][i] = 0.0;
    }
  nc_put_vara_double(cfid, ncw_var_id(cfid, "eta_t"), start, count, d1[0]);
  d_free_2d(d1);
  ncw_close(momgrid->ofile, momgrid->cfid);
}

/* END write_eta()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Prints the .cdl file header                                       */
/*-------------------------------------------------------------------*/
void print_grid_header(momgrid_t *momgrid, FILE *op, int mode)
{

  if (mode) {
    fprintf(op, "netcdf grid_spec {\n");
    fprintf(op, "dimensions:\n");
    fprintf(op, "grid_x_T = %d ;\n", momgrid->grid_x_T);
    fprintf(op, "grid_y_T = %d ;\n", momgrid->grid_y_T);
    fprintf(op, "grid_x_C = %d ;\n", momgrid->grid_x_C);
    fprintf(op, "grid_y_C = %d ;\n", momgrid->grid_y_C);
    fprintf(op, "vertex = 4 ;\n");
    fprintf(op, "zt = %d ;\n", momgrid->nz);
    fprintf(op, "zb = %d ;\n", momgrid->nz);
    fprintf(op, "variables:\n");
    fprintf(op, "        float grid_x_T(grid_x_T) ;\n");
    fprintf(op, "                grid_x_T:long_name = 'Nominal Longitude of T-cell center' ;\n");
    fprintf(op, "                grid_x_T:units = 'degree_east' ;\n");
    fprintf(op, "                grid_x_T:cartesian_axis = 'X' ;\n");
    fprintf(op, "        float grid_y_T(grid_y_T) ;\n");
    fprintf(op, "                grid_y_T:long_name = 'Nominal Latitude of T-cell center' ;\n");
    fprintf(op, "                grid_y_T:units = 'degree_north' ;\n");
    fprintf(op, "                grid_y_T:cartesian_axis = 'Y' ;\n");
    fprintf(op, "        float grid_x_C(grid_x_C) ;\n");
    fprintf(op, "                grid_x_C:long_name = 'Nominal Longitude of C-cell center' ;\n");
    fprintf(op, "                grid_x_C:units = 'degree_east' ;\n");
    fprintf(op, "                grid_x_C:cartesian_axis = 'X' ;\n");
    fprintf(op, "        float grid_y_C(grid_y_C) ;\n");
    fprintf(op, "                grid_y_C:long_name = 'Nominal Latitude of C-cell center' ;\n");
    fprintf(op, "                grid_y_C:units = 'degree_north' ;\n");
    fprintf(op, "                grid_y_C:cartesian_axis = 'Y' ;\n");
    fprintf(op, "        float vertex(vertex) ;\n");
    fprintf(op, "                vertex:long_name = 'Vertex position from southwest couterclockwise' ;\n");
    fprintf(op, "                vertex:units = 'none' ;\n");

  } else {
    fprintf(op, "        float zt(zt) ;\n");
    fprintf(op, "                zt:long_name = 'zt' ;\n");
    fprintf(op, "                zt:units = 'meters' ;\n");
    fprintf(op, "                zt:cartesian_axis = 'z' ;\n");
    fprintf(op, "                zt:positive = 'down' ;\n");
    fprintf(op, "        float zb(zb) ;\n");
    fprintf(op, "                zb:long_name = 'zb' ;\n");
    fprintf(op, "                zb:units = 'meters' ;\n");
    fprintf(op, "                zb:cartesian_axis = 'z' ;\n");
    fprintf(op, "                zb:positive = 'down' ;\n");
    fprintf(op, "        double depth_t(grid_y_T, grid_x_T) ;\n");
    fprintf(op, "                depth_t:long_name = 'topographic depth of T-cell' ;\n");
    fprintf(op, "                depth_t:units = 'meters' ;\n");
    fprintf(op, "        double num_levels(grid_y_T, grid_x_T) ;\n");
    fprintf(op, "                num_levels:long_name = 'number of vertical T-cells' ;\n");
    fprintf(op, "                num_levels:units = 'none' ;\n");
    fprintf(op, "        double wet(grid_y_T, grid_x_T) ;\n");
    fprintf(op, "                wet:long_name = 'land/sea flag (0=land) for T-cell' ;\n");
    fprintf(op, "                wet:units = 'none' ;\n");
    fprintf(op, "\n// global attributes:\n");
    fprintf(op, "                :filename = 'test1_grid.nc' ;\n");
    fprintf(op, "                :xname = 'longitude' ;\n");
    fprintf(op, "                :yname = 'latitude' ;\n");
    fprintf(op, "                :vertex_convention = 'SWCCW' ;\n");
    fprintf(op, "                :y_boundary_type = 'solid_walls' ;\n");
    fprintf(op, "                :x_boundary_type = 'solid_walls' ;\n");
    fprintf(op, "                :topography = 'rectangular_basin' ;\n");
    fprintf(op, "                :input_file = '/archive/fms/mom4/input_data/OCCAM_p5degree.nc' ;\n");
    fprintf(op, "                :input_field = 'topo' ;\n");
    fprintf(op, "                :fill_isolated_cells = 'y' ;\n");
    fprintf(op, "                :adjust_topo = 'y' ;\n");
    fprintf(op, "data:\n\n");

    fprintf(op, "grid_x_T = \n");
    print_farray_1d(op, momgrid->grid_x_T, momgrid->grid_x_Ta);
    fprintf(op, "\ngrid_y_T = \n");
    print_farray_1d(op, momgrid->grid_y_T, momgrid->grid_y_Ta);
    fprintf(op, "\ngrid_x_C = \n");
    print_farray_1d(op, momgrid->grid_x_C, momgrid->grid_x_Ca);
    fprintf(op, "\ngrid_y_C = \n");
    print_farray_1d(op, momgrid->grid_y_C, momgrid->grid_y_Ca);
    fprintf(op, "\nvertex = 1, 2, 3, 4 ;\n");
  }
}

/* END print_grid_header()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Prints the .cdl file attributes                                   */
/*-------------------------------------------------------------------*/
void print_grid_atts(momgrid_t *momgrid, FILE *op, int mode)
{
  char grd = '\0';
  char gx = 0;
  char gy = 0;

  switch(mode) {
  case TCELL:
    grd = 'T';
    gx = 'T';
    gy = 'T';
    break;
  case NCELL:
    grd = 'N';
    gx = 'T';
    gy = 'C';
    break;
  case ECELL:
    grd = 'E';
    gx = 'C';
    gy = 'T';
    break;
  case CCELL:
    grd = 'C';
    gx = 'C';
    gy = 'C';
    break;
  }

  fprintf(op, "        double x_%c(grid_y_%c, grid_x_%c) ;\n", grd, gy, gx);
  fprintf(op, "                x_%c:long_name = 'Geographic longitude of %c_cell centers' ;\n", grd, grd);
  fprintf(op, "                x_%c:units = 'degree_east' ;\n", grd);
  fprintf(op, "        double y_%c(grid_y_%c, grid_x_%c) ;\n", grd, gy, gx);
  fprintf(op, "                y_%c:long_name = 'Geographic latitude of %c_cell centers' ;\n", grd, grd);
  fprintf(op, "                y_%c:units = 'degree_north' ;\n", grd);
  fprintf(op, "        double x_vert_%c(vertex, grid_y_%c, grid_x_%c) ;\n", grd, gy, gx);
  fprintf(op, "                x_vert_%c:long_name = 'Geographic longitude of %c_cell vertices begin southwest counterclockwise' ;\n", grd, grd);
  fprintf(op, "                x_vert_%c:units = 'degree_east' ;\n", grd);
  fprintf(op, "        double y_vert_%c(vertex, grid_y_%c, grid_x_%c) ;\n", grd, gy, gx);
  fprintf(op, "                y_vert_%c:long_name = 'Geographic latitude of %c_cell vertices begin southwest counterclockwise' ;\n", grd, grd);
  fprintf(op, "                y_vert_%c:units = 'degree_north' ;\n", grd);
  fprintf(op, "        double area_%c(grid_y_%c, grid_x_%c) ;\n", grd, gy, gx);
 fprintf(op, "                area_%c:long_name = 'Area of %c_cell' ;\n", grd, grd);
 fprintf(op, "                area_%c:units = 'm2' ;\n", grd);
 fprintf(op, "        double angle_%c(grid_y_%c, grid_x_%c) ;\n", grd, gy, gx);
 fprintf(op, "                angle_%c:long_name = 'Angle clockwise between logical and geographic east of %c_cell' ;\n", grd, grd);
 fprintf(op, "                angle_%c:units = 'degree' ;\n", grd);
 fprintf(op, "        double ds_00_02_%c(grid_y_%c, grid_x_%c) ;\n", grd, gy, gx);
 fprintf(op, "                ds_00_02_%c:long_name = 'Length of western face of %c_cell' ;\n", grd, grd);
 fprintf(op, "                ds_00_02_%c:units = 'm' ;\n", grd);
 fprintf(op, "        double ds_20_22_%c(grid_y_%c, grid_x_%c) ;\n", grd, gy, gx);
 fprintf(op, "                ds_20_22_%c:long_name = 'Length of eastern face of %c_cell' ;\n", grd, grd);
 fprintf(op, "                ds_20_22_%c:units = 'm' ;\n", grd);
 fprintf(op, "        double ds_02_22_%c(grid_y_%c, grid_x_%c) ;\n", grd, gy, gx);
 fprintf(op, "                ds_02_22_%c:long_name = 'Length of northern face of %c_cell' ;\n", grd, grd);
 fprintf(op, "                ds_02_22_%c:units = 'm' ;\n", grd);
 fprintf(op, "        double ds_00_20_%c(grid_y_%c, grid_x_%c) ;\n", grd, gy, gx);
 fprintf(op, "                ds_00_20_%c:long_name = 'Length of southern face of %c_cell' ;\n", grd, grd);
 fprintf(op, "                ds_00_20_%c:units = 'm' ;\n", grd);
 fprintf(op, "        double ds_00_01_%c(grid_y_%c, grid_x_%c) ;\n", grd, gy, gx);
 fprintf(op, "                ds_00_01_%c:long_name = 'Distance from southwest corner to western face center of %c_cell' ;\n", grd, grd);
 fprintf(op, "                ds_00_01_%c:units = 'm' ;\n", grd);
 fprintf(op, "        double ds_01_02_%c(grid_y_%c, grid_x_%c) ;\n", grd, gy, gx);
 fprintf(op, "                ds_01_02_%c:long_name = 'Distance from northwest corner to western face center of %c_cell' ;\n", grd, grd);
 fprintf(op, "                ds_01_02_%c:units = 'm' ;\n", grd);
 fprintf(op, "        double ds_02_12_%c(grid_y_%c, grid_x_%c) ;\n", grd, gy, gx);
 fprintf(op, "                ds_02_12_%c:long_name = 'Distance from northwest corner to northern face center of %c_cell' ;\n", grd, grd);
 fprintf(op, "                ds_02_12_%c:units = 'm' ;\n", grd);
 fprintf(op, "        double ds_12_22_%c(grid_y_%c, grid_x_%c) ;\n", grd, gy, gx);
 fprintf(op, "                ds_12_22_%c:long_name = 'Distance from northeast corner to northern face center of %c_cell' ;\n", grd, grd);
 fprintf(op, "                ds_12_22_%c:units = 'm' ;\n", grd);
 fprintf(op, "        double ds_21_22_%c(grid_y_%c, grid_x_%c) ;\n", grd, gy, gx);
 fprintf(op, "                ds_21_22_%c:long_name = 'Distance from northeast corner to eastern face center of %c_cell' ;\n", grd, grd);
 fprintf(op, "                ds_21_22_%c:units = 'm' ;\n", grd);
 fprintf(op, "        double ds_20_21_%c(grid_y_%c, grid_x_%c) ;\n", grd, gy, gx);
 fprintf(op, "                ds_20_21_%c:long_name = 'Distance from southeast corner to eastern face center of %c_cell' ;\n", grd, grd);
 fprintf(op, "                ds_20_21_%c:units = 'm' ;\n", grd);
 fprintf(op, "        double ds_10_20_%c(grid_y_%c, grid_x_%c) ;\n", grd, gy, gx);
 fprintf(op, "                ds_10_20_%c:long_name = 'Distance from southeast corner to southern face center of %c_cell' ;\n", grd, grd);
 fprintf(op, "                ds_10_20_%c:units = 'm' ;\n", grd);
 fprintf(op, "        double ds_00_10_%c(grid_y_%c, grid_x_%c) ;\n", grd, gy, gx);
 fprintf(op, "                ds_00_10_%c:long_name = 'Distance from southwest corner to southern face center of %c_cell' ;\n", grd, grd);
 fprintf(op, "                ds_00_10_%c:units = 'm' ;\n", grd);
 fprintf(op, "        double ds_01_11_%c(grid_y_%c, grid_x_%c) ;\n", grd, gy, gx);
 fprintf(op, "                ds_01_11_%c:long_name = 'Distance from center to western face of %c_cell' ;\n", grd, grd);
 fprintf(op, "                ds_01_11_%c:units = 'm' ;\n", grd);
 fprintf(op, "        double ds_11_12_%c(grid_y_%c, grid_x_%c) ;\n", grd, gy, gx);
 fprintf(op, "                ds_11_12_%c:long_name = 'Distance from center to northern face of %c_cell' ;\n", grd, grd);
 fprintf(op, "                ds_11_12_%c:units = 'm' ;\n", grd);
 fprintf(op, "        double ds_11_21_%c(grid_y_%c, grid_x_%c) ;\n", grd, gy, gx);
 fprintf(op, "                ds_11_21_%c:long_name = 'Distance from center to eastern face of %c_cell' ;\n", grd, grd);
 fprintf(op, "                ds_11_21_%c:units = 'm' ;\n", grd);
 fprintf(op, "        double ds_10_11_%c(grid_y_%c, grid_x_%c) ;\n", grd, gy, gx);
 fprintf(op, "                ds_10_11_%c:long_name = 'Distance from center to southern face of %c_cell' ;\n", grd, grd);
 fprintf(op, "                ds_10_11_%c:units = 'm' ;\n", grd);
 fprintf(op, "        double ds_01_21_%c(grid_y_%c, grid_x_%c) ;\n", grd, gy, gx);
 fprintf(op, "                ds_01_21_%c:long_name = 'width of %c_cell' ;\n", grd, grd);
 fprintf(op, "                ds_01_21_%c:units = 'm' ;\n", grd);
 fprintf(op, "        double ds_10_12_%c(grid_y_%c, grid_x_%c) ;\n", grd, gy, gx);
 fprintf(op, "                ds_10_12_%c:long_name = 'height of %c_cell' ;\n", grd, grd);
 fprintf(op, "                ds_10_12_%c:units = 'm' ;\n", grd);

}

/* END print_grid_atts()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Prints the grid metrics to file                                   */
/*-------------------------------------------------------------------*/
void print_metrics_cdl(momgrid_t *momgrid, FILE *op, int mode)
{
  char grd = '\0';
  int nx = 0, ny = 0;

  switch(mode) {
  case TCELL:
    grd = 'T';
    nx = momgrid->grid_x_T;
    ny = momgrid->grid_y_T;
    break;
  case NCELL:
    grd = 'N';
    nx = momgrid->grid_x_T;
    ny = momgrid->grid_y_C;
    break;
  case ECELL:
    grd = 'E';
    nx = momgrid->grid_x_C;
    ny = momgrid->grid_y_T;
    break;
  case CCELL:
    grd = 'C';
    nx = momgrid->grid_x_C;
    ny = momgrid->grid_y_C;
    break;
  }

  fprintf(op, "\n\nx_%c =\n",grd);
  print_array_2d(op, nx, ny, momgrid->x_g);

  fprintf(op, "\n\ny_%c =\n",grd);
  print_array_2d(op, nx, ny, momgrid->y_g);

  fprintf(op, "\n\nx_vert_%c =\n",grd);
  print_array_3d(op, nx, ny, momgrid->vertex, momgrid->x_vert);

  fprintf(op, "\n\ny_vert_%c =\n",grd);
  print_array_3d(op, nx, ny, momgrid->vertex, momgrid->y_vert);

  fprintf(op, "\n\narea_%c =\n",grd);
  print_array_2d(op, nx, ny, momgrid->area);

  fprintf(op, "\n\nangle_%c =\n",grd);
  print_array_2d(op, nx, ny, momgrid->angle);

  fprintf(op, "\n\nds_00_02_%c =\n",grd);
  print_array_2d(op, nx, ny, momgrid->ds_00_02);
  fprintf(op, "\n\nds_20_22_%c =\n",grd);
  print_array_2d(op, nx, ny, momgrid->ds_20_22);
  fprintf(op, "\n\nds_02_22_%c =\n",grd);
  print_array_2d(op, nx, ny, momgrid->ds_02_22);
  fprintf(op, "\n\nds_00_20_%c =\n",grd);
  print_array_2d(op, nx, ny, momgrid->ds_00_20);
  fprintf(op, "\n\nds_00_01_%c =\n",grd);
  print_array_2d(op, nx, ny, momgrid->ds_00_01);
  fprintf(op, "\n\nds_01_02_%c =\n",grd);
  print_array_2d(op, nx, ny, momgrid->ds_01_02);
  fprintf(op, "\n\nds_02_12_%c =\n",grd);
  print_array_2d(op, nx, ny, momgrid->ds_02_12);
  fprintf(op, "\n\nds_12_22_%c =\n",grd);
  print_array_2d(op, nx, ny, momgrid->ds_12_22);
  fprintf(op, "\n\nds_21_22_%c =\n",grd);
  print_array_2d(op, nx, ny, momgrid->ds_21_22);
  fprintf(op, "\n\nds_20_21_%c =\n",grd);
  print_array_2d(op, nx, ny, momgrid->ds_20_21);
  fprintf(op, "\n\nds_10_20_%c =\n",grd);
  print_array_2d(op, nx, ny, momgrid->ds_10_20);
  fprintf(op, "\n\nds_00_10_%c =\n",grd);
  print_array_2d(op, nx, ny, momgrid->ds_00_10);
  fprintf(op, "\n\nds_01_11_%c =\n",grd);
  print_array_2d(op, nx, ny, momgrid->ds_01_11);
  fprintf(op, "\n\nds_11_12_%c =\n",grd);
  print_array_2d(op, nx, ny, momgrid->ds_11_12);
  fprintf(op, "\n\nds_11_21_%c =\n",grd);
  print_array_2d(op, nx, ny, momgrid->ds_11_21);
  fprintf(op, "\n\nds_10_11_%c =\n",grd);
  print_array_2d(op, nx, ny, momgrid->ds_10_11);
  fprintf(op, "\n\nds_01_21_%c =\n",grd);
  print_array_2d(op, nx, ny, momgrid->ds_01_21);
  fprintf(op, "\n\nds_10_12_%c =\n",grd);
  print_array_2d(op, nx, ny, momgrid->ds_10_12);
  fprintf(op, "\n\nds_00_02_%c =\n",grd);
  print_array_2d(op, nx, ny, momgrid->ds_00_02);

}

/* END print_metrics()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Prints an array to file                                           */
/*-------------------------------------------------------------------*/
void print_array_2d(FILE *op, int nce1, int nce2, double **a)
{
  int i, j, n;
  int nline = 8;

  for (j = 0; j < nce2; j++) {
    n = 0;
    if (j > 0) fprintf(op,"\n");
    for (i = 0; i < nce1; i++) {
      fprintf(op,"%f",a[j][i]);
      if (i == nce1-1 && j == nce2-1)
	fprintf(op,";");
      else
	fprintf(op,", ");
      if (n == nline) {
	fprintf(op,"\n  ");
	n = 0;
      }
      n++;
    }
  }
}

/* END print_array_2d()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Prints an array to file                                           */
/*-------------------------------------------------------------------*/
void print_array_3d(FILE *op, int nce1, int nce2, int nk, double ***a)
{
  int i, j, k, n;
  int nline = 8;

  for (k = 0; k < nk; k++) {
    for (j = 0; j < nce2; j++) {
      n = 0;
      for (i = 0; i < nce1; i++) {
	fprintf(op,"%f ",a[k][j][i]);
	if (i == nce1-1 && j == nce2-1 && k == nk-1)
	  fprintf(op,";");
	else
	  fprintf(op,", ");
	if (n == nline) {
	  fprintf(op,"\n  ");
	  n = 0;
	}
	n++;
      }
      fprintf(op,"\n");
    }
    fprintf(op,"\n");
  }
}

/* END print_array_3d()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Prints an array to file                                           */
/*-------------------------------------------------------------------*/
void print_array_1d(FILE *op, int nce1, double *a)
{
  int i, n = 0;
  int nline = 8;

  for (i = 0; i < nce1; i++) {
    fprintf(op,"%f",a[i]);
    if (i == nce1-1)
      fprintf(op,";");
    else
      fprintf(op,", ");
    if (n == nline) {
      fprintf(op,"\n  ");
      n = 0;
    }
    n++;
  }
}

/* END print_array_1d()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Prints an array to file                                           */
/*-------------------------------------------------------------------*/
void print_farray_1d(FILE *op, int nce1, float *a)
{
  int i, n = 0;
  int nline = 8;

  for (i = 0; i < nce1; i++) {
    fprintf(op,"%f",a[i]);
    if (i == nce1-1)
      fprintf(op,";");
    else
      fprintf(op,", ");
    if (n == nline) {
      fprintf(op,"\n  ");
      n = 0;
    }
    n++;
  }
}

/* END print_farray_1d()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Closes a .cdl file                                                */
/*-------------------------------------------------------------------*/
void close_cdl(momgrid_t *momgrid, FILE *op)
{

  fprintf(op, "zt = \n");
  print_farray_1d(op, momgrid->nz, momgrid->zt);
  fprintf(op, "\nzb = \n");
  print_farray_1d(op, momgrid->nz, momgrid->zb);
  fprintf(op, "\ndepth_t = \n");
  print_array_2d(op, momgrid->grid_x_T, momgrid->grid_y_T, momgrid->depth);
  fprintf(op, "\nnum_levels = \n");
  print_array_2d(op, momgrid->grid_x_T, momgrid->grid_y_T, momgrid->num_levels);
  fprintf(op, "\nwet = \n");
  print_array_2d(op, momgrid->grid_x_T, momgrid->grid_y_T, momgrid->wet);
  fprintf(op, "\n}");
  fclose(op);
}

/* END close_cdl()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads a MOM grid_spec.nc file and converts to text .prm file      */
/*-------------------------------------------------------------------*/
momgrid_t *read_grid_spec(char *fname, char *ofile)
{
  FILE *op;
  float *zb;
  int fid, n;
  int i, j, ni, nj, mi, mj;
  size_t start[4];
  size_t count[4];
  size_t grid_x_T, grid_y_T;
  size_t grid_x_C, grid_y_C;
  size_t vertex;
  size_t nz;
  momgrid_t *momgrid;
  int ncerr;

  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;
  /* Time independent variables */
  count[0] = 0;
  count[1] = 0;
  count[2] = 0;
  count[3] = 0;

  /*-----------------------------------------------------------------*/
  /* Set up the grid structure                                       */
  momgrid = momgrid_alloc();
  if (ofile != NULL)
    strcpy(momgrid->ofile, ofile);

  /* Open the dump file for reading                                  */
  if ((ncerr = nc_open(fname, NC_NOWRITE, &fid)) != NC_NOERR) {
    printf("Can't find input file %s\n", fname);
    hd_quit((char *)nc_strerror(ncerr));
  }

  /*-----------------------------------------------------------------*/
  /* Get MOM grid dimensions                                         */
  nc_inq_dimlen(fid, ncw_dim_id(fid, "grid_x_T"), &grid_x_T);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "grid_y_T"), &grid_y_T);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "grid_x_C"), &grid_x_C);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "grid_y_C"), &grid_y_C);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "vertex"), &vertex);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "zb"), &nz);
  momgrid->grid_x_T = (int)grid_x_T;
  momgrid->grid_y_T = (int)grid_y_T;
  momgrid->grid_x_C = (int)grid_x_C;
  momgrid->grid_y_C = (int)grid_y_C;
  momgrid->vertex = (int)vertex;
  momgrid->nz = (int)nz;
  momgrid->nce1 = momgrid->grid_x_T + 1;
  momgrid->nce2 = momgrid->grid_y_T + 1;
  momgrid->nfe1 = momgrid->nce1 + 1;
  momgrid->nfe2 = momgrid->nce2 + 1;
  momgrid->cfid = fid;
  momgrid->layers = d_alloc_1d(momgrid->nz + 1);
  momgrid->topo = d_alloc_2d(momgrid->nce1, momgrid->nce2);
  alloc_grid_mem(momgrid, momgrid->nce1, momgrid->nce2);
  alloc_mom_mem(momgrid, momgrid->grid_x_T, momgrid->grid_y_T);
  momgrid->bathy = d_alloc_1d(momgrid->nce1 * momgrid->nce2);

  /*-----------------------------------------------------------------*/
  /* Convert the layers                                              */
  zb = f_alloc_1d(momgrid->nz);
  count[0] = momgrid->nz;
  nc_get_vara_float(fid, ncw_var_id(fid, "zb"), start, count, zb);
  for (n = 0; n < momgrid->nz; n++)
    momgrid->layers[momgrid->nz - n - 1] = -zb[n];
  momgrid->layers[momgrid->nz] = 0.0;
  f_free_1d(zb);

  /*-----------------------------------------------------------------*/
  /* Convert the bathymetry                                          */
  count[0] = momgrid->grid_y_T;
  count[1] = momgrid->grid_x_T;
  nc_get_vara_double(fid, ncw_var_id(fid, "depth_t"), start, count,
                     momgrid->depth[0]);
  nc_get_vara_double(fid, ncw_var_id(fid, "wet"), start, count,
                     momgrid->wet[0]);
  for (j = 0; j < momgrid->grid_y_T; j++)
    for (i = 0; i < momgrid->grid_x_T; i++) {
      momgrid->topo[j][i] = momgrid->depth[j][i];
      if (momgrid->wet[j][i] == 0.0)
	momgrid->topo[j][i] = -LANDCELL;
    }
  j = momgrid->nce2 - 1;
  for (i = 0; i < momgrid->grid_x_T; i++)
    momgrid->topo[j][i] = momgrid->depth[j-1][i];
  i = momgrid->nce1 - 1;
  for (j = 0; j < momgrid->grid_y_T; j++)
    momgrid->topo[j][i] = momgrid->depth[j][i-1];
  i = momgrid->nce1 - 1;
  j = momgrid->nce2 - 1;
  n = 0;
  momgrid->topo[j][i] = momgrid->depth[j-1][i-1];
  for (j = 0; j < momgrid->nce2; j++)
    for (i = 0; i < momgrid->nce1; i++) {
      momgrid->bathy[n] = momgrid->topo[j][i];
      n++;
    }

  /*-----------------------------------------------------------------*/
  /* Read the T cell centered geographic locations                   */
  count[0] = momgrid->grid_y_T;
  count[1] = momgrid->grid_x_T;
  nc_get_vara_double(fid, ncw_var_id(fid, "x_T"), start, count,
                     momgrid->x_g[0]);
  nc_get_vara_double(fid, ncw_var_id(fid, "y_T"), start, count,
                     momgrid->y_g[0]);
  count[0] = momgrid->vertex;
  count[1] = momgrid->grid_y_T;
  count[2] = momgrid->grid_x_T;
  nc_get_vara_double(fid, ncw_var_id(fid, "x_vert_T"), start, count,
                     momgrid->x_vert[0][0]);
  nc_get_vara_double(fid, ncw_var_id(fid, "y_vert_T"), start, count,
                     momgrid->y_vert[0][0]);
  for (j = 0; j < momgrid->grid_y_T; j++)
    for (i = 0; i < momgrid->grid_x_T; i++) {
      momgrid->x[2*j+1][2*i+1] = momgrid->x_g[j][i];
      momgrid->y[2*j+1][2*i+1] = momgrid->y_g[j][i];
      momgrid->x[2*j][2*i] = momgrid->x_vert[0][j][i];
      momgrid->y[2*j][2*i] = momgrid->y_vert[0][j][i];
      momgrid->x[2*j][2*i+2] = momgrid->x_vert[1][j][i];
      momgrid->y[2*j][2*i+2] = momgrid->y_vert[1][j][i];
      momgrid->x[2*j+2][2*i+2] = momgrid->x_vert[2][j][i];
      momgrid->y[2*j+2][2*i+2] = momgrid->y_vert[2][j][i];
      momgrid->x[2*j+2][2*i] = momgrid->x_vert[3][j][i];
      momgrid->y[2*j+2][2*i] = momgrid->y_vert[3][j][i];
    }
  /* Get the locations of the last cell : nce1 edge                    */
  ni = 2 * momgrid->nce1 - 2;
  mi = momgrid->grid_x_C - 1;
  for (j = 0; j < momgrid->grid_y_T; j++) {
    momgrid->x[2*j][ni+2] = 2 * momgrid->x[2*j][ni] - momgrid->x[2*j][ni-2];
    momgrid->y[2*j][ni+2] = 2 * momgrid->y[2*j][ni] - momgrid->y[2*j][ni-2];
    momgrid->x[2*j+2][ni+2] = 2 * momgrid->x[2*j+2][ni] - momgrid->x[2*j+2][ni-2];
    momgrid->y[2*j+2][ni+2] = 2 * momgrid->y[2*j+2][ni] - momgrid->y[2*j+2][ni-2];
  }
  /* Get the locations of the last cell : nce2 edge                    */
  nj = 2 * momgrid->nce2 - 2;
  mj = momgrid->grid_y_C - 1;
  for (i = 0; i < momgrid->grid_x_T; i++) {
    momgrid->x[nj+2][2*i] = 2 * momgrid->x[nj][2*i] - momgrid->x[nj-2][2*i];
    momgrid->y[nj+2][2*i] = 2 * momgrid->y[nj][2*i] - momgrid->y[nj-2][2*1];
    momgrid->x[nj+2][2*i+2] = 2 * momgrid->x[nj][2*i+2] - momgrid->x[nj-2][2*i+2];
    momgrid->y[nj+2][2*i+2] = 2 * momgrid->y[nj][2*i+2] - momgrid->y[nj-2][2*1+2];
  }

  /*-------------------------------------------------------------------*/
  /* Read the E cell centered geographic locations                     */
  count[0] = momgrid->vertex;
  count[1] = momgrid->grid_y_T;
  count[2] = momgrid->grid_x_C;
  nc_get_vara_double(fid, ncw_var_id(fid, "x_vert_E"), start, count,
                     momgrid->x_vert[0][0]);
  nc_get_vara_double(fid, ncw_var_id(fid, "y_vert_E"), start, count,
                     momgrid->y_vert[0][0]);
  for (j = 0; j < momgrid->grid_y_T; j++)
    for (i = 0; i < momgrid->grid_x_C; i++) {
      momgrid->x[2*j][2*i+1] = momgrid->x_vert[0][j][i];
      momgrid->y[2*j][2*i+1] = momgrid->y_vert[0][j][i];
      momgrid->x[2*j+2][2*i+1] = momgrid->x_vert[3][j][i];
      momgrid->y[2*j+2][2*i+1] = momgrid->y_vert[3][j][i];
    }
  /* Get the locations of the last cell : nce2 edge                    */
  nj = 2 * momgrid->nce2 - 2;
  mj = momgrid->grid_y_T - 1;
  for (i = 0; i < momgrid->grid_x_C; i++) {
    momgrid->x[nj+2][2*i+1] = 2 * momgrid->x[nj][2*i+1] - momgrid->x[nj-2][2*i+1];
    momgrid->y[nj+2][2*i+1] = 2 * momgrid->y[nj][2*i+1] - momgrid->y[nj-2][2*1+1];
  }

  /*-------------------------------------------------------------------*/
  /* Read the N cell centered geographic locations                     */
  count[0] = momgrid->vertex;
  count[1] = momgrid->grid_y_C;
  count[2] = momgrid->grid_x_T;
  nc_get_vara_double(fid, ncw_var_id(fid, "x_vert_N"), start, count,
                     momgrid->x_vert[0][0]);
  nc_get_vara_double(fid, ncw_var_id(fid, "y_vert_N"), start, count,
                     momgrid->y_vert[0][0]);
  for (j = 0; j < momgrid->grid_y_C; j++)
    for (i = 0; i < momgrid->grid_x_T; i++) {
      momgrid->x[2*j+1][2*i] = momgrid->x_vert[0][j][i];
      momgrid->y[2*j+1][2*i] = momgrid->y_vert[0][j][i];
      momgrid->x[2*j+1][2*i+2] = momgrid->x_vert[1][j][i];
      momgrid->y[2*j+1][2*i+2] = momgrid->y_vert[1][j][i];
    }
  /* Get the locations of the last cell : nce1 edge                    */
  ni = 2 * momgrid->nce1 - 2;
  mi = momgrid->grid_x_T - 1;
  for (j = 0; j < momgrid->grid_y_C; j++) {
    momgrid->x[2*j+1][ni+2] = 2 * momgrid->x[2*j+1][ni] - momgrid->x[2*j+1][ni-2];
    momgrid->y[2*j+1][ni+2] = 2 * momgrid->y[2*j+1][ni] - momgrid->y[2*j+1][ni-2];
  } 

  /*-------------------------------------------------------------------*/
  /* Get the locations of the last cell : (nce1-1, nce1-1)           */
  ni = 2 * momgrid->nce1 - 2;
  nj = 2 * momgrid->nce2 - 2;
  momgrid->x[nj+2][ni+2] = 2 * momgrid->x[nj+2][ni] - momgrid->x[nj+2][ni-2];
  momgrid->y[nj+2][ni+2] = 2 * momgrid->y[nj][ni+2] - momgrid->y[nj-2][ni+2];
  ni = 2 * momgrid->nce1 - 1;
  nj = 2 * momgrid->nce2;
  momgrid->x[nj][ni] = 0.5 * (momgrid->x[nj][ni+1] + momgrid->x[nj][ni-1]);
  momgrid->y[nj][ni] = 0.5 * (momgrid->y[nj][ni+1] + momgrid->y[nj][ni-1]);
  ni = 2 * momgrid->nce1;
  nj = 2 * momgrid->nce2 - 1;
  momgrid->x[nj][ni] = 0.5 * (momgrid->x[nj+1][ni] + momgrid->x[nj-1][ni]);
  momgrid->y[nj][ni] = 0.5 * (momgrid->y[nj+1][ni] + momgrid->y[nj-1][ni]);

  /*-------------------------------------------------------------------*/
  /* Get the nce1 and nce2 positions by interpolation                  */
  nj = 2 * momgrid->nce2 - 1;
  ni = 2 * momgrid->nce1 - 1;
  for (i = 0; i < ni; i++) {
    momgrid->x[nj][i] = 0.5 * (momgrid->x[nj+1][i] + momgrid->x[nj-1][i]);
    momgrid->y[nj][i] = 0.5 * (momgrid->y[nj+1][i] + momgrid->y[nj-1][i]);
  }
  nj = 2 * momgrid->nce2;
  ni = 2 * momgrid->nce1 - 1;
  for (j = 0; j <= nj; j++) {
    momgrid->x[j][ni] = 0.5 * (momgrid->x[j][ni+1] + momgrid->x[j][ni-1]);
    momgrid->y[j][ni] = 0.5 * (momgrid->y[j][ni+1] + momgrid->y[j][ni-1]);
  }

  /*-----------------------------------------------------------------*/
  /* Open output file and write to file                              */
  if (ofile != NULL) {
    printf("\n** MOM grid generator **\n\n");
    printf("Input parameter file = %s\n", fname);
    printf("Output .prm file = %s\n", momgrid->ofile);
    printf("Grid size = %d x %d x %d\n", momgrid->nce1, momgrid->nce2, momgrid->nz);
    printf("Create ems input file with 'shoc -ag %s'\n", momgrid->ofile);
    printf("\n");
    if((op = fopen(ofile, "w")) == NULL)
      hd_quit("Can't open file %s\n", ofile);
    print_ems_grid(momgrid, op, fname);
    fclose(op);
  }
  ni = (2 * momgrid->nce1 + 1);
  nj = (2 * momgrid->nce2 + 1);
  mom_free(momgrid);
  mom_grid_free(momgrid, UNUSED|ALL);
  ncw_close(fname, fid);
  return(momgrid);
}

/* END read_grid_spec()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Writes EMS grid specification                                     */
/*-------------------------------------------------------------------*/
void print_ems_grid(momgrid_t *momgrid, FILE *op, char *fname)
{
  int i, j, n;
  int ni, nj;
  char buf[MAXSTRLEN];

  fprintf(op, "CODEHEADER           SHOC default version\n");
  fprintf(op, "PARAMETERHEADER      Converted MOM grid_spec file\n");
  fprintf(op, "DESCRIPTION          Converted from MOM %s\n", fname);
  fprintf(op, "NAME                 MOM conversion grid\n");
  fprintf(op, "TIMEUNIT             seconds since 1990-01-01 00:00:00 +10\n");
  fprintf(op, "OUTPUT_TIMEUNIT      seconds since 1990-01-01 00:00:00 +10\n");
  fprintf(op, "LENUNIT              metre\n");
  fprintf(op, "START_TIME           0.0\n");
  fprintf(op, "STOP_TIME            1000\n");
  strncpy(buf, fname, strlen(fname) - 3);
  fprintf(op, "\nINPUT_FILE  %s_ems\n\n", buf);
  fprintf(op, "PROJECTION  GEOGRAPHIC\n");
  fprintf(op, "GRIDTYPE    NUMERICAL\n");
  fprintf(op, "NCE1        %d\n", momgrid->nce1);
  fprintf(op, "NCE2        %d\n", momgrid->nce2);
  fprintf(op, "\n");
  fprintf(op, "LAYERFACES  %d\n", momgrid->nz+1);
  for (i = 0; i <= momgrid->nz; i++)
    fprintf(op, "%f\n", momgrid->layers[i]);
  fprintf(op, "\n");
  n = momgrid->nce1 * momgrid->nce2;
  fprintf(op, "BATHY     %d\n", n);
  for (i = 0; i < n; i++)
    fprintf(op, "%f\n", momgrid->bathy[i]);
  fprintf(op, "\n");
  ni = 2 * momgrid->nce1 + 1;
  nj = 2 * momgrid->nce2 + 1;
  n = ni * nj;
  fprintf(op, "XCOORDS   %d\n", n);
  for (j = 0; j < nj; j++)
    for (i = 0; i < ni; i++)
    fprintf(op, "%f\n", momgrid->x[j][i]);
  fprintf(op, "\n");
  fprintf(op, "YCOORDS   %d\n", n);
  for (j = 0; j < nj; j++)
    for (i = 0; i < ni; i++)
    fprintf(op, "%f\n", momgrid->y[j][i]);
}

/* END print_ems_grid()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Converts MOM grid_spec.nc file to EMS layer structure             */
/*-------------------------------------------------------------------*/
double *read_mom_layers(char *fname, int *nz)
{
  float *zb;
  int fid, n;
  size_t start[4];
  size_t count[4];
  size_t *zn = NULL;
  double *layers = NULL;
  int ncerr;

  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;
  /* Time independent variables */
  count[0] = 0;
  count[1] = 0;
  count[2] = 0;
  count[3] = 0;

  /*-----------------------------------------------------------------*/
  /* Open the dump file for reading                                  */
  if ((ncerr = nc_open(fname, NC_NOWRITE, &fid)) != NC_NOERR) {
    printf("Can't find input file %s\n", fname);
    hd_quit((char *)nc_strerror(ncerr));
  }

  /*-----------------------------------------------------------------*/
  /* Get MOM grid dimensions                                         */
  nc_inq_dimlen(fid, ncw_dim_id(fid, "zb"), zn);
  nz = (int *)zn;

  /*-----------------------------------------------------------------*/
  /* Convert the layers                                              */
  layers = d_alloc_1d(*nz + 1);  
  zb = f_alloc_1d(*nz);
  count[0] = *nz;
  nc_get_vara_float(fid, ncw_var_id(fid, "zb"), start, count, zb);
  for (n = 0; n < *nz; n++)
    layers[*nz - n - 1] = -zb[n];
  layers[*nz] = 0.0;
  f_free_1d(zb);

  ncw_close(fname, fid);
  return(layers);
}

/* END read_mom_layers()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Converts MOM grid_spec.nc file to EMS metrics                     */
/*-------------------------------------------------------------------*/
void read_mom_grid(char *fname, int nx, int ny, double **x, double **y)
{
  int fid;
  int i, j, ni, nj, mi, mj;
  size_t start[4];
  size_t count[4];
  size_t grid_x_T, grid_y_T;
  size_t grid_x_C, grid_y_C;
  int nce1, nce2;
  size_t vertex;
  int gs;
  double **d1, **d2;
  double ***v1, ***v2;
  int ncerr;

  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;
  /* Time independent variables */
  count[0] = 0;
  count[1] = 0;
  count[2] = 0;
  count[3] = 0;

  /*-----------------------------------------------------------------*/
  /* Open the dump file for reading                                  */
  if ((ncerr = nc_open(fname, NC_NOWRITE, &fid)) != NC_NOERR) {
    printf("Can't find input file %s\n", fname);
    hd_quit((char *)nc_strerror(ncerr));
  }

  /*-----------------------------------------------------------------*/
  /* Get MOM grid dimensions                                         */
  gs = (gsize & TCELL) ? 0 : 1;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "grid_x_T"), &grid_x_T);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "grid_y_T"), &grid_y_T);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "grid_x_C"), &grid_x_C);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "grid_y_C"), &grid_y_C);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "vertex"), &vertex);
  nce1 = grid_x_T + gs;
  nce2 = grid_y_T + gs;
  if (nce1 != nx || nce2 != ny)
    hd_quit("Inconsistent grid sizes EMS (%dx%d) : MOM (%dx%d)\n",
	    nx, ny, nce1, nce2);

  /*-----------------------------------------------------------------*/
  /* Read the T cell centered geographic locations                   */
  count[0] = grid_y_T;
  count[1] = grid_x_T;
  d1 = d_alloc_2d(count[1], count[0]);
  d2 = d_alloc_2d(count[1], count[0]);
  nc_get_vara_double(fid, ncw_var_id(fid, "x_T"), start, count, d1[0]);
  nc_get_vara_double(fid, ncw_var_id(fid, "y_T"), start, count, d2[0]);
  count[0] = vertex;
  count[1] = grid_y_T;
  count[2] = grid_x_T;
  v1 = d_alloc_3d(count[2], count[1], count[0]);
  v2 = d_alloc_3d(count[2], count[1], count[0]);
  nc_get_vara_double(fid, ncw_var_id(fid, "x_vert_T"), start, count, v1[0][0]);
  nc_get_vara_double(fid, ncw_var_id(fid, "y_vert_T"), start, count, v2[0][0]);
  for (j = 0; j < grid_y_T; j++)
    for (i = 0; i < grid_x_T; i++) {
      x[2*j+1][2*i+1] = d1[j][i];
      y[2*j+1][2*i+1] = d2[j][i];
      x[2*j][2*i] = v1[0][j][i];
      y[2*j][2*i] = v2[0][j][i];
      x[2*j][2*i+2] = v1[1][j][i];
      y[2*j][2*i+2] = v2[1][j][i];
      x[2*j+2][2*i+2] = v1[2][j][i];
      y[2*j+2][2*i+2] = v2[2][j][i];
      x[2*j+2][2*i] = v1[3][j][i];
      y[2*j+2][2*i] = v2[3][j][i];
    }

  if (gs) {
    /* Get the locations of the last cell : nce1 edge                  */
    ni = 2 * nce1 - 2;
    mi = grid_x_C - 1;
    for (j = 0; j < grid_y_T; j++) {
      x[2*j][ni+2] = 2 * x[2*j][ni] - x[2*j][ni-2];
      y[2*j][ni+2] = 2 * y[2*j][ni] - y[2*j][ni-2];
      x[2*j+2][ni+2] = 2 * x[2*j+2][ni] - x[2*j+2][ni-2];
      y[2*j+2][ni+2] = 2 * y[2*j+2][ni] - y[2*j+2][ni-2];
    }
    /* Get the locations of the last cell : nce2 edge                  */
    nj = 2 * nce2 - 2;
    mj = grid_y_C - 1;
    for (i = 0; i < grid_x_T; i++) {
      x[nj+2][2*i] = 2 * x[nj][2*i] - x[nj-2][2*i];
      y[nj+2][2*i] = 2 * y[nj][2*i] - y[nj-2][2*1];
      x[nj+2][2*i+2] = 2 * x[nj][2*i+2] - x[nj-2][2*i+2];
      y[nj+2][2*i+2] = 2 * y[nj][2*i+2] - y[nj-2][2*1+2];
    }
  }

  /*-------------------------------------------------------------------*/
  /* Read the E cell centered geographic locations                     */
  count[0] = vertex;
  count[1] = grid_y_T;
  count[2] = grid_x_C;
  d_free_3d(v1);
  d_free_3d(v2);
  v1 = d_alloc_3d(count[2], count[1], count[0]);
  v2 = d_alloc_3d(count[2], count[1], count[0]);
  nc_get_vara_double(fid, ncw_var_id(fid, "x_vert_E"), start, count, v1[0][0]);
  nc_get_vara_double(fid, ncw_var_id(fid, "y_vert_E"), start, count, v2[0][0]);
  for (j = 0; j < grid_y_T; j++)
    for (i = 0; i < grid_x_C; i++) {
      x[2*j][2*i+1] = v1[0][j][i];
      y[2*j][2*i+1] = v2[0][j][i];
      x[2*j+2][2*i+1] = v1[3][j][i];
      y[2*j+2][2*i+1] = v2[3][j][i];
    }
  if (gs) {
    /* Get the locations of the last cell : nce2 edge                  */
    nj = 2 * nce2 - 2;
    mj = grid_y_T - 1;
    for (i = 0; i < grid_x_C; i++) {
      x[nj+2][2*i+1] = 2 * x[nj][2*i+1] - x[nj-2][2*i+1];
      y[nj+2][2*i+1] = 2 * y[nj][2*i+1] - y[nj-2][2*1+1];
    }
  }

  /*-------------------------------------------------------------------*/
  /* Read the N cell centered geographic locations                     */
  count[0] = vertex;
  count[1] = grid_y_C;
  count[2] = grid_x_T;
  d_free_3d(v1);
  d_free_3d(v2);
  v1 = d_alloc_3d(count[2], count[1], count[0]);
  v2 = d_alloc_3d(count[2], count[1], count[0]);
  nc_get_vara_double(fid, ncw_var_id(fid, "x_vert_N"), start, count, v1[0][0]);
  nc_get_vara_double(fid, ncw_var_id(fid, "y_vert_N"), start, count, v2[0][0]);
  for (j = 0; j < grid_y_C; j++)
    for (i = 0; i < grid_x_T; i++) {
      x[2*j+1][2*i] = v1[0][j][i];
      y[2*j+1][2*i] = v2[0][j][i];
      x[2*j+1][2*i+2] = v1[1][j][i];
      y[2*j+1][2*i+2] = v2[1][j][i];
    }
  if (gs) {
    /* Get the locations of the last cell : nce1 edge                  */
    ni = 2 * nce1 - 2;
    mi = grid_x_T - 1;
    for (j = 0; j < grid_y_C; j++) {
      x[2*j+1][ni+2] = 2 * x[2*j+1][ni] - x[2*j+1][ni-2];
      y[2*j+1][ni+2] = 2 * y[2*j+1][ni] - y[2*j+1][ni-2];
    }
  }

  /*-------------------------------------------------------------------*/
  /* Get the locations of the last cell : (nce1-1, nce1-1)           */
  if (gs) {
    ni = 2 * nce1 - 2;
    nj = 2 * nce2 - 2;
    x[nj+2][ni+2] = 2 * x[nj+2][ni] - x[nj+2][ni-2];
    y[nj+2][ni+2] = 2 * y[nj][ni+2] - y[nj-2][ni+2];
    ni = 2 * nce1 - 1;
    nj = 2 * nce2;
    x[nj][ni] = 0.5 * (x[nj][ni+1] + x[nj][ni-1]);
    y[nj][ni] = 0.5 * (y[nj][ni+1] + y[nj][ni-1]);
    ni = 2 * nce1;
    nj = 2 * nce2 - 1;
    x[nj][ni] = 0.5 * (x[nj+1][ni] + x[nj-1][ni]);
    y[nj][ni] = 0.5 * (y[nj+1][ni] + y[nj-1][ni]);
  }

  /*-------------------------------------------------------------------*/
  /* Get the nce1 and nce2 positions by interpolation                  */
  if (gs) {
    nj = 2 * nce2 - 1;
    ni = 2 * nce1 - 1;
    for (i = 0; i < ni; i++) {
      x[nj][i] = 0.5 * (x[nj+1][i] + x[nj-1][i]);
      y[nj][i] = 0.5 * (y[nj+1][i] + y[nj-1][i]);
    }
    nj = 2 * nce2;
    ni = 2 * nce1 - 1;
    for (j = 0; j <= nj; j++) {
      x[j][ni] = 0.5 * (x[j][ni+1] + x[j][ni-1]);
      y[j][ni] = 0.5 * (y[j][ni+1] + y[j][ni-1]);
    }
  }

  d_free_2d(d1);
  d_free_2d(d2);
  d_free_3d(v1);
  d_free_3d(v2);
  ncw_close(fname, fid);
}

/* END read_mom_grid()                                               */
/*-------------------------------------------------------------------*/

/* Locate the closest cell to oi, oj that contains valid data.
 * This is not very efficient but so what, we don't run this
 * program very often.
 */
int find_non_nan(double **a, int nx, int ny, int oi, int oj, int *ci, int *cj)
{
  int i, j;
  double mindist;
  int level;

  mindist = 1e38;
  level = 1;
  *ci = -1;
  *cj = -1;
  while (*ci < 0) {
    int finishedLevel = 0;
    int edge = 0;

    /* Scanned the whole grid, must be solid every where. */
    if (level > nx && level > ny)
      return 0;

    while (!finishedLevel) {

      int jfrom = oj - level;
      int jto = oj + level;
      int ifrom = oi - level;
      int ito = oi + level;

      switch (edge) {
      case 0:                  /* Left edge */
        ito = ifrom;
        break;

      case 1:                  /* Bottom edge */
        ++ifrom;
        jto = jfrom;
        break;

      case 2:                  /* Right edge */
        ifrom = ito;
        ++jfrom;
        break;

      case 3:                  /* Top edge */
        ++ifrom;
        jfrom = jto;
        finishedLevel = 1;
        break;
      }
      ++edge;

      for (j = jfrom; j <= jto; ++j) {
        for (i = ifrom; i <= ito; ++i) {
          if ((i < 0) || (j < 0) || (i >= nx) || (j >= ny))
            continue;

          if (!isnan(a[j][i])) {
            double dist = (oi - i) * (oi - i) + (oj - j) * (oj - j);
            if (dist < mindist) {
              *ci = i;
              *cj = j;
              mindist = dist;
            }
          }
        }
      }
    }
    ++level;
  }

  return 1;
}


/* Locate the closest cell to oi, oj that contains valid data.
 * This is not very efficient but so what, we don't run this
 * program very often.
 */
int find_non_nan_c(geometry_t *geom, double *a, int c)
{
  int i, j, cc, ci, cn;
  double mindist;
  int dx, dy;
  int level;

  mindist = 1e38;
  level = 1;
  cn = -1;
  while (cn < 0) {
    int finishedLevel = 0;
    int edge = 0;

    /* Scanned the whole grid, must be solid every where. */
    if (level > geom->nce1 || level > geom->nce2)
      return 0;

    while (!finishedLevel) {

      int jto;
      int ito;
      ci = c;

      switch (edge) {
      case 0:                  /* Left edge */
	for (cc = 0; cc < level; cc++) ci = geom->ym1[ci];
	for (cc = 0; cc < level; cc++) ci = geom->xm1[ci];	
        ito = 0;
	jto = 2 * level;
	dx = dy = -level;
        break;

      case 1:                  /* Bottom edge */
	for (cc = 0; cc < level; cc++) ci = geom->ym1[ci];
	for (cc = 0; cc < level - 1; cc++) ci = geom->xm1[ci];	
        jto = 2 * level - 1;
	dx = 1 - level;
	dy = -level;
        break;

      case 2:                  /* Right edge */
	for (cc = 0; cc < level - 1; cc++) ci = geom->ym1[ci];
	for (cc = 0; cc < level; cc++) ci = geom->xp1[ci];	
        ito = 0;
        jto = 2 * level - 1;
	dx = -level;
	dy = 1 - level;
        break;

      case 3:                  /* Top edge */
	for (cc = 0; cc < level; cc++) ci = geom->yp1[ci];
	for (cc = 0; cc < level - 1; cc++) ci = geom->xm1[ci];	
        ito = 2 * level - 2;
        jto = 0;
	dx = 1 - level;
	dy = -level;
        finishedLevel = 1;
        break;
      }
      ++edge;

      for (j = 0; j < jto; ++j) {
        for (i = 0; i < ito; ++i) {

          if (!isnan(a[c])) {
            double dist = dx * dx + dy * dy;
            if (dist < mindist) {
              cn = ci;
              mindist = dist;
            }
          }
	  if (edge == 0 || edge == 2) {
	    ci = geom->yp1[ci];
	    dy++;
	  }
	  else {
	    ci = geom->xp1[ci];
	    dx++;
	  }
        }
      }
    }
    ++level;
  }
  return(cn);
}

/*
 * For PRE_MARVL, finds name of mom grid file from the output specs
 */
int find_name_marvl(dump_data_t *dumpdata, char *ofile, char *tag, int *p_n)
{
  parameters_t *params = dumpdata->master->params;
  int n;
  
  /* Bail out if not PRE_MARVL */
  if (!(params->runmode & PRE_MARVL))
    return(0);

  /* Loop over all output files and find tag */
  for (n=0; n<params->ndf; n++)
    if (strcmp(params->d_filetype[n], tag) == 0) {
      sprintf(ofile, "%s", params->d_name[n]);
      /* Optionally return the index of this dump file */
      if (p_n != NULL)
	*p_n = n; 
      return(1);
    }

  /* not found */
  return(0);
}
