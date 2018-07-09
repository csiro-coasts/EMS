/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/inputs/readdump.c
 *  
 *  Description:
 *  Routine to read model netCDF dump file
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: readdump.c 5873 2018-07-06 07:23:48Z riz008 $
 *
 */

#include <string.h>
#include <netcdf.h>
#include "hd.h"
#include "tracer.h"


int dumpdata_read_2d_us(dump_data_t *dumpdata, int id, char *name,
			double *p, int dump, int ni, int oset);
int dumpdata_read_3d_us(dump_data_t *dumpdata, int id, char *name, double *p, 
			double **d, int dump, int ni, int ni3, int *map, int **rmap, int oset);
int dumpdata_read_2d0_us(dump_data_t *dumpdata, int id, char *name,
			 double *p, int ni, int oset);
int netCDF_grid_us(parameters_t *params, int cdfid, int *c2cc);
double smoothbathy(double **a, int c, int p, int m, int pp, int mm, int cc,
                   int mode, double bmin);
void intp_undef_s(geometry_t *geom, double *a);
void mean_z_bathy(int i, int j, int zmfe1, int zmfe2, double **bathy);
void set_z_bathy(int i, int j, int zmfe1, int zmfe2, double **bathy, double val);
void adjust_outside_obc(parameters_t * params, double **bathy);
void read_sed_layers(geometry_t *geom, parameters_t *params, master_t *master, dump_data_t *dumpdata);
void remove_channel(parameters_t *params, unsigned long ***flag, double **bathy);
double smooth_bathy(unsigned long ***flag, int nz, int i, int j, int **xm, int **xp, 
		    int **ym, int **yp, double **bathy, double bmin, double bmax);
double bathy_smooth_us(mesh_t *m, double *bathy, int cc);



/*-------------------------------------------------------------------*/
/* Routine to set variables at ghost cells                           */
/*-------------------------------------------------------------------*/
void master_setghosts(geometry_t *geom, /* Sparse global geometry    */
                      master_t *master, /* Master data               */
		      double **cellx,   /* Cell centre x location    */
		      double **celly,   /* Cell centre y location    */
		      double **gridx,   /* Cell corner x location    */
		      double **gridy,   /* Cell corner y location    */
		      int *s2i,         /* Optional sparse to i map  */
		      int *s2j          /* Optional sparse to j map  */
  )
{
  int cb, ci, c, cc, ee, e, ei; /* Sparse couters                    */
  int n;
  double mx = 0;
  double *eta;
  GRID_SPECS *gs;
  int *mask;

  /* Set a no-gradient condition over lateral boundaries             */
  for (ee = 1; ee <= geom->nbpte1S; ee++) {
    cb = geom->bpte1[ee];
    ci = geom->bine1[ee];
    geom->h1au1[cb] = geom->h1au1[ci];
    geom->h1acell[cb] = geom->h1acell[ci];
    geom->h2au1[cb] = geom->h2au1[ci];
    /*if(geom->h2au1[cb]==0.0)printf("h2au1 zero @ %d %d\n", cb, ci);*/
  }

  /* Coriolis and bottom depth                                       */
  for (cc = 1; cc <= geom->nbptS; cc++) {
    cb = geom->bpt[cc];
    ci = geom->bin[cc];
    geom->botz[cb] = geom->botz[ci];
    if (master != NULL)
      master->coriolis[cb] = master->coriolis[ci];
  }

  /* Set the bottom depth on e1 and e2 open boundaries (no gradient) */
  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    for (ee = 1; ee <= open->no2_e1; ee++) {
      e = open->obc_e1[ee];
      ci = open->oi1_e1[ee];
      /*geom->h1au1[e] = geom->h1au1[ci];*/
      geom->h2au1[e] = geom->h2au1[ci];
      /*geom->thetau1[e] = geom->thetau1[ci];*/
      geom->h1acell[e] = geom->h1acell[ci];
    }

    /* Boundary length                                               */
    open->length = 0.0;
    for (ee = 1; ee <= open->no2_e1; ee++) {
      e = open->obc_e1[ee];
      open->length += geom->h1au1[e];
    }

    /* Boundary area                                                 */
    open->area = 0.0;
    for (cc = 1; cc <= open->no3_t; cc++) {
      int cs, zp1;
      double dz;
      c = open->obc_t[cc];
      cs = geom->m2d[c];
      zp1 = geom->zp1[c];
      dz = min(master->eta[cs], geom->gridz[zp1]) -
	max(geom->botz[cs], geom->gridz[c]);
      open->area += dz * geom->h1au1[cs];
    }

    /* Boundary depths                                               */
    open->meandep = open->maxdep = 0.0;
    open->mindep = -1e10;
    for (cc = 1; cc <= open->no2_t; cc++) {
      int cs;
      c = open->obc_t[cc];
      open->maxdep = min(open->maxdep, geom->botz[c]);
      open->mindep = max(open->mindep, geom->botz[c]);
      open->meandep += geom->botz[c];
    }
    open->meandep /= (double)open->no2_t;

    /* Cell location and bottom depth                                */
    for (ee = 1; ee <= open->no2_e1; ee++) {
      int i, j, bi, bj, m, m1, m2;
      c = ci = open->obc_e2[ee];
      e = ei = open->obc_e1[ee];
      for (m = 0; m < open->bgz; m++) {
	m1 = c;
	m2 = open->nmape[ee][c];
	c = open->omape[ee][c];
	e = geom->c2e[open->outi[ee]][c];
	if (open->ceni[ee]) {
	  if(isnan(geom->cellx[c]) && isnan(geom->celly[c])) {
	    if(!isnan(geom->cellx[m1]) && !isnan(geom->cellx[m2]))
	      geom->cellx[c] = intp(geom->cellx[m2], geom->cellx[m1], 0, 1, 2);
	    if(!isnan(geom->celly[m1]) && !isnan(geom->celly[m2]))
	      geom->celly[c] = intp(geom->celly[m2], geom->celly[m1], 0, 1, 2);
	  } else {
	    geom->cellx[c] = geom->cellx[m1];
	    geom->celly[c] = geom->celly[m1];
	  }
	}
	if (!open->ceni[ee]) {
	  if(isnan(geom->cellx[c]) && isnan(geom->celly[c])) {
	    if(!isnan(geom->cellx[m1]) && !isnan(geom->cellx[m2]))
	      geom->cellx[c] = intp(geom->cellx[m1], geom->cellx[m2], 1, 2, 0);
	    if(!isnan(geom->celly[m1]) && !isnan(geom->celly[m2]))
	      geom->celly[c] = intp(geom->celly[m1], geom->celly[m2], 1, 2, 0);
	    else {
	      geom->cellx[c] = geom->cellx[m1];
	      geom->celly[c] = geom->celly[m1];
	    }
	  }
	}
	/* Uniform depth and metrics over the OBC ghosts */
	geom->botz[c] = geom->botz[ci];
	geom->h1au1[e] = geom->h1au1[ei];
	geom->h2au1[e] = geom->h2au1[ei];
	geom->thetau1[e] = geom->thetau1[ei];
	geom->h1acell[e] = geom->h1acell[ei];
      }
    }
    for (ee = 1; ee <= open->no3_e1; ee++) {
      int m;
      c = ci = open->obc_e2[ee];
      e = open->obc_e1[ee];
      for (m = 0; m < open->bgz; m++) {
	c = open->omape[ee][c];
	geom->gridz[c] = geom->gridz[ci];
	geom->cellz[c] = geom->cellz[ci];
      }
      if (!open->bgz) {
	c = open->ogc_t[ee];	
	geom->gridz[c] = geom->gridz[ci];
	geom->cellz[c] = geom->cellz[ci];
      }
    }
  }

  /* Median filter eta if required                                   */
  eta = d_alloc_1d(geom->sgsizS);
  memcpy(eta, master->eta, geom->sgsizS * sizeof(double));
  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    for (ee = 1; ee <= open->no2_e1; ee++) {
      c = open->obc_e2[ee];
      e = open->obc_e1[ee];
      eta[open->omape[ee][c]] = eta[open->nmape[ee][c]];
    }
  }
  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    if (open->options & OP_ETAFIL) {
      for (cc = 1; cc <= open->no2_t; cc++) {
	c = open->obc_t[cc];
	master->eta[c] = con_median(geom, eta, 10.0, c, ST_SQ3);
      }
    }
  }
  d_free_1d(eta);

  /* Set the lateral boundaries again to pick up cells defined on    */
  /* OBCs.                                                           */
  for (cc = 1; cc <= geom->nbptS; cc++) {
    cb = geom->bpt[cc];
    ci = geom->bin[cc];
    geom->botz[cb] = geom->botz[ci];
    if (master != NULL)
      master->coriolis[cb] = master->coriolis[ci];
  }
  for (cc = 1; cc <= geom->nbpte1S; cc++) {
    cb = geom->bpte1[cc];
    ci = geom->bine1[cc];
    geom->h1au1[cb] = geom->h1au1[ci];
    geom->h1acell[cb] = geom->h1acell[ci];
    geom->h2au1[cb] = geom->h2au1[ci];
  }

  /* Structured grids only                                           */
  if (geom->us_type & US_IJ) {
    c = geom->map[geom->nz-1][geom->nce2][geom->nce1];
    ci = geom->c2c[1][geom->c2c[4][c]];
    if (geom->h1au1[c] == 0.0) geom->h1au1[c] = geom->h1au1[ci];

    /* Cell location */
    for (cc = 1; cc <= geom->nbptS; cc++) {
      int i, j, ig, jg;
      c = geom->bpt[cc];
      ci = geom->bin[cc];
      i = ig = geom->s2i[ci];
      j = jg = geom->s2j[ci];
      
      if (geom->c2c[3][ci] == c) ig = i + 1;
      else if (geom->c2c[1][ci] == c) ig = i - 1;
      else if (geom->c2c[2][ci] == c) jg = j + 1;
      else if (geom->c2c[4][ci] == c) jg = j - 1;
      else if (geom->c2c[2][geom->c2c[3][ci]] == c || 
	       geom->c2c[3][geom->c2c[2][ci]] == c) {
	ig = i + 1;
	jg = j + 1;
      } else if (geom->c2c[2][geom->c2c[1][ci]] == c || 
		 geom->c2c[1][geom->c2c[2][ci]] == c) {
	ig = i - 1;
	jg = j + 1;
      } else if (geom->c2c[4][geom->c2c[3][ci]] == c || 
		 geom->c2c[3][geom->c2c[4][ci]] == c) {
	ig = i + 1;
	jg = j - 1;
      } else if (geom->c2c[4][geom->c2c[1][ci]] == c || 
		 geom->c2c[1][geom->c2c[4][ci]] == c) {
	ig = i - 1;
	jg = j - 1;
      }

      ig = min(max(ig, 0), geom->nce1-1);
      jg = min(max(jg, 0), geom->nce2-1);

      geom->cellx[c] = cellx[jg][ig];
      geom->celly[c] = celly[jg][ig];
      geom->gridx[geom->c2v[0][c]] = gridx[jg][ig];
      geom->gridy[geom->c2v[0][c]] = gridy[jg][ig];
      if (s2i != NULL) {
	s2i[c] = ig;
	cellx[jg][ig] = cellx[j][i];
	gridx[jg][ig] = gridx[j][i];
      }
      if (s2j != NULL) {
	s2j[c] = jg;
	celly[jg][ig] = celly[j][i];
	gridy[jg][ig] = gridy[j][i];
      }
    }
  }

  intp_undef_s(geom, geom->cellx);
  intp_undef_s(geom, geom->celly);
  /*
  intp_undef_s(geom, geom->gridx);
  intp_undef_s(geom, geom->gridy);
  */

  /* Set the geographic locations of ghost cells.                    */
  /* Lateral ghost cells.                                            */
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
	ci = geom->c2c[j][c];
	geom->cellx[ci] = geom->u1x[e] + d * cos(r);
	geom->celly[ci] = geom->u1y[e] + d * sin(r);
      }
    }    
  }
  /* Lateral ghost cells again; for safety use work arrays           */
  for(cc = geom->b2_t+1; cc <= geom->n2_t; cc++) {
    double x, y, d, r;
    int j;
    c = geom->w2_t[cc];
    ci = geom->wgst[c];
    for (j = 1; j <= geom->npe[c]; j++) {
      if (c == geom->c2c[j][ci]) {
	e = geom->c2e[j][ci];
	x = geom->u1x[e] - geom->cellx[ci];
	y = geom->u1y[e] - geom->celly[ci];
	d = sqrt(x * x + y * y);
	r = atan2(y, x);
	geom->cellx[c] = geom->u1x[e] + d * cos(r);
	geom->celly[c] = geom->u1y[e] + d * sin(r);
      }
    }
  }
  /* Open boundary ghost cells                                       */
  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    for (ee = 1; ee <= open->no2_e1; ee++) {
      double x, y, d, r;
      int i, j;
      e = open->obc_e1[ee];
      c = open->ogc_t[ee];
      ci = i = open->obc_e2[ee];
      for (j = 1; j <= geom->npe[ci]; j++) {
	if (e == geom->c2e[j][ci]) {
	  x = geom->u1x[e] - geom->cellx[ci];
	  y = geom->u1y[e] - geom->celly[ci];
	  d = 2.0 * sqrt(x * x + y * y);
	  r = atan2(y, x);
	  do {
	    geom->cellx[c] = geom->cellx[ci] + d * cos(r);
	    geom->celly[c] = geom->celly[ci] + d * sin(r);
	    ci = c;
	    c = open->omape[ee][c];
	  } while (c != ci);
	}
      }
    }
  }

  /* Get the surface coordinate corresponding to the maximum bottom */
  /* depth.  */
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    if (geom->botz[c] < mx) {
      geom->cmx = cc;
      mx = geom->botz[c];
    }
  }
}

/* END master_setghosts()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to interpolate grid location at undefined cells           */
/*-------------------------------------------------------------------*/
void intp_undef_s(geometry_t *geom, double *a)
{
  int c, j, npe;
  double *aa = d_alloc_1d(geom->szcS);
  memcpy(aa, a, geom->szcS * sizeof(double));

  if (geom->us_type & US_HEX) {
    for (c = 1; c < geom->szcS; c++) {
      int m1, m2;
      if (isnan(a[c])) {
	npe = geom->npe[c];
	for (j = 1; j <= npe; j++) {
	  m1 = geom->c2c[j][c];
	  m2 = geom->c2c[j][m1];
	  if (!isnan(aa[m1]) && !isnan(aa[m2])) {
	    a[c] = intp(a[m2], a[m1], 0, 1, 2);
	    break;
	  }
	}
      }
    }
    d_free_1d(aa);
  }
}

/* END intp_undef_s()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads in the unstructured input file                              */
/*-------------------------------------------------------------------*/
int dumpdata_read_us(geometry_t *geom, /* Sparse global geometry structure */
		     parameters_t *params, /* Input parameter data structure */
		     master_t *master, /* Master data structure */
		     dump_data_t *dumpdata,  /* Dump data structure */
		     int cdfid, int ti)
{
  int c, cc, cs, e, ee, i, j, k, n, vv, v;
  size_t ndumps;
  size_t start[4];
  size_t count[4];
  size_t nMesh2_node;
  size_t nMesh2_edge;
  size_t nMesh2_face;
  int **k2e, **k2c, **k2v;
  int oset = params->oset;

  /*-----------------------------------------------------------------*/
  /* Dimensions and initialise                                       */
  nc_inq_dimlen(cdfid, ncw_dim_id(cdfid, "nMesh2_face"), &nMesh2_face);
  nc_inq_dimlen(cdfid, ncw_dim_id(cdfid, "nMesh2_edge"), &nMesh2_edge);
  nc_inq_dimlen(cdfid, ncw_dim_id(cdfid, "nMesh2_node"), &nMesh2_node);
  master->dumpdata = dumpdata;

  master->i1 = NULL;
  master->i2 = NULL;
  /* Check that the requested dump exists */
  nc_inq_dimlen(cdfid, ncw_dim_id(cdfid, "record"), &ndumps);
  if (ti > (int)(ndumps - 1)) {
    warn("dump_read: dump %d does not exist in file!\n", ti);
    return (0);
  }

  /*-----------------------------------------------------------------*/
  /* Make the reverse maps                                           */
  k2c = i_alloc_2d(geom->szcS, geom->nz);
  for (cc = 1; cc <= geom->b3_t; cc++) {
    c = geom->w3_t[cc];
    k = geom->s2k[c];
    cs = geom->m2d[c];
    k2c[k][cs-params->oset] = c;
  }
  k2e = i_alloc_2d(geom->szeS, geom->nz);
  for (ee = 1; ee <= geom->n3_e1; ee++) {
    e = geom->w3_e1[ee];
    k = geom->e2k[e];
    cs = geom->m2de[e];
    k2e[k][cs-params->oset] = e;
  }

  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;

  /*-----------------------------------------------------------------*/
  /* Time independent variables                                      */
  count[0] = dumpdata->nz;
  count[1] = 0;
  count[2] = 0;
  count[3] = 0;

  /* Vertical. params->layers is populated from the input netCDF     */
  /* in dump_open_us(), and dumpdata->gridz and dumpdata->cellz are  */
  /* set in dumpdata_build(). Copy these layers structures to geom   */
  /* in 3D.                                                          */
  for (cc = 1; cc <= geom->n3_t; cc++) {
    c = geom->w3_t[cc];
    k = geom->s2k[c];
    geom->gridz[c] = dumpdata->gridz[k];
    geom->cellz[c] = dumpdata->cellz[k];
  }

  /* Horizontal. Note: the following geometry and metric arrays are  */
  /* assigned in build_sparse_grid_us(). They may be re-read from    */
  /* the input file here, either for comparison or overwrite.        */
  /*
  dumpdata_read_2d0_us(dumpdata, cdfid, "Mesh2_node_x", geom->gridx, nMesh2_node);
  dumpdata_read_2d0_us(dumpdata, cdfid, "Mesh2_node_y", geom->gridy, nMesh2_node);
  dumpdata_read_2d0_us(dumpdata, cdfid, "Mesh2_face_x", geom->cellx, nMesh2_face);
  dumpdata_read_2d0_us(dumpdata, cdfid, "Mesh2_face_y", geom->celly, nMesh2_face);
  dumpdata_read_2d0_us(dumpdata, cdfid, "Mesh2_edge_x", geom->u1x, nMesh2_edge);
  dumpdata_read_2d0_us(dumpdata, cdfid, "Mesh2_edge_y", geom->u1y, nMesh2_edge);
  dumpdata_read_2d0_us(dumpdata, cdfid, "h1acell", geom->h1acell, nMesh2_edge);
  dumpdata_read_2d0_us(dumpdata, cdfid, "h1au1", geom->h1au1, nMesh2_edge);
  dumpdata_read_2d0_us(dumpdata, cdfid, "h2au1", geom->h2au1, nMesh2_edge);
  dumpdata_read_2d0_us(dumpdata, cdfid, "thetau1", geom->thetau1, nMesh2_edge);
  */
  dumpdata_read_2d0_us(dumpdata, cdfid, "coriolis", master->coriolis, nMesh2_face, oset);

  /*
  count[0] = dumpdata->nvertex2;
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "Mesh2_node_x"), start, count,
                     geom->gridx);
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "Mesh2_node_y"), start, count,
                     geom->gridy);
  count[0] = dumpdata->nface2;
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "Mesh2_face_x"), start, count,
                     geom->cellx);
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "Mesh2_face_y"), start, count,
                     geom->celly);
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "coriolis"), start, count,
                     master->coriolis);
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "botz"), start, count,
   	             geom->botz);

  count[0] = dumpdata->nedge2;
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "h1acell"), start, count,
                     geom->h1acell);
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "h1au1"), start, count,
                     geom->h1au1);
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "h2au1"), start, count,
                     geom->h2au1);
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "thetau1"), start, count,
                     geom->thetau1);
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "Mesh2_edge_x"), start, count,
                     geom->u1x);
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "Mesh2_edge_y"), start, count,
                     geom->u1y);
  */

  /*-----------------------------------------------------------------*/
  /* Time dependent variables                                        */
  start[0] = ti;
  count[0] = 1L;
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "t"), start, count,
                     &dumpdata->t);

  if (params->timeunit) {
    char timeunits[MAXSTRLEN];
    memset(timeunits, 0, MAXSTRLEN);
    nc_get_att_text(cdfid, ncw_var_id(cdfid, "t"), "units", timeunits);
    tm_change_time_units(timeunits, params->timeunit, &dumpdata->t, 1);
  }
  params->t = dumpdata->t;

  /* Note: use the dumpdata->u1av array to read the backward vel */
  dumpdata_read_2d_us(dumpdata, cdfid, "u1avb", master->u1avb, ti, nMesh2_edge, oset);
  dumpdata_read_2d_us(dumpdata, cdfid, "u1av", master->u1av, ti, nMesh2_edge, oset);
  dumpdata_read_2d_us(dumpdata, cdfid, "u1bot", master->u1bot, ti, nMesh2_edge, oset);
  dumpdata_read_2d_us(dumpdata, cdfid, "wind1", master->wind1, ti, nMesh2_edge, oset);
  dumpdata_read_2d_us(dumpdata, cdfid, "wtop", master->wtop, ti, nMesh2_face, oset);
  dumpdata_read_2d_us(dumpdata, cdfid, "topz", master->topz, ti, nMesh2_face, oset);
  dumpdata_read_2d_us(dumpdata, cdfid, "eta", master->eta, ti, nMesh2_face, oset);
  dumpdata_read_2d_us(dumpdata, cdfid, "patm", master->patm, ti, nMesh2_face, oset);
  dumpdata_read_2d_us(dumpdata, cdfid, "Cd", master->Cd, ti, nMesh2_face, oset);

  dumpdata_read_3d_us(dumpdata, cdfid, "u1b", master->u1b, dumpdata->we, ti, nMesh2_edge, 
		      geom->sze, geom->w3_e1, k2e, oset);
  dumpdata_read_3d_us(dumpdata, cdfid, "u1", master->u1, dumpdata->we, ti, nMesh2_edge, 
		      geom->sze, geom->w3_e1, k2e, oset);
  dumpdata_read_3d_us(dumpdata, cdfid, "w", master->w, dumpdata->wc, ti, nMesh2_face, 
		      geom->szc, geom->w3_t, k2c, oset);
  dumpdata_read_3d_us(dumpdata, cdfid, "dens", master->dens, dumpdata->wc, ti, nMesh2_face, 
		      geom->szc, geom->w3_t, k2c, oset);
  dumpdata_read_3d_us(dumpdata, cdfid, "dens_0", master->dens_0, dumpdata->wc, ti, nMesh2_face, 
		      geom->szc, geom->w3_t, k2c, oset);
  dumpdata_read_3d_us(dumpdata, cdfid, "Kz", master->Kz, dumpdata->wc, ti, nMesh2_face, 
		      geom->szc, geom->w3_t, k2c, oset);
  dumpdata_read_3d_us(dumpdata, cdfid, "Vz", master->Vz, dumpdata->wc, ti, nMesh2_face, 
		      geom->szc, geom->w3_t, k2c, oset);

  /*-----------------------------------------------------------------*/
  /* Check for NaNs and set ghost cells                              */
  for (i = 1; i <= geom->b2_t; i++) {
    int ci, cj;
    c = geom->w2_t[i];
    if (isnan(master->eta[c])) {
      /*
      find_closest_nonnan(dumpdata, dumpdata->eta, geom->s2i[c], geom->s2j[c],
			  geom->nz-1, &ci, &cj);
      */
      master->eta[c] = dumpdata->eta[cj][ci];
      hd_warn("eta_init: Found NaN at (%d %d); replaced with (%d %d).\n", 
	      geom->s2i[c], geom->s2j[c], ci, cj);
    }
  }
  set_lateral_bc_eta(master->eta, geom->nbptS, geom->bpt, geom->bin,
		     geom->bin2, 1);

  /*-----------------------------------------------------------------*/
  /* 3D tracers                                                      */
  for (n = 0; n < dumpdata->ntr; n++) {
    /* If a 3D tracer is included in the input file then initialise  */
    /* with this data. If not use the specification in the parameter */
    /* file tracer list; TRACER#.data                                */
    if (dumpdata_read_3d_us(dumpdata, cdfid, dumpdata->trinfo_3d[n].name, master->tr_wc[n],
			 dumpdata->wc, ti, nMesh2_face, geom->szc, geom->w3_t, k2c, oset)) {
      load_wc_tracer_name(master, params->prmfd, dumpdata->trinfo_3d[n].name, WATER);
    }

    /* Ignore DA variables */
    if (strncasecmp(master->trinfo_3d[n].tag, "DA_OBS", 6) == 0 ) continue;

    /* Check for NaN's in the initial condition */
    for (i = 1; i <= geom->b3_t; i++) {
      c = geom->w3_t[i];
      if (isnan(master->tr_wc[n][c]))
	hd_warn
	  ("dumpdata_read: NaN found in tracer '%s' at (%d %d %d).\n",
	   master->trinfo_3d[n].name, geom->s2i[c], geom->s2j[c], geom->s2k[c]);
      master->tr_wc[n][c] = min(master->tr_wc[n][c], master->trinfo_3d[n].valid_range_wc[1]);
    }
  }

  /*-----------------------------------------------------------------*/
  /* 2D tracers                                                      */
  for (n = 0; n < dumpdata->ntrS; n++) {
    if (dumpdata_read_2d_us(dumpdata, cdfid, dumpdata->trinfo_2d[n].name,
			    master->tr_wcS[n], ti, nMesh2_face, oset)) {
      load_wc_tracer_name(master, params->prmfd, dumpdata->trinfo_2d[n].name, INTER);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Sediment tracers                                                */
  if (dumpdata->sednz > 0) {
    int flg = 0;
    /*
    if (dumpdata_read_3d(dumpdata, cdfid, "z_grid_sed", dumpdata->gridz_sed,
			 ti, dumpdata->sednz + 1, dumpdata->nce2,
			 dumpdata->nce1)) flg = 1;
    for (k = 0; k < params->sednz + 1; k++) {
      c2s_2d(geom, geom->gridz_sed[k], dumpdata->gridz_sed[k], geom->nce1,
             geom->nce2);
    }
    if (dumpdata_read_3d(dumpdata, cdfid, "z_centre_sed", dumpdata->cellz_sed,
			 ti, dumpdata->sednz, dumpdata->nce2, dumpdata->nce1))
      flg = 1;
    for (k = 0; k < params->sednz; k++) {
      c2s_2d(geom, geom->cellz_sed[k], dumpdata->cellz_sed[k], geom->nce1,
             geom->nce2);
    }
    */
    /* If sediment structure does not exist in the input file, then  */
    /* use the structure defined in the parameter file.              */
    if (flg) read_sed_layers(geom, params, master, dumpdata);
    /*
    for (n = 0; n < dumpdata->nsed; n++) {
      char name[MAXSTRLEN];
      strcpy(name, dumpdata->trinfo_sed[n].name);
      strcat(name, "_sed");
      if (dumpdata_read_3d(dumpdata, cdfid, name, dumpdata->tr_sed[n], ti,
			   dumpdata->sednz, dumpdata->nce2, dumpdata->nce1)) {
	load_wc_tracer_name(master, params->prmfd, dumpdata->trinfo_sed[n].name, SEDIM);
      } else {
	for (k = 0; k < params->sednz; k++)
	  c2s_2d(geom, master->tr_sed[n][k], dumpdata->tr_sed[n][k],
		 geom->nce1, geom->nce2);
      }
    }
    */
  }

  master->t = params->t;  /* Required to initialize dumpdata->t     */
  dumpdata_init(dumpdata, geom, master);

  i_free_2d(k2c);
  i_free_2d(k2e);
  return (1);
}

/* END dumpdata_read_us()                                            */
/*-------------------------------------------------------------------*/


/*-----------------------------------------------------------------------*/
/* Reads the time dependent input data from netCDF file                  */
/*-----------------------------------------------------------------------*/
int dump_re_read(master_t *master, /* Master data structure */
		 int cdfid, int ti)
{
  geometry_t *geom = master->geom;
  dump_data_t *dumpdata = master->dumpdata;
  int c, cc, i, j, k, n;
  size_t ndumps;
  size_t start[4];
  size_t count[4];
  size_t nMesh2_node;
  size_t nMesh2_edge;
  size_t nMesh2_face;
  int **k2e, **k2c, **k2v;
  int oset;
  int time_vid;
  char ftimeunits[MAXSTRLEN];
  memset(ftimeunits,0,MAXSTRLEN);
  oset = (dumpdata->start_index) ? 0 : 1;

  /* Check that the requested dump exists */
  nc_inq_dimlen(cdfid, ncw_dim_id(cdfid, "record"), &ndumps);
  if (ti > (int)(ndumps - 1)) {
    hd_warn("dump_read: dump %d does not exist in file!\n", ti);
    return (0);
  }

  if (DEBUG("dump"))
    dlog("dump", "Reading dump %d.\n", ti);

  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;
  /* Time independent variables */
  count[0] = geom->nz;
  count[1] = 0;
  count[2] = 0;
  count[3] = 0;

  /*---------------------------------------------------------------------*/
  /* Time dependent variables                                            */
  start[0] = ti;
  count[0] = 1L;
  time_vid = ncw_var_id(cdfid, "t");
  nc_get_vara_double(cdfid, time_vid, start, count, &master->t);

  /* Make sure time value is consistent with the time units */
  nc_get_att_text(cdfid, time_vid, "units", ftimeunits);
  if (strcmp(ftimeunits, master->timeunit) != 0)
    tm_change_time_units(ftimeunits, master->timeunit, &master->t, 1);
  
  /* Read the 2 dimensional input data into the master sparse array */
  dumpdata_read_2d_us(dumpdata, cdfid, "u1avb", master->u1avb, ti, nMesh2_edge, oset);
  dumpdata_read_2d_us(dumpdata, cdfid, "u1av", master->u1av, ti, nMesh2_edge, oset);
  dumpdata_read_2d_us(dumpdata, cdfid, "u1bot", master->u1bot, ti, nMesh2_edge, oset);
  dumpdata_read_2d_us(dumpdata, cdfid, "wind1", master->wind1, ti, nMesh2_edge, oset);
  dumpdata_read_2d_us(dumpdata, cdfid, "wtop", master->wtop, ti, nMesh2_face, oset);
  dumpdata_read_2d_us(dumpdata, cdfid, "topz", master->topz, ti, nMesh2_face, oset);
  dumpdata_read_2d_us(dumpdata, cdfid, "eta", master->eta, ti, nMesh2_face, oset);
  dumpdata_read_2d_us(dumpdata, cdfid, "patm", master->patm, ti, nMesh2_face, oset);
  dumpdata_read_2d_us(dumpdata, cdfid, "Cd", master->Cd, ti, nMesh2_face, oset);

  dumpdata_read_3d_us(dumpdata, cdfid, "u1b", master->u1b, dumpdata->we, ti, nMesh2_edge, 
		      geom->sze, geom->w3_e1, k2e, oset);
  dumpdata_read_3d_us(dumpdata, cdfid, "u1", master->u1, dumpdata->we, ti, nMesh2_edge, 
		      geom->sze, geom->w3_e1, k2e, oset);
  dumpdata_read_3d_us(dumpdata, cdfid, "w", master->w, dumpdata->wc, ti, nMesh2_face, 
		      geom->szc, geom->w3_t, k2c, oset);
  dumpdata_read_3d_us(dumpdata, cdfid, "dens", master->dens, dumpdata->wc, ti, nMesh2_face, 
		      geom->szc, geom->w3_t, k2c, oset);
  dumpdata_read_3d_us(dumpdata, cdfid, "dens_0", master->dens_0, dumpdata->wc, ti, nMesh2_face, 
		      geom->szc, geom->w3_t, k2c, oset);
  dumpdata_read_3d_us(dumpdata, cdfid, "Kz", master->Kz, dumpdata->wc, ti, nMesh2_face, 
		      geom->szc, geom->w3_t, k2c, oset);
  dumpdata_read_3d_us(dumpdata, cdfid, "Vz", master->Vz, dumpdata->wc, ti, nMesh2_face, 
		      geom->szc, geom->w3_t, k2c, oset);

  /* Check for NaNs                                                  */
  for (cc = 1; cc <= geom->b2_t; cc++) {
    int ci, cj;
    c = geom->w2_t[cc];
    if (isnan(master->eta[c])) {
      hd_warn("dump_read: Found NaN at (%d %d).\n", geom->s2i[c], geom->s2j[c]);
    }
  }
  set_lateral_bc_eta(master->eta, geom->nbptS, geom->bpt, geom->bin,
		     geom->bin2, 1);

  /*-----------------------------------------------------------------*/
  /* 3D tracers                                                      */
  for (n = 0; n < master->ntr; n++) {
    /* Ignore DA variables */
    if (strncmp(master->trinfo_3d[n].tag, "DA_", 3) == 0 ) continue;

    /* If a 3D tracer is included in the input file then initialise  */
    /* with this data. If not use the specification in the parameter */
    /* file tracer list; TRACER#.data                                */
    dumpdata_read_3d_us(dumpdata, cdfid, dumpdata->trinfo_3d[n].name, master->tr_wc[n],
			dumpdata->wc, ti, nMesh2_face, geom->szc, geom->w3_t, k2c, oset);

    /* Check for NaN's in the initial condition */
    for (i = 1; i <= geom->b3_t; i++) {
      c = geom->w3_t[i];
      if (isnan(master->tr_wc[n][c]))
	hd_warn
	  ("dumpdata_read: NaN found in tracer '%s' at (%d %d %d).\n",
	   master->trinfo_3d[n].name, geom->s2i[c], geom->s2j[c], geom->s2k[c]);
    }
  }

  /*-----------------------------------------------------------------*/
  /* 2D tracers                                                      */
  for (n = 0; n < master->ntrS; n++) {
    dumpdata_read_2d_us(dumpdata, cdfid, dumpdata->trinfo_2d[n].name,
			master->tr_wcS[n], ti, nMesh2_face, oset);
  }

  /* Sediment thickness */
  if (geom->sednz > 0) {
    double ***d1 = d_alloc_3d(dumpdata->nce1, dumpdata->nce2, dumpdata->sednz);
    for (n = 0; n < dumpdata->nsed; n++) {
      char name[MAXSTRLEN];
      strcpy(name, dumpdata->trinfo_sed[n].name);
      strcat(name, "_sed");
      /*
      dumpdata_read_3d_us(dumpdata, cdfid, name, d1, ti,
			  dumpdata->sednz, dumpdata->nce2, dumpdata->nce1);
      */
    }
    d_free_3d(d1);
  }

  return (1);
}

/* END dump_re_read()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads in the unstructured grid topology from netCDF               */
/*-------------------------------------------------------------------*/
int netCDF_grid_us(parameters_t *params, int cdfid, int *c2cc)
{
  int c, cc, n;
  size_t start[4];
  size_t count[4];
  int *s2i, *s2j, i1f, i2f;
  double *x, *y;
  double **xbnds, **ybnds;
  size_t ns2, npe;

  nc_inq_dimlen(cdfid, ncw_dim_id(cdfid, "nMesh2_face"), &ns2);
  nc_inq_dimlen(cdfid, ncw_dim_id(cdfid, "nMaxMesh2_face_nodes"), &npe);
  params->npe = (int)npe-1;
  params->ns2 = (int)ns2-1;

  s2i = i_alloc_1d(params->ns2+1);
  s2j = i_alloc_1d(params->ns2+1);
  params->x = d_alloc_2d(params->npe+1, params->ns2+1);
  params->y = d_alloc_2d(params->npe+1, params->ns2+1);
  x = d_alloc_1d(params->ns2 + 1);
  y = d_alloc_1d(params->ns2 + 1);

  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;
  /* time independent variables */
  count[0] = params->ns2+1;
  count[1] = 0;
  count[2] = 0;
  count[3] = 0;

  /* Get the cell centre coordinates */
  if ((i1f = ncw_var_id(cdfid, "Mesh2_iindex")) >= 0)
    nc_get_vara_int(cdfid, ncw_var_id(cdfid, "Mesh2_iindex"), start, count, s2i);  
  if ((i2f = ncw_var_id(cdfid, "Mesh2_jindex")) >= 0)
    nc_get_vara_int(cdfid, ncw_var_id(cdfid, "Mesh2_jindex"), start, count, s2j);
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "Mesh2_face_x"), start, count, x);
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "Mesh2_face_y"), start, count, y);

  params->ic = params->jc = NULL;
  if (i1f >= 0 && i2f >= 0) {
    params->us_type |= US_IJ;
    params->ic = i_alloc_1d(params->ns2+1);
    params->jc = i_alloc_1d(params->ns2+1);
  }
  for (c = 1; c <= params->ns2; c++) {
    cc = c2cc[c];
    params->x[cc][0] = x[c];
    params->y[cc][0] = y[c];
    if (i1f >= 0 && i2f >= 0) {
      params->ic[cc] = s2i[c];
      params->jc[cc] = s2j[c];
    }
  }
  d_free_1d(x);
  d_free_1d(y);

  /* Get the cell vertex coordinates */
  xbnds = d_alloc_2d(params->ns2+1, params->npe+1);
  ybnds = d_alloc_2d(params->ns2+1, params->npe+1);
  count[0] = params->npe+1;
  count[1] = params->ns2+1;
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "Mesh2_face_xbnds"), start, count,
                     xbnds[0]);
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "Mesh2_face_ybnds"), start, count,
                     ybnds[0]);

  for (c = 1; c <= params->ns2; c++) {
    cc = c2cc[c];
    for (n = 1; n <= params->npe; n++) {
      params->x[cc][n] = xbnds[n][c];
      params->y[cc][n] = ybnds[n][c];
    }
  }
  d_free_2d(xbnds);
  d_free_2d(ybnds);
  i_free_1d(s2i);
  i_free_1d(s2j);
}

/* END netCDF_grid_us()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read the bathymetry from the input netCDF file.        */
/* Note: there is generally a difference between the listing of      */
/* bathymetry in an input text file, and the reading of bathymetry   */
/* from a netCDF file. The former is an arbitary arrangement based   */
/* on the grid generator's output, the latter is based on how the    */
/* unstructured indexing is arranged (wet cells surrounded by wet    */
/* cells listed first, then wet cells adjacent to land, then OBC     */
/* cells, then ghost cells). params->bathy resquires the former      */
/* arrangement, hence c2cc is used to convert from netCDF lists to   */
/* grid generator list.                                              */
/* Note: bathymetry is read in create_mesh() using the -p option.    */
/*-------------------------------------------------------------------*/
double *bathy_read_us(parameters_t *params, /* Input parameter data structure */
		      int cdfid, int ti)
{
  double *bathy, *inbath;          /* Bathymetry array */
  int c, cc;
  size_t ns2;
  size_t ndumps;
  size_t start[4];
  size_t count[4];
  int *c2cc;

  /* Check that the requested dump exists */
  nc_inq_dimlen(cdfid, ncw_dim_id(cdfid, "record"), &ndumps);
  if (ti > (int)(ndumps - 1)) {
    warn("dump_read: dump %d does not exist in file!\n", ti);
    return (0);
  }

  /* Allocate memory */
  bathy = d_alloc_1d(params->ns2+1);
  inbath = d_alloc_1d(params->ns2+1);
  c2cc = i_alloc_1d(params->ns2+1);

  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;
  /* time independent variables */
  count[1] = 0;
  count[2] = 0;
  count[3] = 0;

  if (DEBUG("init_m"))
    dlog("init_m", "\nReading unstructured bathymetry from input file\n\n");

  /* Read the bathymetry */
  count[0] = params->ns2+1;
  nc_get_vara_int(cdfid, ncw_var_id(cdfid, "Mesh2_index"), start, count, c2cc);
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "Mesh2_depth"), start, count, inbath);

  for (c = 1; c <= params->ns2; c++) {
    cc = c2cc[c];
    bathy[cc] = inbath[c];
    if(isnan(bathy[cc]))
      bathy[cc] = NOTVALID;
  }

  /* Set the grid topology from the input file data */
  /*netCDF_grid_us(params, cdfid, c2cc); Read in create_mesh()       */

  d_free_1d(inbath);
  i_free_1d(c2cc);
  return (bathy);
}

/* END bathy_read_us()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Read an unstructured double array from id                         */
/*-------------------------------------------------------------------*/
int dumpdata_read_2d_us(dump_data_t *dumpdata, int id, char *name,
			double *p, int dump, int ni, int oset)
{
  size_t start[4];
  size_t count[4];
  int c, cc;
  int vid;
  int istart = (oset) ? 0 : 1;

  start[0] = dump;
  start[1] = 0L;
  start[2] = 0L;
  start[3] = 0L;
  count[0] = 1L;
  count[1] = ni;
  count[2] = 0;
  count[3] = 0;
  
  vid = ncw_var_id(id, name);
  if (vid < 0) {
    memset(p, 0, ni * sizeof(double));
    return(1);
  }
  nc_get_vara_double(id, vid, start, count, dumpdata->w1s);
  for (cc = istart; cc < ni; cc++) {
    c = cc + params->oset;
    p[c] = dumpdata->w1s[cc];
  }
  return(0);
}

/* END dumpdata_read_2d_us()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Read an unstructured double array from id                         */
/*-------------------------------------------------------------------*/
int dumpdata_read_3d_us(dump_data_t *dumpdata, int id, char *name,
			double *p, double **d, int dump, int ni, int ni3, int *map, int **rmap, int oset)
{
  size_t start[4];
  size_t count[4];
  int nz = dumpdata->nz;
  int vid;
  int k, i, ii, c;
  double **in = d_alloc_2d(ni, dumpdata->nz);
  int istart = (oset) ? 0 : 1;

  start[0] = dump;
  start[0] = 0L;
  start[1] = 0L;
  start[2] = 0L;
  start[3] = 0L;
  count[0] = 1L;
  count[1] = dumpdata->nz;
  count[2] = ni;
  count[3] = 0;
  
  vid = ncw_var_id(id, name);
  if (vid < 0) {
    memset((void *)p, 0, sizeof(double) * ni3);
    return(1);
  }

  nc_get_vara_double(id, vid, start, count, in[0]);

  for (k = 0; k < nz; k++) {
    for (ii = istart; ii < ni; ii++) {
      i = map[ii];
      c = rmap[k][ii];
      if (c) p[rmap[k][ii]] = in[k][ii];
    }
  }
  d_free_2d(in);
  return(0);
}

/* END dumpdata_read_3d_us()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Read an unstructured integer array from id                        */
/*-------------------------------------------------------------------*/
int dumpdata_read_2d0_us(dump_data_t *dumpdata, int id, char *name,
			 double *p, int ni, int oset)
{
  size_t start[4];
  size_t count[4];
  int c, cc;
  int vid;
  int istart = (oset) ? 0 : 1;

  start[0] = 0L;
  start[1] = 0L;
  start[2] = 0L;
  start[3] = 0L;
  count[0] = ni;
  count[1] = 0L;
  count[2] = 0;
  count[3] = 0;
  
  vid = ncw_var_id(id, name);
  if (vid < 0) {
    memset(p, 0, ni * sizeof(double));
    return(1);
  }
  nc_get_vara_double(id, vid, start, count, dumpdata->w1s);
  for (cc = istart; cc < ni; cc++) {
    c = cc + params->oset;
    p[c] = dumpdata->w1s[cc];
  }
  return(0);
}

/* END dumpdata_read_2d0_us()                                        */
/*-------------------------------------------------------------------*/

#define DEG2RAD(d) ((d)*M_PI/180.0)
#define RADIUS 6370997.0
#define ECC 0.0
#define GEODESIC(x1, y1, x2, y2) (geod_inv_geod_fwd_sodanos(DEG2RAD(x1), DEG2RAD(y1),\
                                 DEG2RAD(x2), DEG2RAD(y2),\
                                 RADIUS, ECC))

/*-------------------------------------------------------------------*/
/* Modifies the bathymetry for unstructured meshes                   */
/*-------------------------------------------------------------------*/
void set_bathy_us(parameters_t *params)
{
  int cc, c;
  int i, j, k, n, m;            /* Counters                          */
  int is, ie, ip, im;           /* Limits of grid                    */
  int bathylimit = 0;           /* Flag to set min and max bathy     */
  int percent = 0;              /* Flag to set maximum elevation     */
  double min_cell_thickness = 0.0;  /* Minimum cell thickness        */
  double bmin, bmax;            /* Minimum and maximum bathymetry    */
  double etamax;                /* Maximum surface elevation         */
  double *mincelldz;            /* Minimum layer thickness           */
  double *layers;               /* Layer spacing array               */
  double maxgrad;               /* Maximum allowable bottom gradient */
  double mdh;                   /* Max bottom gradient found in file
                                   topo */
  double d1, d2;                /* Dummies                           */
  double dist = 0.0;            /* Grid spacing                      */
  double x, y;                  /* Coordinate distances              */
  int ic = 0;                   /* Grid position at maximum bottom   */
                                /* gradient.                         */
  int nps = 0;                  /* Smoothing passes to reduce bottom */
                                /* gradient.                         */
  double *bathy;                /* New bathymetry array              */
  int ns2 = params->ns2;        /* Size of bathy array               */
  int nz = params->nz;          /* Vertical layer numbers            */
  int maskf = 0, sf;            /* Mask                              */
  mesh_t *mesh;                 /* Mesh structure                    */

  if (params->us_type & US_RS) return;

  /*-----------------------------------------------------------------*/
  /* Allocate memory                                                 */
  mesh = params->mesh;
  bathy = d_alloc_1d(ns2+1);
  layers = params->layers;
  maxgrad = params->maxgrad;

  /*-----------------------------------------------------------------*/
  /* Read parameters required for adjusting the bathymetry           */
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
  /* Read the bathymetry array                                       */
  bmin = bmax = 0.0;
  if (params->bmin)
    bmin = params->bmin;
  if (params->bmax)
    bmax = params->bmax;
  if (bmin && bmax) {
    if (bmin > 0 && bmax > 0 && bmin > bmax)
      hd_quit("ERROR: set_bathy: BATHYMIN (%4.2f) > BATHYMAX (%4.2f)\n", bmin, bmax);
    else
      bathylimit = 1;
  }
  memcpy(bathy, params->bathy, ns2 * sizeof(double));


  /*-----------------------------------------------------------------*/
  /* Adjust for limits                                               */
  for (cc = 1; cc <= mesh->ns2; cc++) {

    /* Adjust botz to min and max bathymetry, checking for outside   */
    /* cells (botz < gridz[0]) and land cells (botz > etamax).       */
    if (bathylimit && bathy[cc] < etamax && bathy[cc] > -bmin)
      bathy[cc] = -bmin;
    if (bathylimit && bathy[cc] >= layers[0] && bathy[cc] < -bmax)
      bathy[cc] = -bmax;

    /* Load the cell values */
    for (k = 0; k < nz; k++) {
      /* Don't allow cells to be less than MIN_CELL_THICKNESS        */
      if (bathy[cc] < layers[k + 1] &&
	  bathy[cc] > layers[k + 1] - mincelldz[k]) {
	double newbotz;
	if ((layers[k + 1] - bathy[cc]) >
	    (layers[k + 1] - layers[k]) / 2)
	  newbotz = layers[k];
	else
	  newbotz = layers[k + 1];
	bathy[cc] = newbotz;
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Average bathymetry from surrounding cells                       */
  if (params->bathyfill & B_AVERAGE) {
    for (cc = 1; cc <= mesh->ns2; cc++) {
      if (bathy[cc] == LANDCELL || bathy[cc] == NOTVALID) {
	double bath = 0.0;
	double nb= 0.0;
	for (j = 1; j <= mesh->npe[cc]; j++) {
	  if ((c = mesh->neic[j][cc])) {
	    if (bathy[c] != LANDCELL && bathy[c] != NOTVALID) {
	      bath += bathy[c];
	      nb += 1.0;
	    }
	  }
	}
	bathy[cc] = (nb) ? bath / nb : -params->bmin;
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Clean up isolated cells                                         */
  for (cc = 1; cc <= mesh->ns2; cc++) {
    int nwc = 0;
    int twc = 0;
    for (j = 1; j <= mesh->npe[cc]; j++) {
      if ((c = mesh->neic[j][cc])) {
	twc++;
	if (bathy[c] != LANDCELL && bathy[c] != NOTVALID)
	  nwc++;
      }
    }
    if (nwc != twc) {
      bathy[cc] = LANDCELL;
      hd_warn("Filling in isolated wet cell at (%d %d)\n",i, j);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set any bathymetry masking                                      */
  d1 = bmin;
  is = ie = -1;
  prm_set_errfn(hd_silent_warn);
  if (prm_read_int(params->prmfd, "BATHY_MASK_IS", &is)) {
    prm_set_errfn(hd_quit);
    prm_read_int(params->prmfd, "BATHY_MASK_IE", &ie);
    prm_read_double(params->prmfd, "BATHY_MASK_DS", &d1);
    prm_read_double(params->prmfd, "BATHY_MASK_DE", &d2);
    d1 *= -1.0; d2 *= -1.0;
    prm_set_errfn(hd_silent_warn);
    maskf = 1;
  }

  /* Single block - gradient or single value */
  mdh = 0.0;
  if (prm_read_double(params->prmfd, "BATHY_MASK_VAL", &mdh)) {
    maskf = 3;
    if (mdh > 0.0) mdh *= -1.0;
  }
  if (maskf) {
    int *bmi, *bmj;
    int mode = read_blocks(params->prmfd, "BATHY_MASK", &n, &bmi, &bmj);

    if (mode & (B_LISTC|B_BLOCKC)) {
      for (k = 1; k <= n; k++) {
	cc = bmi[k];
	if (maskf == 1) {
	  bathy[cc] = (d2 - d1) * (double)(i - is) / (double)(ie - is) + d1;
	}
	else 
	  bathy[cc] = mdh;
      }
    }
  }

  /* Multiple blocks - single value only */
  if (prm_read_int(params->prmfd, "BATHY_NBLOCKS", &im)) {
    char buf[MAXSTRLEN];
    int *bmi, *bmj;
    int mode;
    for (ip = 0; ip < im; ip++) {
      mdh = 0.0;
      sprintf(buf, "BATHY_MASK_VAL%d", ip);
      if (prm_read_double(params->prmfd, buf, &mdh))
	if (mdh > 0.0) mdh *= -1.0;
      sprintf(buf, "BATHY_MASK%d", ip);
      mode = read_blocks(params->prmfd, buf, &n, &bmi, &bmj);
      if (mode & (B_LISTC|B_BLOCKC)) {
	for (k = 1; k <= n; k++) {
	  cc = bmi[k];
	  bathy[cc] = mdh;
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Smooth the bathymetry with a convolution filter if required.    */
  /* Smoothing over a selected area.                                 */
  if(params->smooth) {
    if (prm_read_int(params->prmfd, "SMOOTH_MASK", &m)) {
      double *nbthy = d_alloc_1d(ns2+1);
      double kk;
      n = params->smooth; 
      while(n) {
	prm_flush_line(params->prmfd);
	prm_read_int(params->prmfd, "SMOOTH_MASK", &m);
	for (k = 0; k < m; ++k) {
	  if (fscanf(params->prmfd, "%d", &cc) != 1)
	    hd_quit("set_bathy: Can't read cc in SMOOTH_MASK points list.\n");
	  prm_flush_line(params->prmfd);
	  if (bathy[cc] != LANDCELL && bathy[cc] != NOTVALID) {
	    nbthy[cc] = bathy_smooth_us(mesh, bathy, cc);
	    nbthy[cc] = min(nbthy[cc], -bmin);
	    nbthy[cc] = max(nbthy[cc], -bmax);
	  }
	}

	prm_flush_line(params->prmfd);
	prm_read_int(params->prmfd, "SMOOTH_MASK", &m);
	for (k = 0; k < m; ++k) {
	  if (fscanf(params->prmfd, "%d", &cc) != 1)
	    hd_quit("set_bathy: Can't read cc in SMOOTH_MASK points list.\n");
	  prm_flush_line(params->prmfd);
	  bathy[cc] = nbthy[cc];
	}
	n--;
      }
      d_free_1d(nbthy);
    } else {
      /* Smoothing over a the whole domain.                          */
      double *nbthy = d_alloc_1d(ns2+1);
      double kk;
      n = params->smooth;

      while(n) {
	for (cc = 1; cc <= ns2; cc++) {
	  if (bathy[cc] != LANDCELL && bathy[cc] != NOTVALID) {
	    nbthy[cc] = bathy_smooth_us(mesh, bathy, cc);
	    nbthy[cc] = min(nbthy[cc], -bmin);
	    nbthy[cc] = max(nbthy[cc], -bmax);
	  }
	}
	for (cc = 1; cc <= ns2; cc++) {
	  if (bathy[cc] != LANDCELL && bathy[cc] != NOTVALID) {
	    bathy[cc] = nbthy[cc];
	  }
	}
	n--;
      }
      d_free_1d(nbthy);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Check that the bathymetry gradient is OK                        */
  sf = 1;
  while (maxgrad > 0 && sf) {
    int jj, cn;
    double d;
    int is_geog = (strncasecmp(params->projection,
			       GEOGRAPHIC_TAG, strlen(GEOGRAPHIC_TAG)) == 0) ? 1 : 0;
    mdh = 0.0;
    sf = 0;
    for (cc = 1; cc <= ns2; cc++) {
      if (bathy[cc] != LANDCELL && bathy[cc] != NOTVALID) {
	d1 = d2 = 0.0; 
	for (j = 1; j <= mesh->npe[cc]; j++) {
	  if ((c = mesh->neic[j][cc])) {
	    if (bathy[c] != LANDCELL && bathy[c] != NOTVALID) {
	      if (is_geog) {
		d = GEODESIC(mesh->xloc[mesh->eloc[0][cc][0]], mesh->yloc[mesh->eloc[0][cc][0]],
			     mesh->xloc[mesh->eloc[0][c][0]], mesh->yloc[mesh->eloc[0][c][0]]);
	      } else {
		x = mesh->xloc[mesh->eloc[0][cc][0]] - mesh->xloc[mesh->eloc[0][c][0]];
		y = mesh->yloc[mesh->eloc[0][cc][0]] - mesh->yloc[mesh->eloc[0][c][0]];
		d = sqrt(x * x + y * y);
	      }
	      if (d) {
		d1 += fabs(bathy[c] - bathy[cc]) / d;
		d2 += 1.0;
	      }
	    }
	  }
	}
	if (d2) d1 /= d2;

	if (d1 > maxgrad) {
	  double kk = 0.0;
	  double nbthy = 0.0; 
	  for (j = 1; j <= mesh->npe[cc]; j++) {
	    double kkn = 0.0;
	    double nbthyn = 0.0; 
	    if ((c = mesh->neic[j][cc])) {
	      if (bathy[c] != LANDCELL && bathy[c] != NOTVALID) {
		kk += 1.0;
		nbthy += bathy[c];

		/* Average neighbours of c also                      */
		for (jj = 1; jj <= mesh->npe[c]; jj++) {
		  if ((cn = mesh->neic[j][c])) {
		    if (bathy[cn] != LANDCELL && bathy[cn] != NOTVALID) {
		      kkn += 1.0;
		      nbthyn += bathy[cn];
		    }
		  }
		}
		if (kkn) nbthyn /= kkn;
		if (params->runmode & (AUTO | DUMP)) {
		  bathy[cn] = min(nbthyn, -bmin);
		  bathy[cn] = max(nbthyn, -bmax);
		}
	      }
	    }
	  }
	  if (kk) nbthy /= kk;
	  if (params->runmode & (AUTO | DUMP)) {
	    bathy[cc] = min(nbthy, -bmin);
	    bathy[cc] = max(nbthy, -bmax);
	  }
	  sf = 1;
	  if (d1 > mdh) {
	    mdh = d1;
	    ic = cc;
	    dist = d;
	  }
	}
      }
    }
    nps++;
    if (mdh != 0.0) {
      printf
        ("Warning : Pass %d maximum bottom topography gradient too large : %4.3f\n",
         nps, mdh);
      printf("%6.2fm in %8.2fm at point (%d)\n", mdh * dist, dist, ic);
    }
    if (nps > 10)
      break;
  }

  /*-----------------------------------------------------------------*/
  /* Set uniform bathymetry within bathycon cells from the boundary  */
  for (n = 0; n < params->nobc; n++) {
    open_bdrys_t *open = params->open[n];
    if (open->bathycon) {
      int cn, jj[mesh->mnpe];
      double bathyb;
      /* Loop into the last non-ghost cell bathycon cells into the   */
      /* interior and get the bathymetry at this cell.               */
      for (n = 0; n < mesh->nobc; n++) {
	for (i = 1; i <= mesh->npts[n]; i++) {
	  cc = mesh->loc[n][i];
	  /* Get the direction (median) of the interior              */
	  m = 0;
	  for (j = 1; j <= mesh->npe[cc]; j++) {
	    if ((c = mesh->neic[j][cc])) {
	      if (bathy[c] != LANDCELL && bathy[c] != NOTVALID) {
		jj[m++] = j;
	      }
	    }
	  }
	  j = jj[m/2];
	  /* Map into the interior in this direction bathycon cells  */
	  /* and get the bathymetry at the final cell.               */
	  c = cc;
	  for (k = 0; k < open->bathycon; k++)
	    c = mesh->neic[j][c];
	  bathyb = bathy[c];
	  /* Set bathycon cells in this direction from the boundary  */
	  /* to bathyb.                                              */
	  c = cc;
	  for (k = 0; k < open->bathycon; k++) {
	    bathy[c] = bathyb;
	    c = mesh->neic[j][c];
	  }
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Smooth bathymetry within smooth_z cells from the boundary       */
  for (n = 0; n < params->nobc; n++) {
    open_bdrys_t *open = params->open[n];
    if (open->smooth_z && open->smooth_n) {
      double *nbthy = d_alloc_1d(ns2+1);
      int np;
      int cn, jj[mesh->mnpe];
      double bathyb;
      for (np = 0; np < open->smooth_n; np++) {
	for (n = 0; n < mesh->nobc; n++) {
	  for (i = 1; i <= mesh->npts[n]; i++) {
	    cc = mesh->loc[n][i];
	    /* Get the direction (median) of the interior              */
	    m = 0;
	    for (j = 1; j <= mesh->npe[cc]; j++) {
	      if ((c = mesh->neic[j][cc])) {
		if (bathy[c] != LANDCELL && bathy[c] != NOTVALID) {
		  jj[m++] = j;
		}
	      }
	    }
	    j = jj[m/2];
	    /* Map into the interior in this direction smooth_z cells  */
	    /* and smooth the bathymetry.                              */
	    c = cc;
	    for (k = 0; k < open->smooth_z; k++) {
	      if (c == 0) break;
	      nbthy[c] = bathy_smooth_us(mesh, bathy, c);
	      c = mesh->neic[j][c];
	    }
	    c = cc;
	    for (k = 0; k < open->smooth_z; k++) {
	      bathy[c] = nbthy[c];
	      c = mesh->neic[j][c];
	    }
	  }
	}
      }
      d_free_1d(nbthy);
    }
  }
  memcpy(params->bathy, bathy, ns2 * sizeof(double));
  d_free_1d(bathy);
}

/* END set_bathy_us()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Modifies the bathymetry for structured meshes                     */
/*-------------------------------------------------------------------*/
double **set_bathy(parameters_t *params,  /* Input parameter data
                                             structure */
                   int cdfid    /* Input netCDF file handle */
  )
{
  int i, j, k, n;               /* Counters */
  int is, ie, js, je;           /* Limits of grid */
  int xsize, ysize;             /* Size of original grid */
  int bathylimit = 0;           /* Flag to set min and max bathymetry */
  int percent = 0;              /* Flag to set maximum elevation */
  double min_cell_thickness = 0.0;  /* Minimum cell thickness */
  unsigned long topflag = 0;    /* Surface flag value */
  double bmin, bmax;            /* Minimum and maximum bathymetry */
  double etamax;                /* Maximum surface elevation */
  double *mincelldz;            /* Minimum layer thickness */
  unsigned long ***flag;        /* Flag for Cartesian grid */
  double *layers;               /* Layer spacing array */
  double maxgrad;               /* Maximum allowable bottom gradient */
  double mdh;                   /* Maximum bottom gradient found in file
                                   topo */
  double d1, d2;                /* Dummies */
  double dist = 0.0;            /* Grid spacing over maximum bottom
                                   gradient */
  int ic = 0, jc = 0;           /* Grid position at maximum bottom
                                   gradient */
  int nps = 0;                  /* Smoothing passes to reduce bottom
                                   gradient */
  int **xp, **yp, **xm, **ym;
  int ip, jp, im, jm;
  int ipp, jpp, imm, jmm;
  double **h1au1, **h1au2;
  double **h2au2, **h2au1;
  double **bathy;               /* Bathymetry array */
  int nce1 = params->nce1;
  int nce2 = params->nce2;
  int nz = params->nz;
  int sf = 1;
  int maskf = 0;
  char dir = '\0';

  if (params->us_type & US_RUS) return(NULL);

  /* Allocate memory */
  bathy = d_alloc_2d(nce1, nce2);
  layers = params->layers;
  flag = (unsigned long ***)l_alloc_3d(nce1 + 1, nce2 + 1, nz);
  h1au1 = d_alloc_2d(nce1 + 1, nce2);
  h2au1 = d_alloc_2d(nce1 + 1, nce2);
  h2au2 = d_alloc_2d(nce1, nce2 + 1);
  h1au2 = d_alloc_2d(nce1, nce2 + 1);
  maxgrad = params->maxgrad;

  /*-----------------------------------------------------------------*/
  /* Read parameters required for adjusting the bathymetry */
  is = js = 0;
  ie = xsize = nce1;
  je = ysize = nce2;

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
  /* Read the bathymetry array */
  bmin = bmax = 0.0;
  if (params->bmin)
    bmin = params->bmin;
  if (params->bmax)
    bmax = params->bmax;
  if (bmin && bmax) {
    if (bmin > 0 && bmax > 0 && bmin > bmax)
      hd_quit("ERROR: set_bathy: BATHYMIN (%4.2f) > BATHYMAX (%4.2f)\n", bmin, bmax);
    else
      bathylimit = 1;
  }

  for (j = 0; j < nce2; j++)
    for (i = 0; i < nce1; i++)
      bathy[j][i] = LANDCELL;

  if (params->nvals == 0) {
    /* Use default bathymetry */
    for (j = js; j < je; j++)
      for (i = is; i < ie; i++)
	bathy[j][i] = layers[0];
  } else if (params->nvals == (xsize) * (ysize)) {
    /* Use bathymetry as read */
    n = 0;
    for (j = js; j < je; j++)
      for (i = is; i < ie; i++)
	bathy[j][i] = -params->bathy[n++];
  } else
    hd_quit
      ("set_bathy: bad bathymetry data (neither 0 nor number of cells)\n");

  /*-----------------------------------------------------------------*/
  /* Set boundary OUTSIDE OBC cells if required                      */
  adjust_outside_obc(params, bathy);

  /*-----------------------------------------------------------------*/
  /* Initialise the flag array */
  for (j = 0; j <= nce2; j++)
    for (i = 0; i <= nce1; i++)
      for (k = 0; k < nz; k++) {
        flag[k][j][i] = 0;
      }

  /*-----------------------------------------------------------------*/
  /* Set the flags to solid on the nce1 and nce2 edges */
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

  /* Get the horizontal grid spacings */
  if (params->runmode & (AUTO | DUMP)) {
    for (j = 0; j < params->nce2 + 1; j++) {
      for (i = 0; i < params->nce1; i++) {
        h2au2[j][i] = 2.0 * params->h2[j * 2][i * 2 + 1];
        h1au2[j][i] = 2.0 * params->h1[j * 2][i * 2 + 1];
      }
    }
    for (j = 0; j < params->nce2; j++) {
      for (i = 0; i < params->nce1 + 1; i++) {
        h1au1[j][i] = 2.0 * params->h1[j * 2 + 1][i * 2];
        h2au1[j][i] = 2.0 * params->h2[j * 2 + 1][i * 2];
      }
    }
  } else {
    size_t start[4];
    size_t count[4];

    start[0] = 0;
    start[1] = 0;
    start[2] = 0;
    start[3] = 0;
    count[0] = nce2;
    count[1] = nce1 + 1;
    count[2] = 0;
    count[3] = 0;
    nc_get_vara_double(cdfid, ncw_var_id(cdfid, "h1au1"), start, count,
                       h1au1[0]);
    count[0] = nce2 + 1;
    count[1] = nce1;
    nc_get_vara_double(cdfid, ncw_var_id(cdfid, "h2au2"), start, count,
                       h2au2[0]);
  }

  /* Clean up isolated cells */
    for (j = 0; j < nce2; j++) {
      for (i = 0; i < nce1; i++) {
	if (!(flag[nz - 1][j][i] & (SOLID | OUTSIDE))) {
	  ip = (i < nce1 - 1) ? i+1 : i;
	  im = (i > 0) ? i-1 : i;
	  jp = (j < nce2 - 1) ? j+1 : j;
	  jm = (j > 0) ? j-1 : j;
	  if (flag[nz - 1][j][im] & (SOLID | OUTSIDE) &&
	      flag[nz - 1][j][ip] & (SOLID | OUTSIDE) &&
	      flag[nz - 1][jm][i] & (SOLID | OUTSIDE) &&
	      flag[nz - 1][jp][i] & (SOLID | OUTSIDE)) {
	    bathy[j][i] = LANDCELL;
	    for (k = 0; k < nz; k++) {
	      flag[k][j][i] |= SOLID;
	    }
	    hd_warn("Filling in isolated wet cell at (%d %d)\n",i, j);
	  }
	}
      }
    }

    /*remove_channel(params, flag, bathy);*/

  /*-----------------------------------------------------------------*/
  /* Set any bathymetry masking                                      */
  d1 = bmin;
  is = ie = js = je = -1;
  prm_set_errfn(hd_silent_warn);
  if (prm_read_int(params->prmfd, "BATHY_MASK_IS", &is)) {
    prm_set_errfn(hd_quit);
    prm_read_int(params->prmfd, "BATHY_MASK_IE", &ie);
    prm_read_double(params->prmfd, "BATHY_MASK_DS", &d1);
    prm_read_double(params->prmfd, "BATHY_MASK_DE", &d2);
    d1 *= -1.0; d2 *= -1.0;
    prm_set_errfn(hd_silent_warn);
    maskf = 1;
  } else if (prm_read_int(params->prmfd, "BATHY_MASK_JS", &js)) {
    prm_set_errfn(hd_quit);
    prm_read_int(params->prmfd, "BATHY_MASK_JE", &je);
    prm_read_double(params->prmfd, "BATHY_MASK_DS", &d1);
    prm_read_double(params->prmfd, "BATHY_MASK_DE", &d2);
    d1 *= -1.0; d2 *= -1.0;
    prm_set_errfn(hd_silent_warn);
    maskf = 2;
  }

  /* Single block - gradient or single value */
  mdh = 0.0;
  if (prm_read_double(params->prmfd, "BATHY_MASK_VAL", &mdh)) {
    maskf = 3;
    mdh *= -1.0;
  }
  if (maskf) {
    int *bmi, *bmj;
    int mode = read_blocks(params->prmfd, "BATHY_MASK", &n, &bmi, &bmj);
    if (mode & (B_LISTIJ|B_BLOCKIJ)) {
      for (k = 1; k <= n; k++) {
	i = bmi[k];
	j = bmj[k];
	if (maskf == 1) {
	  bathy[j][i] = (d2 - d1) * (double)(i - is) / (double)(ie - is) + d1;
	}
	else if (maskf == 2)
	  bathy[j][i] = (d2 - d1) * (double)(j - js) / (double)(je - js) + d1;
	else 
	  bathy[j][i] = mdh;
      }
    }
  }

  /* Multiple blocks - single value only */
  if (prm_read_int(params->prmfd, "BATHY_NBLOCKS", &im)) {
    char buf[MAXSTRLEN];
    int *bmi, *bmj;
    for (ip = 0; ip < im; ip++) {
      mdh = 0.0;
      sprintf(buf, "BATHY_MASK_VAL%d", ip);
      if (prm_read_double(params->prmfd, buf, &mdh))
	mdh *= -1.0;
      sprintf(buf, "BATHY_MASK%d", ip);
      read_blocks(params->prmfd, buf, &n, &bmi, &bmj);
      for (k = 1; k <= n; k++) {
	i = bmi[k];
	j = bmj[k];
	bathy[j][i] = mdh;
      }
    }
  }

  /*
  if (prm_read_int(params->prmfd, "BATHY_MASK", &n)) {
    int i1, i2, j1, j2;
    prm_flush_line(params->prmfd);
    if (fscanf(params->prmfd, "(%d,%d)-(%d,%d)", &i1, &i2, &j1, &j2) == 4) {
      prm_flush_line(params->prmfd);
      for (j = 0; j < params->nce2; j++) {
	for (i = 0; i < params->nce1; i++) {
	  if (i >= i1 && i <= i2 && j >= j1 && j <= j2) {
	    bathy[j][i] = mdh;
	    if (i >= is && i <= ie)
	      bathy[j][i] = (d2 - d1) * (double)(i - is) / (double)(ie - is) + d1;
	    if (j >= js && j <= je)
	      bathy[j][i] = (d2 - d1) * (double)(j - js) / (double)(je - js) + d1;
	  }
	}
      }
    }
  }

  if (prm_read_int(params->prmfd, "BATHY_MASK", &n)) {
    prm_flush_line(params->prmfd);
    for (k = 0; k < n; ++k) {
      if (fscanf(params->prmfd, "%d %d", &i, &j) != 2)
	hd_quit("set_bathy: Can't read i j in BATHY_MASK points list.\n");
      prm_flush_line(params->prmfd);
      bathy[j][i] = mdh;
      if (i >= is && i <= ie)
	bathy[j][i] = (d2 - d1) * (double)(i - is) / (double)(ie - is) + d1;
      if (j >= js && j <= je)
	bathy[j][i] = (d2 - d1) * (double)(j - js) / (double)(je - js) + d1;
    }
  }
  */

  /* Set up mapping arrays */
  xp = i_alloc_2d(nce1, nce2);
  xm = i_alloc_2d(nce1, nce2);
  for (j = 0; j < nce2; j++) {
    for (i = 0; i < nce1; i++) {
      if (i < nce1 - 1)
	xp[j][i] =
	  ((flag[nz - 1][j][i + 1] & (SOLID | OUTSIDE) ||
	    i == nce1 - 1) ? i : i + 1);
      else
	xp[j][i] = i;
      if (i > 0)
	xm[j][i] =
	  ((flag[nz - 1][j][i - 1] & (SOLID | OUTSIDE) ||
	    i == 0) ? i : i - 1);
      else
	xm[j][i] = i;
    }
  }
  for (i = 0; i < params->nemx; i++) {
    xp[params->emjsx[i]][params->emisx[i]] = params->emidx[i];
    xm[params->emjdx[i]][params->emidx[i]] = params->emisx[i];
  }

  yp = i_alloc_2d(nce1, nce2);
  ym = i_alloc_2d(nce1, nce2);
  for (i = 0; i < nce1; i++) {
    for (j = 0; j < nce2; j++) {
      if (j < nce2 - 1)
	yp[j][i] =
	  ((flag[nz - 1][j + 1][i] & (SOLID | OUTSIDE) ||
	    j == nce2 - 1) ? j : j + 1);
      else
	yp[j][i] = j;
      if (j > 0)
	ym[j][i] =
	  ((flag[nz - 1][j - 1][i] & (SOLID | OUTSIDE) ||
	    j == 0) ? j : j - 1);
      else
	ym[j][i] = j;
    }
  }
  for (j = 0; j < params->nemy; j++) {
    yp[params->emjsy[j]][params->emisy[j]] = params->emidy[j];
    ym[params->emjdy[j]][params->emidy[j]] = params->emisy[j];
  }

  /*-----------------------------------------------------------------*/
  /* Smooth the bathymetry with a convolution filter if required.    */
  /* Smoothing over a selected area.                                 */
  if(params->smooth) {
    int m;
    if (prm_read_int(params->prmfd, "SMOOTH_MASK", &m)) {
      double **nbthy = d_alloc_2d(nce1, nce2);
      double kk[4][4];
      /* Set the filter weights */
      kk[1][3] = 0.0625;
      kk[2][3] = 0.125;
      kk[3][3] = 0.0625;
      kk[1][2] = 0.125;
      kk[2][2] = 0.25;
      kk[3][2] = 0.125;
      kk[1][1] = 0.0625;
      kk[2][1] = 0.125;
      kk[3][1] = 0.0625;
      n = params->smooth; 
      while(n) {
	prm_flush_line(params->prmfd);
	prm_read_int(params->prmfd, "SMOOTH_MASK", &m);
	for (k = 0; k < m; ++k) {
	  if (fscanf(params->prmfd, "%d %d", &i, &j) != 2)
	    hd_quit("set_bathy: Can't read i j in SMOOTH_MASK points list.\n");
	  prm_flush_line(params->prmfd);
	  if (!(flag[nz - 1][j][i] & (SOLID | OUTSIDE))) {
	    im = xm[j][i];
	    ip = xp[j][i];
	    jm = ym[j][i];
	    jp = yp[j][i];
	    /* Perform the convolution */
	    nbthy[j][i] = kk[1][3] * min(bathy[jp][im], 0.0) + 
	                  kk[2][3] * min(bathy[jp][i], 0.0) + 
	                  kk[3][3] * min(bathy[jp][ip], 0.0) +
	                  kk[1][2] * min(bathy[j][im], 0.0) + 
	                  kk[2][2] * min(bathy[j][i], 0.0) + 
	                  kk[3][2] * min(bathy[j][ip], 0.0) +
                          kk[1][1] * min(bathy[jm][im], 0.0) + 
	                  kk[2][1] * min(bathy[jm][i], 0.0) + 
	                  kk[3][1] * min(bathy[jm][ip], 0.0);
	    nbthy[j][i] = min(nbthy[j][i], -bmin);
	    nbthy[j][i] = max(nbthy[j][i], -bmax);
	  }
	}
	prm_flush_line(params->prmfd);
	prm_read_int(params->prmfd, "SMOOTH_MASK", &m);
	for (k = 0; k < m; ++k) {
	  if (fscanf(params->prmfd, "%d %d", &i, &j) != 2)
	    hd_quit("set_bathy: Can't read i j in SMOOTH_MASK points list.\n");
	  prm_flush_line(params->prmfd);
	  bathy[j][i] = nbthy[j][i];      
	}
	n--;
      }
      d_free_2d(nbthy);
    } else {
      /* Smoothing over a the whole domain.                          */
      double **nbthy = d_alloc_2d(nce1, nce2);
      double kk[4][4];
      n = params->smooth;

      /* Set the filter weights */
      kk[1][3] = 0.0625;
      kk[2][3] = 0.125;
      kk[3][3] = 0.0625;
      kk[1][2] = 0.125;
      kk[2][2] = 0.25;
      kk[3][2] = 0.125;
      kk[1][1] = 0.0625;
      kk[2][1] = 0.125;
      kk[3][1] = 0.0625;

      while(n) {
	for (j = 0; j < nce2; j++) {
	  for (i = 0; i < nce1; i++) {
	    if (!(flag[nz - 1][j][i] & (SOLID | OUTSIDE))) {
	      double k0 = kk[2][2];
	      double k1 = kk[1][2], k2 = kk[3][2];
	      double k3 = kk[1][1], k4 = kk[2][1], k5 = kk[3][1];
	      double k6 = kk[1][3], k7 = kk[2][3], k8 = kk[3][3];
	      if (flag[nz - 1][j][i-1] & (SOLID | OUTSIDE)) {
		k0 += k1;
		k1 = 0.0;
	      }
	      if (flag[nz - 1][j][i+1] & (SOLID | OUTSIDE)) {
		k0 += k2;
		k2 = 0.0;
	      }
	      if (flag[nz - 1][j-1][i-1] & (SOLID | OUTSIDE)) {
		k0 += k3;
		k3 = 0.0;
	      }
	      if (flag[nz - 1][j-1][i] & (SOLID | OUTSIDE)) {
		k0 += k4;
		k4 = 0.0;
	      }
	      if (flag[nz - 1][j-1][i+1] & (SOLID | OUTSIDE)) {
		k0 += k5;
		k5 = 0.0;
	      }
	      if (flag[nz - 1][j+1][i-1] & (SOLID | OUTSIDE)) {
		k0 += k6;
		k6 = 0.0;
	      }
	      if (flag[nz - 1][j+1][i] & (SOLID | OUTSIDE)) {
		k0 += k7;
		k7 = 0.0;
	      }
	      if (flag[nz - 1][j+1][i+1] & (SOLID | OUTSIDE)) {
		k0 += k8;
		k8 = 0.0;
	      }
	      im = xm[j][i];
	      ip = xp[j][i];
	      jm = ym[j][i];
	      jp = yp[j][i];
	      /* Perform the convolution */
	      nbthy[j][i] = k6 * bathy[jp][im] + 
	                    k7 * bathy[jp][i] + 
	                    k8 * bathy[jp][ip] +
	                    k1 * bathy[j][im] + 
	                    k0 * bathy[j][i] + 
	                    k2 * bathy[j][ip] +
                            k3 * bathy[jm][im] + 
	                    k4 * bathy[jm][i] + 
	                    k5 * bathy[jm][ip];
	      nbthy[j][i] = min(nbthy[j][i], -bmin);
	      nbthy[j][i] = max(nbthy[j][i], -bmax);
	    }
	  }
	}
	for (j = 0; j < nce2; j++) {
	  for (i = 0; i < nce1; i++) {
	    if (!(flag[nz - 1][j][i] & (SOLID | OUTSIDE))) {
	      bathy[j][i] = nbthy[j][i];
	    }
	  }
	}
	n--;
      }
      d_free_2d(nbthy);
    }
  }

  /* Check that the bathymetry gradient is OK */
  while (maxgrad > 0 && sf) {
    mdh = 0.0;
    sf = 0;
    for (j = 0; j < nce2; j++) {
      for (i = 0; i < nce1; i++) {
        if (!(flag[nz - 1][j][i] & (SOLID | OUTSIDE))) {
          im = xm[j][i];
          ip = xp[j][i];
          jm = ym[j][i];
          jp = yp[j][i];
          /* Backward x gradient */
          d1 = fabs(bathy[j][i] - bathy[j][im]) / h1au1[j][i];
          d2 = fabs(bathy[j][ip] - bathy[j][i]) / h1au1[j][ip];
          if (d1 > maxgrad || d2 > maxgrad) {
            ipp = xp[j][ip];
            imm = xm[j][im];
            if (params->runmode & (AUTO | DUMP)) {
              bathy[j][i] =
                smoothbathy(bathy, i, ip, im, ipp, imm, j, 1, bmin);
              if (!(flag[nz - 1][j][i + 1] & (SOLID | OUTSIDE)))
                bathy[j][ip] =
                  smoothbathy(bathy, ip, ipp, i, xp[j][ipp], im, j, 1,
                              bmin);
              if (!(flag[nz - 1][j][i - 1] & (SOLID | OUTSIDE)))
                bathy[j][im] =
                  smoothbathy(bathy, im, i, imm, ip, xm[j][imm], j, 1,
                              bmin);
            }
            if (d1 > mdh) {
              mdh = d1;
              ic = i;
              jc = j;
              dist = h1au1[j][i];
              dir = 'x';
            }
          }
          /* Backward y gradient */
          d2 = fabs(bathy[j][i] - bathy[jm][i]) / h2au2[j][i];
          d1 = fabs(bathy[jp][i] - bathy[j][i]) / h2au2[jp][i];
          if (d1 > maxgrad || d2 > maxgrad) {
            jpp = yp[jp][i];
            jmm = ym[jm][i];
            if (params->runmode & (AUTO | DUMP)) {
              bathy[j][i] =
                smoothbathy(bathy, j, jp, jm, jpp, jmm, i, 0, bmin);
              if (!(flag[nz - 1][j + 1][i] & (SOLID | OUTSIDE)))
                bathy[jp][i] =
                  smoothbathy(bathy, jp, jpp, j, yp[jpp][i], jm, i, 0,
                              bmin);
              if (!(flag[nz - 1][j - 1][i] & (SOLID | OUTSIDE)))
                bathy[jm][i] =
                  smoothbathy(bathy, jm, j, jmm, jp, ym[jmm][i], i, 0,
                              bmin);
            }
            sf = 1;
            if (d2 > mdh) {
              mdh = d2;
              ic = i;
              jc = j;
              dist = h2au2[j][i];
              dir = 'y';
            }
          }
        }
      }
    }
    nps++;
    if (mdh != 0.0) {
      printf
        ("Warning : Pass %d %c maximum bottom topography gradient too large : %4.3f\n",
         nps, dir, mdh);
      printf("%6.2fm in %8.2fm at point (%d,%d)\n", mdh * dist, dist, ic,
             jc);
    }
    if (nps > 10)
      break;
  }

  /* Check that uniform flow over topography will not violate the */
  /* vertical velocity Courant number.  */
  if (params->bvel > 0.0) {
    double ***dz, w, w1, w2, ***u1flux, ***u2flux, top, bot;
    double fcbotx, fctopx, fcboty, fctopy;
    short **kb;
    double vel = params->bvel;
    dz = d_alloc_3d(nce1, nce2, nz);
    u1flux = d_alloc_3d(nce1 + 1, nce2, nz);
    u2flux = d_alloc_3d(nce1, nce2 + 1, nz);
    kb = s_alloc_2d(nce1, nce2);
    sf = 1;
    nps = 0;

    while (sf) {
      sf = 0;
      /* Get the cell thickness */
      for (j = 0; j < nce2; j++) {
        for (i = 0; i < nce1; i++) {
          kb[j][i] = nz;
          if (!(flag[nz - 1][j][i] & (SOLID | OUTSIDE))) {
            kb[j][i] = 1;
            while (kb[j][i] < nz && layers[kb[j][i]] < bathy[j][i])
              kb[j][i]++;
            kb[j][i]--;
            for (k = kb[j][i]; k < nz; k++) {
              top = (k == nz - 1) ? 0 : layers[k + 1];
              bot = (k == kb[j][i]) ? bathy[j][i] : layers[k];
              dz[k][j][i] = top - bot;
            }
          }
        }
      }

      /* Get the fluxes */
      for (j = 0; j < nce2; j++) {
        for (i = 1; i < nce1; i++) {
          for (k = 0; k < nz; k++)
            u1flux[k][j][i] = 0.0;
          if (!(flag[nz - 1][j][i] & (SOLID | OUTSIDE))) {
            for (k = kb[j][i]; k < nz; k++) {
              if (!(flag[k][j][i] & (SOLID | OUTSIDE)) &&
                  !(flag[k][j][i - 1] & (SOLID | OUTSIDE))) {
                u1flux[k][j][i] = vel * dz[k][j][i] * h2au1[j][i];
              }
            }
          }
        }
      }
      for (j = 1; j < nce2; j++) {
        for (i = 0; i < nce1; i++) {
          for (k = 0; k < nz; k++)
            u2flux[k][j][i] = 0.0;
          if (!(flag[nz - 1][j][i] & (SOLID | OUTSIDE))) {
            for (k = kb[j][i]; k < nz; k++) {
              if (!(flag[k][j][i] & (SOLID | OUTSIDE)) &&
                  !(flag[k][j - 1][i] & (SOLID | OUTSIDE))) {
                u2flux[k][j][i] = vel * dz[k][j][i] * h1au2[j][i];
              }
            }
          }
        }
      }

      /* Get the vertical velocity */
      for (j = 0; j < nce2; j++) {
        for (i = 0; i < nce1; i++) {
          if (!(flag[nz - 1][j][i] & (SOLID | OUTSIDE))) {
            /* x velocity */
            w1 = w2 = fcbotx = fcboty = 0.0;
            for (k = kb[j][i]; k < nz - 1; k++) {

              fctopx = fcbotx + u1flux[k][j][i] - u1flux[k][j][i + 1];
              w1 = fabs(fctopx / (h2au1[j][i] * h1au2[j][i]));
              fcbotx = fctopx;
              fctopy = fcboty + u2flux[k][j][i] - u2flux[k][j + 1][i];
              w2 = fabs(fctopy / (h2au1[j][i] * h1au2[j][i]));
              fcboty = fctopy;
              w = max(w1, w2);
              d1 = w * params->grid_dt / dz[k][j][i];
              if (d1 > 1.0) {
                sf = 1;
                if (params->runmode == AUTO) {

                  if (!(flag[nz - 1][j][i + 1] & (SOLID | OUTSIDE)) &&
                      !(flag[nz - 1][j][i - 1] & (SOLID | OUTSIDE))) {
                    printf
                      ("Pass %d e1 Courant violation at (%d %d %d) %f %f %f\n",
                       nps, i, j, k, d1, w1, bathy[j][i]);
                    im = xm[j][i];
                    ip = xp[j][i];
                    ipp = xp[j][ip];
                    imm = xm[j][im];
                    bathy[j][i] =
                      smoothbathy(bathy, i, ip, im, ipp, imm, j, 1, bmin);
                    bathy[j][ip] =
                      smoothbathy(bathy, ip, ipp, i, xp[j][ipp], im, j, 1,
                                  bmin);
                    bathy[j][im] =
                      smoothbathy(bathy, im, i, imm, ip, xm[j][imm], j, 1,
                                  bmin);
                  }
                  if (!(flag[nz - 1][j + 1][i] & (SOLID | OUTSIDE)) &&
                      !(flag[nz - 1][j - 1][i] & (SOLID | OUTSIDE))) {
                    printf
                      ("Pass %d e2 Courant violation at (%d %d %d) %f %f %f\n",
                       nps, i, j, k, d1, w2, bathy[j][i]);
                    jm = ym[j][i];
                    jp = yp[j][i];
                    jpp = yp[jp][i];
                    jmm = ym[jm][i];
                    bathy[j][i] =
                      smoothbathy(bathy, j, jp, jm, jpp, jmm, i, 0, bmin);
                    bathy[jp][i] =
                      smoothbathy(bathy, jp, jpp, j, yp[jpp][i], jm, i, 0,
                                  bmin);
                    bathy[jm][i] =
                      smoothbathy(bathy, jm, j, jmm, jp, ym[jmm][i], i, 0,
                                  bmin);
                  }
                }
              }
            }
          }
        }
      }
      nps++;
      if (nps > 10)
        break;
    }

    d_free_3d(dz);
    d_free_3d(u1flux);
    d_free_3d(u2flux);
    s_free_2d(kb);
  }

  /* Set uniform bathymetry within bathycon cells from the boundary  */
  for (n = 0; n < params->nobc; n++) {
    open_bdrys_t *open = params->open[n];
    if (open->bathycon) {
      int m;
      /* Loop into the last non-ghost cell bathycon cells into the   */
      /* interior and get the bathymetry at this cell.               */
      for (m = 0; m < open->npts; m++) {	  
	i = open->iloc[m];
	j = open->jloc[m];
	if (open->type & U2BDRY) {
	  jc = j;
	  if (j == nce2 || (j != nce2 && j == yp[j][i] &&
			   (flag[nz - 1][j][i] & OUTSIDE))) {
	    jc--;
	    for (k = 0; k < open->bathycon; k++) {
	      jc = ym[jc][i];
	    }
	    d1 = bathy[jc][i];
	    jc = j - 1;
	    for (k = 0; k < open->bathycon; k++) {
	      bathy[jc][i] = d1;
	      jc = ym[jc][i];
	    }
	  } else {
	    for (k = 0; k < open->bathycon; k++)
	      jc = yp[jc][i];
	    d1 = bathy[jc][i];
	    jc = j;
	    for (k = 0; k < open->bathycon; k++) {
	      bathy[jc][i] = d1;
	      jc = yp[jc][i];
	    }
	  }
	}
	if (open->type & U1BDRY) {
	  ic = i;
	  if (i == nce1 || (i != nce1 && i == xp[j][i] &&
			   (flag[nz - 1][j][i] & OUTSIDE))) {
	    ic--;
	    for (k = 0; k < open->bathycon; k++)
	      ic = xm[j][ic];
	    d1 = bathy[j][ic];
	    ic = i - 1;
	    for (k = 0; k < open->bathycon; k++) {
	      bathy[j][ic] = d1;
	      ic = xm[j][ic];
	    }
	  } else {
	    for (k = 0; k < open->bathycon; k++)
	      ic = xp[j][ic];
	    d1 = bathy[j][ic];
	    ic = i;
	    for (k = 0; k < open->bathycon; k++) {
	      bathy[j][ic] = d1;
	      ic = xp[j][ic];
	    }
	  }
	}
      }
    }
  }

  /* Smooth bathymetry within smooth_z cells from the boundary  */
  for (n = 0; n < params->nobc; n++) {
    open_bdrys_t *open = params->open[n];
    if (open->smooth_z && open->smooth_n) {
      double **nbthy = d_alloc_2d(nce1, nce2);
      int m, np;
      /* Loop into the last non-ghost cell bathycon cells into the   */
      /* interior and get the bathymetry at this cell.               */
      for (np = 0; np < open->smooth_n; np++) {
	/* Smooth */
	for (m = 0; m < open->npts; m++) {	  
	  i = open->iloc[m];
	  j = open->jloc[m];
	  if (open->type & U2BDRY) {
	    jc = j;
	    if (j == nce2 || (j != nce2 && j == yp[j][i] &&
			      (flag[nz - 1][j][i] & OUTSIDE))) {
	      jc--;
	      for (k = 0; k < open->smooth_z; k++) {
		nbthy[jc][i] = smooth_bathy(flag, nz, i, jc, xm, xp, ym, yp, bathy, bmin, bmax);
		jc = ym[jc][i];
	      }
	    } else {
	      for (k = 0; k < open->smooth_z; k++) {
		nbthy[jc][i] = smooth_bathy(flag, nz, i, jc, xm, xp, ym, yp, bathy, bmin, bmax);
		jc = yp[jc][i];
	      }
	    }
	  }
	  if (open->type & U1BDRY) {
	    ic = i;
	    if (i == nce1 || (i != nce1 && i == xp[j][i] &&
			      (flag[nz - 1][j][i] & OUTSIDE))) {
	      ic--;
	      for (k = 0; k < open->smooth_z; k++) {
		nbthy[j][ic] = smooth_bathy(flag, nz, ic, j, xm, xp, ym, yp, bathy, bmin, bmax);
		ic = xm[j][ic];
	      }
	    } else {
	      for (k = 0; k < open->smooth_z; k++) {
		nbthy[j][ic] = smooth_bathy(flag, nz, ic, j, xm, xp, ym, yp, bathy, bmin, bmax);
		ic = xp[j][ic];
	      }
	    }
	  }
	}
	/* Copy */
	for (m = 0; m < open->npts; m++) {	  
	  i = open->iloc[m];
	  j = open->jloc[m];
	  if (open->type & U2BDRY) {
	    jc = j;
	    if (j == nce2 || (j != nce2 && j == yp[j][i] &&
			      (flag[nz - 1][j][i] & OUTSIDE))) {
	      jc--;
	      for (k = 0; k < open->smooth_z; k++) {
		bathy[jc][i] = nbthy[jc][i];
		jc = ym[jc][i];
	      }
	    } else {
	      for (k = 0; k < open->smooth_z; k++) {
		bathy[jc][i] = nbthy[jc][i];
		jc = yp[jc][i];
	      }
	    }
	  }
	  if (open->type & U1BDRY) {
	    ic = i;
	    if (i == nce1 || (i != nce1 && i == xp[j][i] &&
			      (flag[nz - 1][j][i] & OUTSIDE))) {
	      ic--;
	      for (k = 0; k < open->smooth_z; k++) {
		bathy[j][ic] = nbthy[j][ic];
		ic = xm[j][ic];
	      }
	    } else {
	      for (k = 0; k < open->smooth_z; k++) {
		bathy[j][ic] = nbthy[j][ic];
		ic = xp[j][ic];
	      }
	    }
	  }
	}
      }
      d_free_2d(nbthy);
    }
  }

  i_free_2d(xp);
  i_free_2d(xm);
  i_free_2d(yp);
  i_free_2d(ym);
  d_free_2d(h1au1);
  d_free_2d(h2au2);
  d_free_2d(h2au1);
  d_free_2d(h1au2);
  l_free_3d((long ***)flag);

  return (bathy);
}

/* END set_bathy()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Smooths the bathymetry at a point for unructured meshes           */
/*-------------------------------------------------------------------*/
double bathy_smooth_us(mesh_t *m, double *bathy, int cc) {
  int j, c;
  double kk = 0.0;
  double b = 0.0;

  for (j = 1; j <= m->npe[cc]; j++) {
    if ((c = m->neic[j][cc])) {
      if (bathy[c] != LANDCELL && bathy[c] != NOTVALID) {
	b += bathy[c];
	kk += 1.0;
      }
    }
  }
  if (kk) return(b / kk);
  return(bathy[cc]);
}

/* END bathy_smooth_us()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Smooths the bathymetry at a point for structured grids            */
/*-------------------------------------------------------------------*/
double smooth_bathy(unsigned long ***flag, int nz, int i, int j, int **xm, int **xp, 
		    int **ym, int **yp, double **bathy, double bmin, double bmax)
{
  double nbthy, kk[4][4];
  int im, ip, jm, jp;

  /* Set the filter weights */
  kk[1][3] = 0.0625;
  kk[2][3] = 0.125;
  kk[3][3] = 0.0625;
  kk[1][2] = 0.125;
  kk[2][2] = 0.25;
  kk[3][2] = 0.125;
  kk[1][1] = 0.0625;
  kk[2][1] = 0.125;
  kk[3][1] = 0.0625;


  if (!(flag[nz - 1][j][i] & (SOLID | OUTSIDE))) {
    double k0 = kk[2][2];
    double k1 = kk[1][2], k2 = kk[3][2];
    double k3 = kk[1][1], k4 = kk[2][1], k5 = kk[3][1];
    double k6 = kk[1][3], k7 = kk[2][3], k8 = kk[3][3];
    if (flag[nz - 1][j][i-1] & (SOLID | OUTSIDE)) {
      k0 += k1;
      k1 = 0.0;
    }
    if (flag[nz - 1][j][i+1] & (SOLID | OUTSIDE)) {
      k0 += k2;
      k2 = 0.0;
    }
    if (flag[nz - 1][j-1][i-1] & (SOLID | OUTSIDE)) {
      k0 += k3;
      k3 = 0.0;
    }
    if (flag[nz - 1][j-1][i] & (SOLID | OUTSIDE)) {
      k0 += k4;
      k4 = 0.0;
    }
    if (flag[nz - 1][j-1][i+1] & (SOLID | OUTSIDE)) {
      k0 += k5;
      k5 = 0.0;
    }
    if (flag[nz - 1][j+1][i-1] & (SOLID | OUTSIDE)) {
      k0 += k6;
      k6 = 0.0;
    }
    if (flag[nz - 1][j+1][i] & (SOLID | OUTSIDE)) {
      k0 += k7;
      k7 = 0.0;
    }
    if (flag[nz - 1][j+1][i+1] & (SOLID | OUTSIDE)) {
      k0 += k8;
      k8 = 0.0;
    }
    im = xm[j][i];
    ip = xp[j][i];
    jm = ym[j][i];
    jp = yp[j][i];
    /* Perform the convolution */
    nbthy = k6 * bathy[jp][im] + 
      k7 * bathy[jp][i] + 
      k8 * bathy[jp][ip] +
      k1 * bathy[j][im] + 
      k0 * bathy[j][i] + 
      k2 * bathy[j][ip] +
      k3 * bathy[jm][im] + 
      k4 * bathy[jm][i] + 
      k5 * bathy[jm][ip];
    nbthy = min(nbthy, -bmin);
    nbthy = max(nbthy, -bmax);
  }
  return (nbthy);
}

/* END smooth_bathy()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to apply a convolution smoothing filter to an element of  */
/* a 2-D array.                                                      */
/*-------------------------------------------------------------------*/
double
smoothbathy(double **a, int c, int p, int m, int pp, int mm, int cc,
            int mode, double bmin)
{
  double kk[5];
  double fi;

  /* Set the filter weights */
  kk[0] = 0.1;
  kk[1] = 0.2;
  kk[2] = 0.4;
  kk[3] = 0.2;
  kk[4] = 0.1;

  /* Perform the convolution */
  if (mode)
    fi = kk[0] * a[cc][mm] + kk[1] * a[cc][m] + kk[2] * a[cc][c] +
      kk[3] * a[cc][p] + kk[4] * a[cc][pp];
  else
    fi = kk[0] * a[mm][cc] + kk[1] * a[m][cc] + kk[2] * a[c][cc] +
      kk[3] * a[p][cc] + kk[4] * a[pp][cc];
  if (fi > -bmin)
    fi = -bmin;
  return (fi);
}

/* END smoothbathy()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Moves an open boundary into the interior                          */
/*-------------------------------------------------------------------*/
void adjust_outside_obc(parameters_t * params, double **bathy)
{
  int n, m, i, j, k, ef, m1, m2;
  int *iloc, *jloc, nb;

  for (n = 0; n < params->nobc; n++) {
    open_bdrys_t *open = params->open[n];

    if (open->bout < 0) continue;

    nb = 0;
    iloc = i_alloc_1d(open->npts);
    jloc = i_alloc_1d(open->npts);

    for (m = 0; m < open->npts; m++) {	  
      i = open->iloc[m];
      j = open->jloc[m];
      
      /* U1 boundaries                                               */
      if (open->type & U1BDRY) {
	ef = 0;
	if (i < params->nce1 && i > 0) {
	  if (!iswet(bathy[j][i])) ef = R_EDGE;
	  if (!iswet(bathy[j][i-1])) ef = L_EDGE;
	}
	if (i == params->nce1) ef = R_EDGE;
	if (i == 0) ef = L_EDGE;
	/* R_EDGE                                                    */
	if (ef == R_EDGE) {
	  for (k = 1; k <= open->bout; k++) {
	    i = (i > 0) ? i-1 : i;
	    if (iswet(bathy[j][i])) bathy[j][i] = NOTVALID;		  
	    open->iloc[m] = i;
	  }
	  m1 = (i > 0) ? i-1 : i;
	  if (iswet(bathy[j][m1])) {
	    /* Use below to enforce no OBCs 1 cell wide              */
	    m2 = (m1 > 0) ? m1-1 : m1;
	    if (!iswet(bathy[j][m2]))
	      bathy[j][m1] = NOTVALID;
	    else {
	      iloc[nb] = open->iloc[m];
	      jloc[nb] = open->jloc[m];
	      bathy[j][i] = NOTVALID;
	      nb++;
	    }
	  }
	}
	/* L_EDGE                                                    */
	if (ef == L_EDGE) {
	  for (k = 1; k <= open->bout; k++) {
	    if (iswet(bathy[j][i])) bathy[j][i] = NOTVALID;
	    i = (i < params->nce1) ? i+1 : i;
	    open->iloc[m] = i;
	  }
	  if (iswet(bathy[j][i])) {
	    /* Use below to enforce no OBCs 1 cell wide              */
	    m1 = (i < params->nce1) ? i+1 : i;
	    if (!iswet(bathy[j][m1]))
	      bathy[j][i] = NOTVALID;
	    else {
	      iloc[nb] = open->iloc[m];
	      jloc[nb] = open->jloc[m];
	      bathy[j][i-1] = NOTVALID;
	      nb++;}
	  }
	}
      }
      /* U2 boundaries                                               */
      if (open->type & U2BDRY) {
	ef = 0;
	if (j < params->nce2 && j > 0) {
	  if (!iswet(bathy[j][i])) ef = F_EDGE;
	  if (!iswet(bathy[j-1][i])) ef = B_EDGE;
	}
	if (j == params->nce2) ef = F_EDGE;
	if (j == 0) ef = B_EDGE;
	/* F_EDGE                                                    */
	if (ef == F_EDGE) {
	  for (k = 1; k <= open->bout; k++) {
	    j = (j > 0) ? j-1 : j;
	    if (iswet(bathy[j][i])) bathy[j][i] = NOTVALID;		  
	    open->jloc[m] = j;
	  }
	  m1 = (j > 0) ? j-1 : j;
	  if (iswet(bathy[m1][i])) {
	    m2 = (m1 > 0) ? m1-1 : m1;
	    if (!iswet(bathy[m2][i]))
	      bathy[m1][i] = NOTVALID;
	    else {
	      iloc[nb] = open->iloc[m];
	      jloc[nb] = open->jloc[m];
	      bathy[j][i] = NOTVALID;
	      nb++;
	    }
	  }
	}
	/* B_EDGE                                                    */
	if (ef == B_EDGE) {
	  for (k = 1; k <= open->bout; k++) {
	    if (iswet(bathy[j][i])) bathy[j][i] = NOTVALID;
	    j = (j < params->nce2) ? j+1 : j;
	    open->jloc[m] = j;
	  }
	  if (iswet(bathy[j][i])) {
	    m1 = (j < params->nce2) ? j+1 : j;
	    if (!iswet(bathy[m1][i]))
	      bathy[j][i] = NOTVALID;
	    else {
	      iloc[nb] = open->iloc[m];
	      jloc[nb] = open->jloc[m];
	      bathy[j-1][i] = NOTVALID;
	      nb++;
	    }
	  }
	}
      }
    }
    for (m = 0; m < nb; m++) {
      open->iloc[m] = iloc[m];
      open->jloc[m] = jloc[m];
      open->npts = nb;
    }
    i_free_1d(iloc);
    i_free_1d(jloc);
  }
  /* Omit the OBC points if the OBC cell is defined as NOTVALID above */
  for (n = 0; n < params->nobc; n++) {
    open_bdrys_t *open = params->open[n];
    int m, dof, nb = 0;
    if (open->bout < 0) continue;
    for (m = 0; m < open->npts; m++) {	  
      i = open->iloc[m];
      j = open->jloc[m];
      dof = 1;
      if (open->type & U1BDRY) {
	if (i > 0 && i < params->nce1 && bathy[j][i] == NOTVALID && 
	    bathy[j][i-1] == NOTVALID) 
	  dof = 0;
	else if (i == 0 && bathy[j][i] == NOTVALID)
	  dof = 0;
	else if (i == params->nce1 && bathy[j][i-1] == NOTVALID)
	  dof = 0;
      }
      if (open->type & U2BDRY) {
	if (j > 0 && j < params->nce2 && bathy[j][i] == NOTVALID && 
	    bathy[j-1][i] == NOTVALID) 
	  dof = 0;
	else if (j == 0 && bathy[j][i] == NOTVALID)
	  dof = 0;
	else if (j == params->nce2 && bathy[j-1][i] == NOTVALID)
	  dof = 0;
      }
      if (dof) {
	open->iloc[nb] = i;
	open->jloc[nb] = j;
	nb++;
      }
    }
    open->npts = nb;
    if (open->bout >= 0)
      hd_warn("Boundary%d (%s) new RANGE (%d,%d)-(%d,%d)\n",
	      n,open->name,open->iloc[0],open->jloc[0],
	      open->iloc[nb-1],open->jloc[nb-1]);
  }
}

/* END adjust_outside_obc()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads the sediment layer configuration                            */
/*-------------------------------------------------------------------*/
void read_sed_layers(geometry_t *geom,      /* Global geometry       */
		     parameters_t *params,  /* Input parameters      */
		     master_t *master,      /* Master data           */
		     dump_data_t *dumpdata  /* Dump data             */
		     )
{
  FILE *fp = params->prmfd;
  int i, j, k, m;

  /*-----------------------------------------------------------------*/
  /* Read the sediment layer structure from file                     */
  /* Already read in in read_params()
  if (params->gridz_sed) d_free_1d(params->gridz_sed);
  if (prm_read_darray(fp, "LAYERFACES_SED", &params->gridz_sed,
                      &params->sednz)) {
    if (--params->sednz < 1)
      hd_quit("Number of sediment layers must be 2 or more\n");
  } else {
    double *dz_sed = NULL;
    if (prm_read_int(fp, "NSEDLAYERS", &params->sednz)) {
      dz_sed = d_alloc_1d(params->sednz);
      if (prm_read_darray(fp, "NSEDLAYERS", &dz_sed, &params->sednz)) {
        if (params->sednz) {
          params->gridz_sed = d_alloc_1d(params->sednz + 1);
          params->gridz_sed[params->sednz] = 0.0;
          for (m = params->sednz - 1; m >= 0; m--) {
            params->gridz_sed[m] = params->gridz_sed[m + 1] -
              dz_sed[params->sednz - m - 1];
	  }
        }
      }
      if(dz_sed)
        d_free_1d(dz_sed);
    }
  }
  */

  /*-----------------------------------------------------------------*/
  /* Transfer to the dumpdata variables                              */
  for (j = 0; j < geom->nce2; j++)
    for (i = 0; i < geom->nce1; i++) {
        for (k = 0; k < params->sednz; k++) {
	  dumpdata->gridz_sed[k][j][i] = params->gridz_sed[k];
	  dumpdata->cellz_sed[k][j][i] = 0.5 * (params->gridz_sed[k] +
						params->gridz_sed[k + 1]);
	}
	dumpdata->gridz_sed[params->sednz][j][i] = params->gridz_sed[params->sednz];
    }

  /*-----------------------------------------------------------------*/
  /* Transfer to the window geometry variables                       */
  for (k = 0; k < params->sednz + 1; k++) {
    c2s_2d(geom, geom->gridz_sed[k], dumpdata->gridz_sed[k], geom->nce1,
	   geom->nce2);
  }
  for (k = 0; k < params->sednz; k++) {
    c2s_2d(geom, geom->cellz_sed[k], dumpdata->cellz_sed[k], geom->nce1,
	   geom->nce2);
  }
}

/* END read_sed_layers()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Read a 3d cartesian array from id, zeroing it if an error occurs  */
/* and copying it into the sparse array.                             */
/*-------------------------------------------------------------------*/
int read3d(geometry_t *geom, int id, char *name, double *p, int dump,
            int nk, int nj, int ni)
{
  double ***array = d_alloc_3d(ni, nj, nk); /* Dummy cartesian array */
  int c;                        /* Sparse couters */
  int i, j, k;                  /* Cartesian couters */

  size_t start[4];
  size_t count[4];

  start[0] = dump;
  start[1] = 0L;
  start[2] = 0L;
  start[3] = 0L;
  count[0] = 1L;
  count[1] = nk;
  count[2] = nj;
  count[3] = ni;
  if (nc_get_vara_double
      (id, ncw_var_id(id, name), start, count, array[0][0]) < 0) {
    memset(p, 0, geom->sgsiz * sizeof(double));
    return(1);
  }

  /* Copy into the sparse array */
  for (c = 1; c <= geom->sgnum; c++) {
    i = geom->s2i[c];
    j = geom->s2j[c];
    k = geom->s2k[c];
    if (i == NOTVALID && j == NOTVALID && k == NOTVALID)
      continue;                 /* Continue on ghosts */
    if (i < ni && j < nj)
      p[c] = array[k][j][i];
  }
  d_free_3d(array);
  return(0);
}

/* END read3d()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Removes channels one cell wide from the grid                      */
/*-------------------------------------------------------------------*/
void remove_channel(parameters_t *params, unsigned long ***flag, double **bathy)
{
  int nce1 = params->nce1;
  int nce2 = params->nce2;
  int nz = params->nz;
  int i, j, k;
  int ip, im, jp, jm;

  /* Identify one cell channels or canyons */
  for (k = 0; k < nz - 1; k++) {
    for (j = 0; j < nce2; j++) {
      for (i = 0; i < nce1; i++) {
	if (!(flag[k][j][i] & (SOLID | OUTSIDE))) {
	  ip = (i < nce1 - 1) ? i+1 : i;
	  im = (i > 0) ? i-1 : i;
	  jp = (j < nce2 - 1) ? j+1 : j;
	  jm = (j > 0) ? j-1 : j;

          if (flag[k][j][im] & (SOLID | OUTSIDE) &&
              flag[k][j][ip] & (SOLID | OUTSIDE) &&
              !(flag[k][jm][i] & (SOLID | OUTSIDE)) &&
              !(flag[k][jp][i] & (SOLID | OUTSIDE)))
	    hd_warn("remove_channel: one-grid channel TYPEjj, %d %d %d \n",i, j, k);
          else if (!(flag[k][j][im] & (SOLID | OUTSIDE)) &&
                   !(flag[k][j][ip] & (SOLID | OUTSIDE)) &&
                   flag[k][jm][i] & (SOLID | OUTSIDE) &&
                   flag[k][jp][i] & (SOLID | OUTSIDE))
	    hd_warn("remove_channel: one-grid channel TYPEii, %d %d %d\n",i, j, k);  
          else if (flag[k][j][im] & (SOLID | OUTSIDE) &&
                   flag[k][j][ip] & (SOLID | OUTSIDE) &&
                   !(flag[k][jm][i] & (SOLID | OUTSIDE)) &&
                   flag[k][jp][i] & (SOLID | OUTSIDE))
	    hd_warn("remove_channel: one-grid channel TYPEjm, %d %d %d\n",i,j, k);
          else if (flag[k][j][im] & (SOLID | OUTSIDE) &&
                   flag[k][j][ip] & (SOLID | OUTSIDE) &&
                   flag[k][jm][i] & (SOLID | OUTSIDE) &&
                   !(flag[k][jp][i] & (SOLID | OUTSIDE))) 
	    hd_warn("remove_channel: one-grid channel TYPEjp, %d %d %d\n",i, j, k);
          else if (!(flag[k][j][im] & (SOLID | OUTSIDE)) &&
		   flag[k][j][ip] & (SOLID | OUTSIDE) &&
                   flag[k][jm][i] & (SOLID | OUTSIDE) &&
                   flag[k][jp][i] & (SOLID | OUTSIDE))
	    hd_warn("remove_channel: one-grid channel TYPEim, %d %d %d\n",i, j, k); 
          else if ( flag[k][j][im] & (SOLID | OUTSIDE) &&
		    !(flag[k][j][ip] & (SOLID | OUTSIDE)) &&
		    flag[k][jm][i] & (SOLID | OUTSIDE) &&
		    flag[k][jp][i] & (SOLID | OUTSIDE))
	    hd_warn("remove_channel: one-grid channel TYPEip, %d %d %d\n",i, j, k);
	}
      }
    }
  }

  /* Fill in channels or canyons one cell wide at the surface */
  for (j = 0; j < nce2; j++) {
    for (i = 0; i < nce1; i++) {
      if (!(flag[nz-1][j][i] & (SOLID | OUTSIDE))) {
	ip = (i < nce1 - 1) ? i+1 : i;
	im = (i > 0) ? i-1 : i;
	jp = (j < nce2 - 1) ? j+1 : j;
	jm = (j > 0) ? j-1 : j;
	
	if (flag[nz-1][j][im] & (SOLID | OUTSIDE) &&
	    flag[nz-1][j][ip] & (SOLID | OUTSIDE) &&
	    !(flag[nz-1][jm][i] & (SOLID | OUTSIDE)) &&
	    !(flag[nz-1][jp][i] & (SOLID | OUTSIDE))) {
	  bathy[j][i] = LANDCELL;
	  for (k = 0; k < nz; k++) {
	    flag[k][j][i] |= SOLID;
	  }
	}
	else if (!(flag[nz-1][j][im] & (SOLID | OUTSIDE)) &&
		 !(flag[nz-1][j][ip] & (SOLID | OUTSIDE)) &&
		 flag[nz-1][jm][i] & (SOLID | OUTSIDE) &&
		 flag[nz-1][jp][i] & (SOLID | OUTSIDE)) {
	  bathy[j][i] = LANDCELL;
	  for (k = 0; k < nz; k++) {
	    flag[k][j][i] |= SOLID;
	  }
	}
	else if (flag[nz-1][j][im] & (SOLID | OUTSIDE) &&
		 flag[nz-1][j][ip] & (SOLID | OUTSIDE) &&
		 !(flag[nz-1][jm][i] & (SOLID | OUTSIDE)) &&
		 flag[nz-1][jp][i] & (SOLID | OUTSIDE)) {
	  bathy[j][i] = LANDCELL;
	  for (k = 0; k < nz; k++) {
	    flag[k][j][i] |= SOLID;
	  }
	}
	else if (flag[nz-1][j][im] & (SOLID | OUTSIDE) &&
		 flag[nz-1][j][ip] & (SOLID | OUTSIDE) &&
		 flag[nz-1][jm][i] & (SOLID | OUTSIDE) &&
		 !(flag[nz-1][jp][i] & (SOLID | OUTSIDE))) {
	  bathy[j][i] = LANDCELL;
	  for (k = 0; k < nz; k++) {
	    flag[k][j][i] |= SOLID;
	  }
	}
	else if (!(flag[nz-1][j][im] & (SOLID | OUTSIDE)) &&
		 flag[nz-1][j][ip] & (SOLID | OUTSIDE) &&
		 flag[nz-1][jm][i] & (SOLID | OUTSIDE) &&
		 flag[nz-1][jp][i] & (SOLID | OUTSIDE)) {
	  bathy[j][i] = LANDCELL;
	  for (k = 0; k < nz; k++) {
	    flag[k][j][i] |= SOLID;
	  }
	}
	else if ( flag[nz-1][j][im] & (SOLID | OUTSIDE) &&
		  !(flag[nz-1][j][ip] & (SOLID | OUTSIDE)) &&
		  flag[nz-1][jm][i] & (SOLID | OUTSIDE) &&
		  flag[nz-1][jp][i] & (SOLID | OUTSIDE)) {
	  bathy[j][i] = LANDCELL;
	  for (k = 0; k < nz; k++) {
	    flag[k][j][i] |= SOLID;
	  }
	}
      }
    }
  }
}

/* END remove_channel()                                              */
/*-------------------------------------------------------------------*/


int iswet(int bathy)
{
  if (bathy == NOTVALID || bathy == LANDCELL)
    return 0;
  return 1;
}
