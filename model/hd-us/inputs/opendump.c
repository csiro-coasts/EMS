/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/inputs/opendump.c
 *  
 *  Description:
 *  Routine to open model netCDF dump file for reading
 *  This routine reads the file attributes.
 *  dump_read() should then be used to read the data.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: opendump.c 5873 2018-07-06 07:23:48Z riz008 $
 *
 */

#include <string.h>
#include <netcdf.h>
#include "hd.h"

void create_mesh(parameters_t *params, int fid);
void convert_structured_obc(parameters_t *params);


/*-------------------------------------------------------------------*/
/* Reads in grid information for a run from an unstructured UGRID    */
/* input file.                                                       */
/*-------------------------------------------------------------------*/
int dump_open_us(parameters_t *params, char *name, int in_model)
{
  int i;
  int fid;
  int ncerr;
  char vers[MAXSTRLEN];
  char chead[MAXSTRLEN];
  char phead[MAXSTRLEN];
  char prmname[MAXSTRLEN];
  char buf[MAXSTRLEN];
  size_t kcentresize;
  size_t kgridsize;
  size_t nMesh2_node;
  size_t nMesh2_edge;
  size_t nMesh2_face;
  size_t kcentresize_sed;
  size_t kgridsize_sed;
  int sednz = params->sednz;
  int nface, nedge, nvertex, nmax;

  /* clear the string arrays in case data in file isn't zero terminated */
  for (i = 0; i < MAXSTRLEN; i++) {
    vers[i] = 0;
    chead[i] = 0;
    phead[i] = 0;
    buf[i] = 0;
  }

  /* open the dump file for reading */
  if ((ncerr = nc_open(name, NC_NOWRITE, &fid)) != NC_NOERR) {
    printf("Can't find input file %s\n", name);
    hd_quit((char *)nc_strerror(ncerr));
  }

  /* get dimensions */
  nc_inq_dimlen(fid, ncw_dim_id(fid, "Mesh2_layerfaces"), &kgridsize);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "Mesh2_layers"), &kcentresize);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "nMesh2_node"), &nMesh2_node);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "nMesh2_edge"), &nMesh2_edge);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "nMesh2_face"), &nMesh2_face);
  if (kgridsize != params->nz + 1) {
    hd_warn
      ("dump_open: Different number of layers in dump file (%d) and parameter file (%d)\n", kgridsize, params->nz + 1);
  }
  params->nz = (long)kcentresize;
  nc_get_att_int(fid, NC_GLOBAL, "start_index", &params->oset);
  /*params->oset = atoi(buf);*/
  params->oset = (params->oset) ? 0 : 1;

  if (nc_inq_dimlen(fid, ncw_dim_id(fid, "k_grid_sed"), &kgridsize_sed)
      || nc_inq_dimlen(fid, ncw_dim_id(fid, "k_centre_sed"),
                       &kcentresize_sed))
  {
    /*params->sednz = 0;*/
#if defined(HAVE_SEDIMENT_MODULE)
    if(params->do_sed)
      emstag(LWARN,"hd:opendump:dump_open","Grid var missing for sediment, using %d sediment layers specified in parameter file.\n", sednz);
    /*emstag(LWARN,"hd:opendump:dump_open","Grid var missing for sediment, likely cause is wrong input file, rebuild!");*/
#endif	
  }
  else
    params->sednz = kcentresize_sed;

  if (sednz != params->sednz) {
    hd_warn
      ("dump_open: Different number of sediment layers in dump file and parameter file\n");
    params->sednz = 0;
  }

  nc_get_att_int(fid, NC_GLOBAL, "nface2", &nface);
  nc_get_att_int(fid, NC_GLOBAL, "nedge2", &nedge);
  nc_get_att_int(fid, NC_GLOBAL, "nvertex2", &nvertex);

  if (nvertex != (int)nMesh2_node ||
      nface != (int)nMesh2_face ||
      nedge != (int)nMesh2_edge)
    hd_quit
      ("dump_open: Dimensions not compatible with numbers of cells/faces\n");
  read_grid_atts(params, fid);
  nMesh2_node += params->oset;
  nMesh2_face += params->oset;
  nMesh2_edge += params->oset;
  nface += params->oset;
  nedge += params->oset;
  nvertex += params->oset;

  if (in_model) {
    /* get global attributes */
    nc_get_att_text(fid, NC_GLOBAL, "title", chead);
    nc_get_att_text(fid, NC_GLOBAL, "paramhead", phead);
    nc_get_att_text(fid, NC_GLOBAL, "paramfile", prmname);
    nc_get_att_text(fid, NC_GLOBAL, "version", vers);

    /* check compatibility */
    if (strcmp(version, vers) != 0)
      hd_warn("Input dump file version doesn't match program version\n");
    if (strcmp(codeheader, chead) != 0)
      hd_warn
        ("Parameter file codeheader >%s<\n doesn't match input dump file title >%s<\n",
         codeheader, chead);
    if (strcmp(parameterheader, phead) != 0)
      hd_warn
        ("Parameter file and input dump file parameterheader don't match\n");
  } else {
    /* get global attributes */
    memset(codeheader, 0, MAXSTRLEN);
    memset(parameterheader, 0, MAXSTRLEN);
    memset(prmname, 0, MAXSTRLEN);
    memset(vers,    0, MAXSTRLEN);
    nc_get_att_text(fid, NC_GLOBAL, "title", codeheader);
    nc_get_att_text(fid, NC_GLOBAL, "paramhead", parameterheader);
    nc_get_att_text(fid, NC_GLOBAL, "paramfile", prmname);
    nc_get_att_text(fid, NC_GLOBAL, "version", vers);
    if (strcmp(version, vers) != 0)
      hd_warn("Input dump file version doesn't match program version\n");
  }

  if (strlen(params->timeunit) == 0)
    nc_get_att_text(fid, ncw_var_id(fid, "t"), "units", params->timeunit);
  if (strlen(params->lenunit) == 0)
    nc_get_att_text(fid, ncw_var_id(fid, "Mesh2_layerfaces"), "units",
                    params->lenunit);

  /* Check if geographic. */
  nc_get_att_text(fid, ncw_var_id(fid, "Mesh2_node_x"), "projection", buf);
  if (buf != NULL) {
    if ((strlen(projection) > 0) && (strcasecmp(projection, buf) != 0))
      hd_quit("Input dump file has different projection type");
    strcpy(projection, buf);
    ts_set_default_proj_type(projection);
  }

  /* Suggest the minimum hash table size for evaluation of timeseries
     numeric grids. */
  nmax = max(nface, nedge);
  nmax = max(nmax, nvertex);
  ts_set_default_hashtable_size((nmax + 1) * 4);

  /* time independent variables */
  /* If the vertical layer structure has not been read in, then 
     initialise from the netCDF file. */
  if (params->layers == NULL) {
    size_t start[4];
    size_t count[4];
    params ->layers = d_alloc_1d(params->nz + 1);
    count[0] = params->nz + 1;
    start[0] = 0;
    count[1] = start[1] = 0;
    count[2] = start[2] = 0;
    count[3] = start[3] = 0;
    nc_get_vara_double(fid, ncw_var_id(fid, "Mesh2_layerfaces"), start, count,
		       params->layers);
  }

  /* Read the tracer descriptions from the dump file, but only if called
     from a utility program (as the model itself reads this information
     from the parameter file). */
  if (!in_model) {
    tracer_att_t attr = { "tracer", "" };
    tracer_read_nc(fid, 4, NULL, &attr, &params->ntr, &params->trinfo_3d);
  }
  if (!in_model) {
    tracer_att_t attr = { "tracer2D", "" };
    tracer_read_nc(fid, 3, NULL, &attr, &params->ntrS, &params->trinfo_2d);
  }
  if (!in_model) {
    tracer_att_t attr = { "tracer_sed", "" };
    tracer_read_nc(fid, 4, NULL, &attr, &params->nsed, &params->trinfo_sed);
  }

  create_mesh(params, fid);

  return (fid);
}

/* END dump_open_us()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Creates a mesh structure from information contained in the input  */
/* netCDF file.                                                      */
/*-------------------------------------------------------------------*/
void create_mesh(parameters_t *params, int fid)
{
  int c, cc, v, n, nn;
  size_t start[4];
  size_t count[4];
  size_t nMesh2_face;
  size_t nMesh2_node;
  size_t nMaxMesh2_face_nodes;
  int *index, *c2i = NULL, *c2j = NULL;
  int **ic2v, **c2v, nce1 = 0, nce2 = 0, fill;
  double *cellx, *celly;
  double *gridx, *gridy;
  double *bathy;
  int istart = (params->oset) ? 0 : 1;
  int face_dim = 0;
  int npem;
  int verbose = 0;

  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;

  count[0] = 0;
  count[1] = 0;
  count[2] = 0;
  count[3] = 0;

  /*-----------------------------------------------------------------*/
  /* Get dimensions and allocate                                     */
  nc_inq_dimlen(fid, ncw_dim_id(fid, "nMesh2_face"), &nMesh2_face);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "nMesh2_node"), &nMesh2_node);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "nMaxMesh2_face_nodes"), &nMaxMesh2_face_nodes);
  nc_get_att_int(fid, NC_GLOBAL, "face_dim", &face_dim);

  c2v = i_alloc_2d(nMesh2_face, nMaxMesh2_face_nodes);
  if (face_dim)
    ic2v = i_alloc_2d(nMesh2_face, nMaxMesh2_face_nodes);
  else
    ic2v = i_alloc_2d(nMaxMesh2_face_nodes, nMesh2_face);
  index = i_alloc_1d(nMesh2_face);
  c2i = i_alloc_1d(nMesh2_face);
  c2j = i_alloc_1d(nMesh2_face);
  cellx = d_alloc_1d(nMesh2_face);
  celly = d_alloc_1d(nMesh2_face);
  gridx = d_alloc_1d(nMesh2_node);
  gridy = d_alloc_1d(nMesh2_node);
  bathy = d_alloc_1d(nMesh2_face);

  params->ns2 = nMesh2_face - istart;
  npem = params->npe = nMaxMesh2_face_nodes - istart;
  params->npe2 = i_alloc_1d(params->ns2 + 1);
  params->x = d_alloc_2d(npem + 1, params->ns2 + 1);
  params->y = d_alloc_2d(npem + 1, params->ns2 + 1);
  params->bathy = d_alloc_1d(params->ns2 + 1);

  /*-----------------------------------------------------------------*/
  /* Get the locations of the cell centre                            */
  count[0] = nMesh2_face;
  nc_get_vara_double(fid, ncw_var_id(fid, "Mesh2_face_x"), start, count,
                     cellx);
  nc_get_vara_double(fid, ncw_var_id(fid, "Mesh2_face_y"), start, count,
                     celly);
  nc_get_vara_double(fid, ncw_var_id(fid, "Mesh2_depth"), start, count,
                     bathy);
  nc_get_vara_int(fid, ncw_var_id(fid, "Mesh2_index"), start, count,
		  index);

  /* Get the indices for structured grids                            */
  nc_get_att_int(fid, NC_GLOBAL, "NCE1", &nce1);
  nc_get_att_int(fid, NC_GLOBAL, "NCE2", &nce2);
  if (nce1 > 0 && nce2 > 0) {
    if (autof == 8 || autof == 9) {
      params->nce1 = nce1;
      params->nce2 = nce2;
    } else
      if (nce1 != params->nce1 || nce2 != params->nce2)
	hd_quit("create_mesh: Parameter file dimensions not compatible with netCDF file dimensions: nce1=%d nce2=%d (%d %d)\n",
		nce1, nce2, params->nce1, params->nce2);
    nc_get_vara_int(fid, ncw_var_id(fid, "Mesh2_iindex"), start, count, c2i);
    nc_get_vara_int(fid, ncw_var_id(fid, "Mesh2_jindex"), start, count, c2j);
    params->us_type |= US_IJ;
    params->ic = i_alloc_1d(params->ns2+1);
    params->jc = i_alloc_1d(params->ns2+1);
  }

  /*-----------------------------------------------------------------*/
  /* Get the locations of the cell vertices                          */
  count[0] = nMesh2_node;
  nc_get_vara_double(fid, ncw_var_id(fid, "Mesh2_node_x"), start, count,
                     gridx);
  nc_get_vara_double(fid, ncw_var_id(fid, "Mesh2_node_y"), start, count,
                     gridy);

  /*-----------------------------------------------------------------*/
  /* Get the mapping from centre to vertices and populate the        */
  /* parameters structure.                                           */
  /* Note that the input of the mesh topology to create the input    */
  /* file may list the mesh elements in arbitary order, and this     */
  /* will determine which faces are neighbours in the cell arrays,   */
  /* and which edges are associated with which cells. This in turn   */
  /* will determine on which cells or faces data read in from the    */
  /* input netCDF file is placed. Therefore, make sure to put the    */
  /* grid locations and bathymetry in the correct order via index[]  */
  /* so that the same mappings are created as was used to generate   */
  /* the input netCDF file.                                          */
  set_count(face_dim, count, nMaxMesh2_face_nodes, nMesh2_face);
  nc_get_vara_int(fid, ncw_var_id(fid, "Mesh2_face_nodes"), start, count,
		   ic2v[0]);
  nc_get_att_int(fid, ncw_var_id(fid, "Mesh2_face_nodes"), "_FillValue", &fill);

  for (cc = istart; cc < nMesh2_face; cc++)
    for (n = istart; n < nMaxMesh2_face_nodes; n++) {
      if (face_dim)
	c2v[n][cc] = ic2v[n][cc] + params->oset;
      else
	c2v[n][cc] = ic2v[cc][n] += params->oset;
    }
  for (cc = istart; cc < nMesh2_face; cc++) {
    c = index[cc];
    params->npe2[c] = 0;
    params->x[c][0] = cellx[cc];
    params->y[c][0] = celly[cc];
    params->bathy[c] = bathy[cc];
    if (params->us_type & US_IJ) {
      params->ic[c] = c2i[cc];
      params->jc[c] = c2j[cc];
    }
    if(isnan(params->bathy[c])) params->bathy[c] = NOTVALID;
    for (n = istart; n < nMaxMesh2_face_nodes; n++) {
      nn = n + params->oset;
      if (c2v[n][cc] > 0 && c2v[n][cc] != fill + params->oset) {
	params->npe2[c]++;
	v = c2v[n][cc] - params->oset;
	params->x[c][nn] = gridx[v];
	params->y[c][nn] = gridy[v];
      }
    }
  }

  if (verbose) {
    for (cc = 1; cc <= params->ns2; cc++) {
      if (params->us_type & US_IJ)
	printf("%d %d %f %f : %f : %d(%d %d)\n",cc,params->npe2[cc],params->x[cc][0],params->y[cc][0], 
	       params->bathy[cc], index[cc-params->oset], 
	       params->ic[cc], params->jc[cc]);
      else
	printf("%d %d %f %f : %f : %d\n",cc,params->npe2[cc],params->x[cc][0],params->y[cc][0], 
	       params->bathy[cc], index[cc]);
      for (n = 1; n <= npem; n++) {
	printf(" %d %d %f %f\n",n,params->npe2[cc],params->x[cc][n],params->y[cc][n]);
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Convert open boundary information for structured grids          */
  convert_structured_obc(params);

  /*-----------------------------------------------------------------*/
  /* Create the mesh structure                                       */
  meshstruct_us(params);

  /*-----------------------------------------------------------------*/
  /* Set the indices for structured grids                            */
  if (params->us_type & US_IJ) {
    mesh_t *m = params->mesh;
    for (cc = 1; cc <= m->ns2 ; cc++) {
      m->iloc[cc] = params->ic[cc];
      m->jloc[cc] = params->jc[cc];
    }
  }

  i_free_1d(index);
  i_free_1d(c2i);
  i_free_1d(c2j);
  i_free_2d(c2v);
  i_free_2d(ic2v);
  d_free_1d(cellx);
  d_free_1d(celly);
  d_free_1d(gridx);
  d_free_1d(gridy);
  d_free_1d(bathy);
}

/* END create_mesh()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Converts a structured (i,j) open boundary specification to the    */
/* locations required by the mesh structure to build its OBC         */
/* specification.                                                    */
/*-------------------------------------------------------------------*/
void convert_structured_obc(parameters_t *params)
{
  int nn, n, i, j, cc;
  int nce1 = params->nce1;
  int nce2 = params->nce2;
  double **bathy;

  if (params->us_type & US_IJ) {

    /* Make a structured bathymetry array                            */
    bathy = d_alloc_2d(nce1, nce2);
    for (j = 0; j < params->nce2; j++)
      for (i = 0; i < params->nce1; i++)
	bathy[j][i] = NOTVALID;
    for (cc = 1; cc <= params->ns2; cc++) {
      i = params->ic[cc];
      j = params->jc[cc];
      bathy[j][i] = params->bathy[cc];
    }

    for (n = 0; n < params->nobc; n++) {
      open_bdrys_t *open = params->open[n];

      if (open->locu && open->posx && open->posy) continue;
      open->locu = i_alloc_1d(open->npts);
      open->posx = d_alloc_2d(2 ,open->npts);
      open->posy = d_alloc_2d(2 ,open->npts);

      /* The params->open OBC specification is set here, where the   */
      /* locations are read in with POINTS or RANGE (i.e. intype =   */
      /* O_POI or O_RAN at this stage). These will be converted to   */
      /* the mesh structure in convert_obc_list() called from        */
      /* meshstruct_us(), but since we set the OBC specification to  */
      /* a coordinate type, we must set intype to O_UPC here so that */
      /* the conversion is handled correctly in convert_obc_list().  */
      open->intype = O_UPC;
    
      for (nn = 0; nn < open->npts; nn++) {

	i = open->iloc[nn];
	j = open->jloc[nn];

	/*
	if (open->type & U1BDRY && 
	    (i == 0 || (i > 0 && bathy[j][i-1] < params->layers[0])))
	  open->type |= L_EDGE;
	if (open->type & U2BDRY && 
	    (j == 0 || (j > 0 && bathy[j-1][i] < params->layers[0])))
	  open->type |= B_EDGE;
	if (open->type & U1BDRY && 
	    (i == nce1-1 || (i < nce1-1 && bathy[j][i] < params->layers[0])))
	  open->type |= R_EDGE;
	if (open->type & U2BDRY && 
	    (j == nce2-1 || (j < nce2-1 && bathy[j+1][i] < params->layers[0])))
	  open->type |= F_EDGE;
	*/
	for (cc = 1; cc <= params->ns2; cc++) {
	  if (i == params->ic[cc] && j == params->jc[cc]) {
	    open->locu[nn] = cc;
            /* L_EDGE                                                */
	    if (open->type & U1BDRY && (i == 0 || (i > 0 && bathy[j][i-1] < params->layers[0]))) {
	      open->posx[nn][0] = params->x[cc][1];
	      open->posy[nn][0] = params->y[cc][1];
	      open->posx[nn][1] = params->x[cc][2];
	      open->posy[nn][1] = params->y[cc][2];
	    }
            /* B_EDGE                                                */
	    if (open->type & U2BDRY && (j == 0 || (j > 0 && bathy[j-1][i] < params->layers[0]))) {
	      open->posx[nn][0] = params->x[cc][1];
	      open->posy[nn][0] = params->y[cc][1];
	      open->posx[nn][1] = params->x[cc][4];
	      open->posy[nn][1] = params->y[cc][4];
	    }
	  }

	  if (i - 1 == params->ic[cc] && j == params->jc[cc]) {
	    open->locu[nn] = cc;
            /* R_EDGE                                                */
	    if (open->type & U1BDRY && (i == nce1 || (i < nce1 && bathy[j][i+1] < params->layers[0]))) {
	      open->posx[nn][0] = params->x[cc][3];
	      open->posy[nn][0] = params->y[cc][3];
	      open->posx[nn][1] = params->x[cc][4];
	      open->posy[nn][1] = params->y[cc][4];
	    }
	  }
	  if (i == params->ic[cc] && j - 1 == params->jc[cc]) {
	    open->locu[nn] = cc;
            /* F_EDGE                                                */
	    if (open->type & U2BDRY && (j == nce2 || (j < nce2 && bathy[j+1][i] < params->layers[0]))) {
	      open->posx[nn][0] = params->x[cc][2];
	      open->posy[nn][0] = params->y[cc][2];
	      open->posx[nn][1] = params->x[cc][3];
	      open->posy[nn][1] = params->y[cc][3];
	    }
	  }
	}
      }
    }
    d_free_2d(bathy);
  } else {
  }
}

/* END convert_structured_obc()                                      */
/*-------------------------------------------------------------------*/
