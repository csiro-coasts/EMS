/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/inputs/opendump.c
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
 *  $Id: opendump.c 6960 2021-12-17 03:17:27Z her127 $
 *
 */

#include <string.h>
#include <netcdf.h>
#include "hd.h"


int dump_open(parameters_t *params, char *name, int in_model)
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
  size_t jcentresize;
  size_t icentresize;
  size_t kgridsize;
  size_t jgridsize;
  size_t igridsize;
  size_t jleftsize;
  size_t ileftsize;
  size_t jbacksize;
  size_t ibacksize;
  size_t kcentresize_sed;
  size_t kgridsize_sed;
  int sednz = params->sednz;

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
  nc_inq_dimlen(fid, ncw_dim_id(fid, "k_grid"), &kgridsize);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "j_grid"), &jgridsize);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "i_grid"), &igridsize);

  nc_inq_dimlen(fid, ncw_dim_id(fid, "k_centre"), &kcentresize);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "j_centre"), &jcentresize);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "i_centre"), &icentresize);

  nc_inq_dimlen(fid, ncw_dim_id(fid, "j_left"), &jleftsize);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "i_left"), &ileftsize);

  nc_inq_dimlen(fid, ncw_dim_id(fid, "j_back"), &jbacksize);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "i_back"), &ibacksize);
  if (kgridsize != params->nz + 1) {
    hd_warn
      ("dump_open: Different number of layers in dump file (%d) and parameter file (%d)\n", kgridsize, params->nz + 1);
  }
  params->nz = (long)kcentresize;

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

  nc_get_att_int(fid, NC_GLOBAL, "nce1", &params->nce1);
  nc_get_att_int(fid, NC_GLOBAL, "nfe1", &params->nfe1);
  nc_get_att_int(fid, NC_GLOBAL, "nce2", &params->nce2);
  nc_get_att_int(fid, NC_GLOBAL, "nfe2", &params->nfe2);
  if ((params->nce2 + 1) != (int)jgridsize ||
      (params->nce1 + 1) != (int)igridsize)
    hd_quit
      ("dump_open: Dimensions not compatible with numbers of cells/faces (%d !=%d, %d!=%d)\n");


  read_grid_atts(params, fid);

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
    nc_get_att_text(fid, ncw_var_id(fid, "z_grid"), "units",
                    params->lenunit);

  /* Check if geographic. */
  nc_get_att_text(fid, ncw_var_id(fid, "x_grid"), "projection", buf);
  if (buf != NULL) {
    if ((strlen(projection) > 0) && (strcasecmp(projection, buf) != 0))
      hd_quit("Input dump file has different projection type");
    strcpy(projection, buf);
    ts_set_default_proj_type(projection);
  }

  /* Suggest the minimum hash table size for evaluation of timeseries
     numeric grids. */
  ts_set_default_hashtable_size((params->nce1 + 1) * (params->nce2 +
                                                      1) * 4);

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
    nc_get_vara_double(fid, ncw_var_id(fid, "z_grid"), start, count,
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
  return (fid);
}
