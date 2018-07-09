/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/outputs/createdump.c
 *  
 *  Description:
 *  Create a netCDF dump file for the meco
 *  model
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: createdump.c 5873 2018-07-06 07:23:48Z riz008 $
 *
 */

#include <netcdf.h>
#include <string.h>
#include "hd.h"
#include "tracer.h"


int dump_create(dump_data_t *dumpdata, char *name)
{
  int n;
  int cdfid;
  nc_type fptype = NC_DOUBLE;
  int dims[10];
  int tid;
  int vid;
  /* dimension ids */
  int recdimid;                 /* record dimension id */
  int igridid;                  /* I dimension id at grid corner */
  int jgridid;                  /* J dimension id at grid corner */
  int kgridid;                  /* K dimension id at grid corner */
  int icentreid;                /* I dimension id at grid centre */
  int jcentreid;                /* J dimension id at grid centre */
  int kcentreid;                /* K dimension id at grid centre */
  int ileftid;                  /* I dimension id at left face */
  int jleftid;                  /* J dimension id at left face */
  int ibackid;                  /* I dimension id at back face */
  int jbackid;                  /* J dimension id at back face */
  int kcentreid_sed;            /* K dimension id at grid centre for
                                   sediments */
  int kgridid_sed;              /* K dimension id at grid corner for
                                   sediments */
  int nc_mode;

  /* Check the clobber flag */
  nc_mode = (overwrite_output() ? NC_CLOBBER : NC_NOCLOBBER);

  /* Make 64bit classical the default */
  if (!(params->compatible & V4201))
    nc_mode |= NC_64BIT_OFFSET;

  /* create the netCDF file NC_64BIT_OFFSET */
  if (ncw_create(name, nc_mode, &cdfid) != NC_NOERR)
    hd_quit("Couldn't create output dump file %s\n", name);

  /* define dimensions */
  nc_def_dim(cdfid, "record", NC_UNLIMITED, &recdimid);

  nc_def_dim(cdfid, "k_grid", dumpdata->nz + 1, &kgridid);
  nc_def_dim(cdfid, "j_grid", dumpdata->nce2 + 1, &jgridid);
  nc_def_dim(cdfid, "i_grid", dumpdata->nce1 + 1, &igridid);

  nc_def_dim(cdfid, "k_centre", dumpdata->nz, &kcentreid);
  nc_def_dim(cdfid, "j_centre", dumpdata->nce2, &jcentreid);
  nc_def_dim(cdfid, "i_centre", dumpdata->nce1, &icentreid);

  nc_def_dim(cdfid, "j_left", dumpdata->nce2, &jleftid);
  nc_def_dim(cdfid, "i_left", dumpdata->nce1 + 1, &ileftid);

  nc_def_dim(cdfid, "j_back", dumpdata->nce2 + 1, &jbackid);
  nc_def_dim(cdfid, "i_back", dumpdata->nce1, &ibackid);

  if (dumpdata->sednz > 0) {
    nc_def_dim(cdfid, "k_grid_sed", dumpdata->sednz + 1, &kgridid_sed);
    nc_def_dim(cdfid, "k_centre_sed", dumpdata->sednz, &kcentreid_sed);
  }

  /* time independent variables */
  dims[0] = kgridid;
  nc_def_var(cdfid, "z_grid", fptype, 1, dims, &vid);

  dims[0] = kcentreid;
  nc_def_var(cdfid, "z_centre", fptype, 1, dims, &vid);

  dims[0] = jgridid;
  dims[1] = igridid;
  nc_def_var(cdfid, "x_grid", NC_DOUBLE, 2, dims, &vid);
  nc_def_var(cdfid, "y_grid", NC_DOUBLE, 2, dims, &vid);

  dims[0] = jcentreid;
  dims[1] = icentreid;
  nc_def_var(cdfid, "x_centre", NC_DOUBLE, 2, dims, &vid);
  nc_def_var(cdfid, "y_centre", NC_DOUBLE, 2, dims, &vid);

  dims[0] = jleftid;
  dims[1] = ileftid;
  nc_def_var(cdfid, "x_left", NC_DOUBLE, 2, dims, &vid);
  nc_def_var(cdfid, "y_left", NC_DOUBLE, 2, dims, &vid);

  dims[0] = jbackid;
  dims[1] = ibackid;
  nc_def_var(cdfid, "x_back", NC_DOUBLE, 2, dims, &vid);
  nc_def_var(cdfid, "y_back", NC_DOUBLE, 2, dims, &vid);

  dims[0] = jcentreid;
  dims[1] = icentreid;
  nc_def_var(cdfid, "botz", NC_DOUBLE, 2, dims, &vid);

  dims[0] = jleftid;
  dims[1] = ileftid;
  nc_def_var(cdfid, "h1au1", NC_DOUBLE, 2, dims, &vid);

  dims[0] = jbackid;
  dims[1] = ibackid;
  nc_def_var(cdfid, "h1au2", NC_DOUBLE, 2, dims, &vid);

  dims[0] = jcentreid;
  dims[1] = icentreid;
  nc_def_var(cdfid, "h1acell", NC_DOUBLE, 2, dims, &vid);

  dims[0] = jgridid;
  dims[1] = igridid;
  nc_def_var(cdfid, "h1agrid", NC_DOUBLE, 2, dims, &vid);

  dims[0] = jleftid;
  dims[1] = ileftid;
  nc_def_var(cdfid, "h2au1", NC_DOUBLE, 2, dims, &vid);

  dims[0] = jbackid;
  dims[1] = ibackid;
  nc_def_var(cdfid, "h2au2", NC_DOUBLE, 2, dims, &vid);

  dims[0] = jcentreid;
  dims[1] = icentreid;
  nc_def_var(cdfid, "h2acell", NC_DOUBLE, 2, dims, &vid);

  dims[0] = jgridid;
  dims[1] = igridid;
  nc_def_var(cdfid, "h2agrid", NC_DOUBLE, 2, dims, &vid);

  dims[0] = jleftid;
  dims[1] = ileftid;
  nc_def_var(cdfid, "thetau1", NC_DOUBLE, 2, dims, &vid);

  dims[0] = jbackid;
  dims[1] = ibackid;
  nc_def_var(cdfid, "thetau2", NC_DOUBLE, 2, dims, &vid);

  dims[0] = jcentreid;
  dims[1] = icentreid;
  nc_def_var(cdfid, "coriolis", fptype, 2, dims, &vid);

  dims[0] = icentreid;
  nc_def_var(cdfid, "crci", NC_SHORT, 1, dims, &vid);
  nc_def_var(cdfid, "clci", NC_SHORT, 1, dims, &vid);

  dims[0] = igridid;
  nc_def_var(cdfid, "crfi", NC_SHORT, 1, dims, &vid);
  nc_def_var(cdfid, "clfi", NC_SHORT, 1, dims, &vid);

  dims[0] = icentreid;
  nc_def_var(cdfid, "frci", NC_SHORT, 1, dims, &vid);
  nc_def_var(cdfid, "flci", NC_SHORT, 1, dims, &vid);

  dims[0] = igridid;
  nc_def_var(cdfid, "frfi", NC_SHORT, 1, dims, &vid);
  nc_def_var(cdfid, "flfi", NC_SHORT, 1, dims, &vid);

  dims[0] = jcentreid;
  nc_def_var(cdfid, "cfcj", NC_SHORT, 1, dims, &vid);
  nc_def_var(cdfid, "cbcj", NC_SHORT, 1, dims, &vid);

  dims[0] = jgridid;
  nc_def_var(cdfid, "cffj", NC_SHORT, 1, dims, &vid);
  nc_def_var(cdfid, "cbfj", NC_SHORT, 1, dims, &vid);

  dims[0] = jcentreid;
  nc_def_var(cdfid, "ffcj", NC_SHORT, 1, dims, &vid);
  nc_def_var(cdfid, "fbcj", NC_SHORT, 1, dims, &vid);

  dims[0] = jgridid;
  nc_def_var(cdfid, "fffj", NC_SHORT, 1, dims, &vid);
  nc_def_var(cdfid, "fbfj", NC_SHORT, 1, dims, &vid);

  /* time dependent variables */
  dims[0] = recdimid;
  tid = nc_def_var(cdfid, "t", NC_DOUBLE, 1, dims, &vid);

  if (dumpdata->sednz > 0) {
    dims[0] = recdimid;
    dims[1] = kgridid_sed;
    dims[2] = jcentreid;
    dims[3] = icentreid;
    nc_def_var(cdfid, "z_grid_sed", fptype, 4, dims, &vid);

    dims[1] = kcentreid_sed;
    nc_def_var(cdfid, "z_centre_sed", fptype, 4, dims, &vid);
  }

  if (dumpdata->large) {
    dims[1] = jleftid;
    dims[2] = ileftid;
    nc_def_var(cdfid, "u1av", fptype, 3, dims, &vid);
    
    dims[1] = jbackid;
    dims[2] = ibackid;
    nc_def_var(cdfid, "u2av", fptype, 3, dims, &vid);
    
    dims[1] = jcentreid;
    dims[2] = icentreid;
    nc_def_var(cdfid, "wtop", fptype, 3, dims, &vid);
    nc_def_var(cdfid, "topz", fptype, 3, dims, &vid);
  }
  dims[1] = jcentreid;
  dims[2] = icentreid;
  nc_def_var(cdfid, "eta", fptype, 3, dims, &vid);

  /* 2D tracers */
  if (dumpdata->large) {
    for (n = 0; n < dumpdata->ntrS; n++)
      nc_def_var(cdfid, dumpdata->trinfo_2d[n].name, fptype, 3, dims, &vid);
  }

  dims[1] = jleftid;
  dims[2] = ileftid;
  nc_def_var(cdfid, "wind1", fptype, 3, dims, &vid);

  dims[1] = jbackid;
  dims[2] = ibackid;
  nc_def_var(cdfid, "wind2", fptype, 3, dims, &vid);

  dims[1] = jcentreid;
  dims[2] = icentreid;
  nc_def_var(cdfid, "patm", fptype, 3, dims, &vid);

  if (dumpdata->large) {
    dims[1] = kcentreid;
    dims[2] = jleftid;
    dims[3] = ileftid;
    nc_def_var(cdfid, "u1", fptype, 4, dims, &vid);

    dims[1] = kcentreid;
    dims[2] = jbackid;
    dims[3] = ibackid;
    nc_def_var(cdfid, "u2", fptype, 4, dims, &vid);

    dims[1] = kcentreid;
    dims[2] = jcentreid;
    dims[3] = icentreid;
    nc_def_var(cdfid, "w", fptype, 4, dims, &vid);
  }

  /* list of tracers */
  dims[2] = jcentreid;
  dims[3] = icentreid;
  /* watercolumn */
  dims[1] = kcentreid;
  for (n = 0; n < dumpdata->ntr; n++) {
    if (!dumpdata->large && 
	(strcmp(dumpdata->trinfo_3d[n].name, "salt") != 0) &&
	(strcmp(dumpdata->trinfo_3d[n].name, "temp") != 0)) continue;
      nc_def_var(cdfid, dumpdata->trinfo_3d[n].name, fptype, 4, dims, &vid);
  }

  /* sediment */
  if (dumpdata->large) {
    dims[1] = kcentreid_sed;
    for (n = 0; n < dumpdata->nsed; n++) {
      char name[MAXSTRLEN];
      strcpy(name, dumpdata->trinfo_sed[n].name);
      strcat(name, "_sed");
      nc_def_var(cdfid, name, fptype, 4, dims, &vid);
    }
  }

  dims[1] = kcentreid;
  dims[2] = jgridid;
  dims[3] = igridid;
  nc_def_var(cdfid, "flag", NC_INT, 4, dims, &vid);

  if (dumpdata->large) {
    dims[1] = kcentreid;
    dims[2] = jcentreid;
    dims[3] = icentreid;
    nc_def_var(cdfid, "dens", fptype, 4, dims, &vid);
    nc_def_var(cdfid, "dens_0", fptype, 4, dims, &vid);
    nc_def_var(cdfid, "Kz", fptype, 4, dims, &vid);
    nc_def_var(cdfid, "Vz", fptype, 4, dims, &vid);

    dims[1] = jleftid;
    dims[2] = ileftid;
    nc_def_var(cdfid, "u1bot", fptype, 3, dims, &vid);

    dims[1] = jbackid;
    dims[2] = ibackid;
    nc_def_var(cdfid, "u2bot", fptype, 3, dims, &vid);
  }
  write_dump_attributes(dumpdata, cdfid, fptype, 0, dumpdata->nce1,
                        dumpdata->nfe1, 0, dumpdata->nce2, dumpdata->nfe2, "",
			dumpdata->output_tunit);

  nc_sync(cdfid);
  nc_enddef(cdfid);
  return (cdfid);
}

void create_df(dump_data_t *dumpdata, char *name, int mode)
{
  dump_file_t *df = NULL;
  df = (dump_file_t *)malloc(sizeof(dump_file_t));
  memset(df, 0, sizeof(dump_file_t));
  strcpy(df->name, name);
  df->tout = dumpdata->t;
  df->tinc = 1.0;
  df->tstop = df->tout + df->tinc;
  df->bpv = 8;

  df->nvars = parseline(strdup("ALL"), df->vars, MAXNUMVARS);   

  df->landfill = locate_landfill_function("default");

  df->ilower = 0;
  df->jlower = 0;
  df->klower = 0;
  df->nce1 = dumpdata->nce1;
  df->nfe1 = dumpdata->nfe1;
  df->nce2 = dumpdata->nce2;
  df->nfe2 = dumpdata->nfe2;
  df->nz = dumpdata->nz;
  df->nz_sed = dumpdata->sednz;
  df->ns2 = dumpdata->ns2;
  df->ns3 = dumpdata->ns3;
  df->nface2 = dumpdata->nface2;
  df->nedge2 = dumpdata->nedge2;
  df->nvertex2 = dumpdata->nvertex2;
  df->nface3 = dumpdata->nface3;
  df->nedge3 = dumpdata->nedge3;
  df->nvertex3 = dumpdata->nvertex3;
  df->da_cycle = (NO_DA|DO_DA|NONE); // write in all cases
  df->finished = 0;
  df->compress = 0;
  df->long0_360 = 0;
  strcpy(df->tunit, dumpdata->output_tunit);

  if (mode & US_WUS) {
    strcpy(df->name, name);
    df_ugrid_create(dumpdata, df);
    df_ugrid_write(dumpdata, df, dumpdata->t);
    df_ugrid_close(dumpdata, df);
  } else {
    dumpdata_fill(geom, master, dumpdata);
    df_std_create(dumpdata, df);
    df_std_write(dumpdata, df, dumpdata->t);
    df_std_close(dumpdata, df);
  }

  free((dump_file_t *) df);
}
