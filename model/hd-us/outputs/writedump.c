/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/outputs/writedump.c
 *  
 *  Description:
 *  Routine to write netCDF model dump file
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: writedump.c 5873 2018-07-06 07:23:48Z riz008 $
 *
 */

#include <netcdf.h>
#include "hd.h"
#include "tracer.h"

void dump_write(dump_data_t *dumpdata, int cdfid, int ti)
{
  int n;
  size_t start[4];
  size_t count[4];
  double newt = dumpdata->t;

  if (DEBUG("dump"))
    dlog("dump", "Writing dump %ld at time %f days.\n", ti,
         dumpdata->t / 86400.0);

  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;
  /* time independent variables */
  count[0] = dumpdata->nz + 1;
  count[1] = 0;
  count[2] = 0;
  count[3] = 0;
  nc_put_vara_double(cdfid, ncw_var_id(cdfid, "z_grid"), start, count,
                     dumpdata->gridz);

  start[0] = 0;
  /* time independent variables */
  count[0] = dumpdata->nz;
  nc_put_vara_double(cdfid, ncw_var_id(cdfid, "z_centre"), start, count,
                     dumpdata->cellz);
  count[0] = dumpdata->nce2 + 1;
  count[1] = dumpdata->nce1 + 1;
  nc_put_vara_double(cdfid, ncw_var_id(cdfid, "x_grid"), start, count,
                     dumpdata->gridx[0]);
  nc_put_vara_double(cdfid, ncw_var_id(cdfid, "y_grid"), start, count,
                     dumpdata->gridy[0]);
  nc_put_vara_double(cdfid, ncw_var_id(cdfid, "h1agrid"), start, count,
                     dumpdata->h1agrid[0]);
  nc_put_vara_double(cdfid, ncw_var_id(cdfid, "h2agrid"), start, count,
                     dumpdata->h2agrid[0]);

  count[0] = dumpdata->nce2;
  count[1] = dumpdata->nce1;
  nc_put_vara_double(cdfid, ncw_var_id(cdfid, "x_centre"), start, count,
                     dumpdata->cellx[0]);
  nc_put_vara_double(cdfid, ncw_var_id(cdfid, "y_centre"), start, count,
                     dumpdata->celly[0]);
  nc_put_vara_double(cdfid, ncw_var_id(cdfid, "h1acell"), start, count,
                     dumpdata->h1acell[0]);
  nc_put_vara_double(cdfid, ncw_var_id(cdfid, "h2acell"), start, count,
                     dumpdata->h2acell[0]);
  nc_put_vara_double(cdfid, ncw_var_id(cdfid, "coriolis"), start, count,
                     dumpdata->coriolis[0]);
  nc_put_vara_double(cdfid, ncw_var_id(cdfid, "botz"), start, count,
                     dumpdata->botz[0]);

  count[0] = dumpdata->nce2;
  count[1] = dumpdata->nce1 + 1;
  nc_put_vara_double(cdfid, ncw_var_id(cdfid, "x_left"), start, count,
                     dumpdata->u1x[0]);
  nc_put_vara_double(cdfid, ncw_var_id(cdfid, "y_left"), start, count,
                     dumpdata->u1y[0]);
  nc_put_vara_double(cdfid, ncw_var_id(cdfid, "h1au1"), start, count,
                     dumpdata->h1au1[0]);
  nc_put_vara_double(cdfid, ncw_var_id(cdfid, "h2au1"), start, count,
                     dumpdata->h2au1[0]);
  nc_put_vara_double(cdfid, ncw_var_id(cdfid, "thetau1"), start, count,
                     dumpdata->thetau1[0]);

  count[0] = dumpdata->nce2 + 1;
  count[1] = dumpdata->nce1;
  nc_put_vara_double(cdfid, ncw_var_id(cdfid, "x_back"), start, count,
                     dumpdata->u2x[0]);
  nc_put_vara_double(cdfid, ncw_var_id(cdfid, "y_back"), start, count,
                     dumpdata->u2y[0]);
  nc_put_vara_double(cdfid, ncw_var_id(cdfid, "h1au2"), start, count,
                     dumpdata->h1au2[0]);
  nc_put_vara_double(cdfid, ncw_var_id(cdfid, "h2au2"), start, count,
                     dumpdata->h2au2[0]);
  nc_put_vara_double(cdfid, ncw_var_id(cdfid, "thetau2"), start, count,
                     dumpdata->thetau2[0]);

  count[0] = dumpdata->nce1;
  nc_put_vara_short(cdfid, ncw_var_id(cdfid, "crci"), start, count,
                    dumpdata->crci);
  nc_put_vara_short(cdfid, ncw_var_id(cdfid, "clci"), start, count,
                    dumpdata->clci);
  nc_put_vara_short(cdfid, ncw_var_id(cdfid, "frci"), start, count,
                    dumpdata->frci);
  nc_put_vara_short(cdfid, ncw_var_id(cdfid, "flci"), start, count,
                    dumpdata->flci);
  count[0] = dumpdata->nfe1;
  nc_put_vara_short(cdfid, ncw_var_id(cdfid, "crfi"), start, count,
                    dumpdata->crfi);
  nc_put_vara_short(cdfid, ncw_var_id(cdfid, "clfi"), start, count,
                    dumpdata->clfi);
  nc_put_vara_short(cdfid, ncw_var_id(cdfid, "frfi"), start, count,
                    dumpdata->frfi);
  nc_put_vara_short(cdfid, ncw_var_id(cdfid, "flfi"), start, count,
                    dumpdata->flfi);
  count[0] = dumpdata->nce2;
  nc_put_vara_short(cdfid, ncw_var_id(cdfid, "cfcj"), start, count,
                    dumpdata->cfcj);
  nc_put_vara_short(cdfid, ncw_var_id(cdfid, "cbcj"), start, count,
                    dumpdata->cbcj);
  nc_put_vara_short(cdfid, ncw_var_id(cdfid, "ffcj"), start, count,
                    dumpdata->ffcj);
  nc_put_vara_short(cdfid, ncw_var_id(cdfid, "fbcj"), start, count,
                    dumpdata->fbcj);
  count[0] = dumpdata->nfe2;
  nc_put_vara_short(cdfid, ncw_var_id(cdfid, "cffj"), start, count,
                    dumpdata->cffj);
  nc_put_vara_short(cdfid, ncw_var_id(cdfid, "cbfj"), start, count,
                    dumpdata->cbfj);
  nc_put_vara_short(cdfid, ncw_var_id(cdfid, "fffj"), start, count,
                    dumpdata->fffj);
  nc_put_vara_short(cdfid, ncw_var_id(cdfid, "fbfj"), start, count,
                    dumpdata->fbfj);

  /* time dependent variables */
  tm_change_time_units(dumpdata->timeunit, dumpdata->output_tunit, &newt,
                         1);
  start[0] = ti;
  count[0] = 1L;
  nc_put_vara_double(cdfid, ncw_var_id(cdfid, "t"), start, count, &newt);

  if (dumpdata->sednz > 0) {
    count[1] = dumpdata->sednz;
    count[2] = dumpdata->nce2;
    count[3] = dumpdata->nce1;
    nc_put_vara_double(cdfid, ncw_var_id(cdfid, "z_centre_sed"), start,
                       count, dumpdata->cellz_sed[0][0]);
    count[1] = dumpdata->sednz + 1;
    nc_put_vara_double(cdfid, ncw_var_id(cdfid, "z_grid_sed"), start,
                       count, dumpdata->gridz_sed[0][0]);
  }

  count[1] = dumpdata->nce2;
  count[2] = dumpdata->nfe1;
  if (dumpdata->large)
    nc_put_vara_double(cdfid, ncw_var_id(cdfid, "u1av"), start, count,
		       dumpdata->u1av[0]);
  nc_put_vara_double(cdfid, ncw_var_id(cdfid, "wind1"), start, count,
                     dumpdata->wind1[0]);
  count[1] = dumpdata->nfe2;
  count[2] = dumpdata->nce1;
  if (dumpdata->large)
    nc_put_vara_double(cdfid, ncw_var_id(cdfid, "u2av"), start, count,
		       dumpdata->u2av[0]);
  nc_put_vara_double(cdfid, ncw_var_id(cdfid, "wind2"), start, count,
                     dumpdata->wind2[0]);

  count[1] = dumpdata->nce2;
  count[2] = dumpdata->nce1;
  if (dumpdata->large) {
    nc_put_vara_double(cdfid, ncw_var_id(cdfid, "wtop"), start, count,
		       dumpdata->wtop[0]);
    nc_put_vara_double(cdfid, ncw_var_id(cdfid, "topz"), start, count,
		       dumpdata->topz[0]);
  }
  nc_put_vara_double(cdfid, ncw_var_id(cdfid, "eta"), start, count,
                     dumpdata->eta[0]);
  nc_put_vara_double(cdfid, ncw_var_id(cdfid, "patm"), start, count,
                     dumpdata->patm[0]);


  count[1] = dumpdata->nz;
  count[2] = dumpdata->nfe2;
  count[3] = dumpdata->nfe1;
  nc_put_vara_long(cdfid, ncw_var_id(cdfid, "flag"), start, count,
                   (long *)dumpdata->flag[0][0]);

  if (dumpdata->large) {
    count[2] = dumpdata->nce2;
    count[3] = dumpdata->nfe1;
    nc_put_vara_double(cdfid, ncw_var_id(cdfid, "u1"), start, count,
		       dumpdata->u1[0][0]);
    count[2] = dumpdata->nfe2;
    count[3] = dumpdata->nce1;
    nc_put_vara_double(cdfid, ncw_var_id(cdfid, "u2"), start, count,
		       dumpdata->u2[0][0]);
    count[2] = dumpdata->nce2;
    count[3] = dumpdata->nce1;
    nc_put_vara_double(cdfid, ncw_var_id(cdfid, "u"), start, count,
		       dumpdata->u[0][0]);
    count[2] = dumpdata->nce2;
    count[3] = dumpdata->nce1;
    nc_put_vara_double(cdfid, ncw_var_id(cdfid, "v"), start, count,
		       dumpdata->v[0][0]);
    count[2] = dumpdata->nce2;
    count[3] = dumpdata->nce1;
    nc_put_vara_double(cdfid, ncw_var_id(cdfid, "w"), start, count,
		       dumpdata->w[0][0]);
  }

  /* tracer list */
  /* 3D watercolumn tracers */
  count[2] = dumpdata->nce2;
  count[3] = dumpdata->nce1;
  for (n = 0; n < dumpdata->ntr; n++) {
    if (!dumpdata->large && 
	(strcmp(dumpdata->trinfo_3d[n].name, "salt") != 0) &&
	(strcmp(dumpdata->trinfo_3d[n].name, "temp") != 0)) continue;
    nc_put_vara_double(cdfid, ncw_var_id(cdfid, dumpdata->trinfo_3d[n].name),
		       start, count, (dumpdata->tr_wc[n][0][0]));
  }

  if (dumpdata->large) {
    /* 2D tracers */
    count[1] = dumpdata->nce2;
    count[2] = dumpdata->nce1;
    for (n = 0; n < dumpdata->ntrS; n++) {
      nc_put_vara_double(cdfid, ncw_var_id(cdfid, dumpdata->trinfo_2d[n].name),
			 start, count, (dumpdata->tr_wcS[n][0]));
    }

    /* sediments */
    count[1] = dumpdata->sednz;
    count[2] = dumpdata->nce2;
    count[3] = dumpdata->nce1;
    if (dumpdata->sednz > 0)
      for (n = 0; n < dumpdata->nsed; n++) {
	char name[MAXSTRLEN];
	strcpy(name, dumpdata->trinfo_sed[n].name);
	strcat(name, "_sed");
	nc_put_vara_double(cdfid, ncw_var_id(cdfid, name), start, count,
			   (dumpdata->tr_sed[n][0][0]));
      }

    /* diagnostic variables temporarily in dump file */
    count[1] = dumpdata->nz;
    nc_put_vara_double(cdfid, ncw_var_id(cdfid, "dens"), start, count,
			 dumpdata->dens[0][0]);
  }

  /* do a flush so that if the model stops abnormally we can examine the
     dump file */
  nc_sync(cdfid);
}
