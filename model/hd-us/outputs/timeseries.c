/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/outputs/timeseries.c
 *  
 *  Description:
 *  Routines to deal with time series output
 *  from the meco model. These are output data
 *  at a small number of points which are written
 *  to a file at regular intervals (much more
 *  frequently than whole model dumps are
 *  written)
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: timeseries.c 6446 2020-01-12 22:48:46Z riz008 $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hd.h"
#include "tracer.h"


// local functions
static int depth2c(master_t *master, ts_point_t ts, double depth, 
		   int mode, int cs);
void read_ts_data_init(FILE *fp, int n, master_t *master, ts_point_t *ts);
int stencil_size(geometry_t *geom, int kernal, int c);
static int d2c(master_t *master, double depth, int mode, int cs);
int read_ts_data(master_t *master, ts_point_t ts, double t, int c);
double average_glider_data(master_t *master, ts_point_t *tslist, double t, int c, int tn, 
			   int *cu, int *cd, double *tu, double *td, double *nvals);

/*
 * Logs a single timepoint to the timeseries file
 */
static void ts_add_point(double t, int n, int mode)
{
  master_t *master = tslist[n].master;
  geometry_t *geom = master->geom;
  int i, j, m, c, cs, e, es, t_index;
  double sinth, costh;
  int e1, e2;
  double u1val, u2val;
  double u1av, u2av;
  double xcpt, ycpt;
  double xcpta, ycpta;
  timeseries_t *loc_ts = &tslist[n].ts;
 
  if (!mode) {
    tm_change_time_units(master->timeunit, tslist[n].tsunits, &t, 1);
    fprintf(tslist[n].fp, "%f NaN NaN NaN NaN NaN", t);
  
    if (loc_ts->t != NULL) fprintf(tslist[n].fp, " NaN NaN NaN 0 0 0");
  
    for (m = 0; m < tslist[n].nvars; m++) {
      fprintf(tslist[n].fp, " NaN");
    }
    if (master->ptconc) fprintf(tslist[n].fp, " NaN");
    if (master->light)  fprintf(tslist[n].fp, " NaN");  
    if (tslist[n].ndata) {
      for (m = 0; m < tslist[n].dnvars; m++) {
	fprintf(tslist[n].fp, " NaN NaN");
	if (tslist[n].minv != NULL && tslist[n].maxv != NULL)
	  fprintf(tslist[n].fp, " NaN NaN");
	if (tslist[n].nvals != NULL && tslist[n].maxv != NULL)
	  fprintf(tslist[n].fp, " NaN");
      }
    }
    fprintf(tslist[n].fp, "\n");
    return;
  }

  c = depth2c(master, tslist[n], tslist[n].z, tslist[n].v_offset, tslist[n].c);
  cs = geom->m2d[c];
  if (loc_ts->t != NULL) {
    // Need to recalculate these trigs for the new position
    tslist[n].sinth = 0.0;
    tslist[n].costh = 0.0;
    for (m = 1; m <= geom->npe[cs]; m++) {
      e = geom->c2e[m][c];
      es = geom->m2de[e];
      tslist[n].sinth += sin(geom->thetau1[es]);
      tslist[n].costh += cos(geom->thetau1[es]);
    }
    tslist[n].sinth /= (double)geom->npe[cs];
    tslist[n].costh /= (double)geom->npe[cs];
  }      
  sinth = tslist[n].sinth;
  costh = tslist[n].costh;

  /*
  e1 = geom->c2e[6][c];
  e2 = geom->c2e[3][c];
  u1val = (master->u1[e1] + master->u1[e2]) / 2.0;
  u1av =
    (master->u1av[geom->m2de[e1]] + master->u1av[geom->m2de[e2]]) / 2.0;
  e1 = geom->c2e[2][c];
  e2 = geom->c2e[4][c];
  u2val = (master->u1[e1] + master->u1[e2]) / 2.0;
  u2av =
    (master->u1av[geom->m2de[e1]] + master->u1av[geom->m2de[e2]]) / 2.0;

  xcpt = u1val * costh - u2val * sinth;
  ycpt = u1val * sinth + u2val * costh;
  xcpta = u1av * costh - u2av * sinth;
  ycpta = u1av * sinth + u2av * costh;
  */

  xcpt = master->u[c];
  ycpt = master->v[c];
  xcpta = master->uav[cs];
  ycpta = master->vav[cs];

  /* Read in data */
  if (tslist[n].ndata) {
    if (read_ts_data(master, tslist[n], t, c))
      return;
  }

  // Make sure we log in proper time units
  tm_change_time_units(master->timeunit, tslist[n].tsunits, &t, 1);

  fprintf(tslist[n].fp, "%f %f %f %f %f %f", t, master->eta[cs],
	  xcpta, ycpta, xcpt, ycpt);
  
  if (loc_ts->t != NULL) {
    // write out the x,y,z,i,j & k
    fprintf(tslist[n].fp, " %f %f %f %d %d %d",
	    tslist[n].x, tslist[n].y, tslist[n].z,
	    tslist[n].i, tslist[n].j, geom->s2k[c]);
  }
  
  for (m = 0; m < tslist[n].nvars; m++) {
    int tn = tslist[n].vars[m];
    if(tslist[n].var_type[m] == WATER)
      fprintf(tslist[n].fp, " %f", master->tr_wc[tn][c]);
    else
      fprintf(tslist[n].fp, " %f", master->tr_wcS[tn][cs]);
  }
  if (master->ptconc)
    fprintf(tslist[n].fp, " %f", master->ptconc[c]);

  if (master->light)
    fprintf(tslist[n].fp, " %f", master->light[cs]);
  
  if (tslist[n].ndata) {
    for (m = 0; m < tslist[n].dnvars; m++) {
	fprintf(tslist[n].fp, " %f %f", tslist[n].val[m], tslist[n].obs[m]);
      if (tslist[n].minv != NULL && tslist[n].maxv != NULL)
	fprintf(tslist[n].fp, " %f %f", tslist[n].minv[m], tslist[n].maxv[m]);
      if (tslist[n].nvals != NULL && tslist[n].maxv != NULL)
	  fprintf(tslist[n].fp, " %f", tslist[n].nvals[m]);
    }
  }
  fprintf(tslist[n].fp, "\n");
}

/* 
 * Time scheduler_t functions for outputing spot point timeseries.
 */
static double ts_event(sched_event_t *event, double t)
{
  int n;
  timeseries_t *loc_ts = NULL;
  double tsout = schedule->stop_time;
  int fill_all = 0;
  
  /*
   * If using locations from file, we're forced to fill the whole
   * domain as we don't know the location ahead of time
   */
  for (n = 0; n < nts; ++n) {
    if (tslist[n].ts.t != NULL) {
      fill_all = 1;
      break;
    }
  }

  if (fill_all) {
    /*
     * Guard against dereferencing the function pointer before its had
     * a chance to be initialised, as at registration time
     */
    /*master->master_fill(master, window, windat, wincon);*/
    master_fill_ts(master, window, windat, wincon);
  } else {
    for (n = 0; n < nts; ++n) {
      if (tslist[n].metric & TS_GLIDER)
	master_fill_glider(master, window, windat, wincon, &tslist[n], t);
    }
    master_fill_ts(master, window, windat, wincon);
  }

  /* 
   * Write the next record for each timeseries file that has reached its
   * output time
   */
  for (n = 0; n < nts; ++n) {
    /* Skip if the DA flag doesn't match for this ts file */
    if (master->da && !(tslist[n].da_cycle & master->da)) continue;

    loc_ts = &tslist[n].ts;
    if (t >= (tslist[n].tsout - DT_EPS)) {
      // Keep reading from file, if necessary
      if (loc_ts->t != NULL) {
	double ts_t = t;

	// Make time consistent. i.e. convert time from input file to
	//                                                  model time
	tm_change_time_units(master->timeunit, loc_ts->t_units, &ts_t, 1);
	
	/* 
	 * See if we're in the range of the file
	 */
	if (ts_has_time(loc_ts, ts_t)) {
	  // Get the new i,j,k

	  tslist[n].x = ts_eval(loc_ts, tslist[n].varids[0], ts_t);
	  tslist[n].y = ts_eval(loc_ts, tslist[n].varids[1], ts_t);
	  tslist[n].z = ts_eval(loc_ts, tslist[n].varids[2], ts_t);
	  
	  if (loc_ts->df->variables[tslist[n].varids[2]].z_is_depth)
	    tslist[n].z = -tslist[n].z;
	  
	  /*
	   * Skip over this point if not in the domain
	   *
	   * We could issue a warning here but that can fill up the runlog
	   * file really quickly
	   */
	  if ((tslist[n].c = hd_grid_xytoij(master, tslist[n].x,  tslist[n].y, 
					    &tslist[n].i, &tslist[n].j)))
	    ts_add_point(t, n, 1);
	} else
	  ts_add_point(t, n, 0);
      } else {
	if (tslist[n].c > 0)
	  ts_add_point(t, n, 1); // time should already be in master units
      }
      tslist[n].tsout += tslist[n].tsdt;
      fflush(tslist[n].fp);
    }
  }

  /* Run through the timeseries list and find the * the one that next
     fires. */
  for (n = 0; n < nts; ++n) {
    if (master->da && !(tslist[n].da_cycle & master->da)) continue;
    if ((tslist[n].c > 0) && (tslist[n].tsout < tsout))
      tsout = tslist[n].tsout;
  }

  return tsout;
}

/*-------------------------------------------------------------------*/
/* Create the time series file for each element in the tslist, and   */
/* and dump the first time.                                          */
/*-------------------------------------------------------------------*/
int ts_init(sched_event_t *event)
{
  int i, j, c, cs, ncols;
  int coli; // column index

  /* Create the time series file for each element in the tslist */
  for (i = 0; i < nts; ++i) {
    timeseries_t *ots = NULL;
    master_t *master = tslist[i].master;
    geometry_t *geom = master->geom;

    /* If restarting an existing run, then read the old timeseries file
     * into memory */
    if (forced_restart || nrt_restart) {
      FILE *fp = fopen(tslist[i].pname, "rb");
      if (fp != NULL) {
	fclose(fp);
	ots = (timeseries_t *)malloc(sizeof(timeseries_t));
	memset(ots, 0, sizeof(timeseries_t));
	ts_read(tslist[i].pname, ots);
      }
    }
    
    /* Count the number of columns; */
    ncols = 6;  /* Time, ETA, UAV, VAV, U, V */
    coli  = 1;  /* Reset */

    /* Add spatial dimension, if called for */
    if (tslist[i].ts.t != NULL)
      ncols += (2*3); // add the x,y,z,i,j,k
    
    /* Add variables */
    ncols += tslist[i].nvars;
    
    /* Add particle concentration */
    if (master->ptconc)
      ++ncols;

    /* Add light */
    if (master->light)
      ++ncols;


    /* Add fuzzy verification comparisons */
    if (tslist[i].ndata) {
      ncols += (2 * tslist[i].dnvars);
    
      if (tslist[i].minv != NULL && tslist[i].maxv != NULL) 
	ncols += (2 * tslist[i].dnvars);

      if (tslist[i].nvals != NULL)
	ncols += tslist[i].dnvars;
    }
    
    /* Create the time series and write the header */
    if ((tslist[i].fp = fopen(tslist[i].pname, "w")) == NULL)
      hd_quit_and_dump("ts_init: Can't create timeseries file '%s'.\n",
		       tslist[i].pname);

    /* Write the comment header */
    fprintf(tslist[i].fp, "# Time series data from model run\n");
    fprintf(tslist[i].fp, "#\n");
    
    if (tslist[i].ts.t == NULL) {
      fprintf(tslist[i].fp,
              "# Points: x y z i j k cellx celly cellz botz\n");
      c = depth2c(master, tslist[i], tslist[i].z, tslist[i].v_offset,
		  tslist[i].c);
      if (tslist[i].c > 0)
	cs = geom->m2d[c];
      else
	c = cs = 0;
      fprintf(tslist[i].fp,
              "# %s %.8g %.8g %.8g %d %d %d %.8g %.8g %.8g %.8g\n",
              tslist[i].pname, tslist[i].x, tslist[i].y, tslist[i].z,
              tslist[i].i, tslist[i].j, geom->s2k[c], geom->cellx[cs],
              geom->celly[cs], geom->cellz[c] * master->Ds[cs],
	      geom->botz[cs] * master->Ds[cs]
	      );
    } else {
      fprintf(tslist[i].fp,
              "# Position data read from file '%s'\n", tslist[i].ts.name);
    }
    fprintf(tslist[i].fp, "#\n");
    
    /* Write out the number of columns */
    fprintf(tslist[i].fp, "## COLUMNS %d\n##\n", ncols);
    
    /* Write non-variable headers */
    fprintf(tslist[i].fp, "## COLUMN1.name time\n\
## COLUMN1.long_name Time\n\
## COLUMN1.units %s\n\
## COLUMN1.coordinate_type TIME\n\
## COLUMN1.fill_value 0\n\
## COLUMN1.missing_value -999999999\n\
##\n", tslist[i].tsunits);
      fprintf(tslist[i].fp, "## COLUMN2.name eta\n\
## COLUMN2.long_name Surface Elevation\n\
## COLUMN2.units metres\n\
## COLUMN2.fill_value 0\n\
## COLUMN2.missing_value -999\n\
##\n");
      fprintf(tslist[i].fp, "## COLUMN3.name uav\n\
## COLUMN3.long_name Depth averaged east current\n\
## COLUMN3.units ms-1\n\
## COLUMN3.fill_value 0\n\
## COLUMN3.missing_value -999\n\
##\n");
      fprintf(tslist[i].fp, "## COLUMN4.name vav\n\
## COLUMN4.long_name Depth averaged north current\n\
## COLUMN4.units ms-1\n\
## COLUMN4.fill_value 0\n\
## COLUMN4.missing_value -999\n\
##\n");
      fprintf(tslist[i].fp, "## COLUMN5.name u\n\
## COLUMN5.long_name East Current\n\
## COLUMN5.units m s-1\n\
## COLUMN5.fill_value 0\n\
## COLUMN5.missing_value -999\n\
##\n");
      fprintf(tslist[i].fp, "## COLUMN6.name v\n\
## COLUMN6.long_name North Current\n\
## COLUMN6.units m s-1\n\
## COLUMN6.fill_value 0\n\
## COLUMN6.missing_value -999\n\
##\n");
    coli += 6;
    if (tslist[i].ts.t != NULL) {
      fprintf(tslist[i].fp, "## COLUMN7.name x\n\
## COLUMN7.long_name Latitude\n\
## COLUMN7.units degrees\n\
## COLUMN7.fill_value 0\n\
## COLUMN7.missing_value -999\n\
##\n");
      fprintf(tslist[i].fp, "## COLUMN8.name y\n\
## COLUMN8.long_name Longitutde\n\
## COLUMN8.units degrees\n\
## COLUMN8.fill_value 0\n\
## COLUMN8.missing_value -999\n\
##\n");
      fprintf(tslist[i].fp, "## COLUMN9.name z\n\
## COLUMN9.long_name depth\n\
## COLUMN9.units m\n\
## COLUMN9.fill_value 0\n\
## COLUMN9.missing_value -999\n\
##\n");
      fprintf(tslist[i].fp, "## COLUMN10.name i\n\
## COLUMN10.long_name i\n\
## COLUMN10.units index\n\
## COLUMN10.fill_value 0\n\
## COLUMN10.missing_value -999\n\
##\n");
      fprintf(tslist[i].fp, "## COLUMN11.name j\n\
## COLUMN11.long_name i\n\
## COLUMN11.units index\n\
## COLUMN11.fill_value 0\n\
## COLUMN11.missing_value -999\n\
##\n");
      fprintf(tslist[i].fp, "## COLUMN12.name k\n\
## COLUMN12.long_name i\n\
## COLUMN12.units index\n\
## COLUMN12.fill_value 0\n\
## COLUMN12.missing_value -999\n\
##\n");
      coli += 6;
    }
    
    /* Write the tracer header */
    for (j = 0; j < tslist[i].nvars; ++j) {
      int tn = tslist[i].vars[j];
      tracer_info_t *trinfo;
      
      trinfo = (tslist[i].var_type[j] == WATER) ? master->trinfo_3d :
	master->trinfo_2d;
      fprintf(tslist[i].fp, "## COLUMN%1d.name %s\n", coli,trinfo[tn].name);
      fprintf(tslist[i].fp, "## COLUMN%1d.long_name %s\n", coli,
	      trinfo[tn].long_name);
      fprintf(tslist[i].fp, "## COLUMN%1d.units %s\n", coli,
	      trinfo[tn].units);
      fprintf(tslist[i].fp, "## COLUMN%1d.fill_value_wc %g\n", coli,
	      trinfo[tn].fill_value_wc);
      fprintf(tslist[i].fp, "## COLUMN%1d.missing_value -99999999\n", coli);
      fprintf(tslist[i].fp, "##\n");
      coli++;
    }
    /* Header for particle concentrations */
    if (master->ptconc) {
      fprintf(tslist[i].fp, "## COLUMN%1d.name ptconc\n", coli);
      fprintf(tslist[i].fp,
	      "## COLUMN%1d.long_name Particle_Concentration\n", coli);
      fprintf(tslist[i].fp, "## COLUMN%1d.units 1\n", coli);
      fprintf(tslist[i].fp, "## COLUMN%1d.fill_value 0\n", coli);
      fprintf(tslist[i].fp, "## COLUMN%1d.missing_value -99999999\n",
	      coli);
      fprintf(tslist[i].fp, "##\n");
      coli++;
    }
    /* Header for light */
    if (master->light) {
      fprintf(tslist[i].fp, "## COLUMN%1d.name light\n", coli);
      fprintf(tslist[i].fp,
	      "## COLUMN%1d.long_name Light\n", coli);
      fprintf(tslist[i].fp, "## COLUMN%1d.units Wm-2\n", coli);
      fprintf(tslist[i].fp, "## COLUMN%1d.fill_value 0\n", coli);
      fprintf(tslist[i].fp, "## COLUMN%1d.missing_value -99999999\n",
	      coli);
      fprintf(tslist[i].fp, "##\n");
      coli++;
    }

    /* Header for model-data fuzzy verification */
    if (tslist[i].ndata) {
      for (j = 0; j < tslist[i].dnvars; ++j) {
	char name[MAXSTRLEN], units[MAXSTRLEN];

	if (strcmp(tslist[i].dvars[j], "PRES") == 0) {
	  strcpy(name, "N2");
	  strcpy(units, "s-2");
	} else {
	  strcpy(name, tslist[i].dvars[j]);
	  strcpy(units, tslist[i].dunits[j]);
	}
	fprintf(tslist[i].fp, "## COLUMN%1d.name %s_comp\n", coli, name);
	fprintf(tslist[i].fp, "## COLUMN%1d.long_name %s_comp\n", coli, name);
	fprintf(tslist[i].fp, "## COLUMN%1d.units %s\n", coli, units);
	fprintf(tslist[i].fp, "## COLUMN%1d.fill_value 0\n", coli);
	fprintf(tslist[i].fp, "## COLUMN%1d.missing_value -99999999\n",
		coli);
	fprintf(tslist[i].fp, "##\n");
	coli++;
	fprintf(tslist[i].fp, "## COLUMN%1d.name %s_obs\n", coli, name);
	fprintf(tslist[i].fp, "## COLUMN%1d.long_name %s_obs\n", coli, name);
	fprintf(tslist[i].fp, "## COLUMN%1d.units %s\n", coli, units);
	fprintf(tslist[i].fp, "## COLUMN%1d.fill_value 0\n", coli);
	fprintf(tslist[i].fp, "## COLUMN%1d.missing_value -99999999\n",
		coli);
	fprintf(tslist[i].fp, "##\n");
	coli++;
	if (tslist[i].minv != NULL && tslist[i].maxv != NULL) {
	  fprintf(tslist[i].fp, "## COLUMN%1d.name %s_minv\n", coli, name);
	  fprintf(tslist[i].fp, "## COLUMN%1d.long_name %s_minv\n", coli, name);
	  fprintf(tslist[i].fp, "## COLUMN%1d.units %s\n", coli, units);
	  fprintf(tslist[i].fp, "## COLUMN%1d.fill_value 0\n", coli);
	  fprintf(tslist[i].fp, "## COLUMN%1d.missing_value -99999999\n",
		  coli);
	  fprintf(tslist[i].fp, "##\n");
	  coli++;
	  fprintf(tslist[i].fp, "## COLUMN%1d.name %s_maxv\n", coli, name);
	  fprintf(tslist[i].fp, "## COLUMN%1d.long_name %s_maxv\n", coli, name);
	  fprintf(tslist[i].fp, "## COLUMN%1d.units %s\n", coli, units);
	  fprintf(tslist[i].fp, "## COLUMN%1d.fill_value 0\n", coli);
	  fprintf(tslist[i].fp, "## COLUMN%1d.missing_value -99999999\n",
		  coli);
	  fprintf(tslist[i].fp, "##\n");
	  coli++;
	}
	if (tslist[i].nvals != NULL) {
	  fprintf(tslist[i].fp, "## COLUMN%1d.name %s_nvals\n", coli, name);
	  fprintf(tslist[i].fp, "## COLUMN%1d.long_name %s_nvals\n", coli, name);
	  fprintf(tslist[i].fp, "## COLUMN%1d.units \n", coli);
	  fprintf(tslist[i].fp, "## COLUMN%1d.fill_value 0\n", coli);
	  fprintf(tslist[i].fp, "## COLUMN%1d.missing_value -99999999\n",
		  coli);
	  fprintf(tslist[i].fp, "##\n");
	  coli++;
	}
      }
    }

    /* If there was a file prior to restart, then output records upto
     * the current time. */
    if (ots != NULL) {
      int r=0;
      int v=0;
      for (r=0; r<ots->nt; ++r) {
	double newt = ots->t[r];
	tm_change_time_units(tslist[i].tsunits, master->timeunit, &newt, 1);
	if (newt >= master->t) break;
	
	for (v=0; v<ots->nv; ++v)
	  fprintf(tslist[i].fp, "%f ", ots->df->variables[v].data[r]);
	fprintf(tslist[i].fp, "\n");
      }
      
      free(ots);
      ots = NULL;
    }

    fflush(tslist[i].fp);
  }
  return 1;
}


/* Close all of the time series files.
 */
static void ts_cleanup(sched_event_t *event, double t)
{
  int i;

  for (i = 0; i < nts; ++i) {
    if (tslist[i].c > 0) {
      fclose(tslist[i].fp);
      tslist[i].i = -1;
      tslist[i].c = 0;
    }
  }
}



void timeseries_init(FILE * prmfd, master_t *master,geometry_t *geom,
		     geometry_t **window, dump_data_t *dumpdata)
{
  int n, m, wn;
  int i;
  int j;
  int c, cs = 0, e, es;
  char key[MAXSTRLEN];
  char buf[MAXSTRLEN];
  char fname[MAXSTRLEN];
  char *vars[MAXNUMVARS];       /* Tracer output variable names */
  FILE * tsfd = prmfd;
  FILE *locF = NULL;

  if (master->runmode & DUMP) return;
	
  /* Get number of time series points from parameter file */
  prm_set_errfn(hd_silent_warn);
  
  if(!prm_read_char(tsfd, "TSPOINTS", buf)) {
    if(prm_read_char(tsfd, "TSFILE", buf))
      tsfd = fopen(buf,"r");
    else {
      nts = 0;
      return;
    }
  }

  if ((prm_read_int(tsfd, "TSPOINTS", &nts)) <= 0 || (nts <= 0)) {
    nts = 0;
    return;                     /* nothing to do - no points */
  }
  if (master->tmode & SP_CHECK) {
    nts = 0;
    return;                     /* nothing to do - no points */
  }

  prm_set_errfn(hd_silent_warn);
  /* allocate storage for the list of points */
  if ((tslist = (ts_point_t *)malloc(sizeof(ts_point_t) * nts)) == NULL)
    hd_quit_and_dump
      ("timeseries_init: Not enough memory for point list\n");
  /* read point list */
  prm_set_errfn(hd_quit);

  for (n = 0; n < nts; n++) {
    timeseries_t *loc_ts = &tslist[n].ts;

    sprintf(key, "TS%1d.name", n);
    prm_read_char(tsfd, key, tslist[n].pname);

    if (strlen(master->opath)) {
      sprintf(buf, "%s%s", master->opath, tslist[n].pname);
      strcpy(tslist[n].pname, buf);
    }

    sprintf(key, "TS%1d.location", n);
    prm_read_char(tsfd, key, buf);
    memset(loc_ts, 0, sizeof(timeseries_t));

    // See if location is a timeseries file
    loc_ts->t = NULL;
    if ( (locF = fopen(fv_get_filename(buf, fname), "r")) != NULL ) {
      fclose(locF);
      ts_read(fname, loc_ts);
    }

    if (loc_ts->t == NULL) {
      // Assume numeric
      if (sscanf(buf, "%lf %lf %lf",
		 &tslist[n].x, &tslist[n].y, &tslist[n].z) != 3)
	hd_quit_and_dump("timeseries_init: Can't read point list\n");
    } else {
      // Cache the variable id's
      char *varname;
      char buf2[MAXSTRLEN];
      // Make sure buf points to the the full name string
      varname = fv_get_varname(buf, "lon", buf2);
      if ((tslist[n].varids[0] = ts_get_index(loc_ts, varname)) < 0)
	hd_quit("ts_init: variable for 'lon' not found in '%s'\n", buf);

      varname = fv_get_varname(buf, "lat", buf2);
      if ((tslist[n].varids[1] = ts_get_index(loc_ts, varname)) < 0)
	hd_quit("ts_init: variable for 'lat' not found in '%s'\n", buf);

      varname = fv_get_varname(buf, "depth", buf2);
      if ((tslist[n].varids[2] = ts_get_index(loc_ts, varname)) < 0)
	hd_quit("ts_init: variable for 'depth' not found in '%s'\n", buf);
    }
    sprintf(key, "TS%1d.dt", n);

    prm_get_time_in_secs(tsfd, key, &tslist[n].tsdt);
    if (tslist[n].tsdt < 0.0)
      hd_quit("ts_init: Invalid TS step.\n");

    prm_set_errfn(hd_silent_warn);
    strcpy(tslist[n].tsunits, "default");

  /*---------------------------------------------------------*/
    /* Include variable vertical references for the */
    /* time series point.  */
    /* 0 = referenced to mean sea level */
    /* 1 = referenced to the bottom */
    /* 2 = referenced to the free surface */
    tslist[n].v_offset = 0;
    sprintf(key, "TS%1d.reference", n);
    if (prm_read_char(tsfd, key, buf)) {
      if (strcmp(buf, "msl") == 0 || strcmp(buf, "MSL") == 0)
        tslist[n].v_offset = 0;
      if (strcmp(buf, "bottom") == 0 || strcmp(buf, "BOTTOM") == 0)
        tslist[n].v_offset = 1;
      if (strcmp(buf, "surface") == 0 || strcmp(buf, "SURFACE") == 0)
        tslist[n].v_offset = 2;
    }

    tslist[n].type = SIMPLE;
    sprintf(key, "TS%1d.type", n);
    if (prm_read_char(tsfd, key, buf)) {
      if (strcmp(buf, "simple") == 0 || strcmp(buf, "SIMPLE") == 0)
        tslist[n].type = SIMPLE;
      if (strcmp(buf, "standard") == 0 || strcmp(buf, "STANDARD") == 0)
        tslist[n].type = STANDARD;
    }

    /*
     * DA flags
     * Please make sure any changes are in sync with dumpfile_init
     */
    tslist[n].da_cycle = NONE;
    if (master->params->da) {
      /* Make forecast mode the default */
      tslist[n].da_cycle = NO_DA;

      sprintf(key, "TS%d.da", n);
      if (prm_read_char(tsfd, key, buf)) {
	if (contains_token(buf, "forecast"))
	  tslist[n].da_cycle = NO_DA;
	else if (contains_token(buf, "reanalysis"))
	  tslist[n].da_cycle = DO_DA;
	else
	  hd_quit("ts_init: unknown da_cycle '%s' specified\n", buf);
      }
    }

    /* Tracer variables */
    sprintf(key, "TS%1d.vars", n);
    prm_skip_to_end_of_key(tsfd, key);
    if (fgets(buf, MAXSTRLEN, tsfd) == NULL) {
      tslist[n].nvars = master->ntr + master->ntrS;
      for(i = 0; i < master->ntr; i++) {
				tslist[n].vars[i] = i;
				tslist[n].var_type[i] = WATER;
      }
      for(i = 0; i < master->ntrS; i++) {
				tslist[n].vars[i + master->ntr] = i;
				tslist[n].var_type[i + master->ntr] = INTER;
      }
    }
    else {
      tslist[n].nvars = parseline(strdup(buf), vars, MAXNUMVARS);
      for(i = 0; i < tslist[n].nvars; i++)
	if((tslist[n].vars[i] = tracer_find_index(vars[i], master->ntr,
						  master->trinfo_3d)) >= 0)
	  tslist[n].var_type[i] = WATER;
	else if((tslist[n].vars[i] = tracer_find_index(vars[i], master->ntrS,
						       master->trinfo_2d)) >= 0)
	  tslist[n].var_type[i] = INTER;
	else {
	  hd_warn("timeseries_init: can't find tracer %s\n",vars[i]);
	  tslist[n].nvars--;
	  i--;
	}
    }

    /*
     * units field is optional
     */
    strcpy(tslist[n].tsunits, "default");
    sprintf(key, "TS%1d.units", n);
    if (prm_read_char(tsfd, key, buf)) {
      // copy over the time units
      strcpy(tslist[n].tsunits, buf);
    }

    prm_set_errfn(hd_quit);

    // Some defaults
    tslist[n].i = -1;           /* flag point as bad */
    tslist[n].c = 0;
    tslist[n].fp = NULL;
    tslist[n].master = master;
    tslist[n].tsout = schedule->start_time;

    /* SIGMA : Multiply botz by total depth - ?? FR */
    if (loc_ts->t == NULL) {

      if ((c = hd_grid_xytoij_r(master, tslist[n].x, tslist[n].y, &i, &j))) {
	cs = geom->m2d[c];

	if (cs == 0) {
	  hd_warn("timeseries_init: point %s ignored: not in water\n",
		  tslist[n].pname);
	  continue;
	}
	if (tslist[n].v_offset == 1) {
	  double depth = geom->botz[cs] * master->Ds[cs] + fabs(tslist[n].z);
	  c = depth2c(master, tslist[n], depth, 2, cs);
	} else
	  c = depth2c(master, tslist[n], tslist[n].z, 2, cs);
	
	if (c >= 0 && geom->m2d[c] <= geom->ewetS) {
	  tslist[n].c = cs;
	  tslist[n].i = i;
	  tslist[n].j = j;
	  tslist[n].k = c;
	  
	  /* calculate sin and cos of the angle between the e1 direction and
	     the x axis at this point */
	  tslist[n].sinth = 0.0;
	  tslist[n].costh = 0.0;
	  for (m = 1; m <= geom->npe[cs]; m++) {
	    e = geom->c2e[m][c];
	    es = geom->m2de[e];
	    tslist[n].sinth += sin(geom->thetau1[es]);
	    tslist[n].costh += cos(geom->thetau1[cs]);
	  }
	  tslist[n].sinth /= (double)geom->npe[cs];
	  tslist[n].costh /= (double)geom->npe[cs];
	}
      } else
	hd_warn("timeseries_init: point %s ignored: not in water\n",
		tslist[n].pname);
      if (tslist[n].i < 0)
	hd_warn("timeseries_init: point %s ignored: not in water\n",
		tslist[n].pname);
    } else {
      /*
       * Need this >= 0 for the case when we start the model earlier
       * than the first timestamp in the file
       */
      tslist[n].i = 0;
      tslist[n].c = 1;
    }

    /* Check if we should use the default timeunit */
    if (strcasecmp(tslist[n].tsunits, "default") == 0)
      strcpy(tslist[n].tsunits, master->output_tunit);

    /* Input data */
    read_ts_data_init(tsfd, n, master, &tslist[n]);    

  } // foreach n=nts

  if(tsfd != prmfd)
    fclose(tsfd);

  /* Register an interest for output with the time scheduler. */
  sched_register(schedule, "timeseries", ts_init,
                 ts_event, ts_cleanup, NULL, NULL, NULL);

  /*-----------------------------------------------------------------*/
  /* Save the locations of surface time series points. The whole     */
  /* column is transferred since the exact location in the water     */
  /* column of the time series point may vary according to the       */
  /* reference point.                                                */
  for (wn = 1; wn <= master->nwindows; wn++) {
    int cc = 0, m;
    for (n = 0; n < nts; ++n) {
      int *st = NULL, ssize;
      if (tslist[n].c <= 0) continue;
      c = tslist[n].c;
      if (c <=0 || c > geom->sgnum) continue;
      if (geom->fm[c].wn == wn) {
	ssize = stencil_size(geom, tslist[n].kernal, c);
	st = stencil(geom, c, &ssize, ST_SIZED, 0);
	for (m = 0; m < ssize; m++) {
	  c = st[m];
	  while (c != geom->zm1[c]) {
	    cc++;
	    c = geom->zm1[c];
	  }
	}
	i_free_1d(st);
      }
    }
    window[wn]->ns2m_ts = cc;
    if (cc) window[wn]->s2m_ts = i_alloc_1d(cc + 1);
    cc = 1;
    for (n = 0; n < nts; ++n) {
      int *st = NULL, ssize;
      if (tslist[n].c <= 0) continue;
      c = tslist[n].c;
      if (c <=0 || c > geom->sgnum) continue;
      if (geom->fm[c].wn == wn) {
	ssize = stencil_size(geom, tslist[n].kernal, c);
	st = stencil(geom, c, &ssize, ST_SIZED, 0);
	for (m = 0; m < ssize; m++) {
	  c = st[m];
	  while (c != geom->zm1[c]) {
	    window[wn]->s2m_ts[cc] = geom->fm[c].sc;
	    cc++;
	    c = geom->zm1[c];
	  }
	}
	i_free_1d(st);
      }
    }
    if (DEBUG("init_w"))
      dlog("init_w",
	   "  Window %d : %d slave to master timeseries transfer cells\n", wn,
	   window[wn]->ns2m_ts);
  }
}


void timeseries_end(void)
{
  sched_deregister(schedule, "timeseries");
  if (tslist)
    free(tslist);
  tslist = NULL;
}

/*-------------------------------------------------------------------*/
/* Returns the vertical cell location given a depth.                 */
/*-------------------------------------------------------------------*/
static int depth2c(master_t *master,   /* Master data structure */
		   ts_point_t ts,      /* Time series point data structure */
		   double depth,       /* Depth to convert to vertical index */
		   int mode,           /* Reference level code */
		   int cs              /* Surface coordinate */
  )
{
  geometry_t *geom = master->geom;
  int c = 0;
  if (mode == 0 || mode == 1)
    return (ts.k);              /* ts.k contains the sparse coordinate */
  else {
    c = cs;
    if (c > 0 && c <= geom->sgnum) {
      double dc = 0.0;
      while (c != geom->zm1[c]) {
        dc += master->dz[c] * master->Ds[cs];
        c = geom->zm1[c];
        if (dc >= fabs(depth))
          break;
      }
    }
    return (geom->zp1[c]);
  }
}

static int d2c(master_t *master,   /* Master data structure */
	       double depth,       /* Depth to convert to vertical index */
	       int mode,           /* Reference level code */
	       int cs              /* Surface coordinate */
  )
{
  geometry_t *geom = master->geom;
  int c = 0;

  c = cs;
  if (c > 0 && c <= geom->sgnum) {
    double dc = 0.0;
    while (c != geom->zm1[c]) {
      dc += master->dz[c] * master->Ds[cs];
      c = geom->zm1[c];
      if (dc >= fabs(depth))
	break;
    }
  }
  return (geom->zp1[c]);
}

/* END depth2c()                                                     */
/*-------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Timeseries file initialisation on the master                     */
/*------------------------------------------------------------------*/
void timeseries_init_m(FILE * prmfd, master_t *master,geometry_t *geom,
		       geometry_t **window, dump_data_t *dumpdata)
{
  int n, wn, m;
  int i;
  int j;
  int c, cs = 0, e, es;
  char key[MAXSTRLEN];
  char buf[MAXSTRLEN];
  char fname[MAXSTRLEN];
  char *vars[MAXNUMVARS];       /* Tracer output variable names */
  FILE * tsfd = prmfd;
  FILE *locF = NULL;

	
  /* Get number of time series points from parameter file */
  prm_set_errfn(hd_silent_warn);
  
  if(!prm_read_char(tsfd, "TSPOINTS", buf)) {
    if(prm_read_char(tsfd, "TSFILE", buf))
      tsfd = fopen(buf,"r");
    else {
      nts = 0;
      return;
    }
  }

  if ((prm_read_int(tsfd, "TSPOINTS", &nts)) <= 0 || (nts <= 0)) {
    nts = 0;
    return;                     /* nothing to do - no points */
  }

  prm_set_errfn(hd_silent_warn);
  /* allocate storage for the list of points */
  if ((tslist = (ts_point_t *)malloc(sizeof(ts_point_t) * nts)) == NULL)
    hd_quit_and_dump
      ("timeseries_init: Not enough memory for point list\n");
  /* read point list */
  prm_set_errfn(hd_quit);

  for (n = 0; n < nts; n++) {
    timeseries_t *loc_ts = &tslist[n].ts;

    sprintf(key, "TS%1d.name", n);
    prm_read_char(tsfd, key, tslist[n].pname);

    if (strlen(master->opath)) {
      sprintf(buf, "%s%s", master->opath, tslist[n].pname);
      strcpy(tslist[n].pname, buf);
    }

    sprintf(key, "TS%1d.location", n);
    prm_read_char(tsfd, key, buf);
    memset(loc_ts, 0, sizeof(timeseries_t));

    // See if location is a timeseries file
    loc_ts->t = NULL;
    if ( (locF = fopen(fv_get_filename(buf, fname), "r")) != NULL ) {
      fclose(locF);
      ts_read(fname, loc_ts);
    }

    if (loc_ts->t == NULL) {
      // Assume numeric

      if (sscanf(buf, "%lf %lf %lf",
		 &tslist[n].x, &tslist[n].y, &tslist[n].z) != 3)
	hd_quit_and_dump("timeseries_init: Can't read point list\n");

    } else {
      // Cache the variable id's
      char *varname;

      varname = fv_get_varname(fname, "lon", buf);
      if ((tslist[n].varids[0] = ts_get_index(loc_ts, varname)) < 0)
	hd_quit("ts_init: variable for 'lon' not found in '%s'\n", buf);

      varname = fv_get_varname(fname, "lat", buf);
      if ((tslist[n].varids[1] = ts_get_index(loc_ts, varname)) < 0)
	hd_quit("ts_init: variable for 'lat' not found in '%s'\n", buf);

      varname = fv_get_varname(fname, "depth", buf);
      if ((tslist[n].varids[2] = ts_get_index(loc_ts, varname)) < 0)
	hd_quit("ts_init: variable for 'depth' not found in '%s'\n", buf);
    }
    sprintf(key, "TS%1d.dt", n);

    prm_get_time_in_secs(tsfd, key, &tslist[n].tsdt);
    if (tslist[n].tsdt < 0.0)
      hd_quit("ts_init: Invalid TS step.\n");

    prm_set_errfn(hd_silent_warn);
    strcpy(tslist[n].tsunits, "default");

    /*---------------------------------------------------------*/
    /* Include variable vertical references for the */
    /* time series point.  */
    /* 0 = referenced to mean sea level */
    /* 1 = referenced to the bottom */
    /* 2 = referenced to the free surface */
    tslist[n].v_offset = 0;
    sprintf(key, "TS%1d.reference", n);
    if (prm_read_char(tsfd, key, buf)) {
      if (strcmp(buf, "msl") == 0 || strcmp(buf, "MSL") == 0)
        tslist[n].v_offset = 0;
      if (strcmp(buf, "bottom") == 0 || strcmp(buf, "BOTTOM") == 0)
        tslist[n].v_offset = 1;
      if (strcmp(buf, "surface") == 0 || strcmp(buf, "SURFACE") == 0)
        tslist[n].v_offset = 2;
    }

    tslist[n].type = SIMPLE;
    sprintf(key, "TS%1d.type", n);
    if (prm_read_char(tsfd, key, buf)) {
      if (strcmp(buf, "simple") == 0 || strcmp(buf, "SIMPLE") == 0)
        tslist[n].type = SIMPLE;
      if (strcmp(buf, "standard") == 0 || strcmp(buf, "STANDARD") == 0)
        tslist[n].type = STANDARD;
    }

    /*
     * DA flags
     * Please make sure any changes are in sync with dumpfile_init
     */
    tslist[n].da_cycle = NONE;
    if (master->params->da) {
      /* Make forecast mode the default */
      tslist[n].da_cycle = NO_DA;

      sprintf(key, "TS%d.da", n);
      if (prm_read_char(tsfd, key, buf)) {
	if (contains_token(buf, "forecast"))
	  tslist[n].da_cycle = NO_DA;
	else if (contains_token(buf, "reanalysis"))
	  tslist[n].da_cycle = DO_DA;
	else
	  hd_quit("ts_init: unknown da_cycle '%s' specified\n", buf);
      }
    }

    /* Tracer variables */
    sprintf(key, "TS%1d.vars", n);
    prm_skip_to_end_of_key(tsfd, key);
    if (fgets(buf, MAXSTRLEN, tsfd) == NULL) {
      tslist[n].nvars = master->ntr + master->ntrS;
      for(i = 0; i < master->ntr; i++) {
				tslist[n].vars[i] = i;
				tslist[n].var_type[i] = WATER;
      }
      for(i = 0; i < master->ntrS; i++) {
				tslist[n].vars[i + master->ntr] = i;
				tslist[n].var_type[i + master->ntr] = INTER;
      }
    }
    else {
      tslist[n].nvars = parseline(strdup(buf), vars, MAXNUMVARS);
      for(i = 0; i < tslist[n].nvars; i++)
	if((tslist[n].vars[i] = tracer_find_index(vars[i], master->ntr,
						  master->trinfo_3d)) >= 0)
	  tslist[n].var_type[i] = WATER;
	else if((tslist[n].vars[i] = tracer_find_index(vars[i], master->ntrS,
						       master->trinfo_2d)) >= 0)
	  tslist[n].var_type[i] = INTER;
	else {
	  hd_warn("timeseries_init: can't find tracer %s\n",vars[i]);
	  tslist[n].nvars--;
	  i--;
	}
    }
    
    /*
     * units field is optional
     */
    strcpy(tslist[n].tsunits, "default");
    sprintf(key, "TS%1d.units", n);
    if (prm_read_char(tsfd, key, buf)) {
      // copy over the time units
      strcpy(tslist[n].tsunits, buf);
    }

    prm_set_errfn(hd_quit);

    // Some defaults
    tslist[n].i = -1;           /* flag point as bad */
    tslist[n].c = 0;
    tslist[n].fp = NULL;
    tslist[n].master = master;
    tslist[n].tsout = schedule->start_time;
    
    /* SIGMA : Multiply botz by total depth - ?? FR */
    if (loc_ts->t == NULL) {
      if ((c = hd_grid_xytoij(master, tslist[n].x, tslist[n].y, &i, &j))) {
	cs = geom->m2d[c];

	if (cs == 0) {
	  hd_warn("timeseries_init: point %s ignored: not in water\n",
		  tslist[n].pname);
	  continue;
	}
	if (tslist[n].v_offset == 1) {
	  double depth = geom->botz[cs] * master->Ds[cs] + fabs(tslist[n].z);
	  c = depth2c(master, tslist[n], depth, 2, cs);
	} else
	  c = depth2c(master, tslist[n], tslist[n].z, 2, cs);
	
	if (c >= 0 && geom->m2d[c] <= geom->ewetS) {
	  tslist[n].c = cs;
	  tslist[n].i = i;
	  tslist[n].j = j;
	  tslist[n].k = c;
	  
	  /* calculate sin and cos of the angle between the e1 direction and
	     the x axis at this point */
	  tslist[n].sinth = 0.0;
	  tslist[n].costh = 0.0;
	  for (m = 1; m <= geom->npe[cs]; m++) {
	    e = geom->c2e[m][c];
	    es = geom->m2de[e];
	    tslist[n].sinth += sin(geom->thetau1[es]);
	    tslist[n].costh += cos(geom->thetau1[cs]);
	  }
	  tslist[n].sinth /= (double)geom->npe[cs];
	  tslist[n].costh /= (double)geom->npe[cs];
	}
      } else
	hd_warn("timeseries_init: point %s ignored: not in water\n",
		tslist[n].pname);
      if (tslist[n].i < 0)
	hd_warn("timeseries_init: point %s ignored: not in water\n",
		tslist[n].pname);
    } else {
      /*
       * Need this >= 0 for the case when we start the model earlier
       * than the first timestamp in the file
       */
      tslist[n].i = 0;
      tslist[n].c = 1;
    }

    /* Check if we should use the default timeunit */
    if (strcasecmp(tslist[n].tsunits, "default") == 0)
      strcpy(tslist[n].tsunits, master->output_tunit);
    
  } // foreach n=nts
  
  if(tsfd != prmfd)
    fclose(tsfd);
}

/* END timeseries_init_m()                                          */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Timeseries file initialisation on the windows                    */
/*------------------------------------------------------------------*/
void timeseries_init_w(master_t *master, geometry_t **window)
		       
{
  geometry_t *geom = master->geom;
  int wn, n, m;
  int c, cc = 0;

  for (wn = 1; wn <= master->nwindows; wn++) {

    /*---------------------------------------------------------------*/
    /* Save the locations of surface time series points. The whole   */
    /* column is transferred since the exact location in the water   */
    /* column of the time series point may vary according to the     */
    /* reference point.                                              */
    for (n = 0; n < nts; ++n) {
      int *st = NULL, ssize;
      if (tslist[n].c <= 0) continue;
      c = tslist[n].c;
      if (c <=0 || c > geom->sgnum) continue;
      if (geom->fm[c].wn == wn) {
	ssize = stencil_size(geom, tslist[n].kernal, c);
	st = stencil(geom, c, &ssize, ST_SIZED, 0);
	for (m = 0; m < ssize; m++) {
	  c = st[m];
	  while (c != geom->zm1[c]) {
	    cc++;
	    c = geom->zm1[c];
	  }
	}
	i_free_1d(st);
      }
    }
    window[wn]->ns2m_ts = cc;
    if (cc) {
      if (window[wn]->s2m_ts) i_free_1d(window[wn]->s2m_ts);
      window[wn]->s2m_ts = i_alloc_1d(cc + 1);
    }
    cc = 1;
    for (n = 0; n < nts; ++n) {
      int *st = NULL, ssize;
      if (tslist[n].c <= 0) continue;
      c = tslist[n].c;
      if (c <=0 || c > geom->sgnum) continue;
      if (geom->fm[c].wn == wn) {
	ssize = stencil_size(geom, tslist[n].kernal, c);
	st = stencil(geom, c, &ssize, ST_SIZED, 0);
	for (m = 0; m < ssize; m++) {
	  c = st[m];
	  while (c != geom->zm1[c]) {
	    window[wn]->s2m_ts[cc] = geom->fm[c].sc;
	    cc++;
	    c = geom->zm1[c];
	  }
	}
	i_free_1d(st);
      }
    }
  }
  master->regf &= ~RS_TSSET;
}

/* END timeseries_init_w()                                          */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Close timeseries, re-initialise and populate the timeseries      */
/* structure.                                                       */
/*------------------------------------------------------------------*/
void ts_resetup(master_t *master, FILE *fp)
{
  geometry_t *geom = master->geom;
  dump_data_t *dumpdata = master->dumpdata;
  int n, i, onts, fr = forced_restart;
  char **dname;
  int *next;

  onts = nts;
  dname = (char **)malloc(onts * sizeof(char *));
  next = i_alloc_1d(onts);

  /*----------------------------------------------------------------*/
  /* Close the dumpfiles and free the dumplist                      */
  for (i = 0; i < onts; ++i) {
    if (tslist[i].c > 0) {
      dname[i] = (char *)malloc(sizeof(char)*MAXSTRLEN);
      strcpy(dname[i], tslist[i].pname);
      next[i] = tslist[i].tsout;
      fclose(tslist[i].fp);
      tslist[i].i = -1;
      tslist[i].c = 0;
    }
  }
  free((ts_point_t *) tslist);

  /*----------------------------------------------------------------*/
  /* Reread and initialise the new dumplist                         */
  forced_restart = 1;
  timeseries_init_m(fp, master, geom, window, dumpdata);
  ts_init(sched_get_even_by_name(schedule, "timeseries"));

  /*----------------------------------------------------------------*/
  /* Reset the next dump events                                     */
  for (i = 0; i < nts; i++) {

    /* Skip dump if the DA flag doesn't match for this dumpfile     */
    if (master->da && !(tslist[i].da_cycle & master->da)) continue;
    
    tslist[i].tsout = master->t;
    for (n = 0; n < onts; n++) {
      if (strcmp(dname[n], tslist[i].pname) == 0) {
	tslist[i].tsout = next[n];
      }
    }
  }

  forced_restart = fr;
  free((char **)dname);
  i_free_1d(next);
}

/* END ts_resetup()                                                 */
/*------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Initialised external data to be read into the timeseries          */
/*-------------------------------------------------------------------*/
void read_ts_data_init(FILE *fp,          /* Input parameters data   */
		       int n,             /* Timeseries number       */
		       master_t *master,  /* Master data             */
		       ts_point_t *ts     /* Timeseries data         */
		       )
{
  int i, tn;                 /* Counters                             */
  char buf[MAXSTRLEN];       /* Dummy buffer                         */
  char key[MAXSTRLEN];       /* Dummy buffer                         */
  char v_name[MAXSTRLEN];    /* Scaling variable name                */
  char f_name[MAXSTRLEN];    /* Forcing file name                    */
  geometry_t *geom = master->geom;

  /*-----------------------------------------------------------------*/
  /* Read the scaling variable and units                             */
  ts->ndata = 0;
  ts->dnvars = 0;
  ts->kernal = 1;
  sprintf(key, "TS%1d.data_file", n);
  if (prm_read_char(fp, key, f_name)) {
    char files[MAXNUMTSFILES][MAXSTRLEN];
    cstring *filenames;
    char *names[MAXSTRLEN * MAXNUMARGS];
    char *thresholds[MAXSTRLEN * MAXNUMARGS];
    int nvars;

    sprintf(v_name, "%c", '\0');
    sprintf(key, "TS%1d.data_name", n);
    prm_read_char(fp, key, v_name);
    
    sprintf(key, "TS%1d.data_kernal", n);
    if (prm_read_char(fp, key, buf)) {
      if (parseline(buf, names, MAXNUMARGS) == 2) {
	if (strcmp(names[1], "km") == 0)
	  ts->kernal = -(int)(atof(names[0]) * 1000.);
	else
	  ts->kernal = -atoi(names[0]);
      } else {
	ts->kernal = atoi(buf);
	/* Odd kernals only */
	if (ts->kernal > 1 && ts->kernal%2 == 0) ts->kernal--;
      }
    } else
      ts->kernal = 1;

    ts->metric = TS_CLOS;
    sprintf(key, "TS%1d.data_metric", n);
    if (prm_read_char(fp, key, buf)) {
      if (strcmp(buf, "CLOSEST") == 0)
	ts->metric = TS_CLOS;      
      if (strcmp(buf, "DIFF") == 0)
	ts->metric = TS_DIFF;      
      if (strcmp(buf, "MEAN") == 0)
	ts->metric = TS_MEAN;      
      if (strcmp(buf, "RMSE") == 0)
	ts->metric = TS_RMSE;      
      if (strcmp(buf, "CATEGORICAL") == 0)
	ts->metric = TS_CAT;
      if (strcmp(buf, "TRUE_SKILL") == 0)
	ts->metric = TS_HK;
      if (strcmp(buf, "CRITICAL_SUCCESS") == 0)
	ts->metric = TS_TS;
      if (strcmp(buf, "HIT_RATE") == 0)
	ts->metric = TS_H;
      if (strcmp(buf, "FALSE_ALARM_RATE") == 0)
	ts->metric = TS_F;
      if (strcmp(buf, "PREDICTION_RATE") == 0)
	ts->metric = TS_PRED;
      if (strcmp(buf, "GLIDER") == 0) {
	ts->metric = TS_GLIDER;
	if (strlen(v_name))
	  strcat(v_name, " temp salt dens");
	else
	  strcpy(v_name, "temp salt dens");
	ts->last = 0;
	if (!(master->numbers & GLIDER))
	  hd_quit("Set NUMBERS GLIDER to invoke the glider metric.\n");
      }
    }

    ts->ndata = parseline(f_name, (char **)files, MAXNUMTSFILES);
    filenames = (cstring *)malloc(sizeof(cstring) * ts->ndata);
    if ((ts->dnvars = parseline(v_name, names, MAXNUMARGS)) == 0)
      hd_quit("Can't file 'TS%1d.data_name'\n", n);
    ts->dvarids = i_alloc_2d(ts->ndata, ts->dnvars);

    ts->thresh = NULL;
    sprintf(key, "TS%1d.data_threshold", n);
    if (!(ts->metric & (TS_CLOS|TS_MEAN|TS_RMSE)) &&
	prm_read_char(fp, key, buf)) {
      ts->thresh = d_alloc_2d(2, ts->dnvars);
      tn = parseline(buf, thresholds, MAXNUMARGS);
      for (tn = 0; tn < ts->dnvars; tn++) {
	int found = 0;
	for (i = 0; i < master->ntr; i++) {
	  if (strcmp(thresholds[tn], master->trinfo_3d[i].name) == 0) {
	    ts->thresh[tn][0] = 0.0;
	    ts->thresh[tn][1] = (double)i;
	    found = 1;
	  }
	}
	if (!found) {
	  for (i = 0; i < master->ntrS; i++) {
	    if (strcmp(thresholds[tn], master->trinfo_2d[i].name) == 0) {
	      ts->thresh[tn][0] = 1.0;
	      ts->thresh[tn][1] = (double)i;
	      found = 1;
	    }
	  }
	}
	if (!found) {
	  ts->thresh[tn][0] = 2.0;
	  ts->thresh[tn][1] = atof(thresholds[tn]);
	}
      }
    }

    for (i = 0; i < ts->ndata; ++i)
     strcpy(filenames[i], ((char **)files)[i]);

    ts->tsdata = hd_ts_multifile_read(master, ts->ndata, filenames);
    if (ts->tsdata == NULL) {
      free(filenames);
      free(ts->tsdata);
      hd_warn("read_ts_data: Can't open forcing file %s.\n", f_name);
      return;
    }

    /* Check if density is in a glider file; if not read pressure    */
    /*  compute.                                                     */
    if(ts->metric & TS_GLIDER) {
      int found = 0;
      for (tn = 0; tn < ts->dnvars; tn++) {
	if (strcmp(names[tn], "dens") == 0) {
	  found = tn;
	  for (i = 0; i < ts->ndata; ++i) {
	    int varid;
	    if ((varid = ts_get_index(ts->tsdata[i], fv_get_varname(filenames[i], names[tn], buf)) >= 0))
	      found = 0;
	  }
	}
      }
      if (found) {
	hd_warn("read_ts_data: Can't find '%s' in glider file : computing density.\n", 
		names[found]);
	strcpy(names[found], "PRES");
      }
    }

    for (tn = 0; tn < ts->dnvars; tn++) {
      if (hd_ts_multifile_get_index(ts->ndata, ts->tsdata, filenames, names[tn], ts->dvarids[tn]) > 0) {
	for(i = 0; i < ts->ndata; i++) {
	  if (ts->dvarids[tn][i] < 0)
	    hd_warn("read_ts_data: Can't find variable '%s' in timeseries file %s.\n", 
		    names[tn], filenames[i]);
	}
      } else
	hd_warn("read_ts_data: Can't find variable '%s' in any timeseries file %s.\n", 
		names[tn], f_name);
    }

    ts->val = d_alloc_1d(ts->dnvars);
    ts->obs = d_alloc_1d(ts->dnvars);
    ts->minv = ts->maxv = ts->nvals = NULL;
    if (ts->metric & TS_GLIDER) {
      ts->minv = d_alloc_1d(ts->dnvars);
      ts->maxv = d_alloc_1d(ts->dnvars);
      ts->nvals = d_alloc_1d(ts->dnvars);
    }
    ts->data = (double **)p_alloc_1d(ts->dnvars);
    ts->ddim = i_alloc_1d(ts->dnvars);
    ts->dvars = (char **)malloc(ts->dnvars * sizeof(char *));
    ts->dunits = (char **)malloc(ts->dnvars * sizeof(char *));
    for (tn = 0; tn < ts->dnvars; tn++) {
      ts->dvars[tn] = (char *)malloc(sizeof(char)*MAXSTRLEN);
      ts->dunits[tn] = (char *)malloc(sizeof(char)*MAXSTRLEN);
    }

    for (tn = 0; tn < ts->dnvars; tn++) {
      int found = 0;
      ts->ddim[tn] = 1;
      if (strcmp(names[tn], "eta") == 0) {
	ts->data[tn] = master->eta;
	strcpy(ts->dvars[tn], "eta");
	strcpy(ts->dunits[tn], "m");
	ts->ddim[tn] = 0;
      } else if (strcmp(names[tn], "u1") == 0) {
	ts->data[tn] = master->u1;
	strcpy(ts->dvars[tn], "u1");
	strcpy(ts->dunits[tn], "ms-1");
      } else if (strcmp(names[tn], "u2") == 0) {
	ts->data[tn] = master->u2;
	strcpy(ts->dvars[tn], "u2");
	strcpy(ts->dunits[tn], "ms-1");
      } else if (strcmp(names[tn], "dens") == 0) {
	ts->data[tn] = master->dens_0;
	strcpy(ts->dvars[tn], "N2");	
	strcpy(ts->dunits[tn], "s-2");
      } else if (strcmp(names[tn], "PRES") == 0) {
	ts->data[tn] = master->dens_0;
	strcpy(ts->dvars[tn], "PRES");	
	strcpy(ts->dunits[tn], "dbar");
      } else if (strcmp(names[tn], "depth") == 0) {
	ts->data[tn] = geom->cellz;
	strcpy(ts->dvars[tn], "depth");
	strcpy(ts->dunits[tn], "m");
      } else {
	for (i = 0; i < master->ntr; i++) {
	  if (strcmp(names[tn], master->trinfo_3d[i].name) == 0) {
	    ts->data[tn] = master->tr_wc[i];
	    sprintf(ts->dvars[tn], "%s", names[tn]);
	    strcpy(ts->dunits[tn], master->trinfo_3d[i].units);
	    found = 1;
	  }
	}
	if (found == 0) {
	  for (i = 0; i < master->ntrS; i++) {
	    if (strcmp(names[tn], master->trinfo_2d[i].name) == 0) {
	      ts->data[tn] = master->tr_wcS[i];
	      ts->ddim[tn] = 0;
	      sprintf(ts->dvars[tn], "%s", names[tn]);
	      strcpy(ts->dunits[tn], master->trinfo_2d[i].units);
	      found = 1;
	    }
	  }
	}
      }
    }
  }
}

/* END read_ts_data_init()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Returns a stencil size                                            */
/*-------------------------------------------------------------------*/
int stencil_size(geometry_t *geom, int kernal, int c)
{
  int ssize;

  if (kernal > 0)
    ssize = kernal;
  else {   /* Convert from a length scale to number of cells         */
    ssize = -kernal / sqrt(geom->cellarea[geom->m2d[c]]);
    if (ssize > 1 && ssize%2 == 0) ssize--;
    ssize = max(ssize, 1);
  }
  return(ssize);
}

/* END stencil_size()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Read external data and compare to model data                      */
/*-------------------------------------------------------------------*/
int read_ts_data(master_t *master, ts_point_t tslist, double t, int c)
{
  geometry_t *geom = master->geom;
  int cc, cs, c2, cg, tn, tid, sid;
  double val, d1, d2;
  double tu, td;
  double temp, sal, p, dens;
  int *st = NULL, ssize;
  int posf;

  /* Get the kernal                                                  */
  ssize = stencil_size(geom, tslist.kernal, c);
  c2 = geom->m2d[c];
  st = stencil(geom, c, &ssize, ST_SIZED, 0);

  /* Initialize glider data if it enters a new column                */
  if (master->glider && tslist.metric & TS_GLIDER) {
    if (c == tslist.last)
      return(1);
    if (c2 != geom->m2d[tslist.last]) {
      cs = c2;
      while (cs != geom->zm1[cs]) {
	master->glider[cs] = 0.0;
	cs = geom->zm1[cs];
      }
    }
    tslist.last = c ;
  }

  /* Loop though the variables specified for comparison              */
  for (tn = 0; tn < tslist.dnvars; tn++) {
    int n2f = (strcmp(tslist.dvars[tn], "N2") == 0 || strcmp(tslist.dvars[tn], "PRES") == 0) ? 1 : 0;
    /* Get the observation value                                     */
    if (tslist.ddim[tn]) {

      /* For glider metrics we compute the average denity in a cell  */
      if (tslist.metric & TS_GLIDER) {
	/* Add the density to the glider array                       */
	char var[MAXSTRLEN];
	double tzm1, nvals;
	int cu, cd, cb, cud, cdd;
	int zp1 = geom->zp1[c];
	int zm1 = geom->zm1[c];

	/* If the glider has not entered this cell, then average its */
	/* data while in the cell.                                   */
	if (master->glider[c] == 0.0) {
	  strcpy(var, tslist.dvars[tn]);
	  /* Compute buoyancy frequency                              */
	  if (n2f) {
	    int densf = 0;
	    if (strcmp(var, "N2") == 0) 
	      strcpy(var, "dens");
	    else
	      densf = 1;
	    master->glider[c] = average_glider_data(master, &tslist, t, c, tn, &cu, &cd, &tu, &td, &nvals);
	    tslist.nvals[tn] = nvals;

	    /* Compute glider density if not supplied in the file    */
	    if (densf) {
	      p = master->glider[c] * 1e6;
	      /* eos2 expects pressure as N/m2 and converts to bars  */
	      /* (divides by 1e5), so multiply p by 1e6 to retrieve  */
	      /* bars in eos2().                                     */
	      temp = average_glider_data(master, &tslist, tzm1, zm1, tid, &cud, &cdd, &d1, &d2, &nvals);
	      sal = average_glider_data(master, &tslist, tzm1, zm1, sid, &cud, &cdd, &d1, &d2, &nvals);
	      eos2(sal, temp, p, &dens, &master->glider[c]);
	    }

	    /* Get the average density in cells below c              */
	    /* Get the orientation of the cell below or above c      */
	    if (cd == geom->zm1[c]) {
	      cg = c;
	      zm1 = cd;
	      tzm1 = td;
	      posf = 1;
	    } else if (cu == geom->zm1[c]) {
	      cg = c;
	      zm1 = cu;
	      tzm1 = tu;
	      posf = 1;
	    } else if (c == geom->zm1[cu]) {
	      zm1 = cg = cu;
	      tzm1 = tu;
	      posf = 0;
	    } else if (c == geom->zm1[cd]) {
	      zm1 = cg = cd;
	      tzm1 = td;
	      posf = 0;
	    } else {
	      /*printf("Can't find c\n");*/
	      val = 0.0;
	      break;
	    }
	    /* Get the density in this adjacent cell                 */
	    master->glider[zm1] = average_glider_data(master, &tslist, tzm1, zm1, tn, &cud, &cdd, &d1, &d2, &nvals);
	    if (densf) {
	      p = master->glider[zm1] * 1e6;
	      temp = average_glider_data(master, &tslist, tzm1, zm1, tid, &cud, &cdd, &d1, &d2, &nvals);
	      sal = average_glider_data(master, &tslist, tzm1, zm1, sid, &cud, &cdd, &d1, &d2, &nvals);
	      eos2(sal, temp, p, &dens, &master->glider[zm1]);
	    }
	    val = buoyancy_frequency2_m(master, geom, master->glider, cg);
	  } else {
	    /* Average the variable in the cell. Save temp and salt  */
	    /* for use in computing density if required.             */
	    val = average_glider_data(master, &tslist, t, c, tn, &cu, &cd, &d1, &d2, &nvals);
	    tslist.nvals[tn] = nvals;
	    if (strcmp(tslist.dvars[tn], "temp") == 0) {
	      tid = tn;
	      temp = val;
	    }
	    if (strcmp(tslist.dvars[tn], "salt") == 0) {
	      sid = tn;
	      sal = val;
	    }
	  }
	}
      } else {	
	val = hd_ts_multifile_eval_xyz(tslist.ndata, tslist.tsdata, tslist.dvarids[tn],
				       t, tslist.x, tslist.y, tslist.z);
      }
      tslist.obs[tn] = val;

      /* Set ghost cells in the model (on the master)                */
      for (cc = 1; cc <= geom->nbpt; cc++) {
	int c1, c2;
	c1 = geom->bpt[cc];
	c2 = geom->bin[cc];
	tslist.data[tn][c1] = tslist.data[tn][c2];
      }
    } else {
      val = hd_ts_multifile_eval_xy(tslist.ndata, tslist.tsdata, tslist.dvarids[tn],
				    t, tslist.x, tslist.y);
      /* Set ghost cells in the model (on the master)                */
      for (cc = 1; cc <= geom->nbptS; cc++) {
	int c1, c2;
	c1 = geom->bpt[cc];
	c2 = geom->bin[cc];
	tslist.data[tn][c1] = tslist.data[tn][c2];
      }
    }

    /* Compare to the model within the kernal                        */
    if (tslist.metric & (TS_CLOS)) {
      d1 = HUGE;
      for (cc = 0; cc < ssize; cc++) {
	cs = st[cc];
	if (fabs(val - tslist.data[tn][cs]) < d1) {
	  d1 = fabs(val - tslist.data[tn][cs]);
	  tslist.val[tn] = (tslist.metric & TS_DIFF) ? val - tslist.data[tn][cs] : tslist.data[tn][cs];
	}
      }
    } else if (tslist.metric & TS_DIFF && tslist.thresh != NULL) {
      double thresh, d1 = 0.0;
      tslist.val[tn] = 0.0;
      for (cc = 0; cc < ssize; cc++) {
	cs = st[cc];
	if (tslist.thresh[tn][0] == 0)
	  thresh = master->tr_wc[(int)tslist.thresh[tn][1]][cs];
	else if (tslist.thresh[tn][0] == 1)
	  thresh = master->tr_wcS[(int)tslist.thresh[tn][1]][cs];
	else
	  thresh = tslist.thresh[tn][1];
	if (fabs(val - tslist.data[tn][cs]) < thresh) tslist.val[tn] += 1.0;
      }
      tslist.val[tn] /= (double)ssize;
    } else if (tslist.metric & TS_MEAN) {
      d1 = 0.0;
      for (cc = 0; cc < ssize; cc++) {
	cs = st[cc];
	d1 += tslist.data[tn][cs];
      }
      d1 /= (double)(ssize);
      tslist.val[tn] = d1;
    } else if (tslist.metric & TS_RMSE) {
      d1 = 0.0;
      for (cc = 0; cc < ssize; cc++) {
	cs = st[cc];
	d2 = (val - tslist.data[tn][cs]);
	d1 += (d2 * d2);
      }
      tslist.val[tn] = sqrt(d1 / (double)(ssize));
    } else if (tslist.metric & TS_GLIDER) {
      if (n2f) {
	/* Buoyancy frequency computation; use the same cells as the */
	/* glider used.                                              */
	cs = (posf) ? c : geom->zp1[c];
	tslist.val[tn] = buoyancy_frequency2_m(master, geom, master->dens, cs);
      }
      else
	tslist.val[tn] = tslist.data[tn][c];

      tslist.minv[tn] = HUGE;
      tslist.maxv[tn] = 0.0;
      for (cc = 0; cc < ssize; cc++) {
	cs = st[cc];
	if (n2f) {
	  if (!posf) cs = geom->zp1[cs];
	  d2 = buoyancy_frequency2_m(master, geom, master->dens, cs);
	} else
	  d2 = tslist.data[tn][cs];
	tslist.minv[tn] = min(tslist.minv[tn], d2);
	tslist.maxv[tn] = max(tslist.maxv[tn], d2);
      }
    } else if (tslist.metric & (TS_CAT|TS_HK|TS_TS|TS_H|TS_F|TS_PRED) && tslist.thresh != NULL) {
      /* Categorical contingency table                               */
      /* Hit = 0, miss = 1, false alarm = 2, reject = 3              */
      /* Hit = observation and prediction above the threshold        */
      /* Reject = observation and prediction below the threshold     */
      /* Miss = observation above and and prediction below threshold */
      /* False = observation below and prediction above threshold    */
      double cat[4], thresh;
      tslist.val[tn] = -1.0;
      for (cc = 0; cc < 4; cc++) cat[cc] = 0.0;
      for (cc = 0; cc < ssize; cc++) {
	cs = st[cc];
	if (tslist.thresh[tn][0] == 0)
	  thresh = master->tr_wc[(int)tslist.thresh[tn][1]][cs];
	else if (tslist.thresh[tn][0] == 1)
	  thresh = master->tr_wcS[(int)tslist.thresh[tn][1]][cs];
	else
	  thresh = tslist.thresh[tn][1];
	if (val >= thresh && tslist.data[tn][cs] >= thresh) cat[0] += 1.0;
	if (val >= thresh && tslist.data[tn][cs] < thresh) cat[1] += 1.0;
	if (val < thresh && tslist.data[tn][cs] >= thresh) cat[2] += 1.0;
	if (val < thresh && tslist.data[tn][cs] < thresh) cat[3] += 1.0;
      }
      if (tslist.metric & TS_CAT) {
	if (cat[0] >= cat[1] && cat[0] >= cat[2] && cat[0] >= cat[3]) tslist.val[tn] = 0;
	if (cat[1] >= cat[0] && cat[1] >= cat[2] && cat[1] >= cat[3]) tslist.val[tn] = 1;
	if (cat[2] >= cat[1] && cat[2] >= cat[0] && cat[2] >= cat[3]) tslist.val[tn] = 2;
	if (cat[3] >= cat[1] && cat[3] >= cat[2] && cat[3] >= cat[0]) tslist.val[tn] = 3;
      } else if (tslist.metric & TS_HK) {    
	d1 = cat[0] + cat[1];
	d2 = cat[2] + cat[3];
	if (d1 && d2)
	  tslist.val[tn] = cat[0] / d1 - cat[2] / d2;
	else 
	  tslist.val[tn] = -1.0;
      } else if (tslist.metric & TS_TS) {
	d1 = cat[0] + cat[1] + cat[2];
	tslist.val[tn] = (d1) ? cat[0] / d1 : -1.0;
      } else if (tslist.metric & TS_H) {
	d1 = cat[0] + cat[1];
	tslist.val[tn] = (d1) ? cat[0] / d1 : -1.0;
      } else if (tslist.metric & TS_F) {
	d1 = cat[2] + cat[3];
	tslist.val[tn] = (d1) ? cat[2] / d1 : -1.0;
      } else if (tslist.metric & TS_PRED) {
	d1 = cat[0] + cat[1] + cat[2] + cat[3];
	tslist.val[tn] = (d1) ? (cat[0] + cat[3]) / d1 : -1.0;
      }
    }

    /*printf("%f %d %f : %f %f %f : %f %f\n",t,tn,val,tslist.x, tslist.y, tslist.z, tslist.data[tn][c], tslist.val[tn]);*/
  }

  i_free_1d(st);
  return(0);
}

/* END read_ts_data()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to average glider data onto a model cell                  */
/*-------------------------------------------------------------------*/
double average_glider_data(master_t *master,   /* Master data        */
			   ts_point_t *tslist, /* Time series        */
			   double t,           /* Time               */
			   int c,              /* Glider location    */
			   int tn,             /* Number of variable */
			   int *cu,            /* Cell above c       */
			   int *cd,            /* Cell below c       */
			   double *tu,         /* Time when in cu    */ 
			   double *td,         /* Time when in cd    */ 
			   double *nvals       /* Number of samples  */
			   )
{
  int m, n = 0;
  int cp, cg = c;
  int r1, r2;
  int i, j;
  double val = 0.0, v1;
  double tg = t;
  double frac;
  double x, y, z;
  timeseries_t **tsdata = tslist->tsdata;
  timeseries_t *loc_ts = &tslist->ts;
  datafile_t *df;
  df_variable_t *v;
  int fv = 1, firstval = 0;

  for (m = 0; m < tslist->ndata; ++m) {
    timeseries_t *ts = tsdata[m];

    if (tslist->dvarids[tn][i] < 0) continue;
    df = ts->df;
    v = df_get_variable(df, tslist->dvarids[tn][m]);

    /* Get the record number bracketing time t                       */
    df_find_record(df, tg, &r1, &r2, &frac);

    /* Read in glider records previous to time t until a new cell is */
    /* entered.                                                      */
    while (cg == c) {
      df_read_records(df, v, r1, 1);        /* Read the record       */
      v1 = v->data[r1 - v->start_record];   /* Record value          */
      if (fv && !isnan(v1)) {
	val += v1;                          /* Average value         */
	n++;
	if (firstval) fv = 0;
      }
      r1 = max(r1 - 1, 0);                  /* Increment the record  */
      tg = df->records[r1];                 /* Glider time of record */
      /* Get the new location in the grid                            */
      cp = cg;
      if (!(cg = get_glider_loc(master, tslist, loc_ts, tg))) {
	break;
      }
    }
    *cu = cg;
    *tu = tg;

    /* Read in glider records subsequrnt to time t until a new cell  */
    /* is entered.                                                   */
    cg = c;
    tg = t;
    fv = 1;
    while (cg == c) {
      df_read_records(df, v, r2, 1);        /* Read the record       */
      v1 = v->data[r2 - v->start_record];   /* Record value          */
      if (fv && !isnan(v1)) {
	val += v1;                          /* Average value         */
	n++;
	if (firstval) fv = 0;
      }
      r2 = min(r2 + 1, v->nrecords - 1);    /* Increment the record  */
      tg = df->records[r2];                 /* Glider time of record */
      cp = cg;
      if (!(cg = get_glider_loc(master, tslist, loc_ts, tg))) {
	break;
      }
    }
    *cd = cg;
    *td = tg;
  }
  val = (n) ? val / (double)n : NaN;
  *nvals = (double)n;
  return(val);
}

/* END average_glider_data()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Gets a new glider location in the grid                            */
/*-------------------------------------------------------------------*/
int get_glider_loc(master_t *master, ts_point_t *tslist, timeseries_t *loc_ts, double t)
{
  int c, cs, i, j;
  double x, y, z, ts_t = t;

  tm_change_time_units(master->timeunit, loc_ts->t_units, &ts_t, 1);

  if (ts_has_time(loc_ts, ts_t)) {
    x = ts_eval(loc_ts, tslist->varids[0], ts_t);
    y = ts_eval(loc_ts, tslist->varids[1], ts_t);
    z = ts_eval(loc_ts, tslist->varids[2], ts_t);
    if (loc_ts->df->variables[tslist->varids[2]].z_is_depth)
      z = -z;

    if ((c = hd_grid_xytoij(master, x, y, &i, &j))) {
      cs = geom->m2d[c];
      c = d2c(master, z, tslist->v_offset, cs);
    } else
      return(0);
  } else
    return(0);
  return(c);
}

/* END get_glider_loc()                                              */
/*-------------------------------------------------------------------*/
