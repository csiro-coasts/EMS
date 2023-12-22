/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/forcings/swan.c
 *  
 *  Description:
 *  Event scheduler routines for Atmospheric cloud.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: swan.c 5873 2018-07-06 07:23:48Z riz008 $
 *
 */

#include <math.h>
#include <string.h>
#include "hd.h"

#define SWAN_INIT    0x01
#define SWAN_STEP    0x02

static int swan_init(sched_event_t *event);
double swan_event(sched_event_t *event, double t);
static void swan_cleanup(sched_event_t *event, double t);
void write_swan_params(master_t *master, char *fname, char *mname);
void bdry_wave(master_t *master, geometry_t *geom, swan_data_t *data);


/* Functions for reading the schedule the swan forcings. */
void swan_couple_init(master_t *master)
{
  geometry_t *geom = master->geom;
  swan_data_t *data = NULL;
  int c, cc, cw, nb, n;

  if (!(master->do_wave & W_SWAN)) return;

  /* Read parameters */
  prm_set_errfn(hd_silent_warn);

  data = (swan_data_t *)malloc(sizeof(swan_data_t));

  data->master = master;
  data->dt = master->wavedt;
  data->nw2c = geom->nw2c;
  data->swan_hs = master->swan_hs;
  data->swan_eta = d_alloc_1d(geom->szcS);
  data->swan_uav = d_alloc_1d(geom->szcS);
  data->swan_vav = d_alloc_1d(geom->szcS);
  data->swan_wx = d_alloc_1d(geom->szcS);
  data->swan_wy = d_alloc_1d(geom->szcS);
  data->swan_dep = d_alloc_1d(geom->szcS);
  data->swan_amp = d_alloc_1d(geom->szcS);
  data->swan_per = d_alloc_1d(geom->szcS);
  data->swan_dir = d_alloc_1d(geom->szcS);
  data->swan_ub = d_alloc_1d(geom->szcS);
  data->swan_Fx = d_alloc_1d(geom->szcS);
  data->swan_Fy = d_alloc_1d(geom->szcS);
  data->swan_ste1 = d_alloc_1d(geom->szcS);
  data->swan_ste2 = d_alloc_1d(geom->szcS);
  data->swan_Kb = d_alloc_1d(geom->szcS);
  data->swan_k = d_alloc_1d(geom->szcS);
  data->swan_fwcapx = d_alloc_1d(geom->szcS);
  data->swan_fwcapy = d_alloc_1d(geom->szcS);
  data->swan_fbrex = d_alloc_1d(geom->szcS);
  data->swan_fbrey = d_alloc_1d(geom->szcS);
  data->swan_fbotx = d_alloc_1d(geom->szcS);
  data->swan_fboty = d_alloc_1d(geom->szcS);
  data->swan_fsurx = d_alloc_1d(geom->szcS);
  data->swan_fsury = d_alloc_1d(geom->szcS);
  data->swan_wfdx = d_alloc_1d(geom->szcS);
  data->swan_wfdy = d_alloc_1d(geom->szcS);
  data->swan_wovsx = d_alloc_1d(geom->szcS);
  data->swan_wovsy = d_alloc_1d(geom->szcS);
  /*
    data->swan_frolx = d_alloc_1d(geom->szcS);
    data->swan_froly = d_alloc_1d(geom->szcS);
  */

  /* Set the default wavenumber (must be non-zero).                  */
  /* The centres updated with wave data are vertices in SWAN, and    */
  /* these must map to three centres. This means land boundary cells */
  /* are excluded, and will be assigned zero wavenumber, which leads */
  /* to errors in wave_av_sbc(). To avoid this we set the wavenumber */
  /* here to the same default used in SWAN.                          */
  if (master->wave_k) {
    if ((n = tracer_find_index("wave_k", master->ntrS, master->trinfo_2d)) >= 0) {
      double per, kmin;
      /* Using wave period fill value                                */
      if ((nb = tracer_find_index("wave_period", master->ntrS, 
				  master->trinfo_2d)) >= 0) {
	per = master->trinfo_2d[nb].fill_value_wc;
	kmin = 4.0 * PI/ (master->g * per * per);
      }
      /* Using wavenumber fill value                                 */
      for (cc = 1; cc <= geom->b2_t; cc++) {
	c = geom->w2_t[cc];
	master->wave_k[c] = master->trinfo_2d[n].fill_value_wc;
	master->wave_k[c] = kmin;	
      }
    }
  }
  /* Set initial wave conditions. These are only used if OBCs are    */
  /* passive.                                                        */
  if ((n = tracer_find_index("wave_amp", master->ntrS, master->trinfo_2d)) >= 0) {
    tracer_re_read(&master->trinfo_2d[n], params->prmfd, INTER);
    for (cc = 1; cc <= geom->nw2c; cc++) {
      c = geom->w2c[cc];
      data->swan_amp[cc] = master->trinfo_2d[n].fill_value_wc;
    }
  }
  if ((n = tracer_find_index("wave_period", master->ntrS, master->trinfo_2d)) >= 0) {
    tracer_re_read(&master->trinfo_2d[n], params->prmfd, INTER);
    for (cc = 1; cc <= geom->nw2c; cc++) {
      c = geom->w2c[cc];
      data->swan_per[cc] = master->trinfo_2d[n].fill_value_wc;
    }
  }
  if ((n = tracer_find_index("wave_dir", master->ntrS, master->trinfo_2d)) >= 0) {
    tracer_re_read(&master->trinfo_2d[n], params->prmfd, INTER);
    for (cc = 1; cc <= geom->nw2c; cc++) {
      c = geom->w2c[cc];
      data->swan_dir[cc] = master->trinfo_2d[n].fill_value_wc * PI / 180.0;
    }
  }
  /*
  for (cc = 1; cc <= geom->nw2c; cc++) {
    c = geom->w2c[cc];
    data->swan_amp[cc] = 0.02;
    data->swan_per[cc] = 2.0;
    data->swan_dir[cc] = 270.0 * PI / 180.0;
  }
  */
  /* Get the open boundaries                                         */
  data->no2_t = 0;
  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    if (open->options & OP_WAVES) {
      for (cc = 1; cc <= open->no2_t; cc++) {
	c = open->obc_t[cc];
	if (geom->c2w[c]) data->no2_t++;
      }
    }
  }

  data->obc_t = i_alloc_1d(data->no2_t+1);
  nb = 1;
  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    if (open->options & OP_WAVES) {
      for (cc = 1; cc <= open->no2_t; cc++) {
	c = open->obc_t[cc];
	cw = geom->c2w[c];
	if (cw) data->obc_t[nb++] = cw;
      }
    }
  }

#if defined(HAVE_WAVE_MODULE)
  data->wave = wave_build_m(data, master->prmfd);
#endif

  data->flag = SWAN_INIT;


  /* Check the open boundaries for wave data                         */
  for (n = 0; n < geom->nobc; ++n) {
    open_bdrys_t *open = geom->open[n];
    if (open->bcond_wav & FILEIN) {
      hd_ts_multifile_check(open->ntsfiles, open->tsfiles, open->filenames, 
			    "wave_amp", schedule->start_time, schedule->stop_time);
      hd_ts_multifile_check(open->ntsfiles, open->tsfiles, open->filenames, 
			    "wave_period", schedule->start_time, schedule->stop_time);
      hd_ts_multifile_check(open->ntsfiles, open->tsfiles, open->filenames, 
			    "wave_dir", schedule->start_time, schedule->stop_time);
    }
  }
  sched_register(schedule, "swan", swan_init,
                 swan_event, swan_cleanup, data, NULL, NULL);
}

static int swan_init(sched_event_t *event)
{
  return 1;
}

double swan_event(sched_event_t *event, double t)
{
  swan_data_t *data = (swan_data_t *)schedGetPublicData(event);
  master_t *master = data->master;
  geometry_t *geom = master->geom;
  int cc, c, c1, cw, n;

  /* Don't call the swan routine if this event is called from the    */
  /* scheduler setup.                                                */
  if (data->flag & SWAN_INIT) {
    data->flag = SWAN_STEP;
    return event->next_event;
  }

  if (t >= (event->next_event - SEPS)) {
#if defined(HAVE_WAVE_MODULE)

    wave_t *wave = data->wave;
    /* Map the variables required by SWAN into dummy arrays          */
    for (cc = 1; cc <= geom->nw2c; cc++) {
      c = geom->w2c[cc];
      data->swan_eta[cc] = master->eta[c];
      data->swan_uav[cc] = master->u[c];
      data->swan_vav[cc] = master->v[c];
      if (master->swind1 && master->swind2) {
	data->swan_wx[cc] = master->swind1[c];
	data->swan_wy[cc] = master->swind2[c];
      }
      data->swan_dep[cc] = -geom->botz[c];
    }
    data->swan_hs = master->swan_hs;

    /* Update the open bounries if required                          */
    bdry_wave(master, geom, data);

    /* Invoke the SWAN wave coupling                                 */
    swan_step(data->wave);

    /* Get the wave current friction                               
       for (cc = 1; cc <= wincon->vca2; cc++) {
       c = wincon->s2[cc];
       if (wave->do_Cd) {
       ustrcw_estim(data->wave, c);
       }
       }
    */
    /* Map the SWAN variables back to arrays                         */
    for (cc = 1; cc <= geom->nw2c; cc++) {
      c = geom->w2c[cc];
      if (wave->do_amp) master->wave_amp[c] = data->swan_amp[cc];
      if (wave->do_per) master->wave_period[c] = data->swan_per[cc];
      if (wave->do_dir) master->wave_dir[c] = data->swan_dir[cc];
      if (wave->do_ub) master->wave_ub[c] = data->swan_ub[cc];
      if (wave->do_wif) {
	master->wave_Fx[c] = data->swan_Fx[cc];
	master->wave_Fy[c] = data->swan_Fy[cc];
      }
      if (wave->do_stokes) {
	master->wave_ste1[c] = data->swan_ste1[cc];
	master->wave_ste2[c] = data->swan_ste2[cc];
      }
      if (wave->do_Kb)
	master->wave_Kb[c] = data->swan_Kb[cc];
      if (wave->do_k)
      	master->wave_k[c] = data->swan_k[cc];
      if (wave->do_noncon) {
	master->wave_fwcapx[c] = data->swan_fwcapx[cc];
	master->wave_fwcapy[c] = data->swan_fwcapy[cc];
	master->wave_fbrex[c] = data->swan_fbrex[cc];
	master->wave_fbrey[c] = data->swan_fbrey[cc];
	master->wave_fbotx[c] = data->swan_fbotx[cc];
	master->wave_fboty[c] = data->swan_fboty[cc];
	master->wave_fsurx[c] = data->swan_fsurx[cc];
	master->wave_fsury[c] = data->swan_fsury[cc];
	master->wave_wfdx[c] = data->swan_wfdx[cc];
	master->wave_wfdy[c] = data->swan_wfdy[cc];
	master->wave_wovsx[c] = data->swan_wovsx[cc];
	master->wave_wovsy[c] = data->swan_wovsy[cc];
	/*
 	master->wave_fwcapx[c] = 0.0;
 	master->wave_fwcapy[c] = 0.0;
 	master->wave_fbrex[c] = 0.0;
 	master->wave_fbrey[c] = 0.0;
 	master->wave_fbotx[c] = 0.0;
 	master->wave_fboty[c] = 0.0;
 	master->wave_fsurx[c] = 0.0;
 	master->wave_fsury[c] = 0.0;
 	master->wave_wfdx[c] = 0.0;
 	master->wave_wfdy[c] = 0.0;
 	master->wave_wovsx[c] = 0.0;
 	master->wave_wovsy[c] = 0.0;
	*/
	/*
	master->wave_frolx[c] = data->swan_frolx[cc];
	master->wave_froly[c] = data->swan_froly[cc];
	*/
	master->wave_frolx[c] = 0.0;
	master->wave_froly[c] = 0.0;
      }
    }
    event->next_event += data->dt;
  }

#endif
  return event->next_event;
}

void swan_cleanup(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  swan_data_t *data = (swan_data_t *)schedGetPrivateData(event);

  if (data != NULL) {
    d_free_1d(data->swan_eta);
    d_free_1d(data->swan_uav);
    d_free_1d(data->swan_vav);
    d_free_1d(data->swan_wx);
    d_free_1d(data->swan_wy);
    d_free_1d(data->swan_dep);
    d_free_1d(data->swan_amp);
    d_free_1d(data->swan_per);
    d_free_1d(data->swan_dir);
    d_free_1d(data->swan_ub);
    d_free_1d(data->swan_Fx);
    d_free_1d(data->swan_Fy);
    d_free_1d(data->swan_ste1);
    d_free_1d(data->swan_ste2);
    d_free_1d(data->swan_Kb);
    d_free_1d(data->swan_k);
    d_free_1d(data->swan_fwcapx);
    d_free_1d(data->swan_fwcapy);
    d_free_1d(data->swan_fbrex);
    d_free_1d(data->swan_fbrey);
    d_free_1d(data->swan_fbotx);
    d_free_1d(data->swan_fboty);
    d_free_1d(data->swan_fsurx);
    d_free_1d(data->swan_fsury);
    d_free_1d(data->swan_wfdx);
    d_free_1d(data->swan_wfdy);
    d_free_1d(data->swan_wovsx);
    d_free_1d(data->swan_wovsy);
    /*
    d_free_1d(data->swan_frolx);
    d_free_1d(data->swan_froly);
    */
    free(data);
  }
}

/*-------------------------------------------------------------------*/
/* Prints the dual of the COMPAS mesh (i.e. the triangulation) to    */
/* files with format compatible for input into SWAN.                 */
/* Typically a vertex will map to <= 2 centres for land boundaries.  */
/* In some cases (e.g. where quad grids are integrated in the mesh)  */
/* there may be 3 centre maps however. Interior vertices should      */
/* always map to 3 centres.                                          */
/* Note: SWAN will fail in SwanCheckGrid.f90 if the number of non-   */
/* boundary cells (vmark=0) around a vertex is smaller than 4 or     */
/* larger than 10. We force a cell to be a boundary cell (vmark=1)   */
/* if cells are < 4, but the > 10 criterion holds.                   */
/* vmark is the last entry in the .node file:                        */
/* vmark = 0 : interior cell                                         */
/* vmark = 1 : land boundary cell                                    */
/* vmark > 1 : open boundary cell                                    */
/* Files output:                                                     */
/* .node : Vertex coordinates and vmark                              */
/* .ele : Indices in .node of triangle vertices                      */
/* .bth : bathymetry                                                 */
/* _swan.txt : SWAN triangulation for plotting                       */
/*-------------------------------------------------------------------*/
void write_swan_mesh(master_t *master,
		     geometry_t **window,
		     window_t **windat,
		     win_priv_t **wincon
		     )
{
  geometry_t *geom = master->geom;
  parameters_t *params = master->params;
  FILE *np, *ep, *bp, *op;
  char buf[MAXSTRLEN], key[MAXSTRLEN];
  int c, cc, c1, cg, v, vv, nc, nt, n;
  int *mask;
  int mk;
  int mode = 0;
  int checkf = 1;

  if (!params->doswan && !(params->do_wave & W_SWAN)) return;
  if (params->doswan) mode |= 1;
  if (params->do_wave & W_SWAN) mode |= 2;
  geom->c2w = i_alloc_1d(geom->szcS);
  memset(geom->c2w, 0, geom->szcS * sizeof(int));
  mask = i_alloc_1d(geom->szcS);
  memset(mask, 0, geom->szcS * sizeof(int));
  if (mode & 2) {
    geom->nw2c = 0;
    geom->w2c = i_alloc_1d(geom->szcS);
    memset(geom->w2c, 0, geom->szcS * sizeof(int));
  }
  if (mode & 1) mesh_ofile_name(params, key);

  /* Get the mask for valid centres. Note: we do not include         */
  /* boundary vertices, as these will map to less than 3 centres.    */
  nc = nt = 0;
  for (vv = 1; vv <= geom->v2_e2; vv++) {
    v = geom->w2_e2[vv];
    c1 = 0;
    for (cc = 1; cc <= geom->nvc[v]; cc++) {
      c = geom->v2c[v][cc];
      if (!mask[c] && c && geom->nvc[v] == 3) {
	nc++;
	mask[c] = 1;
      }
      if (!c || geom->nvc[v] != 3) c1 = 1;
    }
    if (!c1) nt++;
  }

  if (nc != geom->b2_t) 
    hd_warn("SWAN convert: swan mesh size(%d) != COMPAS size(%d)\n",
	    nc, geom->b2_t);

  /* Print the triangulation of non-boundary cells                   */
  if (mode & 1) {
    sprintf(buf, "%s.node", key);
    np = fopen(buf, "w");
    sprintf(buf, "%s.bth", key);
    bp = fopen(buf, "w");
    fprintf(np, "%d 2 0 1\n", nc);
  }
  geom->nw2c = nc;
  nc = 1;
  for (cc = 1; cc <= geom->v2_t; cc++) {
    c = geom->w2_t[cc];

    if (mask[c]) {
      mk = 0;
      for (n = 1; n <= geom->npe[c]; n++) {
	c1 = geom->c2c[n][c];
	if (geom->wgst[c1]) {
	  mk = 1;
	  break;
	}
      }
      if (mode & 1) {
	fprintf(np, "%d %f %f %d\n", nc, geom->cellx[c], geom->celly[c], mk);
	fprintf(bp, "%f\n", geom->botz[c]);
      }
      if (mode & 2) geom->w2c[nc] = c;
      geom->c2w[c] = nc++;
    }
  }

  /* Print the triangulation of boundary cells                       */
  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    for (cc = 1; cc <= open->no2_t; cc++) {
      c = open->obc_t[cc];
      if (mask[c]) {
	if (mode & 1) {
	  fprintf(np, "%d %f %f %d\n", nc, geom->cellx[c], geom->celly[c], open->id+2);
	  fprintf(bp, "%f\n", geom->botz[c]);
	}
	if (mode & 2) geom->w2c[nc] = c;
	geom->c2w[c] = nc++;
      }
    }
  }

  /* Print the triangle connectivity                                 */
  if (mode & 1) {
    fclose(np);
    fclose(bp);
    sprintf(buf, "%s.ele", key);
    ep = fopen(buf, "w");
    fprintf(ep, "%d 3 0\n", nt);
    if (checkf) {
      sprintf(buf, "%s_swan.txt", key);
      op = fopen(buf, "w");
    }
  }
  nc = 1;
  for (vv = 1; vv <= geom->v2_e2; vv++) {
    v = geom->w2_e2[vv];
    if (geom->nvc[v] != 3) 
      hd_warn("swan_mesh(): found more than 3 v2c maps at v=%d[%f %f]\n", v,
	      geom->gridx[v], geom->gridy[v]);
    /* Check for valid trianges                                      */
    c1 = 0;
    for (cc = 1; cc <= geom->nvc[v]; cc++) {
      c = geom->v2c[v][cc];
      if (!c || geom->nvc[v] != 3) c1 = 1;
    }
    if (c1) continue;

    /* Print the index of triange nodes to file                      */
    if (mode & 1) fprintf(ep, "%d ", nc++);
    for (cc = 1; cc <= geom->nvc[v]; cc++) {
      c = geom->v2c[v][cc];
      if (c == 0) hd_warn("swan_mesh(): zero centre found at v=%d, cc=%d\n", v, cc);
      if (mode & 1) {
	fprintf(ep, "%d ", geom->c2w[c]);
	if (checkf) fprintf(op, "%f %f\n", geom->cellx[c], geom->celly[c]);
      }
    }
    if (mode & 1) fprintf(ep, "\n");

    /* Check if required                                             */
    if (checkf) {
      c = geom->v2c[v][1];
      if (mode & 1) {
	fprintf(op, "%f %f\n", geom->cellx[c], geom->celly[c]);
	fprintf(op, "NaN NaN\n");
      }
    }
  }
  i_free_1d(mask);

  if (mode & 1) {
    fclose(ep);
    /* Check: Read the .ele and .node files back in, and dump the    */
    /* triangulation.                                                */
    if (checkf) {
      int v1, v2, v3, *vmark;
      double *x, *y;
      char *fields[MAXSTRLEN * MAXNUMARGS];
      fclose(op);

      sprintf(buf, "%s.node", key);
      if ((np = fopen(buf, "r")) == NULL)
	hd_quit("Can't open file %s\n", buf);
      c = 0;
      while (fgets(buf, 256, np) != NULL) {
	cc = parseline(buf, fields, MAXNUMARGS);
	if (cc != 4) 
	  hd_warn("write_swan_mesh(): Error at line %d: %s. Not enough entries (%d).\n",c, cc);
	if (c == 0) {
	  n = atoi(fields[0]);
	  x = d_alloc_1d(n+1);
	  y = d_alloc_1d(n+1);
	  vmark = i_alloc_1d(n+1);
	} else {
	  if (c != atoi(fields[0])) 
	    hd_warn("write_swan_mesh(): .node entries not sequential at %d\n",c);
	  x[c] = atof(fields[1]);
	  y[c] = atof(fields[2]);
	  vmark[c] = atoi(fields[3]);
	}
	c++;
      }

      if (c-1 != n) 
	hd_warn("write_swn_mesh(): Incorrect number of entries in .node: %d vs %d\n",n, c);
      fclose(np);

      sprintf(buf, "%s.ele", key);
      if ((ep = fopen(buf, "r")) == NULL)
	hd_quit("Can't open file %s\n", buf);

      sprintf(buf, "%s_swan.txt", key);
      op = fopen(buf, "w");
      c = 0;
      while (fgets(buf, 256, ep) != NULL) {
	cc = parseline(buf, fields, MAXNUMARGS);
	if (c == 0) {
	  n = atoi(fields[0]);
	} else {
	  if (cc != 4) hd_quit("write_swan_mesh(): Error at line %d: triangles != 3 (%d)\n",c, cc);
	  if (c != atoi(fields[0])) 
	    hd_quit("write_swan_mesh(): .ele entries not sequential at %d\n",c);
	  v1 = atoi(fields[1]);
	  v2 = atoi(fields[2]);
	  v3 = atoi(fields[3]);
	  /* Interior cells only
	     if (vmark[v1] == 0 && vmark[v2] == 0 && vmark[v3] == 0) {*/
	  fprintf(op, "%f %f\n",x[v1], y[v1]);
	  fprintf(op, "%f %f\n",x[v2], y[v2]);
	  fprintf(op, "%f %f\n",x[v3], y[v3]);
	  fprintf(op, "%f %f\n",x[v1], y[v1]);
	  fprintf(op, "NaN NaN\n");
	}
	c++;
      }
      if (c-1 != n) 
	hd_warn("write_swn_mesh(): Incorrect number of entries in .ele: %d vs %d\n",n, c);
      d_free_1d(x);
      d_free_1d(y);
      i_free_1d(vmark);
      fclose(ep);
      fclose(op);
    }
  }

  if (mode & 2) {
    /* Get the mapping from the sparse vector to SWAN nodes          */
    mask = i_alloc_1d(geom->szcS);
    memset(mask, 0, geom->szcS * sizeof(int));
    for (cc = 1; cc <= geom->nw2c; cc++) {
      cg = geom->w2c[cc];
      mask[cg] = geom->fm[cg].wn;
    }
    for (n = 1; n <= geom->nwindows; n++) {
      window[n]->nw2c = 0;
      for (cc = 1; cc <= window[n]->b2_t; cc++) {
	c = window[n]->w2_t[cc];
	cg = window[n]->wsa[c];
	if (mask[cg] == window[n]->wn) window[n]->nw2c++;
      }
      window[n]->w2c = i_alloc_1d(window[n]->nw2c+1);
      nc = 1;
      for (cc = 1; cc <= window[n]->b2_t; cc++) {
	c = window[n]->w2_t[cc];
	cg = window[n]->wsa[c];
	if (mask[cg] == window[n]->wn) window[n]->w2c[nc++] = c;
      }
    }
    i_free_1d(mask);

    /* Write a SWAN parameter file                                   */
    mesh_ofile_name(params, key);
    sprintf(buf, "%s.swn", key);
    write_swan_params(master, buf, key);
    if (params->runmode & (AUTO | DUMP))
      write_swan_params(master, "INPUT", key);
  }
}

/* END write_swan_mesh()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Writes a SWAN parameter file for coupling with COMPAS             */
/*-------------------------------------------------------------------*/
void write_swan_params(master_t *master, char *fname, char *mname)
{
  geometry_t *geom = master->geom;
  parameters_t *params = master->params;
  FILE *op;
  char buf[MAXSTRLEN], key[MAXSTRLEN], dt[MAXSTRLEN];
  char restart[MAXSTRLEN];
  int n, i, bf = 0, rs = 0;
  double d1;
  int y, mo, d, h, mi, s;
  char *sdate, *edate;
  char ed[MAXSTRLEN], sd[MAXSTRLEN];

  op = fopen(fname, "w");

  /* Get the restart name                                            */
  if (strlen(params->restart_name)) {
    rs = 1;
    if (endswith(params->restart_name, ".nc")) {
      n = strlen(params->restart_name);
      for (i = 0; i < n-3; i++)
	key[i] = params->restart_name[i];
      sprintf(restart, "%s_swan", key);
    } else
      sprintf(restart, "%s_swan", params->restart_name);
  }

  fprintf(op, "$------------------------\n");
  fprintf(op, "$ %s input file for UNSWAN\n", params->grid_name);
  fprintf(op, "$------------------------\n\n");

  fprintf(op, "$ Start-up\n");
  strcpy(buf, "001");
  if (strlen(params->runnoc)) strcpy(buf, params->runnoc);
  /* PROJECT name must be < 16 characters. SWAN does not operate     */
  /* properly if it isn't.                                           */
  if (strlen(params->grid_name) > 16) {
    for (i = 0; i < 16; i++)
      key[i] = params->grid_name[i];
    key[i] = '\0';
    hd_warn("SWAN PROJECT name > 16 characters. Changing to '%s'.\n", key);
  } else 
    strcpy(key, params->grid_name);
  fprintf(op, "PROJECT '%s' '%s'\n", key, buf);
  /*fprintf(op, "SET 0 90 0.5 200 1 9.81 1025 inrhog=1 CARTesian\n");*/
  fprintf(op, "SET 0 90 0.5 200 1 9.81 1025 inrhog=0 CARTesian\n");
  fprintf(op, "MODE NONSTATIONARY TWODIMENSIONAL\n");
  fprintf(op, "COORDINATES SPHERICAL\n\n");

  /*fprintf(op, "CGRID UNSTRUCTURED CIRCLE 24 0.05 1.0 24\n");*/
  fprintf(op, "CGRID UNSTRUCTURED CIRCLE 36 0.0412 1. 35 \n");
  fprintf(op, "READgrid UNSTRUCTURED TRIANGLE '%s'\n\n", mname);
 
  fprintf(op, "INPgrid  BOTTOM UNSTRUCtured\n");
  fprintf(op, "READinp  COMBOT\n\n");

  fprintf(op, "INPgrid  WLEVel UNSTRUCtured\n");
  fprintf(op, "READinp  COMWL\n\n");

  fprintf(op, "INPgrid  CURrent UNSTRUCtured\n");
  fprintf(op, "READinp  COMCUR\n\n");

  fprintf(op, "INPgrid  WInd UNSTRUCtured\n");
  fprintf(op, "READinp  COMWIND\n\n");

  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    if (open->options & OP_WAVES) {
      bf = 1;
      break;
    }
  }
  if (bf) {
    fprintf(op, "BOU SHAPE JON PEAK DSPR POWER\n");
    fprintf(op, "BOUNDSPEC COMBOUND\n\n");
  }

  fprintf(op, "$ Wave model processes and parameter\n");
  if (forced_restart && rs) {
    fprintf(op, "INITial HOTStart '%s'\n", restart);
  } else {
    if (rs) fprintf(op, "$ INITial HOTStart '%s'\n", restart);
    fprintf(op, "INITIAL DEFAULT\n");
  }
  fprintf(op, "GEN3 AGROW\n");
  fprintf(op, "FRIC\n");
  fprintf(op, "BREAKING\n");
  fprintf(op, "TRIAD\n");
  fprintf(op, "PROP BSBT\n");
  fprintf(op, "NUM ACCUR 0.02 0.02 0.02 95 NONSTAT 7\n\n");

  if (rs) {
    fprintf(op, "RESTart '%s'\n\n", restart);
  }


  sdate = tm_time_to_datestr(schedule->start_time, master->timeunit);
  strcpy(sd, sdate);
  edate = tm_time_to_datestr(schedule->stop_time, master->timeunit);
  strcpy(ed, edate);
  /* Convert to ISO format for SWAN                                  */
  sscanf(sd, "%4d-%02d-%02d %02d:%02d:%02d", &y, &mo, &d, &h, &mi, &s);
  sprintf(sd, "%4d%02d%02d.%02d%02d%02d", y, mo, d, h, mi, s);
  sscanf(ed, "%4d-%02d-%02d %02d:%02d:%02d", &y, &mo, &d, &h, &mi, &s);
  sprintf(ed, "%4d%02d%02d.%02d%02d%02d", y, mo, d, h, mi, s);

  d1 = params->wavedt;
  if (d1 >= 86400.0)
    sprintf(dt, "%5.1f DAY", d1 / 86400.0);
  else if (d1 >= 3600.0)
    sprintf(dt, "%5.1f HR", d1 / 3600.0);
  else if (d1 >= 60.0)
    sprintf(dt, "%5.1f MIN", d1 / 60.0);
  else
    sprintf(dt, "%5.1f SEC", d1);
  fprintf(op, "BLOCK  'COMPGRID' HEADER 'noncon.nc' HSIGN GENWIND DISBOT DISSURF DISWCAP DISTUR OUTput %s %s \n",sd, dt);
  fprintf(op, "$ Run simulation\n");
  fprintf(op, "COMPUTE NONSTAT %s %s %s\n", sd, dt, ed);

  /*
  if (rs) {
    fprintf(op, "\nHOTFile '%s'\n", restart);
  }
  */

  fflush(op);
  fclose(op);
}

/* END write_swan_params()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to set wave variables on all open boundaries.             */
/*-------------------------------------------------------------------*/
void bdry_wave(master_t *master, /* Master data                      */
	       geometry_t *geom, /* Master geometry                  */
	       swan_data_t *data
	       )
{
  int n, cc, c, cw;              /* Counters */
  double x, y, z;

  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    double x, y, z;

    if (open->bcond_wav & NOTHIN) continue;

    for (cc = 1; cc <= open->no2_t; cc++) {
      c = open->obc_t[cc];
      cw = geom->c2w[c];
      
      if (open->bcond_wav & FILEIN) {
	x = geom->cellx[c];
	y = geom->celly[c];
	z = geom->cellz[c] * master->Ds[c];
	data->swan_amp[cw] = 
	  hd_ts_multifile_eval_xyz_by_name(open->ntsfiles, 
					   open->tsfiles,
					   open->filenames, 
					   "wave_amp",
					   master->t3d,
					   x, y, z);
	data->swan_per[cw] = 
	  hd_ts_multifile_eval_xyz_by_name(open->ntsfiles, 
					   open->tsfiles,
					   open->filenames, 
					   "wave_period",
					   master->t3d,
					   x, y, z);
	/*if(master->t3d==0.0)data->swan_per[cw] = 2.0;*/
	data->swan_dir[cw] = 
	  hd_ts_multifile_eval_xyz_by_name(open->ntsfiles, 
					   open->tsfiles,
					   open->filenames, 
					   "wave_dir",
					   master->t3d,
					   x, y, z);
      }
      if (open->bcond_wav & CLAMPD) {
	data->swan_amp[cw] = 0.0;
	data->swan_per[cw] = 0.0;
	data->swan_dir[cw] = 0.0;
      }
      if (open->bcond_wav & NOGRAD) {
	int c1 = geom->c2w[open->oi1_t[cc]];
	data->swan_amp[cw] = data->swan_amp[c1];
	data->swan_per[cw] = data->swan_per[c1];
	data->swan_dir[cw] = data->swan_dir[c1];
      }
      if (open->bcond_wav & (CYCLIC|CYCLED)) {
	int c1 = geom->c2w[open->cyc_t[cc]];
	data->swan_amp[cw] = data->swan_amp[c1];
	data->swan_per[cw] = data->swan_per[c1];
	data->swan_dir[cw] = data->swan_dir[c1];
      }
      if (open->bcond_wav & LINEXT) {
	int *imap = open->nmap;
	int c1 = open->oi1_t[cc];
	int c2 = open->oi2_t[cc];
	int c3 = imap[c2];
	int c4 = imap[c3];
	data->swan_amp[cw] = bc_leastsq(data->swan_amp[c4], 
				       data->swan_amp[c3], 
				       data->swan_amp[c2], 
				       data->swan_amp[c1]);
	data->swan_per[cw] = bc_leastsq(data->swan_per[c4],
				       data->swan_per[c3],
				       data->swan_per[c2],
				       data->swan_per[c1]);
	data->swan_dir[cw] = bc_leastsq(data->swan_dir[c4], 
				       data->swan_dir[c3], 
				       data->swan_dir[c2], 
				       data->swan_dir[c1]);
      }
      if (open->bcond_wav & POLEXT) {
	int *imap = open->nmap;
	int c1 = open->oi1_t[cc];
	int c2 = open->oi2_t[cc];
	int c3 = imap[c2];
	data->swan_amp[cw] = bc_polint(data->swan_amp[c3], 
				      data->swan_amp[c2], 
				      data->swan_amp[c1]);
	data->swan_per[cw] = bc_polint(data->swan_per[c3],
				      data->swan_per[c2],
				      data->swan_per[c1]);
	data->swan_dir[cw] = bc_polint(data->swan_dir[c3], 
				      data->swan_dir[c2], 
				      data->swan_dir[c1]);
      }
      /* Wave direction is input to SWAN as degrees but saved        */
      /* internally in radians.                                      */
      if (!(open->bcond_wav & FILEIN))
	data->swan_dir[cw] *= 180.0 / PI;

    }
  }
}

/* END bdry_wave()                                                   */
/*-------------------------------------------------------------------*/
