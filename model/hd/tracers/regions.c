/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/tracers/regions.c
 *  
 *  Description:
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: regions.c 5992 2018-10-17 04:54:39Z her127 $
 *
 */

#include <stdlib.h>
#include <string.h>
#include "hd.h"
#include "tracer.h"

#define ALL_TRACERS       1
#define TRACERS_WC        2
#define TRACERS_DIAGN_WC  4
#define MAXREGIONS        1000

region_t *region_alloc(void);
int read_region(master_t *master, parameters_t *params);
void get_regions(geometry_t *win, region_t *region, double *mask, int nregions);
void region_schedule(region_t *region, double dt, double trem, char *timeunit);
double get_totflux(region_t *region, int tr);
int find_trindex(region_t *region, int trn);

/*-------------------------------------------------------------------*/
/* Reads the region infromation from the parameter file              */
/*-------------------------------------------------------------------*/
void read_region_info(parameters_t *params, FILE *fp)
{

  if (prm_read_char(fp, "REGION", params->regions)) {
    strcpy(params->region_dt, params->stop_time);
    if (!prm_read_char(fp, "REGION_DT", params->region_dt))
      sprintf(params->region_dt, "%f days", get_run_length(fp) + 1);
    prm_read_char(fp, "REGION_VARS", params->region_vars);
    prm_read_char(fp, "REGION_MODE", params->region_mode);
    prm_read_int(fp, "REGION_OBC_ZONE", &params->region_obc);
    params->ntr += 2;
  }
}

/* END read_region_info()                                            */ 
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets up the global regions data structure                         */
/*-------------------------------------------------------------------*/
void init_regions_g(master_t *master, parameters_t *params)
{
  geometry_t *geom = master->geom;
  int nregions;
  int n, m, ntr, rn, c, cc, cs, c1;
  region_t *region;
  char *fields[MAXSTRLEN * MAXNUMARGS];
  int varf = 0, *nobc;
  int **rmap;

  if (!strlen(params->regions)) {
    geom->nregions = 0;
    return;
  }

  /*-----------------------------------------------------------------*/
  /* Read the region mask from file                                  */
  geom->nregions = nregions = read_region(master, params);
  if (nregions == 0) return;

  /*-----------------------------------------------------------------*/
  /* Get any region attributes                                       */
  if (strlen(params->region_mode)) {
    char buf[MAXSTRLEN];
    strcpy(buf, params->region_mode);
    if (contains_token(buf, "NONE") != NULL) {
      master->region_mode = NONE;
    } else {
      master->region_mode = 0;
      if (contains_token(buf, "OBC_BDRY") != NULL)
	master->region_mode |= RG_BDRY;
      if (contains_token(buf, "OBC_AREA") != NULL)
	master->region_mode |= RG_AREA;
      if (contains_token(buf, "ALL_TRANSFERS") != NULL)
	master->region_mode |= RG_ALLT;
    }
  }
  master->region_obcz = params->region_obc;
  if (params->runmode & TRANS && params->trasc != FFSL) {
    if (!(master->region_mode & RG_AREA)) {
      if (master->region_mode & NONE)
	master->region_mode = RG_AREA;
      else
	master->region_mode |= RG_AREA;
      master->region_obcz = 3;
    } else {
      if (master->region_obcz < 2) {
	master->region_obcz = 2;
	hd_warn("regions_init: minimum OBC_AREA size = 3\n");
      }
    }
  } else {
    if (!(master->region_mode & RG_BDRY) && master->region_obcz == 0)
      hd_warn("regions_init: OBC fluxes non-conservative - use OBC_BDRY with REGION_OBC_ZONE >= 1\n");
  }

  /*-----------------------------------------------------------------*/
  /* Get the region interval                                         */
  if (strlen(params->region_dt)) {
    if (tm_scale_to_secs(params->region_dt, &master->region_dt)) {
      master->region_dt = (master->region_dt < master->dt) ? master->dt :
	master->region_dt;
      master->region_next = master->t + master->region_dt;
    } else if (contains_token(params->region_dt, "YEARLY")) {
      master->region_dt = YEARLY;
      master->region_next = master->t + next_year(master->t, master->timeunit);
    } else if (contains_token(params->region_dt, "SEASONAL")) {
      master->region_dt = SEASONAL;
      master->region_next = master->t + next_season(master->t, 
						    master->timeunit, &c);
    } else if (contains_token(params->region_dt, "MONTHLY")) {
      master->region_dt = MONTHLY;
      master->region_next = master->t + next_month(master->t, 
						   master->timeunit, &c);
    } else {
      hd_warn("Incorrect REGION_DT specified.\n");
      geom->nregions = 0;
      return;
    }
  } else {
    hd_warn("No REGION_DT specified.\n");
    geom->nregions = 0;
    return;
  }

  /*-----------------------------------------------------------------*/
  /* Get the region variables                                        */
  master->region_nvar = 3;     /* volume, temp, salt always included */
  master->region_nvarS = 0;
  if (strlen(params->region_vars)) {
    n = parseline(params->region_vars, fields, MAXNUMARGS);
    if (n) {
      if (strcmp(fields[0], "ALL") == 0) {
	master->region_nvar = master->ntr + 1;
	varf = ALL_TRACERS;
      } else if (strcmp(fields[0], "TRACERS_WC") == 0) {
	master->region_nvar = 1;
	for (ntr = 0; ntr < master->ntr; ntr++)
	  if (!master->trinfo_3d[ntr].diagn) master->region_nvar++;
	varf = TRACERS_WC;
      } else if (strcmp(fields[0], "TRACERS_DIAGN_WC") == 0) {
	for (ntr = 0; ntr < master->ntr; ntr++)
	  if (master->trinfo_3d[ntr].diagn) master->region_nvar++;
	varf = TRACERS_DIAGN_WC;
      } else {
	for (rn = 0; rn < n; rn++) {
	  if ((tracer_find_index(fields[rn], params->ntr, params->trinfo_3d)) >= 0) {
	    if ((strcmp(fields[rn], "temp") != 0) && (strcmp(fields[rn], "salt") != 0))
	      master->region_nvar++;
	  }
	  if ((tracer_find_index(fields[rn], params->ntrS, params->trinfo_2d)) >= 0) {
	    master->region_nvarS++;
	  }
	}
      }
    }
  } else
    strcpy(params->region_vars, "temp salt");

  master->region_var = i_alloc_1d(master->region_nvar);
  if (master->region_nvarS)
    master->region_varS = i_alloc_1d(master->region_nvarS);
  /* Index 0 = volume                                                */
  if (varf & ALL_TRACERS)
    for (ntr = 0; ntr < master->ntr; ntr++) {
      master->region_var[ntr] = ntr;
    }
  else if (varf & TRACERS_WC) {
    m = 1;
    for (ntr = 0; ntr < master->ntr; ntr++) {
      if (!master->trinfo_3d[ntr].diagn) {
	master->region_var[m] = ntr;
	m++;
      }
    }
  }
  else if (varf & TRACERS_DIAGN_WC) {
    m = 3;
    for (ntr = 0; ntr < master->ntr; ntr++) {
      if (master->trinfo_3d[ntr].diagn) {
	master->region_var[m] = ntr;
	m++;
      }
    }
  } else {
    master->region_var[1] = tracer_find_index("temp", params->ntr, params->trinfo_3d);
    master->region_var[2] = tracer_find_index("salt", params->ntr, params->trinfo_3d);
    m = 3;
    for (c = 0; c < n; c++) {
      ntr = tracer_find_index(fields[c], params->ntr, params->trinfo_3d);
      if (ntr >= 0) {
	if ((strcmp(fields[c], "temp") != 0) && (strcmp(fields[c], "salt") != 0)) {
	  master->region_var[m] = ntr;
	  m++;
	}
      }
    }
  }
  m = 0;
  for (c = 0; c < n; c++) {
    ntr = tracer_find_index(fields[c], params->ntrS, params->trinfo_2d);
    if (ntr >= 0) {
      master->region_varS[m] = ntr;
      m++;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set up the open boundary mask = 1 master->region_obcz cells     */
  /* into the interior.                                              */
  for (cc = 1; cc <= geom->nbpt; cc++) {
    c = geom->bpt[cc];
    master->regionid[c] = NOTVALID;
  }
  if (geom->nobc) {
    nobc = i_alloc_1d(geom->nobc);
    memset(nobc, 0, geom->nobc * sizeof(int));
  }
  for (m = 0; m < geom->nobc; m++) {
    int bb;
    open_bdrys_t *open = geom->open[m];
    for (cc = 1; cc <= open->no3_t; cc++) {
      c = open->obc_t[cc];
      if (open->stagger & OUTFACE && open->ocodex & R_EDGE && c != geom->xp1[c])
	master->regionid[geom->xp1[c]] = NOTVALID;
      if (open->stagger & OUTFACE && open->ocodex & L_EDGE && c != geom->xm1[c])
	master->regionid[geom->xm1[c]] = NOTVALID;
      if (open->stagger & OUTFACE && open->ocodey & F_EDGE && c != geom->yp1[c])
	master->regionid[geom->yp1[c]] = NOTVALID;
      if (open->stagger & OUTFACE && open->ocodey & B_EDGE && c != geom->ym1[c])
	master->regionid[geom->ym1[c]] = NOTVALID;
      if (master->regionid[c] != NOTVALID) {
	if (!nobc[m] && master->region_mode & RG_AREA) {
	  nobc[m] = 1;
	  nregions++;
	}
	for (bb = 0; bb < master->region_obcz; bb++) {
	  master->regionid[c] = (master->region_mode & RG_AREA) ? nregions-1 : NOTVALID;
	  c = open->nmap[c];
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Populate the regions structures                                 */
  geom->nregions = nregions;
  geom->region = (region_t **)malloc(sizeof(region_t *) * (geom->nregions));
  for (n = 0; n < geom->nregions; n++) {
    geom->region[n] = region_alloc();
    region = geom->region[n];
    region->nregions = geom->nregions;
    region->t = master->t;
    region->dt = master->region_dt;
    region->next = master->region_next;
    region->nfaces = nregions + geom->nobc;
    region->id = n;
    region->mode = master->region_mode;
    region->obcz = master->region_obcz;
    sprintf(region->name, "Region%d", n);

    /* Set the output variables */
    if (master->Vi) region->mode |= RG_VERR;
    if (master->fillf & MONOTONIC) region->mode |= RG_GLOB;
    if (master->runmode & TRANS && master->npss) region->mode |= RG_PSS;
#ifdef HAVE_SEDIMENT_MODULE
    if (master->do_sed) region->mode |= RG_SED;
#endif
#ifdef HAVE_ECOLOGY_MODULE
    if (master->do_eco) region->mode |= RG_ECO;
#endif

    /* Save the variables                                            */
    region->nvar = master->region_nvar;
    region->var = i_alloc_1d(master->region_nvar);
    region->nvarS = master->region_nvarS;
    if (master->region_nvarS)
      region->varS = i_alloc_1d(master->region_nvarS);
    for (m = 0; m < region->nvar; m++)
      region->var[m] = master->region_var[m];
    for (m = 0; m < region->nvarS; m++)
      region->varS[m] = master->region_varS[m];

    /* Get the region locations                                      */
    get_regions(geom, region, master->regionid, geom->nregions);

    /* Allocate the mass arrays                                      */
    region->mmass = d_alloc_1d(region->nvar);
    region->stdev = d_alloc_1d(region->nvar);
    region->smass = d_alloc_1d(region->nvar);
    region->emass = d_alloc_1d(region->nvar);
    if (region->mode & RG_PSS) region->pssf = d_alloc_1d(region->nvar);
    if (region->mode & RG_VERR) region->merr = d_alloc_1d(region->nvar);
    if (region->mode & RG_GLOB) region->gfer = d_alloc_1d(region->nvar);
    if (region->mode & RG_SED) region->sedtr = d_alloc_1d(region->nvar);
    if (region->mode & RG_ECO) region->eco = d_alloc_1d(region->nvar);
    if (region->nboundaries) {
      region->trflux = d_alloc_2d(region->nboundaries, region->nvar);
      region->w1 = d_alloc_2d(region->nboundaries, region->nvar);
    }
    memset(region->smass, 0, region->nvar * sizeof(double));
    memset(region->emass, 0, region->nvar * sizeof(double));
    memset(region->mmass, 0, region->nvar * sizeof(double));
    memset(region->stdev, 0, region->nvar * sizeof(double));
    if (region->pssf) memset(region->pssf, 0, region->nvar * sizeof(double));
    if (region->merr) memset(region->merr, 0, region->nvar * sizeof(double));
    if (region->nboundaries) {
      for (m = 0; m < region->nvar; m++) {
	memset(region->trflux[m], 0, region->nboundaries * sizeof(double));
	memset(region->w1[m], 0, region->nboundaries * sizeof(double));
      }
    }
    region->transferf = RG_BEGIN;
    region->pdt = 0.0;
    region->vfluxp = region->vfluxn = 0.0;
  }

  /* Find the nearest boundary in the reverse maps                 */
  if (!(region->mode & RG_ALLT)) {
    rmap = i_alloc_2d(nregions, nregions);
    for (n = 0; n < nregions; n++) {
      region = geom->region[n];
      for (m = 0; m < region->nregions; m++) {
	rmap[n][m] = region->rmap[m];
	if (m != n && region->rmap[m] == NOTVALID) {
	  region_t *reg = geom->region[m];
	  int n1, n2;
	  for (n1 = 0; n1 < reg->nregions; n1++) {
	    for (n2 = 0; n2 < region->nregions; n2++) {
	      if (region->rmap[n2] != NOTVALID && reg->rmap[n1] != NOTVALID && n1 == n2)
		rmap[n][m] = region->rmap[n2];
	    }
	  }
	}
      }
    }
    for (n = 0; n < nregions; n++) {
      region = geom->region[n];
      for (m = 0; m < region->nregions; m++) {
	region->rmap[m] = rmap[n][m];
      }
    }
    i_free_2d(rmap);
  }

  /* Set the region mask for this window                           */
  /* Set ghost points                                              */
  for (cc = 1; cc <= geom->nbpt; cc++) {
    c = geom->bpt[cc];
    c1 = geom->bin[cc];
    master->regionid[c] = master->regionid[c1];
  }
  /* Set sediment points                                           */
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->bot_t[cc];
    master->regionid[geom->zm1[c]] = master->regionid[c];
  }
  /* Ghost cells for sediments                                     */
  for (cc = 1; cc <= geom->ngsed; cc++) {
    c = geom->gsed_t[cc];
    c1 = geom->ised_t[cc];
    master->regionid[c] = master->regionid[c1];
  }
}

/* END init_regions_g()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets up the regions data structure on the windows                 */
/*-------------------------------------------------------------------*/
void init_regions_w(master_t *master, geometry_t **window)
{
  geometry_t *geom = master->geom;
  region_t *region;
  int n, m, nr, wn, wr;
  int c, cc, c1, bb;
  int nregions, *nobc, **rmap;

  if (geom->nregions == 0)
    return;
  if (geom->nwindows == 1) {
    window[1]->nregions = geom->nregions;
    window[1]->region = geom->region;
    return;
  }

  /*-----------------------------------------------------------------*/
  /* Loop through the windows and find regions                       */
  for (wn = 1; wn <= geom->nwindows; wn++) {
    window_t *windat = window[wn]->windat;
    win_priv_t *wincon = window[wn]->wincon;
    int wr[geom->nregions];
    double *mask = windat->regionid;

    for (m = 0; m < geom->nregions; m++) wr[m] = NOTVALID;
    for (c = 1; c <= window[wn]->enon; c++)
      mask[c] = NOTVALID;

    /* Get the number of regions in this window                      */
    for (cc = 1; cc <= window[wn]->b3_t; cc++) {
      c = window[wn]->w3_t[cc];
      c1 = window[wn]->wsa[c];
      mask[c] = master->regionid[c1];
      nr = (int)mask[c];
      if (nr >= 0 && nr < geom->nregions) wr[nr] = nr;
    }
    for (cc = 1; cc <= window[wn]->a3_e1; cc++) {
      c = window[wn]->w3_e1[cc];
      c1 = window[wn]->wsa[c];
      if (mask[c] == NOTVALID && master->regionid[c1] != NOTVALID) {
	mask[c] = master->regionid[c1];
	nr = (int)mask[c];
      }
    }
    for (cc = 1; cc <= window[wn]->a3_e2; cc++) {
      c = window[wn]->w3_e2[cc];
      c1 = window[wn]->wsa[c];
      if (mask[c] == NOTVALID && master->regionid[c1] != NOTVALID) {
	mask[c] = master->regionid[c1];
	nr = (int)mask[c];
      }
    }

    /* Set the region mask for this window                           */
    /* Set ghost points                                              */
    for (cc = 1; cc <= window[wn]->nbpt; cc++) {
      c = window[wn]->bpt[cc];
      c1 = window[wn]->bin[cc];
      mask[c] = mask[c1];
    }

    /* Set sediment points                                           */
    for (cc = 1; cc <= window[wn]->b2_t; cc++) {
      c = window[wn]->bot_t[cc];
      mask[window[wn]->zm1[c]] = mask[c];
    }

    /* Ghost cells for sediments                                     */
    /*
    for (cc = 1; cc <= window[wn]->ngsed; cc++) {
      c = window[wn]->gsed_t[cc];
      c1 = window[wn]->ised_t[cc];
      mask[c] = mask[c1];
    }
    */

    /* Count the number of regions                                   */
    nregions = 0;
    for (nr = 0; nr < geom->nregions; nr++) if (wr[nr] != NOTVALID) nregions++;

    /* Open boundaries                                               */
    if (window[wn]->nobc) {
      nobc = i_alloc_1d(window[wn]->nobc);
      memset(nobc, 0, window[wn]->nobc * sizeof(int));
    }
    for (m = 0; m < window[wn]->nobc; m++) {
      open_bdrys_t *open = window[wn]->open[m];
      /* Set up the open boundary mask = 1 master->region_obcz cells */
      /* into the interior.                                          */
      for (cc = 1; cc <= open->no3_t; cc++) {
	c = open->obc_t[cc];
	if (open->ocodex & R_EDGE && open->stagger & OUTFACE && c != window[wn]->xp1[c])
	  mask[window[wn]->xp1[c]] = NOTVALID;
	if (open->ocodey & F_EDGE && open->stagger & OUTFACE && c != window[wn]->yp1[c])
	  mask[window[wn]->yp1[c]] = NOTVALID;
	if (mask[c] != NOTVALID) {
	  if (!nobc[m] && master->region_mode & RG_AREA) {
	    nobc[m] = 1;
	    nregions++;
	  }
	  for (bb = 0; bb < master->region_obcz; bb++) {
	    mask[c] = (master->region_mode & RG_AREA) ? nregions - 1 : NOTVALID;
	    c = open->nmap[c];
	  }
	}
      }
    }

    /*---------------------------------------------------------------*/
    /* Allocate the region structure                                 */
    nr = 0;
    for (m = 0; m < window[wn]->nobc; m++) nr += nobc[m];
    if (nobc) i_free_1d(nobc);
    if (master->region_obcz && master->region_mode & RG_AREA) 
      nregions += nr;

    window[wn]->region = (region_t **)malloc(sizeof(region_t *) * (nregions));    
    /* Populate the region structure                                 */
    for (nr = 0; nr < nregions; nr++) {
      window[wn]->nregions = nregions;
      window[wn]->region[nr] = region_alloc();
      region = window[wn]->region[nr];
      region->nregions = nregions;
      region->t = master->t;
      region->dt = master->region_dt;
      region->next = master->region_next;
      region->nfaces = geom->nregions + window[wn]->nobc;
      region->obcz = master->region_obcz;
      region->id = NOTVALID;
      for (m = 0; m < geom->nregions; m++) {
	if (wr[m] != NOTVALID) {
	  region->id = n = wr[m];
	  wr[m] = NOTVALID;
	  break;
	}
      }
      region->mode = geom->region[n]->mode;

      /* Set the output variables */
      if (windat->Vi) region->mode |= RG_VERR;
      if (wincon->fillf & MONOTONIC) region->mode |= RG_GLOB;
      if (master->runmode & TRANS && wincon->npss) region->mode |= RG_PSS;
#ifdef HAVE_SEDIMENT_MODULE
      if (wincon->do_sed) region->mode |= RG_SED;
#endif
#ifdef HAVE_ECOLOGY_MODULE
      if (wincon->do_eco) region->mode |= RG_ECO;
#endif
      /* Save the variables                                          */
      region->nvar = master->region_nvar;
      region->var = i_alloc_1d(master->region_nvar);
      region->nvarS = master->region_nvarS;
      if (master->region_nvarS)
	region->varS = i_alloc_1d(master->region_nvarS);
      for (m = 0; m < region->nvar; m++)
	region->var[m] = master->region_var[m];
      for (m = 0; m < region->nvarS; m++)
	region->varS[m] = master->region_varS[m];

      /* Populate the regions                                        */
      get_regions(window[wn], region, mask, geom->nregions);

      /* Get the map from local to global boundaries */
      if (region->gmap) i_free_1d(region->gmap);
      if (region->nboundaries)
	region->gmap = i_alloc_1d(region->nboundaries);
      for (m = 0; m < region->nboundaries; m++) {
	for (bb = 0; bb < geom->region[n]->nboundaries; bb++) {
	  if (strcmp(region->fluxname[m], geom->region[n]->fluxname[bb]) == 0) {
	    region->gmap[m] = bb;
	    break;
	  }
	}
      }

      /* Allocate the mass arrays                                    */
      region->mmass = d_alloc_1d(region->nvar);
      region->stdev = d_alloc_1d(region->nvar);
      region->smass = d_alloc_1d(region->nvar);
      region->emass = d_alloc_1d(region->nvar);
      if (region->mode & RG_PSS) region->pssf = d_alloc_1d(region->nvar);
      if (region->mode & RG_VERR) region->merr = d_alloc_1d(region->nvar);
      if (region->mode & RG_GLOB) region->gfer = d_alloc_1d(region->nvar);
      if (region->mode & RG_SED) region->sedtr = d_alloc_1d(region->nvar);
      if (region->mode & RG_ECO) region->eco = d_alloc_1d(region->nvar);
      if (region->nboundaries) {
	region->trflux = d_alloc_2d(region->nboundaries, region->nvar);
	region->w1 = d_alloc_2d(region->nboundaries, region->nvar);
      }
      region->pdt = 0.0;
      region->transferf = RG_BEGIN;
      region->vfluxp = region->vfluxn = 0.0;
    }
    /* Find the nearest boundary in the reverse maps                 */
    rmap = i_alloc_2d(nregions, nregions);
    for (nr = 0; nr < nregions; nr++) {
      region = window[wn]->region[nr];
      for (m = 0; m < region->nregions; m++) {
	rmap[nr][m] = region->rmap[m];
	if (m != nr && region->rmap[m] == NOTVALID) {
	  region_t *reg = window[wn]->region[m];
	  int n1, n2;
	  for (n1 = 0; n1 < reg->nregions; n1++) {
	    for (n2 = 0; n2 < region->nregions; n2++) {
	      if (region->rmap[n2] != NOTVALID && reg->rmap[n1] != NOTVALID && n1 == n2)
		rmap[nr][m] = region->rmap[n2];
	    }
	  }
	}
      }
    }
    for (nr = 0; nr < nregions; nr++) {
      region = window[wn]->region[nr];
      for (m = 0; m < region->nregions; m++) {
	region->rmap[m] = rmap[nr][m];
      }
    }
    i_free_2d(rmap);
  }
}

/* END init_regions_w()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to get the region locations                               */
/*-------------------------------------------------------------------*/
void get_regions(geometry_t *win, region_t *region, double *mask, int nregions)
{
  int c, cc, cr, bb, m, rn, hr;
  int n = region->id;
  int dotan = 1;

  /*-----------------------------------------------------------------*/
  /* Get the region locations                                        */
  region->nvec = 0;
  for (cc = 1; cc <= win->b3_t; cc++) {
    c = win->w3_t[cc];
    if (n == (int)mask[c])
      region->nvec++;
  }
  if (region->vec) i_free_1d(region->vec);
  if (region->nvec) region->vec = i_alloc_1d(region->nvec + 1);
  m = 1;
  for (cc = 1; cc <= win->b3_t; cc++) {
    c = win->w3_t[cc];
    if (n == (int)mask[c])
      region->vec[m++] = c;
  }

  /*-----------------------------------------------------------------*/
  /* Get the region boundaries in the e1 direction. The convention   */
  /* that a positive flux implies a flux into the region. Boundaries */
  /* refer to exchanges between regions and an exchange between a    */
  /* region and an OBC.                                              */
  if (region->nbe1) i_free_1d(region->nbe1);
  region->nbe1 = i_alloc_1d(region->nfaces);
  memset(region->nbe1, 0, region->nfaces * sizeof(int));
  /* Region to region boundaries                                     */
  for (cc = 1; cc <= region->nvec; cc++) {
    c = region->vec[cc];
    rn = (int)mask[win->xp1[c]];
    if (rn != NOTVALID && n != rn) {
      region->nbe1[rn]++;
    }
    rn = (int)mask[win->xm1[c]];
    if (rn != NOTVALID && n != rn) {
      region->nbe1[rn]++;
    }
  }
  /* Region to OBC boundaries                                        */
  for (m = 0; m < win->nobc; m++) {
    open_bdrys_t *open = win->open[m];
    region->nbe1[m + nregions] = 0;
    if (open->ocodex & (L_EDGE|R_EDGE)) {
      for (cc = 1; cc <= open->no3_t; cc++) {
	c = open->obc_t[cc];
	if (region->mode & RG_BDRY)
	  for (bb = 0; bb < region->obcz; bb++)
	    c = open->nmap[c];
	if (n == (int)mask[c])
	  region->nbe1[m + nregions]++;
      }
    }
    /* Tangential for RG_BDRY                                        */
    if (dotan && open->ocodey & (B_EDGE|F_EDGE)) {
      for (cc = 1; cc <= open->no3_t; cc++) {
	c = open->obc_t[cc];
	if (region->mode & RG_BDRY) {
	  for (bb = 0; bb < region->obcz; bb++) {
	    if (n == (int)mask[win->xp1[c]])
	      region->nbe1[m + nregions]++;
	    if (n == (int)mask[win->xm1[c]])
	      region->nbe1[m + nregions]++;
	    c = open->nmap[c];
	  }
	}
      }
    }
  }
  /* Get the maximum face size                                       */
  rn = 0;
  for (m = 0; m < region->nfaces; m++)
    if (region->nbe1[m] > rn) rn = region->nbe1[m];
  if (region->be1) i_free_2d(region->be1);
  if (region->sgne1) d_free_2d(region->sgne1);
  if (rn) {
    region->be1 = i_alloc_2d(rn, region->nfaces);
    region->sgne1 = d_alloc_2d(rn, region->nfaces);

    /* Populate the face arrays                                      */
    memset(region->nbe1, 0, region->nfaces * sizeof(int));
    for (cc = 1; cc <= region->nvec; cc++) {
      c = region->vec[cc];
      rn = (int)mask[win->xp1[c]];
      if (rn != NOTVALID && n != rn) {
	region->sgne1[rn][region->nbe1[rn]] = -1.0;
	region->be1[rn][region->nbe1[rn]++] = win->xp1[c];
      }
      rn = (int)mask[win->xm1[c]];
      if (rn != NOTVALID && n != rn) {
	region->sgne1[rn][region->nbe1[rn]] = 1.0;
	region->be1[rn][region->nbe1[rn]++] = c;
      }
    }
    for (m = 0; m < win->nobc; m++) {
      open_bdrys_t *open = win->open[m];
      region->nbe1[m + nregions] = 0;
      if (open->ocodex & (L_EDGE|R_EDGE)) {
	for (cc = 1; cc <= open->no3_t; cc++) {
	  c = open->obc_t[cc];
	  region->sgne1[m + nregions][region->nbe1[m + nregions]] = 1.0;
	  if (open->ocodex & R_EDGE) 
	    region->sgne1[m + nregions][region->nbe1[m + nregions]] = -1.0;
	  if (region->mode & RG_BDRY)
	    for (bb = 0; bb < region->obcz; bb++)
	      c = open->nmap[c];
	  cr = (open->ocodex & R_EDGE && open->stagger & OUTFACE ) ? win->xp1[c] : c;
	  if (n == (int)mask[c]) {
	    region->be1[m + nregions][region->nbe1[m + nregions]] = cr;
	    region->nbe1[m + nregions]++;
	  }
	}
      }
      /* Tangential for RG_BDRY                                      */
      if (dotan && open->ocodey & (B_EDGE|F_EDGE)) {
	for (cc = 1; cc <= open->no3_t; cc++) {
	  c = open->obc_t[cc];
	  if (region->mode & RG_BDRY) {
	    for (bb = 0; bb < region->obcz; bb++) {
	      cr = win->xp1[c];
	      if (n == (int)mask[cr]) {
		region->be1[m + nregions][region->nbe1[m + nregions]] = cr;
		region->sgne1[m + nregions][region->nbe1[m + nregions]] = 1.0;
		region->nbe1[m + nregions]++;
	      }
	      cr = win->xm1[c];
	      if (n == (int)mask[cr]) {
		region->be1[m + nregions][region->nbe1[m + nregions]] = c;
		region->sgne1[m + nregions][region->nbe1[m + nregions]] = -1.0;
		region->nbe1[m + nregions]++;
	      }
	      c = open->nmap[c];
	    }
	  }
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the region boundaries in the e2 direction                   */
  if (region->nbe2) i_free_1d(region->nbe2);
  region->nbe2 = i_alloc_1d(region->nfaces);
  memset(region->nbe2, 0, region->nfaces * sizeof(int));
  /* Region to region boundaries                                     */
  for (cc = 1; cc <= region->nvec; cc++) {
    c = region->vec[cc];
    rn = (int)mask[win->yp1[c]];
    if (rn != NOTVALID && n != rn) {
      region->nbe2[rn]++;
    }
    rn = (int)mask[win->ym1[c]];
    if (rn != NOTVALID && n != rn) {
      region->nbe2[rn]++;
    }
  }
  /* Region to OBC boundaries                                        */
  for (m = 0; m < win->nobc; m++) {
    open_bdrys_t *open = win->open[m];
    region->nbe2[m + nregions] = 0;
    if (open->ocodey & (B_EDGE|F_EDGE)) {
      for (cc = 1; cc <= open->no3_t; cc++) {
	c = open->obc_t[cc];
	if (region->mode & RG_BDRY)
	  for (bb = 0; bb < region->obcz; bb++)
	    c = open->nmap[c];
	if (n == (int)mask[c])
	  region->nbe2[m + nregions]++;
      }
    }
    /* Tangential for RG_BDRY                                        */
    if (dotan && open->ocodex & (L_EDGE|R_EDGE)) {
      for (cc = 1; cc <= open->no3_t; cc++) {
	c = open->obc_t[cc];
	if (region->mode & RG_BDRY) {
	  for (bb = 0; bb < region->obcz; bb++) {
	    if (n == (int)mask[win->yp1[c]])
	      region->nbe2[m + nregions]++;
	    if (n == (int)mask[win->ym1[c]])
	      region->nbe2[m + nregions]++;
	    c = open->nmap[c];
	  }
	}
      }
    }
  }
  /* Get the maximum face size                                       */
  rn = 0;
  for (m = 0; m < region->nfaces; m++)
    if (region->nbe2[m] > rn) rn = region->nbe2[m];
  if (region->be2) i_free_2d(region->be2);
  if (region->sgne2) d_free_2d(region->sgne2);
  if (rn) {
    region->be2 = i_alloc_2d(rn, region->nfaces);
    region->sgne2 = d_alloc_2d(rn, region->nfaces);

    /* Populate the face arrays                                      */
    memset(region->nbe2, 0, region->nfaces * sizeof(int));
    for (cc = 1; cc <= region->nvec; cc++) {
      c = region->vec[cc];
      rn = (int)mask[win->yp1[c]];
      if (rn != NOTVALID && n != rn) {
	region->sgne2[rn][region->nbe2[rn]] = -1.0;
	region->be2[rn][region->nbe2[rn]++] = win->yp1[c];
      }
      rn = (int)mask[win->ym1[c]];
      if (rn != NOTVALID && n != rn) {
	region->sgne2[rn][region->nbe2[rn]] = 1.0;
	region->be2[rn][region->nbe2[rn]++] = c;
      }
    }
    for (m = 0; m < win->nobc; m++) {
      open_bdrys_t *open = win->open[m];
      region->nbe2[m + nregions] = 0;
      if (open->ocodey & (B_EDGE|F_EDGE)) {
	for (cc = 1; cc <= open->no3_t; cc++) {
	  c = open->obc_t[cc];
	  region->sgne2[m + nregions][region->nbe2[m + nregions]] = 1.0;
	  if (open->ocodey & F_EDGE) 
	    region->sgne2[m + nregions][region->nbe2[m + nregions]] = -1.0;
	  if (region->mode & RG_BDRY)
	    for (bb = 0; bb < region->obcz; bb++)
	      c = open->nmap[c];
	  cr = (open->ocodey & F_EDGE && open->stagger & OUTFACE) ? win->yp1[c] : c;
	  if (n == (int)mask[c]) {
	    region->be2[m + nregions][region->nbe2[m + nregions]] = cr;
	    region->nbe2[m + nregions]++;
	  }
	}
      }
      /* Tangential for RG_BDRY                                      */
      if (dotan && open->ocodex & (L_EDGE|R_EDGE)) {
	for (cc = 1; cc <= open->no3_t; cc++) {
	  c = open->obc_t[cc];
	  if (region->mode & RG_BDRY) {
	    for (bb = 0; bb < region->obcz; bb++) {
	      cr = win->yp1[c];
	      if (n == (int)mask[cr]) {
		region->be2[m + nregions][region->nbe2[m + nregions]] = cr;
		region->sgne2[m + nregions][region->nbe2[m + nregions]] = 1.0;
		region->nbe2[m + nregions]++;
	      }
	      cr = win->ym1[c];
	      if (n == (int)mask[cr]) {
		region->be2[m + nregions][region->nbe2[m + nregions]] = c;
		region->sgne2[m + nregions][region->nbe2[m + nregions]] = -1.0;
		region->nbe2[m + nregions]++;
	      }
	      c = open->nmap[c];
	    }
	  }
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the region boundaries in the z direction                    */
  if (region->nbz) i_free_1d(region->nbz);
  region->nbz = i_alloc_1d(region->nfaces);
  memset(region->nbz, 0, region->nfaces * sizeof(int));
  if (region->sgnz) d_free_1d(region->sgnz);
  region->sgnz = d_alloc_1d(region->nfaces);
  /* Region to region boundaries                                     */
  for (cc = 1; cc <= region->nvec; cc++) {
    c = region->vec[cc];
    rn = (int)mask[win->zp1[c]];
    if (rn != NOTVALID && n != rn) {
      region->nbz[rn]++;
      region->sgnz[rn] = -1.0;
    }
    rn = (int)mask[win->zm1[c]];
    if (rn != NOTVALID && n != rn) {
      region->nbz[rn]++;
      region->sgnz[rn] = 1.0;
    }
  }
  /* Get the maximum face size                                       */
  rn = 0;
  for (m = 0; m < region->nfaces; m++)
    if (region->nbz[m] > rn) rn = region->nbz[m];
  if (region->bz) i_free_2d(region->bz);
  if (rn) region->bz = i_alloc_2d(rn, region->nfaces);
  /* Populate the face arrays                                        */
  memset(region->nbz, 0, region->nfaces * sizeof(int));
  for (cc = 1; cc <= region->nvec; cc++) {
    c = region->vec[cc];
    rn = (int)mask[win->zp1[c]];
    if (rn != NOTVALID && n != rn) 
      region->bz[rn][region->nbz[rn]++] = win->zp1[c];
    rn = (int)mask[win->zm1[c]];
    if (rn != NOTVALID && n != rn) 
      region->bz[rn][region->nbz[rn]++] = c;
  }

  /*-----------------------------------------------------------------*/
  /* Get a map for all boundaries of the region                      */
  region->nboundaries = 0;
  if (region->bmap) i_free_1d(region->bmap);
  region->bmap = i_alloc_1d(region->nfaces);
  if (region->rmap) i_free_1d(region->rmap);
  region->rmap = i_alloc_1d(region->nfaces);
  for (m = 0; m < region->nfaces; m++) {
    region->rmap[m] = NOTVALID;
    if (region->mode & RG_ALLT) {
      region->bmap[region->nboundaries] = m;
      region->rmap[m] = region->nboundaries;
      region->nboundaries++;
    } else {
      if (m != region->id  && (region->nbe1[m] || region->nbe2[m] || region->nbz[m])) {
	region->bmap[region->nboundaries] = m;
	region->rmap[m] = region->nboundaries;
	region->nboundaries++;
      }
    }
  }

  if (region->dz) d_free_1d(region->dz);
  if (region->nvec) region->dz = d_alloc_1d(region->nvec + 1);

  if (region->fluxname) free((char **)region->fluxname);
  if (region->nboundaries) 
    region->fluxname = (char **)malloc(region->nboundaries * sizeof(char *));
  for (m = 0; m < region->nboundaries; m++) {
    region->fluxname[m] = (char *)malloc(sizeof(char)*MAXSTRLEN);
    memset(region->fluxname[m], 0, sizeof(char)*MAXSTRLEN);
    rn = region->bmap[m];
    if (rn < nregions) {
      if (region->nbe1[rn] || region->nbe2[rn])
	sprintf(region->fluxname[m], "horizontal flux from region %d to %d", region->id, rn);
      else
	sprintf(region->fluxname[m], "vertical flux from region %d to %d", region->id, rn);
    } else {
      sprintf(region->fluxname[m], "flux from region %d through OBC %s", region->id, 
	      win->open[rn - nregions]->name);
    }
  }
}

/* END get_regions()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to transfer region information to the master              */
/*-------------------------------------------------------------------*/
void region_transfer(master_t *master, geometry_t *window)
{
  geometry_t *geom = master->geom;
  int n, m, wn, c, cc, tt, bb, bg;
  region_t *region, *region_w;

  if (geom->nwindows > 1) {
    for (n = 0; n < geom->nregions; n++) {
      region = geom->region[n];
      for (m = 0; m < window->nregions; m++) {
	region_w = window->region[m];
	if (region->id == region_w->id) {
	  region->transferf = region_w->transferf;
	  region->pt = region_w->pt;
	  region->step = region_w->step;
	  region->vfluxp += region_w->vfluxp;
	  region->vfluxn += region_w->vfluxn;
	  for (tt = 0; tt < region->nvar; tt++) {
	    region->emass[tt] += region_w->emass[tt];
	    region->mmass[tt] += region_w->mmass[tt];
	    region->stdev[tt] += region_w->stdev[tt];
	    if (region_w->transferf & RG_START)
	      region->smass[tt] += region_w->smass[tt];
	    if (region->mode & RG_PSS)
	      region->pssf[tt] += region_w->pssf[tt];
	    if (region->mode & RG_VERR)
	      region->merr[tt] += region_w->merr[tt];
	    if (region->mode & RG_GLOB)
	      region->gfer[tt] += region_w->gfer[tt];
	    if (region->mode & RG_SED)
	      region->sedtr[tt] += region_w->sedtr[tt];
	    if (region->mode & RG_ECO)
	      region->eco[tt] += region_w->eco[tt];
	    for (bb = 0; bb < region_w->nboundaries; bb++) {
	      bg = region_w->gmap[bb];
	      region->trflux[tt][bg] += region_w->trflux[tt][bb];
	      region->w1[tt][bg] += region_w->w1[tt][bb];
	    }
	  }
	}
      }
    }
  }
}

/* END region_transfer()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to compute the mass in each region                        */
/*-------------------------------------------------------------------*/
void region_mass(geometry_t *window, 
		 window_t *windat, 
		 win_priv_t *wincon)
{
  int n, tn, tt, c, cc, cs, bb, bn;
  double dtr, mass;
  double dt = windat->dttr;
  int pf = -1;
  double scale = 1.0 / 86400.0;

  for (n = 0; n < window->nregions; n++) {
    region_t *region = window->region[n];

    /*---------------------------------------------------------------*/
    /* Get the cell thickness for wet cells                          */
    for (cc = 1; cc <= region->nvec; cc++) {
      c = region->vec[cc];
      region->dz[cc] = (wincon->c1[c]) ? wincon->dz[c]: 0.0;
    }
  }

  for (n = 0; n < window->nregions; n++) {
    region_t *region = window->region[n];

    /*---------------------------------------------------------------*/
    /* Get the volumes and volume fluxes store in index 0. Fluxes    */
    /* are computed here since the region_flux computes fluxes for   */
    /* a specific tracer in the tracer list (i.e. volume is not      */
    /* included in the tracer list).                                 */
    region->emass[0] = 0.0;
    for (cc = 1; cc <= region->nvec; cc++) {
      c = region->vec[cc];
      cs = window->m2d[c];
      region->emass[0] += window->cellarea[cs] * region->dz[cc];
    }

    /*---------------------------------------------------------------*/
    /* Get the remaining tracers total mass                          */
    for (tt = 1; tt < region->nvar; tt++) {
      tn = region->var[tt];
      region->emass[tt] = 0.0;
      for (cc = 1; cc <= region->nvec; cc++) {
	c = region->vec[cc];
	cs = window->m2d[c];
	mass = window->cellarea[cs] * region->dz[cc] * windat->tr_wc[tn][c];
	region->emass[tt] += mass;
	if (region->mode & RG_GLOB)
	  region->gfer[tt] += mass;
	if (region->mode & RG_SED)
	  region->sedtr[tt] -= mass;
      }
    }

    /*---------------------------------------------------------------*/
    /* Get the mean                                                  */
    for (tt = 0; tt < region->nvar; tt++) {
      region->mmass[tt] = (region->mmass[tt] * region->pdt +
			   region->emass[tt] * dt) / (region->pdt + dt);
    }

    /*---------------------------------------------------------------*/
    /* Get the residence time                                        */
    /*
    dtr = scale * region->emass[0] / (min(fabs(region->vfluxp), fabs(region->vfluxn)));
    region->residence = (region->residence * region->pdt + dtr * dt) / (region->pdt + dt);
    */
    region->residence = (region->pdt + dt) * scale * region->mmass[0] / 
      (min(fabs(region->vfluxp), fabs(region->vfluxn)));
    region->residencen = (region->pdt + dt) * scale * region->mmass[0] / 
      fabs(region->vfluxp + region->vfluxn);
    for (cc = 1; cc <= region->nvec; cc++) {
      c = region->vec[cc];
      windat->regres[c] = region->residence;
    }

    /*---------------------------------------------------------------*/
    /* Update the time                                               */
    region->pdt += dt;

    /*---------------------------------------------------------------*/
    /* Get the standard deviation                                    */
    region->step += 1.0;
    for (tt = 0; tt < region->nvar; tt++) {
      double var, sump, sum, sum2, nstep;
      nstep = region->step;
      /* Get the sum of tracer at the current nstep                  */
      sump = region->mmass[tt] * nstep;
      /* Get the sum of tracer trn+1 up to the previous nstep        */
      sum = sump - region->emass[tt];
      /* Get the sum of tracer squared up to the previous nstep      */
      sum = sum * sum / (nstep - 1.0);
      sum2 = region->stdev[tt] * (nstep - 2.0) + sum;
      /* Increment the sum of tracer squared (i.e. get the sum of    */
      /* tracer squared up to the current nstep).                    */
      sum2 += region->emass[tt] * region->emass[tt];
      /* Get the variance                                            */
      sum =  (sump * sump) / nstep;
      var = (sum2 - sum) / (nstep - 1.0);
      /* Get the standard deviation                                  */
      region->stdev[tt] = (nstep > 1.0 && var > 0.0) ? sqrt(var) : 0.0;
    }

    /* Initialise                                                    */
    if (region->transferf & RG_BEGIN) {
      memcpy(region->smass, region->emass, region->nvar * sizeof(double));
    }

    /* Diagnostics                                                   */
    if(region->id == pf) {
      tt = 0;
      printf("%f %f : %f %f %f : %f ",master->days,region->next/86400,region->smass[tt],region->emass[tt],get_totflux(region, tt),region->merr[tt]/region->smass[tt]);
      mass = region->smass[tt]-region->emass[tt]+get_totflux(region,tt)+region->merr[tt];
      if (region->mode & RG_PSS) {
	printf("%f ",region->pssf[tt]);
	mass += region->pssf[tt];
      }
      if (region->mode & RG_GLOB) {
	printf("%f ",region->gfer[tt]);
	mass += region->gfer[tt];
      }
      printf(": %f\n",mass);
    }
  }
}

/* END region_mass()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets the timing flags for the regions                             */
/*-------------------------------------------------------------------*/
void region_schedule(region_t *region, double dt, double trem, char *timeunit)
{
  int n, tt, bb;
  double dtr;
  int trfin = region->transferf;

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */
  if (region->transferf & RG_START) {
    memcpy(region->smass, region->emass, region->nvar * sizeof(double));
    memcpy(region->mmass, region->emass, region->nvar * sizeof(double));
    /*region->residence = 0.0;*/
    region->vfluxp = region->vfluxn = 0.0;
    region->step = 0.0;
    region->pdt = 0.0;
    if (region->mode & RG_VERR)
      memset(region->merr, 0, region->nvar * sizeof(double));
    if (region->mode & RG_PSS)
      memset(region->pssf, 0, region->nvar * sizeof(double));
    if (region->mode & RG_GLOB)
      memset(region->gfer, 0, region->nvar * sizeof(double));
    if (region->mode & RG_SED)
      memset(region->sedtr, 0, region->nvar * sizeof(double));
    if (region->mode & RG_ECO)
      memset(region->eco, 0, region->nvar * sizeof(double));
    for (tt = 0; tt < region->nvar; tt++) {
      for (bb = 0; bb < region->nboundaries; bb++)
	region->trflux[tt][bb] = 0.0;
    }
  }
  /*region->transferf = RG_NONE; pre v5373 code */
  region->transferf = (trfin & RG_BEGIN && trem > dt) ? RG_NONE|RG_BEGIN : RG_NONE;

  /*---------------------------------------------------------------*/
  /* Set the flag for dumping when region_transfer() is called.    */
  /* The masses and fluxes are computed in blocks of region->dt;   */
  /* end mass is at region->next-dt, and the next computational    */
  /* block starts at region->next. Fluxes are computed from the    */
  /* start of the block to region->next-dt.                        */
  /* region->t + dt added post v5373                               */
  if (trem - dt <= 0 && region->t + dt >= region->next - DT_EPS) {
    region->transferf |= RG_WRITE;
    /*region->pt = region->t;  pre v5373 code */
    region->pt = region->t + dt;
  }

  /*---------------------------------------------------------------*/
  /* Set the flag for the start of a new block. The mean array and */
  /* fluxes are re-initialised to zero here, and in the transfer   */
  /* routine the start mass is set.                                */
  /* region->t + dt added post v5373                               */
  if (trem - dt <= 0 && region->t + dt >= region->next - DT_EPS) {
    if (region->dt == YEARLY)
      dtr = next_year(region->t, timeunit);
    if (region->dt == SEASONAL)
      dtr = next_season(region->t, timeunit, &n); 
    else if (region->dt == MONTHLY)
      dtr = next_month(region->t, timeunit, &n); 
    else
      dtr = region->dt;
    /*region->next = region->t + dtr;  pre v5373 code */
    region->next = region->next + dtr;
    
    region->transferf |= RG_START;
    /*region->pt = region->t; pre v5373 code (not required - set above) */
  }
  region->t += dt;

  /* Set the flag for the first interval start                       */
  /*if (trfin & RG_BEGIN) region->transferf |= RG_START; pre v5373 code */
  if (trfin & RG_BEGIN && trem - dt <= 0) region->transferf |= RG_START;
}

/* END region_schedule()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to compute the mass fluxes across each region boundary    */
/* for the hydro model.                                              */
/*-------------------------------------------------------------------*/
void region_flux_coup(geometry_t *window, window_t *windat, win_priv_t *wincon,
		      double *Fx, double *Fy, double *Fz, double dt, int trn)
{
  int n, tt, tn, c, cc, cs, bb, bn;

  for (n = 0; n < window->nregions; n++) {
    region_t *region = window->region[n];

    /*---------------------------------------------------------------*/
    /* Is this tracer a variable within the region?                  */
    if (!(tt = find_trindex(region, trn))) continue;

    /*---------------------------------------------------------------*/
    /* Compute the fluxes                                            */
    for (bb = 0; bb < region->nboundaries; bb++) {
      bn = region->bmap[bb];
      if (region->nbe1[bn]) {
	for (cc = 0; cc < region->nbe1[bn]; cc++) {
	  c = region->be1[bn][cc];
	  region->trflux[tt][bb] += (Fx[c] * dt * region->sgne1[bn][cc]);
	}
      }
      if (region->nbe2[bn]) {
	for (cc = 0; cc < region->nbe2[bn]; cc++) {
	  c = region->be2[bn][cc];
	  region->trflux[tt][bb] += (Fy[c] * dt * region->sgne2[bn][cc]);
	}
      }
      if (region->nbz[bn]) {
	for (cc = 0; cc < region->nbz[bn]; cc++) {
	  c = region->bz[bn][cc];
	  region->trflux[tt][bb] += (Fz[c] * region->sgnz[bn]);
	}
      }
    }
  }
}

/* END region_flux_coup()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Volume fluxes for the hydro model. The first computation in the   */
/* block of region->dt the fluxes are not computed.                  */
/*-------------------------------------------------------------------*/
void region_volume_flux_coup(geometry_t *window, 
			     window_t *windat, 
			     win_priv_t *wincon,
			     double dt, double trem)
{
  int n, cc, c, cs, bb, bn;

  for (n = 0; n < window->nregions; n++) {
    region_t *region = window->region[n];
    double flux;

    /*---------------------------------------------------------------*/
    /* Set schedule flags                                            */
    region_schedule(region, dt, trem, wincon->timeunit);

    /*---------------------------------------------------------------*/
    /* Get the volume fluxes                                         */
    for (bb = 0; bb < region->nboundaries; bb++) {
      bn = region->bmap[bb];
      if (region->nbe1[bn]) {
	for (cc = 0; cc < region->nbe1[bn]; cc++) {
	  c = region->be1[bn][cc];
	  flux = (windat->u1flux3d[c] * dt * region->sgne1[bn][cc]);
	  /*if(region->id==2&&c==window->m2d[c])printf("u1 %d %d(%d %d) %f\n",bb,c,window->s2i[c],window->s2j[c],flux);*/
	  region->trflux[0][bb] += flux;
	  if (flux >= 0.0)
	    region->vfluxp += flux;
	  else
	    region->vfluxn += flux;
	}
      }
      if (region->nbe2[bn]) {
	for (cc = 0; cc < region->nbe2[bn]; cc++) {
	  c = region->be2[bn][cc];
	  flux = (windat->u2flux3d[c] * dt * region->sgne2[bn][cc]);
	  /*if(region->id==2&&c==window->m2d[c])printf("u2 %d %d(%d %d) %f %f\n",bb,c,window->s2i[c],window->s2j[c],flux,windat->u2[c]);*/
	  region->trflux[0][bb] += flux;
	  if (flux >= 0.0)
	    region->vfluxp += flux;
	  else
	    region->vfluxn += flux;
	}
      }
      if (region->nbz[bn]) {
	for (cc = 0; cc < region->nbz[bn]; cc++) {
	  c = region->bz[bn][cc];
	  cs = window->m2d[c];
	  flux = (window->cellarea[cs] * windat->w[c] * dt * region->sgnz[bn]);
	  region->trflux[0][bb] += flux;
	  if (flux >= 0.0)
	    region->vfluxp += flux;
	  else
	    region->vfluxn += flux;
	}
      }
    }
  }    
}

/* END region_volume_flux_coup()                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads the region mask from file and sets ghost cells              */
/*-------------------------------------------------------------------*/
int read_region(master_t *master, parameters_t *params)
{
  geometry_t *geom = master->geom;
  timeseries_t *ts;
  int timeindex = 0;
  int fid, ncerr, id, nregions;
  int n, c, cc, cs, c1;
  double d1;
  size_t nr, nce1, nce2, nz, record;
  char buf[MAXSTRLEN];
  int rg[MAXREGIONS], found;

  /* Open the file                                                   */
  if ((ncerr = nc_open(params->regions, NC_NOWRITE, &fid)) != NC_NOERR) {
    hd_warn("read_region: Can't find regions file %s\n", params->regions);
    return(0);
  }

  /* Get the dimensions                                              */
  nc_inq_dimlen(fid, ncw_dim_id(fid, "nbox"), &nr);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "i_centre"), &nce1);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "j_centre"), &nce2);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "k_centre"), &nz);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "record"), &record);

  /* Exact dimesions; read data 'as is'.                             */
  if ((int)nce1 == geom->nce1 && (int)nce2 == geom->nce2 && (int)nz == geom->nz) {
    if (record > 1)
      timeindex = dump_choose_by_time(params, fid, master->t);
    if (timeindex > (int)(record - 1)) {
      hd_warn("read_region: dump %d does not exist in file %s!\n", timeindex, params->regions);
      return(0);
    }

    /* Without this initialisation we'll leave a valid 0 in places where it shouldn't be */
    for (c = 1; c <= geom->sgnum; c++) master->regionid[c] = NOTVALID;

    /* Get the region id from file */
    read3d(geom, fid, "boxnos", master->regionid, timeindex, geom->nz, geom->nce2, geom->nce1);
    nc_close(fid);
  } else {  /* Interpolate onto the grid                             */
    hd_warn("read_region: interpolating regions from file %s onto the grid.\n", 
	    params->regions);
    ts = hd_ts_read(master, params->regions, 0);

    id = ts_get_index(ts, fv_get_varname(buf, "boxnos", buf));
    if (id < 0) {
      hd_warn("read_region: The file '%s' does not contain the variable 'boxnos'.\n",
	      params->regions);
      return(0);
    }
    for (c = 1; c <= geom->sgnum; c++) master->regionid[c] = NOTVALID;
    for (cc = 1; cc <= geom->b3_t; cc++) {
      c = geom->w3_t[cc];
      cs = geom->m2d[c];
      d1 = ts_eval_xyz(ts, id, master->t, geom->cellx[cs], geom->celly[cs], geom->cellz[c]);
      if (fmod(d1, 1) <= 0.5)		     
	master->regionid[c] = floor(d1);
      else
	master->regionid[c] = ceil(d1);
    }
    nc_close(fid);
  }

  /* Count the number of regions                                     */
  nregions = 0;
  for (cc = 1; cc <= geom->b3_t; cc++) {
    c = geom->w3_t[cc];
    if (isnan(master->regionid[c])) master->regionid[c] = NOTVALID;
    found = 1;
    for (n = 0; n < nregions; n++)
      if (rg[n] == (int)master->regionid[c]) found = 0;
    if (master->regionid[c] != NOTVALID && found) {
      rg[nregions] = (int)master->regionid[c];
      nregions++;
    }
  }

  /* Set ghost points                                                */
  for (cc = 1; cc <= geom->nbpt; cc++) {
    c = geom->bpt[cc];
    c1 = geom->bin[cc];
    master->regionid[c] = master->regionid[c1];
  }
  /* Set sediment points                                             */
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->bot_t[cc];
    master->regionid[geom->zm1[c]] = master->regionid[c];
  }
  /* Set RIGHT and FRONT OBCs                                        */
  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    if (open->ocodex & R_EDGE && open->stagger & OUTFACE ) {
      for (cc = 1; cc <= open->no3_e1; cc++) {
        c = open->obc_e1[cc];
        c1 = open->oi1_e1[cc];
        master->regionid[c] = master->regionid[c1];
      }
    }
    if (open->ocodey & F_EDGE && open->stagger & OUTFACE) {
      for (cc = 1; cc <= open->no3_e2; cc++) {
        c = open->obc_e2[cc];
        c1 = open->oi1_e2[cc];
        master->regionid[c] = master->regionid[c1];
      }
    }
  }
  return(nregions);
}

/* END read_region()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to allocate memory for the regions structure              */
/*-------------------------------------------------------------------*/
region_t *region_alloc(void)
{
  region_t *region = (region_t *)malloc(sizeof(region_t));
  memset(region, 0, sizeof(region_t));
  return region;
}

/* END region_alloc()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to allocate memory for the regions structure              */
/*-------------------------------------------------------------------*/
void region_clean(geometry_t *win)
{
  int n;

  for (n = 0; n < win->nregions; n++) {
    region_t *region = win->region[n];
    if (region->nbe1) i_free_1d(region->nbe1);
    if (region->nbe2) i_free_1d(region->nbe2);
    if (region->be1) i_free_2d(region->be1);
    if (region->be2) i_free_2d(region->be2);
    if (region->bmap) i_free_1d(region->bmap);
    if (region->rmap) i_free_1d(region->rmap);
    if (region->gmap) i_free_1d(region->gmap);
    if (region->sgne1) d_free_2d(region->sgne1);
    if (region->sgne2) d_free_2d(region->sgne2);
    if (region->trflux) d_free_2d(region->trflux);
    if (region->w1) d_free_2d(region->w1);
    if (region->dz) d_free_1d(region->dz);
    if (region->fluxname) free((char **)region->fluxname);
    d_free_1d(region->mmass);
    d_free_1d(region->stdev);
    d_free_1d(region->smass);
    d_free_1d(region->emass);
    if (region->pssf) d_free_1d(region->pssf);
    if (region->merr) d_free_1d(region->merr);
    if (region->gfer) d_free_1d(region->gfer);
    if (region->sedtr) d_free_1d(region->sedtr);
    if (region->eco) d_free_1d(region->eco);
    if (region->fp) fclose(region->fp);
    free(region);
  }
  free(win->region);
}

/* END region_clean()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads the region mask from file and sets ghost cells              */
/*-------------------------------------------------------------------*/
int read_regioni(master_t *master, char *rname, double *regionid)
{
  geometry_t *geom = master->geom;
  timeseries_t *ts;
  int timeindex = 0;
  int fid, ncerr, id, nregions;
  int n, c, cc, cs, c1;
  double d1;
  size_t nr, nce1, nce2, nz, record;
  char buf[MAXSTRLEN];
  int rg[MAXREGIONS], found;

  /* Open the file                                                   */
  if ((ncerr = nc_open(rname, NC_NOWRITE, &fid)) != NC_NOERR) {
    hd_warn("read_region: Can't find regions file %s\n", rname);
    return(0);
  }

  /* Get the dimensions                                              */
  nc_inq_dimlen(fid, ncw_dim_id(fid, "nbox"), &nr);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "i_centre"), &nce1);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "j_centre"), &nce2);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "k_centre"), &nz);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "record"), &record);

  /* Exact dimesions; read data 'as is'.                             */
  if ((int)nce1 == geom->nce1 && (int)nce2 == geom->nce2 && (int)nz == geom->nz) {
    if (record > 1)
      timeindex = dump_choose_by_time_s(fid, master->t);
    if (timeindex > (int)(record - 1)) {
      hd_warn("get_region: dump %d does not exist in file %s!\n", timeindex, rname);
      return(0);
    }
    /* Get the region id from file */
    read3d(geom, fid, "boxnos", regionid, timeindex, geom->nz, geom->nce2, geom->nce1);
    nc_close(fid);
  } else {  /* Interpolate onto the grid                             */
    hd_warn("read_region: interpolating regions from file %s onto the grid.\n", 
	    rname);
    ts = hd_ts_read(master, rname, 0);

    id = ts_get_index(ts, fv_get_varname(buf, "boxnos", buf));
    if (id < 0) {
      hd_warn("read_region: The file '%s' does not contain the variable 'boxnos'.\n",
	      rname);
      return(0);
    }
    for (c = 1; c <= geom->sgnum; c++) regionid[c] = NOTVALID;
    for (cc = 1; cc <= geom->b3_t; cc++) {
      c = geom->w3_t[cc];
      cs = geom->m2d[c];
      d1 = ts_eval_xyz(ts, id, master->t, geom->cellx[cs], geom->celly[cs], geom->cellz[c]);
      if (fmod(d1, 1) <= 0.5)		     
	regionid[c] = floor(d1);
      else
	regionid[c] = ceil(d1);
    }
    nc_close(fid);
  }

  /* Count the number of regions                                     */
  nregions = 0;
  for (cc = 1; cc <= geom->b3_t; cc++) {
    c = geom->w3_t[cc];
    if (isnan(regionid[c])) regionid[c] = NOTVALID;
    found = 1;
    for (n = 0; n < nregions; n++)
      if (rg[n] == (int)regionid[c]) found = 0;
    if (regionid[c] != NOTVALID && found) {
      rg[nregions] = (int)regionid[c];
      nregions++;
    }
  }

  /* Set ghost points                                                */
  for (cc = 1; cc <= geom->nbpt; cc++) {
    c = geom->bpt[cc];
    c1 = geom->bin[cc];
    regionid[c] = regionid[c1];
  }
  /* Set sediment points                                             */
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->bot_t[cc];
    regionid[geom->zm1[c]] = regionid[c];
  }
  /* Set RIGHT and FRONT OBCs                                        */
  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    if (open->ocodex & R_EDGE && open->stagger & OUTFACE ) {
      for (cc = 1; cc <= open->no3_e1; cc++) {
        c = open->obc_e1[cc];
        c1 = open->oi1_e1[cc];
        regionid[c] = regionid[c1];
      }
    }
    if (open->ocodey & F_EDGE && open->stagger & OUTFACE) {
      for (cc = 1; cc <= open->no3_e2; cc++) {
        c = open->obc_e2[cc];
        c1 = open->oi1_e2[cc];
        regionid[c] = regionid[c1];
      }
    }
  }
  return(nregions);
}

/* END read_regioni()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Controls the dumping of regions                                   */
/*-------------------------------------------------------------------*/
void dump_regions(master_t *master)
{
  geometry_t *geom = master->geom;
  region_t *region;
  int n, tt, prf = 2;

  if (prf) {
    for (n = 0; n < geom->nregions; n++) {
      region = geom->region[n];
      if (region->transferf & RG_WRITE) {
	if (prf == 1)
	  region_print(master, region);
	else if (prf == 2)
	  region_write(master, region);
      }
    }
  }

  /* Re-initialise global regions for the next timestep */
  if (geom->nwindows > 1) {
    for (n = 0; n < geom->nregions; n++) {
      region = geom->region[n];
      if (region->transferf & RG_WRITE)
	memset(region->smass, 0, region->nvar * sizeof(double));
      memset(region->emass, 0, region->nvar * sizeof(double));
      memset(region->mmass, 0, region->nvar * sizeof(double));
      memset(region->stdev, 0, region->nvar * sizeof(double));
      region->vfluxp = region->vfluxn = 0.0;
      region->step = 0.0;
      if (region->mode & RG_PSS)
	memset(region->pssf, 0, region->nvar * sizeof(double));
      if (region->mode & RG_VERR) {
	memset(region->mmass, 0, region->nvar * sizeof(double));
	memset(region->stdev, 0, region->nvar * sizeof(double));
      }
      if (region->mode & RG_GLOB)
	memset(region->gfer, 0, region->nvar * sizeof(double));
      if (region->mode & RG_SED)
	memset(region->sedtr, 0, region->nvar * sizeof(double));
      if (region->mode & RG_ECO)
	memset(region->eco, 0, region->nvar * sizeof(double));
      for (tt = 0; tt < region->nvar; tt++)
	memset(region->trflux[tt], 0, region->nboundaries * sizeof(double));
    }
  }
}

/* END dump_regions()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Prints the mass and flux infromation in a region                  */
/*-------------------------------------------------------------------*/
void region_print(master_t *master, region_t *region)
{
  geometry_t *geom = master->geom;
  int tt, tn, bb;
  char key[MAXSTRLEN];
  double time;
  printf("\n-----------------------\n");
  time = region->pt;
  tm_change_time_units(master->timeunit, master->output_tunit, &time, 1);
  printf("Time %f : budget for %s\n", time, region->name);
  for (tt = 0; tt < region->nvar; tt++) {
    tn = region->var[tt];
    if (tt) {
      printf("Tracer %d (%s)\n", tn, master->trinfo_3d[tn].name);
      strcpy(key, "Mass");
    } else {
      printf("Volume\n");
      strcpy(key, "Volume");
    }
    printf("  %s at start of interval = %6.3e\n",key, region->smass[tt]);
    printf("  %s at end of interval = %6.3e\n",key, region->emass[tt]);
    printf("  %s average over interval = %6.3e\n",key, region->mmass[tt]);
    printf("  %s standard deviation over interval = %6.3e\n",key, region->stdev[tt]);
    for (bb = 0; bb < region->nboundaries; bb++) {
      printf("  %s = %6.3e\n", region->fluxname[bb], region->trflux[tt][bb]);
    }
  }
}

/* End region_print()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Prints the mass and flux infromation in a region                  */
/*-------------------------------------------------------------------*/
void region_write(master_t *master, region_t *region)
{
  geometry_t *geom = master->geom;
  int n, tt, tn, bb, bn;
  double d1, time;

  /* Create the file if required */
  if (region->fp == NULL) {
    timeseries_t *ots = NULL;
    char buf[MAXSTRLEN];

    /* Construct regions file name */
    if (strlen(master->opath))
      sprintf(buf, "%sregion%d.ts", master->opath, region->id);
    else
      sprintf(buf, "region%d.ts", region->id);

    /* For restarts, read in the old regions file */
    if (forced_restart) {
      FILE *fp = fopen(buf, "r");
      if (fp != NULL) {
	fclose(fp);
	ots = (timeseries_t *)malloc(sizeof(timeseries_t));
	memset(ots, 0, sizeof(timeseries_t));
	ts_read(buf, ots);
      }
    }
    region->fp = fopen(buf, "w");
    n = 1 + 2; // +2 added for residence
    for (tt = 0; tt < region->nvar; tt++) {
      n += 5; // tt=0 corresponds to volume
      if (region->mode & RG_VERR) n++;
      if (tt && region->mode & RG_GLOB) n++;
      if (region->mode & RG_PSS) n++;
      if (tt && region->mode & RG_SED) n++;
      if (tt && region->mode & RG_ECO) n++;
      n += region->nboundaries;
    }
    fprintf(region->fp, "## COLUMNS %d\n", n);
    fprintf(region->fp, "##\n");
    fprintf(region->fp, "## COLUMN1.name  Time\n");
    fprintf(region->fp, "## COLUMN1.long_name  Time\n");
    fprintf(region->fp,
	    "## COLUMN1.units  %s\n", master->output_tunit);
    fprintf(region->fp, "## COLUMN1.missing_value -999\n");
    fprintf(region->fp, "##\n");
    n = 2;
    fprintf(region->fp, "## COLUMN%d.name  residence_time\n", n);
    fprintf(region->fp, "## COLUMN%d.long_name  Residence time\n", n);
    fprintf(region->fp, "## COLUMN%d.units  days\n", n);
    fprintf(region->fp, "## COLUMN%d.missing_value -999\n",n);
    fprintf(region->fp, "##\n");
    n++;
    fprintf(region->fp, "## COLUMN%d.name  residence_time_net\n", n);
    fprintf(region->fp, "## COLUMN%d.long_name  Residence time using net flow\n", n);
    fprintf(region->fp, "## COLUMN%d.units  days\n", n);
    fprintf(region->fp, "## COLUMN%d.missing_value -999\n",n);
    fprintf(region->fp, "##\n");
    n++;
    fprintf(region->fp, "## COLUMN%d.name  start_volume\n", n);
    fprintf(region->fp, "## COLUMN%d.long_name  Start volume\n", n);
    fprintf(region->fp, "## COLUMN%d.units  m^3\n", n);
    fprintf(region->fp, "## COLUMN%d.missing_value -999\n",n);
    fprintf(region->fp, "##\n");
    n++;
    fprintf(region->fp, "## COLUMN%d.name  end_volume\n", n);
    fprintf(region->fp, "## COLUMN%d.long_name  End volume\n", n);
    fprintf(region->fp, "## COLUMN%d.units  m^3\n", n);
    fprintf(region->fp, "## COLUMN%d.missing_value -999\n",n);
    fprintf(region->fp, "##\n");
    n++;
    fprintf(region->fp, "## COLUMN%d.name  mean_volume\n", n);
    fprintf(region->fp, "## COLUMN%d.long_name  Mean volume\n", n);
    fprintf(region->fp, "## COLUMN%d.units  m^3\n", n);
    fprintf(region->fp, "## COLUMN%d.missing_value -999\n",n);
    fprintf(region->fp, "##\n");
    n++;
    fprintf(region->fp, "## COLUMN%d.name  stdev_volume\n", n);
    fprintf(region->fp, "## COLUMN%d.long_name  Standard deviation volume\n", n);
    fprintf(region->fp, "## COLUMN%d.units  m^3\n", n);
    fprintf(region->fp, "## COLUMN%d.missing_value -999\n",n);
    fprintf(region->fp, "##\n");
    n++;
    fprintf(region->fp, "## COLUMN%d.name  volume_budget\n", n);
    fprintf(region->fp, "## COLUMN%d.long_name  Volume budget\n", n);
    fprintf(region->fp, "## COLUMN%d.units  m^3\n", n);
    fprintf(region->fp, "## COLUMN%d.missing_value -999\n",n);
    fprintf(region->fp, "##\n");
    n++;
    if (region->mode & RG_VERR) {
      fprintf(region->fp, "## COLUMN%d.name  volume_error\n", n);
      fprintf(region->fp, "## COLUMN%d.long_name  Volume error\n", n);
      fprintf(region->fp, "## COLUMN%d.units  m^3\n", n);
      fprintf(region->fp, "## COLUMN%d.missing_value -999\n",n);
      fprintf(region->fp, "##\n");
      n++;
    }
    if (region->mode & RG_PSS) {
      fprintf(region->fp, "## COLUMN%d.name  pss_flux\n", n);
      fprintf(region->fp, "## COLUMN%d.long_name  Point source input\n", n);
      fprintf(region->fp, "## COLUMN%d.units  m^3\n", n);
      fprintf(region->fp, "## COLUMN%d.missing_value -999\n",n);
      fprintf(region->fp, "##\n");
      n++;
    }
    for (bb = 0; bb < region->nboundaries; bb++) {
      bn = region->bmap[bb];
      fprintf(region->fp, "## COLUMN%d.name  volume_flux%d-%d\n", n, region->id, bn);
      fprintf(region->fp, "## COLUMN%d.long_name  Volume %s\n", n, region->fluxname[bb]);
      fprintf(region->fp, "## COLUMN%d.units  m^3\n", n);
      fprintf(region->fp, "## COLUMN%d.missing_value -999\n",n);
      fprintf(region->fp, "##\n");
      n++;
    }
    for (tt = 1; tt < region->nvar; tt++) {
      tn = region->var[tt];
      fprintf(region->fp, "## COLUMN%d.name  start_mass_%s\n", n, master->trinfo_3d[tn].name);
      fprintf(region->fp, "## COLUMN%d.long_name  Tracer %s start mass\n", n,
	      master->trinfo_3d[tn].name);
      fprintf(region->fp, "## COLUMN%d.units  %s.m^3\n", n, master->trinfo_3d[tn].units);
      fprintf(region->fp, "## COLUMN%d.missing_value -999\n",n);
      fprintf(region->fp, "##\n");
      n++;
      fprintf(region->fp, "## COLUMN%d.name  end_mass_%s\n", n, master->trinfo_3d[tn].name);
      fprintf(region->fp, "## COLUMN%d.long_name  Tracer %s end mass\n", n,
	      master->trinfo_3d[tn].name);
      fprintf(region->fp, "## COLUMN%d.units  %s.m^3\n", n, master->trinfo_3d[tn].units);
      fprintf(region->fp, "## COLUMN%d.missing_value -999\n",n);
      fprintf(region->fp, "##\n");
      n++;
      fprintf(region->fp, "## COLUMN%d.name  mean_mass_%s\n", n, master->trinfo_3d[tn].name);
      fprintf(region->fp, "## COLUMN%d.long_name  Tracer %s mean mass\n", n,
	      master->trinfo_3d[tn].name);
      fprintf(region->fp, "## COLUMN%d.units  %s.m^3\n", n, master->trinfo_3d[tn].units);
      fprintf(region->fp, "## COLUMN%d.missing_value -999\n",n);
      fprintf(region->fp, "##\n");
      n++;
      fprintf(region->fp, "## COLUMN%d.name  stdev_mass_%s\n", n, master->trinfo_3d[tn].name);
      fprintf(region->fp, "## COLUMN%d.long_name  Tracer %s standard deviation\n", n,
	      master->trinfo_3d[tn].name);
      fprintf(region->fp, "## COLUMN%d.units  %s.m^3\n", n, master->trinfo_3d[tn].units);
      fprintf(region->fp, "## COLUMN%d.missing_value -999\n",n);
      fprintf(region->fp, "##\n");
      n++;
      fprintf(region->fp, "## COLUMN%d.name  mass_budget_%s\n", n, master->trinfo_3d[tn].name);
      fprintf(region->fp, "## COLUMN%d.long_name  Tracer %s mass balance\n", n,
	      master->trinfo_3d[tn].name);
      fprintf(region->fp, "## COLUMN%d.units  %s.m^3\n", n, master->trinfo_3d[tn].units);
      fprintf(region->fp, "## COLUMN%d.missing_value -999\n",n);
      fprintf(region->fp, "##\n");
      n++;
      if (region->mode & RG_VERR) {
	fprintf(region->fp, "## COLUMN%d.name  mass_error_%s\n", n, master->trinfo_3d[tn].name);
	fprintf(region->fp, "## COLUMN%d.long_name  Tracer %s mass error\n", n,
		master->trinfo_3d[tn].name);
	fprintf(region->fp, "## COLUMN%d.units  %s.m^3\n", n, master->trinfo_3d[tn].units);
	fprintf(region->fp, "## COLUMN%d.missing_value -999\n",n);
	fprintf(region->fp, "##\n");
	n++;
      }
      if (region->mode & RG_GLOB) {
	fprintf(region->fp, "## COLUMN%d.name  global_fill_%s\n", n, master->trinfo_3d[tn].name);
	fprintf(region->fp, "## COLUMN%d.long_name  Tracer %s global fill error\n", n,
		master->trinfo_3d[tn].name);
	fprintf(region->fp, "## COLUMN%d.units  %s.m^3\n", n, master->trinfo_3d[tn].units);
	fprintf(region->fp, "## COLUMN%d.missing_value -999\n",n);
	fprintf(region->fp, "##\n");
	n++;
      }
      if (region->mode & RG_PSS) {
	fprintf(region->fp, "## COLUMN%d.name  pss_flux_%s\n", n, master->trinfo_3d[tn].name);
	fprintf(region->fp, "## COLUMN%d.long_name  Point source %s input\n", n,
		master->trinfo_3d[tn].name);
	fprintf(region->fp, "## COLUMN%d.units  %s.m^3\n", n, master->trinfo_3d[tn].units);
	fprintf(region->fp, "## COLUMN%d.missing_value -999\n",n);
	fprintf(region->fp, "##\n");
	n++;
      }
      if (region->mode & RG_SED) {
	fprintf(region->fp, "## COLUMN%d.name  sed_contribution_%s\n", n, master->trinfo_3d[tn].name);
	fprintf(region->fp, "## COLUMN%d.long_name  %s sediment transport contribution\n", n,
		master->trinfo_3d[tn].name);
	fprintf(region->fp, "## COLUMN%d.units  %s.m^3\n", n, master->trinfo_3d[tn].units);
	fprintf(region->fp, "## COLUMN%d.missing_value -999\n",n);
	fprintf(region->fp, "##\n");
	n++;
      }
      if (region->mode & RG_ECO) {
	fprintf(region->fp, "## COLUMN%d.name  eco_contribution_%s\n", n, master->trinfo_3d[tn].name);
	fprintf(region->fp, "## COLUMN%d.long_name  %s ecology contribution\n", n,
		master->trinfo_3d[tn].name);
	fprintf(region->fp, "## COLUMN%d.units  %s.m^3\n", n, master->trinfo_3d[tn].units);
	fprintf(region->fp, "## COLUMN%d.missing_value -999\n",n);
	fprintf(region->fp, "##\n");
	n++;
      }
      for (bb = 0; bb < region->nboundaries; bb++) {
	bn = region->bmap[bb];
	fprintf(region->fp, "## COLUMN%d.name  %s_mass_flux%d-%d\n", n, 
		master->trinfo_3d[tn].name, region->id, bn);
	fprintf(region->fp, "## COLUMN%d.long_name  Tracer %s mass %s\n", n, 
		master->trinfo_3d[tn].name, region->fluxname[bb]);
	fprintf(region->fp, "## COLUMN%d.units  %s.m^3\n", n,
		master->trinfo_3d[tn].units);
	fprintf(region->fp, "## COLUMN%d.missing_value -999\n",n);
	fprintf(region->fp, "##\n");
	n++;
      }
    }
    /* If there was a file prior to restart, write data up to current time */
    if (ots != NULL) {
      int r=0;
      int v=0;
      for (r=0; r<ots->nt; ++r) {
	double newt = ots->t[r];
	tm_change_time_units(master->output_tunit, master->timeunit, &newt, 1);
	if (newt >= master->t) break;
	for (v=0; v<ots->nv; ++v)
	  fprintf(region->fp, "%f ", ots->df->variables[v].data[r]);
	fprintf(region->fp, "\n");
      }
      free(ots);
      ots = NULL;
    }
  }
  time = region->pt;
  tm_change_time_units(master->timeunit, master->output_tunit, &time, 1);
  fprintf(region->fp, "%f ",time);
  fprintf(region->fp, "%f %f ",region->residence, region->residencen);
  for (tt = 0; tt < region->nvar; tt++) {
    d1 = region->smass[tt] - region->emass[tt];
    for (bb = 0; bb < region->nboundaries; bb++)
      d1 += region->trflux[tt][bb];
    fprintf(region->fp, "%f %f %f %f %f ",region->smass[tt], region->emass[tt], 
	    region->mmass[tt], region->stdev[tt], d1);
    if (region->mode & RG_VERR)
      fprintf(region->fp, "%f ",region->merr[tt]);
    if (tt && region->mode & RG_GLOB)
      fprintf(region->fp, "%f ",region->gfer[tt]);
    if (region->mode & RG_PSS)
      fprintf(region->fp, "%f ",region->pssf[tt]);
    if (tt && region->mode & RG_SED)
      fprintf(region->fp, "%f ",region->sedtr[tt]);
    if (tt && region->mode & RG_ECO)
      fprintf(region->fp, "%f ",region->eco[tt]);
    for (bb = 0; bb < region->nboundaries; bb++)
      fprintf(region->fp, "%f ",region->trflux[tt][bb]);
  }
  fprintf(region->fp, "\n");
  fflush(region->fp);
}

/* End region_write()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Get the total flux from a region                                  */
/*-------------------------------------------------------------------*/
double get_totflux(region_t *region, int tr)
{
  int bb;
  double sum = 0;

  for (bb = 0; bb < region->nboundaries; bb++) {
    sum += region->trflux[tr][bb];
  }
  return(sum);
}

/* END get_totflux()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Volume fluxes. The first computation in the block of region->dt   */
/* the fluxes are not computed.                                      */
/*-------------------------------------------------------------------*/
void region_volume_flux_trans(geometry_t *window, window_t *windat, win_priv_t *wincon)
{
  int c, cs, cc, c1, c2, j;
  int dr, sr, m;
  int *S = wincon->s2;
  double **A= wincon->wgt;
  double *obc = wincon->w5;
  double dt = windat->dttr;
  double cellvol, dvol;
  region_t **region=window->region;

  if (window->nregions == 0) return;

  /* Schedule and get the volume error                               */
  for (m = 0; m < window->nregions; m++) {

    region_schedule(region[m], dt, dt, wincon->timeunit);

    /* Set the dz array                                              */
    for (cc = 1; cc <= region[m]->nvec; cc++) {
      c = region[m]->vec[cc];
      region[m]->dz[cc] = (wincon->c1[c]) ? wincon->dz[c]: 0.0;
    }
    if (region[m]->mode & RG_VERR) {
      for (cc = 1; cc <= region[m]->nvec; cc++) {
	c = region[m]->vec[cc];
	region[m]->merr[0] += windat->Vi[c];
      }
    }
    if (region[m]->mode & RG_PSS) {
      for (cc = 1; cc <= region[m]->nvec; cc++) {
	c = region[m]->vec[cc];
	region[m]->pssf[0] += windat->waterss[c] * dt;
      }
    }
  }

  /* Open boundary mask                                              */
  memset(obc, 0, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->v3_t; cc++) {
    c = window->w3_t[cc];
    obc[c] = 1.0;
  }

  /*-----------------------------------------------------------------*/
  /* Sum (aijVi(t))                                                  */
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    c2 = window->m2d[c];
    cs = S[cc];
    cellvol = window->cellarea[c2] * wincon->dz[c];
    dr = master->regionid[c];
    if (dr == NOTVALID) continue;
    for (j = 0; j < 8; j++) {
      if (A[c][j]) {
	c1 = wincon->lmap[cs][j];
	sr = master->regionid[c1];
	if (sr != NOTVALID && dr != sr) {
	  m = region[dr]->rmap[sr];
	  dvol = (A[c][j] * cellvol * obc[c]);
	  if (m != NOTVALID) {
	    region[dr]->trflux[0][m] += dvol;
	    region[dr]->vfluxp += dvol;
	  }
	  m = region[sr]->rmap[dr];
	  if (m != NOTVALID) {
	    region[sr]->trflux[0][m] -= dvol;
	    region[sr]->vfluxn -= dvol;
	  }
	}
      }
    }
  }
}

/* END region_volume_flux_trans()                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Mass fluxes for the transport model                               */
/*-------------------------------------------------------------------*/
void region_flux_trans(geometry_t *window,  /* Window geometry       */
		       window_t *windat,    /* Window data           */
		       win_priv_t *wincon,  /* Window constants      */
		       int trn,             /* Tracer number         */
		       double *dtracer      /* pss flux              */
		       )
{
  int c, cs, cc, j, c1, c2;
  int tt, tn, dr, sr, m;
  int *S = wincon->s2;
  double **A= wincon->wgt;
  double *obc = wincon->w5;
  double *tr = windat->tr_wc[trn];
  double cellvol;
  region_t **region=window->region;
  int found;

  if (window->nregions == 0) return;

  /* Get the volume error                                            */
  for (m = 0; m < window->nregions; m++) {
    if ((tt = find_trindex(region[m], trn))) {
      if (region[m]->mode & RG_VERR) {
	for (cc = 1; cc <= region[m]->nvec; cc++) {
	  c = region[m]->vec[cc];
	  region[m]->merr[tt] += (tr[c] * windat->Vi[c]);
	}
      }
      if (region[m]->mode & RG_PSS) {
	for (cc = 1; cc <= region[m]->nvec; cc++) {
	  c = region[m]->vec[cc];
	  region[m]->pssf[tt] -= (dtracer[c] * (double)wincon->c1[c]);
	}
      }
    }
  }

  /* Open boundary mask                                              */
  memset(obc, 0, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->v3_t; cc++) {
    c = window->w3_t[cc];
    obc[c] = 1.0;
  }

  /*-----------------------------------------------------------------*/
  /* Sum (aijVi(t))                                                  */
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    c2 = window->m2d[c];
    cs = S[cc];
    cellvol = window->cellarea[c2] * wincon->dz[c];
    dr = master->regionid[c];
    if (dr == NOTVALID || dr < 0 || dr > window->nregions) continue;
    for (j = 0; j < 8; j++) {
      if (A[c][j]) {
	c1 = wincon->lmap[cs][j];
	sr = master->regionid[c1];
	if (sr == NOTVALID || sr < 0 || sr > window->nregions) continue;
	if (sr != NOTVALID && dr != sr) {
	  tt = find_trindex(region[dr], trn);
	  m = region[dr]->rmap[sr];
	  if (tt && m != NOTVALID)
	    region[dr]->trflux[tt][m] += (A[c][j] * cellvol) * obc[c] * tr[c1];
	  tt = find_trindex(region[sr], trn);
	  m = region[sr]->rmap[dr];
	  if (tt && m != NOTVALID)
	    region[sr]->trflux[tt][m] -= (A[c][j] * cellvol) * obc[c] * tr[c1];
	}
      }
    }
  }
}

/* END region flux_trans()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Finds a tracer index in the variable vector                       */
/*-------------------------------------------------------------------*/
int find_trindex(region_t *region, int trn)
{
  int tt, tn;

  for (tt = 1; tt < region->nvar; tt++) {
    tn = region->var[tt];
    if (tn == trn)
      return(tt);
  }
  return(0);
}

/* END find_trindex()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Mass fluxes for the transport model                               */
/*-------------------------------------------------------------------*/
void region_mass_tr(geometry_t *window,  /* Window geometry     */
		    window_t *windat,    /* Window data         */
		    win_priv_t *wincon,  /* Window constants    */
		    double sign,         /* Add or subtract */
		    int trn,             /* Tracer number */
		    int mode             /* Array to increment */
		    )
{
  int m, tt, tn, c, cc, cs;
  double *tr;
  region_t **region=window->region;
  double *vec;

  if (window->nregions == 0) return;

  /* Get the volume error                                            */
  for (m = 0; m < window->nregions; m++) {    
    /* Get the mass array to fill                                    */
    if (mode & RG_GLOB && region[m]->mode & RG_GLOB)
      vec = region[m]->gfer;
    else if (mode & RG_SED && region[m]->mode & RG_SED)
      vec = region[m]->sedtr;
    else if (mode & RG_ECO && region[m]->mode & RG_ECO)
      vec = region[m]->eco;
    else
      continue;
    if (trn == 0) {
      for (cc = 1; cc <= region[m]->nvec; cc++) {
	c = region[m]->vec[cc];
	cs = window->m2d[c];
	vec[0] += sign * window->cellarea[cs] * region[m]->dz[cc];
      }
    } else {
      /* Get the tracer index                                        */
      tr = windat->tr_wc[trn];
      tt = find_trindex(region[m], trn);
      /* Increment the array                                         */
      if (tt) {
	for (cc = 1; cc <= region[m]->nvec; cc++) {
	  c = region[m]->vec[cc];
	  cs = window->m2d[c];
	  vec[tt] += sign * tr[c] * window->cellarea[cs] * region[m]->dz[cc];
	}
      }
    }
  }
}

/* END region_mass_tr()                                              */
/*-------------------------------------------------------------------*/
