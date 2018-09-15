/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/boundary/boundaryio.c
 *  
 *  Description:
 *  Routines to read (and write) lists of specified
 *  boundary conditions from a text file
 *  (usually the parameter file)
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: boundaryio.c 5943 2018-09-13 04:39:09Z her127 $
 *
 */

#include <string.h>
#include <stdio.h>
#include <time.h>
#include "hd.h"

int read_bdry_custom(open_bdrys_t *open, int bnum, char *name,
                     double fill, bdry_details_t *data, char *cusname);
void locate_boundary_function(open_bdrys_t *open, const char *tag,
                              bdrycustom_init_m * init_m,
                              bdrycustom_init_w * init_w,
                              bdrycustom_m * func_m,
                              bdrycustom_w * func_w, bdrytransfer * trans,
			      bdryfree_m * free_m, bdryfree_w * free_w );
tidal_memory_t **tide_alloc_2d(long n1, long n2);
int cbc(int num);
tide_details_t *sreadtideforce(FILE * fp, int bnum, int n);
int get_bcond_no(char *list);
int bcond_no(char *list);
void scale_data_copy(scale_details_t *data, scale_details_t *din);
void read_bdry_scale(open_bdrys_t *io, open_bdrys_t *open, char *trname[]);
void read_bdry_eta_scale(master_t *master, open_bdrys_t *io, open_bdrys_t *open);
void get_relax_time(open_bdrys_t *open, FILE *fp, char *bname, int n, 
		    int obc1, int obc2);
void custom_free_m(master_t *master, bdry_details_t *data);
void custom_free_w(geometry_t *window, bdry_details_t *data);
void bdry_custom_init(master_t *master, open_bdrys_t *open, tracer_info_t *tracer, int type);
void std_bdry_init(open_bdrys_t *open, char *stdbdry,
		   double dt, double dt2d, char *bdrypath, tracer_info_t *tracer);


/*-------------------------------------------------------------------*/
/* Routine to read the open boundary condition types and save to the */
/* open boundary data structure.                                     */
/* Note: any parameters read in here must be copied to the master in */
/* copy_OBC_conds() and from the master to each window open boundary */
/* in OBC_build() in windows.c.                                      */
/*-------------------------------------------------------------------*/
void get_OBC_conds(parameters_t *params,   /*      Input parameters        */
		   open_bdrys_t *open,     /* Open boundary data structure */
		   FILE *fp,        /* File handle to input parameter file */
		   int n,                         /* Open boundary number  */
		   tracer_info_t *tracers         /* Tracer info structure */
		   )
{
  char buf[MAXSTRLEN];          /* Dummy string for text */
  int i;                        /* Counter */
  char keyword[MAXSTRLEN];      /* Buffer for I/O */
  char key2[MAXSTRLEN];         /* Buffer for I/O */
  char bname[MAXSTRLEN];        /* Boundary condation name */
  int bdf = 0;                  /* Check for custom routine success */
  double d1, d2, d3, d4;        /* Dummies */

  prm_set_errfn(hd_quit);

  /* Get the open boundary condition type */
  open->id = n;
  sprintf(keyword, "BOUNDARY%1d.TYPE", n);
  prm_read_char(fp, keyword, buf);
  if (strcasecmp(buf, "u1") == 0) {
    open->type |= U1BDRY;
  } else if (strcasecmp(buf, "u2") == 0) {
    open->type |= U2BDRY;
  } else
    hd_quit
      ("Unsupported boundary type '%s'.\nMust be 'u1', 'u2' or 'velocity'.\n",
       buf);

  /* Get the name of the open boundary */
  prm_set_errfn(hd_silent_warn);
  sprintf(keyword, "BOUNDARY%1d.NAME", n);
  prm_read_char(fp, keyword, open->name);

  /*-----------------------------------------------------------------*/
  /* Initialise and allocate memory */
  init_OBC_conds(params, open);

  /*-----------------------------------------------------------------*/
  /* Read general boundary information */
  sprintf(keyword, "BOUNDARY%1d.NSPONGE_VERT", n);
  prm_read_int(fp, keyword, &open->sponge_zone);
  sprintf(keyword, "BOUNDARY%1d.NSPONGE_HORZ", n);
  prm_read_double(fp, keyword, &open->sponge_zone_h);
  sprintf(keyword, "BOUNDARY%1d.SPONGE_FACT", n);
  prm_read_double(fp, keyword, &open->sponge_f);
  sprintf(keyword, "BOUNDARY%1d.SMOOTH_PHASE", n);
  prm_read_double(fp, keyword, &open->spf);
  sprintf(keyword, "BOUNDARY%1d.BATHY_CON", n);
  prm_read_int(fp, keyword, &open->bathycon);
  sprintf(keyword, "BOUNDARY%1d.BATHY_SMOOTH", n);
  if (prm_read_char(fp, keyword, buf))
    sscanf(buf,"%d %d", &open->smooth_z, &open->smooth_n);
  sprintf(keyword, "BOUNDARY%1d.OUTSIDE_ZONE", n);
  prm_read_int(fp, keyword, &open->bout);
  sprintf(keyword, "BOUNDARY%1d.ADJUST_FLUX", n);
  if (prm_read_char(fp, keyword, buf) > 0)
    tm_scale_to_secs(buf, &open->adjust_flux);
  sprintf(keyword, "BOUNDARY%1d.ADJUST_TIDE", n);
  if (prm_read_char(fp, keyword, buf) > 0)
    tm_scale_to_secs(buf, &open->adjust_flux_s);
  if(strlen(params->patm)) {
    sprintf(keyword, "BOUNDARY%1d.INVERSE_BAROMETER", n);
    if (prm_read_char(fp, keyword, buf))
      open->inverse_barometer = is_true(buf);
    sprintf(keyword, "BOUNDARY%1d.INV_BAR", n);
    if (prm_read_char(fp, keyword, buf))
      open->inverse_barometer = is_true(buf);
  } else
    open->inverse_barometer = 0;
  sprintf(keyword, "BOUNDARY%1d.RELAX_ELE", n);
  if (prm_read_char(fp, keyword, buf)) {
    if (sscanf(buf, "%lf %lf %lf", &d1, &d2, &d3) == 3) {
      open->relax_ele = (int)d1; 
      d4 = params->grid_dt / (double)params->iratio;
      open->rele_b = d2 * d4;
      open->rele_i = d3 * d4;
      if (open->rele_b == 0.0) open->rele_b = d4;
      if (open->rele_i == 0.0) open->rele_i = d4;
    } else
      hd_warn("Boundary %d elevation relaxation format: RELAX_ELE zone boundary_tc interior_tc\n",n);
  }
  if (params->runmode & (MANUAL | RE_ROAM | TRANS | DUMP)) {
    sprintf(keyword, "BOUNDARY%1d.STAGGER", n);
    prm_read_char(fp, keyword, buf);
    if (strcmp(buf,"INFACE") == 0) open->stagger = INFACE|INFACET;
  }

  /*-----------------------------------------------------------------*/
  /* Standard conditions */
  for (i = 0; i < NBSTD; i++) {
    sprintf(keyword, "BOUNDARY%1d.BCOND%1d", n, i);
    if (prm_read_char(fp, keyword, open->bstd[i]))
      open->nbstd++;
  }

  /*-----------------------------------------------------------------*/
  /* Normal velocity components */
  open->bcond_nor = NOTHIN;
  sprintf(keyword, "BOUNDARY%1d.BCOND_NOR", n);
  if (prm_read_char(fp, keyword, buf)) {
    open->bcond_nor = get_bcond_no(buf);
    bcname(open->bcond_nor, bname);
    if (open->bcond_nor <= 0 || open->bcond_nor > maxbc) {
      hd_warn("Boundary %d: Normal boundary condition type unspecified.\n", n);
      open->bcond_nor = FILEIN;
    } else if (open->bcond_nor & (UPSTRM | TIDALM | TIDEBC | TRFLUX | TRCONC | TRCONF)) {
      hd_quit("Boundary %d: Unsupported normal boundary condition type %s.\n", n, bname);
    }
  }
  open->bcond_nor2d = open->bcond_nor;
  if (open->bcond_nor & (FILEIN | CUSTOM))
    open->bcond_nor2d = VERTIN;
  sprintf(buf, "NULL");
  sprintf(keyword, "BOUNDARY%1d.BCOND_NOR2D", n);
  prm_read_char(fp, keyword, buf);
  if (strcmp(buf, "NULL") != 0)
    open->bcond_nor2d = get_bcond_no(buf);
  get_relax_time(open, fp, bname, n, open->bcond_nor, open->bcond_nor2d);
  open->relax_zone_nor = 0;
  open->rnor_b = open->rnor_i = 0.0;
  sprintf(keyword, "BOUNDARY%1d.RELAX_ZONE_NOR", n);
  if (prm_read_char(fp, keyword, buf)) {
    char *fields[MAXSTRLEN * MAXNUMARGS];
    i = parseline(buf, fields, MAXNUMARGS);
    open->relax_zone_nor = atoi(fields[0]);
    if (i == 3) {
      open->rnor_b = atof(fields[1]);
      open->rnor_i = atof(fields[2]);
    }
  }
  open->linear_zone_nor = 0;
  sprintf(keyword, "BOUNDARY%1d.LINEAR_ZONE_NOR", n);
  prm_read_int(fp, keyword, &open->linear_zone_nor);

  /*-----------------------------------------------------------------*/
  /* Tracer components */
  open->tidemem_depth = 0.0;
  sprintf(buf, "NULL");
  /* Set all the automatically generated tracers */
  for (i = 0; i < open->atr; i++) {
    open->bcond_tra[i] = NOGRAD;
    /* Set any auto sediment tracers */
#if defined(HAVE_SEDIMENT_MODULE)
    if (is_sed_var(tracers[i].name))
	open->bcond_tra[i] = sed_get_obc(&tracers[i]);
#endif
#if defined(HAVE_ECOLOGY_MODULE)
    if (is_eco_var(tracers[i].name))
      open->bcond_tra[i] = eco_get_obc(&tracers[i]);
#endif
  }

  /* Set all the tracers */
  sprintf(keyword, "BOUNDARY%1d.BCOND_TRA_ALL", n);
  prm_read_char(fp, keyword, buf);
  if (strcmp(buf, "NULL") != 0) {
    int bc = get_bcond_no(buf);
    for (i = open->atr; i < open->ntr; i++) {
      open->bcond_tra[i] = bc;
      /* Set the default UPSTRM condition to use file data */
      if (open->bcond_tra[i] == UPSTRM)
        open->bcond_tra[i] = UPSTRM | FILEIN;
      if (open->bcond_tra[i] == TRFLUX)
        open->bcond_tra[i] = TRFLUX | FILEIN;
      if (open->bcond_tra[i] == TRCONC)
        open->bcond_tra[i] = TRCONC | FILEIN;
      if (open->bcond_tra[i] == TRCONF)
        open->bcond_tra[i] = TRCONF | FILEIN;
      open->relax_zone_tra[i] = 0;
      open->rtra_b[i] = 0;
      open->rtra_i[i] = 0;
      sprintf(keyword, "BOUNDARY%1d.RELAX_ZONE_ALL", n);
      if (prm_read_char(fp, keyword, buf)) {
	char *fields[MAXSTRLEN * MAXNUMARGS];
	i = parseline(buf, fields, MAXNUMARGS);
	open->relax_zone_tra[i] = atoi(fields[0]);
	if (i == 3) {
	  open->rtra_b[i] = atof(fields[1]);
	  open->rtra_i[i] = atof(fields[2]);
	}
      }
      sprintf(keyword, "BOUNDARY%1d.TRPC_ALL", n);
      prm_read_double(fp, keyword, &open->trpc[i]);
    }
  } else {
    for (i = open->atr; i < open->ntr; i++) {
      open->bcond_tra[i] = NOGRAD;
      open->relax_zone_tra[i] = 0;
    }
  }

  /* Set the individual tracers by tracer number */
  for (i = open->atr; i < open->ntr; i++) {
    int tm = tracers[i].m;
    sprintf(buf, "NULL");
    sprintf(keyword, "BOUNDARY%1d.BCOND_TRA%d", n, tm);
    prm_read_char(fp, keyword, buf);
    if (strcmp(buf, "NULL") != 0) {
      open->bcond_tra[i] = get_bcond_no(buf);
      /* Set the default UPSTRM condition to use file data */
      if (open->bcond_tra[i] == UPSTRM)
        open->bcond_tra[i] = UPSTRM | FILEIN;
      if (open->bcond_tra[i] == TRFLUX)
        open->bcond_tra[i] = TRFLUX | FILEIN;
      if (open->bcond_tra[i] == TRCONC)
        open->bcond_tra[i] = TRCONC | FILEIN;
      if (open->bcond_tra[i] == TRCONF)
        open->bcond_tra[i] = TRCONF | FILEIN;
      if (!(open->bcond_tra[i] & (CLAMPD | CYCLIC | CYCLED | UPSTRM | FILEIN | TRCONF |
				  CUSTOM | PROFIL | DEPROF | LINEXT | POLEXT | DESCAL |
				  NOGRAD | TIDALM | STATIS | TRFLUX | TRCONC | NOTHIN)))
	hd_quit("Boundary %d: Unsupported tracer%d boundary condition type %s.\n", 
		n, i, buf);
      open->relax_zone_tra[i] = 0;
      open->rtra_b[i] = 0;
      open->rtra_i[i] = 0;
      sprintf(keyword, "BOUNDARY%1d.RELAX_ZONE_TRA%d", n, tm);
      if (prm_read_char(fp, keyword, buf)) {
	char *fields[MAXSTRLEN * MAXNUMARGS];
	i = parseline(buf, fields, MAXNUMARGS);
	open->relax_zone_tra[i] = atoi(fields[0]);
	if (i == 3) {
	  open->rtra_b[i] = atof(fields[1]);
	  open->rtra_i[i] = atof(fields[2]);
	}
      }
      sprintf(keyword, "BOUNDARY%1d.TRPC_TRA%d", n, tm);
      prm_read_double(fp, keyword, &open->trpc[i]);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the individual tracers by tracer name. Look for autotracers */
  /* also.                                                           */
  /*for (i = open->atr; i < open->ntr; i++) {*/
  for (i = 0; i < open->ntr; i++) {
    int tm = tracers[i].m;
    sprintf(buf, "NULL");
    sprintf(keyword, "BOUNDARY%1d.BCOND_%s", n, tracers[i].name);
    prm_read_char(fp, keyword, buf);
    if (strcmp(buf, "NULL") != 0) {
      open->bcond_tra[i] = get_bcond_no(buf);
      /* Set the default UPSTRM condition to use file data */
      if (open->bcond_tra[i] == UPSTRM)
        open->bcond_tra[i] = UPSTRM | FILEIN;
      if (open->bcond_tra[i] == TRFLUX)
        open->bcond_tra[i] = TRFLUX | FILEIN;
      if (open->bcond_tra[i] == TRCONC)
        open->bcond_tra[i] = TRCONC | FILEIN;
      if (open->bcond_tra[i] == TRCONF)
        open->bcond_tra[i] = TRCONF | FILEIN;
      if (!(open->bcond_tra[i] & (CLAMPD | CYCLIC | CYCLED | UPSTRM | FILEIN | TRCONF |
				  CUSTOM | PROFIL | DEPROF | LINEXT | POLEXT | DESCAL |
				  NOGRAD | TIDALM | STATIS | TRFLUX | TRCONC | NOTHIN)))
	hd_quit("Boundary %d: Unsupported tracer%d boundary condition type %s.\n", 
		n, i, buf);
      open->relax_zone_tra[i] = 0;
      open->rtra_b[i] = 0;
      open->rtra_i[i] = 0;
      sprintf(keyword, "BOUNDARY%1d.RELAX_ZONE_%s", n, tracers[i].name);
      if (prm_read_char(fp, keyword, buf)) {
	char *fields[MAXSTRLEN * MAXNUMARGS];
	i = parseline(buf, fields, MAXNUMARGS);
	open->relax_zone_tra[i] = atoi(fields[0]);
	if (i == 3) {
	  open->rtra_b[i] = atof(fields[1]);
	  open->rtra_i[i] = atof(fields[2]);
	}
      }
      sprintf(keyword, "BOUNDARY%1d.TRPC_%s", n, tracers[i].name);
      prm_read_double(fp, keyword, &open->trpc[i]);
    }
  }

  for (i = 0; i < open->ntr; i++) {
    bcname(open->bcond_tra[i], bname);
    if (open->bcond_tra[i] < 0 || open->bcond_tra[i] > maxbc) {
      hd_warn("Boundary %d: Tracer%d boundary condition type unspecified.\n", n, i);
      open->bcond_tra[i] = FILEIN;
      open->relax_zone_tra[i] = 0;
    }

    /* Get the depth above which to impose tidal memory */
    if (open->bcond_tra[i] & TIDALM) {
      sprintf(keyword, "BOUNDARY%1d.TM_DEPTH", n);
      prm_read_double(fp, keyword, &open->tidemem_depth);
    }
    /* Get the fill value for clamped open boundary conditions */
    if (open->bcond_tra[i] & CLAMPD) {
      sprintf(keyword, "BOUNDARY%1d.CLAMP_VAL", n);
      if (!prm_read_double(fp, keyword, &open->clampv[i]))
	open->clampv[i] = tracers[i].fill_value_wc;
    }
  }

  /* Set the velocity face / center to be used with UPSTRM           */
  open->upmeth = INTERIOR;
  sprintf(keyword, "BOUNDARY%1d.UPSTRM_METHOD", n);
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "CENTER") == 0)
      open->upmeth = CENTER;
    if (strcmp(buf, "FACE") == 0)
      open->upmeth = FACE;
    if (strcmp(buf, "INTERIOR") == 0)
      open->upmeth = INTERIOR;
    if (strcmp(buf, "ADAPTIVE") == 0)
      open->upmeth = ADAPTIVE;
  }

  /*-----------------------------------------------------------------*/
  /* Elevation */
  open->bcond_ele = NOTHIN;
  sprintf(keyword, "BOUNDARY%1d.BCOND_ELE", n);
  if (prm_read_char(fp, keyword, buf)) {
    open->bcond_ele = get_bcond_no(buf);
    bcname(open->bcond_ele, bname);
    if (open->bcond_ele <= 0 || open->bcond_ele > maxbc) {
      hd_warn("Boundary %d: Elevation boundary condition type unspecified.\n", n);
      open->bcond_ele = NOTHIN;
    } else if (open->bcond_ele & (UPSTRM | TIDALM | VERTIN | TRFLUX | TRCONC | TRCONF)) {
      hd_quit("Boundary %d: Unsupported elevation boundary condition type %s.\n", n, bname);
    }
  }
  if (open->adjust_flux == 0.0 && open->nbstd == 0)
    get_relax_time(open, fp, bname, n, open->bcond_ele, NOTHIN);
  if (open->bcond_ele & TIDEBC) {
    prm_set_errfn(hd_quit);
    sprintf(keyword, "BOUNDARY%1d.T_CONSTITUENTS", n);
    prm_read_int(fp, keyword, &open->ntide);
    prm_flush_line(fp);
    open->tideforce = sreadtideforce(fp, n, open->ntide);
    prm_set_errfn(hd_silent_warn);
  }
  if(open->bcond_ele & TIDALC ||
     open->bcond_nor & TIDALC ||
     open->bcond_nor2d & TIDALC ||
     open->bcond_tan & TIDALC ||
     open->bcond_tan2d & TIDALC) {
    prm_set_errfn(hd_quit);
    sprintf(keyword, "BOUNDARY%1d.T_CONSTITUENTS", n);
    prm_read_char(fp, keyword, open->tide_con);
    prm_set_errfn(hd_silent_warn);
  }
  if (open->bcond_ele & FLATHR) {
    open->bcond_ele &= ~FLATHR;
    open->bcond_ele |= FLATHE;
  }
  if (!(params->compatible & V1670) && open->adjust_flux && 
      !(open->bcond_ele & (FILEIN|TIDALC)) && !strlen(open->bstd[0]))
	hd_quit("Boundary%d elevation DATA must be supplied using ADJUST_FLUX.\n", open->id);

  open->relax_zone_ele = 0;
  sprintf(keyword, "BOUNDARY%1d.RELAX_ZONE_ELE", n);
  prm_read_int(fp, keyword, &open->relax_zone_ele);
  if (open->bcond_ele & FILEIN || strlen(open->bstd[0])) {
    sprintf(keyword, "BOUNDARY%1d.FILEIN_DT", n);
    if (prm_read_char(fp, keyword, buf))
      tm_scale_to_secs(buf, &open->file_dt);
  }

  /*-----------------------------------------------------------------*/
  /* Tangential velocity components */
  open->bcond_tan = NOTHIN;
  sprintf(keyword, "BOUNDARY%1d.BCOND_TAN", n);
  if (prm_read_char(fp, keyword, buf)) {
    open->bcond_tan = get_bcond_no(buf);
    bcname(open->bcond_tan, bname);
    if (open->bcond_tan <= 0 || open->bcond_tan > maxbc) {
      hd_warn("Boundary %d: Tangential boundary condition type unspecified.\n", n);
      open->bcond_tan2d = open->bcond_tan = CLAMPD;
    } else if (open->bcond_tan & (UPSTRM | TIDALM | TIDEBC | TRFLUX | TRCONC | TRCONF))
      hd_quit("Boundary %d: Unsupported tangential boundary condition type %s.\n", n, bname);
  }
  open->bcond_tan2d = open->bcond_tan;
  if (open->bcond_tan & (FILEIN | CUSTOM))
    open->bcond_tan2d = VERTIN;
  sprintf(buf, "NULL");
  sprintf(keyword, "BOUNDARY%1d.BCOND_TAN2D", n);
  prm_read_char(fp, keyword, buf);
  if (strcmp(buf, "NULL") != 0)
    open->bcond_tan2d = get_bcond_no(buf);
  get_relax_time(open, fp, bname, n, open->bcond_tan, open->bcond_tan2d);
  open->relax_zone_tan = 0;
  open->rnor_b = open->rnor_i = 0.0;
  sprintf(keyword, "BOUNDARY%1d.RELAX_ZONE_TAN", n);
  if (prm_read_char(fp, keyword, buf)) {
    char *fields[MAXSTRLEN * MAXNUMARGS];
    i = parseline(buf, fields, MAXNUMARGS);
    open->relax_zone_tan = atoi(fields[0]);
    if (i == 3) {
      open->rtan_b = atof(fields[1]);
      open->rtan_i = atof(fields[2]);
    }
  }
  open->linear_zone_tan = 0;
  sprintf(keyword, "BOUNDARY%1d.LINEAR_ZONE_TAN", n);
  prm_read_int(fp, keyword, &open->linear_zone_tan);

  if (open->adjust_flux == 0.0) {
    if ((open->bcond_nor & (FILEIN | CUSTOM)) && 
	open->bcond_ele & (FILEIN | CUSTOM))
      hd_warn("Multiple forcing (nor. velocity and eta) specified for %s.\n",
	      open->name);
    if (open->bcond_tan & (FILEIN | CUSTOM) &&
	open->bcond_ele & (FILEIN | CUSTOM))
      hd_warn("Multiple forcing (tan. velocity and eta) specified for %s.\n",
	      open->name);
  }

  /*-----------------------------------------------------------------*/
  /* Vertical velocity coefficients */
  open->bcond_w = NOTHIN;
  if (params->runmode & TRANS) open->bcond_w = NOGRAD;
  sprintf(keyword, "BOUNDARY%1d.BCOND_W", n);
  if (prm_read_char(fp, keyword, buf)) {
    open->bcond_w = get_bcond_no(buf);
    bcname(open->bcond_w, bname);
  }
  if (open->bcond_w <= 0 || open->bcond_w > maxbc) {
    hd_warn("Boundary %d: Vertical velocity boundary condition type unspecified.\n", n);
    open->bcond_w = NOTHIN;
  } else if (open->bcond_w & (UPSTRM | TIDALM | TIDEBC | TRFLUX | TRCONC | TRCONF | FILEIN | CUSTOM)) {
    hd_quit("Boundary %d: Unsupported vertical velocity boundary condition type %s.\n", n, bname);
  }

  /*-----------------------------------------------------------------*/
  /* Vertical mixing coefficients */
  open->bcond_Vz = open->bcond_Kz = NOTHIN;
  sprintf(keyword, "BOUNDARY%1d.BCOND_VZ", n);
  if (prm_read_char(fp, keyword, buf)) {
    open->bcond_Vz = get_bcond_no(buf);
    bcname(open->bcond_Vz, bname);
  }
  if (open->bcond_Vz <= 0 || open->bcond_Vz > maxbc) {
    hd_warn("Boundary %d: Vertical viscosity boundary condition type unspecified.\n", n);
    open->bcond_Vz = NOTHIN;
  } else
    if (!
        (open->
         bcond_Vz & (CLAMPD | CYCLIC | CYCLED | UPSTRM | FILEIN | CUSTOM | 
                     LINEXT | POLEXT | NOGRAD | NOTHIN | TRFLUX | TRCONC | TRCONF))) {
      hd_quit("Boundary %d: Unsupported vertical viscosity boundary condition type %s.\n", n, bname);
  }
  sprintf(keyword, "BOUNDARY%1d.BCOND_KZ", n);
  if (prm_read_char(fp, keyword, buf)) {
    open->bcond_Kz = get_bcond_no(buf);
    bcname(open->bcond_Vz, bname);
  }
  if (open->bcond_Kz <= 0 || open->bcond_Kz > maxbc) {
    hd_warn("Boundary %d: Vertical diffusion boundary condition type unspecified.\n", n);
    open->bcond_Kz = NOTHIN;
  } else
    if (!
        (open->
         bcond_Kz & (CLAMPD | CYCLIC | CYCLED | UPSTRM | FILEIN | CUSTOM | 
                     LINEXT | POLEXT | NOGRAD | NOTHIN | TRFLUX | TRCONC | TRCONF))) {
      hd_quit("Boundary %d: Unsupported vertical diffusion boundary condition type %s.\n", n, bname);
  }

  prm_set_errfn(hd_silent_warn);

  /*-----------------------------------------------------------------*/
  /* Read boundary forcing data if present */
  sprintf(open->tsfn, "%c", '\0');
  sprintf(keyword, "BOUNDARY%1d.DATA", n);
  prm_read_char(fp, keyword, open->tsfn);
  if (strlen(params->bdrypath) && strlen(open->tsfn)) {
    sprintf(buf, "%s%s", params->bdrypath, open->tsfn);
    strcpy(open->tsfn, buf);
  }

  /*-----------------------------------------------------------------*/
  /* Check if a custom boundary condition has been specified.        */
  /* Normal velocity custom routines are associated with cusname_u1  */
  /* and tangential custom routines with cusname_u2.                 */
  sprintf(open->cusname_u1, "%c", '\0');
  sprintf(open->cusname_u2, "%c", '\0');
  sprintf(open->cusname_u1av, "%c", '\0');
  sprintf(open->cusname_u2av, "%c", '\0');
  sprintf(open->cusname_eta, "%c", '\0');
  open->custype = U1BDRY;
  sprintf(keyword, "BOUNDARY%1d.CUSTOM.%s", n, "eta");
  sprintf(keyword, "BOUNDARY%1d.CUSTOM.%s", n, "eta");
  if (prm_read_char(fp, keyword, open->cusname_eta))
    bdf = 1;
  sprintf(keyword, "BOUNDARY%1d.CUSTOM.%s", n, "u1");
  sprintf(key2, "BOUNDARY%1d.CUSTOM.%s", n, "nor");
  if ((prm_read_char(fp, keyword, open->cusname_u1) ||
       prm_read_char(fp, key2, open->cusname_u1)) && bdf)
    hd_warn
      ("Multiple custom routines (eta and nor) specified for BOUNDARY%1d.\n",
       n);
  sprintf(keyword, "BOUNDARY%1d.CUSTOM.%s", n, "u2");
  sprintf(key2, "BOUNDARY%1d.CUSTOM.%s", n, "tan");
  if ((prm_read_char(fp, keyword, open->cusname_u2) ||
       prm_read_char(fp, key2, open->cusname_u2)) && bdf)
    hd_quit
      ("Multiple custom routines (eta and tan) specified for BOUNDARY%1d.\n",
       n);
  
  sprintf(keyword, "BOUNDARY%1d.CUSTOM.%s", n, "u1av");
  prm_read_char(fp, keyword, open->cusname_u1av);
  sprintf(keyword, "BOUNDARY%1d.CUSTOM.%s", n, "nor2d");
  prm_read_char(fp, keyword, open->cusname_u1av);
  sprintf(keyword, "BOUNDARY%1d.CUSTOM.%s", n, "u2av");
  prm_read_char(fp, keyword, open->cusname_u2av);
  sprintf(keyword, "BOUNDARY%1d.CUSTOM.%s", n, "tan2d");
  prm_read_char(fp, keyword, open->cusname_u2av);
  if (strcmp(open->cusname_u1, "u1flowbdry") == 0 || strlen(open->bstd[0])) {
    sprintf(keyword, "BOUNDARY%1d.U1_FLOW", n);
    prm_read_char(fp, keyword, open->bflow);
    sprintf(keyword, "BOUNDARY%1d.U1_HC", n);
    prm_read_double(fp, keyword, &open->bhc);
    sprintf(keyword, "BOUNDARY%1d.U1_LENGTH", n);
    prm_read_double(fp, keyword, &open->rlen);
  }

  /* eta scaling */
  sprintf(open->scale_e, "%c", '\0');
  sprintf(keyword, "BOUNDARY%1d.SCALE_ETA", n);
  prm_read_char(fp, keyword, open->scale_e);

  /* Tracers */
  for (i = 0; i < open->ntr; i++) {
    open->cusname_t[i] = (char *)malloc(MAXSTRLEN);
    sprintf(open->cusname_t[i], "%c", '\0');
    sprintf(keyword, "BOUNDARY%1d.CUSTOM.%s", n, tracers[i].name);
    prm_read_char(fp, keyword, open->cusname_t[i]);

    open->scale_s[i] = (char *)malloc(MAXSTRLEN);
    sprintf(open->scale_s[i], "%c", '\0');
    sprintf(keyword, "BOUNDARY%1d.SCALE_S.%s", n, tracers[i].name);
    prm_read_char(fp, keyword, open->scale_s[i]);

    open->scale_p[i] = (char *)malloc(MAXSTRLEN);
    sprintf(open->scale_p[i], "%c", '\0');
    sprintf(keyword, "BOUNDARY%1d.SCALE_P.%s", n, tracers[i].name);
    prm_read_char(fp, keyword, open->scale_p[i]);

    open->scale_d[i] = (char *)malloc(MAXSTRLEN);
    sprintf(open->scale_d[i], "%c", '\0');
    sprintf(keyword, "BOUNDARY%1d.SCALE_D.%s", n, tracers[i].name);
    prm_read_char(fp, keyword, open->scale_d[i]);
  }

  /*-----------------------------------------------------------------*/
  /* Options                                                         */
  sprintf(keyword, "BOUNDARY%1d.OPTIONS", n);
  if (prm_read_char(fp, keyword, buf)) {
    if (contains_token(buf, "NONE") != NULL) {
      open->options = NONE;
    } else {
      open->options = 0;
      if (contains_token(buf, "UPSTRM") != NULL)
	open->options |= OP_UPSTRM;
      if (contains_token(buf, "GEOSTR") != NULL)
	open->options |= OP_GEOSTR;
      if (contains_token(buf, "RLEN") != NULL)
	open->options |= OP_RLEN;
      if (contains_token(buf, "DYNAMIC_HC") != NULL)
	open->options |= OP_DYNAHC;
      if (contains_token(buf, "NO_HDIFF") != NULL)
	open->options |= OP_NOHDIF;
      if (contains_token(buf, "YANKOVSKY") != NULL)
	open->options |= OP_YANKOVSKY;
      if (contains_token(buf, "NO_SALT") != NULL)
	open->options |= OP_NOSALT;
      if (contains_token(buf, "FULL_DEPTH") != NULL)
	open->options |= OP_FDEPTH;
      if (contains_token(buf, "FRESH_FLOW") != NULL)
	open->options |= OP_IFRESH;
      if (contains_token(buf, "TRUNC_LAYER") != NULL)
	open->options |= OP_TRUNCL;
      if (contains_token(buf, "SCALE_MULT") != NULL)
	open->options |= OP_MULTF;
      if (contains_token(buf, "ETA_MEDFIL") != NULL)
	open->options |= OP_ETAFIL;
      if (contains_token(buf, "CORNER_MEANS") != NULL)
	open->options |= OP_OBCCM;
      if (contains_token(buf, "ISO_SPONGE") != NULL)
	open->options |= OP_ISPNG;
      if (contains_token(buf, "NO_OUTFLOW") != NULL)
	open->options |= OP_PARFLOW;
      if (contains_token(buf, "NEST_BARO") != NULL)
	open->options |= OP_NESTBARO;
      if (contains_token(buf, "MACREADY") != NULL)
	open->options |= OP_MACREADY;
      if (contains_token(buf, "OVERWRITE") != NULL)
	open->options |= OP_OWRITE;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Standard conditions */
  if (strlen(open->bstd[0])) {
    open->bstdf = 0;
    std_bdry_init(open, open->bstd[0], params->grid_dt, 
		  params->grid_dt / (double)params->iratio,
		  params->bdrypath, params->trinfo_3d);
		  
  }
  if (params->runmode & DUMP && open->stagger & INFACE)
    open->stagger = OUTFACE|INFACET;

  /*-----------------------------------------------------------------*/
  /* Boundary ghost zone */
  sprintf(keyword, "BOUNDARY%1d.GHOST_CELLS", n);
  if (!(prm_read_int(fp, keyword, &open->bgz))) {
    for (i = 0; i < open->ntr; i++) {
      if (open->bcond_tra[i] & (TRCONC|TRCONF))
	open->bgz = laux;
    }
  }
}

/* END get_OBC_conds()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Initialise the open boundary structure                            */
/*-------------------------------------------------------------------*/
void init_OBC_conds(parameters_t *params, open_bdrys_t *open)
{
  int i;

  open->bathycon = 0;
  open->smooth_z = 0;
  open->smooth_n = 0;
  open->bout = -1;
  open->relax_time = 0;
  open->relax_timei = 0;
  open->sponge_zone = 0;
  open->sponge_zone_h = 0.0;
  open->sponge_f = 0.0;
  open->relax_zone_nor = 0;
  open->rnor_b = 0.0;
  open->rnor_i = 0.0;
  open->relax_zone_tan = 0;
  open->rtan_b = 0.0;
  open->rtan_i = 0.0;
  open->inverse_barometer = 1;
  open->meanc = 0.0;
  open->bflux_2d = 0.0;
  open->bflux_3d = 0.0;
  open->adjust_flux = 0.0;
  open->adjust_flux_s = 0.0;
  open->spf = 0.0;
  open->relax_ele = 0.0;
  open->file_dt = 0.0;
  open->file_next = params->t;
  open->rele_b = open->rele_i = 0.0; 
  open->stagger = OUTFACE;
  sprintf(open->bflow, "%c", '\0');
  open->rlen = 0.0;
  open->bhc = NOTVALID;
  open->bgz = 0;
  open->v1 = open->v2 = open->v3 = 0.0;
  open->sbcond = 0;
  open->options = NONE;
  open->minlat = NOTVALID;
  open->minlon = NOTVALID;
  open->maxlat = NOTVALID;
  open->maxlon = NOTVALID;
  open->elon = open->elat = NOTVALID;
  open->slon = open->slat = NOTVALID;
  open->mlon = open->mlat = NOTVALID;
  open->intype = 0;
  open->nedges = 0;
  if (open->ntr > 0) {
    /*UR-FIX assign extra element  bcond_tra - otherwise read 
     * from uninitialised value 
     * org
     * open->bcond_tra = i_alloc_1d(open->ntr);
     */
    open->bcond_tra = i_alloc_1d(open->ntr+1);
    open->clampv = d_alloc_1d(open->ntr);
    open->relax_zone_tra = i_alloc_1d(open->ntr);
    open->rtra_b = d_alloc_1d(open->ntr);
    open->rtra_i = d_alloc_1d(open->ntr);
    open->trpc = d_alloc_1d(open->ntr);
    /*UR-FIX assign last element */
    open->bcond_tra[open->ntr] = -1;
  }
  memset(open->relax_zone_tra, 0, open->ntr * sizeof(int));
  memset(open->rtra_b, 0, open->ntr * sizeof(double));
  memset(open->rtra_i, 0, open->ntr * sizeof(double));
  memset(open->clampv, 0, open->ntr * sizeof(int));
  memset(open->trpc, 0, open->ntr * sizeof(double));
  open->bstdf = -1;
  open->nbstd = 0;
  open->bstd = (char **)malloc(NBSTD * sizeof(char *));
  for (i = 0; i < NBSTD; i++) {
    open->bstd[i] = (char *)malloc(sizeof(char)*MAXSTRLEN);
    sprintf(open->bstd[i], "%c", '\0');
  }
}

/* END init_OBC_conds()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read the open boundary condition types and save to the */
/* open boundary data structure.                                     */
/*-------------------------------------------------------------------*/
void copy_OBC_conds(open_bdrys_t *io, /* ParamStruct open boundary data
                                         structure */
                    open_bdrys_t *open, /* master_t open boundary data
                                           structure */
                    int n,      /* Open boundary number */
                    tracer_info_t *tracer /* Tracer info data structure */
  )
{
  int i;                        /* Counters */

  /* Get the open boundary condition type */
  open->type = io->type;
  open->id = io->id;

  /* Get the name of the open boundary */
  strcpy(open->name, io->name);

  /* Allocate tracer memory */
  open->ntr = io->ntr;
  open->atr = io->atr;
  if (open->ntr > 0) {
    open->bcond_tra = i_alloc_1d(open->ntr+1);
    open->clampv = d_alloc_1d(open->ntr);
    open->relax_zone_tra = i_alloc_1d(open->ntr);
    open->rtra_b = d_alloc_1d(open->ntr);
    open->rtra_i = d_alloc_1d(open->ntr);
    open->trpc = d_alloc_1d(open->ntr);
    open->bcond_tra[open->ntr] = -1;
  }

  /* Get the open boundary conditions for each component */
  open->bcond_nor = io->bcond_nor;
  open->bcond_nor2d = io->bcond_nor2d;
  open->relax_zone_nor = io->relax_zone_nor;
  open->rnor_b = io->rnor_b;
  open->rnor_i = io->rnor_i;
  open->linear_zone_nor = io->linear_zone_nor;
  open->bcond_tan = io->bcond_tan;
  open->bcond_tan2d = io->bcond_tan2d;
  open->relax_zone_tan = io->relax_zone_tan;
  open->rtan_b = io->rtan_b;
  open->rtan_i = io->rtan_i;
  open->linear_zone_tan = io->linear_zone_tan;
  open->bcond_ele = io->bcond_ele;
  open->bcond_w = io->bcond_w;
  open->bcond_Vz = io->bcond_Vz;
  open->bcond_Kz = io->bcond_Kz;
  open->relax_zone_ele = io->relax_zone_ele;
  open->relax_ele = io->relax_ele;
  open->rele_b = io->rele_b;
  open->rele_i = io->rele_i;
  open->adjust_flux = io->adjust_flux;
  open->adjust_flux_s = io->adjust_flux_s;
  open->spf = io->spf;
  open->stagger = io->stagger;
  open->meanc = io->meanc;
  open->bflux_2d = io->bflux_2d;
  open->bflux_3d = io->bflux_3d;
  open->intype = io->intype;
  for (i = 0; i < open->ntr; i++) {
    open->relax_zone_tra[i] = io->relax_zone_tra[i];
    open->rtra_b[i] = io->rtra_b[i];
    open->rtra_i[i] = io->rtra_i[i];
    open->bcond_tra[i] = io->bcond_tra[i];
    open->clampv[i] = io->clampv[i];
    open->trpc[i] = io->trpc[i];
  }
  open->relax_time = io->relax_time;
  open->relax_timei = io->relax_timei;
  open->tidemem_depth = io->tidemem_depth;
  open->ncyc = io->ncyc;
  if (open->ncyc > 0) {
    open->ilocc = i_alloc_1d(open->ncyc);
    memcpy(open->ilocc, io->ilocc, open->ncyc * sizeof(int));
    open->jlocc = i_alloc_1d(open->ncyc);
    memcpy(open->jlocc, io->jlocc, open->ncyc * sizeof(int));
  }

  /*-----------------------------------------------------------------*/
  /* Read general boundary information */
  open->sponge_zone = io->sponge_zone;
  open->sponge_zone_h = io->sponge_zone_h;
  open->sponge_f = io->sponge_f;
  open->inverse_barometer = io->inverse_barometer;
  open->upmeth = io->upmeth;
  open->bgz = io->bgz;
  open->rlen = io->rlen;
  open->options = io->options;
  strcpy(open->tsfn, io->tsfn);
  open->sbcond = io->sbcond;
  open->bstdf = io->bstdf;
  open->nbstd = io->nbstd;
  if (open->nbstd) {
    open->bstd = (char **)malloc(open->nbstd * sizeof(char *));
    for (i = 0; i < open->nbstd; i++) {
      open->bstd[i] = (char *)malloc(sizeof(char)*MAXSTRLEN);
      strcpy(open->bstd[i], io->bstd[i]);
    }
  }
  open->maxlat = io->maxlat;
  open->minlat = io->minlat;
  open->maxlon = io->maxlon;
  open->minlon = io->minlon;
  open->nedges = io->nedges;
  /*
  if (io->nedges) {
    io->edges = i_alloc_1d(open->nedges);
    memcpy(io->edges, open->edges, io->nedges * sizeof(int));
  }
  */
  /* Tidal forcing for elevation */
  strcpy(open->tide_con, io->tide_con);
  open->ntide = io->ntide;
  if (open->ntide > 0)
    open->tideforce =
      (tide_details_t *)malloc(sizeof(tide_details_t) * open->ntide);
  for (i = 0; i < open->ntide; i++) {
    open->tideforce[i].ntides = io->tideforce[i].ntides;
    strcpy(open->tideforce[i].tname, io->tideforce[i].tname);
    open->tideforce[i].ic = io->tideforce[i].ic;
    open->tideforce[i].jc = io->tideforce[i].jc;
    open->tideforce[i].amp = io->tideforce[i].amp;
    open->tideforce[i].per = io->tideforce[i].per;
    open->tideforce[i].mta = io->tideforce[i].mta;
    open->tideforce[i].dta = io->tideforce[i].dta;
    open->tideforce[i].mtp = io->tideforce[i].mtp;
    open->tideforce[i].dtp = io->tideforce[i].dtp;
  }
}

/* END copy_OBC_conds()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read in the eta boundary relaxation                    */
/*-------------------------------------------------------------------*/
void get_OBC_relax(parameters_t *params, open_bdrys_t *open, FILE *fp, int n)
{
  char buf[MAXSTRLEN];

  if (open->relax_ele || (open->adjust_flux && params->compatible & V1670)) {
    if(strlen(params->etarlxn)) {
      hd_warn("Eta relaxation file for boundary %d = %s\n", n,
	      params->etarlxn);
      params->etarlx |= BOUNDARY;
    } else {
      if (prm_read_char(fp, "eta_relaxation_file", params->etarlxn) > 0) {
	if (prm_read_char(fp, "eta_relaxation_input_dt", buf) > 0)
	  tm_scale_to_secs(buf, &params->etarlxdt);
	else
	  params->etarlxdt = 86400;
	params->etarlx = BOUNDARY;
	params->ntrS++;
      } else {
	hd_warn("Must supply an eta_relaxation_file for boundary %d (%s)\n", n, open->name);
	open->relax_ele = 0;
	open->rele_b = open->rele_i = 0.0; 
	open->adjust_flux = 0.0;
      }
    }
  }
}

/* END get_OBC_relax()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Returns 1 if a passive OBC is used in conjunction with an active  */
/* OBC and a relaxation time constant is required.                   */
/*-------------------------------------------------------------------*/
int sforce_relax(int bcond      /* OBC code */
  )
{
  if (bcond == (NOGRAD | FILEIN))
    return (1);
  else if (bcond == (ORLANS | FILEIN))
    return (1);
  else if (bcond == (GRAVTY | FILEIN))
    return (1);
  else if (bcond == (CAMOBR | FILEIN))
    return (1);
  else if (bcond == (RAYMND | FILEIN))
    return (1);
  else if (bcond == (MILLER | FILEIN))
    return (1);
  else if (bcond == (CLAMPD | FILEIN))
    return (1);
  else if (bcond == (LINEXT | FILEIN))
    return (1);
  else if (bcond == (POLEXT | FILEIN))
    return (1);
  else if (bcond == (NOTHIN | FILEIN))
    return (1);
  else if (bcond == (NOTHIN | FILEIN | RAYMND))
    return (1);
  else if (bcond == (NOTHIN | FILEIN | GRAVTY))
    return (1);
  else if (bcond == (NOTHIN | FILEIN | ORLANS))
    return (1);
  else if (bcond == (NOTHIN | FILEIN | MILLER))
    return (1);
  else if (bcond == (NOTHIN | FILEIN | CAMOBR))
    return (1);
  else if (bcond == (NOGRAD | CUSTOM))
    return (1);
  else if (bcond == (ORLANS | CUSTOM))
    return (1);
  else if (bcond == (GRAVTY | CUSTOM))
    return (1);
  else if (bcond == (CAMOBR | CUSTOM))
    return (1);
  else if (bcond == (RAYMND | CUSTOM))
    return (1);
  else if (bcond == (MILLER | CUSTOM))
    return (1);
  else if (bcond == (CLAMPD | CUSTOM))
    return (1);
  else if (bcond == (LINEXT | CUSTOM))
    return (1);
  else if (bcond == (POLEXT | CUSTOM))
    return (1);
  else if (bcond == (NOTHIN | CUSTOM))
    return (1);
  else if (bcond == (NOGRAD | LINEAR))
    return (1);
  else if (bcond == (ORLANS | LINEAR))
    return (1);
  else if (bcond == (GRAVTY | LINEAR))
    return (1);
  else if (bcond == (CAMOBR | LINEAR))
    return (1);
  else if (bcond == (RAYMND | LINEAR))
    return (1);
  else if (bcond == (MILLER | LINEAR))
    return (1);
  else if (bcond == (CLAMPD | LINEAR))
    return (1);
  else if (bcond == (LINEXT | LINEAR))
    return (1);
  else if (bcond == (POLEXT | LINEAR))
    return (1);
  else if (bcond == (NOGRAD | NOTHIN))
    return (1);
  else if (bcond == (ORLANS | NOTHIN))
    return (1);
  else if (bcond == (GRAVTY | NOTHIN))
    return (1);
  else if (bcond == (CAMOBR | NOTHIN))
    return (1);
  else if (bcond == (RAYMND | NOTHIN))
    return (1);
  else if (bcond == (MILLER | NOTHIN))
    return (1);
  else if (bcond == (CLAMPD | NOTHIN))
    return (1);
  else if (bcond == (LINEXT | NOTHIN))
    return (1);
  else if (bcond == (POLEXT | NOTHIN))
    return (1);
  else if (bcond == (NOGRAD | TIDEBC))
    return (1);
  else if (bcond == (ORLANS | TIDEBC))
    return (1);
  else if (bcond == (GRAVTY | TIDEBC))
    return (1);
  else if (bcond == (CAMOBR | TIDEBC))
    return (1);
  else if (bcond == (RAYMND | TIDEBC))
    return (1);
  else if (bcond == (MILLER | TIDEBC))
    return (1);
  else if (bcond == (CLAMPD | TIDEBC))
    return (1);
  else if (bcond == (LINEXT | TIDEBC))
    return (1);
  else if (bcond == (POLEXT | TIDEBC))
    return (1);
  else if (bcond == (NOGRAD | TIDALH))
    return (1);
  else if (bcond == (ORLANS | TIDALH))
    return (1);
  else if (bcond == (GRAVTY | TIDALH))
    return (1);
  else if (bcond == (CAMOBR | TIDALH))
    return (1);
  else if (bcond == (RAYMND | TIDALH))
    return (1);
  else if (bcond == (MILLER | TIDALH))
    return (1);
  else if (bcond == (CLAMPD | TIDALH))
    return (1);
  else if (bcond == (LINEXT | TIDALH))
    return (1);
  else if (bcond == (POLEXT | TIDALH))
    return (1);
  else if (bcond == (FILEIN | NOGRAD | TIDALH))
    return (1);
  else if (bcond == (FILEIN | ORLANS | TIDALH))
    return (1);
  else if (bcond == (FILEIN | GRAVTY | TIDALH))
    return (1);
  else if (bcond == (FILEIN | CAMOBR | TIDALH))
    return (1);
  else if (bcond == (FILEIN | RAYMND | TIDALH))
    return (1);
  else if (bcond == (FILEIN | MILLER | TIDALH))
    return (1);
  else if (bcond == (FILEIN | CLAMPD | TIDALH))
    return (1);
  else if (bcond == (FILEIN | LINEXT | TIDALH))
    return (1);
  else if (bcond == (FILEIN | POLEXT | TIDALH))
    return (1);
  else if (bcond == (NOGRAD | TIDALC))
    return (1);
  else if (bcond == (ORLANS | TIDALC))
    return (1);
  else if (bcond == (GRAVTY | TIDALC))
    return (1);
  else if (bcond == (CAMOBR | TIDALC))
    return (1);
  else if (bcond == (RAYMND | TIDALC))
    return (1);
  else if (bcond == (MILLER | TIDALC))
    return (1);
  else if (bcond == (CLAMPD | TIDALC))
    return (1);
  else if (bcond == (LINEXT | TIDALC))
    return (1);
  else if (bcond == (POLEXT | TIDALC))
    return (1);
  else if (bcond == (FILEIN | NOGRAD | TIDALC))
    return (1);
  else if (bcond == (FILEIN | ORLANS | TIDALC))
    return (1);
  else if (bcond == (FILEIN | GRAVTY | TIDALC))
    return (1);
  else if (bcond == (FILEIN | CAMOBR | TIDALC))
    return (1);
  else if (bcond == (FILEIN | RAYMND | TIDALC))
    return (1);
  else if (bcond == (FILEIN | MILLER | TIDALC))
    return (1);
  else if (bcond == (FILEIN | CLAMPD | TIDALC))
    return (1);
  else if (bcond == (FILEIN | LINEXT | TIDALC))
    return (1);
  else if (bcond == (FILEIN | POLEXT | TIDALC))
    return (1);

  else
    return (0);
}

/* END sforce_relax()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Returns the boundary condition code number for a given boundary   */
/* condition text string.                                            */
/*-------------------------------------------------------------------*/
int get_bcond_no(char *list)
{
  int code = 0;
  char *tok;
  tok = strtok(list, "|");
  code |= bcond_no(tok);
  while (tok != NULL) {
    tok = strtok(NULL, "|");
    code |= bcond_no(tok);
  }
  if (code == 0)
    hd_warn("Cannot find valid boundary condition code %d for %s\n", code, list);
  return (code);
}

/* END get_bcond_no()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Compares text with boundary condition names                       */
/*-------------------------------------------------------------------*/
int bcond_no(char *list)
{
  int code = 0;
  if (list == NULL)
    return (code);
  if (strcmp(list, "FILEIN") == 0)
    code = FILEIN;
  if (strcmp(list, "ORLANS") == 0)
    code = ORLANS;
  if (strcmp(list, "NOGRAD") == 0)
    code = NOGRAD;
  if (strcmp(list, "CLAMPD") == 0)
    code = CLAMPD;
  if (strcmp(list, "LINEXT") == 0)
    code = LINEXT;
  if (strcmp(list, "POLEXT") == 0)
    code = POLEXT;
  if (strcmp(list, "GRAVTY") == 0)
    code = GRAVTY;
  if (strcmp(list, "UPSTRM") == 0)
    code = UPSTRM;
  if (strcmp(list, "TRFLUX") == 0)
    code = TRFLUX;
  if (strcmp(list, "TRCONC") == 0)
    code = TRCONC;
  if (strcmp(list, "TRCONF") == 0)
    code = TRCONF;
  if (strcmp(list, "PROFIL") == 0)
    code = PROFIL;
  if (strcmp(list, "DEPROF") == 0)
    code = DEPROF;
  if (strcmp(list, "DESCAL") == 0)
    code = DESCAL;
  if (strcmp(list, "TIDALM") == 0)
    code = TIDALM;
  if (strcmp(list, "TIDALH") == 0)
    code = TIDALH;
  if (strcmp(list, "TIDALC") == 0)
    code = TIDALC;
  if (strcmp(list, "TIDEBC") == 0)
    code = TIDEBC;
  if (strcmp(list, "NOTHIN") == 0)
    code = NOTHIN;
  if (strcmp(list, "CYCLIC") == 0)
    code = CYCLIC;
  if (strcmp(list, "CYCLED") == 0)
    code = CYCLED;
  if (strcmp(list, "CAMOBR") == 0)
    code = CAMOBR;
  if (strcmp(list, "RAYMND") == 0)
    code = RAYMND;
  if (strcmp(list, "MILLER") == 0)
    code = MILLER;
  if (strcmp(list, "VERTIN") == 0)
    code = VERTIN;
  if (strcmp(list, "CUSTOM") == 0)
    code = CUSTOM;
  if (strcmp(list, "LINEAR") == 0)
    code = LINEAR;
  if (strcmp(list, "STATIS") == 0)
    code = STATIS;
  if (strcmp(list, "FLATHR") == 0)
    code = FLATHR;
  if (strcmp(list, "FLATHE") == 0)
    code = FLATHE;
  if (strcmp(list, "LOCALN") == 0)
    code = LOCALN;
  if (strcmp(list, "LOCALT") == 0)
    code = LOCALT;
  if (strcmp(list, "LOCALE") == 0)
    code = LOCALE;
  return (code);
}

/* END bcond_no()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
void bcname(int code, char *name)
{
  /* 
     name=malloc(MAXSTRLEN); strcpy(name,"NULL"); */
  if (code & STATIS)
    strcpy(name, "STATIS");
  if (code & LINEAR)
    strcpy(name, "LINEAR");
  if (code & CUSTOM)
    strcpy(name, "CUSTOM");
  if (code & VERTIN)
    strcpy(name, "VERTIN");
  if (code & MILLER)
    strcpy(name, "MILLER");
  if (code & RAYMND)
    strcpy(name, "RAYMND");
  if (code & CAMOBR)
    strcpy(name, "CAMOBR");
  if (code & CYCLIC)
    strcpy(name, "CYCLIC");
  if (code & CYCLED)
    strcpy(name, "CYCLED");
  if (code & TIDEBC)
    strcpy(name, "TIDEBC");
  if (code & NOTHIN)
    strcpy(name, "NOTHIN");
  if (code & TIDALH)
    strcpy(name, "TIDALH");
  if (code & TIDALC)
    strcpy(name, "TIDALC");
  /*
  if (code & TIDALM)
    strcpy(name, "TIDALM");
  */
  if (code & UPSTRM)
    strcpy(name, "UPSTRM");
  if (code & TRFLUX)
    strcpy(name, "TRFLUX");
  if (code & TRCONC)
    strcpy(name, "TRCONC");
  if (code & TRCONF)
    strcpy(name, "TRCONF");
  if (code & PROFIL)
    strcpy(name, "PROFIL");
  if (code & DEPROF)
    strcpy(name, "DEPROF");
  if (code & DESCAL)
    strcpy(name, "DESCAL");
  if (code & GRAVTY)
    strcpy(name, "GRAVTY");
  if (code & POLEXT)
    strcpy(name, "POLEXT");
  if (code & LINEXT)
    strcpy(name, "LINEXT");
  if (code & NOGRAD)
    strcpy(name, "NOGRAD");
  if (code & CLAMPD)
    strcpy(name, "CLAMPD");
  if (code & ORLANS)
    strcpy(name, "ORLANS");
  if (code & FILEIN)
    strcpy(name, "FILEIN");
  if (code & LOCALN)
    strcpy(name, "LOCALN");
  if (code & LOCALT)
    strcpy(name, "LOCALT");
  if (code & LOCALE)
    strcpy(name, "LOCALE");
  if (code & FLATHR)
    strcpy(name, "FLATHR");
  if (code & FLATHE)
    strcpy(name, "FLATHR");

  if (code & FILEIN && strcmp(name, "FILEIN") != 0)
    sprintf(name, "%s|FILEIN", name);
  if (code & DESCAL && strcmp(name, "DESCAL") != 0)
    sprintf(name, "%s|DESCAL", name);
  if (code & ORLANS && strcmp(name, "ORLANS") != 0)
    sprintf(name, "%s|ORLANS", name);
  if (code & CLAMPD && strcmp(name, "CLAMPD") != 0)
    sprintf(name, "%s|CLAMPD", name);
  if (code & NOGRAD && strcmp(name, "NOGRAD") != 0)
    sprintf(name, "%s|NOGRAD", name);
  if (code & LINEXT && strcmp(name, "LINEXT") != 0)
    sprintf(name, "%s|LINEXT", name);
  if (code & POLEXT && strcmp(name, "POLEXT") != 0)
    sprintf(name, "%s|POLEXT", name);
  if (code & GRAVTY && strcmp(name, "GRAVTY") != 0)
    sprintf(name, "%s|GRAVTY", name);
  if (code & UPSTRM && strcmp(name, "UPSTRM") != 0)
    sprintf(name, "%s|UPSTRM", name);
  if (code & TRFLUX && strcmp(name, "TRFLUX") != 0)
    sprintf(name, "%s|TRFLUX", name);
  if (code & TRCONC && strcmp(name, "TRCONC") != 0)
    sprintf(name, "%s|TRCONC", name);
  if (code & TRCONF && strcmp(name, "TRCONF") != 0)
    sprintf(name, "%s|TRCONF", name);
  if (code & PROFIL && strcmp(name, "PROFIL") != 0)
    sprintf(name, "%s|PROFIL", name);
  if (code & DEPROF && strcmp(name, "DEPROF") != 0)
    sprintf(name, "%s|DEPROF", name);
  /*
  if (code & TIDALM && strcmp(name, "TIDALM") != 0)
    sprintf(name, "%s|TIDALM", name);
  */
  if (code & TIDALH && strcmp(name, "TIDALH") != 0)
    sprintf(name, "%s|TIDALH", name);
  if (code & TIDALC && strcmp(name, "TIDALC") != 0)
    sprintf(name, "%s|TIDALC", name);
  if (code & NOTHIN && strcmp(name, "NOTHIN") != 0)
    sprintf(name, "%s|NOTHIN", name);
  if (code & TIDEBC && strcmp(name, "TIDEBC") != 0)
    sprintf(name, "%s|TIDEBC", name);
  if (code & CYCLIC && strcmp(name, "CYCLIC") != 0)
    sprintf(name, "%s|CYCLIC", name);
  if (code & CYCLED && strcmp(name, "CYCLED") != 0)
    sprintf(name, "%s|CYCLED", name);
  if (code & CAMOBR && strcmp(name, "CAMOBR") != 0)
    sprintf(name, "%s|CAMOBR", name);
  if (code & RAYMND && strcmp(name, "RAYMND") != 0)
    sprintf(name, "%s|RAYMND", name);
  if (code & MILLER && strcmp(name, "MILLER") != 0)
    sprintf(name, "%s|MILLER", name);
  if (code & VERTIN && strcmp(name, "VERTIN") != 0)
    sprintf(name, "%s|VERTIN", name);
  if (code & CUSTOM && strcmp(name, "CUSTOM") != 0)
    sprintf(name, "%s|CUSTOM", name);
  if (code & LINEAR && strcmp(name, "LINEAR") != 0)
    sprintf(name, "%s|LINEAR", name);
  if (code & STATIS && strcmp(name, "STATIS") != 0)
    sprintf(name, "%s|STATIS", name);
  if (code & LOCALN && strcmp(name, "LOCALN") != 0)
    sprintf(name, "%s|LOCALN", name);
  if (code & LOCALT && strcmp(name, "LOCALT") != 0)
    sprintf(name, "%s|LOCALT", name);
  if (code & LOCALE && strcmp(name, "LOCALE") != 0)
    sprintf(name, "%s|LOCALE", name);
  /* 
     free(name); return(name); */
}

/* END bcname()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads the relaxation time scales                                  */
/*-------------------------------------------------------------------*/
void get_relax_time(open_bdrys_t *open,  /* Open boundary strcuture  */
		    FILE *fp,            /* Parameter file handle    */
		    char *bname,         /* OBC name                 */
		    int n,               /* OBC number               */
		    int obc1,            /* OBC type                 */
		    int obc2             /* OBC type                 */
		    )
{
  char key1[MAXSTRLEN];
  char key2[MAXSTRLEN];

  if (open->relax_time <= 0 && (sforce_relax(obc1) ||
				sforce_relax(obc2))) {
    sprintf(key1, "BOUNDARY%1d.RELAX_TIME", n);
    sprintf(key2, "BOUNDARY%1d.RELAX_OUT", n);
    if (!(prm_get_time_in_secs(fp, key1, &open->relax_time)) &&
	!(prm_get_time_in_secs(fp, key2, &open->relax_time)))
      hd_quit("Boundary %d condition %s requires RELAX_TIME\n", n, bname);
    sprintf(key1, "BOUNDARY%1d.RELAX_IN", n);
    if (!(prm_get_time_in_secs(fp, key1, &open->relax_timei)))
      open->relax_timei = open->relax_time;
  }
}
/* END get_relax_time()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read in tidal constituent data to the tide_details_t   */
/* structure for each edge.                                          */
/*-------------------------------------------------------------------*/
tide_details_t *sreadtideforce(FILE * fp, int bnum, int n)
{
  int i;
  tide_details_t *tide = NULL;
  char keyword[MAXSTRLEN];      /* Buffer for I/O */

  if (n > 0) {
    tide = (tide_details_t *)malloc(sizeof(tide_details_t) * n);
    memset(tide, 0, sizeof(tide_details_t) * n);
    for (i = 0; i < n; i++) {
      tide[i].ntides = n;
      sprintf(keyword, "BOUNDARY%1d.T_NAME", bnum);
      prm_read_char(fp, keyword, tide[i].tname);
      sprintf(keyword, "BOUNDARY%1d.T_XLOCATION", bnum);
      prm_read_double(fp, keyword, &tide[i].ic);
      sprintf(keyword, "BOUNDARY%1d.T_YLOCATION", bnum);
      prm_read_double(fp, keyword, &tide[i].jc);
      sprintf(keyword, "BOUNDARY%1d.T_AMPLITUDE", bnum);
      prm_read_double(fp, keyword, &tide[i].amp);
      sprintf(keyword, "BOUNDARY%1d.T_PERIOD", bnum);
      prm_get_time_in_secs(fp, keyword, &tide[i].per);
      /*
      prm_read_double(fp, keyword, &tide[i].per);
      tide[i].per *= 3600.0;
      */
      sprintf(keyword, "BOUNDARY%1d.T_MOD_AMP", bnum);
      prm_read_double(fp, keyword, &tide[i].mta);
      tide[i].mta *= 1e-5;
      sprintf(keyword, "BOUNDARY%1d.T_DIR_AMP", bnum);
      prm_read_double(fp, keyword, &tide[i].dta);
      tide[i].dta = 270.0 - tide[i].dta;
      sprintf(keyword, "BOUNDARY%1d.T_MOD_PSE", bnum);
      prm_read_double(fp, keyword, &tide[i].mtp);
      tide[i].mtp *= (1e-3 * PI / 180);
      sprintf(keyword, "BOUNDARY%1d.T_DIR_PSE", bnum);
      prm_read_double(fp, keyword, &tide[i].dtp);
      tide[i].dtp = 270.0 - tide[i].dtp;
    }
  } else
    hd_quit
      ("readtideforce: No tidal forcing specified for boundary edge %d.\n",
       bnum);
  prm_set_errfn(hd_silent_warn);

  return tide;
}

/* END sreadtideforce()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set up the custom boundary routines on the master.     */
/* Note, this can only be done after the master data structure is    */
/* created and initialised.                                          */
/*-------------------------------------------------------------------*/
void bdry_custom_m(parameters_t *params,  /* Input parameter data    */
                   geometry_t *geom,      /* Global geometry         */
                   master_t *master       /* Master data             */
  )
{
  int n, t;
  open_bdrys_t *open, *io;
  tracer_info_t *tracer;
  char files[MAXNUMTSFILES][MAXSTRLEN];

  for (n = 0; n < geom->nobc; n++) {
    open = geom->open[n];
    io = params->open[n];
    tracer = params->trinfo_3d;

    /*---------------------------------------------------------------*/
    /* Check if a custom boundary condition has been specified       */
    strcpy(open->cusname_u1, io->cusname_u1);
    strcpy(open->cusname_u2, io->cusname_u2);
    strcpy(open->cusname_u1av, io->cusname_u1av);
    strcpy(open->cusname_u2av, io->cusname_u2av);
    strcpy(open->cusname_eta, io->cusname_eta);
    for (t = 0; t < open->ntr; ++t) {
      open->cusname_t[t] = (char *)malloc(MAXSTRLEN);
      strcpy(open->cusname_t[t], io->cusname_t[t]);
    }
    strcpy(open->bflow, io->bflow);
    if (io->bhc == NOTVALID) io->bhc = open->meandep;
    open->bhc = io->bhc;
    bdry_custom_init(master, open, tracer, io->custype);

    /*---------------------------------------------------------------*/
    /* Read the boundary scaling                                     */
    read_bdry_eta_scale(master, io, open);
    if (open->ntr > 0)
      read_bdry_scale(io, open, master->trname);
  }

  /*-----------------------------------------------------------------*/
  /* Initialise the custom boundary routine                          */
  bdry_init_m(master);
}

/* END bdry_custom_m()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set up custom boundaries                               */
/*-------------------------------------------------------------------*/
void bdry_custom_init(master_t *master, open_bdrys_t *open, 
		      tracer_info_t *tracer, int type)
{
  int t;
  char files[MAXNUMTSFILES][MAXSTRLEN];

  /*-----------------------------------------------------------------*/
  /* Copy boundary forcing data if present                           */
  if (strlen(open->tsfn)) {
    open->ntsfiles =
      parseline(open->tsfn, (char **)files, MAXNUMTSFILES);
    
    open->filenames =
      (cstring *) malloc(sizeof(cstring) * open->ntsfiles);
    for (t = 0; t < open->ntsfiles; ++t)
      strcpy(open->filenames[t], ((char **)files)[t]);
      open->tsfiles = hd_ts_multifile_read(master, open->ntsfiles,
                                           open->filenames);
  }

  /*-----------------------------------------------------------------*/
  /* Check if a custom boundary condition has been specified         */
  sprintf(open->datau1.name, "%c", '\0');
  sprintf(open->datau2.name, "%c", '\0');
  sprintf(open->datau1av.name, "%c", '\0');
  sprintf(open->datau2av.name, "%c", '\0');
  sprintf(open->etadata.name, "%c", '\0');
  if (open->bhc == NOTVALID) open->bhc = open->meandep;

  /* Check for u1av FILEIN or CUSTOM input                         */
  read_bdry_custom(open, open->id, "u1av", 0.0, &open->datau1av,
		   open->cusname_u1av);
  /* Check for u2av FILEIN or CUSTOM input                         */
  read_bdry_custom(open, open->id, "u2av", 0.0, &open->datau2av,
		   open->cusname_u2av);
  /* Check for eta FILEIN or CUSTOM input                          */
  read_bdry_custom(open, open->id, "eta", 0.0, &open->etadata,
		   open->cusname_eta);
  /* Check for u1 or u2 FILEIN or CUSTOM input.                    */
  read_bdry_custom(open, open->id, "u1", 0.0, &open->datau1, 
		   open->cusname_u1);
  read_bdry_custom(open, open->id, "u2", 0.0, &open->datau2, 
		   open->cusname_u2);

  /*---------------------------------------------------------------*/
  /* Read the tracer custom routines                               */
  if (open->ntr > 0) {
    open->bdata_t = (bdry_details_t *)malloc(sizeof(bdry_details_t) *
					     open->ntr);
    memset(open->bdata_t, 0, sizeof(bdry_details_t) * open->ntr);
    for (t = 0; t < open->ntr; ++t) {
      read_bdry_custom(open, open->id, tracer[t].name, tracer[t].fill_value_wc,
		       &open->bdata_t[t], open->cusname_t[t]);
    }
  }
}

/* END bdry_custom_init()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to allocate custom boudary functions                      */
/*-------------------------------------------------------------------*/
int read_bdry_custom(open_bdrys_t *open,  /* Open boundary structure */
                     int bnum,            /* Boundary number         */
                     char *name,  /* Name of boundary type           */
                     double fill, /* Fill value for boundary         */
                     bdry_details_t *data,  /* Custom boundary data  */
                     char *cusname  /* Name of the custom boundary
                                       function */
  )
{
  int i;
  int ret = 0;
  prm_set_errfn(hd_silent_warn);
  data->explct = 0;
  data->nargs = 0;
  data->args = NULL;

  /*-----------------------------------------------------------------*/
  /* Set the edge name according to the specified forcing            */
  if (open->type & U1BDRY && strcmp(name, "u1") == 0 &&       
      open->bcond_nor & (FILEIN | CUSTOM)) {            /* u1 NOR    */
    strcpy(data->name, "u1");
    data->type = NOR;
  }
  else if (open->type & U1BDRY && strcmp(name, "u2") == 0 && 
	   open->bcond_tan & (FILEIN | CUSTOM)) {       /* u2 TAN    */
    strcpy(data->name, "u2");
    data->type = TAN;
  }
  else if (open->type & U1BDRY && strcmp(name, "u1av") == 0 && 
	   open->bcond_nor2d & (FILEIN | CUSTOM)) {     /* u1av NOR  */
    strcpy(data->name, "u1av");
    data->type = NOR;
  }
  else if (open->type & U1BDRY && strcmp(name, "u2av") == 0 && 
	   open->bcond_tan2d & (FILEIN | CUSTOM)) {     /* u2av TAN  */
    strcpy(data->name, "u2av");
    data->type = TAN;
  }
  else if (strcmp(name, "eta") == 0 && open->bcond_ele & (FILEIN | CUSTOM))
    strcpy(data->name, "eta");
  else if (strcmp(name, "eta") != 0 && strcmp(name, "u1") != 0 &&
           strcmp(name, "u2") != 0)
    strcpy(data->name, name);

  if (strlen(cusname)) {
    char *fields[MAXSTRLEN * MAXNUMARGS];
    int nf = parseline(cusname, fields, MAXNUMARGS);
    data->explct = 1;

    /*---------------------------------------------------------------*/
    /* Check if it is a function                                     */
    if (nf <= 0)
      hd_quit
        ("readbdrycustom : '%s' custom boundary incorrectly specified.\n",
         name);

    locate_boundary_function(open, fields[0], &data->init_m, &data->init_w,
                             &data->custom_m, &data->custom_w,
                             &data->trans, &data->free_m, &data->free_w);

    if (data->custom_m == NULL) {
      if (nf == 1) {
        if (strcasecmp(fields[0], "default") == 0)
          data->fill_value = fill;
        else if (sscanf(fields[0], "%lf", &data->fill_value) != 1)
          hd_quit
            ("readbdrycustom : '%s' custom boundary incorrectly specified.\n",
             name);
      } else
        hd_quit
          ("readbdrycustom : '%s' custom boundary incorrectly specified.\n",
           name);
    }


    /*---------------------------------------------------------------*/
    /* Copy the custom information over.                             */
    else {
      strcpy(data->custom_tag, fields[0]);
      data->nargs = nf - 1;
      if (data->nargs) {
        data->args = (cstring *) malloc(sizeof(cstring) * data->nargs);
        memset(data->args, 0, sizeof(cstring) * data->nargs);
      }
      for (i = 0; i < data->nargs; ++i) {
        strcpy(data->args[i], fields[i + 1]);
      }
    }
    ret = 1;
  }
  prm_set_errfn(hd_silent_warn);

  return (ret);
}

/* END read_bdry_custom()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to initialise the custom boundary routines in each        */
/* window.                                                           */
/*-------------------------------------------------------------------*/
void bdry_custom_w(geometry_t *geom,  /* Global sparse geometry */
                   geometry_t **window  /* Window data structure */
  )
{
  open_bdrys_t **open;          /* Global open boundary structure */
  int n, nn, tn;                /* Counters */
  int c1;                       /* Sparse coordinates / counters */
  int nwindows = geom->nwindows;

  open = geom->open;
  /*-----------------------------------------------------------------*/
  /* Fill the window OBC structues with custom information */
  for (nn = 1; nn <= nwindows; nn++) {
    for (n = 0; n < geom->nobc; n++) {
      if ((c1 = geom->owc[n][nn]) != -1) {

        /* File forcing info : direct pointing. These variables */
        /* should never be used by the windows since the master */
        /* reads all file information. Include for consistency.  */
        window[nn]->open[c1]->ntsfiles = open[n]->ntsfiles;
        window[nn]->open[c1]->tsfiles = open[n]->tsfiles;
        window[nn]->open[c1]->filenames = open[n]->filenames;

        /* u1 or u2 custom boundary routines */
        if (strlen(open[n]->datau1.name))
          bdry_data_copy(&window[nn]->open[c1]->datau1, &open[n]->datau1);
        if (strlen(open[n]->datau2.name))
          bdry_data_copy(&window[nn]->open[c1]->datau2, &open[n]->datau2);
        if (strlen(open[n]->datau1av.name))
          bdry_data_copy(&window[nn]->open[c1]->datau1av, &open[n]->datau1av);
        if (strlen(open[n]->datau2av.name))
          bdry_data_copy(&window[nn]->open[c1]->datau2av, &open[n]->datau2av);
        if (strlen(open[n]->etadata.name))
          bdry_data_copy(&window[nn]->open[c1]->etadata,
                         &open[n]->etadata);

        /* Tracer custom boundary routines */
        window[nn]->open[c1]->bdata_t = (bdry_details_t *)
          malloc(sizeof(bdry_details_t) * open[n]->ntr);
        for (tn = 0; tn < open[n]->ntr; tn++)
          bdry_data_copy(&window[nn]->open[c1]->bdata_t[tn],
                         &open[n]->bdata_t[tn]);

        /* Initialise the window custom boundary function */
        bdry_init_w(window[nn], window[nn]->open[c1], open[n]);

	/*-----------------------------------------------------------*/
	/* Copy the boundary scaling                                 */
	window[nn]->open[c1]->sdata_e = 
	  (scale_details_t *)malloc(sizeof(scale_details_t));
	scale_data_copy(window[nn]->open[c1]->sdata_e,
			open[n]->sdata_e);
	if (open[n]->ntr > 0) {
	  window[nn]->open[c1]->sdata_t = 
	    (scale_details_t *)malloc(sizeof(scale_details_t) *
				      open[n]->ntr);
	  for (tn = 0; tn < open[n]->ntr; tn++)
	    scale_data_copy(&window[nn]->open[c1]->sdata_t[tn],
			    &open[n]->sdata_t[tn]);
	}
      }
    }
  }
}

/* END bdry_custom_w()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to copy boundary data structures                          */
/*-------------------------------------------------------------------*/
void bdry_data_copy(bdry_details_t *data, bdry_details_t *din)
{
  int n;

  memset(data, 0, sizeof(bdry_details_t));
  data->explct = din->explct;
  strcpy(data->name, din->name);
  data->fill_value = din->fill_value;
  data->type = din->type;
  strcpy(data->custom_tag, din->custom_tag);
  data->nargs = din->nargs;
  if (data->nargs) {
    data->args = (cstring *) malloc(sizeof(cstring) * data->nargs);
    for (n = 0; n < data->nargs; ++n)
      strcpy(data->args[n], din->args[n]);
  }
  /* Set the custom function for the window the same as the master */
  data->custom_w = din->custom_w;
  data->init_w = din->init_w;
  data->trans = din->trans;
  data->free_w = din->free_w;
}

/* END bdry_data_copy()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to copy boundary scale structures                         */
/*-------------------------------------------------------------------*/
void scale_data_copy(scale_details_t *data, scale_details_t *din)
{
  memset(data, 0, sizeof(scale_details_t));
  strcpy(data->name, din->name);
  data->type = din->type;
  data->ntr = din->ntr;
  data->fact = din->fact;
  data->val = din->val;
  data->flag = din->flag;
}

/* END scale_data_copy()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Checks if specified custom routines are standard custom routines  */
/*-------------------------------------------------------------------*/
void locate_boundary_function(open_bdrys_t *open, const char *tag,
                              bdrycustom_init_m * init_m,
                              bdrycustom_init_w * init_w,
                              bdrycustom_m * func_m,
                              bdrycustom_w * func_w, bdrytransfer * trans,
			      bdryfree_m * free_m, bdryfree_w * free_w )
{

  int i = 0;
  static bdryfn_map_t bdryfn_standard[] = {
    {"hdstd_to_u1", bf_hdstd_to_u1_init_m, bf_hdstd_to_u1_init_w,
     bf_hdstd_to_u1_m, bf_hdstd_to_u1_w, bf_hdstd_to_u1_t,
     bf_ts_free, bf_void_free},
    {"hdstd_to_u2", bf_hdstd_to_u2_init_m, bf_hdstd_to_u2_init_w,
     bf_hdstd_to_u2_m, bf_hdstd_to_u2_w, bf_hdstd_to_u2_t,
     bf_ts_free, bf_void_free},

    {"hdstd_to_u1av", bf_hdstd_to_u1_init_m, bf_uv_to_uav_init_w,
     bf_hdstd_to_u1av_m, bf_uv_to_u1av_w, bf_uvav_to_u1av_t,
     bf_ts_free, bf_void_free},
    {"hdstd_to_u2av", bf_hdstd_to_u2_init_m, bf_uv_to_uav_init_w,
     bf_hdstd_to_u2av_m, bf_uv_to_u2av_w, bf_uvav_to_u2av_t,
     bf_ts_free, bf_void_free},

    {"uv_to_u1", bf_uv_to_u1_init_m, bf_uv_to_u1_init_w,
     bf_uv_to_u1_m, bf_uv_to_u1_w, bf_uv_to_u1_t,
     bf_ts_free, bf_void_free},
    {"uv_to_u2", bf_uv_to_u2_init_m, bf_uv_to_u2_init_w,
     bf_uv_to_u2_m, bf_uv_to_u2_w, bf_uv_to_u2_t,
     bf_ts_free, bf_void_free},

    {"uv_to_u1av", bf_uv_to_u1_init_m, bf_uv_to_uav_init_w,
     bf_uv_to_u1av_m, bf_uv_to_u1av_w, bf_uvav_to_u1av_t,
     bf_ts_free, bf_c2cc_free},
    {"uv_to_u2av", bf_uv_to_u2_init_m, bf_uv_to_uav_init_w,
     bf_uv_to_u2av_m, bf_uv_to_u2av_w, bf_uvav_to_u2av_t,
     bf_ts_free, bf_c2cc_free},

    {"uvav_to_u1av", bf_uvav_to_u1av_init_m, bf_uvav_to_u1av_init_w,
     bf_uvav_to_u1av_m, bf_uvav_to_u1av_w, bf_uvav_to_u1av_t},
    {"uvav_to_u2av", bf_uvav_to_u2av_init_m, bf_uvav_to_u2av_init_w,
     bf_uvav_to_u2av_m, bf_uvav_to_u2av_w, bf_uvav_to_u2av_t,
     bf_ts_free, bf_void_free},
    /*
    {"uv_adj_u2av", bf_uv_adj_u2_init_m, bf_uv_adj_uav_init_w,
     bf_uv_adj_u2av_m, bf_uv_adj_u2av_w, bf_uv_adj_u2_t},
    */
    {"u1flowbdry", bf_u1_flow_init_m, bf_u1_flow_init_w,
     bf_flow_m, bf_u1_flow_w, bf_flow_t,
     bf_flow_free, bf_void_free},
    {"use_eqn", bf_use_eqn_init_m, NULL, bf_use_eqn_m, bf_use_eqn_w, 
     bf_use_eqn_trans, bf_eqn_free, bf_void_free},
    {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL}
  };

  while (bdryfn_standard[i].name != NULL) {
    if (strcmp(bdryfn_standard[i].name, tag) == 0) {
      *init_m = bdryfn_standard[i].init_m;
      *init_w = bdryfn_standard[i].init_w;
      *func_m = bdryfn_standard[i].func_m;
      *func_w = bdryfn_standard[i].func_w;
      *trans = bdryfn_standard[i].trans;
      *free_m = bdryfn_standard[i].free_m;
      *free_w = bdryfn_standard[i].free_w;
      return;
    }
    ++i;
  }
}

/* END locate_boudary_function()                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Frees boundary data on the master                                 */
/*-------------------------------------------------------------------*/
void bdry_custom_free(geometry_t *geom,    /* Global geometry        */
		      master_t *master,    /* Master data            */		 
		      geometry_t **window, /* Window geometry        */
		      open_bdrys_t *open  /* Open boundary structure */
		      )
{
  int n, bn, t;

  tracer_info_t *tracer;
  char files[MAXNUMTSFILES][MAXSTRLEN];

  /* Free ts files */
  if (open->ntsfiles) {
    for (t = 0; t < open->ntsfiles; ++t) {
      hd_ts_free(master, open->tsfiles[t]);
    }
    free(open->filenames);
  }    

  /* Reset the custom specification */
  sprintf(open->cusname_u1, "%c", '\0');
  sprintf(open->cusname_u2, "%c", '\0');
  sprintf(open->cusname_u1av, "%c", '\0');
  sprintf(open->cusname_u2av, "%c", '\0');
  sprintf(open->cusname_eta, "%c", '\0');
  for (t = 0; t < open->ntr; t++) {
    if (open->cusname_t[t])
      free((char *)open->cusname_t[t]);
    if (open->bdata_t)
      custom_free_m(master, &open->bdata_t[t]);
  }
  if (open->bdata_t)
    free((bdry_details_t *)open->bdata_t);
  custom_free_m(master, &open->datau1);
  custom_free_m(master, &open->datau2);
  custom_free_m(master, &open->datau1av);
  custom_free_m(master, &open->datau2av);
  custom_free_m(master, &open->etadata);

  /* Find the corresponding window boundaries */
  for (n = 1; n <= geom->nwindows; n++) {
    for (bn = 0; bn < window[n]->nobc; bn++) {
      open_bdrys_t *open_w = window[n]->open[bn];
      if (open->id == open_w->id) {
	for (t = 0; t < open_w->ntr; t++) {
	  if (open_w->bdata_t)
	    custom_free_w(window[n], &open_w->bdata_t[t]);
	}
	if (open_w->bdata_t)
	  free((bdry_details_t *)open_w->bdata_t);	
	custom_free_w(window[n], &open_w->datau1);
	custom_free_w(window[n], &open_w->datau2);
	custom_free_w(window[n], &open_w->datau1av);
	custom_free_w(window[n], &open_w->datau2av);
	custom_free_w(window[n], &open_w->etadata);
      }
    }
  }
}

void custom_free_m(master_t *master, bdry_details_t *data)
{
  int i;

  sprintf(data->name, "%c", '\0');
  if (data->free_m)
    data->free_m(master, data);
  data->free_m = NULL;
  data->init_m = NULL;
  data->custom_m = NULL;
  data->explct = 0;
  if (data->nargs)
    free(data->args);
  data->nargs = 0;
  data->args = NULL;
}

void custom_free_w(geometry_t *window, bdry_details_t *data)
{
  int i;

  sprintf(data->name, "%c", '\0');
  if (data->free_w)
    data->free_w(window, data);
  data->free_w = NULL;
  data->init_w = NULL;
  data->custom_w = NULL;
  data->trans = NULL;
  data->explct = 0;
  if (data->nargs)
    free(data->args);
  data->nargs = 0;
  data->args = NULL;
}

/* END bdry_custom_free_m()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to reconfigure the master if the model fails              */
/*-------------------------------------------------------------------*/
void bdry_reconfigure(master_t *master, geometry_t **window)
{
  int i, bn, m, n;
  char *files[MAXSTRLEN * MAXNUMARGS];
  char bstd[MAXSTRLEN], buf[MAXSTRLEN], buf1[MAXSTRLEN];
  geometry_t *geom = master->geom;

  strcpy(buf, master->runerror);
  m = parseline(buf, files, MAXNUMARGS);

  /* Open boundary received a DF_FAIL error code : chage OBCs        */
  if (strcmp(files[0], "OBC") == 0) {
    open_bdrys_t *open;
    n = atoi(files[1]);
    open = geom->open[n];
    if (open->bstdf == -1 || open->nbstd <= 1) {
      hd_warn("Reconfigure: standard boundaries not in use; can't reconfigure OBCs\n");
      return;
    }
    sprintf(bstd, "%c", '\0');
    for (i = 2; i < m; i++) {
      /* Find any 2-way nesting standard boundary conditions         */
      if (strcmp(files[i], "NEST2WAY") == 0) {
	int nn;
	char *list[MAXSTRLEN * MAXNUMARGS];
	for (bn = 0; bn < open->nbstd; bn++) {
	  strcpy(buf1, open->bstd[bn]);
	  nn = parseline(buf1, list, MAXNUMARGS);
	  if (strcmp(list[0], "NEST2WAY") == 0) {
	    strcpy(bstd, open->bstd[bn]);
	    open->bstdf = bn;
	    break;
	  }
	}
      }
      /* Find any 1-way nesting standard boundary conditions         */
      else if (strcmp(files[i], "NEST1WAY") == 0) {
	int nn;
	char *list[MAXSTRLEN * MAXNUMARGS];
	for (bn = 0; bn < open->nbstd; bn++) {
	  strcpy(buf1, open->bstd[bn]);
	  nn = parseline(buf1, list, MAXNUMARGS);
	  if (strcmp(list[0], "NEST1WAY") == 0) {
	    strcpy(bstd, open->bstd[bn]);
	    open->bstdf = bn;
	    break;
	  }
	}
      }
      /* Find any river standard boundary conditions                 */
      if (strcmp(files[i], "RIVER") == 0) {
	int nn;
	char *list[MAXSTRLEN * MAXNUMARGS];
	for (bn = 0; bn < open->nbstd; bn++) {
	  strcpy(buf1, open->bstd[bn]);
	  nn = parseline(buf1, list, MAXNUMARGS);
	  if (strcmp(list[0], "RIVER") == 0) {
	    strcpy(bstd, open->bstd[bn]);
	    open->bstdf = bn;
	    break;
	  }
	}
      }
    }
    /* Switch between standard conditions 0 and 1                    */
    if (!strlen(bstd)) {
      if (open->bstdf == 0) {
	strcpy(bstd, open->bstd[1]);
	open->bstdf = 1;
      } else if (open->bstdf == 1) {
	strcpy(bstd, open->bstd[0]);
	open->bstdf = 0;
      }
    }
    strcpy(buf, bstd);
    m = parseline(buf, files, MAXNUMARGS);
    hd_warn("Reconfigure: reconfiguring OBC %s to %s at %.2f days\n", open->name, 
	    files[0], master->days);
    if (bdry_reinit(geom, master, window, open, bstd))
      hd_quit("OBC reinit: no alternative OBC supplied for bdry %d\n", n);
  }
  master->regf &= ~RS_OBCSET;
}

/* END bdry_reconfigure()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to re-initialise an open boundary configuration           */
/*-------------------------------------------------------------------*/
int bdry_reinit(geometry_t *geom,        /* Global geometry          */
		master_t *master,        /* Master data              */		 
		geometry_t **window,     /* Window geometry          */
		open_bdrys_t *open,      /* Open boundary structure  */
		char *bstd               /* Standard OBC             */
		)
{
  int n, bn, t;
  char buf[MAXSTRLEN];

  /* Free memory associated with this boundary                       */
  bdry_custom_free(geom, master, window, open);

  /* Reset the open boundary condition                               */
  if (open->bstdf >= 0) {

    strcpy(buf, bstd);
    std_bdry_init(open, buf, master->grid_dt, master->dt2d, 
		  master->bdrypath, master->trinfo_3d);

    /* Set up the custom specification on the master                 */
    for (t = 0; t < open->ntr; ++t) {
      parameters_t *params = master->params;
      open_bdrys_t *io = params->open[open->id];
      open->cusname_t[t] = (char *)malloc(MAXSTRLEN);
      strcpy(open->cusname_t[t], io->cusname_t[t]);
    }
    bdry_custom_init(master, open, master->trinfo_3d, open->custype);
    bdry_init_m(master);

    /* Set up the custom routines on the windows                     */
    for (n = 1; n <= geom->nwindows; n++) {
      for (bn = 0; bn < window[n]->nobc; bn++) {
	open_bdrys_t *open_w = window[n]->open[bn];
	if (open->id == open_w->id) {
	  /* Copy the boundary specification to the windows          */
	  open_w->bcond_ele = open->bcond_ele;
	  open_w->bcond_nor = open->bcond_nor;
	  open_w->bcond_nor2d = open->bcond_nor2d;
	  open_w->bcond_tan = open->bcond_tan;
	  open_w->bcond_tan2d = open->bcond_tan2d;
	  open_w->bcond_Vz = open->bcond_Vz;
	  open_w->bcond_Kz = open->bcond_Kz;
	  open_w->adjust_flux = open->adjust_flux;
	  open_w->adjust_flux_s = open->adjust_flux_s;
	  open_w->inverse_barometer = open->inverse_barometer;
	  strcpy(open_w->bflow, open->bflow);
	  open_w->bhc = open->bhc;
	  open_w->rlen = open->rlen;
	  open_w->options = open->options;
	  strcpy(open_w->tsfn, open->tsfn);
	  strcpy(open_w->cusname_u1, open->cusname_u1);
	  strcpy(open_w->cusname_u2, open->cusname_u2);
	  for (t = 0; t < open_w->ntr; t++)
	    open_w->bcond_tra[t] = open->bcond_tra[t]; 	  
	  /* These are included for consistency                      */
	  open_w->ntsfiles = open->ntsfiles;
	  open_w->tsfiles = open->tsfiles;
	  open_w->filenames = open->filenames;
	  /* u1 or u2 custom boundary routines                       */
	  if (strlen(open->datau1.name))
	    bdry_data_copy(&open_w->datau1, &open->datau1);
	  if (strlen(open->datau2.name))
	    bdry_data_copy(&open_w->datau2, &open->datau2);
	  if (strlen(open->datau1av.name))
	    bdry_data_copy(&open_w->datau1av, &open->datau1av);
	  if (strlen(open->datau2av.name))
	    bdry_data_copy(&open_w->datau2av, &open->datau2av);
	  if (strlen(open->etadata.name))
	    bdry_data_copy(&open_w->etadata, &open->etadata);
	  /* Tracer custom boundary routines                         */
	  open_w->bdata_t = (bdry_details_t *)
	    malloc(sizeof(bdry_details_t) * open->ntr);
	  for (t = 0; t < open->ntr; t++)
	    bdry_data_copy(&open_w->bdata_t[t], &open->bdata_t[t]);
	  /* Initialise the window custom boundary function */
	  bdry_init_w(window[n], open_w, open);
	}
      }
    }
  } else
    return 1;
  return 0;
}

/* END bdry_reinit()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to assign tracer numbers of scaling functions             */
/*-------------------------------------------------------------------*/
void read_bdry_scale(open_bdrys_t *io,    /* Open boundary strcuture */
		     open_bdrys_t *open,  /* Open boundary strcuture */
		     char *trname[]       /* Tracer names            */
		     )
{
  int t, tn;
  int ret;

  open->sdata_t = (scale_details_t *)malloc(sizeof(scale_details_t) *
					    open->ntr);
  memset(open->sdata_t, 0, sizeof(scale_details_t) * open->ntr);

  for (t = 0; t < open->ntr; t++) {
    scale_details_t *scale = &open->sdata_t[t];
    sprintf(scale->name, "%c", '\0');

    /* These are all stack variables of MAXSTRLEN length, check first index */
    if (*io->scale_s == NULL || *io->scale_p == NULL || *io->scale_d == NULL)
      continue;

    /* Scaling to the density gradient */
    if (strlen(io->scale_d[t])) {
      char *fields[MAXSTRLEN * MAXNUMARGS];
      int nf = parseline(io->scale_d[t], fields, MAXNUMARGS);
      scale->type |= TRSC_DEN;
      scale->flag = 0;
      scale->val = atof(fields[0]);
      scale->fact = atof(fields[1]);
      if (nf == 3) {
	if (strcmp(fields[2], "c") == 0) {
	    scale->flag |= TRSC_CPY;
	    sprintf(scale->name, "Density gradient scaled: FILEIN at depth %5.2f and below, gradient multiplier=%5.2f",
		    scale->val, scale->fact);
	}
	if (strcmp(fields[2], "t") == 0) {
	  scale->flag |= TRSC_TRN;
	  sprintf(scale->name, "Density gradient scaled: endpoint=FILEIN at depth %5.2f (truncated below), gradient multiplier=%5.2f",
		  scale->val, scale->fact);
	}
      } else
	sprintf(scale->name, "Density gradient scaled: FILEIN at depth %5.2f, gradient multiplier=%5.2f",
		scale->val, scale->fact);
      scale->ntr = -1;
      continue;
    }

    if (strlen(io->scale_s[t])) {
      scale->type |= TRSC_SUM;
      strcpy(scale->name, io->scale_s[t]);
    } else if (strlen(io->scale_p[t])) {
      scale->type |= TRSC_PCT;
      strcpy(scale->name, io->scale_p[t]);
    } else 
      continue;

    if (strlen(scale->name)) {
      char *fields[MAXSTRLEN * MAXNUMARGS];
      int nf = parseline(scale->name, fields, MAXNUMARGS);

      /*-------------------------------------------------------------*/
      /* Check if it is a function                                   */
      if (nf <= 0)
	hd_quit
	  ("readbdryscale : bdry %s tracer %d scaling incorrectly specified.\n",
	   open->name, t);

      /* Not implemented
      locate_scale_function(open, fields[0], &scale->custom_w);
      scale->type |= TRSC_CUS;
      */
      scale->ntr = -1;
      scale->fact = 0.0;
      /* Check if it is a tracer name or constant                    */

      if (nf == 1) {
	scale->type |= TRSC_NUM;  /* Assume it is a constant         */	
	ret = 0;
	for (tn = 0; tn < open->ntr; ++tn) {
	  if (strcmp(scale->name, trname[tn]) == 0) {
	    ret = 1;
	    break;
	  }
	}
	if (ret) {
	  scale->type &= ~TRSC_NUM;
	  scale->type |= TRSC_TRA;
	  scale->ntr = tn;
	} else if (sscanf(fields[0], "%lf", &scale->fact) != 1)
          hd_quit
	    ("readbdryscale : bdry %s tracer %d scaling incorrectly specified.\n",
	     open->name, t);
      } else
	hd_quit
	  ("readbdryscale : bdry %s tracer %d scaling incorrectly specified.\n",
	   open->name, t);
    }
  }
}

/* END read_bdry_scale()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to assign tracer numbers of scaling functions             */
/*-------------------------------------------------------------------*/
void read_bdry_eta_scale(master_t *master,   /* Master data             */
			 open_bdrys_t *io,   /* Open boundary strcuture */
			 open_bdrys_t *open  /* Open boundary strcuture */
			 )
{
  int ret, tn;
  scale_details_t *scale;

  open->sdata_e = (scale_details_t *)malloc(sizeof(scale_details_t));
  memset(open->sdata_e, 0, sizeof(scale_details_t));

  scale = open->sdata_e;

  if (io->scale_e == NULL)
    return;

  if (strlen(io->scale_e)) {
    scale->type |= TRSC_SUM;
    strcpy(scale->name, io->scale_e);
  } else 
    return;

  if (strlen(scale->name)) {
    char *fields[MAXSTRLEN * MAXNUMARGS];
    int nf = parseline(scale->name, fields, MAXNUMARGS);

    /*-------------------------------------------------------------*/
    /* Check if it is a function                                   */
    if (nf <= 0)
      hd_quit
	("readbdryscale : bdry %s eta scaling incorrectly specified.\n",
	 open->name);

    /* Not implemented
       locate_scale_function(open, fields[0], &scale->custom_w);
       scale->type |= TRSC_CUS;
    */
    scale->ntr = -1;
    scale->fact = 0.0;
    /* Check if it is a tracer name or constant                    */
    if (nf == 1) {
      scale->type |= TRSC_NUM;  /* Assume it is a constant         */	
      ret = 0;
      for (tn = 0; tn < master->ntrS; ++tn) {
	if (strcmp(scale->name, master->trinfo_2d[tn].name) == 0) {
	  ret = 1;
	  break;
	}
      }
      if (ret) {
	scale->type &= ~TRSC_NUM;
	scale->type |= TRSC_TRA;
	scale->ntr = tn;
      } else if (sscanf(fields[0], "%lf", &scale->fact) != 1)
	hd_quit
	  ("readbdryscale : bdry %s eta scaling incorrectly specified.\n",
	   open->name);
    } else if (nf == 2 && strcmp(fields[0], "*") == 0) {
      scale->type &= ~TRSC_SUM;
      scale->type |= TRSC_PCT;
      scale->type |= TRSC_NUM;  /* Assume it is a constant         */	
      ret = 0;
      for (tn = 0; tn < master->ntrS; ++tn) {
	if (strcmp(scale->name, master->trinfo_2d[tn].name) == 0) {
	  ret = 1;
	  break;
	}
      }
      if (ret) {
	scale->type &= ~TRSC_NUM;
	scale->type |= TRSC_TRA;
	scale->ntr = tn;
      } else if (sscanf(fields[1], "%lf", &scale->fact) != 1)
	hd_quit
	  ("readbdryscale : bdry %s eta scaling incorrectly specified.\n",
	   open->name);
    } else
      hd_quit
	("readbdryscale : bdry %s eta scaling incorrectly specified.\n",
	 open->name);
  }
}

/* END read_bdry_eta_scale()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read in boundary cell lists into the params open       */
/* boundary structure. These are subsequently converted to the OBC   */
/* specification in the mesh structure.                              */
/*-------------------------------------------------------------------*/
void get_obc_list(open_bdrys_t *open, FILE *fp, int n, char *key) {
  int nn, m;
  int is, ie, js, je;
  char key1[MAXSTRLEN], key2[MAXSTRLEN], key3[MAXSTRLEN], key4[MAXSTRLEN];
  char buf[MAXSTRLEN], key5[MAXSTRLEN];

  prm_set_errfn(hd_warn);
  sprintf(key1, "%s%1d.POINTS", key, n);
  sprintf(key2, "%s%1d.RANGE", key, n);
  sprintf(key3, "%s%1d.LIMITS", key, n);
  sprintf(key4, "%s%1d.UPOINTS", key, n);
  sprintf(key5, "%s%1d.START_LOC", key, n);

  if (prm_read_int(fp, key1, &open->npts)) {
    open->intype = O_POI;
    open->iloc = i_alloc_1d(open->npts);
    open->jloc = i_alloc_1d(open->npts);
    for (nn = 0; nn < open->npts; nn++) {
      if (fscanf(fp, "%d %d", &open->iloc[nn],
		 &open->jloc[nn]) != 2)
	hd_quit("params_read: Can't read i j in boundary  points list.\n");
      /* Flush the remainder of the line */
      prm_flush_line(fp);
    }
  } else if (prm_read_char(fp, key2, buf)) {
    open->intype = O_RAN;
    sscanf(buf, "(%d,%d)-(%d,%d)", &is, &js, &ie, &je);
    if (ie-is+1 <= 0)
      hd_quit("params_read: start > end i coordinate for boundary %d (%d != %d).\n", n, is, ie);
    else if (je-js+1 <= 0)
      hd_quit("params_read: start > end j coordinate for boundary %d (%d != %d).\n", n, js, je);
    else {
      int dir = (ie - is + 1 > je - js + 1) ? 1 : 0;
      open->npts = max(ie - is + 1, je - js + 1);
      open->iloc = i_alloc_1d(open->npts);
      open->jloc = i_alloc_1d(open->npts);
      for (nn = 0; nn < open->npts; nn++) {
	if (dir) {
	  open->iloc[nn] = is + nn;
	  open->jloc[nn] = js;
	} else {
	  open->iloc[nn] = is;
	  open->jloc[nn] = js + nn;
	}
      }
    }
  } else if (prm_read_char(fp, key3, buf)) {
    double d1;
    open->intype = O_LIM;
    sscanf(buf, "(%lf,%lf)-(%lf,%lf)", &open->minlon, &open->minlat, 
	   &open->maxlon, &open->maxlat);
    if (open->minlat > open->maxlat) {
      d1 = open->minlat;
      open->minlat = open->maxlat;
      open->maxlat = d1;
    }
    if (open->minlon > open->maxlon) {
	d1 = open->minlon;
	open->minlon = open->maxlon;
	open->maxlon = d1;
    }
    sprintf(key1, "%s%1d.EDGES", key, n);
    if (prm_read_char(fp, key1, buf)) {
      char *fields[MAXSTRLEN * MAXNUMARGS];
      open->nedges = parseline(buf, fields, MAXNUMARGS);
      open->edges = i_alloc_1d(open->nedges);
      for (nn = 0; nn < open->nedges; nn++)
	open->edges[nn] = atoi(fields[nn]);
    }
  } else if (prm_read_int(fp, key4, &open->npts)) {
    /* If an unstructured mesh specification is read in read_mesh_us */
    /* (US_IUS is set) then OBCs in the mesh structure have been     */
    /* directly read in.                                             */
    char *fields[MAXSTRLEN * MAXNUMARGS];
    open->locu = i_alloc_1d(open->npts);
    open->posx = d_alloc_2d(2 ,open->npts);
    open->posy = d_alloc_2d(2 ,open->npts);
    prm_flush_line(fp);
    fgets(buf,56,fp);
    m = parseline(buf, fields, MAXNUMARGS);
    prm_skip_to_end_of_key(fp, key4);
    prm_flush_line(fp);
    for (nn = 0; nn < open->npts; nn++) {
      open->posx[nn][1] = open->posy[nn][1] = NOTVALID;
      if (m == 3) {
	open->intype = O_UPI;
	if (fscanf(fp, "%d (%lf %lf)", 
		   &open->locu[nn], &open->posx[nn][0], &open->posy[nn][0]) != 3)
	  hd_quit("params_read: Can't read edge coordinates in boundary points list. Format '1 (2 3)'\n");
      }
      if (m == 2) {
	open->intype = O_UPC;
	if (fscanf(fp, "%d (%lf,%lf)-(%lf,%lf)", 
		   &open->locu[nn], &open->posx[nn][0], &open->posy[nn][0],
		   &open->posx[nn][1], &open->posy[nn][1]) != 5)
	  hd_quit("params_read: Can't read edge coordinates in boundary points list. Format '1 (1.0,2.0)-(3.0,4.0)\n");
      }
      /* Flush the remainder of the line */
      prm_flush_line(fp);
    }
  } else if (prm_read_char(fp, key5, buf)) {
    open->intype = O_COR;
    sscanf(buf, "%lf %lf", &open->slon, &open->slat);
    sprintf(key5, "%s%1d.END_LOC", key, n);
    if (prm_read_char(fp, key5, buf)) {
      sscanf(buf, "%lf %lf", &open->elon, &open->elat);
      sprintf(key5, "%s%1d.MID_LOC", key, n);
      if (prm_read_char(fp, key5, buf)) {
	sscanf(buf, "%lf %lf", &open->mlon, &open->mlat);
      }
    } else {
      hd_quit("params_read: cannot find %s in OBC%d (%s)\n", key5, n, open->name);
    }
  } else {
    hd_quit("params_read: cannot find %s%d (%s) , UPOINTS, POINTS, RANGE or START_LOC.\n", key, n, open->name);
  }

  /* Read the cell locations to map CYCLED boundaries to */
  prm_set_errfn(hd_silent_warn);
  sprintf(key1, "%s%1d.CYCLED_POINTS", key, n);
  if (prm_read_int(fp, key1, &open->ncyc)) {
    if (open->ncyc != open->npts) {
      hd_quit("params_read: CYCLED boundary  points must = %d for boundary %d list.\n", open->npts, n);
    }
    open->ilocc = i_alloc_1d(open->ncyc);
    open->jlocc = i_alloc_1d(open->ncyc);
    for (nn = 0; nn < open->ncyc; nn++) {
      if (fscanf(fp, "%d %d", &open->ilocc[nn],
		 &open->jlocc[nn]) != 2)
	hd_quit("params_read: Can't read i j in CYCLED boundary  points list.\n");
      /* Flush the remainder of the line */
      prm_flush_line(fp);
    }
  }

  sprintf(key1, "%s%1d.CYCLED_RANGE", key, n);
  if (prm_read_char(fp, key1, buf)) {
    sscanf(buf, "(%d,%d)-(%d,%d)", &is, &js, &ie, &je);
    if (ie-is+1 <= 0)
      hd_quit("params_read: start > end i coordinate for cycled boundary %d (%d != %d).\n", n, js, je);
    else if (je-js+1 <= 0)
      hd_quit("params_read: start > end j coordinate for cycled boundary %d (%d != %d).\n", n, js, je);
    else {
      int dir = (ie - is + 1 > je - js + 1) ? 1 : 0;
      open->ncyc = max(ie - is + 1, je - js + 1);
      open->ilocc = i_alloc_1d(open->ncyc);
      open->jlocc = i_alloc_1d(open->ncyc);
      for (nn = 0; nn < open->ncyc; nn++) {
	if (dir) {
	  open->ilocc[nn] = is + nn;
	  open->jlocc[nn] = js;
	} else {
	  open->ilocc[nn] = is;
	  open->jlocc[nn] = js + nn;
	}
      }
    }
    if (open->ncyc != open->npts) {
      hd_quit("params_read: CYCLED boundary  points must = %d for boundary %d list.\n", open->npts, n);
    }
  }
  prm_set_errfn(hd_quit);
}

/* END get_obc_list()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Converts an open boundary specification in the params->open       */
/* structure to the mesh OBC specification.                          */
/* Conversion takes place according to the following:                */
/* For -g mode:                                                      */
/* O_POI and O_RAN are converted in this routine called from         */
/*           meshstruct_s().                                         */
/* O_UPI and O_UPC are converted in this routine called from         */
/*           meshstruct_us().                                        */
/* O_COR are converted in get_mesh_obc() (called within              */
/*           build_sparse_grid_us().                                 */
/* O_LIM are converted in get_mesh_obc_limit() and has been          */
/*           superseded with O_COR.                                  */
/* For -p mode:                                                      */
/* O_POI and O_RAN are specified in the params->open structure in    */
/*           convert_structured_obc() (called from create_mesh())    */
/*           in a format corresponding to O_UPC. These are then      */
/*           converted in this routine called from meshstruct_us()   */
/*           (which is also called from create_mesh()).              */
/* O_UPI and O_UPC are converted in this routine, called from        */
/*           meshstruct_us() within create_mesh()).                  */
/* O_COR are converted in get_mesh_obc() (called within              */
/*           build_sparse_grid_us().                                 */
/*-------------------------------------------------------------------*/
void convert_obc_list(parameters_t *params, /* Parameters info       */
		      open_bdrys_t *open,   /* Open boundary info    */
		      int n,                /* Open boundary number  */
		      geometry_t *geom,     /* Structured geometry   */
		      int *cmap             /* Mapping function      */
		      )
		      
{
  FILE *op, *sp;
  mesh_t *mesh = params->mesh;
  int m, c, cc, cco;
  double xb1, yb1, xb2, yb2;
  double eps = 1e-6;    /* Precision for OBC comparisons             */
  int verbose = 0;
  int filef = 1;
  int obcf = 1;

  if (open->intype & O_COR) {
    /* OBC specification from START & END coordinates is computed in */
    /* get_mesh_obc() (called within build_sparse_grid_us().         */
    return;
  }

  if (filef) {
    if (n == 0) {
      if ((op = fopen("boundary.txt", "w")) == NULL)
	filef = 0;
    } else {
      if ((op = fopen("boundary.txt", "a")) == NULL)
	filef = 0;
    }
  }
  if (obcf) {
    if (n == 0) {
      if ((sp = fopen("obc_spec.txt", "w")) == NULL)
	obcf = 0;
    } else {
      if ((sp = fopen("obc_spec.txt", "a")) == NULL)
	obcf = 0;
    }
    fprintf(sp, "\nBOUNDARY%d.UPOINTS     %d\n", n, mesh->npts[n]);
  }

  if (open->intype & (O_UPI|O_UPC)) {
    if (open->posx == NULL || open->posy == NULL) return;
    for (cc = 0; cc < open->npts; cc++) {
      /* Note: open->npts loops from 0:npts-1, mesh->npts loops from */
      /* 1:npts.                                                     */
      xb1 = open->posx[cc][0];
      yb1 = open->posy[cc][0];
      xb2 = open->posx[cc][1];
      yb2 = open->posy[cc][1];
      cco = open->locu[cc];	
      c = cmap[cco];
      if (open->intype & O_UPI) {
	/* OBC indices do not need to be remapped in                 */
	/* convert_mesh_input(); set mesh->loc[n][0] = NOTVALID to   */
	/* indicate this.                                            */
	mesh->loc[n][cc+1] = c;
	mesh->obc[n][cc+1][0] = (int)xb1;
	mesh->obc[n][cc+1][1] = (int)yb1;
	mesh->loc[n][0] = NOTVALID;
	if (verbose) printf("OBC%d UPI %d : cc=%d (%d %d)\n",n, cc+1, cco, 
			    mesh->obc[n][cc+1][0], mesh->obc[n][cc+1][1]);
	if (filef) {
	  fprintf(op, "%f %f\n", mesh->xloc[(int)xb1], mesh->yloc[(int)xb1]);
	  fprintf(op, "%f %f\n", mesh->xloc[(int)yb1], mesh->yloc[(int)yb1]);
	  fprintf(op, "NaN NaN\n");
	}
	if (obcf) {
	  fprintf(sp, "%d (%lf,%lf)-(%lf,%lf)\n", mesh->loc[n][cc+1], 
		  mesh->xloc[(int)xb1], mesh->yloc[(int)xb1],
		  mesh->xloc[(int)yb1], mesh->yloc[(int)yb1]);
	}
      }
      if (open->intype & O_UPC) {
	double x1, y1, x2, y2;
	int j, jj;
	for (j = 1; j <= mesh->npe[c]; j++) {
	  double x1, y1, x2, y2;
	  x1 = params->x[c][j];
	  y1 = params->y[c][j];
	  jj = (j == mesh->npe[c]) ? 1 : j + 1;
	  x2 = params->x[c][jj];
	  y2 = params->y[c][jj];
	  if ((fabs(xb1-x1) < eps && fabs(yb1-y1) < eps && 
	       fabs(xb2-x2) < eps && fabs(yb2-y2) < eps) ||
	      (fabs(xb1-x2) < eps && fabs(yb1-y2) < eps && 
	       fabs(xb2-x1) < eps && fabs(yb2-y1) < eps)) {
	    mesh->loc[n][cc+1] = c;
	    mesh->obc[n][cc+1][0] = j;
	    mesh->obc[n][cc+1][1] = jj;
	    if (verbose) printf("OBC%d %d : cc=%d (%d %d)\n",n, cc+1, cco, j, jj); 
	    if (filef) {
	      fprintf(op, "%f %f\n", xb1, yb1);
	      fprintf(op, "%f %f\n", xb2, yb2);
	      fprintf(op, "NaN NaN\n");
	    }
	  }
	}
      }
    }
  }

  if (open->intype & (O_POI|O_RAN)) {
    mesh->npts[n] = open->no2_t;
    for (cc = 1; cc <= open->no2_t; cc++) {
      c = open->obc_t[cc];
      cco = geom->c2cc[c];
      mesh->loc[n][cc] = cco;
      if (open->type & U1BDRY && open->ocodex == L_EDGE) {
	mesh->obc[n][cc][0] = 1;
	mesh->obc[n][cc][1] = 2;
      }
      if (open->type & U2BDRY && open->ocodey == B_EDGE) {
	mesh->obc[n][cc][0] = 4;
	  mesh->obc[n][cc][1] = 1;
      }
      if (open->type & U1BDRY && open->ocodex == R_EDGE) {
	mesh->obc[n][cc][0] = 3;
	mesh->obc[n][cc][1] = 4;
      }
      if (open->type & U2BDRY && open->ocodey == F_EDGE) {
	mesh->obc[n][cc][0] = 2;
	mesh->obc[n][cc][1] = 3;
      }
    }
  }
  if (filef) fclose(op);
  if (obcf) fclose(sp);
}

/* END convert_obc_list()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets the flux adjustment timescale                                */
/*-------------------------------------------------------------------*/
void get_flux_adjust(FILE *fp, open_bdrys_t *open, int n) {

  char buf[MAXSTRLEN], buf1[MAXSTRLEN];
  char keyword[MAXSTRLEN];
  char tu0[MAXSTRLEN], tu1[MAXSTRLEN];
  double r0, r1, d0, d1;

  sprintf(keyword, "BOUNDARY%1d.ADJUST_FLUX", n);
  if (prm_read_char(fp, keyword, buf) > 0) {
    open->fas = (fa_info_t *)malloc(sizeof(relax_info_t));
    memset(open->fas, 0, sizeof(relax_info_t));
    if (sscanf(buf, "%s %lf %lf %s %lf %lf %s", 
	       buf1, &d0, &r0, tu0, &d1, &r1, tu1) == 7) {
      
      if(strcmp(buf1, "linear") == 0) {
	open->fas->dv0 = d0;
	open->fas->dv1 = d1;
	sprintf(buf, "%f %s", r0, tu0);
	tm_scale_to_secs(buf, &r0);
	sprintf(buf, "%f %s", r1, tu1);
	tm_scale_to_secs(buf, &r1);
	open->fas->tc0 = r0;
	open->fas->tc1 = r1;
	open->fas->slope = (d1 - d0) ? (r1 - r0) / (d1 - d0) : 0.0;
	open->fas->tctype = (RLX_ADPT|RLX_LINR);
      } else if(strcmp(buf, "temporal") == 0) {
	open->fas->dv0 = d0;
	open->fas->dv1 = d1;
	sprintf(buf, "%f %s", r0, tu0);
	tm_scale_to_secs(buf, &r0);
	sprintf(buf, "%f %s", r1, tu1);
	tm_scale_to_secs(buf, &r1);
	open->fas->tc0 = r0;
	open->fas->tc1 = r1;
	open->fas->slope = (d1 - d0) ? (r1 - r0) / (d1 - d0) : 0.0;
	open->fas->tctype = (RLX_ADPT|RLX_TIME);
      }
    } else if (sscanf(buf, "%s %lf %lf %s", buf1, &d0, &r0, tu0) == 4) {
      if (strcmp(buf1, "exponential") == 0) {
	sprintf(buf, "%f %s", r0, tu0);
	tm_scale_to_secs(buf, &r0);
	open->fas->dv0 = d0 * log(r0);
	open->fas->dv1 = d0;
	open->fas->slope = d0 * log(r0 / 86400.0);
	open->fas->tc0 = r0;
	open->fas->tctype = (RLX_ADPT|RLX_EXP);
      }
    } else {
      tm_scale_to_secs(buf, &open->adjust_flux);
      open->fas->tctype = RLX_CONS;
      open->fas->tc0 = open->adjust_flux;
    }
  }
}

/* END get_flux_adjust()                                             */
/*-------------------------------------------------------------------*/

tidal_memory_t **tide_alloc_2d(long n1, long n2)
{
  tidal_memory_t *p;
  tidal_memory_t **pp;
  long i, size;

  /* first allocate main storage */
  size = n1 * n2 * sizeof(tidal_memory_t);
  if ((p = (tidal_memory_t *)malloc(size)) == NULL)
    hd_quit("tide_alloc_2d: Not enough memory\n");

  /* now allocate row pointers */
  size = n2 * sizeof(tidal_memory_t *);
  if ((pp = (tidal_memory_t **)malloc(size)) == NULL)
    hd_quit("tide_alloc_2d: Can't allocate row pointers\n");
  /* point to the rows */
  for (i = 0; i < n2; i++)
    pp[i] = p + i * n1;
  return (pp);
}

/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to initialise standard boundaries                         */
/*-------------------------------------------------------------------*/
void std_bdry_init(open_bdrys_t *open, char *stdbdry, 
		   double dt, double dt2d, char *bdrypath, tracer_info_t *tracer)
{
  int m;
  char buf[MAXSTRLEN];
  char *files[MAXSTRLEN * MAXNUMARGS];

  strcpy(buf, stdbdry);
  m = parseline(buf, files, MAXNUMARGS);
  if (contains_token(files[0], "NEST1WAY") != NULL) {
    if (m < 4) hd_quit("1-WAY_NEST format : 'NEST1WAY data_ets.nc data_uv_nor.nc data_uv_tan.nc'\n");
    std_bdry(open, m, files, tracer, dt, dt2d, bdrypath, NEST1WAY);
  }
  if (contains_token(files[0], "NEST1WAY|STD") != NULL) {
    if (m < 4) hd_quit("1-WAY_NEST format : 'NEST1WAY|STD data_ets.nc data_uv_nor.nc data_uv_tan.nc'\n");
    std_bdry(open, m, files, tracer, dt, dt2d, bdrypath, NEST1WAY|NEST_STD);
  }
  if (contains_token(files[0], "NEST2WAY") != NULL) {
    if (m < 4) hd_quit("2-WAY_NEST format : 'NEST2WAY data_ets.mpk data_uv_nor.mpk data_uv_tan.mpk'\n");
    std_bdry(open, m, files, tracer, dt, dt2d, bdrypath, NEST2WAY);
  }
  if (contains_token(files[0], "NEST_FLA") != NULL) {
    if (m < 3) hd_quit("Flather radiation format : 'NEST_FLA data_ets.nc data_uv_nor.nc data_uvav_tan.nc'\n");
    std_bdry(open, m, files, tracer, dt, dt2d, bdrypath, NEST_FLA);
  }
  if (contains_token(files[0], "NEST_FLA|STD") != NULL) {
    if (m < 3) hd_quit("Flather radiation format : 'FLATHR data_ets.nc data_uvav_nor.nc data_uvav_tan.nc'\n");
    std_bdry(open, m, files, tracer, dt, dt2d, bdrypath, NEST_FLA|NEST_STD);
  }
  if (contains_token(files[0], "NEST_RAD") != NULL) {
    if (m < 2) hd_quit("Sommerfeld radiation format : 'NEST_RAD data_ets.nc'\n");
    std_bdry(open, m, files, tracer, dt, dt2d, bdrypath, NEST_RAD);
  }
  if (contains_token(files[0], "NEST_CPD") != NULL) {
    if (m < 4) hd_quit("CLAMPED nesting format : 'NEST_CPD data_ets.mpk data_uv_nor.mpk data_uv_tan.mpk'\n");
    std_bdry(open, m, files, tracer, dt, dt2d, bdrypath, NEST_CPD);
  }
  if (contains_token(files[0], "NEST_CPD|STD") != NULL) {
    if (m < 4) hd_quit("CLAMPED nesting format : 'NEST_CPD data_ets.nc data_uv_nor.nc data_uv_tan.nc'\n");
    std_bdry(open, m, files, tracer, dt, dt2d, bdrypath, NEST_CPD|NEST_STD);
  }
  if (contains_token(files[0], "RIVER") != NULL) {
    if (m < 3)
      hd_quit("RIVER format : 'RIVER flow.ts temp.ts\n");
    std_bdry(open, m, files, tracer, dt, dt2d, bdrypath, RIVER);
  }
  if (contains_token(files[0], "SOLID") != NULL) {
    std_bdry(open, m, files, tracer, dt, dt2d, bdrypath, SOLID);
  }
  if (contains_token(files[0], "NOTHIN") != NULL) {
    std_bdry(open, m, files, tracer, dt, dt2d, bdrypath, NOTHIN);
  }
}

/* END std_bdry_init()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets up standard boundaries                                       */
/*-------------------------------------------------------------------*/
void std_bdry(open_bdrys_t *open,	       
	      int nfiles,
	      char **file,
	      tracer_info_t *tracers, 
	      double dt,
	      double dt2d,
	      char *bdrypath,
	      int mode)
{
  int i, tid, sid;
  char buf[MAXSTRLEN];
  char keyword[MAXSTRLEN];
  char cusu1[MAXSTRLEN];
  char cusu2[MAXSTRLEN];

  open->custype = open->type;

  /*-----------------------------------------------------------------*/
  /* Solid boindaries or no action taken                             */
  if (mode & (SOLID|NOTHIN)) {
    open->sbcond = NOTHIN;
    open->bcond_ele = NOTHIN;
    open->bcond_nor = NOTHIN;
    open->bcond_nor2d = NOTHIN;
    if (mode & SOLID) {
      open->sbcond = SOLID;
      open->bcond_nor = CLAMPD;
      open->bcond_nor2d = CLAMPD;
    }
    open->bcond_tan = NOTHIN;
    open->bcond_tan2d = NOTHIN;
    if (open->bcond_Vz == NOTHIN) open->bcond_Vz = NOGRAD;
    if (open->bcond_Kz == NOTHIN) open->bcond_Kz = NOGRAD;
    if (nfiles > 1)
      strcpy(open->tsfn, file[1]);
    for (i = 2; i < nfiles; i++)
      sprintf(open->tsfn, "%s %s",open->tsfn, file[i]);      
    if (strlen(bdrypath)) {
      sprintf(buf, "%s%s", bdrypath, open->tsfn);
      strcpy(open->tsfn, buf);
    }
    for (i = 0; i < open->ntr; i++) {
      if (open->bcond_tra[i] == NOGRAD)
	open->bcond_tra[i] = TRCONC|NOTHIN;
    }
    return;
  }

  /*-----------------------------------------------------------------*/
  /* 1 way nesting                                                   */
  if (mode & NEST1WAY) {
    int do_baro = 0;

    open->sbcond = NEST1WAY;
    if (do_baro)
      open->sbcond |= NEST_BARO;

    /* Sea level                                                     */
    if (open->bcond_ele == NOTHIN)
      open->bcond_ele = NOTHIN|FILEIN;

    /* Temperature and salinity                                      */
    for (i = open->atr; i < open->ntr; i++) {
      if (strcmp(tracers[i].name, "temp") == 0) {
	tid = i;
	if (open->bcond_tra[i] == NOGRAD)
	  open->bcond_tra[i] = TRCONC|FILEIN;
      } else if (strcmp(tracers[i].name, "salt") == 0) {
	sid = i;
	if (open->bcond_tra[i] == NOGRAD)
	  open->bcond_tra[i] = TRCONC|FILEIN;
      } else
	open->bcond_tra[i] = NOGRAD;
    }

    /* Velocity                                                      */
    open->bcond_nor = CUSTOM;
    open->bcond_nor2d = VERTIN;
    open->bcond_tan = CUSTOM;
    open->bcond_tan2d = VERTIN;
    open->bcond_Vz = NOGRAD;
    open->bcond_Kz = NOGRAD;
    if (open->adjust_flux == 0.0)
      open->adjust_flux = -0.2 * dt;
    /* Set a default tidal relaxation
    if (open->adjust_flux_s == 0.0)
      open->adjust_flux_s = -dt2d;
    */
    strcpy(open->tsfn, file[1]);
    for (i = 2; i < nfiles - 2; i++)
      sprintf(open->tsfn, "%s %s",open->tsfn, file[i]);
    if (strlen(bdrypath)) {
      sprintf(buf, "%s%s", bdrypath, open->tsfn);
      strcpy(open->tsfn, buf);
    }
    strcpy(cusu1, "uv_to_u1");
    strcpy(cusu2, "uv_to_u2");
    if (mode & NEST_STD) {
      strcpy(cusu1, "hdstd_to_u1");
      strcpy(cusu2, "hdstd_to_u2");
    }
    sprintf(open->cusname_u1, "%s %s", cusu1, file[nfiles-2]);
    sprintf(open->cusname_u2, "%s %s", cusu2, file[nfiles-1]);
  }

  /*-----------------------------------------------------------------*/
  /* 2 way nesting                                                   */
  if (mode & NEST2WAY) {
    int do_baro = 0;

    open->sbcond = NEST2WAY;
    if (do_baro)
      open->sbcond |= NEST_BARO;

    for (i = 1; i < nfiles; i++)
      if (!endswith(file[i], ".mpk"))
	hd_warn("nest2way: file %s is not the memory packet format required for 2-way nesting.\n", file[i]);

    /* Sea level                                                     */
    if (open->bcond_ele == NOTHIN)
      open->bcond_ele = NOTHIN|FILEIN;

    /* Temperature and salinity                                      */
    for (i = open->atr; i < open->ntr; i++) {
      if (strcmp(tracers[i].name, "temp") == 0) {
	tid = i;
	if (open->bcond_tra[i] == NOGRAD)
	  open->bcond_tra[i] = TRCONC|FILEIN;
	/*open->bcond_tra[i] = FILEIN;*/
      }
      if (strcmp(tracers[i].name, "salt") == 0) {
	sid = i;
	if (open->bcond_tra[i] == NOGRAD)
	  open->bcond_tra[i] = TRCONC|FILEIN;
	/*open->bcond_tra[i] = FILEIN;*/
      }
    }

    /* Velocity                                                      */
    open->bcond_nor = CUSTOM;
    open->bcond_nor2d = CUSTOM;
    open->bcond_tan = CUSTOM;
    open->bcond_tan2d = CUSTOM;
    open->bcond_Vz = NOGRAD;
    open->bcond_Kz = NOGRAD;

    if (open->adjust_flux == 0.0)
      open->adjust_flux = -0.2 * dt;
    /*
    if (open->adjust_flux_s == 0.0)
      open->adjust_flux_s = -dt2d;
    */
    if (strlen(bdrypath))
      sprintf(open->tsfn, "%s%s", bdrypath, file[1]);
    else
      strcpy(open->tsfn, file[1]);
    for (i = 2; i < nfiles - 2; i++) {
      if (strlen(bdrypath))
	sprintf(open->tsfn, "%s %s%s",open->tsfn, bdrypath, file[i]);
      else
	sprintf(open->tsfn, "%s %s",open->tsfn, file[i]);
    }
    strcpy(cusu1, "uv_to_u1");
    strcpy(cusu2, "uv_to_u2");
    if (mode & NEST_STD) {
      strcpy(cusu1, "hdstd_to_u1");
      strcpy(cusu2, "hdstd_to_u2");
    }
    sprintf(open->cusname_u1, "%s %s_uv_nor.mpk", cusu1, file[nfiles-2]);
    sprintf(open->cusname_u2, "%s %s_uv_tan.mpk", cusu2, file[nfiles-1]);
    if (do_baro) {
      strcpy(cusu1, "uvav_to_u1av");
      strcpy(cusu2, "uvav_to_u2av");
      if (mode & NEST_STD) {
	strcpy(cusu1, "hdstd_to_u1av");
	strcpy(cusu2, "hdstd_to_u2av");
      }
      sprintf(open->cusname_u1av, "%s %s_uvav_nor.mpk", cusu1, file[nfiles-2]);
      sprintf(open->cusname_u2av, "%s %s_uvav_tan.mpk", cusu2, file[nfiles-1]);
    } else {
      open->bcond_nor2d = VERTIN;
      open->bcond_tan2d = VERTIN;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Clamped nesting                                                 */
  if (mode & NEST_CPD) {
    int do_baro = 0;
    int suff = 0;

    if (open->options & OP_NESTBARO) do_baro = 1;
    open->sbcond = NEST_CPD;
    if (do_baro)
      open->sbcond |= NEST_BARO;

    for (i = 1; i < nfiles; i++) {
      if (!endswith(file[i], ".mpk"))
	hd_warn("clamped nesting: file %s is not the memory packet format required for clamped nesting.\n", file[i]);
      if (endswith(file[i], ".nc") || endswith(file[i], ".ts")) suff = 1;
    }

    /* Sea level                                                     */
    if (open->bcond_ele == NOTHIN) {
      open->bcond_ele = FILEIN|GRAVTY;

      open->relax_time = 900.0;
      open->relax_timei = 900.0;

    }

    /* Temperature and salinity                                      */
    for (i = open->atr; i < open->ntr; i++) {
      if (strcmp(tracers[i].name, "temp") == 0) {
	tid = i;
	if (open->bcond_tra[i] == NOGRAD)
	  open->bcond_tra[i] = FILEIN;
	/*open->bcond_tra[i] = FILEIN|UPSTRM;*/

      }
      if (strcmp(tracers[i].name, "salt") == 0) {
	sid = i;
	if (open->bcond_tra[i] == NOGRAD)
	  open->bcond_tra[i] = FILEIN;
	/*open->bcond_tra[i] = FILEIN|UPSTRM;*/
      }
    }

    /* Velocity                                                      */
    open->bcond_nor = CUSTOM;
    open->bcond_nor2d = CUSTOM;
    open->bcond_tan = CUSTOM;
    open->bcond_tan2d = CUSTOM;
    open->bcond_Vz = NOGRAD;
    open->bcond_Kz = NOGRAD;
    /* Should set INFACE for L_EDGE and B_EDGE; set in preprocess.c  */

    if (strlen(bdrypath))
      sprintf(open->tsfn, "%s%s", bdrypath, file[1]);
    else
      strcpy(open->tsfn, file[1]);
    for (i = 2; i < nfiles - 2; i++) {
      if (strlen(bdrypath))
	sprintf(open->tsfn, "%s %s%s",open->tsfn, bdrypath, file[i]);
      else
	sprintf(open->tsfn, "%s %s",open->tsfn, file[i]);
    }
    strcpy(cusu1, "uv_to_u1");
    strcpy(cusu2, "uv_to_u2");
    if (mode & NEST_STD) {
      strcpy(cusu1, "hdstd_to_u1");
      strcpy(cusu2, "hdstd_to_u2");
    }
    if (suff) {
      sprintf(open->cusname_u1, "%s %s", cusu1, file[nfiles-2]);
      sprintf(open->cusname_u2, "%s %s", cusu2, file[nfiles-1]);
    } else {
      sprintf(open->cusname_u1, "%s %s_uv_nor.mpk", cusu1, file[nfiles-2]);
      sprintf(open->cusname_u2, "%s %s_uv_tan.mpk", cusu2, file[nfiles-1]);
    }
    if (do_baro) {
      strcpy(cusu1, "uvav_to_u1av");
      strcpy(cusu2, "uvav_to_u2av");
      if (mode & NEST_STD) {
	strcpy(cusu1, "hdstd_to_u1av");
	strcpy(cusu2, "hdstd_to_u2av");
      }
      if (suff) {
	sprintf(open->cusname_u1, "%s %s", cusu1, file[nfiles-2]);
	sprintf(open->cusname_u2, "%s %s", cusu2, file[nfiles-1]);
      } else {
	sprintf(open->cusname_u1av, "%s %s_uv_nor.mpk", cusu1, file[nfiles-2]);
	sprintf(open->cusname_u2av, "%s %s_uv_tan.mpk", cusu2, file[nfiles-1]);
      }
    } else {
      open->bcond_nor2d = VERTIN;
      open->bcond_tan2d = VERTIN;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Flather radiation                                               */
  if (mode & NEST_FLA) {
    int velf = 0;
    open->sbcond = NEST_FLA|NEST_BARO;

    /* Sea level                                                     */
    open->bcond_ele = FLATHE|FILEIN|GRAVTY;
    if (velf) open->bcond_ele = FLATHR|FILEIN;

    /* Temperature and salinity                                      */
    for (i = open->atr; i < open->ntr; i++) {
      if (strcmp(tracers[i].name, "temp") == 0) {
	tid = i;
	if (open->bcond_tra[i] == NOGRAD)
	  open->bcond_tra[i] = FILEIN|UPSTRM;
      }
      if (strcmp(tracers[i].name, "salt") == 0) {
	sid = i;
	if (open->bcond_tra[i] == NOGRAD)
	  open->bcond_tra[i] = FILEIN|UPSTRM;
      }
    }

    /* Velocity                                                      */
    open->bcond_nor = NOGRAD;
    if (velf) open->bcond_nor = CUSTOM;
    open->bcond_nor2d = FLATHR|CUSTOM;
    open->bcond_tan = GRAVTY;
    open->bcond_tan2d = GRAVTY;
    if (velf) {
      open->bcond_tan = CUSTOM;
      open->bcond_tan2d = VERTIN;
    }
    open->stagger = INFACE;
    open->bcond_Vz = NOGRAD;
    open->bcond_Kz = NOGRAD;

    if (strlen(bdrypath))
      sprintf(open->tsfn, "%s%s", bdrypath, file[1]);
    else
      strcpy(open->tsfn, file[1]);
    for (i = 2; i < nfiles - 2; i++) {
      if (strlen(bdrypath))
	sprintf(open->tsfn, "%s %s%s",open->tsfn, bdrypath, file[i]);
      else
	sprintf(open->tsfn, "%s %s",open->tsfn, file[i]);
    }
    strcpy(cusu1, "uvav_to_u1av");
    if (mode & NEST_STD)
      strcpy(cusu1, "hdstd_to_u1av");
    sprintf(open->cusname_u1av, "%s %s", cusu1, file[nfiles-2]);
    if (velf) {
      strcpy(cusu1, "uv_to_u1");
      sprintf(open->cusname_u1, "%s %s", cusu1, file[nfiles-2]);
      strcpy(cusu2, "uv_to_u2");
      sprintf(open->cusname_u2, "%s %s", cusu2, file[nfiles-1]);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Sommerfeld radiation                                            */
  if (mode & NEST_RAD) {
    open->sbcond = NEST_RAD|NEST_BARO;

    /* Sea level                                                     */
    open->bcond_ele = FILEIN|RAYMND;

    /* Temperature and salinity                                      */
    for (i = open->atr; i < open->ntr; i++) {
      if (strcmp(tracers[i].name, "temp") == 0) {
	tid = i;
	open->bcond_tra[i] = UPSTRM;
      }
      if (strcmp(tracers[i].name, "salt") == 0) {
	sid = i;
	open->bcond_tra[i] = UPSTRM;
      }
    }

    /* Velocity                                                      */
    open->bcond_nor = NOGRAD;
    open->bcond_nor2d = NOGRAD;
    open->bcond_tan = RAYMND;
    open->bcond_tan2d = CLAMPD;
    open->bcond_Vz = NOGRAD;
    open->bcond_Kz = NOGRAD;
    open->relax_timei = 600.0;
    open->relax_time = 864000.0;

    if (strlen(bdrypath))
      sprintf(open->tsfn, "%s%s", bdrypath, file[1]);
    else
      strcpy(open->tsfn, file[1]);
    for (i = 2; i < nfiles; i++) {
      if (strlen(bdrypath))
	sprintf(open->tsfn, "%s %s%s",open->tsfn, bdrypath, file[i]);
      else
	sprintf(open->tsfn, "%s %s",open->tsfn, file[i]);
    }
  }

  /*-----------------------------------------------------------------*/
  /* River forcing                                                   */
  if (mode & RIVER) {
    open->sbcond = RIVER;
    open->bcond_ele = NOTHIN;
    open->bcond_nor = CUSTOM;
    open->bcond_nor2d = VERTIN;
    open->bcond_tan = CLAMPD;
    open->bcond_tan2d = CLAMPD;
    open->bcond_Vz = NOGRAD;
    open->bcond_Kz = NOGRAD;
    open->inverse_barometer = 0;
    for (i = open->atr; i < open->ntr; i++) {
      if (strcmp(tracers[i].name, "temp") == 0) {
	tid = i;
	if (open->bcond_tra[i] == NOGRAD)
	  open->bcond_tra[i] = TRCONC|FILEIN;
      }
      if (strcmp(tracers[i].name, "salt") == 0) {
	sid = i;
	if (open->bcond_tra[i] == NOGRAD) {
	  open->bcond_tra[i] = TRCONC|CUSTOM;
	  if (!strlen(open->cusname_t[i]))
	    sprintf(open->cusname_t[i],"0.0");
	}
      }
    }
    strcpy(open->tsfn, file[2]);
    for (i = 3; i < nfiles; i++)
      sprintf(open->tsfn, "%s %s",open->tsfn, file[i]);      
    if (strlen(bdrypath)) {
      sprintf(buf, "%s%s", bdrypath, open->tsfn);
      strcpy(open->tsfn, buf);
    }

    sprintf(open->cusname_u1, "u1flowbdry");
    sprintf(open->bflow, "%s", file[1]);
    /*open->bhc = NOTVALID;*/
  }
}

/* END std_bdry()                                                    */
/*-------------------------------------------------------------------*/

