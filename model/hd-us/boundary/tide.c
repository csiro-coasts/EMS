/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/boundary/tide.c
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
 *  $Id: tide.c 6308 2019-09-13 04:29:41Z her127 $
 *
 */

#include <string.h>
#include <stdio.h>
#include <time.h>
#include "hd.h"

#define MAXPOSITIONS 1001
#define TRUE  1
#define FALSE 0
#define MAXCONSTIT 30
#define UVSIZE 3119040
#define PCON 115                /* Maximum no. of constituents used
                                   (including Z0) */

typedef struct {
  double latmin;
  double latmax;
  double lonmin;
  double lonmax;
  double *uv;
  int nx;
  int ny;
  short pseudo;
  short radiation;
  short idiff;
} ortho_cons_t;

int lpeqmt(double *ts, double *dlat, double *tlp);
double zoneshift(char *name, double oldphase, int oldzone, int newzone);
char *strlwr(char *in);
void tptide_eot_init(ortho_cons_t *oc);
int tptide_eot(double *dlat, double *dlon, short *isdata, double *u,
               double *v, ortho_cons_t *oc);
int csrtide_eot(double *dlat, double *dlon, short *isdata, double *u,
                double *v, ortho_cons_t *oc);
int admit2(double *dlat, double *dlon, short *pseudo, short *radiation,
           double *u, double *v, double *x, double *y, double *h1,
           double *h2, double *amp, double *pha);
void tide_bdry_init(master_t *master, geometry_t *window,
                    open_bdrys_t *open);
void tide_vec_init(master_t *master, tidal_consts_t *tc, double *cellx,
		   double *celly, int *vec, int nvc);
int time_to_yr(double t, char *u);
void set_tmaps_3d(geometry_t *window, int bcond, int cs, int vs, int *map, 
                  int is, int ie, int *ivec);
int check_tidefile_type(char *name, char *fname, int mode);
double mjd2gmst(double mjd);
double mjd2gmst_d(double mjd);
void lunar_angles(double mjd, double *A, double *dec);
void moon_pos(double mjd, double *L, double *B);
double frac(double a);
double r2r(double x);
double sgn(double a);
double atan3 (double a, double b);
void moonvars(double jdate, double *rasc, double *decl, double *plon, double *plat, double *dist);

double c_b13 = -999.0;
double c_b129 = 360.0;
double c_b153 = 24.0;
int c_n1 = -1;

/* Hardwired variables                                               */
int yrs = 62;                   /* Number of years for nodal corrections */
double z0 = 0.0;                /* Mean sea level above datum */
int ms = 1;                     /* Start constituent number in i_1[] (0,1
                                   or 2) */
int me = 2;                     /* End constituent number in i_2[] (0,1 or 
                                   2) */
int i_1[3] = { 1, 9, 20 };      /* Start choice of constituent range */
int i_2[3] = { 8, 19, 30 };     /* End choice of constituent range */
int tmode = 1;                  /* tmode=0 : calculate tide from
                                   orthoweights */
         /* tmode=1 : calculate tide from harmonics */

/*-------------------------------------------------------------------*/
/* Routine to initialise custom tide forcing along open boundaries.  */
/* The constituent list is read from the parameter file, then phase  */
/* and amplitude information are read from file and populate the     */
/* tide data structures.                                             */
/*-------------------------------------------------------------------*/
void custom_tide_init(master_t *master, 
		      geometry_t **window,
		      int mode              /* TD_ETA = elevation    */
		                            /* TD_U = u velocity     */
		                            /* TD_V = v velocity     */
		      )
{
  geometry_t *geom = master->geom;
  open_bdrys_t *open;
  cstring *con_name;
  tidal_consts_t *tc = NULL;
  FILE *fp2;
  int i, k, t, c, cc, m, n, nn;
  int nc, ncmax;
  int size;
  int *mask;
  char buf[MAXSTRLEN], name[5];
  int have_tide_files = 0;
  timeseries_t *tsa, *tsp;
  int tsp_id, tsa_id;
  int fid, ncerr;
  int y1;
  double stime = master->t;
  double tzi = 0.0, tz = tm_tz_offset(master->timeunit);
  double d1, d2 = 0.0;
  int bcond = 0;
  char units[MAXSTRLEN], ampunits[MAXSTRLEN], type[MAXSTRLEN];
  int vs;
  int bs, be, *bec;
  int botzf, fw;
  double *xloc, *yloc;
  char files[MAXNUMTSFILES][MAXSTRLEN];

  datafile_t *df;
  df_variable_t *v;
  df_attribute_t *a;

  char *tname[MAXCONSTIT + 1] = {
    /* long period */
    "LP   55.565",              /* 1 */
    "Sa   56.554",              /* 2 */
    "Ssa  57.555",              /* 3 */
    "TERa 58.554",              /* 4 */
    "Mm   65.455",              /* 5 */
    "Mf   75.555",              /* 6 */
    "TERm 85.455",              /* 7 */
    "93a  93.555",              /* 8 */
    /* diurnal */
    "2Q1 125.755",              /* 9 */
    "Q1  135.655",              /* 10 */
    "O1  145.555",              /* 11 */
    "M1  155.655",              /* 12 */
    "P1  163.555",              /* 13 */
    "S1  164.556",              /* 14 */
    "K1  165.555",              /* 15 */
    "PHI1 167.555",             /* 16 */
    "J1  175.455",              /* 17 */
    "OO1 185.555",              /* 18 */
    "NU1 195.455",              /* 19 */
    /* semi-diurnal */
    "227 227.655",              /* 20 */
    "2N2 235.755",              /* 21 */
    "MU2 237.555",              /* 22 */
    "N2  245.655",              /* 23 */
    "NU2 247.455",              /* 24 */
    "M2  255.555",              /* 25 */
    "L2  265.455",              /* 26 */
    "T2  272.556",              /* 27 */
    "S2  273.555",              /* 28 */
    "K2  275.555",              /* 29 */
    "285 285.455"
  };                            /* 30 */

  /* Check to see if the CSR tides have been initialised.            */
  if(strlen(master->nodal_dir) > 0)
    have_tide_files = 1;

  /* Get the timezone                                                */
  if (strlen(master->tide_con_file)) {
    tsp = (timeseries_t *)malloc(sizeof(timeseries_t));
    memset(tsp, 0, sizeof(timeseries_t));
    ts_read(master->tide_con_file, tsp);
    tzi = tm_tz_offset(tsp->t_units);
    ts_free(tsp);
    free(tsp);
    tsa = hd_ts_read(master, master->tide_con_file, 0);
  }

  /*-----------------------------------------------------------------*/
  /* Parse windows and boundaries                                    */
  for (n = 1; n <= geom->nwindows; n++) {
    ncmax = 0;
    fw = 0;
    for (nn = 0; nn < window[n]->nobc; nn++) {
      open = window[n]->open[nn];
      vs = 0;

      /* Set the pointers                                            */
      if (mode == TD_NU) {
	tc = &open->tun;
	tc->type = TD_NU;
	strcpy(type, "normal u");
      } else if (mode == TD_NV) {
	tc = &open->tvn;
	tc->type = TD_NV;
	strcpy(type, "normal v");
      } else if (mode == TD_TU) {
	tc = &open->tut;
	tc->type = TD_TU;
	strcpy(type, "tangential u");
      } else if (mode == TD_TV) {
	tc = &open->tvt;
	tc->type = TD_TV;
	strcpy(type, "tangential v");
      } else {
	tc = &open->tc;
	tc->type = TD_ETA;
	strcpy(type, "elevation");
      }

      /* Get the cells to prescribe tidal harmonics for              */
      if (mode & (TD_NU|TD_NV)) {
	if (open->bcond_nor & TIDALC || open->bcond_nor2d & TIDALC) {
	  vs = open->no2_e1;
	  bs = 1;
	  be = open->no2_e1;
	  bec = open->obc_e1;
	  size = open->no3_e1;
	  bcond = open->bcond_nor;
	}
	strcpy(units, "ms-1");
	xloc = window[n]->u1x;
	yloc = window[n]->u1y;
      } else if (mode & (TD_TU|TD_TV)) {
	if (open->bcond_tan & TIDALC || open->bcond_tan2d & TIDALC) {
	  vs = open->to2_e1 - open->no3_e1;
	  bs = open->no3_e1 + 1;
	  be = open->to2_e1;
	  bec = open->obc_e1;
	  size = open->to3_e1;
	  bcond = open->bcond_tan;
	}
	strcpy(units, "ms-1");
	xloc = window[n]->u1x;
	yloc = window[n]->u1y;
      } else {
	if (open->bcond_ele & TIDALC) {
	  vs = open->no2_t;
	  bs = 1;
	  be = open->no2_t;
	  bec = open->obc_t;
	  size = open->no2_t;
	}
	strcpy(units, "m");
	xloc = window[n]->cellx;
	yloc = window[n]->celly;
      }
      strcpy(buf, open->tide_con);
      ncmax = parseline(buf, (char **)files, MAXNUMTSFILES);

      if (vs) {

	/*-----------------------------------------------------------*/
	/* Allocate the cells for the tidal structures               */
	tc->amp = d_alloc_2d(ncmax + 1, vs + 1);
	tc->pha = d_alloc_2d(ncmax + 1, vs + 1);
	tc->map = i_alloc_1d(size + 1);
	memset(tc->map, 0, (size + 1) * sizeof(int));

	/*-----------------------------------------------------------*/
	/* Set the harmonics if TIDALC boundaries are found          */
	if (!have_tide_files)
	  hd_quit
	    ("custom_tide_init: Custom tide specified, but no constituent files.");
	
	/* Read in the desired constituent names                     */
	strcpy(buf, open->tide_con);
	nc = parseline(buf, (char **)files, MAXNUMTSFILES);

	con_name = (cstring *)malloc(sizeof(cstring) * nc);
	for (t = 0; t < nc; t++) {
	  strcpy(con_name[t], ((char **)files)[t]);
	}
	mask = i_alloc_1d(nc+1);
	for (t = 0; t <= nc; t++) mask[t] = -1;

	/* Get the number of constituents                            */
	tc->nt = 0;
	for (t = 0; t < nc; t++) {
	  for (i = i_1[0]; i <= i_2[2]; i++) {
	    size = strlen(tname[i-1]) - strlen(strpbrk(tname[i-1], " "));
	    strncpy(name, tname[i-1], size);
	    name[size] = '\0';
	    if (strcmp(con_name[t], name) == 0) {
	      tc->nt++;
	      mask[tc->nt] = i-1;
	      break;
	    }
	  }
	  if (mask[tc->nt] == -1)
	    hd_warn("custom_tide_init: Can't recognize tidal constituent %s. Ignoring this constituent.\n", con_name[t]);
	}
	if (!tc->nt && !fw) {
	  hd_warn("custom_tide_init: No tide constitudents specified for %s OBC%d %s\n",type, nn, open->name);
	  fw = 1;
	}

	/* Allocate memory */
	tc->z0 = z0;
	tc->yr = d_alloc_1d(yrs + 1);
	size = tc->nt + 1;
	tc->sigma = d_alloc_1d(size);
	tc->j = d_alloc_2d(size, yrs + 1);
	tc->v = d_alloc_2d(size, yrs + 1);
	tc->vpv = d_alloc_2d(size, yrs + 1);

	/* Get the tidal constituent names */
	tc->tname = (char **)malloc((tc->nt + 1) * sizeof(char *));
	for (i = 1; i <= tc->nt; i++) {
	  tc->tname[i] = (char *)malloc(MAXCONSTIT * sizeof(char));
	  strcpy(tc->tname[i], tname[mask[i]]);
	}
	i_free_1d(mask);
	free(con_name);

	/* Read in astronomical arguments */
	for (i = 1; i <= tc->nt; i++) {
	  size = strlen(tc->tname[i]) - strlen(strpbrk(tc->tname[i], " "));
	  strncpy(name, tc->tname[i], size);
	  name[size] = '\0';
	  sprintf(buf, "%s%s.con", master->nodal_dir, strlwr(name));
	  if (strcasecmp(name, "z0") != 0) {
	    if ((fp2 = fopen(buf, "r")) == NULL) {
	      hd_quit("custom_bdry_init: Can't open file '%s'.", buf);
	      continue;
	    }
	    fscanf(fp2, "%lf", &d1);
	    tc->sigma[i] = d1 / 3600.0; /* Convert to degrees per second */
	    k = 1;
	    while (k <= yrs) {
	      if ((y1 = feof(fp2)))
		break;
	      fscanf(fp2, "%lf %lf %lf %lf", &tc->yr[k], &tc->j[k][i],
		     &tc->v[k][i], &tc->vpv[k][i]);
	      k++;
	    }
	    fclose(fp2);
	    if (k - 1 != yrs)
	      hd_quit("custom_bdry_init: Could not read %d years.", yrs);
	  }
	}

	/* Check that the model timeunit is in seconds */
	sprintf(buf, "%s", strlwr(master->timeunit));
	if (strcasecmp(strtok(buf, " "), "seconds") == 0)
	  d2 = 86400.0;
	else if (strcasecmp(strtok(buf, " "), "minutes") == 0)
	  d2 = 1440.0;
	else if (strcasecmp(strtok(buf, " "), "minutes") == 0)
	  d2 = 1440.0;
	else if (strcasecmp(strtok(buf, " "), "hours") == 0)
	  d2 = 24.0;
	else if (strcasecmp(strtok(buf, " "), "days") == 0)
	  d2 = 1.0;
	else
	  hd_quit("%s is an invalid timeunit for tidal forcing\n",
		  master->timeunit);

	/* Get the Julian days of valid years */
	sprintf(buf, "days since 1990-01-01 00:00:00");
	m = (int)(tz / 3600);
	if (m >= 0 && m < 10)
	  sprintf(buf, "%s +0%d", buf, m);
	else if (m >= 10)
	  sprintf(buf, "%s +%d", buf, m);
	else if (m < 0 && m > -10)
	  sprintf(buf, "%s -0%d", buf, -m);
	else if (m <= -10)
	  sprintf(buf, "%s -%d", buf, -m);
	tc->yr[1] = -7305.0;          /* 1 Jan 1970 relative to 1 Jan 1990 */
	tm_change_time_units(buf, master->timeunit, &tc->yr[1], 1);
	y1 = 1970;
	for (k = 2; k <= yrs; k++) {
	  if (y1 % 4 == 0 && (y1 % 100 != 0 || y1 % 400 == 0))
	    d1 = 366.0;
	  else
	    d1 = 365.0;
	  y1++;
	  /* Nodal correction years converted to seconds relative to the */
	  /* model timeunit.  */
	  tc->yr[k] = tc->yr[k - 1] + d1 * d2;
	}

	/* Find the first year less than Julian day */
	k = 1;
	while (stime - tz > tc->yr[k] && k <= yrs)
	  k++;
	tc->k = k - 1;
	if (k == yrs) {
	  hd_warn("custom_tide_init: Simulation period lies outside nodel correction range.\n");
	}

	/* Read the amplitude and phase files and check they contain */
	/* the specified constituents. */
	for (i = 1; i <= tc->nt; i++) {

	  size = strlen(tc->tname[i]) - strlen(strpbrk(tc->tname[i], " "));
	  strncpy(name, tc->tname[i], size);
	  name[size] = '\0';
	  if (mode &  (TD_NU|TD_TU)) {
	    sprintf(buf,"%s_amp",name);
	    botzf = check_tidefile_type(buf, master->tide_con_file, mode);
	    if (botzf == -1) 
	      hd_quit("Can't find units ms-1 or m2s-1 for %s in %s\n", buf, master->tide_con_file);
	  } else if (mode &  (TD_NV|TD_TV)) {
	    sprintf(buf,"%s_amp",name);
	    botzf = check_tidefile_type(buf, master->tide_con_file, mode);
	    if (botzf == -1) 
	      hd_quit("Can't find units ms-1 or m2s-1 for %s in %s\n", buf, master->tide_con_file);
	  } else {
	   sprintf(buf,"%s_amp",name);
	   botzf = 0;
	  }
	  if (botzf) {
	    master->tidef |= TD_TRAN;
	  } else {
	    master->tidef |= TD_VEL;
	  }

	  /*
	  tsa = frcw_read_cell_ts(master, master->tide_con_file, stime,
				  buf, units, &d2, &tsa_id, NULL);
	  */
	  
	  if ((tsa_id = ts_get_index(tsa, buf)) < 0)
	    hd_quit("Can't find variable %s in file %s\n", buf, master->tide_con_file);
	  
	  if (mode &  (TD_NU|TD_TU)) {
	    if (botzf)
	      sprintf(buf,"%s_phase_U",name);
	    else
	      sprintf(buf,"%s_phase_u",name);
	  } else if (mode &  (TD_NV|TD_TV)) {
	    if (botzf)
	      sprintf(buf,"%s_phase_V",name);
	    else
	      sprintf(buf,"%s_phase_v",name);
	  } else
	    sprintf(buf,"%s_phase",name);	  

	  /*
	  tsp = frc_read_cell_ts(master, master->tide_con_file, stime,
				 buf, "degrees", &d2, &tsp_id, NULL);
	  */
	  if ((tsp_id = ts_get_index(tsa, buf)) < 0)
	    hd_quit("Can't find variable %s in file %s\n", buf, master->tide_con_file);
	  
	  /* Get the tidal harmonics for each cell on the boundary */
	  vs = 1;
	  for (cc = bs; cc <= be; cc++) {
	    c = bec[cc];
	    tc->amp[vs][i] = ts_eval_xy(tsa, tsa_id, stime, xloc[c], yloc[c]);
	    tc->pha[vs][i] = ts_eval_xy(tsa, tsp_id, stime, xloc[c], yloc[c]);
	    tc->pha[vs][i] = zoneshift(name, tc->pha[vs][i], (int)tzi/3600, 0); 
	    tc->map[cc] = vs;
	    set_tmaps_3d(window[n], bcond, c, vs, tc->map, bs, be, bec);
	    /* Velocoty amplitudes from file are transports; divide  */
	    /* by model depth to get velocity.                       */
	    if (botzf) {
	      int c1 = window[n]->e2c[c][0];
	      int c2 = window[n]->e2c[c][1];
	      /*double depth = 0.5 * fabs(window[n]->botz[c1] + window[n]->botz[c2])*/
	      double depth = fabs(window[n]->botzu1[c]);
	      tc->amp[vs][i] /= depth;
	    }
	    vs++;
	  }
	  
	  /*hd_ts_free(master, tsa);*/
	  /*hd_ts_free(master, tsp);*/
	  
	}
	if (DEBUG("init_m")) {
	  cc = 1;
	  dlog("init_m", "\n\n%s (window %d) has %d tidal constituents\n",
	       open->name, n, open->tc.nt);
	  c = open->obc_t[cc];
	  dlog("init_m", "Values at %5.2f %5.2f (%d %d)\n",
	       window[n]->cellx[c], window[n]->celly[c], 
	       geom->s2i[window[n]->wsa[c]], geom->s2j[window[n]->wsa[c]]);
	  /* Freq. in deg/sec = 24*(360/freq)/86400 hours */
	  dlog("init_m", " Name Doodson   Freq(o/s) Amp(m)  Phase(o)\n");
	  for(i = 1; i <= open->tc.nt; i++) {
	    dlog("init_m", "  %s %f %f %f\n",
		 open->tc.tname[i], open->tc.sigma[i],
		 open->tc.amp[cc][i], open->tc.pha[cc][i]);
	  }
	}
      }
    }
  }
  hd_ts_free(master, tsa);
  /*hd_ts_free(master, tsp);*/
}

/* END custom_tide_init()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to check if a tide file contains transports (with units   */
/* m2s-1) for the amplitute of velocity constituents, or speed (with */
/* units ms-1). Returns 1 for the former and 0 for the latter.       */
/*-------------------------------------------------------------------*/
int check_tidefile_type(char *name, char *fname, int mode)
{
  char buf[MAXSTRLEN];
  char buf1[MAXSTRLEN];
  char buf2[MAXSTRLEN];
  char units[MAXSTRLEN];
  char units1[MAXSTRLEN];
  int fid;

  nc_open(fname, NC_NOWRITE, &fid);

  /* Check for velocity amplitudes in the file */
  if (mode & (TD_NU|TD_TU))
    sprintf(buf, "%s_u",name);
  else if (mode &  (TD_NV|TD_TV))
    sprintf(buf, "%s_v",name);
  else
    sprintf(buf, "%s", name);
  sprintf(units, "%c", '\0');
  if (nc_get_att_text(fid, ncw_var_id(fid, buf), "units", units) >= 0) {
    strcpy(buf1, units);
    if (strcmp(buf1, "ms-1") == 0) {
      strcpy(name, buf);
      nc_close(fid);
      return(0);
    }
  }
  if (mode &  (TD_NU|TD_TU))
    sprintf(buf, "%s_U",name);
  else if (mode &  (TD_NV|TD_TV))
    sprintf(buf, "%s_V",name);
  else
    sprintf(buf, "%s",name);
  sprintf(buf1, "%c", '\0');
  nc_get_att_text(fid, ncw_var_id(fid, buf), "units", units1);
  strcpy(buf2, units1);
  if (strcmp(buf2, "m2s-1") == 0) {
    strcpy(name, buf);
    nc_close(fid);
    return(1);
  }
  nc_close(fid);
  return(-1);
}

/* END check_tidefile_type()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to initialise custom tide forcing within the mesh.        */
/* The constituent list is read from the parameter file, then phase  */
/* and amplitude information are read from file and populate the     */
/* tide data structures.                                             */
/*-------------------------------------------------------------------*/
void custom_tide_grid_init(master_t *master, 
			   geometry_t **window
			   )
{
  geometry_t *geom = master->geom;
  open_bdrys_t *open;
  cstring *con_name;
  tidal_consts_t *tc = NULL;
  FILE *fp2;
  int i, k, t, c, cc, m, n, nn;
  int nc, ncmax;
  int size;
  int *mask;
  char buf[MAXSTRLEN], name[5];
  int have_tide_files = 0;
  timeseries_t *tsa, *tsp;
  int tsp_id, tsa_id;
  int y1;
  double stime = master->t;
  double tzi = 0.0, tz = tm_tz_offset(master->timeunit);
  double d1, d2 = 0.0;
  char files[MAXNUMTSFILES][MAXSTRLEN];
  char tide_con[MAXSTRLEN];

  char *tname[MAXCONSTIT + 1] = {
    /* long period */
    "LP   55.565",              /* 1 */
    "Sa   56.554",              /* 2 */
    "Ssa  57.555",              /* 3 */
    "TERa 58.554",              /* 4 */
    "Mm   65.455",              /* 5 */
    "Mf   75.555",              /* 6 */
    "TERm 85.455",              /* 7 */
    "93a  93.555",              /* 8 */
    /* diurnal */
    "2Q1 125.755",              /* 9 */
    "Q1  135.655",              /* 10 */
    "O1  145.555",              /* 11 */
    "M1  155.655",              /* 12 */
    "P1  163.555",              /* 13 */
    "S1  164.556",              /* 14 */
    "K1  165.555",              /* 15 */
    "PHI1 167.555",             /* 16 */
    "J1  175.455",              /* 17 */
    "OO1 185.555",              /* 18 */
    "NU1 195.455",              /* 19 */
    /* semi-diurnal */
    "227 227.655",              /* 20 */
    "2N2 235.755",              /* 21 */
    "MU2 237.555",              /* 22 */
    "N2  245.655",              /* 23 */
    "NU2 247.455",              /* 24 */
    "M2  255.555",              /* 25 */
    "L2  265.455",              /* 26 */
    "T2  272.556",              /* 27 */
    "S2  273.555",              /* 28 */
    "K2  275.555",              /* 29 */
    "285 285.455"
  };                            /* 30 */

  strcpy(tide_con, "M2 S2 N2 K2 K1 O1 P1 Q1");

  /* Check to see if the CSR tides have been initialised.            */
  if(strlen(master->nodal_dir) > 0)
    have_tide_files = 1;

  /* Get the timezone                                                */
  if (strlen(master->tide_con_file)) {
    tsp = (timeseries_t *)malloc(sizeof(timeseries_t));
    memset(tsp, 0, sizeof(timeseries_t));
    ts_read(master->tide_con_file, tsp);
    tzi = tm_tz_offset(tsp->t_units);
    ts_free(tsp);
    free(tsp);
  } else {
    hd_warn("Can't find TIDE_CONSTITUENTS file for custom tide computation.\n");
    master->numbers1 &= ~TPXO;
    for (n = 1; n <= geom->nwindows; n++) {
      win_priv_t *wincon = window[n]->wincon;
      wincon->numbers1 &= ~TPXO;
    }
    return;
  }

  /*-----------------------------------------------------------------*/
  /* Parse windows and boundaries                                    */
  for (n = 1; n <= geom->nwindows; n++) {
    win_priv_t *wincon = window[n]->wincon;

    tc = &wincon->tc;
    ncmax = 0;

    if (wincon->csr_tide == 1)
      hd_quit("Cannot allocate tidal structure: CSR TIDAL_REMOVAL already set.\n");

    tc->map = NULL;
    size = window[n]->szcS;

    strcpy(buf, tide_con);
    ncmax = parseline(buf, (char **)files, MAXNUMTSFILES);

    if (size) {

      /*-------------------------------------------------------------*/
      /* Allocate the cells for the tidal structures                 */
      tc->amp = d_alloc_2d(ncmax + 1, size);
      tc->pha = d_alloc_2d(ncmax + 1, size);
      
      /*-------------------------------------------------------------*/
      /* Set the harmonics if TIDALC boundaries are found            */
      if (!have_tide_files)
	hd_quit
	  ("custom_tide_grid_init: Custom tide specified, but no nodal constituent files.");
	
      /* Read in the desired constituent names                       */
      strcpy(buf, tide_con);
      nc = parseline(buf, (char **)files, MAXNUMTSFILES);

      con_name = (cstring *)malloc(sizeof(cstring) * nc);
      for (t = 0; t < nc; t++) {
	strcpy(con_name[t], ((char **)files)[t]);
      }
      mask = i_alloc_1d(nc+1);
      for (t = 0; t <= nc; t++) mask[t] = -1;

      /* Get the number of constituents                              */
      tc->nt = 0;
      for (t = 0; t < nc; t++) {
	for (i = i_1[0]; i <= i_2[2]; i++) {
	  size = strlen(tname[i-1]) - strlen(strpbrk(tname[i-1], " "));
	  strncpy(name, tname[i-1], size);
	  name[size] = '\0';
	  if (strcmp(con_name[t], name) == 0) {
	    tc->nt++;
	    mask[tc->nt] = i-1;
	    break;
	  }
	}
	if (mask[tc->nt] == -1)
	  hd_warn("custom_tide_grid_init: Can't recognize tidal constituent %s. Ignoring this constituent.\n", con_name[t]);
      }
      
      /* Allocate memory                                             */
      tc->z0 = z0;
      tc->yr = d_alloc_1d(yrs + 1);
      size = tc->nt + 1;
      tc->sigma = d_alloc_1d(size);
      tc->j = d_alloc_2d(size, yrs + 1);
      tc->v = d_alloc_2d(size, yrs + 1);
      tc->vpv = d_alloc_2d(size, yrs + 1);

      /* Get the tidal constituent names                             */
      tc->tname = (char **)malloc((tc->nt + 1) * sizeof(char *));
      for (i = 1; i <= tc->nt; i++) {
	tc->tname[i] = (char *)malloc(MAXCONSTIT * sizeof(char));
	strcpy(tc->tname[i], tname[mask[i]]);
      }
      i_free_1d(mask);
      free(con_name);
      
      /* Read in astronomical arguments                              */
      for (i = 1; i <= tc->nt; i++) {
	size = strlen(tc->tname[i]) - strlen(strpbrk(tc->tname[i], " "));
	strncpy(name, tc->tname[i], size);
	name[size] = '\0';
	sprintf(buf, "%s%s.con", master->nodal_dir, strlwr(name));
	if (strcasecmp(name, "z0") != 0) {
	  if ((fp2 = fopen(buf, "r")) == NULL) {
	    hd_quit("custom_bdry_init: Can't open file '%s'.", buf);
	    continue;
	  }
	  fscanf(fp2, "%lf", &d1);
	  tc->sigma[i] = d1 / 3600.0; /* Convert to degrees per second */
	  k = 1;
	  while (k <= yrs) {
	    if ((y1 = feof(fp2)))
	      break;
	    fscanf(fp2, "%lf %lf %lf %lf", &tc->yr[k], &tc->j[k][i],
		   &tc->v[k][i], &tc->vpv[k][i]);
	    k++;
	  }
	  fclose(fp2);
	  if (k - 1 != yrs)
	    hd_quit("custom_bdry_init: Could not read %d years.", yrs);
	}
      }

      /* Check that the model timeunit is in seconds                 */
      sprintf(buf, "%s", strlwr(master->timeunit));
      if (strcasecmp(strtok(buf, " "), "seconds") == 0)
	d2 = 86400.0;
      else if (strcasecmp(strtok(buf, " "), "minutes") == 0)
	d2 = 1440.0;
      else if (strcasecmp(strtok(buf, " "), "minutes") == 0)
	d2 = 1440.0;
      else if (strcasecmp(strtok(buf, " "), "hours") == 0)
	d2 = 24.0;
      else if (strcasecmp(strtok(buf, " "), "days") == 0)
	d2 = 1.0;
      else
	hd_quit("%s is an invalid timeunit for tidal forcing\n",
		master->timeunit);

      /* Get the Julian days of valid years                          */
      sprintf(buf, "days since 1990-01-01 00:00:00");
      m = (int)(tz / 3600);
      if (m >= 0 && m < 10)
	sprintf(buf, "%s +0%d", buf, m);
      else if (m >= 10)
	sprintf(buf, "%s +%d", buf, m);
      else if (m < 0 && m > -10)
	sprintf(buf, "%s -0%d", buf, -m);
      else if (m <= -10)
	sprintf(buf, "%s -%d", buf, -m);
      tc->yr[1] = -7305.0;          /* 1 Jan 1970 relative to 1 Jan 1990 */
      tm_change_time_units(buf, master->timeunit, &tc->yr[1], 1);
      y1 = 1970;
      for (k = 2; k <= yrs; k++) {
	if (y1 % 4 == 0 && (y1 % 100 != 0 || y1 % 400 == 0))
	  d1 = 366.0;
	else
	  d1 = 365.0;
	y1++;
	/* Nodal correction years converted to seconds relative to   */
	/* the model timeunit.                                       */
	tc->yr[k] = tc->yr[k - 1] + d1 * d2;
      }

      /* Find the first year less than Julian day                    */
      k = 1;
      while (stime - tz > tc->yr[k] && k <= yrs)
	k++;
      tc->k = k - 1;
      if (k == yrs) {
	hd_warn("custom_tide_init: Simulation period lies outside nodel correction range.\n");
      }
      
      /* Read the amplitude and phase files and check they contain   */
      /* the specified constituents.                                 */
      for (i = 1; i <= tc->nt; i++) {
	
	size = strlen(tc->tname[i]) - strlen(strpbrk(tc->tname[i], " "));
	strncpy(name, tc->tname[i], size);
	name[size] = '\0';
	sprintf(buf,"%s_amp",name);
	tsa = frc_read_cell_ts(master, master->tide_con_file, stime,
			       buf, "m", &d2, &tsa_id, NULL);
	sprintf(buf,"%s_phase",name);	  
	tsp = frc_read_cell_ts(master, master->tide_con_file, stime,
			       buf, "degrees", &d2, &tsp_id, NULL);
	
	/* Get the tidal harmonics for each cell in the domain       */
	for (cc = 1; cc <= window[n]->b2_t; cc++) {
	  c = window[n]->w2_t[cc];
	  tc->amp[c][i] = ts_eval_xy(tsa, tsa_id, stime, window[n]->cellx[c], window[n]->celly[c]);
	  tc->pha[c][i] = ts_eval_xy(tsp, tsp_id, stime, window[n]->cellx[c], window[n]->celly[c]);
	  tc->pha[c][i] = zoneshift(name, tc->pha[c][i], (int)tzi/3600, 0);
	}
	/* Open boundary ghosts
	for (nn = 0; nn < window[n]->nobc; nn++) {
	  open_bdrys_t *open = window[n]->open[nn];
	  for (cc = 1; cc <= open->no2_e1; cc++) {
	    c = open->ogc_t[cc];
	    tc->amp[c][i] = ts_eval_xy(tsa, tsa_id, stime, window[n]->cellx[c], window[n]->celly[c]);
	    tc->pha[c][i] = ts_eval_xy(tsp, tsp_id, stime, window[n]->cellx[c], window[n]->celly[c]);
	    tc->pha[c][i] = zoneshift(name, tc->pha[c][i], (int)tzi/3600, 0);
	  }
	}
	*/
	hd_ts_free(master, tsa);
	hd_ts_free(master, tsp);
      }
      if (DEBUG("init_m")) {
	cc = 1;
	dlog("init_m", "\n\nwindow %d) has %d tidal constituents\n", n, tc->nt);	     
	/* Freq. in deg/sec = 24*(360/freq)/86400 hours              */
	dlog("init_m", " Name Doodson   Freq(o/s) Amp(m)  Phase(o)\n");
	for(i = 1; i <= tc->nt; i++) {
	  dlog("init_m", "  %s %f %f %f\n",tc->tname[i], tc->sigma[i],
	       tc->amp[cc][i], tc->pha[cc][i]);
	}
      }
    }
  }

  hd_ts_free(master, tsa);
  hd_ts_free(master, tsp);
}

/* END custom_tide_grid_init()                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to initialise custom tide forcing within the mesh.        */
/* The constituent list is read from the parameter file, then phase  */
/* and amplitude information are read from file and populate the     */
/* tide data structures.                                             */
/*-------------------------------------------------------------------*/
void custom_tide_grid_init_uv(master_t *master, 
			      geometry_t **window,
			      int mode
			      )
{
  geometry_t *geom = master->geom;
  cstring *con_name;
  tidal_consts_t *tc = NULL;
  FILE *fp2;
  int i, k, t, c, cc, m, n, nn;
  int nc, ncmax;
  int size;
  int *mask;
  char buf[MAXSTRLEN], name[5];
  int have_tide_files = 0;
  int botzf;
  timeseries_t *tsa;
  int tsp_id, tsa_id;
  int y1;
  double stime = master->t;
  double tzi = 0.0, tz = tm_tz_offset(master->timeunit);
  double d1, d2 = 0.0;
  char files[MAXNUMTSFILES][MAXSTRLEN];
  char tide_con[MAXSTRLEN];

  char *tname[MAXCONSTIT + 1] = {
    /* long period */
    "LP   55.565",              /* 1 */
    "Sa   56.554",              /* 2 */
    "Ssa  57.555",              /* 3 */
    "TERa 58.554",              /* 4 */
    "Mm   65.455",              /* 5 */
    "Mf   75.555",              /* 6 */
    "TERm 85.455",              /* 7 */
    "93a  93.555",              /* 8 */
    /* diurnal */
    "2Q1 125.755",              /* 9 */
    "Q1  135.655",              /* 10 */
    "O1  145.555",              /* 11 */
    "M1  155.655",              /* 12 */
    "P1  163.555",              /* 13 */
    "S1  164.556",              /* 14 */
    "K1  165.555",              /* 15 */
    "PHI1 167.555",             /* 16 */
    "J1  175.455",              /* 17 */
    "OO1 185.555",              /* 18 */
    "NU1 195.455",              /* 19 */
    /* semi-diurnal */
    "227 227.655",              /* 20 */
    "2N2 235.755",              /* 21 */
    "MU2 237.555",              /* 22 */
    "N2  245.655",              /* 23 */
    "NU2 247.455",              /* 24 */
    "M2  255.555",              /* 25 */
    "L2  265.455",              /* 26 */
    "T2  272.556",              /* 27 */
    "S2  273.555",              /* 28 */
    "K2  275.555",              /* 29 */
    "285 285.455"
  };                            /* 30 */

  strcpy(tide_con, "M2 S2 N2 K2 K1 O1 P1 Q1");

  /* Check to see if the CSR tides have been initialised.            */
  if(strlen(master->nodal_dir) > 0)
    have_tide_files = 1;

  /* Get the timezone                                                */
  if (strlen(master->tide_con_file)) {
    tsa = hd_ts_read(master, master->tide_con_file, 0);
    tzi = tm_tz_offset(tsa->t_units);
  } else {
    hd_warn("Can't find TIDE_CONSTITUENTS file for custom tide computation.\n");
    return;
  }

  /*-----------------------------------------------------------------*/
  /* Parse windows and boundaries                                    */
  for (n = 1; n <= geom->nwindows; n++) {
    win_priv_t *wincon = window[n]->wincon;

    if (mode & TD_NU)
      tc = &wincon->tcu;
    else if (mode & TD_NV)
      tc = &wincon->tcv;
    ncmax = 0;

    if (wincon->csr_tide == 1)
      hd_quit("Cannot allocate tidal structure: CSR TIDAL_REMOVAL already set.\n");

    tc->map = NULL;
    size = window[n]->szeS;

    strcpy(buf, tide_con);
    ncmax = parseline(buf, (char **)files, MAXNUMTSFILES);

    if (size) {

      /*-------------------------------------------------------------*/
      /* Allocate the cells for the tidal structures                 */
      tc->amp = d_alloc_2d(ncmax + 1, size);
      tc->pha = d_alloc_2d(ncmax + 1, size);
      
      /*-------------------------------------------------------------*/
      /* Set the harmonics if TIDALC boundaries are found            */
      if (!have_tide_files)
	hd_quit
	  ("custom_tide_grid_init: Custom tide specified, but no nodal constituent files.");
	
      /* Read in the desired constituent names                       */
      strcpy(buf, tide_con);
      nc = parseline(buf, (char **)files, MAXNUMTSFILES);

      con_name = (cstring *)malloc(sizeof(cstring) * nc);
      for (t = 0; t < nc; t++) {
	strcpy(con_name[t], ((char **)files)[t]);
      }
      mask = i_alloc_1d(nc+1);
      for (t = 0; t <= nc; t++) mask[t] = -1;

      /* Get the number of constituents                              */
      tc->nt = 0;
      for (t = 0; t < nc; t++) {
	for (i = i_1[0]; i <= i_2[2]; i++) {
	  size = strlen(tname[i-1]) - strlen(strpbrk(tname[i-1], " "));
	  strncpy(name, tname[i-1], size);
	  name[size] = '\0';
	  if (strcmp(con_name[t], name) == 0) {
	    tc->nt++;
	    mask[tc->nt] = i-1;
	    break;
	  }
	}
	if (mask[tc->nt] == -1)
	  hd_warn("custom_tide_grid_init: Can't recognize tidal constituent %s. Ignoring this constituent.\n", con_name[t]);
      }
      
      /* Allocate memory                                             */
      tc->z0 = z0;
      tc->yr = d_alloc_1d(yrs + 1);
      size = tc->nt + 1;
      tc->sigma = d_alloc_1d(size);
      tc->j = d_alloc_2d(size, yrs + 1);
      tc->v = d_alloc_2d(size, yrs + 1);
      tc->vpv = d_alloc_2d(size, yrs + 1);

      /* Get the tidal constituent names                             */
      tc->tname = (char **)malloc((tc->nt + 1) * sizeof(char *));
      for (i = 1; i <= tc->nt; i++) {
	tc->tname[i] = (char *)malloc(MAXCONSTIT * sizeof(char));
	strcpy(tc->tname[i], tname[mask[i]]);
      }
      i_free_1d(mask);
      free(con_name);

      /* Read in astronomical arguments                              */
      for (i = 1; i <= tc->nt; i++) {
	size = strlen(tc->tname[i]) - strlen(strpbrk(tc->tname[i], " "));
	strncpy(name, tc->tname[i], size);
	name[size] = '\0';
	sprintf(buf, "%s%s.con", master->nodal_dir, strlwr(name));
	if (strcasecmp(name, "z0") != 0) {
	  if ((fp2 = fopen(buf, "r")) == NULL) {
	    hd_quit("custom_bdry_init: Can't open file '%s'.", buf);
	    continue;
	  }
	  fscanf(fp2, "%lf", &d1);
	  tc->sigma[i] = d1 / 3600.0; /* Convert to degrees per second */
	  k = 1;
	  while (k <= yrs) {
	    if ((y1 = feof(fp2)))
	      break;
	    fscanf(fp2, "%lf %lf %lf %lf", &tc->yr[k], &tc->j[k][i],
		   &tc->v[k][i], &tc->vpv[k][i]);
	    k++;
	  }
	  fclose(fp2);
	  if (k - 1 != yrs)
	    hd_quit("custom_bdry_init: Could not read %d years.", yrs);
	}
      }

      /* Check that the model timeunit is in seconds                 */
      sprintf(buf, "%s", strlwr(master->timeunit));
      if (strcasecmp(strtok(buf, " "), "seconds") == 0)
	d2 = 86400.0;
      else if (strcasecmp(strtok(buf, " "), "minutes") == 0)
	d2 = 1440.0;
      else if (strcasecmp(strtok(buf, " "), "minutes") == 0)
	d2 = 1440.0;
      else if (strcasecmp(strtok(buf, " "), "hours") == 0)
	d2 = 24.0;
      else if (strcasecmp(strtok(buf, " "), "days") == 0)
	d2 = 1.0;
      else
	hd_quit("%s is an invalid timeunit for tidal forcing\n",
		master->timeunit);

      /* Get the Julian days of valid years                          */
      sprintf(buf, "days since 1990-01-01 00:00:00");
      m = (int)(tz / 3600);
      if (m >= 0 && m < 10)
	sprintf(buf, "%s +0%d", buf, m);
      else if (m >= 10)
	sprintf(buf, "%s +%d", buf, m);
      else if (m < 0 && m > -10)
	sprintf(buf, "%s -0%d", buf, -m);
      else if (m <= -10)
	sprintf(buf, "%s -%d", buf, -m);
      tc->yr[1] = -7305.0;          /* 1 Jan 1970 relative to 1 Jan 1990 */
      tm_change_time_units(buf, master->timeunit, &tc->yr[1], 1);
      y1 = 1970;
      for (k = 2; k <= yrs; k++) {
	if (y1 % 4 == 0 && (y1 % 100 != 0 || y1 % 400 == 0))
	  d1 = 366.0;
	else
	  d1 = 365.0;
	y1++;
	/* Nodal correction years converted to seconds relative to   */
	/* the model timeunit.                                       */
	tc->yr[k] = tc->yr[k - 1] + d1 * d2;
      }

      /* Find the first year less than Julian day                    */
      k = 1;
      while (stime - tz > tc->yr[k] && k <= yrs)
	k++;
      tc->k = k - 1;
      if (k == yrs) {
	hd_warn("custom_tide_init: Simulation period lies outside nodel correction range.\n");
      }
      
      /* Read the amplitude and phase files and check they contain   */
      /* the specified constituents.                                 */
      for (i = 1; i <= tc->nt; i++) {
	
	size = strlen(tc->tname[i]) - strlen(strpbrk(tc->tname[i], " "));
	strncpy(name, tc->tname[i], size);
	name[size] = '\0';

	sprintf(buf,"%s_amp",name);
	botzf = check_tidefile_type(buf, master->tide_con_file, mode);
	/*
	tsa = frcw_read_cell_ts(master, master->tide_con_file, stime,
				buf, "ms-1", &d2, &tsa_id, NULL);
	*/
	if (botzf) {
	  master->tidef |= TD_TRAN;
	} else {
	  master->tidef |= TD_VEL;
	}

	if ((tsa_id = ts_get_index(tsa, buf)) < 0)
	  hd_quit("Can't find variable %s in file %s\n", buf, master->tide_con_file);

	if (mode & TD_NU) {
	  if (botzf)
	    sprintf(buf,"%s_phase_U",name);
	  else
	    sprintf(buf,"%s_phase_u",name);
	} else {
	  if (botzf)
	    sprintf(buf,"%s_phase_V",name);
	  else
	    sprintf(buf,"%s_phase_v",name);
	}
	if ((tsp_id = ts_get_index(tsa, buf)) < 0)
	  hd_quit("Can't find variable %s in file %s\n", buf, master->tide_con_file);
	/*
	tsp = frc_read_cell_ts(master, master->tide_con_file, stime,
			       buf, "degrees", &d2, &tsp_id, NULL);
	*/
	/* Get the tidal harmonics for each cell in the domain       */
	for (cc = 1; cc <= window[n]->b2_e1; cc++) {
	  c = window[n]->w2_e1[cc];
	  tc->amp[c][i] = ts_eval_xy(tsa, tsa_id, stime, window[n]->u1x[c], window[n]->u1y[c]);
	  tc->pha[c][i] = ts_eval_xy(tsa, tsp_id, stime, window[n]->u1x[c], window[n]->u1y[c]);
	  tc->pha[c][i] = zoneshift(name, tc->pha[c][i], (int)tzi/3600, 0);
	  /* Velocoty amplitudes from file are transports; divide    */
	  /* by model depth to get velocity.                         */
	  if (botzf) {
	    int c1 = window[n]->e2c[c][0];
	    int c2 = window[n]->e2c[c][1];
	    /*double depth = 0.5 * fabs(window[n]->botz[c1] + window[n]->botz[c2]);*/
	    double depth = fabs(window[n]->botzu1[c]);
	    tc->amp[c][i] /= depth;
	  }
	}
	/*
	hd_ts_free(master, tsa);
	hd_ts_free(master, tsp);
	*/
      }
    }
  }
  hd_ts_free(master, tsa);
}

/* END custom_tide_grid_init_uv()                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets the tidal mapping function for 3D boundaries                 */
/*-------------------------------------------------------------------*/
void set_tmaps_3d(geometry_t *window, 
                  int bcond, 
                  int cs, 
                  int vs, 
                  int *map, 
                  int is, 
                  int ie, 
                  int *ivec
                  )
{
  int c, cc;

  if (bcond & TIDALC) {
    for (cc = is; cc <= ie; cc++) {
      c = ivec[cc];
      if (window->m2d[c] == cs)
        map[cc] = vs;
    }
  }
}

/* END set_tmaps_3d()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to initialise tidal forcing. This can only be done after  */
/* the window is initialised with cellx and celly data, and the cell */
/* centered locations for the open boundaries within windows have    */
/* been defined.                                                     */
/*-------------------------------------------------------------------*/
void csr_tide_init(master_t *master, geometry_t **window)
{
  int n, nn;
  geometry_t *geom = master->geom;
  open_bdrys_t *open;
  int have_tide_files = 0;

  /* Check to see if the CSR tides have been initialised. */
  if( strlen(master->orthoweights) > 0 && strlen(master->nodal_dir) > 0)
    have_tide_files = 1;

  /* Parse boundaries */
  for (n = 1; n <= geom->nwindows; n++) {
    for (nn = 0; nn < window[n]->nobc; nn++) {
      open = window[n]->open[nn];
      if (open->bcond_ele & TIDALH || open->bcond_nor2d == (FLATHR|TIDALH)) {

        if (!have_tide_files)
          hd_quit
            ("csr_tide_init: Tidal boundary nodes specified, but no constituent files.");

        /*tide_bdry_init(master, window[n], open);*/
	tide_vec_init(master, &open->tc, window[n]->cellx, window[n]->celly, 
		      open->obc_t, open->no2_t);

	if (DEBUG("init_m")) {
	  int i, cc = 1, c;
	  dlog("init_m", "\n\n%s (window %d) has %d tidal constituents\n",
	       open->name, n, open->tc.nt);
	  c = open->obc_t[cc];
	  dlog("init_m", "Values at %5.2f %5.2f (%d %d)\n",
	       window[n]->cellx[c], window[n]->celly[c], 
	       geom->s2i[window[n]->wsa[c]], geom->s2j[window[n]->wsa[c]]);
	  /* Freq. in deg/sec = 24*(360/freq)/86400 hours */
	  dlog("init_m", " Name Doodson   Freq(o/s) Amp(m)  Phase(o)\n");
	  for(i = 1; i <= open->tc.nt; i++) {
	    dlog("init_m", "  %s %f %f %f\n",
		 open->tc.tname[i], open->tc.sigma[i],
		 open->tc.amp[cc][i], open->tc.pha[cc][i]);
	  }
	}
      }
    }
  }
}

/* END csr_tide_init()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to initialise tidal forcing. This can only be done after  */
/* the window is initialised with cellx and celly data, and the cell */
/* centered locations for the open boundaries within windows have    */
/* been defined.                                                     */
/*-------------------------------------------------------------------*/
void csr_tide_grid_init(master_t *master, geometry_t **window)
{
  int n, do_tide;

  /* Check to see if csr tidal constituents are required for ACTIVE  */
  /* alert checking, using eta relaxation.                           */
  if (master->alertf & ACTIVE && (master->div2d_f || master->eta_f) &&
      master->etarlx & ALERT)
    do_tide = 1;

  /* Check to see if the csr tides have been initialised. If not     */
  /* boundaries are not forced with csr tides (TIDALH) and removing  */
  /* an estimated tide from the sea level is not required.           */
  if( strlen(master->orthoweights) > 0 && strlen(master->nodal_dir) > 0)
    do_tide = 1;
  else
    do_tide = 0;

  /* Loop through the windows and initialise the csr tide constants  */
  /* if required.                                                    */
  for (n = 1; n <= master->nwindows; n++) {
    win_priv_t *wincon = window[n]->wincon;    /* Window constants   */
    if (do_tide && wincon->tide_r & CSR_R) {
      tide_vec_init(master, &wincon->tc, window[n]->cellx, window[n]->celly, 
		    window[n]->w2_t, window[n]->b2_t);
      wincon->csr_tide = 1;
    } else {
      wincon->csr_tide = 0;
    }
  }
}

/* END csr_tide_grid_init()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to create and initialise the tidal harmonics data         */
/* structure for each open boundary.                                 */
/*-------------------------------------------------------------------*/
void tide_bdry_init(master_t *master, geometry_t *window, /* Window
                                                             geometry */
                    open_bdrys_t *open  /* Open boundary data structure */
  )
{
  tidal_consts_t *tc = NULL;
  ortho_cons_t *oc = NULL;
  FILE *fp2 = NULL;
  short isdata;
  int i, k, m, cc, c;
  int size;
  int y1;
  double stime = master->t;
  double tz = tm_tz_offset(master->timeunit);
  double d1, d2 = 0.0;

  double u[6], v[6];
  double *uv = NULL;
  double x[MAXPOSITIONS], y[MAXPOSITIONS];
  double h1[MAXPOSITIONS], h2[MAXPOSITIONS];
  double amp[MAXPOSITIONS], pha[MAXPOSITIONS];
  double dlat, dlon;
  char buf[MAXSTRLEN], name[5];
  map_proj_t *mp = NULL;

  char *tname[MAXCONSTIT + 1] = {
    /* long period */
    "LP   55.565",              /* 1 */
    "Sa   56.554",              /* 2 */
    "Ssa  57.555",              /* 3 */
    "TERa 58.554",              /* 4 */
    "Mm   65.455",              /* 5 */
    "Mf   75.555",              /* 6 */
    "TERm 85.455",              /* 7 */
    "93a  93.555",              /* 8 */
    /* diurnal */
    "2Q1 125.755",              /* 9 */
    "Q1  135.655",              /* 10 */
    "O1  145.555",              /* 11 */
    "M1  155.655",              /* 12 */
    "P1  163.555",              /* 13 */
    "S1  164.556",              /* 14 */
    "K1  165.555",              /* 15 */
    "PHI1 167.555",             /* 16 */
    "J1  175.455",              /* 17 */
    "OO1 185.555",              /* 18 */
    "NU1 195.455",              /* 19 */
    /* semi-diurnal */
    "227 227.655",              /* 20 */
    "2N2 235.755",              /* 21 */
    "MU2 237.555",              /* 22 */
    "N2  245.655",              /* 23 */
    "NU2 247.455",              /* 24 */
    "M2  255.555",              /* 25 */
    "L2  265.455",              /* 26 */
    "T2  272.556",              /* 27 */
    "S2  273.555",              /* 28 */
    "K2  275.555",              /* 29 */
    "285 285.455"
  };                            /* 30 */
/*  mask[i]=1 : tidal constituent included; mask[i]=0 : constituent omitted 
                0 1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 20 1 2 3 4 5 6 7 8 9 0*/
  /* int mask[31]={0,0,0,0,0,0,0,0,0,1,1, 1,1,1,1,1,1,1,1,0,0,
     1,1,1,1,1,1,1,1,1,0}; */
  int mask[31] =
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 0
  };

  /*-----------------------------------------------------------------*/
  /* Open and initialize the tide data structure */
  /* Allocate memory for the tide structure */
  tc = &open->tc;

  oc = (ortho_cons_t *)malloc(sizeof(ortho_cons_t));
  memset(oc, 0, sizeof(ortho_cons_t));
  oc->uv = (double *)malloc(UVSIZE * sizeof(double));

  /* Get the number of constituents and names */
  tc->nt = 0;
  for (m = ms; m <= me; ++m)
    for (i = i_1[m]; i <= i_2[m]; i++)
      if (mask[i])
        tc->nt++;

  /* Allocate memory */
  tc->z0 = z0;
  tc->yr = d_alloc_1d(yrs + 1);
  size = tc->nt + 1;
  tc->sigma = d_alloc_1d(size);
  tc->j = d_alloc_2d(size, yrs + 1);
  tc->v = d_alloc_2d(size, yrs + 1);
  tc->vpv = d_alloc_2d(size, yrs + 1);
  tc->amp = d_alloc_2d(size, open->no2_t + 1);
  tc->pha = d_alloc_2d(size, open->no2_t + 1);

  /* Get the tidal constituent names */
  tc->tname = (char **)malloc((tc->nt + 1) * sizeof(char *));
  for (i = 1; i <= tc->nt; i++) {
    tc->tname[i] = (char *)malloc(MAXCONSTIT * sizeof(char));
  }
  k = 1;
  for (m = ms; m <= me; ++m)
    for (i = i_1[m]; i <= i_2[m]; i++) {
      if (mask[i]) {
        strcpy(tc->tname[k], tname[i - 1]);
        k++;
      }
    }

  /* Read in astronomical arguments */
  for (i = 1; i <= tc->nt; i++) {
    size = strlen(tc->tname[i]) - strlen(strpbrk(tc->tname[i], " "));
    strncpy(name, tc->tname[i], size);
    name[size] = '\0';
    sprintf(buf, "%s%s.con", master->nodal_dir, strlwr(name));
    if (strcasecmp(name, "z0") != 0) {
      if ((fp2 = fopen(buf, "r")) == NULL) {
        hd_quit("tide_bdry_init: Can't open file '%s'.", buf);
        continue;
      }
      fscanf(fp2, "%lf", &d1);
      tc->sigma[i] = d1 / 3600.0; /* Convert to degrees per second */
      k = 1;
      while (k <= yrs) {
        if ((y1 = feof(fp2)))
          break;
        fscanf(fp2, "%lf %lf %lf %lf", &tc->yr[k], &tc->j[k][i],
               &tc->v[k][i], &tc->vpv[k][i]);
        k++;
      }
      fclose(fp2);
      if (k - 1 != yrs)
        hd_quit("tide_bdry_init: Could not read %d years.", yrs);
    }
  }

  /* Check that the model timeunit is in seconds */
  sprintf(buf, "%s", strlwr(master->timeunit));
  if (strcasecmp(strtok(buf, " "), "seconds") == 0)
    d2 = 86400.0;
  else if (strcasecmp(strtok(buf, " "), "minutes") == 0)
    d2 = 1440.0;
  else if (strcasecmp(strtok(buf, " "), "minutes") == 0)
    d2 = 1440.0;
  else if (strcasecmp(strtok(buf, " "), "hours") == 0)
    d2 = 24.0;
  else if (strcasecmp(strtok(buf, " "), "days") == 0)
    d2 = 1.0;
  else
    hd_quit("%s is an invalid timeunit for tidal forcing\n",
            master->timeunit);

  /* Get the Julian days of valid years */
  sprintf(buf, "days since 1990-01-01 00:00:00");
  m = (int)(tz / 3600);
  if (m >= 0 && m < 10)
    sprintf(buf, "%s +0%d", buf, m);
  else if (m >= 10)
    sprintf(buf, "%s +%d", buf, m);
  else if (m < 0 && m > -10)
    sprintf(buf, "%s -0%d", buf, -m);
  else if (m <= -10)
    sprintf(buf, "%s -%d", buf, -m);
  tc->yr[1] = -7305.0;          /* 1 Jan 1970 relative to 1 Jan 1990 */
  tm_change_time_units(buf, master->timeunit, &tc->yr[1], 1);
  y1 = 1970;
  for (k = 2; k <= yrs; k++) {
    if (y1 % 4 == 0 && (y1 % 100 != 0 || y1 % 400 == 0))
      d1 = 366.0;
    else
      d1 = 365.0;
    y1++;
    /* Nodal correction years converted to seconds relative to the */
    /* model timeunit.  */
    tc->yr[k] = tc->yr[k - 1] + d1 * d2;
  }

  /* Find the first year less than Julian day */
  k = 1;
  while (stime - tz > tc->yr[k] && k <= yrs)
    k++;
  tc->k = k - 1;
  if (k == yrs) hd_warn("tide_bdry_init: Simulation period lies outside nodel correction range.\n");

  /*-----------------------------------------------------------------*/
  /* Get the tidal harmonics from the orthoweights.  */

  /* Initialise the orthoweight structure.  */
  tptide_eot_init(oc);

  if (strlen(projection) > 0) {
    if (strncasecmp(projection, GEOGRAPHIC_TAG, strlen(GEOGRAPHIC_TAG)) != 0) {
       char *args[256];
       int nargs = parseline(strdup(projection), args, 256);
       mp = mp_init(nargs, args);
    }
  }

  /* Get the tidal harmonics for each cell on the boundary */
  for (cc = 1; cc <= open->no2_t; cc++) {
    c = open->obc_t[cc];

    dlat = window->celly[c];
    dlon = window->cellx[c];

    if (mp != NULL)
       mp_inverse(mp, dlon, dlat, &dlat, &dlon);

    /* Get the orthoweights */
    csrtide_eot(&dlat, &dlon, &isdata, u, v, oc);

    /* Get the amplitude and phase */
    if (isdata) {
      admit2(&dlat, &dlon, &oc->pseudo, &oc->radiation, u, v, x, y, h1, h2,
             amp, pha);

      /* Copy the valid constituents to the data structure */
      k = 1;
      for (m = ms; m <= me; ++m)
        for (i = i_1[m]; i <= i_2[m]; i++) {
          if (mask[i]) {
            tc->amp[cc][k] = amp[i - 1] / 1e3;
            tc->pha[cc][k] = pha[i - 1];
            k++;
          }
        }
    } else {
      hd_warn("Invalid harmonics for (%f,%f)\n", dlat, dlon);
    }
  }
  free((void *)uv);

  if (mp != NULL)
     mp_cleanup(mp);
}

/* END tide_bdry_init()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to create and initialise the tidal harmonics data         */
/* structure for each cell in vector vec[].                          */
/*-------------------------------------------------------------------*/
void tide_vec_init(master_t *master,    /* Master data               */
		   tidal_consts_t *tc,  /* Tidal constants           */
		   double *cellx,       /* e1 cell center locations  */
		   double *celly,       /* e2 cell center locations  */
		   int *vec,            /* Array of tidal cells      */
		   int nvc              /* Size of vec               */
  )
{
  ortho_cons_t *oc = NULL;
  FILE *fp2 = NULL;
  short isdata;
  int i, k, m, cc, c;
  int size;
  int y1;
  double stime = master->t;
  double tz = tm_tz_offset(master->timeunit);
  double d1, d2 = 0.0;

  double u[6], v[6];
  double *uv = NULL;
  double x[MAXPOSITIONS], y[MAXPOSITIONS];
  double h1[MAXPOSITIONS], h2[MAXPOSITIONS];
  double amp[MAXPOSITIONS], pha[MAXPOSITIONS];
  double dlat, dlon;
  char buf[MAXSTRLEN], name[5];
  map_proj_t *mp = NULL;

  char *tname[MAXCONSTIT + 1] = {
    /* long period */
    "LP   55.565",              /* 1 */
    "Sa   56.554",              /* 2 */
    "Ssa  57.555",              /* 3 */
    "TERa 58.554",              /* 4 */
    "Mm   65.455",              /* 5 */
    "Mf   75.555",              /* 6 */
    "TERm 85.455",              /* 7 */
    "93a  93.555",              /* 8 */
    /* diurnal */
    "2Q1 125.755",              /* 9 */
    "Q1  135.655",              /* 10 */
    "O1  145.555",              /* 11 */
    "M1  155.655",              /* 12 */
    "P1  163.555",              /* 13 */
    "S1  164.556",              /* 14 */
    "K1  165.555",              /* 15 */
    "PHI1 167.555",             /* 16 */
    "J1  175.455",              /* 17 */
    "OO1 185.555",              /* 18 */
    "NU1 195.455",              /* 19 */
    /* semi-diurnal */
    "227 227.655",              /* 20 */
    "2N2 235.755",              /* 21 */
    "MU2 237.555",              /* 22 */
    "N2  245.655",              /* 23 */
    "NU2 247.455",              /* 24 */
    "M2  255.555",              /* 25 */
    "L2  265.455",              /* 26 */
    "T2  272.556",              /* 27 */
    "S2  273.555",              /* 28 */
    "K2  275.555",              /* 29 */
    "285 285.455"
  };                            /* 30 */
/*  mask[i]=1 : tidal constituent included; mask[i]=0 : constituent omitted 
                0 1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 20 1 2 3 4 5 6 7 8 9 0*/
  /* int mask[31]={0,0,0,0,0,0,0,0,0,1,1, 1,1,1,1,1,1,1,1,0,0,
     1,1,1,1,1,1,1,1,1,0}; */
  int mask[31] =
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 0
  };

  /*-----------------------------------------------------------------*/
  /* Open and initialize the tide data structure */
  oc = (ortho_cons_t *)malloc(sizeof(ortho_cons_t));
  memset(oc, 0, sizeof(ortho_cons_t));
  oc->uv = (double *)malloc(UVSIZE * sizeof(double));

  /* Get the number of constituents and names */
  tc->nt = 0;
  for (m = ms; m <= me; ++m)
    for (i = i_1[m]; i <= i_2[m]; i++)
      if (mask[i])
        tc->nt++;

  /* Allocate memory */
  tc->z0 = z0;
  tc->yr = d_alloc_1d(yrs + 1);
  size = tc->nt + 1;
  tc->sigma = d_alloc_1d(size);
  tc->j = d_alloc_2d(size, yrs + 1);
  tc->v = d_alloc_2d(size, yrs + 1);
  tc->vpv = d_alloc_2d(size, yrs + 1);
  tc->amp = d_alloc_2d(size, nvc + 1);
  tc->pha = d_alloc_2d(size, nvc + 1);

  /* Get the tidal constituent names */
  tc->tname = (char **)malloc((tc->nt + 1) * sizeof(char *));
  for (i = 1; i <= tc->nt; i++) {
    tc->tname[i] = (char *)malloc(MAXCONSTIT * sizeof(char));
  }
  k = 1;
  for (m = ms; m <= me; ++m)
    for (i = i_1[m]; i <= i_2[m]; i++) {
      if (mask[i]) {
        strcpy(tc->tname[k], tname[i - 1]);
        k++;
      }
    }

  /* Read in astronomical arguments */
  for (i = 1; i <= tc->nt; i++) {
    size = strlen(tc->tname[i]) - strlen(strpbrk(tc->tname[i], " "));
    strncpy(name, tc->tname[i], size);
    name[size] = '\0';
    sprintf(buf, "%s%s.con", master->nodal_dir, strlwr(name));
    if (strcasecmp(name, "z0") != 0) {
      if ((fp2 = fopen(buf, "r")) == NULL) {
        hd_quit("tide_vec_init: Can't open file '%s'.", buf);
        continue;
      }
      fscanf(fp2, "%lf", &d1);
      tc->sigma[i] = d1 / 3600.0; /* Convert to degrees per second */
      k = 1;
      while (k <= yrs) {
        if ((y1 = feof(fp2)))
          break;
        fscanf(fp2, "%lf %lf %lf %lf", &tc->yr[k], &tc->j[k][i],
               &tc->v[k][i], &tc->vpv[k][i]);
        k++;
      }
      fclose(fp2);
      if (k - 1 != yrs)
        hd_quit("tide_vec_init: Could not read %d years.", yrs);
    }
  }

  /* Check that the model timeunit is in seconds */
  sprintf(buf, "%s", strlwr(master->timeunit));
  if (strcasecmp(strtok(buf, " "), "seconds") == 0)
    d2 = 86400.0;
  else if (strcasecmp(strtok(buf, " "), "minutes") == 0)
    d2 = 1440.0;
  else if (strcasecmp(strtok(buf, " "), "minutes") == 0)
    d2 = 1440.0;
  else if (strcasecmp(strtok(buf, " "), "hours") == 0)
    d2 = 24.0;
  else if (strcasecmp(strtok(buf, " "), "days") == 0)
    d2 = 1.0;
  else
    hd_quit("%s is an invalid timeunit for tidal forcing\n",
            master->timeunit);

  /* Get the Julian days of valid years */
  sprintf(buf, "days since 1990-01-01 00:00:00");
  m = (int)(tz / 3600);
  if (m >= 0 && m < 10)
    sprintf(buf, "%s +0%d", buf, m);
  else if (m >= 10)
    sprintf(buf, "%s +%d", buf, m);
  else if (m < 0 && m > -10)
    sprintf(buf, "%s -0%d", buf, -m);
  else if (m <= -10)
    sprintf(buf, "%s -%d", buf, -m);
  tc->yr[1] = -7305.0;          /* 1 Jan 1970 relative to 1 Jan 1990 */
  tm_change_time_units(buf, master->timeunit, &tc->yr[1], 1);
  y1 = 1970;
  for (k = 2; k <= yrs; k++) {
    if (y1 % 4 == 0 && (y1 % 100 != 0 || y1 % 400 == 0))
      d1 = 366.0;
    else
      d1 = 365.0;
    y1++;
    /* Nodal correction years converted to seconds relative to the */
    /* model timeunit.  */
    tc->yr[k] = tc->yr[k - 1] + d1 * d2;
  }

  /* Find the first year less than Julian day */
  k = 1;
  while (stime - tz > tc->yr[k] && k <= yrs)
    k++;
  tc->k = k - 1;
  if (k == yrs) hd_warn("tide_vec_init: Simulation period lies outside nodel correction range.\n");

  /*-----------------------------------------------------------------*/
  /* Get the tidal harmonics from the orthoweights.  */

  /* Initialise the orthoweight structure.  */
  tptide_eot_init(oc);

  if (strlen(projection) > 0) {
    if (strncasecmp(projection, GEOGRAPHIC_TAG, strlen(GEOGRAPHIC_TAG)) != 0) {
       char *args[256];
       int nargs = parseline(strdup(projection), args, 256);
       mp = mp_init(nargs, args);
    }
  }

  /* Get the tidal harmonics for each cell on the boundary */
  for (cc = 1; cc <= nvc; cc++) {
    c = vec[cc];

    dlat = celly[c];
    dlon = cellx[c];

    if (mp != NULL)
       mp_inverse(mp, dlon, dlat, &dlat, &dlon);

    /* Get the orthoweights */
    csrtide_eot(&dlat, &dlon, &isdata, u, v, oc);

    /* Get the amplitude and phase */
    if (isdata) {
      admit2(&dlat, &dlon, &oc->pseudo, &oc->radiation, u, v, x, y, h1, h2,
             amp, pha);

      /* Copy the valid constituents to the data structure */
      k = 1;
      for (m = ms; m <= me; ++m)
        for (i = i_1[m]; i <= i_2[m]; i++) {
          if (mask[i]) {
            tc->amp[cc][k] = amp[i - 1] / 1e3;
            tc->pha[cc][k] = pha[i - 1];
            k++;
          }
        }
    } else {
      hd_warn("Invalid harmonics for (%f,%f)\n", dlat, dlon);
    }
  }
  free((void *)uv);

  if (mp != NULL)
     mp_cleanup(mp);
}

/* END tide_vec_init()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates the tide from tidal harmonics.                         */
/* eta = z0 + sum(fn.Hn.cos(s.t+un-gn)                               */
/* z0 = mean sea level above datum                                   */
/* gn = phase of constituent of equilibrium tide at Greenwhich       */
/* s = speed of constituent (deg/s)                                  */
/* Hn = amplitude of constituent (m)                                 */
/* fn and un = nodal corrections relative to Greenwhich              */
/* Note : time input to this function should be in seconds.          */
/*-------------------------------------------------------------------*/
double csr_tide_eval(tidal_consts_t *tc,  /* Tidal constants         */
                     int cc,              /* Boundary index          */
                     double jd  /* Julian time at Greenwhich (s)     */
  )
{
  int i, k1, k2, ci;
  double j1, j2, v1, v2, w1, w2;
  double vpv, v0, nj, nv;
  double eta, time;
  double d2r = atan(1.0) * 4.0 / 180.0; /* Degrees to radians        */

  ci = (tc->map != NULL) ? tc->map[cc] : cc;

  if (jd <= tc->yr[yrs]) {
    /* Get the indicies for the years bracketing jd                  */
    k1 = tc->k;
    k2 = (k1 == yrs) ? k1 : k1 + 1;
    if (jd > tc->yr[k2]) {
      k1 = tc->k = k1 + 1;
      k2 = (k1 == yrs) ? k1 : k1 + 1;
    }
    time = jd - tc->yr[k1];
    w2 = time / (tc->yr[k2] - tc->yr[k1]);
    w1 = 1.0 - w2;
  } else {
    k1 = k2 = tc->k;
    time = jd - tc->yr[k1];
    w2 = 1.0;
    w1 = 0.0;
  }

  /* Sum the tide. Note tidal phases and nodal phase correcions are  */
  /* given relative to Greenwhich: the Julian day input must be      */
  /* corrected for the local time zone.                              */
  eta = tc->z0;
  for (i = 1; i <= tc->nt; i++) {
    j1 = tc->j[k1][i];
    j2 = tc->j[k2][i];
    v1 = tc->v[k1][i];
    v2 = tc->v[k2][i];
    vpv = tc->vpv[k1][i];
    v0 = vpv - v1;
    nj = j1 * w1 + j2 * w2;
    nv = v0 + v1 * w1 + v2 * w2;
    eta +=
      tc->amp[ci][i] * nj *
      cos((tc->sigma[i] * time + nv - tc->pha[ci][i]) * d2r);
  }
  return (eta);
}

/* END tide_eval()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* To compute diurnal and semidiurnal ocean tides from an            */
/*      "orthoweight" file and long period tides from equilibrium    */
/*      theory for a given time, latitude and longitude.             */
/* Note: An "orthoweight" file can contain geocentric tide, ocean    */
/*         tide, load tide, or tide model difference orthoweights.   */
/*         If the orthoweight file in use contains differences or    */
/*         loads then the user should beware of using the "tide"     */
/*         return variable and should use "tds" instead.             */
/*  This version of the driver presumes that the orthoweight file is */
/*     that for elastic ocean tide, or geocentric tide.  The file    */
/*     should be linked to fortran unit 51.                          */
/* Coded by Richard Eanes, University of Texas Center for Space      */
/*     Research eanes@csr.utexas.edu                                 */
/* Version 2.0 - 03 Nov 94                                           */
/*             - 13 Nov 94, added better documentation on djmd.      */
/* Version 3.0 - March 1995, new version of tptide allows compact    */
/*               integer format for orthoweight file                 */
/*                                                                   */
/* Inputs: (double precision)                                        */
/*   dmjd    Modified Julian Date in Days, number of days (including */
/*           fraction) since 17 Nov 1858 00:00.  The standard        */
/*           astronomical julian date (jd) is related to dmjd by     */
/*           jd = 2400000.5 + dmjd.  The modified julian date of     */
/*           1 Jan 1958 00:00 is 36204.  The modified julian date of */
/*           1 Jan 1994 12:00 is 49353.5.                            */
/*   dlat    North Latitude (in degrees)                             */
/*   dlon    East Longitude (in degrees)                             */
/*                                                                   */
/* Outputs (double precision)                                        */
/*   tds     Diurnal and semidiurnal tide value in millimeters       */
/*   tlp     Long period tide value in millimeters (equilibrium      */
/*           theory)                                                 */
/*   tide    Total tide value in millimeters                         */
/*   isdata  Logical, true if location is within the domain of the   */
/*           tide model, false otherwise.  When false the tide       */
/*           values are returned as zero.                            */
/* The following outputs can be used for constructing maps of the    */
/* tidal fields at the major diurnal and semidiurnal constituents    */
/* using subroutine admit2.                                          */
/*   u,v     Diurnal and semidiurnal orthoweights for input location */
/*           (can be used in subroutine admit to compute harmonic    */
/*           information).                                           */
/*   pseudo  Logical, switch passed to subroutine admit determining  */
/*           what variety of orthoweights are on the tide file.      */
/*   radiation  Logical, true if orthoweights modified for radiation */
/*           potential.                                              */
/*-------------------------------------------------------------------*/
int csrtide_eot(double *dlat, double *dlon, short *isdata, double *u,
                double *v, ortho_cons_t *oc)
{

  /* Parameter adjustments */
  --v;
  --u;

  /* lpeqmt(&smjd, dlat, tlp); */
  tptide_eot(dlat, dlon, isdata, &u[1], &v[1], oc);

  return 0;
}

/* END csrtide_eot()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read the orthoweight file and store in uv              */
/*-------------------------------------------------------------------*/
void tptide_eot_init(ortho_cons_t *oc)
{
  FILE *fp1;
  double undef = 99999.0;
  int ml, val;
  int itype;
  int irad, ifmt;
  int i, j, n, i1 = 0, idiff;
  int nx, ny;
  int iuv[14];
  double u1, v1, u2, v2, u3, v3, u4, v4, u5, v5, u6, v6;
  char title1[80], title2[80], buf[100];
  double latmin, latmax, lonmin, lonmax;
  double *uv = oc->uv;

  /*-----------------------------------------------------------------*/
  /* Open the orthoweight file */
  if((fp1 = fopen(master->orthoweights, "r")) == NULL)
    hd_quit("Can't open file %s\n", master->orthoweights);

  /*-----------------------------------------------------------------*/
  /* Read t/p -derived ortho-weights */
  for (j = 1; j <= 361; ++j) {
    for (i = 1; i <= 720; ++i) {
      uv[(((((i + j * 720) << 1) + 1) * 3 + 1) << 1) - 8660] = undef;
    }
  }
  fgets(title1, 80, fp1);
  fgets(title2, 80, fp1);
  fscanf(fp1, "%d %d %lf %lf %lf %lf %d %d %d %d", &nx, &ny,
         &latmin, &latmax, &lonmin, &lonmax, &idiff, &itype, &irad, &ifmt);
  oc->nx = nx;
  oc->ny = ny;
  oc->latmin = latmin;
  oc->latmax = latmax;
  oc->lonmin = lonmin;
  oc->lonmax = lonmax;
  oc->idiff = idiff;

  /* Orthoweight file format options */
  /* idiff = 1 : zero tide returned for points outside of domain */
  /* and isdata set to true.  */
  /* itype = 1 : cartwright&ray 1991 pseudo-orthoweights are */
  /* assumed.  */
  /* irad = 1 : cartwright&ray radiation potential modification is */
  /* used.  */
  /* ifmt = 1 : compact integer format used for orthoweights */

  if (itype == 1)
    oc->pseudo = TRUE;          /* c&r style orthoweights */
  else
    oc->pseudo = FALSE;         /* corrected orthoweights */
  if (irad == 1)
    oc->radiation = TRUE;       /* corrected for radiation potential */
  else
    oc->radiation = FALSE;      /* not corrected for radiation potential */
  if (nx > 720 || ny > 361) {
    hd_quit("tptide_eot_init: must change dimensions\n");
  }

  do {
    if (ifmt == 0) {            /* original c&r format */

      if ((i1 = fscanf(fp1, "%d%d%lf%lf%lf%lf%lf%lf",
                       &i, &j, &u1, &v1, &u2, &v2, &u3, &v3)) != 0)
        break;
      if ((i1 = fscanf(fp1, "%lf%lf%lf%lf%lf%lf",
                       &u4, &v4, &u5, &v5, &u6, &v6) != 0))
        break;

      oc->uv[(((((i + j * 720) << 1) + 1) * 3 + 1) << 1) - 8660] = u1;
      oc->uv[(((((i + j * 720) << 1) + 1) * 3 + 1) << 1) - 8659] = v1;
      oc->uv[(((((i + j * 720) << 1) + 1) * 3 + 2) << 1) - 8660] = u2;
      oc->uv[(((((i + j * 720) << 1) + 1) * 3 + 2) << 1) - 8659] = v2;
      oc->uv[(((((i + j * 720) << 1) + 1) * 3 + 3) << 1) - 8660] = u3;
      oc->uv[(((((i + j * 720) << 1) + 1) * 3 + 3) << 1) - 8659] = v3;
      oc->uv[(((((i + j * 720) << 1) + 2) * 3 + 1) << 1) - 8660] = u4;
      oc->uv[(((((i + j * 720) << 1) + 2) * 3 + 1) << 1) - 8659] = v4;
      oc->uv[(((((i + j * 720) << 1) + 2) * 3 + 2) << 1) - 8660] = u5;
      oc->uv[(((((i + j * 720) << 1) + 2) * 3 + 2) << 1) - 8659] = v5;
      oc->uv[(((((i + j * 720) << 1) + 2) * 3 + 3) << 1) - 8660] = u6;
      oc->uv[(((((i + j * 720) << 1) + 2) * 3 + 3) << 1) - 8659] = v6;

    }
    if (ifmt == 1) {            /* compact integer format */
      n = 0;
      if ((i1 = feof(fp1)))
        break;
      while (n < 14) {
        fscanf(fp1, "%s", buf);
        if (strchr(buf, '*') != NULL) {
          ml = n + atoi(strtok(buf, "*"));
          val = atoi(strtok(NULL, "*"));
          for (; n < ml; n++)
            iuv[n] = val;
        } else {
          iuv[n] = atoi(buf);
          n++;
        }
      }
      i = iuv[0];
      j = iuv[1];
      oc->uv[(((((i + j * 720) << 1) + 1) * 3 + 1) << 1) - 8660] =
        iuv[2] / 1e3;
      oc->uv[(((((i + j * 720) << 1) + 1) * 3 + 1) << 1) - 8659] =
        iuv[3] / 1e3;
      oc->uv[(((((i + j * 720) << 1) + 1) * 3 + 2) << 1) - 8660] =
        iuv[4] / 1e3;
      oc->uv[(((((i + j * 720) << 1) + 1) * 3 + 2) << 1) - 8659] =
        iuv[5] / 1e3;
      oc->uv[(((((i + j * 720) << 1) + 1) * 3 + 3) << 1) - 8660] =
        iuv[6] / 1e3;
      oc->uv[(((((i + j * 720) << 1) + 1) * 3 + 3) << 1) - 8659] =
        iuv[7] / 1e3;
      oc->uv[(((((i + j * 720) << 1) + 2) * 3 + 1) << 1) - 8660] =
        iuv[8] / 1e3;
      oc->uv[(((((i + j * 720) << 1) + 2) * 3 + 1) << 1) - 8659] =
        iuv[9] / 1e3;
      oc->uv[(((((i + j * 720) << 1) + 2) * 3 + 2) << 1) - 8660] =
        iuv[10] / 1e3;
      oc->uv[(((((i + j * 720) << 1) + 2) * 3 + 2) << 1) - 8659] =
        iuv[11] / 1e3;
      oc->uv[(((((i + j * 720) << 1) + 2) * 3 + 3) << 1) - 8660] =
        iuv[12] / 1e3;
      oc->uv[(((((i + j * 720) << 1) + 2) * 3 + 3) << 1) - 8659] =
        iuv[13] / 1e3;
    }
  } while (!i1);
  fclose(fp1);
}

/* END tptide_eot_init()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* This code is based upon the gstide program used to compute the    */
/* cartwright&ray geosat tide solution.                              */
/*  name - tptide_eot       name derivation - topex/poseidon elastic */
/*                          ocean tide */
/*  function -  to compute the ocean tidal height at a given time    */
/*              and location.                                        */
/*  arguments -                                                      */
/*     name      type  i/o               description                 */
/*     ----      ----  ---               -----------                 */
/*     dlat       d     i    north latitude (in degrees, -90 to  90) */
/*     dlon       d     i    east longitude (in degrees,   0 to 360) */
/*     time       d     i    desired time, in seconds of mjd.  e.g., */
/*                           jan 1, 1986 00:00 is 4011638400.        */
/*     tide       d     o    computed tidal height, in cm.           */
/*     isdata     l     o    logical denoting whether tide data      */
/*                           exist at desired location. If false,    */
/*                           then tide is not modified.              */
/*     u,v        d     o    the orthoweights interpolated to dlat,  */
/*                           dlon, dimensioned (3,2)                 */
/*     pseudo     l     o    true for c&r style "orthoweights"       */
/*                           this can be passed to admit.f for       */
/*                           proper interpretation of the            */
/*                           orthoweights. Set per integer on        */
/*                           orthoweight file.                       */
/*     radiation  l     o    true if the correction for radiation    */
/*                           potential is activated.                 */
/*                           set per integer on orthoweight file.    */
/*  Usage notes -                                                    */
/*  Using the user supplied orthoweight file, this routine computes  */
/*  the height of the (ocean + load) tide at the desired location.   */
/*  This quantity is appropriate for use in satellite altimetry.     */
/*  The computed tide is composed of the 30 largest spectral lines   */
/*  within both the diurnal and semidiurnal bands, sufficient for    */
/*  representing all major constituents including nodal modulations. */
/*                                                                   */
/*  Processing logic -                                               */
/*  The tidal height is computed at four grid points surrounding     */
/*  the desired location, with the final result being a bilinear     */
/*  interpolation.                                                   */
/*                                                                   */
/*  File references -                                                */
/*  The orthoweight data are read on initial call to routine:        */
/*  file iunit - global array of orthoweights defining tidal         */
/*  constants at all valid locations. (see parameter statement)      */
/*                                                                   */
/*  Important local variables -                                      */
/*  Array uv is loaded with orthoweights as follows:                 */
/*         uv(i,j,l,m,n), where                                      */
/*           i = 1 for u; 2 for v                                    */
/*           j = 1,2 or 3  (3 complex coeffs per species.)           */
/*           l = tidal species (1 or 2)                              */
/*           m = 1 to nx, longitude index number                     */
/*           n = 1 to ny, latitude index                             */
/*                                                                   */
/*  For machines like the cray, the compiler option to ignore double */
/*  precision should be selected.                                    */
/*                                                                   */
/*  Error processing - none.                                         */
/*                                                                   */
/*  Technical references -                                           */
/*  D. Cartwright and R. Ray, Oceanic tides from geosat altimetry,   */
/*        Journal of Geophysical Research, 95, 3069-3090, 1990.      */
/*        (particularly appendix a).                                 */
/*                                                                   */
/*  History -                                                        */
/*   version   date       programmer        change description       */
/*   -------   ----       ----------        ------------------       */
/*     1.0    8/17/90   ray/cartwright  Initial version.             */
/*     1.1   12/18/90   ray/cartwright  Enhanced documentation;      */
/*                                      use ascii input for          */
/*                                      portability.                 */
/*     1.2   12/31/91   r. ray        Changed min w to 0.5 (old      */
/*                                    value of 0.2 allowed too much  */
/*                                    extrapolation).                */
/*     2.0    2/25/94   r. eanes      Changed to use t/p derived     */
/*                                    orthoweights.                  */
/*     2.1    3/02/94   r. eanes      Generalized grid size up to    */
/*                                    1 x 1 degree and added data    */
/*                                    statements for the harmonic    */
/*                                    amplitudes and phases (a la    */
/*                                    the jpl mods). Iinterpolate    */
/*                                    orthoweights instead of tides. */
/*     2.2    3/03/94   r. eanes      Return orthoweights (u,v) for  */
/*                                    other use. Added itype to      */
/*                                    switch to c&r style            */
/*                                    orthoweights. Change min w     */
/*                                    (wmin) back to 0.2.            */
/*     2.3    4/18/94   r. eanes      Take out t/p kludge for        */
/*                                    latitudes > 66.2               */
/*     2.4    5/20/94   r. eanes      Change treatment of            */
/*                                    extrapolating past the         */
/*                                    southern and northernmost grid */
/*                                    points. Set wmin to 0.1.       */
/*     2.5                            Added r.ray's radiation tide   */
/*                                    correction. wmin changed to    */
/*                                    0.5 (again).                   */
/*     3.0    1/24/95   r. eanes      Parameter statements for       */
/*                                    dimensions and unit number.    */
/*                                    Compact format option for      */
/*                                    orthoweight file. Increased    */
/*                                    dimensions to allow 0.5 x 0.5  */
/*                                    deg grid. Changed wmin to 0.2  */
/*                                    to make rsc94b work at gauges. */
/*-------------------------------------------------------------------*/
int tptide_eot(double *dlat, double *dlon, short *isdata, double *u,
               double *v, ortho_cons_t *oc)
{
  double atan(), cos(), sin();

  /* Initialized data */
  double undef = 99999.0;
  double eps = 1e-10;
  double wmin = 0.2;

  /* Local variables */
  int jlat1, jlat2, ilon1, ilon2, i, j;
  int nx = oc->nx, ny = oc->ny;
  double xlat, xlon;
  double xlats, xlatn;
  double w;
  double w1, w2, w3, w4;
  double dx, dy;
  double wx1, wx2, wy1, wy2;
  double *uv = oc->uv;
  double latmin = oc->latmin, latmax = oc->latmax;
  double lonmin = oc->lonmin, lonmax = oc->lonmax;

  /* Parameter adjustments */
  v -= 4;
  u -= 4;

  *isdata = TRUE;
  dx = (lonmax - lonmin) / (nx - 1);
  dy = (latmax - latmin) / (ny - 1);
  xlats = latmin - dy;          /* southern limit of extr */
  xlatn = latmax + dy;          /* northern limit of extr */

  for (i = 1; i <= 3; ++i) {
    for (j = 1; j <= 2; ++j) {
      u[i + j * 3] = 0.0;
      v[i + j * 3] = 0.0;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Compute indices for desired position */
  xlat = *dlat;
  jlat1 = (int)((xlat - latmin) / dy) + 1;
  if (jlat1 == 0 && xlat > xlats) {
    jlat1 = 1;
  }
  /* allow extr */
  if (jlat1 == ny && xlat < xlatn) {
    jlat1 = ny - 1;
  }
  /* allow extr */
  jlat2 = jlat1 + 1;
  if (jlat1 < 1 || jlat2 > ny) {
    *isdata = FALSE;
    /* If the orthoweight file repesents differences return zero */
    /* instead of saying there is no data.  */
    if (oc->idiff == 1) {
      *isdata = TRUE;
    }
    return 0;
  }
  xlon = *dlon;
  if (xlon < lonmin) {
    xlon = xlon + 360. - eps;
  }
  /* eps to avoid pr */
  ilon1 = (int)((xlon - lonmin) / dx) + 1;
  if (ilon1 > nx) {
    --ilon1;
  }
  /* avoid problem f */
  ilon2 = ilon1 + 1;
  if (ilon2 > nx) {
    ilon2 = 1;
  }
  /* longitude grid */
  if (ilon1 < 1 || ilon1 > nx) {
    hd_quit("tptide_eot: error in tptide: incorrect longitude input\n");
    exit(0);
  }
  if (uv[(((((ilon1 + jlat1 * 720) << 1) + 1) * 3 + 1) << 1) - 8660] ==
      undef &&
      uv[(((((ilon2 + jlat1 * 720) << 1) + 1) * 3 + 1) << 1) - 8660] ==
      undef &&
      uv[(((((ilon1 + jlat2 * 720) << 1) + 1) * 3 + 1) << 1) - 8660] ==
      undef &&
      uv[(((((ilon2 + jlat2 * 720) << 1) + 1) * 3 + 1) << 1) - 8660] ==
      undef) {
    *isdata = FALSE;
    if (oc->idiff == 1) {
      *isdata = TRUE;
    }
    return 0;
  }

  /*-----------------------------------------------------------------*/
  /* distances from corners in units of cell size */
  wx1 = (dx - (xlon - (double)(ilon1 - 1) * dx - lonmin)) / dx;
  wx2 = 1.0 - wx1;
  wy1 = (dy - (xlat - (double)(jlat1 - 1) * dy - latmin)) / dy;
  wy2 = 1.0 - wy1;
  if (wy2 < 0.0) {              /* north corner weights zero in south */
    wy1 += wy2;
    wy2 = 0.0;
  }
  if (wy1 < 0.0) {              /* south corner weights zero in north */
    wy2 = wy1 + wy2;
    wy1 = 0.0;
  }
  /* interpolation weights: */
  /* w1,w2,w3,w4 are for northwest, northeast, southeast, southwest */
  /* corners.  */
  w1 = wx1 * wy2;
  w2 = wx2 * wy2;
  w3 = wx2 * wy1;
  w4 = wx1 * wy1;
  /* initialize sum of weights and orthoweights */
  w = 0.;
  for (i = 1; i <= 3; ++i) {
    for (j = 1; j <= 2; ++j) {
      u[i + j * 3] = 0.;
      v[i + j * 3] = 0.;
    }
  }
  /*-----------------------------------------------------------------*/
  /* interpolate the orthoweights */
  if (uv[(((((ilon1 + jlat1 * 720) << 1) + 1) * 3 + 1) << 1) - 8660] !=
      undef) {
    for (i = 1; i <= 3; ++i) {
      u[i + 3] = uv[((i + ((((ilon1 + jlat1 * 720) << 1) + 1) * 3)) << 1)
                    - 8660] * w4;
      v[i + 3] = uv[((i + ((((ilon1 + jlat1 * 720) << 1) + 1) * 3)) << 1)
                    - 8659] * w4;
      u[i + 6] = uv[((i + ((((ilon1 + jlat1 * 720) << 1) + 2) * 3)) << 1)
                    - 8660] * w4;
      v[i + 6] = uv[((i + ((((ilon1 + jlat1 * 720) << 1) + 2) * 3)) << 1)
                    - 8659] * w4;
    }
    w = w4;
  }

  if (uv[(((((ilon1 + jlat2 * 720) << 1) + 1) * 3 + 1) << 1) - 8660] !=
      undef) {
    for (i = 1; i <= 3; ++i) {
      u[i + 3] = uv[((i + ((((ilon1 + jlat2 * 720) << 1) + 1) * 3)) << 1)
                    - 8660] * w1 + u[i + 3];
      v[i + 3] = uv[((i + ((((ilon1 + jlat2 * 720) << 1) + 1) * 3)) << 1)
                    - 8659] * w1 + v[i + 3];
      u[i + 6] = uv[((i + ((((ilon1 + jlat2 * 720) << 1) + 2) * 3)) << 1)
                    - 8660] * w1 + u[i + 6];
      v[i + 6] = uv[((i + ((((ilon1 + jlat2 * 720) << 1) + 2) * 3)) << 1)
                    - 8659] * w1 + v[i + 6];
    }
    w += w1;
  }

  if (uv[(((((ilon2 + jlat2 * 720) << 1) + 1) * 3 + 1) << 1) - 8660] !=
      undef) {
    for (i = 1; i <= 3; ++i) {
      u[i + 3] = uv[((i + ((((ilon2 + jlat2 * 720) << 1) + 1) * 3)) << 1)
                    - 8660] * w2 + u[i + 3];
      v[i + 3] = uv[((i + ((((ilon2 + jlat2 * 720) << 1) + 1) * 3)) << 1)
                    - 8659] * w2 + v[i + 3];
      u[i + 6] = uv[((i + ((((ilon2 + jlat2 * 720) << 1) + 2) * 3)) << 1)
                    - 8660] * w2 + u[i + 6];
      v[i + 6] = uv[((i + ((((ilon2 + jlat2 * 720) << 1) + 2) * 3)) << 1)
                    - 8659] * w2 + v[i + 6];
    }
    w += w2;
  }

  if (uv[(((((ilon2 + jlat1 * 720) << 1) + 1) * 3 + 1) << 1) - 8660] !=
      undef) {
    for (i = 1; i <= 3; ++i) {
      u[i + 3] = uv[((i + ((((ilon2 + jlat1 * 720) << 1) + 1) * 3)) << 1)
                    - 8660] * w3 + u[i + 3];
      v[i + 3] = uv[((i + ((((ilon2 + jlat1 * 720) << 1) + 1) * 3)) << 1)
                    - 8659] * w3 + v[i + 3];
      u[i + 6] = uv[((i + ((((ilon2 + jlat1 * 720) << 1) + 2) * 3)) << 1)
                    - 8660] * w3 + u[i + 6];
      v[i + 6] = uv[((i + ((((ilon2 + jlat1 * 720) << 1) + 2) * 3)) << 1)
                    - 8659] * w3 + v[i + 6];
    }
    w += w3;
  }

  if (w > wmin) {
    for (i = 1; i <= 3; ++i) {
      for (j = 1; j <= 2; ++j) {
        u[i + j * 3] /= w;
        v[i + j * 3] /= w;
      }
    }
  } else {
    *isdata = FALSE;
    if (oc->idiff == 1) {
      *isdata = TRUE;
      for (i = 1; i <= 3; ++i) {
        for (j = 1; j <= 2; ++j) {
          u[i + j * 3] = 0.0;
          v[i + j * 3] = 0.0;
        }
      }
    }
  }
  return 0;
}

/* END tptide_eot()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computes the tide at time ts from the orthoweights u,v.           */
/*  This is a kernel routine, meant to be called only by tptide.     */
/*-------------------------------------------------------------------*/
int
tptid1(double *u, double *v, double *ts, short *pseudo, short *radiation,
       double *tide)
{
  /* Initialized data */
  double phc[4] = { 290.21, 280.12, 274.35, 343.51 };
  double dpd[4] = { 13.1763965, .9856473, .1114041, .0529539 };
  double u00[2] = { .0298, .02 };
  double u20[2] = { .1408, .0905 };
  double u21[2] = { .0805, .0638 };
  double u40[2] = { .6002, .3476 };
  double u41[2] = { .3025, .1645 };
  double v41[2] = { .1517, .0923 };
  double tc = 4043174400.;
  double tslast = -9.99e10;
  short init = TRUE;
  int indx[240] = { -3, -3, -2, -2, -2, -2, -1,
    -1, -1, -1, -1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3,
    3, -3, -3, -2, -2, -2,
    -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2,
    3, 3, 3, 0, 2, 0, 0, 2,
    2, 0, 0, 0, 0, 2, 0, 0, 0, 2, -3, -2, -2, 0, 0, 0, 1, 2, -2, 0, 0, -2,
    0, 0, 0, 0, 2, 0, 2, 3,
    -1, 0, 0, 1, 2, 3, -2, -1, 0, 0, 1, -2, 0, 0, 0, -3, -2, -1, 0, 0, 0,
    0, -2, 0, 0, 2, 0, 1, 1,
    -1, -1, 0, 0, 0, 2, 0, -1, 1, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1,
    0, -2, 0, 0, 3, 1, 2,
    0, 0, 1, 1, 1, 1, -1, -1, 2, 0, 0, 0, 0, 1, -1, 1, 1, 0, 0, 0, 0, 0, 0,
    0, 1, -1, -1, 0, 0, -1,
    0, -1, 0, -2, -1, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, -1, 0, 1, 0, 0, 0, 0,
    1, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, -1, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1, 0, 1,
    2, 0, 0, 1
  };
  double amp[60] = { .00663, .00802, .00947,
    .05019, .0018, .00954, .00152, .04946, .26218, .00171, .00343, .00741,
    .02062, .00414, .00394, .00713, .00137, .122, .0073, .36874, .05002,
    .00293, .00524, .00395, .02061, .00409, .00342, .00169, .01128, .00723,
    .0018, .00467, .01601, .01932, .0013, .00102, .00451, .121, .00113,
    .02298,
    .00106, .0019, .00218, .02358, .63194, .00193, .00466, .01786, .00447,
    .00197, .01718, .29402, .00305, .00102, .07994, .02382, .00259, 8.6e-4,
    .00447, .00195
  };
  double pha[60] = { 90., 90., 90., 90., 90., 90.,
    270., 90., 90., 270., 270., 270., 270., 270., 270., 13., 270., 90.,
    90., 270.,
    270., -13., 270., 270., 270., 270., 270., 270., 270., 270., 0., 0., 0.,
    0., 77.,
    103., 180., 0., 77., 0., 77., 180., 103., 180., 0., 77., 180., 180.,
    0., 0.,
    283., 0., 262., 180., 0., 0., 0., 0., 0., 0.
  };
  /* 
     data indx !( l ,m,n) -3,-3,-2,-2,-2,-2,-1,-1,-1,-1,-1, 0, 0, 0, 0, !( 
     1:15,1,1) 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, !(16:30,1,1)
     -3,-3,-2,-2,-2,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, !( 1:15,2,1) 0, 1, 1,
     1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, !(16:30,2,1) 0, 2, 0, 0, 2, 2, 0, 
     0, 0, 0, 2, 0, 0, 0, 2, !( 1:15,1,2) -3,-2,-2, 0, 0, 0, 1, 2,-2, 0,
     0,-2, 0, 0, 0, !(16:30,1,2) 0, 2, 0, 2, 3,-1, 0, 0, 1, 2, 3,-2,-1, 0, 
     0, !( 1:15,2,2) 1,-2, 0, 0, 0,-3,-2,-1, 0, 0, 0, 0,-2, 0, 0,
     !(16:30,2,2) 2, 0, 1, 1,-1,-1, 0, 0, 0, 2, 0,-1, 1, 1,-1, !(
     1:15,1,3) 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 0,-2, 0, 0, !(16:30,1,3)
     3, 1, 2, 0, 0, 1, 1, 1, 1,-1,-1, 2, 0, 0, 0, !( 1:15,2,3) 0, 1,-1, 1, 
     1, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, !(16:30,2,3) 0, 0,-1, 0,-1, 0,-2,-1, 
     0, 0, 0, 0, 0, 1, 0, !( 1:15,1,4) 0,-1, 0,-1, 0, 1, 0, 0, 0, 0, 1, 0, 
     0, 0, 1, !(16:30,1,4) 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0,-1, 0, !( 
     1:15,2,4) 0, 0, 0, 0, 1, 0, 0, 0,-1, 0, 1, 2, 0, 0, 1/ !(16:30,2,4)
     data amp !( l ,m) 0.00663,0.00802,0.00947,0.05019,0.0018 ,0.00954, !( 
     1: 6,1) 0.00152,0.04946,0.26218,0.00171,0.00343,0.00741, !( 7:12,1)
     0.02062,0.00414,0.00394,0.00713,0.00137,0.122 , !(13:18,1) 0.0073
     ,0.36874,0.05002,0.00293,0.00524,0.00395, !(19:24,1)
     0.02061,0.00409,0.00342,0.00169,0.01128,0.00723, !(25:30,1) 0.0018
     ,0.00467,0.01601,0.01932,0.0013 ,0.00102, !( 1: 6,2) 0.00451,0.121
     ,0.00113,0.02298,0.00106,0.0019 , !( 7:12,2)
     0.00218,0.02358,0.63194,0.00193,0.00466,0.01786, !(13:18,2)
     0.00447,0.00197,0.01718,0.29402,0.00305,0.00102, !(19:24,2)
     0.07994,0.02382,0.00259,0.00086,0.00447,0.00195/ !(25:30,2) data pha
     !( l ,m) 90., 90., 90., 90., 90., 90.,270., 90., 90.,270., !( 1:10,1)
     270.,270.,270.,270.,270., 13.,270., 90., 90.,270., !(11:20,1)
     270.,-13.,270.,270.,270.,270.,270.,270.,270.,270., !(21:30,1) 0., 0.,
     0., 0., 77.,103.,180., 0., 77., 0., !( 1:10,2) 77.,180.,103.,180., 0., 
     77.,180.,180., 0., 0., !(11:20,2) 283., 0.,262.,180., 0., 0., 0., 0.,
     0., 0./ !(21:30,2)

     long period LP 55.565 1 Sa 56.554 2 Ssa 57.555 3 TERa 58.554 4 Mm
     65.455 5 Mf 75.555 6 TERm 85.455 7 93a 93.555 8 diurnal 2Q1 125.755 9
     Q1 135.655 10 * O1 145.555 11 * M1 155.655 12 P1 163.555 13 * S1
     164.556 14 * K1 165.555 15 * PHI1167.555 16 J1 175.455 17 OO1 185.555
     18 NU1 195.455 19 semi-diurnal 227 227.655 20 2N2 235.755 21 * MU2
     237.555 22 * N2 245.655 23 * NU2 247.455 24 * M2 255.555 25 * L2
     265.455 26 * T2 272.556 27 * S2 273.555 28 * K2 275.555 29 * 285
     285.455 30 mask[i]=1 : tidal constituent included; mask[i]=0 :
     constituent omitted 0 1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 20 1 2 3 
     4 5 6 7 8 9 0 */
  int mask[31] =
    { 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1
  };

  /* Local variables */
  int i, k, l, m, n;
  double fdpd, dpld, crct, freq[60], shpn[4], crst, srct, srst, a, e = 0;
  double p[6], q[6], otide, theta, twopi, fd, at[3], bt[3], ct[60];
  double td, ph, cr, sr, st[60], omegat, rad;

  /* Parameter adjustments */
  v -= 4;
  u -= 4;

  if (init) {
    init = FALSE;
    rad = atan(1.0) / 45.0;
    dpld = 360.0 - dpd[0] + dpd[1];
    twopi = atan(1.0) * 8.0;
    fdpd = 0.;
    /* printf("semidiurnal & diurnal tidal potential terms:\n"); */
    for (m = 1; m <= 2; ++m) {
      fdpd += dpld;
      for (l = 1; l <= 30; ++l) {
        if (!mask[l])
          continue;
        theta = fdpd;
        for (n = 1; n <= 4; ++n) {
          theta += (double)indx[l + (m + (n << 1)) * 30 - 91] * dpd[n - 1];
        }
        freq[l + m * 30 - 31] = theta;
        omegat = theta * 2. * rad;
        ct[l + m * 30 - 31] = cos(omegat) * 2.;
        st[l + m * 30 - 31] = sin(omegat) * 2.;
        if (*radiation && m == 2 && (l == 21 || l == 22 || l == 23)) {
          amp[l + m * 30 - 31] *= .97;
          pha[l + m * 30 - 31] += -5.9;
          /* write(6,503) l,m,(indx(l,m,n),n=1,4),freq(l,m),amp(l,m), */
/*     *                      pha(l,m) */
        } else {
          /* write(6,502)
             l,m,(indx(l,m,n),n=1,4),freq(l,m),amp(l,m),pha(l,m) */
        }
      }
    }
  }
  if (*ts != tslast) {
    td = (*ts - tc) / 86400.;
    /* Compute 4 principal mean longitudes for a given time td */
    for (n = 1; n <= 4; ++n) {
      ph = phc[n - 1] + td * dpd[n - 1];
      shpn[n - 1] = fmod(ph, c_b129);
    }
    fd = td - floor(td);
    e = fd * 360. - shpn[0] + shpn[1];
  }

  otide = 0.;
  for (m = 1; m <= 2; ++m) {
    if (*ts != tslast) {
      /* Compute orthotides p(3),q(3) for given time ts in terms */
      /* of potential amplitudes at(3), bt(3) */
      for (k = 1; k <= 3; ++k) {
        at[k - 1] = 0.;
        bt[k - 1] = 0.;
      }
      for (l = 1; l <= 30; ++l) {
        if (!mask[l])
          continue;
        ph = m * e + pha[l + m * 30 - 31];
        for (n = 1; n <= 4; ++n) {
          i = indx[l + (m + (n << 1)) * 30 - 91];
          if (i == 0)
            break;
          ph += (double)i *shpn[n - 1];
        }
        theta = ph * rad;
        a = amp[l + m * 30 - 31] * 100.;
        cr = a * cos(theta);
        sr = a * sin(theta);
        crct = cr * ct[l + m * 30 - 31];
        crst = cr * st[l + m * 30 - 31];
        srct = sr * ct[l + m * 30 - 31];
        srst = sr * st[l + m * 30 - 31];
        at[0] += cr;
        bt[0] -= sr;
        at[1] += crct;
        bt[1] -= srct;
        at[2] += crst;
        bt[2] -= srst;
      }
      p[m * 3 - 3] = u00[m - 1] * at[0];
      q[m * 3 - 3] = u00[m - 1] * bt[0];
      p[m * 3 - 2] = u20[m - 1] * at[0] - u21[m - 1] * at[1];
      q[m * 3 - 2] = u20[m - 1] * bt[0] - u21[m - 1] * bt[1];
      if (*pseudo) {
        /* Cartwright and ray pseudo-orthoweights */
        p[m * 3 - 1] = u40[m - 1] * at[0] - u41[m - 1] * at[1] -
          v41[m - 1] * at[2];
        q[m * 3 - 1] = u40[m - 1] * bt[0] - u41[m - 1] * bt[1] -
          v41[m - 1] * bt[2];
      } else {
        /* Correct orthoweight definition */
        p[m * 3 - 1] = u40[m - 1] * at[0] - u41[m - 1] * at[1] +
          v41[m - 1] * at[2];
        q[m * 3 - 1] = u40[m - 1] * bt[0] - u41[m - 1] * bt[1] +
          v41[m - 1] * bt[2];
      }
    }
    for (k = 1; k <= 3; ++k) {
      otide = otide + u[k + m * 3] * p[k + m * 3 - 4] + v[k + m * 3] *
        q[k + m * 3 - 4];
    }
  }
  *tide = otide;
  tslast = *ts;
  return 0;
}

/* END tptid1()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computes admittance and in-phase and quadrature components from   */
/* orthoweights.                                                     */
/* Coded by Richard Eanes, UT/CSR, 1993 and 1994                     */
/*                                                                   */
/* Inputs:                                                           */
/*  pseudo   -- Logical variable, true for CR91 style pseudo-        */
/*              orthoweights.                                        */
/*              Use output from tptide for this.                     */
/*  radiation - Logical variable, true if the orthoweights are       */
/*              corrected for the radiation potential at T2,S2,R2    */
/*              (R2 not in this list).                               */
/*  u,v       - The ortwhoweights for the diurnal and semidiurnal    */
/*              bands.                                               */
/*                                                                   */
/* Outputs:                                                          */
/*  x,y       - The tidal admittances at a list of frequencies       */
/*              spanning the diurnal and semidiurnal bands (see data */
/*              below).                                              */
/*  h1,h2     - In-phase and quadrature amplitudes at the same       */
/*              frequencies (cm).                                    */
/*  amp,pha   - Amplitude and phase at the same list of frequencies. */
/*              (cm,deg).                                            */
/* The orthotide weight factors (signs with + were wrong in C&R 90   */
/*             Table A1).                                            */
/*-------------------------------------------------------------------*/
int
admit2(double *dlat, double *dlon, short *pseudo, short *radiation,
       double *u, double *v, double *x, double *y, double *h1, double *h2,
       double *amp, double *pha)
{
  double atan(), cos(), sin(), sqrt(), atan2();

  /* Initialized data */
  double sp[12] = { .0298, .1408, .0805, .6002, .3025, .1517, .02, .0905,
    .0638, .3476, .1645, .0923
  };
  double freq[30] = { 1.47094e-4, .002737779, .005475819,
    .008213597, .036291647, .073202203, .10949385, .140928587, .856952384,
    .893244048, .929535695, .966446233, .997262079, 1., 1.002737898,
    1.008213699, 1.039029527, 1.075940083, 1.11223173, 1.828255527,
    1.859690264, 1.864547173, 1.895981946, 1.90083882, 1.932273593,
    1.968565204, 1.997262163, 2., 2.005475795, 2.041767407
  };
  double hsm[30] = { .0279, -.0049, -.031, -.0018, -.0352, -.0666,
    -.0128, -.002, -.0066, -.0502, -.2622, .0206, -.122, .0029, .3688,
    .0053,
    .0206, .0113, .0022, .0047, .016, .0193, .121, .023, .6319, -.0179,
    .0172,
    .294, .08, .0045
  };
  /* 
     char *tname[MAXCONSTIT+1]= { "LP", "Sa", "Ssa", "TERa", "Mm", "Mf",
     "TERm", "93a", "2Q1", "Q1", "O1", "M1", "P1", "S1", "K1", "PHI1",
     "J1", "OO1", "NU1", "227", "2N2", "MU2", "N2", "NU2", "M2", "L2",
     "T2", "S2", "K2", "285"}; */
  double dt = 2.;
  double rad_fac = .97;
  double rad_phase = 5.9;
  short first = TRUE;
  short mn = 1;                 /* 0 for all constituents, 1 for
                                   diurnal/semi-diurnal */

  /* System generated locals */
  int i1;
  double d__1, d__2;

  /* Local variables */
  int i, m;
  double theta, twopi;
  int mf;
  double cw, sw, xp1[30], xp2[30], deg, fxy[30];

  /* Parameter adjustments */
  --pha;
  --amp;
  --h2;
  --h1;
  --y;
  --x;
  v -= 4;
  u -= 4;

  /* Function Body */
  /* name f(cycle/d) period hs (m) phase */
  /* long period */
  /* diurnal */
  /* semi-diurnal */
  /* Ray-Sanchez-Cartrwigh */
  /* ...compute constants depending only on frequency */
  first = TRUE;
  if (first) {
    first = FALSE;
    twopi = atan(1.0) * 8.0;
    deg = 360.0 / twopi;
    for (m = mn; m <= 2; ++m) {
      mf = (int)pow((double)c_n1, (double)m);
      i1 = i_2[m];
      for (i = i_1[m]; i <= i1; ++i) {
        theta = twopi * freq[i - 1] * dt;
        cw = cos(theta) * (float)2.;
        sw = sin(theta) * (float)2.;
        xp1[i - 1] = sp[m * 6 - 5] - sp[m * 6 - 4] * cw;
        if (*pseudo) {
          /* Cartwright&Ray pseudo-orthoweights */
          xp2[i - 1] = sp[m * 6 - 3] - sp[m * 6 - 2] * cw -
            sp[m * 6 - 1] * sw;
        } else {
          /* True orthoweights */
          xp2[i - 1] = sp[m * 6 - 3] - sp[m * 6 - 2] * cw +
            sp[m * 6 - 1] * sw;
        }
        fxy[i - 1] = mf * (d__1 = hsm[i - 1], fabs(d__1)) * 1e3;
      }
    }
  }
  for (m = mn; m <= 2; ++m) {
    i1 = i_2[m];
    for (i = i_1[m]; i <= i1; ++i) {
      x[i] = sp[m * 6 - 6] * u[m * 3 + 1] + xp1[i - 1] *
        u[m * 3 + 2] + xp2[i - 1] * u[m * 3 + 3];
      y[i] = sp[m * 6 - 6] * v[m * 3 + 1] + xp1[i - 1] *
        v[m * 3 + 2] + xp2[i - 1] * v[m * 3 + 3];
      h1[i] = fxy[i - 1] * x[i];
      h2[i] = -fxy[i - 1] * y[i];
      /* Computing 2nd power */
      d__1 = h1[i];
      /* Computing 2nd power */
      d__2 = h2[i];
      amp[i] = sqrt(d__1 * d__1 + d__2 * d__2);
      if (amp[i] != 0.0) {
        pha[i] = atan2(h2[i], h1[i]) * deg;
        if (pha[i] < 0.0) {
          pha[i] += 360.0;
        }
      } else {
        pha[i] = 0.0;
      }
      /* Correct T2 and S2 for radiation potential */
      if (*radiation && (i == 27 || i == 28)) {
        amp[i] *= rad_fac;
        pha[i] += rad_phase;
        h1[i] = amp[i] * cos(pha[i] / deg);
        h2[i] = amp[i] * sin(pha[i] / deg);
      }
      d__1 = amp[i] / 1e3;

      /*
      pha[i]=zoneshift(tname[i-1],pha[i],0,10); 
      printf("%5s %8.2f %8.2f %10.7f %8.2f\n", tname[i-1], *dlat, *dlon, d__1, 
	     zoneshift(tname[i-1], pha[i], 0, 10));
      printf("%5s %8.2f %8.2f %10.7f %8.2f\n", tname[i-1], *dlat,* dlon, d__1, pha[i]);
      */
    }
  }
  return 0;
}

/* END admit2()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computes the long-period equilibrium ocean tides.                 */
/*  arguments -                                                      */
/*     name      type  i/o               description                 */
/*     ----      ----  ---               -----------                 */
/*     ts         d     i   modified julian day, in seconds,         */
/*                          denoting time at which tide is to be     */
/*                          computed.                                */
/*     dlat       d     i   latitude in degrees (positive north) for */
/*                          the position at which tide is computed.  */
/*     tlp        d     o   computed long-period tide (centimeters)  */
/*                                                                   */
/*  processing logic -                                               */
/*  Fifteen tidal spectral lines from the cartwright-tayler-edden    */
/*  tables are summed over to compute the long-period tide.          */
/*                                                                   */
/*  technical references -                                           */
/*     Cartwright & Tayler, Geophys. j. r.a.s., 23, 45, 1971.        */
/*     Cartwright & Edden, Geophys. j. r.a.s., 33, 253, 1973.        */
/*                                                                   */
/*  history -                                                        */
/*    version   date       programmer         change description     */
/*    -------   ----       ----------         ------------------     */
/*      1.0   11/27/90     d. cartwright    documented by r. ray     */
/*                                                                   */
/* solar perigee held constant                                       */
/*-------------------------------------------------------------------*/
int lpeqmt(double *ts, double *dlat, double *tlp)
{
  double cos(), sin();

  /* Initialized data */
  double phc[4] = { 290.21, 280.12, 274.35, 343.51 };
  double dpd[4] = { 13.1763965, .9856473, .1114041, .0529539 };
  double tc = 4043174400.;
  /* Local variables */
  double shpn[4];
  int n;
  double s, td, ph, zlp;

  /* Compute 4 principal mean longitudes in radians at time td */
  td = (*ts - tc) / 86400.;
  /* days past 1 jan 87 00:00 */
  for (n = 1; n <= 4; ++n) {
    ph = phc[n - 1] + td * dpd[n - 1];
    shpn[n - 1] = fmod(ph, c_b129) * .017453292519943;
  }

  /* Assemble long-period tide potential from 15 cte terms > 1 mm.  */
  /* Nodal term is included but not the constant term.  */
  zlp = cos(shpn[3]) * 2.79 - cos(shpn[1] - 4.9392817831438691) * .49 -
    cos(shpn[1] * 2.) * 3.1;
  ph = shpn[0];
  zlp = zlp - cos(ph - shpn[1] * 2. + shpn[2]) * .67 -
    (3.52 - cos(shpn[3]) * .46) * cos(ph - shpn[2]);
  ph += shpn[0];
  zlp = zlp - cos(ph) * 6.66 - cos(ph + shpn[3]) * 2.76 -
    cos(ph + shpn[3] * 2.) * .26 - cos(ph - shpn[1] * 2.) * .58 -
    cos(ph - shpn[2] * 2.) * .29;
  ph += shpn[0];
  zlp =
    zlp - cos(ph - shpn[2]) * 1.27 - cos(ph - shpn[2] + shpn[3]) * .53 -
    cos(ph - shpn[1] * 2. + shpn[2]) * .24;

  /* Multiply by gamma_2 * sqrt(5/4 pi) * p20(lat) */
  s = sin(*dlat * .017453292519943);
  *tlp = zlp * .437 * (s * 1.5 * s - .5);
  return 0;
}

/* lpeqmt()                                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Compute the hour, minute, and second.                             */
/*-------------------------------------------------------------------*/
int kalday(double *utc, int *mdyhms)
{
  /* System generated locals */
  double d__1;

  /* Local variables */
  int i;
  double t, xtime;
  int jd, lx, nx, iptime[5];
  double ptimes;

  /* Parameter adjustments */
  --mdyhms;

  /* Function Body */
  t = *utc;
L5:
  jd = (int)t;
  d__1 = (t - jd) * 24. + 12.;
  xtime = fmod(d__1, c_b153);
  iptime[3] = (int)xtime;
  xtime = (xtime - iptime[3]) * 60.;
  iptime[4] = (int)xtime;
  ptimes = (xtime - iptime[4]) * 60.;

  /* Compute the year, month, and day.  */
  if (iptime[3] < 12) {
    ++jd;
  }
  lx = jd + 68569;
  nx = (lx << 2) / 146097;
  lx -= (nx * 146097 + 3) / 4;
  iptime[2] = (lx + 1) * 4000 / 1461001;
  lx = lx - iptime[2] * 1461 / 4 + 31;
  iptime[0] = lx * 80 / 2447;
  iptime[1] = lx - iptime[0] * 2447 / 80;
  lx = iptime[0] / 11;
  iptime[0] = iptime[0] + 2 - lx * 12;
  iptime[2] = (nx - 49) * 100 + iptime[2] + lx;
  for (i = 1; i <= 5; ++i) {
    mdyhms[i] = iptime[i - 1];
  }
  mdyhms[3] += -1900;
  mdyhms[6] = (int)(ptimes + .5);

  if (mdyhms[4] == 23 && mdyhms[5] == 59 && mdyhms[6] == 60) {
    t += 1e-6;
    goto L5;
  }
  return 0;
}

/* END kalday()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
double zoneshift(char *name, double oldphase, int oldzone, int newzone)
{
  double newphase, sigma = 0.0;
  if (strcasecmp(name, "M2") == 0)
    sigma = 28.9841042;
  else if (strcasecmp(name, "S2") == 0)
    sigma = 30.0000000;
  else if (strcasecmp(name, "K1") == 0)
    sigma = 15.0410686;
  else if (strcasecmp(name, "O1") == 0)
    sigma = 13.9430356;
  else if (strcasecmp(name, "SA") == 0)
    sigma = 0.0410686;
  else if (strcasecmp(name, "SS") == 0)
    sigma = 0.0821373;
  else if ((strcasecmp(name, "MM") == 0) || (strcasecmp(name, "Mm") == 0))
    sigma = 0.5443747;
  else if (strcasecmp(name, "MS") == 0)
    sigma = 1.0158958;
  else if (strcasecmp(name, "MF") == 0)
    sigma = 1.0980331;
  else if (strcasecmp(name, "S1") == 0)
    sigma = 15.0000000;
  else if (strcasecmp(name, "Q1") == 0)
    sigma = 13.3986609;
  else if (strcasecmp(name, "P1") == 0)
    sigma = 14.9589314;
  else if (strcasecmp(name, "N2") == 0)
    sigma = 28.4397295;
  else if (strcasecmp(name, "NU2") == 0)
    sigma = 28.5125831;
  else if (strcasecmp(name, "K2") == 0)
    sigma = 30.0821373;
  else if (strcasecmp(name, "L2") == 0)
    sigma = 29.5284789;
  else if (strcasecmp(name, "2N2") == 0)
    sigma = 27.8953548;
  else if (strcasecmp(name, "MU2") == 0)
    sigma = 27.9682084;
  else if (strcasecmp(name, "T2") == 0)
    sigma = 29.9589333;
  else if (strcasecmp(name, "M4") == 0)
    sigma = 57.9682084;
  else if (strcasecmp(name, "MS4") == 0)
    sigma = 58.9841042;
  else if (strcasecmp(name, "2MS6") == 0)
    sigma = 87.9682084;
  else if (strcasecmp(name, "2Q1") == 0)
    sigma=12.8542862000;
  else if (strcasecmp(name, "PHI1") == 0)
    sigma=15.123205900;
  else if (strcasecmp(name, "J1") == 0)
    sigma=15.5854433000;
  else if ((strcasecmp(name, "OO1") == 0) || (strcasecmp(name, "001") == 0))
    sigma=16.1391017000;

  newphase = oldphase + sigma * (newzone - oldzone);
  if (newphase < 0)
    newphase = newphase + 360;
  if (newphase > 360)
    newphase = newphase - 360;

  return (newphase);
}

char *strlwr(char *in)
{
  int i;
  char c;
  for (i = 0; i < (int)strlen(in); i++) {
    c = in[i];
    if (c >= 65 && c <= 90) {
      c = c + 32;
      in[i] = (char)c;
    }
  }
  return (in);
}


/* Extracts a year from a time and units */
int time_to_yr(double t, char *u)
{
  char units[] = "1 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";
  double mult;
  int y, mo, d;
  int h, mi, s;
  double j = tm_time_to_julsecs(u);



  /* Read units, convert to seconds */
  if (sscanf(u, "%10s since %*d-%*d-%*d", &units[2]) != 1)
    quit("time_to_yr: Can't understand %s\n", u);
  if (!tm_scale_to_secs(units, &mult))
    quit("time_to_yr: Can't convert %s to seconds\n", units);

  /* Add time value to epoch */
  j += t * mult / 86400.0;

  /* Convert back to date */
  tm_to_julsecs(j, &y, &mo, &d, &h, &mi, &s);

  return(y);
}

/*-------------------------------------------------------------------*/
/* Computes the equilibrium tide                                     */
/*-------------------------------------------------------------------*/
void equ_tide_eval(geometry_t *window,    /* Window geometry         */
		   window_t *windat,      /* Window data             */
		   win_priv_t *wincon,    /* Window constants        */
		   double *equitide       /* Equilibrium tide array  */
		   )
{
  int n, c, cc, cs;
  int yr, mon, nday, jd; /* Year, month, day, Julian date            */
  double day;            /* Day of year                              */
  double hrs;            /* Hours at longitude lon                   */
  double ml = 7.35e22;   /* Mass of the moon                         */
  double ms = 1.99e30;   /* Mass of the sun                          */
  double me = 5.97e24;   /* Mass of the earth                        */
  double a = 6370997.0;  /* Radius of the earth                      */
  double Rl = 3844e5;    /* Mean distance of moon from earth         */
  double Rs = 1496e8;    /* Mean distance of sun from earth          */
  double es = 0.0167086; /* Solar eccentricity                       */
  double el = 0.0549;    /* Lunar eccentricity                       */
  double as = 149.60e9;  /* Sun semi-maor axis (m)                   */
  double lat;            /* Latitude (radians)                       */ 
  double lon;            /* Longitude (degrees)                      */
  double gmst;           /* Greenwhich mean sidereal time (hours)    */
  double Al;             /* Lunar right ascension of ascending node  */
  double elon, elat;     /* Lunar ecliptic longitude and latitude    */
  double mass[2];        /* Lunar and solar mass (kg)                */
  double radius[2];      /* Lunar and solar distance (m)             */
  double dec[2];         /* Lunar and solar declination (rad)        */
  double hrang[2];       /* Lunar and solar hour angle (rad)         */
  double jt;             /* Julian date                              */
  double d2r = PI/180.0; /* Degrees to radians                       */
  double C0, C1, C2;
  double d1, sins, coss;

  /* Set up arrays                                                   */
  mass[0] = ml; mass[1] = ms;
  radius[0] = Rl; radius[1] = Rs;
  el *= d2r;
  for (n = 0; n <= 1; n++)
    mass[n] = a * (mass[n] / me);

  /* Get the Julian date                                             */
  tm_time_to_ymd(windat->t, master->timeunit, &yr, &mon, &nday);
  jd = date_to_jul(mon, nday, yr);
  /* Julian date starts at midday - adjust half a day                */
  jt = (double)jd - 0.5;

  /* Get the Greenwhich mean sidereal time in hours                  */
  /* see https://aa.usno.navy.mil/faq/docs/GAST.php                  */
  jt += windat->days - floor(windat->days);
  gmst = fmod(18.697374558 + 24.06570982441908 * (jt - 2451545.0), 24.0) * PI / 12.0;

  /* Get the lunar right ascension, radius and declination           */
  moonvars(jt, &Al, &dec[0], &elon, &elat, &radius[0]);
  radius[0] *= 1e3;

  /* Get the equilibrium tide                                        */
  for (cc = 1; cc <= window->a2_t; cc++) {
    double ramp = (wincon->rampf & (TIDALH|TIDALC)) ? windat->rampval : 1.0;
    c = window->w2_t[cc];
    cs = window->m2d[c];
    lon = window->cellx[cs];
    lat = window->celly[cs] * d2r;

    /* Get the day of the year                                       */
    if (window->is_geog)
      dtime(NULL, master->timeunit, windat->t, &yr, &day, &lon);
    else
      continue;

    /* Get the solar hour angle                                      */
    nday = (int)day;
    hrs = 24 * (day - (double)nday);
    hrang[1] = (hrs - 12.0) * PI / 12.0;
    radius[1] = d2r * day * 360.0 / 365.25;
    radius[1] = as * (1.0 - es * es) / (1.0 + es *cos(radius[1]));
    
    /* Get the solar declination                                     */
    d1 = nday * 2 * PI / 365.0;
    dec[1] = 0.006918 + 0.070257 * sin(d1) - 0.399912 * cos(d1)
      + 0.000907 * sin(2 * d1) - 0.006758 * cos(2 * d1)
      + 0.00148 * sin(3 * d1) - 0.002697 * cos(3 * d1);

    /* Get the lunar hour angle (Pugh Eq. 3.20a)                     */
    hrang[0] = lon * d2r + gmst - Al;

    /* Compute the equibibrium tide contribution                     */
    equitide[cs] = 0.0;
    for (n = 0; n <= 1; n++) {
      /* Get the time dependent coefficients                         */
      d1 = pow(a / radius[n], 3.0);
      C0 = d1 * (1.5 * sin(dec[n]) * sin(dec[n]) - 0.5);
      C1 = d1 * (0.75 * sin(2.0*dec[n]) * sin(2.0*dec[n])*cos(hrang[n]));
      C2 = d1 * (0.75 * cos(dec[n]) * cos(dec[n])*cos(2.0*hrang[n]));

      sins = sin(lat)*sin(lat);
      coss = cos(lat)*cos(lat);
      equitide[cs] += mass[n] * (C0 * (1.5 * sins - 0.5) +
				 C1 * sin(2.0 * lat)+
				 C2 * coss);
      /*printf("%d %d %f %f %f : %f %f\n",c,n,dec[n],lat*180/PI,hrang[n]*180/PI,windat->eta[cs],eta);*/
    }
    equitide[c] *= ramp;
  }
  /* Set the ghost cells                                             */
  set_lateral_bc_eta(equitide, window->nbptS, window->bpt,
		     window->bin, window->bin2, 0);
}


double equ_tide_evalo(geometry_t *window, window_t *windat, int c)
{
  int n, cs = window->m2d[c];
  int yr, mon, nday, jd; /* Year, month, day, Julian date            */
  double day;            /* Day of year                              */
  double hrs;            /* Hours at longitude lon                   */
  double ml = 7.35e22;   /* Mass of the moon                         */
  double ms = 1.99e30;   /* Mass of the sun                          */
  double me = 5.97e24;   /* Mass of the earth                        */
  double a = 6370997.0;  /* Radius of the earth                      */
  double Rl = 3844e5;    /* Mean distance of moon from earth         */
  double Rs = 1496e8;    /* Mean distance of sun from earth          */
  double es = 0.0167086; /* Solar eccentricity                       */
  double el = 0.0549;    /* Lunar eccentricity                       */
  double as = 149.60e9;  /* Sun semi-maor axis (m)                   */
  double lat;            /* Latitude (radians)                       */ 
  double lon = window->cellx[cs];  /* Longitude (degrees)            */
  double gmst;           /* Greenwhich mean sidereal time (hours)    */
  double Al;             /* Lunar right ascension of ascending node  */
  double elon, elat;     /* Lunar ecliptic longitude and latitude    */
  double mass[2];        /* Lunar and solar mass (kg)                */
  double radius[2];      /* Lunar and solar distance (m)             */
  double dec[2];         /* Lunar and solar declination (rad)        */
  double hrang[2];       /* Lunar and solar hour angle (rad)         */
  double jt;             /* Julian date                              */
  double eta;            /* Equilibrium tide height (m)              */
  double w0 = 1.0;        /* Mean solar day, cycles/mean solar day   */
  double w1 = 0.9661369;  /* Mean lunar day                          */
  double w2 = 0.0366009;  /* Sidereal month                          */
  double w3 = 0.0027379;  /* Tropical year                           */
  double w4 = 0.0030937;  /* Lunar perigee                           */
  double w5 = 0.0001471;  /* Regression of lunar node                */
  double w0d = 15.0;      /* Mean solar day, deg/mean solar hour     */
  double w1d = 14.4921;   /* Mean lunar day                          */
  double w2d = 0.5490;    /* Sidereal month                          */
  double w3d = 0.0411;    /* Tropical year                           */
  double w4d = 0.0046;    /* Lunar perigee                           */
  double w5d = 0.0022;    /* Regression of lunar node                */
  double d2r = PI/180.0;  /* Degrees to radians                      */
  double C0, C1, C2;
  double d1, sins, coss;

  /* Set up arrays                                                   */
  mass[0] = ml; mass[1] = ms;
  radius[0] = Rl; radius[1] = Rs;
  lat = window->celly[cs] * d2r;
  el *= d2r;

  /* Get the Julian date                                             */
  tm_time_to_ymd(windat->t, master->timeunit, &yr, &mon, &nday);
  jd = date_to_jul(mon, nday, yr);
  /* Julian date starts at midday - adjust half a day                */
  jt = (double)jd - 0.5;

  /* Get the day of the year                                         */
  if (window->is_geog) {
    /*dtime(master->params->output_tunit, master->timeunit, windat->t, &yr, &day, NULL);*/
    dtime(NULL, master->timeunit, windat->t, &yr, &day, &lon);
  } else {
    hd_warn("Can't compute equilibrium tide from a non-geographic mesh.\n");
    return(0.0);
  }

  /* Get the solar hour angle                                        */
  nday = (int)day;
  hrs = 24 * (day - (double)nday);
  hrang[1] = (hrs - 12.0) * PI / 12.0;
  radius[1] = d2r * day * 360.0 / 365.25;
  radius[1] = as * (1.0 - es * es) / (1.0 + es *cos(radius[1]));

  /* Get the solar declination                                       */
  d1 = nday * 2 * PI / 365.0;
  dec[1] = 0.006918 + 0.070257 * sin(d1) - 0.399912 * cos(d1)
    + 0.000907 * sin(2 * d1) - 0.006758 * cos(2 * d1)
    + 0.00148 * sin(3 * d1) - 0.002697 * cos(3 * d1);

  /* Get the Greenwhich mean sidereal time in hours                  */
  /* see https://aa.usno.navy.mil/faq/docs/GAST.php                  */
  jt += windat->days - floor(windat->days);
  gmst = fmod(18.697374558 + 24.06570982441908 * (jt - 2451545.0), 24.0) * PI / 12.0;
  /*gmst = mjd2gmst(jt) * PI / 12.0;*/

  /* Get the lunar right ascension, radius and declination           */
  moonvars(jt, &Al, &dec[0], &elon, &elat, &radius[0]);
  radius[0] *= 1e3;

  /* Get the lunar hour angle (Pugh Eq. 3.20a)                       */
  /*hrang[0] = lon * d2r + d2r * (w0d + w3d) * gmst - Al;*/
  hrang[0] = lon * d2r + gmst - Al;

  /*if(window->wn==1&&c==100)printf("t=%f gmst=%f lon=%f lat=%f A=%f hrang[0]=%f %f %f %f\n",windat->days,gmst*12./PI,elon/d2r,elat/d2r,Al/d2r,hrang[0]/d2r,dec[0]/d2r,radius[0]/1e3,jt);*/

  eta = 0.0;
  for (n = 0; n <= 1; n++) {
    d1 = pow(a/radius[n],3.0);
    C0 = d1 * (1.5 * sin(dec[n]) * sin(dec[n]) - 0.5);
    C1 = d1 * (0.75 * sin(2.0*dec[n]) * sin(2.0*dec[n])*cos(hrang[n]));
    C2 = d1 * (0.75 * cos(dec[n]) * cos(dec[n])*cos(2.0*hrang[n]));

    /* Compute the equibibrium tide contribution                     */
    sins = sin(lat)*sin(lat);
    coss = cos(lat)*cos(lat);
    eta += a * (mass[n] / me) * (C0*(1.5*sins-0.5) +
				 C1*sin(2.0*lat)+
				 C2*coss);
    /*printf("%d %d %f %f %f : %f %f\n",c,n,dec[n],lat*180/PI,hrang[n]*180/PI,windat->eta[cs],eta);*/
  }
  return(eta);
}


#define r2r(x)  (2.0 * PI * (x - floor(x)))

/*-------------------------------------------------------------------*/
/* lunar ephemeris                                                   */
/* input jdate = julian date                                         */
/* output                                                            */
/*  rasc  = right ascension of the moon (radians)                    */
/*          (0 <= rasc <= 2 pi)                                      */
/*  decl  = declination of the moon (radians)                        */
/*          (-pi/2 <= decl <= pi/2)                                  */
/* https://au.mathworks.com/matlabcentral/fileexchange/              */
/*       39191-low-precision-ephemeris?focused=3789856&tab=function  */
/*-------------------------------------------------------------------*/
void moonvars(double jdate, 
	      double *rasc, 
	      double *decl, 
	      double *plon, 
	      double *plat,
	      double *dist
	      )
{
double atr = PI / 648000.;
/* time arguments                                                    */
double djd = jdate - 2451545.;
double t = (djd / 3652.5) + 1;
/* fundamental trig arguments (radians)                              */
double gm = r2r(0.374897 + 0.03629164709 * djd);
double gm2 = 2. * gm;
double gm3 = 3. * gm;
double fm = r2r(0.259091 + 0.0367481952 * djd);
double fm2 = 2. * fm;
double em = r2r(0.827362 + 0.03386319198 * djd);
double em2 = 2. * em;
double em4 = 4. * em;
double gs = r2r(0.993126 + 0.0027377785 * djd);
double lv = r2r(0.505498 + 0.00445046867 * djd);
double lm = r2r(0.606434 + 0.03660110129 * djd);
double ls = r2r(0.779072 + 0.00273790931 * djd);
double rm = r2r(0.347343 - 0.00014709391 * djd);
/* geocentric, ecliptic longitude of the moon (radians)              */
 double a, b, obliq, r;
double l = 22640. * sin(gm) - 4586. * sin(gm - em2) + 2370. * sin(em2);
 l = l + 769. * sin(gm2) - 668. * sin(gs) - 412. * sin(fm2);
 l = l - 212. * sin(gm2 - em2) - 206. * sin(gm - em2 + gs);
 l = l + 192. * sin(gm + em2) + 165. * sin(em2 - gs);
 l = l + 148. * sin(gm - gs) - 125. * sin(em) - 110. * sin(gm + gs);
 l = l - 55. * sin(fm2 - em2) - 45. * sin(gm + fm2) + 40. * sin(gm - fm2);
 l = l - 38. * sin(gm - em4) + 36. * sin(gm3) - 31. * sin(gm2 - em4);
 l = l + 28. * sin(gm - em2 - gs) - 24. * sin(em2 + gs) + 19. * sin(gm - em);
 l = l + 18. * sin(em + gs) + 15. * sin(gm + em2 - gs) + 14. * sin(gm2 + em2);
 l = l + 14. * sin(em4) - 13. * sin(gm3 - em2) - 17. * sin(rm);
 l = l - 11. * sin(gm + 16. * ls - 18. * lv) + 10. * sin(gm2 - gs) + 
   9. * sin(gm - fm2 - em2);
 l = l + 9. * (cos(gm + 16. * ls - 18. * lv) - sin(gm2 - em2 + gs)) -
   8. * sin(gm + em);
 l = l + 8. * (sin(2. * (em - gs)) - sin(gm2 + gs)) - 7. * (sin(2. * gs) + 
   sin(gm - 2. * (em - gs)) - sin(rm));
 l = l - 6. * (sin(gm - fm2 + em2) + sin(fm2 + em2)) -
   4. * (sin(gm - em4 + gs) - t * cos(gm + 16. * ls - 18 * lv));
 l = l - 4. * (sin(gm2 + fm2) - t * sin(gm + 16. * ls - 18. * lv));
 l = l + 3. * (sin(gm - 3. * em) - sin(gm + em2 + gs) -
    sin(gm2 - em4 + gs) + sin(gm - 2. * gs) + sin(gm - em2 - 2. * gs));
 l = l - 2. * (sin(gm2 - em2 - gs) + sin(fm2 - em2 + gs) - sin(gm + em4));
l = l + 2. * (sin(4. * gm) + sin(em4 - gs) + sin(gm2 - em));
 *plon = lm + atr * l;
 /* geocentric, ecliptic latitude of the moon (radians)              */
 b = 18461. * sin(fm) + 1010. * sin(gm + fm) + 1000. * sin(gm - fm);
 b = b - 624. * sin(fm - em2) - 199. * sin(gm - fm - em2) -
    167. * sin(gm + fm - em2);
 b = b + 117. * sin(fm + em2) + 62. * sin(gm2 + fm) + 33. * sin(gm - fm + em2);
 b = b + 32. * sin(gm2 - fm) - 30. * sin(fm - em2 + gs) -
    16. * sin(gm2 - em2 + fm);
 b = b + 15. * sin(gm + fm + em2) + 12. * sin(fm - em2 - gs) -
   9. * sin(gm - fm - em2 + gs);
 b = b - 8. * (sin(fm + rm) - sin(fm + em2 - gs)) -
    7. * sin(gm + fm - em2 + gs);
 b = b + 7. * (sin(gm + fm - gs) - sin(gm + fm - em4));
 b = b - 6. * (sin(fm + gs) + sin(3. * fm) - sin(gm - fm - gs));
 b = b - 5. * (sin(fm + em) + sin(gm + fm + gs) + sin(gm - fm + gs) -
	       sin(fm - gs) - sin(fm - em));
 b = b + 4. * (sin(gm3 + fm) - sin(fm - em4)) - 3. * (sin(gm - fm - em4) -
						      sin(gm - 3. * fm));
 b = b - 2. * (sin(gm2 - fm - em4) + sin(3. * fm - em2) - sin(gm2 - fm + em2) -
	       sin(gm - fm + em2 - gs));
 *plat = atr * (b + 2. * (sin(gm2 - fm - em2) + sin(gm3 - fm)));
 /* obliquity of the ecliptic (radians)                              */
 obliq = atr * (84428 - 47 * t + 9 * cos(rm));
 /* geocentric distance (kilometers)                                 */
 r = 60.36298 - 3.27746 * cos(gm) - .57994 * cos(gm - em2);
 r = r - .46357 * cos(em2) - .08904 * cos(gm2) + .03865 * cos(gm2 - em2);
 r = r - .03237 * cos(em2 - gs) - .02688 * cos(gm + em2) - 
   .02358 * cos(gm - em2 + gs);
 r = r - .0203 * cos(gm - gs) + .01719 * cos(em) + .01671 * cos(gm + gs);
 r = r + .01247 * cos(gm - fm2) + .00704 * cos(gs) + .00529 * cos(em2 + gs);
 r = r - .00524 * cos(gm - em4) + .00398 * cos(gm - em2 - gs) -
   .00366 * cos(gm3);
 r = r - .00295 * cos(gm2 - em4) - .00263 * cos(em + gs) +
   .00249 * cos(gm3 - em2);
 r = r - .00221 * cos(gm + em2 - gs) + .00185 * cos(fm2 - em2) -
    .00161 * cos(2 * (em - gs));
 r = r + 0.00147 * cos(gm + fm2 - em2) - 0.00142 * cos(em4) +
   0.00139 * cos(gm2 - em2 + gs);
 *dist = 6378.14 * (r - 0.00118 * cos(gm - em4 + gs) - 0.00116 * cos(gm2 + em2) -
		    0.0011 * cos(gm2 - gs));
 /* geocentric, equatorial right ascension and declination (radians) */
 a = sin(*plon) * cos(obliq) - tan(*plat) * sin(obliq);
 b = cos(*plon);
 *rasc = atan3(a, b);
 *decl = asin(sin(*plat) * cos(obliq) + cos(*plat) * sin(obliq) * sin(*plon));
}

/* Four quadrant inverse tangent */
/*  a = sine of angle */
/*  b = cosine of angle */
/* y = angle (radians; 0 =< c <= 2 * pi) */
double atan3 (double a, double b)
{
double epsilon = 0.0000000001;
double pidiv2 = 0.5 * PI;
 double y, c;

 if (fabs(a) < epsilon) {
    y = (1. - sgn(b)) * pidiv2;
    return(y);
 } else {
    c = (2. - sgn(a)) * pidiv2;
 }

 if (fabs(b) < epsilon) {
    y = c;
    return(y);
 } else
    y = c + sgn(a) * sgn(b) * (fabs(atan(a / b)) - pidiv2);
 return(y);
}

double sgn(double a)
{
  double ret;
  if (a > 0.0) ret = 1.0;
  if (a== 0.0) ret = 0.0;
  if (a < 0.0) ret = -1.0;
  return(ret);
}


/* mjd2gmst modified julian date to greenwich mean sidereal time.    */
/* gmst = mjd2gmst(mjd) converts modified julian date to greenwich   */
/* mean sidereal time using the algorithm from the Astronomical      */
/* Almanac 2002, pg. B6.                                             */
double mjd2gmst(double mjd)
{
  /* Calculate the greenwich mean sidereal time at midnight          */
  double mjd2000 = 51544.5; /* Modified Julian Date of Epoch J2000.0 */
  double int_mjd = floor(mjd);
  double frac_mjd = mjd-int_mjd;
  double Tu = (int_mjd - mjd2000)/36525.0;
  double gmst = 24110.54841 + Tu * (8640184.812866 +
				    Tu * (0.093104 -  Tu * 6.2e-6));
  /* Add the mean sidereal time interval from midnight to time       */
  gmst = fmod(gmst + frac_mjd*86400*1.00273790934,86400);
  /* Convert to hours                                                */
  gmst = gmst/3600;
  return(gmst);
}

/* Right ascension and declination of Moon pole in EME2000, rad      */
/* https://www.hq.nasa.gov/alsj/lunar_cmd_2005_jpl_d32296.pdf        */
void lunar_angles(double mjd, double *A, double *dec)
{
  double J2000 = 2451545.0; /* Reference epoch of J2000, Julian days */
  double D = mjd-J2000;     /* Days past epoch of J2000              */
  double T = D / 36525.0;   /* Julian centuries past J2000           */
  double E1 = 125.045 - 0.0529921 * D;
  double E2 = 250.089 - 0.1059842 * D;
  double E3 = 260.008 + 13.0120009 * D;
  double E4 = 176.625 + 13.3407154 * D;
  double E6 = 311.589 + 26.4057084 * D;
  double E7 = 134.963 + 13.0649930 * D;
  double E10 = 15.134 - 0.1589763 * D;
  double E13 = 25.053 + 12.9590088 * D;
  double d2r = PI / 180.0;
  double d0 = 66.5392;
  /*d0 = -10.218;*/

  *A = (269.9949 + 0.0031 * T - 3.8787*sin(E1*d2r) - 0.1204*sin(E2*d2r) + 
	0.0700*sin(E3*d2r) - 0.0172*sin(E4*d2r) + 0.0072*sin(E6*d2r) - 0.0052*sin(E10*d2r)
	+ 0.0043*sin(E13*d2r)) * d2r;

  *dec = (d0 + 0.0130 * T + 1.5419*cos(E1*d2r) + 0.0239*cos(E2*d2r) - 0.0278*cos(E3*d2r) +
	  0.0068*cos(E4*d2r) - 0.0029*cos(E6*d2r) + 0.0009*cos(E7*d2r) +
	  0.0008*cos(E10*d2r) - 0.0009*cos(E13*d2r)) * d2r;

}


/*-------------------------------------------------------------------*/
/* Moon: Computes the Moon's geocentric position using a low         */
/* precision analytical series.                                      */
/* Input: Mjd_TT=Terrestrial Time (Modified Julian Date)             */
//* Output: ecliptic longitude (L) and latitude (B)                  */
/* Modified from matlab code by M. Mahooti.                          */
/* https://www.mathworks.com/matlabcentral/fileexchange/             */
/*                        56041-moon-position?s_tid=prof_contriblnk  */
/*-------------------------------------------------------------------*/
void moon_pos(double mjd, double *L, double *B)
{
  double pi2 = 2*PI;              /* 2pi                             */
  double Rad = PI/180.0;          /* Radians per degree              */
  double Arcs = 3600.0*180.0/PI;  /* Arcseconds per radian           */
  double MJD_J2000 = 51544.5;     /* Modified Julian Date of J2000   */
  /* Constants                                                       */
  double ep = 23.43929111*Rad;    /* Obliquity of J2000 ecliptic     */
  double T = (mjd-MJD_J2000)/36525.0; /* Julian cent. since J2000    */
  /* Mean elements of lunar orbit                                    */
  double L_0 = frac(0.606433 + 1336.851344*T);    /* Mean longitude [rev] w.r.t. J2000 equinox */
  double l = pi2*frac(0.374897 + 1325.552410*T);  /* Moon's mean anomaly [rad] */
  double lp = pi2*frac(0.993133 + 99.997361*T);   /* Sun's mean anomaly [rad] */
  double D = pi2*frac(0.827361 + 1236.853086*T);  /* Diff. long. Moon-Sun [rad] */
  double F = pi2*frac(0.259086 + 1342.227825*T);  /* Argument of latitude */
  /* Ecliptic longitude (w.r.t. equinox of J2000)                    */
  double dL = +22640.0*sin(l) - 4586.0*sin(l-2*D) + 2370.0*sin(2*D) +  769.0*sin(2*l) -
    668.9*sin(lp) - 412.0*sin(2*F) - 212.0*sin(2*l-2*D) - 206.0*sin(l+lp-2*D) +
    192.0*sin(l+2*D) - 165.0*sin(lp-2*D) - 125.0*sin(D) - 110.0*sin(l+lp) +
    148.0*sin(l-lp) - 55.0*sin(2*F-2*D);
  *L = pi2 * frac(L_0 + dL/1296.0e3);         /* [rad]                */
  /* Ecliptic latitude                                                */
  double S = F + (dL+412*sin(2*F)+541*sin(lp)) / Arcs; 
  double h = F-2*D;
  double N = -526.0*sin(h) + 44.0*sin(l+h) - 31.0*sin(-l+h) - 23.0*sin(lp+h) +
    11.0*sin(-lp+h) - 25.0*sin(-2*l+F) + 21.0*sin(-l+F);
  *B = ( 18520.0*sin(S) + N ) / Arcs;         /* [rad]                */
}

double frac(double a)
{
  return(a-floor(a));
}
