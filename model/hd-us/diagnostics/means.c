/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/diagnostics/means.c
 *  
 *  Description:
 *  Routines relating to computation of averages
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: means.c 6413 2019-11-22 00:19:06Z her127 $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hd.h"

#define START_EPS (0.001)
#define nm3d 22
#define nm2d 5
char *mn_3d[nm3d] = {
  "u1mean", "u2mean", "wmean", "Kzmean", "temp_mean", "salt_mean",
  "flux_e1", "flux_e2", "flux_w", "flux_kz",
  "u1_adv", "u1_hdif", "u1_vdif", "u1_cor", "u1_btp", "u1_bcp",
  "u2_adv", "u2_hdif", "u2_vdif", "u2_cor", "u2_btp", "u2_bcp"
};
int mf_3d[nm3d] = {
  VEL3D, VEL3D, VEL3D, KZ_M, TS, TS,
  FLUX, FLUX, FLUX, FLUX,
  TENDENCY, TENDENCY, TENDENCY, TENDENCY, TENDENCY, TENDENCY, 
  TENDENCY, TENDENCY, TENDENCY, TENDENCY, TENDENCY, TENDENCY
};
char *mn_2d[nm2d] = {
  "eta_mean",
  "u1av_mean", "u2av_mean",
  "w1mean", "w2mean"
};
int mf_2d[nm2d] = {
  ETA_M,
  VEL2D, VEL2D,
  WIND, WIND
};

int read_mean_3d(master_t *master, double t, int sf, char *var, double *v, 
		 double ***vc);
void reset_means3d_m(master_t *master, int trm);
int read_mean_2d(master_t *master, double t, int sf, char *var, double *v, 
		 double **vc);
void reset_means2d_m(master_t *master, int trm);

/* Local functions */
static int means_init(sched_event_t *event);
double means_event(sched_event_t *event, double t);
static void means_cleanup(sched_event_t *event, double t);

typedef struct {
  master_t *master;             /* Grid associated with              */
  int means;                    /* Mean velocity diagnostic          */
  double means_dt;              /* Mean velocity averaging interval  */
  double means_next;            /* Next time for zeroing means       */
  double means_os;              /* Offset for restarts               */
  int means_tra;                /* Mean tracer number                */
  int tstep;                    /* Counter for TRANSPORT             */
  int tnext;                    /* Reinitialisation for TRANSPORT    */
} means_data_t;


/*-------------------------------------------------------------------*/
/* Initialize the mean quantities                                    */
/*-------------------------------------------------------------------*/
void init_means(master_t *master, parameters_t *params)
{
  means_data_t *mean = NULL;
  int c, cc;
  int ns, tn;

  master->means = params->means;

  /* Allocate memory for means structure, populate and register      */
  /* the scheduler events.                                           */
  mean = (means_data_t *)malloc(sizeof(means_data_t));
  memset(mean, 0, sizeof(means_data_t));

  mean->master = master;
  mean->means = master->means;

  /* Allocate                                                        */
  if (master->means & VEL3D) {
    master->ume = d_alloc_1d(geom->sze);
    memset(master->ume, 0, geom->sze * sizeof(double));
  }
  if (master->means & VEL2D) {
    master->uame = d_alloc_1d(geom->szeS);
    memset(master->uame, 0, geom->szeS * sizeof(double));
  }
  if (master->means & VOLFLUX || master->tmode & SP_U1VM) {
    master->u1vm = d_alloc_1d(geom->sze);
    memset(master->u1vm, 0, geom->sze * sizeof(double));
  }
  if (master->means & TRANSPORT) {
    master->means_dt = master->tratio * master->dt;
    mean->tnext = master->tratio;
    mean->tstep = 0;
  }
  if (master->means & NONE) return;

  /* Set the mean step                                               */
  if (strlen(params->means_dt)) {
    if (tm_scale_to_secs(params->means_dt, &master->means_dt))
      master->means_next = master->t + master->means_dt;
    else if (contains_token(params->means_dt, "YEARLY")) {
      master->means_dt = YEARLY;
      master->means_next = master->t + next_year(master->t, master->timeunit);
    } else if (contains_token(params->means_dt, "SEASONAL")) {
      master->means_dt = SEASONAL;
      master->meancs = d_alloc_1d(13);
      memset(master->meancs, 0, 13 * sizeof(double));
      master->means_next = master->t + next_season(master->t, 
						   master->timeunit, &c);
    } else if (contains_token(params->means_dt, "MONTHLY")) {
      master->means_dt = MONTHLY;
      master->meancs = d_alloc_1d(13);
      memset(master->meancs, 0, 13 * sizeof(double));
      master->means_next = master->t + next_month(master->t, 
						  master->timeunit, &c);
    } else if (contains_token(params->means_dt, "DAILY")) {
      master->means_dt = DAILY;
      master->meancs = d_alloc_1d(366);
      memset(master->meancs, 0, 366 * sizeof(double));
      master->means_next = master->t + next_day(master->t, 
						master->timeunit, &c);
    } else
      master->means_dt = 0.0;

    if (strlen(params->means_mc) && master->means_dt == SEASONAL || 
	       master->means_dt == MONTHLY || master->means_dt == DAILY) {
      char *fields[MAXSTRLEN * MAXNUMARGS];
      cc = parseline(params->means_mc, fields, MAXNUMARGS);
      for (c = 1; c <= cc; c++) {
	master->meancs[c] = atof(fields[c-1]);
      }
    }
    if (strlen(params->means_os)) {
      tm_scale_to_secs(params->stop_time, &master->means_os);
      for (cc = 1; cc <= geom->enonS; cc++) {
        master->meanc[cc] = master->means_os;
      }
    } else
      master->means_os = 0.0;
  } else
    master->means_dt = 0.0;
  mean->means_dt = master->means_dt;
  mean->means_next = master->means_next;
  mean->means_os = master->means_os;

  /* Set the mean transport step                                     */
  if (master->means & TRANSPORT) {
    if (master->means_dt == 0.0)
      hd_warn("constants : Must set MEAN_DT using TRANSPORT with MEAN\n");
    else {
      if (fmod(master->means_dt, master->dt) != 0.0) {
        hd_warn
          ("constants : MEAN_DT is not integrally divisable by DT : resetting.\n");
        c = (int)master->means_dt / master->dt;
        master->means_dt = (double)c *master->dt;
      }
    }
    /*master->tratio = (int)master->means_dt / master->dt;*/
  }

  /* Store all 3d mean tracers in tm_3d                              */
  master->ntm_3d = 0;
  if (master->means & VEL3D) master->ntm_3d += 3;
  if (master->means & TS) master->ntm_3d += 2;
  if (master->means & KZ_M) master->ntm_3d += 1;
  if (master->means & FLUX) master->ntm_3d += 4;
  if (master->means & TENDENCY) master->ntm_3d += 12;
  /*if (master->means & VOLFLUX) master->ntm_3d += 2;*/
  if (master->means & MTRA3D) {
    master->ntm_3d += 1;
    if ((master->means_tra = tracer_find_index(params->means_tra, master->ntr, master->trinfo_3d)) < 0) {
      hd_warn("compute_constants: Can't find 3D tracer %s for MEAN tracer.\n", params->means_tra);
      master->means &= ~MTRA3D;
    }
  }
  if (master->ntm_3d) {
    master->tm_3d = i_alloc_1d(master->ntm_3d);
    ns = 0;
    for (tn = 0; tn < nm3d; tn++) {
      c = tracer_find_index(mn_3d[tn], master->ntr, master->trinfo_3d);
      if (master->means & mf_3d[tn] && c >= 0) {
	master->tm_3d[ns] = c;
	ns++;
      }
    }
  }

  /* Store all 2d mean tracers in tm_2d                              */
  master->ntm_2d = 0;
  if (master->means & VEL2D) master->ntm_2d += 2;
  if (master->means & ETA_M) master->ntm_2d += 1;
  if (master->means & WIND) master->ntm_2d += 2;
  if (master->means & MTRA2D) {
    master->ntm_2d += 1;
    if ((master->means_tra = tracer_find_index(params->means_tra, master->ntrS, master->trinfo_2d)) < 0) {
      hd_warn("compute_constants: Can't find 2D tracer %s for MEAN tracer.\n", params->means_tra);
      master->means &= ~MTRA2D;
    }
  }
  if (master->ntm_2d) {
    master->tm_2d = i_alloc_1d(master->ntm_2d);
    ns = 0;
    for (tn = 0; tn < nm2d; tn++) {
      c = tracer_find_index(mn_2d[tn], master->ntrS, master->trinfo_2d);
      if (master->means & mf_2d[tn] && c >= 0) {
	master->tm_2d[ns] = c;
	ns++;
      }
    }
  }

  /* Register the scheduled function                                 */
  if (!(master->means & TRANSPORT))
    sched_register(schedule, "means", means_init, means_event, means_cleanup,
		   mean, NULL, NULL);
}

/* END init_means()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Scheduled initialsation function                                  */
/*-------------------------------------------------------------------*/
static int means_init(sched_event_t *event)
{
  return 1;
}

/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Re-initializes means on the master                                */
/*-------------------------------------------------------------------*/
double means_event(sched_event_t *event, double t)
{
  means_data_t *mean = (means_data_t *)schedGetPublicData(event);
  master_t *master = mean->master;

  if (master->means & NONE) return;

  if (t >= (event->next_event - SEPS)) {
    reset_means_m(master);
    event->next_event = master->means_next;
  }

  return event->next_event;
}

/* END means_event()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void means_cleanup(sched_event_t *event, double t)
{
  means_data_t *mean = (means_data_t *)schedGetPrivateData(event);

  if (mean != NULL) {
    free(mean);
  }
}

/* END mean_cleanup()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Re-initialises the master mean counters                           */
/*-------------------------------------------------------------------*/
void reset_means_m(master_t *master) 
{
  geometry_t *geom = master->geom;;
  int mon, day;
  double itc;
  int n;

  master->means |= RESET;    /* Set the RESET flag for transfers     */
  itc =  master->meanc[1];   /* Save the time counter                */
  /* Re-initialize the time counter.                                 */
  memset(master->meanc, 0, geom->szcS * sizeof(double));
  /* Set the time counter flag; this is set to 0 when meanc is       */
  /* re-initialized in reset_means_3d_m() or reset_means_2d_m()      */
  /* for SEASONAL or MONTHLY means below.                            */
  master->meanc[0] = 1.0;  
  /* Re-initialize 3D mean tracers                                   */
  for (n = 0; n < master->ntm_3d; n++)
    reset_means3d_m(master, master->tm_3d[n]);
  if (master->means & VEL3D) memset(master->ume, 0, geom->sze * sizeof(double));
  if (master->means & VOLFLUX) memset(master->u1vm, 0, geom->sze * sizeof(double));
  /* Re-initialize 2D mean tracers                                   */
  for (n = 0; n < master->ntm_2d; n++)
    reset_means2d_m(master, master->tm_2d[n]);
  /* Set the next re-initialization time                             */
  if (master->means_dt == SEASONAL) {
    master->means_next = master->t + next_season(master->t, 
						 master->timeunit, &mon);
    master->meancs[mon] = itc;
  } else if (master->means_dt == MONTHLY) {
    master->means_next = master->t + next_month(master->t, 
						master->timeunit, &mon);
    master->meancs[mon] = itc;
  } else if (master->means_dt == DAILY) {
    master->means_next = master->t + next_day(master->t, 
					      master->timeunit, &day);
    master->meancs[day] = itc;
  }
  else
    master->means_next = master->t + master->means_dt;
}

/* END reset_means_m()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to zero mean arrays at interval means_dt                  */
/*-------------------------------------------------------------------*/
void reset_means(geometry_t *window,  /* Window geometry             */
                 window_t *windat,    /* Window data                 */
                 win_priv_t *wincon,  /* Window constants            */
                 int mode             /* Operation mode              */
  )
{
  int c, cc;                    /* Counters                          */

  double ttime = 24;            /* Predominant tidal period (hours)  */
  double time;                  /* Cumulative time counter           */
  double *detadt = wincon->d4;  /* Mean forward eta difference       */
  double *w1 = wincon->d1;      /* EW component of wind stress       */
  double *w2 = wincon->d2;      /* NS component of wind stress       */

  if (wincon->means & NONE || !wincon->means_dt)
    return;
  if (mode & ALL) return;                  /* Legacy code            */

  /*-----------------------------------------------------------------*/
  /* Get the mean wind if required                                   */
  if (mode & WIND && wincon->means & WIND) {
    int cc, c;
    double ns = windat->dtf;
    vel_cen(window, windat, wincon, windat->wind1, NULL, w1, w2, NULL, NULL, 1);
    if (windat->w1m) {
      for (cc = 1; cc <= window->b2_t; cc++) {
        c = window->w2_t[cc];
	windat->w1m[c] = (windat->w1m[c] * windat->meanc[c] + 
			  w1[c] * ns) / (windat->meanc[c] + ns);
      }
    }
    if (windat->w2m) {
      for (cc = 1; cc <= window->b2_t; cc++) {
        c = window->w2_t[cc];
	windat->w2m[c] = (windat->w2m[c] * windat->meanc[c] + 
			  w2[c] * ns) / (windat->meanc[c] + ns);
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the mean TS if required                                     */
  if (mode & TS && wincon->means & (TS|MTRA3D|MTRA2D)) {
    int cc, c, cs;
    double ns = windat->dtf;
    if (windat->tempm && windat->saltm) {
      if (wincon->means & MMM) {
	for (cc = 1; cc <= window->b3_t; cc++) {
	  c = window->w3_t[cc];
	  cs = window->m2d[c];
	  windat->tempm[c] = max(windat->tempm[c], windat->temp[c]);
	  windat->saltm[c] = max(windat->saltm[c],windat->sal[c]);
	}
      } else {
	for (cc = 1; cc <= window->b3_t; cc++) {
	  c = window->w3_t[cc];
	  cs = window->m2d[c];
	  windat->tempm[c] = (windat->tempm[c] * windat->meanc[cs] + 
			      windat->temp[c] * ns) / (windat->meanc[cs] + ns);
	  windat->saltm[c] = (windat->saltm[c] * windat->meanc[cs] + 
			      windat->sal[c] * ns) / (windat->meanc[cs] + ns);			      
	}
      }
    }
    if (windat->tram) {
      int vs = (wincon->means & MTRA3D) ? window->b3_t : window->b2_t;
      int *vec = (wincon->means & MTRA3D) ? window->w3_t : window->w2_t;
      if (wincon->means & MMM) {
	for (cc = 1; cc <= vs; cc++) {
	  c = vec[cc];
	  windat->tram[c] = max(windat->tram[c], windat->tr_wc[wincon->means_tra][c]);
	}
      } else {
	for (cc = 1; cc <= vs; cc++) {
	  c = vec[cc];
	  cs = window->m2d[c];
	  windat->tram[c] = (windat->tram[c] * windat->meanc[cs] + 
			      windat->tr_wc[wincon->means_tra][c] * ns) / 
	    (windat->meanc[cs] + ns);
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Increment the mean time counter and zero arrays                 */
  if (mode & RESET) {
    /* Zero the mean arrays                                          */
    if (windat->t >= wincon->means_next) {
      if (wincon->means_dt == SEASONAL)
	wincon->means_next = windat->t + next_season(windat->t, 
						     wincon->timeunit, &c);
      else if (wincon->means_dt == MONTHLY)
	wincon->means_next = windat->t + next_month(windat->t, 
						    wincon->timeunit, &c);
      else if (wincon->means_dt == DAILY)
	wincon->means_next = windat->t + next_day(windat->t, 
						  wincon->timeunit, &c);
      else
	wincon->means_next = windat->t + wincon->means_dt;
    }

    /* Increment the time counter                                    */
    if (wincon->means & TIDAL) {
      ttime *= (0.8 * 3600.0);
      /* Average the elevation differences over the 2D loop          */
      for (cc = 1; cc <= window->b2_t; cc++) {
        c = window->w2_t[cc];
        time = windat->meanc[c];
        detadt[c] = (windat->eta[c] - wincon->oldeta[c]) / windat->dt;

        /* Find the start of the next tidal cycle                    */
        if (windat->meanc[c] == 0 && detadt[c] <= 0 &&
            windat->odeta[c] > 0)
          windat->meanc[c] += windat->dt;

        /* Increment and check for the tidal cycle end               */
        if (windat->meanc[0] >= 0 && windat->meanc[c] > 0) {

          /* Increment the iteration counter                         */
          windat->meanc[c] += windat->dt;

          /* Find when a tidal cycle is complete. When this occurs   */
          /* meanc contains the negative of the iterations performed */
          /* over the tidal cycle and meanc[0] is incremented. When  */
          /* meanc[0]=b2_t then all cells in the wet domain have     */
          /* undergone a full tidal cycle.                           */
          if (time > ttime && detadt[c] <= 0 && windat->odeta[c] > 0) {
            windat->meanc[c] *= -1;
            windat->meanc[0]++;
          }
        }
      }
      memcpy(windat->odeta, detadt, window->szcS * sizeof(double));
      /* When all cells have undergone a tidal cycle set meanc[0]<0  */
      /* and set the iterations performed over the cycle > 0.        */
      if (windat->meanc[0] == (double)window->b2_t) {
        windat->meanc[0] = -1;
      }
    } else {
      /*-------------------------------------------------------------*/
      /* Increment the iteration counter                             */
      windat->meanc[0] = 1.0;
      for (cc = 1; cc <= window->enonS; cc++) {
        c = window->wsa[cc];
        windat->meanc[cc] += windat->dtf;
      }
    }
  }
}

/* END reset_means()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Re-initializes means on the master. Legacy code as this is now    */
/* handed in the scheduler via means_event().                        */
/*-------------------------------------------------------------------*/
void reset_means_m_o(master_t *master)
{
  geometry_t *geom = master->geom;;
  int mon, day;
  double itc;
  int n;

  if (master->means & NONE) return;

  if (master->t >= master->means_next) {
    master->means |= RESET;  /* Set the RESET flag for transfers     */
    itc =  master->meanc[1]; /* Save the time counter                */
    /* Re-initialize the time counter.                               */
    memset(master->meanc, 0, geom->szcS * sizeof(double));
    /* Set the time counter flag; this is set to 0 when meanc is     */
    /* re-initialized in reset_means_3d_m() or reset_means_2d_m()    */
    /* for SEASONAL or MONTHLY means below.                          */
    master->meanc[0] = 1.0;  
    /* Re-initialize 3D mean tracers                                 */
    for (n = 0; n < master->ntm_3d; n++)
      reset_means3d_m(master, master->tm_3d[n]);
    if (master->means & VEL3D) memset(master->ume, 0, geom->sze * sizeof(double));
    if (master->means & VOLFLUX) memset(master->u1vm, 0, geom->sze * sizeof(double));
    /* Re-initialize 2D mean tracers                                 */
    for (n = 0; n < master->ntm_2d; n++)
      reset_means2d_m(master, master->tm_2d[n]);
    /* Set the next re-initialization time                           */
    if (master->means_dt == SEASONAL) {
      master->means_next = master->t + next_season(master->t, 
						   master->timeunit, &mon);
      master->meancs[mon] = itc;
    } else if (master->means_dt == MONTHLY) {
      master->means_next = master->t + next_month(master->t, 
						  master->timeunit, &mon);
      master->meancs[mon] = itc;
    } else if (master->means_dt == DAILY) {
      master->means_next = master->t + next_day(master->t, 
						  master->timeunit, &day);
      master->meancs[day] = itc;
    }
    else
      master->means_next = master->t + master->means_dt;
  }
}

/* END reset_means_m_o()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to zero 3d mean variables, or initialize from file.       */
/*-------------------------------------------------------------------*/
void reset_means3d_m(master_t *master, int trm)
{
  geometry_t *geom = master->geom;
  dump_data_t *dumpdata = master->dumpdata;
  double t;
  int cc, trd;
  int season, mon, sf;
  char *vname;

  vname = master->trinfo_3d[trm].name;
  if (master->means_dt == SEASONAL || master->means_dt == MONTHLY) {
    trd = tracer_find_index(vname, dumpdata->ntr, dumpdata->trinfo_3d);
    if (master->means_dt == SEASONAL) {
      t = prev_season(master->t, master->timeunit, &season, &mon);
      sf = SEASONAL;
    }
    if (master->means_dt == MONTHLY) {
      t = prev_month(master->t, master->timeunit, &mon);
      sf = MONTHLY;
    }
    if (master->means_dt == DAILY) {
      t = prev_day(master->t, master->timeunit, &mon);
      sf = DAILY;
    }
    if (read_mean_3d(master, t, sf, vname, master->tr_wc[trm], 
		     dumpdata->tr_wc[trd])) {
      if (master->meanc[0]) {
	for (cc = 1; cc <= geom->enonS; cc++)
	  master->meanc[cc] = master->meancs[mon];
	master->meanc[0] = 0.0;
      }
    }
  } else {
    memset(master->tr_wc[trm], 0, geom->szc * sizeof(double));
  }
}

/* END reset_means3d_m()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads in variable 'var' from a file contained in the dumplist     */
/* with SEASONAL|MONTHLY tinc, and places the data in v.             */
/*-------------------------------------------------------------------*/
int read_mean_3d(master_t *master, double t, int sf, char *var, 
		 double *v, double ***vc)
{
  geometry_t *geom = master->geom;
  dump_data_t *dumpdata = master->dumpdata;
  size_t nce1, nce2, nz, nMesh2_face, ndumps;
  int ti;
  int n, m, nn;
  size_t start[4];
  size_t count[4];
  int fid;
  double tvals;
  char timeunits[MAXSTRLEN];
  char incname[MAXSTRLEN];

  strcpy(incname, "SEASONAL");
  if (sf == MONTHLY)
    strcpy(incname, "MONTHLY");
  if (sf == DAILY)
    strcpy(incname, "DAILY");

  /*-----------------------------------------------------------------*/
  /* Loop through the dump files                                     */
  fid = -1;
  for (n = 0; n < dumpdata->ndf; ++n) {
    dump_file_t *dump = &dumpdata->dumplist[n];
    /* Find a file with SEASONAL|MONTHLY output                      */
    if (dump->incf == sf) {
      /* Check if 'var' is dumped in this file                       */
      for (m = 0; m < dump->nvars; m++) {
	if (strcmp(dump->vars[m], var) == 0) {
	  /* Open the file                                           */
	  if ((ncerr = nc_open(dump->name, NC_NOWRITE, &fid)) != NC_NOERR) {
	    hd_warn("Can't find %s dump file %s containing %s.\n", 
		    incname, dump->name, var);
	    hd_quit((char *)nc_strerror(ncerr));
	  }
	  /* Read and check dimensions                               */
	  if (nc_inq_dimlen(fid, ncw_dim_id(fid, "k_centre"), &nz) == 0 &&
	      nc_inq_dimlen(fid, ncw_dim_id(fid, "j_centre"), &nce2) == 0 &&
	      nc_inq_dimlen(fid, ncw_dim_id(fid, "i_centre"), &nce1) == 0) {
	    if (nce1 != geom->nce1 || nce2 != geom->nce2 || nz != geom->nz)
	      continue;
	    nc_inq_dimlen(fid, ncw_dim_id(fid, "record"), &ndumps);
	    /* Get time units                                        */
	    memset(timeunits, 0, MAXSTRLEN);
	    nc_get_att_text(fid, ncw_var_id(fid, "t"), "units", timeunits);
	    /* Find the dump number corresponding to t               */
	    ti = -1;
	    for (nn = 0; nn < ndumps; ++nn) {
	      size_t st = nn;
	      size_t co = 1;
	      nc_get_vara_double(fid, ncw_var_id(fid, "t"), &st, &co, &tvals);
	      tm_change_time_units(timeunits, master->timeunit, &tvals, 1);
	      if (fabs(tvals - t) < START_EPS) {
		ti = nn;
		break;
	      }
	    }
	    /* Dump file doesn't contain the time; zero and return   */
	    if (ti < 0) {
	      memset(v, 0, geom->szc * sizeof(double));
	      return(0);
	    }
	    /* Read in the data and map to the sparse variable       */
	    start[0] = ti;
	    start[1] = 0;
	    start[2] = 0;
	    start[3] = 0;
	    count[0] = 1L;
	    count[1] = nz;
	    count[2] = nce2;
	    count[3] = nce1;
	    nc_get_vara_double(fid, ncw_var_id(fid, var), start, count, vc[0][0]);
	    c2s_3d(geom, v, vc, geom->nce1, geom->nce2, geom->nz);
	    nc_close(fid);
	    return(1);
	  } else if (nc_inq_dimlen(fid, ncw_dim_id(fid, "Mesh2_layers"), &nz) == 0 &&
		     nc_inq_dimlen(fid, ncw_dim_id(fid, "nMesh2_face"), &nMesh2_face) == 0) {
	    int oset, cc, c, cs, k, **k2c;
	    if (nce1 != geom->b2_t || nz != geom->nz)
	      continue;
	    nc_get_att_int(fid, NC_GLOBAL, "start_index", &oset);
	    oset = (oset) ? 0 : 1;
	    nc_inq_dimlen(fid, ncw_dim_id(fid, "record"), &ndumps);
	    /* Get time units                                        */
	    memset(timeunits, 0, MAXSTRLEN);
	    nc_get_att_text(fid, ncw_var_id(fid, "t"), "units", timeunits);
	    /* Find the dump number corresponding to t               */
	    ti = -1;
	    for (nn = 0; nn < ndumps; ++nn) {
	      size_t st = nn;
	      size_t co = 1;
	      nc_get_vara_double(fid, ncw_var_id(fid, "t"), &st, &co, &tvals);
	      tm_change_time_units(timeunits, master->timeunit, &tvals, 1);
	      if (fabs(tvals - t) < START_EPS) {
		ti = nn;
		break;
	      }
	    }
	    /* Dump file doesn't contain the time; zero and return   */
	    if (ti < 0) {
	      memset(v, 0, geom->szc * sizeof(double));
	      return(0);
	    }
	    /* Make the reverse maps                                 */
	    k2c = i_alloc_2d(geom->szcS, geom->nz);
	    for (cc = 1; cc <= geom->b3_t; cc++) {
	      c = geom->w3_t[cc];
	      k = geom->s2k[c];
	      cs = geom->m2d[c];
	      k2c[k][cs-oset] = c;
	    }

	    /* Read in the data and map to the sparse variable       */
	    dumpdata_read_3d_us(dumpdata, fid, var, v, dumpdata->wc, ti, 
				nMesh2_face, geom->szc, geom->w3_t, k2c, oset);
	    nc_close(fid);
	    i_free_2d(k2c);
	    return(1);
	  }
	}
      }
    }
  }
  memset(v, 0, geom->szc * sizeof(double));
  return(0);
}

/* END reset_means_3d_m()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to zero 2d mean variables, or initialize from file.       */
/*-------------------------------------------------------------------*/
void reset_means2d_m(master_t *master, int trm)
{
  geometry_t *geom = master->geom;
  dump_data_t *dumpdata = master->dumpdata;
  double t;
  int cc, trd;
  int season, mon, sf;
  char *vname;

  vname = master->trinfo_2d[trm].name;
  if (master->means_dt == SEASONAL || master->means_dt == MONTHLY) {
    trd = tracer_find_index(vname, dumpdata->ntrS, dumpdata->trinfo_2d);
    if (master->means_dt == SEASONAL) {
      t = prev_season(master->t, master->timeunit, &season, &mon);
      sf = SEASONAL;
    }
    if (master->means_dt == MONTHLY) {
      t = prev_month(master->t, master->timeunit, &mon);
      sf = MONTHLY;
    }
    if (master->means_dt == DAILY) {
      t = prev_day(master->t, master->timeunit, &mon);
      sf = DAILY;
    }
    if (read_mean_2d(master, t, sf, vname, master->tr_wcS[trm], 
		     dumpdata->tr_wcS[trd])) {
      if (master->meanc[0]) {
	for (cc = 1; cc <= geom->enonS; cc++)
	  master->meanc[cc] = master->meancs[mon];
	master->meanc[0] = 0.0;
      }
    }
  } else
    memset(master->tr_wcS[trm], 0, geom->szcS * sizeof(double));
}

/* END reset_means2d_m()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads in variable 'var' from a file contained in the dumplist     */
/* with SEASONAL|MONTHLY tinc, and places the data in v.             */
/*-------------------------------------------------------------------*/
int read_mean_2d(master_t *master, double t, int sf, char *var, 
		 double *v, double **vc)
{
  geometry_t *geom = master->geom;
  dump_data_t *dumpdata = master->dumpdata;
  size_t nce1, nce2, nMesh2_face, ndumps;
  int ti;
  int n, m, nn;
  size_t start[4];
  size_t count[4];
  int fid;
  double tvals;
  char timeunits[MAXSTRLEN];
  char incname[MAXSTRLEN];

  strcpy(incname, "SEASONAL");
  if (sf == MONTHLY)
    strcpy(incname, "MONTHLY");
  if (sf == DAILY)
    strcpy(incname, "DAILY");

  /*-----------------------------------------------------------------*/
  /* Loop through the dump files                                     */
  fid = -1;
  for (n = 0; n < dumpdata->ndf; ++n) {
    dump_file_t *dump = &dumpdata->dumplist[n];
    /* Find a file with SEASONAL|MONTHLY output                      */
    if (dump->incf == sf) {
      /* Check if 'var' is dumped in this file                       */
      for (m = 0; m < dump->nvars; m++) {
	if (strcmp(dump->vars[m], var) == 0) {
	  /* Open the file                                           */
	  if ((ncerr = nc_open(dump->name, NC_NOWRITE, &fid)) != NC_NOERR) {
	    hd_warn("Can't find %s dump file %s containing %s.\n", 
		    incname, dump->name, var);
	    hd_quit((char *)nc_strerror(ncerr));
	  }
	  /* Read and check dimensions                               */
	  if (nc_inq_dimlen(fid, ncw_dim_id(fid, "j_centre"), &nce2) == 0 &&
	      nc_inq_dimlen(fid, ncw_dim_id(fid, "i_centre"), &nce1) == 0) {
	    if (nce1 != geom->nce1 || nce2 != geom->nce2)
	      continue;
	    nc_inq_dimlen(fid, ncw_dim_id(fid, "record"), &ndumps);
	    /* Get time units                                        */
	    memset(timeunits, 0, MAXSTRLEN);
	    nc_get_att_text(fid, ncw_var_id(fid, "t"), "units", timeunits);
	    /* Find the dump number corresponding to t               */
	    ti = -1;
	    for (nn = 0; nn < ndumps; ++nn) {
	      size_t st = nn;
	      size_t co = 1;
	      nc_get_vara_double(fid, ncw_var_id(fid, "t"), &st, &co, &tvals);
	      tm_change_time_units(timeunits, master->timeunit, &tvals, 1);
	      if (fabs(tvals - t) < START_EPS) {
		ti = nn;
		break;
	      }
	    }
	    /* Dump file doesn't contain the time; zero and return   */
	    if (ti < 0) {
	      memset(v, 0, geom->szcS * sizeof(double));
	      return(0);
	    }
	    /* Read in the data and map to the sparse variable       */
	    start[0] = ti;
	    start[1] = 0;
	    start[2] = 0;
	    start[3] = 0;
	    count[0] = 1L;
	    count[1] = nce2;
	    count[2] = nce1;
	    count[3] = 0;
	    nc_get_vara_double(fid, ncw_var_id(fid, var), start, count, vc[0]);
	    c2s_2d(geom, v, vc, geom->nce1, geom->nce2);
	    nc_close(fid);
	    return(1);
	  } else if (nc_inq_dimlen(fid, ncw_dim_id(fid, "nMesh2_face"), &nMesh2_face) == 0) {
	    int oset;
	    if (nce1 != geom->b2_t)
	      continue;
	    nc_get_att_int(fid, NC_GLOBAL, "start_index", &oset);
	    oset = (oset) ? 0 : 1;
	    nc_inq_dimlen(fid, ncw_dim_id(fid, "record"), &ndumps);
	    /* Get time units                                        */
	    memset(timeunits, 0, MAXSTRLEN);
	    nc_get_att_text(fid, ncw_var_id(fid, "t"), "units", timeunits);
	    /* Find the dump number corresponding to t               */
	    ti = -1;
	    for (nn = 0; nn < ndumps; ++nn) {
	      size_t st = nn;
	      size_t co = 1;
	      nc_get_vara_double(fid, ncw_var_id(fid, "t"), &st, &co, &tvals);
	      tm_change_time_units(timeunits, master->timeunit, &tvals, 1);
	      if (fabs(tvals - t) < START_EPS) {
		ti = nn;
		break;
	      }
	    }
	    /* Dump file doesn't contain the time; zero and return   */
	    if (ti < 0) {
	      memset(v, 0, geom->szc * sizeof(double));
	      return(0);
	    }

	    /* Read in the data and map to the sparse variable       */
	    dumpdata_read_2d_us(dumpdata, fid, var, v, ti, nMesh2_face, oset);
	    nc_close(fid);
	    return(1);
	  }
	}
      }
    }
  }
     memset(v, 0, geom->szcS * sizeof(double));
  return(0);
}

/* END read_mean_2d_m()                                              */
/*-------------------------------------------------------------------*/
