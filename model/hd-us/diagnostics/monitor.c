/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/diagnostics/monitor.c
 *  
 *  Description:
 *  Prints run time diagnostics to file
 *  from a file.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: monitor.c 7427 2023-10-25 01:26:35Z her127 $
 *
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include "hd.h"

#define START_EPS (0.001)
#define RHO_0 (1024.0)           /* standard density */
#define STAT_INCREMENT 10        /* Statistics output every seconds */

double fsw[10][11] =
  { {2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {10.0, 4.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {44.0, 15.0, -6.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {186.0, 56.0, -28.0, 9.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {772.0, 210.0, -120.0, 45.0, -10.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {3172.0, 792.0, -495.0, 220.0, -66.0, 12.0, -1.0, 0.0, 0.0, 0.0, 0.0},
    {12952.0, 3003.0, -2002.0, 1001.0, -364.0, -91.0, -14.0, 1.0, 0.0, 0.0, 0.0},
    {52666.0, 11440.0, -8008.0, 4368.0, -1820.0, 560.0, -120.0, 16.0, -1.0, 0.0, 0.0},
    {213324.0, 43758.0, -31824.0, 18564.0, -8568.0, 3060.0, -816.0, 153.0, -18.0, 1.0, 0.0},
    {863820.0, 167960.0, -125970.0, 77520.0, -38760.0, 15504.0, -4845.0, 1140.0, -190.0, 20.0, -1.0}
};


void wbdrycustom(parameters_t *params, FILE * fp, int bnum,
                 bdry_details_t *data);
void reset_cfl(master_t *master);
void reset_cfl_event(master_t *master, double mcfl, int mode);
double get_max_div_2d(geometry_t *window, window_t *windat,
		      win_priv_t *wincon, int *cl);
double get_max_div_3d(geometry_t *window, window_t *windat,
		      win_priv_t *wincon, int *cl);
double get_ddt(geometry_t *window, double *T, double *TB, int *cl);
double buoyancy_frequency2(geometry_t *window, window_t *windat,
			   win_priv_t *wincon, int c);
double calc_min_cfl_3d(geometry_t *window, window_t *windat,
		       win_priv_t *wincon, int *cl);
void remove_tide(geometry_t *window, double *eta);
void print_total_mass(master_t *master);
void psorts(double *a, int *c, int n);
void porders(double *p, double *q, int *i, int *j);
void sound_channel(geometry_t *window, window_t *windat, win_priv_t *wincon);
void quicksort(double *p, int *q, int l, int r);
int partition(double *p, int *q, int l, int r);
int mergeSort(int a[], int start, int end );
double vertex_mean(geometry_t *window, double *a, int c);
dump_file_t *create_output(master_t *master, int *fid, char *filename);
double percentiles[] = {0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40,
			0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85,
			0.90, 0.95, 1.00};
#define Ns (sizeof(percentiles) / sizeof(double))

/*-------------------------------------------------------------------*/
/* Routine to print simulation time details to file                  */
/*-------------------------------------------------------------------*/
void monitor(master_t *master,     /* Master data                    */
	     geometry_t **window,  /* Window geometry                */
             int mode           /* mode=0, initialise, mode=1 output */
	     )
{
  static int prev_elaps;
  double d1 = 0.0, d2, elapsed = 0.0;        /* Dummies              */
  int n;
  FILE *dfp, *fopen();

  if (autof == 9) return;

  if (mode == 0) {
    master->iclk = dp_clock();
    master->tvol = master->tmass = 0.0;
    memset(master->trtot, 0, master->ntot * sizeof(double));
    if (master->vf)
      memset(master->vf, 0, master->geom->nobc * sizeof(double));
  } else if (mode == -1) {
    master->iclk = 0.0;
    master->mclk = 0.0;
    master->mrtr = 0.0;
    master->mcfl2d = master->mcfl3d = HUGE;
  } else {
    if (master->totals)
      print_total_mass(master);
    if (master->errornorm & N_ERRN)
      print_error_norms(master);
    if (!(master->decf & NONE)) {
      int dim = (master->decf == DEC_ETA) ? 2 : 3;
      calc_decorr(geom, master->decv, master->decv1, master->decorr, 
		  dim, master->decs);
    }
    if (diag_log) {
      dfp = fopen(diag_logfile, "w");
      if (dfp != NULL) {
        fprintf(dfp, "\nSimulation start = %9.4f (days) : %s\n",
                schedule->start_time / 86400.0,
		tm_time_to_datestr(schedule->start_time, master->timeunit));
        fprintf(dfp, "Simulation stop  = %9.4f (days) : %s\n",
                schedule->stop_time / 86400.0,
		tm_time_to_datestr(schedule->stop_time, master->timeunit));
        fprintf(dfp, "Simulation time  = %9.4f (days) : %s\n\n",
                master->t / 86400.0, 
		tm_time_to_datestr(master->t, master->timeunit));

	/* DA stuff */
	if (master->da == NO_DA) 
	  fprintf(dfp, "Data Assimilation : forecast cycle\n\n");
	if (master->da == DO_DA)
	  fprintf(dfp, "Data Assimilation : reanalysis cycle\n\n");
	if (master->da & FCST_DA)
	  fprintf(dfp, "Data Assimilation : extended forecast cycle\n\n");

        master->iclk = dp_clock() - master->iclk;
        if (master->iclk < 0) {
          fprintf(dfp, "Negative CPU time = %5.3f (sec)\n", master->iclk);
          master->iclk = master->mclk / master->nstep;
        }
        master->mclk += master->iclk;
        if (master->nstep != 0)
          d1 = master->mclk / master->nstep;
        fprintf(dfp, "CPU time used this iteration = %5.3f (sec)\n",
                master->iclk);
        fprintf(dfp, "Mean CPU time used / iteration = %5.3f (sec)\n", d1);
        fprintf(dfp, "CPU run time ratio = %f\n", master->dt / d1);

        if (master->nwindows > 1) {
          d2 = 0.0;
          for (n = 1; n <= master->nwindows; n++)
            d2 += master->wclk[n];
          for (n = 1; n <= master->nwindows; n++)
            fprintf(dfp, "Window %d integration time = %3.1f %%\n", n,
                    100 * master->wclk[n] / d2);
          /* Get the new window sizes                                */
          if (master->geom->win_size) {
            geometry_t *geom = master->geom;
            double d3 = 0.0;
            double size[geom->nwindows+1];
            for (n = 1; n <= master->nwindows; n++) {
              size[n] = geom->win_size[n] *
                (d2 / (master->wclk[n] * master->nwindows));
              d3 += size[n];
            }
            fprintf(dfp, "   Actual load balance = ");
            for (n = 1; n <= master->nwindows; n++) {
              fprintf(dfp, "%4.3f ", geom->win_size[n]);
            }
            fprintf(dfp, "\nPredicted load balance = ");
            for (n = 1; n <= master->nwindows; n++) {
              size[n] /= d3;
              fprintf(dfp, "%4.3f ", size[n]);
            }
            fprintf(dfp, "\n");
          }
        }

        /* Elapsed time                                              */
        elapsed = difftime(time(NULL), schedule->exec_start_time) / 3600.0;
        fprintf(dfp, "Elapsed time = %d day(s) %2.2d:%2.2d:%2.2d\n",
              (int)(elapsed/24.0), (int)fmod(elapsed, 24.0), (int)fmod(elapsed*60, 60.0),
              (int)fmod(elapsed*3600, 60.0));
        d1 = ((master->t - schedule->start_time)/3600) / elapsed;
        d2 = (schedule->stop_time - master->t)/ 3600.0;
	if (master->nstep > 1 && d1 < HUGE)
	  master->mrtr = (master->mrtr * (double)(master->nstep-1) + d1) / (double)master->nstep;
        fprintf(dfp, "Total time ratio = %.2f (%d:1)\n", d1, (int)master->mrtr);
        /* Estimated time to completion                              */
	if (d1 > 0)
	  d2 = d2 /d1;
        fprintf(dfp, "Time to completion = %d day(s) %2.2d:%2.2d:%2.2d\n",
              (int)(d2/24.0), (int)fmod(d2, 24.0), (int)fmod(d2*60, 60.0),
              (int)fmod(d2*3600, 60.0));
        /* Percent complete                                          */
        fprintf(dfp, "Percent complete = %2.1f%%\n",
            (master->t - schedule->start_time)*100.0
              / (schedule->stop_time - schedule->start_time));

        if (!(master->cfl & NONE) && !(master->mode2d)) {
          fprintf(dfp, "CFL Min: 2D=%4.1f 3D=%4.1f\n", master->mcfl2d,
                  master->mcfl3d);
          if (master->cfl & ACTIVE)
            reset_cfl(master);
          if (master->cfl & ACTIVE3D)
            reset_cfl_event(master, master->mcfl3d, CFL_DIAG);
        }
	if (master->magec > 0.0) {
          fprintf(dfp, "Average age of %d particles = %4.1f (days)\n",
		  (int)master->magec, master->mage /
		  (86400.0 * master->magec));
	}
#if defined(HAVE_ECOLOGY_MODULE)
	if (DEBUG("ecology")) {
	  for(n = 1; n <= master->nwindows; n++) {
	    win_priv_t *wincon = window[n]->wincon;
	    if(wincon->do_eco) {
	      fprintf(dfp, "ecology : window %d\n", n);
	      ecology_printstepstats(wincon->e, dfp);
	    }
	  }
	}
#endif
	if (mode == 3) {
	  fprintf(dfp, "** CRASHED **\n");
	  history_log(master, HST_NOK);
	} else {
	  if (master->t == schedule->stop_time) {
	    fprintf(dfp, "Run successful.\n");
	    history_log(master, HST_OK);
	  } else
	    fprintf(dfp, "Running...\n");
	}
        n = fclose(dfp);
	FLUSH_TIMING
      }
    }
    /* Runtime statistics log, this should only be enabled during    */
    /* development.                                                    */
    if(stat_log) {
      int elaps = (int)(elapsed*3600);
      if(master->t - master->dt - schedule->start_time < 1) {
        dfp = fopen(stat_logfile, "w");
        if (dfp != NULL)
            fprintf(dfp, "Run time (sec) \t Simulation time (sec) \t complete \n");

      } else if(elaps % STAT_INCREMENT == 0 && elaps > prev_elaps) /* only every 10 seconds */
        dfp = fopen(stat_logfile, "a");
      else
        dfp = NULL;

      if (dfp != NULL) {
        fprintf(dfp, " %d %d %3.1f\n",
                            elaps,
                            (int)(master->t - schedule->start_time),
                            ((master->t - schedule->start_time)*100.0 / (schedule->stop_time - schedule->start_time)));
        fclose(dfp);
        prev_elaps = elaps;
      }
    }
  }
}

/* END monitor()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* print debugging information for a nominated location              */
/*-------------------------------------------------------------------*/
void debug_c(geometry_t *window, int var, int mode)
{
  window_t *windat = window->windat;    /* Window data               */
  win_priv_t *wincon = window->wincon;  /* Window constants          */
  int wn, e, es, cc, c, cs, cm, j, jd;
  int ee, e1, e1s;
  double max;

  if (wincon->dbc == 0) return;
  if (window->wn != wincon->dbw) return;
  if (wincon->dbgf & D_AFTER && windat->t < wincon->dbgtime) return;
  c = wincon->dbc;
  jd = wincon->dbj;
  e1 = window->c2e[jd][c];
  cs = window->m2d[c];
  e1s = window->m2de[e1];
  cc = window->c2cc[c];
  wn = window->wn;

  if (var == D_INIT) {
    if (!(wincon->dbgf & D_OPEN)) {
      wincon->dbf = fopen("debug.txt", "w");
      wincon->dbgf |= D_OPEN;
    }
    if (mode == 0 && wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "Start HD step\n");
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "\nDebugging for %d = (%d %d %d), time %f\n\n",
	   c, window->s2i[c], window->s2j[c], window->s2k[c], windat->days);
    fprintf(wincon->dbf, "nstep = %d\n", windat->nstep);
    fprintf(wincon->dbf, "Window number = %d\n", wn);
    fprintf(wincon->dbf, "edges = %d\n", window->npe[cs]);
    fprintf(wincon->dbf, "cc = %d\n", cc);
    fprintf(wincon->dbf, "cs = %d\n", cs);
    fprintf(wincon->dbf, "cb = %d\n", window->bot_t[cc]);
    fprintf(wincon->dbf, "cg = %d\n", window->wsa[c]);
    for (j = 1; j <= window->npe[cs]; j++) {
      e = window->c2e[j][c];
      es = window->m2de[e];
      fprintf(wincon->dbf, "e = %d : e1s = %d, dir=%d\n", e, es, j);
    }
    for (j = 1; j <= window->npe[cs]; j++) {
      fprintf(wincon->dbf, "c2c[%d] = %d (cg=%d)\n", j, window->c2c[j][c],
	      window->wsa[window->c2c[j][c]]);
    }
    fprintf(wincon->dbf, "depth = %5.2f\n", window->botz[cs]);
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = wincon->s1[cc];
      if (cs == window->m2d[c]) break;
    }
    if (c == 0)c = cs;
    e1 = window->c2e[1][c];
    e1s = window->m2de[e];
    fprintf(wincon->dbf,"\nc       k   eta    u1av   u1     temp  salt   dz    dzu1  gridz\n");
    fprintf(wincon->dbf,"%-7d %-3d %-5.3f  %-6.3f %-6.3f %5.2f %5.2f %5.2f %5.2f %5.2f\n",
	    c, window->s2k[c], windat->eta[cs], windat->u1av[e1s], windat->u1[e1], 
	    windat->temp[c], windat->sal[c], wincon->dz[c], windat->dzu1[e1], 
	    window->gridz[c]);    
    c = window->zm1[c];
    if (c != window->zm1[c]) {
      do {
	fprintf(wincon->dbf,"%-7d %-3d               %-6.3f %5.2f %5.2f %5.2f %5.2f %5.2f\n",
		c, window->s2k[c], windat->u1[e1], windat->temp[c], 
		windat->sal[c], wincon->dz[c], windat->dzu1[e1],
		window->gridz[c]);
	c = window->zm1[c];
      } while (c != window->zm1[c]);
    }
    fflush(wincon->dbf);
  }

  if (var == D_U && mode == D_ADVECT) {
    fprintf(wincon->dbf, "\n3D mode\n");
    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "u1 advection\n");
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "Post-advection U1 = %f\n", windat->nu1[e1]);
    fprintf(wincon->dbf, "  u1inter(3D)  = %f\n", wincon->u1inter[e1s]);
  }
  if (var == D_U && mode == D_HDIFF) {
    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "u1 horizontal diffusion\n");
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "Post-diffusion U1 = %f\n", windat->nu1[e1]);
    fprintf(wincon->dbf, "  u1inter (+2D) = %f\n", wincon->u1inter[e1s]);
  }
  if (var == D_U && mode == D_PRESSURE) {
    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "u1 pressure\n");
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "Post-pressure term U1 = %f\n", windat->nu1[e1]);
    fprintf(wincon->dbf, "  u1inter (+integral rho) = %f\n", wincon->u1inter[e1s]);
  }
  if (var == D_U && mode == D_VZ) {
    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "u1 vertical mixing\n");
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "Post-vertical mixing U1 = %f\n", windat->nu1[e1]);
    fprintf(wincon->dbf, "  u1inter (+wind stress) = %f\n", wincon->u1inter[e1s]);
  }
  if (var == D_U && mode == D_POST) {
    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "u1 complete\n");
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "New U1 = %f\n", windat->nu1[e1]);
    fprintf(wincon->dbf, "  u1inter = %f\n\n", wincon->u1inter[e1s]);
  }


  if (var == D_UA && mode == D_ADVECT) {
    fprintf(wincon->dbf, "\n2D mode step %d\n", wincon->ic);
    fprintf(wincon->dbf, "---------------\n");
    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "u1av advection\n");
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "Post-advection U1AV = %f\n", windat->nu1av[e1s]);
  }
  if (var == D_UA && mode == D_HDIFF) {
    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "u1av horizontal diffusion\n");
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "Post-diffusion U1AV = %f\n", windat->nu1av[e1s] + 
	    wincon->tend2d[T_HDF][e1s]);
  }
  if (var == D_UA && mode == D_POST) {
    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "u1av update\n");
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "Barotropic pressure U1AV = %e\n", wincon->b1);
    fprintf(wincon->dbf, "Coriolis U1AV = %e\n", wincon->b2);
    fprintf(wincon->dbf, "Bottom friction U1AV = %e\n", wincon->b3);
    fprintf(wincon->dbf, "Wind stress/baroclinic pressure U1AV = %e\n", wincon->u1inter[e1s]);
    fprintf(wincon->dbf, "U1AV + tendencies = %f\n", windat->nu1av[e1s]);
  }
  if (var == D_UA && mode == D_BDRY) {
    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "u1av complete\n");
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "New U1AV = %f\n", windat->nu1av[e1s]);
    max = 0;
    cm = 1;
    for (ee = 1; ee <= window->b2_e1; ee++) {
      e = window->w2_e1[ee];
      if (fabs(windat->nu1av[e]) > max) {
	max = fabs(windat->nu1av[e]);
	cm = window->e2ijk[e];
	e1 = e;
      }
    }
    for (ee = 1; ee <= window->npe[cm]; ee++) {
      if (window->c2e[ee][cm] == e1) break;
    }
    fprintf(wincon->dbf, "Maximum |u1av| in window%d = %f @ (e=%d ee=%d c=%d cg=%d)\n\n",
	    window->wn, max, e1, ee, window->s2i[cm], window->s2j[cm]);
    fflush(wincon->dbf);
  }


  if (var == D_ETA && mode == D_POST) {
    double colflux;
    double neweta;

    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "eta update\n");
      fflush(wincon->dbf);
      return;
    }
    colflux = 0.0;
    for (j = 1; j <= window->npe[cs]; j++) {
      fprintf(wincon->dbf, "edge%d u1av = %-6.3f : flux = %6.3e\n", j, windat->u1av[window->c2e[j][c]],wincon->ba[j]);
      colflux += wincon->ba[j];
    }
    fprintf(wincon->dbf, "colflux = %f\n", colflux);
    neweta = max(windat->etab[c] - colflux / window->cellarea[c],
		 window->botz[c]);
    fprintf(wincon->dbf, "eta(t+1) = %f\n\n", neweta);
    fflush(wincon->dbf);
  }

  if (var == D_TS && mode == D_INIT) {
    fprintf(wincon->dbf, "Start tracer step\n");
    if (wincon->dbgf & D_STEP) {
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "  Temperature = %f\n", windat->temp[c]);
    fprintf(wincon->dbf, "  Salinity = %f\n", windat->sal[c]);
  }
  if (var == D_TS && mode == D_ADVECT) {
    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "Tracer advection complete\n");
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "After tracer advection\n");
    fprintf(wincon->dbf, "  Temperature = %f\n", windat->temp[c]);
    fprintf(wincon->dbf, "  Salinity = %f\n", windat->sal[c]);
  }
  if (var == D_TS && mode == D_STRML) {
    cm = wincon->s2[wincon->s4[c]];
    fprintf(wincon->dbf, "Source cell for %d = %d (%d %d %d)\n",c, cm,
	    window->s2i[cm], window->s2j[cm], window->s2k[cm]);	    
  }
  if (var == D_TS && mode == D_HDIFF) {
    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "Tracer vertical diffusion complete\n");
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "After tracer vertical diffusion\n");
    fprintf(wincon->dbf, "  Temperature = %f\n", windat->temp[c]);
    fprintf(wincon->dbf, "  Salinity = %f\n", windat->sal[c]);
  }
  if (var == D_TS && mode == D_BDRY) {
    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "Tracer OBCs complete\n");
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "After tracer OBCs\n");
    fprintf(wincon->dbf, "  Temperature = %f\n", windat->temp[c]);
    fprintf(wincon->dbf, "  Salinity = %f\n", windat->sal[c]);
  }
  if (var == D_TS && mode == D_POST && wincon->dbgf & D_OPEN) {
    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "Tracer update complete\n");
      fflush(wincon->dbf);
    } else {
      fprintf(wincon->dbf, "Temperature = %f\n", windat->temp[c]);
      fprintf(wincon->dbf, "Salinity = %f\n", windat->sal[c]);
      fprintf(wincon->dbf, "Density = %f\n", windat->dens[c]);
    }
    fprintf(wincon->dbf,"\n-----------------------------------------------\n");
    fprintf(wincon->dbf,"-----------------------------------------------\n\n\n");
    if (!(wincon->dbgf & D_APPEND)) {
      fclose(wincon->dbf);
      wincon->dbgf &= ~D_OPEN;
    } else
      fflush(wincon->dbf);
  }
}

/* END debug_c()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* print debugging information for model crashes                     */
/*-------------------------------------------------------------------*/
void crash_c(geometry_t *window, int c, int e, char *text)
{
  FILE *fp;
  window_t *windat = window->windat;    /* Window data               */
  win_priv_t *wincon = window->wincon;  /* Window constants          */
  int wn, es, cc, cs, cm, j;
  int ee, e1, e1s;
  int eed = 0;
  double max;

  if (c == 0) c = window->e2c[e][0];
  cs = window->m2d[c];
  cc = window->c2cc[c];
  wn = window->wn;
  if (e == 0) {
    max = 0.0;
    for (ee = 1; ee <= window->npe[cs]; ee++) {
      e1 = window->c2e[ee][cs];
      if (isnan(windat->u1av[e1]) || fabs(windat->u1av[e1]) > max) {
	e = e1;
	eed = ee;
	max = fabs(windat->u1av[e1]);
	if (isnan(windat->u1av[e1])) break;
      }
    }
  }
  es = window->m2de[e];

  fp = fopen("crash.txt", "w");
  fprintf(fp, "Crash diagnostics : %s instability\n\n", text);
  fprintf(fp, "Cell location %d = (%d %d %d), time %f\n\n",
	  c, window->s2i[c], window->s2j[c], window->s2k[c], windat->days);
  fprintf(fp, "nstep = %d\n", windat->nstep);
  fprintf(fp, "Window number = %d\n", wn);
  fprintf(fp, "edges = %d\n", window->npe[cs]);
  fprintf(fp, "cc = %d\n", cc);
  fprintf(fp, "cs = %d\n", cs);
  fprintf(fp, "cb = %d\n", window->bot_t[cc]);
  fprintf(fp, "cg = %d\n", window->wsa[c]);
  fprintf(fp, "Cell location = %f %f\n", window->cellx[cs], window->celly[cs]);
  fprintf(fp, "Edge location = %f %f\n", window->u1x[es], window->u1y[es]);
  if (eed)
    fprintf(fp, "e = %d, eg = %d, ee = %d\n", e, window->wse[e], eed);
  else
    fprintf(fp, "e = %d, eg = %d\n", e, window->wse[e]);
  for (j = 1; j <= window->npe[cs]; j++) {
    e = window->c2e[j][c];
    es = window->m2de[e];
    fprintf(fp, "e = %d : e1s = %d, dir=%d\n", e, es, j);
  }
  for (j = 1; j <= window->npe[cs]; j++) {
    fprintf(fp, "c2c[%d] = %d\n", j, window->c2c[j][c]);
  }
  fprintf(fp, "depth = %5.2f\n", window->botz[cs]);
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = wincon->s1[cc];
    if (cs == window->m2d[c]) break;
  }
  if (c == 0)c = cs;
  fprintf(fp,"\nc       k   eta    u1av   u1     temp  salt   dz    dzu1  gridz\n");
  fprintf(fp,"%-7d %-3d %-5.3f  %-6.3f %-6.3f %5.2f %5.2f %5.2f %5.2f %5.2f\n",
	  c, window->s2k[c], windat->eta[cs], windat->u1av[es], windat->u1[e], 
	  windat->temp[c], windat->sal[c], wincon->dz[c], windat->dzu1[e], 
	  window->gridz[c]);    
  c = window->zm1[c];
  if (c != window->zm1[c]) {
    do {
      fprintf(fp,"%-7d %-3d               %-6.3f %5.2f %5.2f %5.2f %5.2f %5.2f\n",
	      c, window->s2k[c], windat->u1[e], windat->temp[c], 
	      windat->sal[c], wincon->dz[c], windat->dzu1[e],
	      window->gridz[c]);
      c = window->zm1[c];
    } while (c != window->zm1[c]);
  }

  fprintf(fp, "\n3D velocity tendencies\n");
  fprintf(fp, "Advection tendency = %f\n", wincon->tend3d[T_ADV][e]);
  fprintf(fp, "Horizontal diffusion tendency = %f\n", wincon->tend3d[T_HDF][e]);
  fprintf(fp, "Vertical diffusion tendency = %f\n", wincon->tend3d[T_VDF][e]);
  fprintf(fp, "Barotropic pressure tendency = %f\n", wincon->tend3d[T_BTP][e]);
  fprintf(fp, "Baroclinic pressure tendency = %f\n", wincon->tend3d[T_BCP][e]);
  if(wincon->waves & STOKES_DRIFT)
    fprintf(fp, "Stokes tendency = %f\n", wincon->tend3d[T_STK][e]);

  fprintf(fp, "\n2D velocity tendencies\n");
  fprintf(fp, "Advection tendency = %f\n", wincon->tend2d[T_ADV][es]);
  fprintf(fp, "Horizontal diffusion tendency = %f\n", wincon->tend2d[T_HDF][es]);
  fprintf(fp, "Barotropic pressure tendency = %f\n", wincon->tend2d[T_BTP][es]);
  fprintf(fp, "Bottom stress tendency = %f\n", wincon->tend2d[T_BOT][es]);
  fprintf(fp, "3D contributions (wind, density) tendency = %f\n", wincon->u1inter[es]);
  fprintf(fp, "Non-linear contributions = %f\n", wincon->tend2d[T_NLI][es]);

  fprintf(fp, "\n2D velocity\n");
  for (ee = 1; ee <= window->npe[cs]; ee++) {
    e = window->c2e[ee][cs];
    fprintf(fp, "%d %d %f\n",ee, e, windat->u1av[e]);
  }
  fclose(fp);
}

/* END crash_c()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to get the tendency of a particular process and store in  */
/* an array.                                                         */
/*-------------------------------------------------------------------*/
void get_tend(geometry_t *window, /* Window geometry                 */
              int *vec,           /* Vector of cells to process      */
              int eb,             /* End index for vec               */
              double *vel,        /* Velocity array after process    */
              double *ovel,       /* Velocity array before process   */
              double *tendency    /* Tendency of process             */
  )
{
  int c, cs, cc;                  /* Sparse coordinate/counter       */
  double t, tend;
  win_priv_t *wincon = window->wincon;
  window_t *windat = window->windat;

  /* Set the tendency if required                                    */
  if (tendency) {
    if (!(wincon->means & TENDENCY)) {
      memset(tendency, 0, window->sze * sizeof(double));
      for (cc = 1; cc <= eb; cc++) {
        c = vec[cc];
        tendency[c] += (vel[c] - ovel[c]);
      }
    } else {
      t = windat->dtf;
      for (cc = 1; cc <= eb; cc++) {
        c = vec[cc];
	cs = window->m2d[c];
	tend = (vel[c] - ovel[c]);
	if (windat->meanc[cs] > 0)
	  tendency[c] = (tendency[c] * windat->meanc[cs] + tend * t) / 
	    (windat->meanc[cs] + t);
      }
    }
  }
  /* Reset the velocity after this process                           */
  memcpy(ovel, vel, window->sze * sizeof(double));
}

void get_tendv(geometry_t *window, /* Window geometry                */
	       int *vec,           /* Vector of cells to process     */
	       int eb,             /* End index for vec              */
	       double *tend,       /* Tendency                       */
	       double *tend1,      /* East tendency of process       */
	       double *tend2       /* North tendency of process      */
  )
{
  int c, cs, cc, e, es, ee;       /* Counters                        */
  double t, area;
  win_priv_t *wincon = window->wincon;
  window_t *windat = window->windat;

  /* Set the tendency if required                                    */
  if (tend1 && tend2) {
    if (!(wincon->means & TENDENCY)) {
      /* Rotate into east and north cell centered components         */
      vel_cen(window, windat, wincon, tend, NULL, tend1, tend2, NULL, NULL, 0);
    } else {
      double *ut = wincon->w4;
      double *vt = wincon->w5;
      vel_cen(window, windat, wincon, tend, NULL, ut, vt, NULL, NULL, 0);
      t = windat->dtf;
      for (cc = 1; cc <= window->b3_t; cc++) {
        c = window->w3_t[cc];
	cs = window->m2d[c];

	/* Update the cell centered mean                             */
	if (windat->meanc[cs] > 0)
	  tend1[c] = (tend1[c] * windat->meanc[cs] + ut[c] * t) / 
	    (windat->meanc[cs] + t);
	  tend2[c] = (tend2[c] * windat->meanc[cs] + vt[c] * t) / 
	    (windat->meanc[cs] + t);
      }
    }
  }
}

/* END get_tend()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Initialise the alert output file and alert variables.             */
/*-------------------------------------------------------------------*/
void alerts_init(master_t *master,     /* Master data                */
		 geometry_t **window   /* Window geometry            */
  )
{
  geometry_t *geom = master->geom;     /* Master geometry            */
  int nb, c;                           /* Counters                   */

  if (!(master->alertf & NONE)) {
    /* Initialise the maximums                                       */
    master->mcfl = master->scfl = HUGE;
    master->ncfl = -1;
    master->meta = 0.0;
    master->mu1 = 0.0; master->mu1a = 0.0;
    master->mu2 = 0.0; master->mu2a = 0.0;
    master->mdeta = 0.0; master->mdw = 0.0;
    master->mdt = master->mds = 0.0;
    master->msmin = master->mtmin = HUGE;
    master->msmax = master->mtmax = 0.0;
    master->mshear = 0.0;
    master->ma1 = master->ma2 = 0.0;
    master->mh1 = master->mh2 = 0.0;
    master->mv1 = master->mv2 = 0.0;
    master->mb1 = master->mb2 = 0.0;
    master->md1 = master->md2 = 0.0;
    master->mc1 = master->mc2 = 0.0;
    if (geom->nobc)
      master->ef = d_alloc_1d(geom->nobc);
    master->nalert = i_alloc_1d(nna);
    memset(master->nalert, 0, nna * sizeof(int));
    /* Window alerts                                                 */
    for (nb = 1; nb <= master->nwindows; nb++) {
      window_t *windat = window[nb]->windat;
      if (window[nb]->nobc)
	windat->ef = d_alloc_1d(window[nb]->nobc);
      windat->nalert = i_alloc_1d(nna);
      memset(windat->nalert, 0, nna * sizeof(int));
      windat->meta = 0.0;
      windat->mu1 = 0.0; windat->mu1a = 0.0;
      windat->mu2 = 0.0; windat->mu2a = 0.0;
      windat->mdeta = 0.0; windat->mdw = 0.0;
      windat->ma1 = 0.0; windat->mh1 = 0.0; windat->mv1 = 0.0;
      windat->mb1 = 0.0; windat->md1 = 0.0; windat->mc1 = 0.0;
      windat->ma2 = 0.0; windat->mh2 = 0.0; windat->mv2 = 0.0;
      windat->mb2 = 0.0; windat->md2 = 0.0; windat->mc2 = 0.0;
      windat->mdt = windat->mds = 0.0;
      windat->msmin = windat->mtmin = HUGE;
      windat->msmax = windat->mtmax = 0.0;
      windat->mshear = 0.0;
      windat->ceta = i_alloc_1d(window[nb]->szcS);
      windat->cu1 = i_alloc_1d(window[nb]->sze);
      windat->cu1a = i_alloc_1d(window[nb]->szeS);
      windat->tempb = d_alloc_1d(window[nb]->szc);
      windat->salb = d_alloc_1d(window[nb]->szc);
    }
    if (master->alert_dt) {
      tsalert.fp = NULL;
      if (strlen(master->opath))
	sprintf(tsalert.pname,"%s%s.ts", master->opath, master->alertname);
      else
	sprintf(tsalert.pname,"%s.ts",master->alertname);
      tsalert.fp = fopen(tsalert.pname, "w");
      tsalert.master = master;
      tsalert.i = 0;
      tsalert.j = 0;
      tsalert.k = 0;
      tsalert.tsout = master->t;
      tsalert.tsdt = master->alert_dt;
      c = 15 + geom->nobc;
      fprintf(tsalert.fp, "## COLUMNS %d\n", 1 + c);
      fprintf(tsalert.fp, "##\n");
      fprintf(tsalert.fp, "## COLUMN1.name  Time\n");
      fprintf(tsalert.fp, "## COLUMN1.long_name  Time\n");
      fprintf(tsalert.fp, "## COLUMN1.units  %s\n", master->output_tunit);
      fprintf(tsalert.fp, "##\n");
      c = 1;
      fprintf(tsalert.fp, "## COLUMN%1.1d.name  eta\n", ++c);
      fprintf(tsalert.fp, "## COLUMN%1.1d.long_name Maximum surface elevation\n",c);
      fprintf(tsalert.fp, "## COLUMN%1.1d.units m\n",c);
      fprintf(tsalert.fp, "##\n");
      fprintf(tsalert.fp, "## COLUMN%1.1d.name  u1\n", ++c);
      fprintf(tsalert.fp, "## COLUMN%1.1d.long_name Maximum u1 velocity\n",c);
      fprintf(tsalert.fp, "## COLUMN%1.1d.units ms-1\n",c);
      fprintf(tsalert.fp, "##\n");
      fprintf(tsalert.fp, "## COLUMN%1.1d.name  u1av\n", ++c);
      fprintf(tsalert.fp, "## COLUMN%1.1d.long_name Maximum u1av velocity\n",c);
      fprintf(tsalert.fp, "## COLUMN%1.1d.units ms-1\n",c);
      fprintf(tsalert.fp, "##\n");
      fprintf(tsalert.fp, "## COLUMN%1.1d.name  u2\n", ++c);
      fprintf(tsalert.fp, "## COLUMN%1.1d.long_name Maximum u2 velocity\n",c);
      fprintf(tsalert.fp, "## COLUMN%1.1d.units ms-1\n",c);
      fprintf(tsalert.fp, "##\n");
      fprintf(tsalert.fp, "## COLUMN%1.1d.name  u2av\n", ++c);
      fprintf(tsalert.fp, "## COLUMN%1.1d.long_name Maximum u2av velocity\n",c);
      fprintf(tsalert.fp, "## COLUMN%1.1d.units ms-1\n",c);
      fprintf(tsalert.fp, "##\n");
      fprintf(tsalert.fp, "## COLUMN%1.1d.name  w\n", ++c);
      fprintf(tsalert.fp, "## COLUMN%1.1d.long_name Maximum vertical velocity\n",c);
      fprintf(tsalert.fp, "## COLUMN%1.1d.units ms-1\n",c);
      fprintf(tsalert.fp, "##\n");
      fprintf(tsalert.fp, "## COLUMN%1.1d.name  detadt\n", ++c);
      fprintf(tsalert.fp, "## COLUMN%1.1d.long_name Maximum rate of change of eta\n",c);
      fprintf(tsalert.fp, "## COLUMN%1.1d.units ms-1\n",c);
      fprintf(tsalert.fp, "##\n");
      fprintf(tsalert.fp, "## COLUMN%1.1d.name  dwdz\n", ++c);
      fprintf(tsalert.fp, "## COLUMN%1.1d.long_name Maximum vertical velocity gradient\n",c);
      fprintf(tsalert.fp, "## COLUMN%1.1d.units s-1\n",c);
      fprintf(tsalert.fp, "##\n");
      fprintf(tsalert.fp, "## COLUMN%1.1d.name  shear\n", ++c);
      fprintf(tsalert.fp, "## COLUMN%1.1d.long_name Maximum horizontal shear\n",c);
      fprintf(tsalert.fp, "## COLUMN%1.1d.units s-1\n",c);
      fprintf(tsalert.fp, "##\n");
      fprintf(tsalert.fp, "## COLUMN%1.1d.name  cfl\n", ++c);
      fprintf(tsalert.fp, "## COLUMN%1.1d.long_name Minimum 3D cfl\n",c);
      fprintf(tsalert.fp, "## COLUMN%1.1d.units s\n",c);
      fprintf(tsalert.fp, "##\n");
      fprintf(tsalert.fp, "## COLUMN%1.1d.name  dT/dt\n", ++c);
      fprintf(tsalert.fp, "## COLUMN%1.1d.long_name Maximum rate of change of T\n",c);
      fprintf(tsalert.fp, "## COLUMN%1.1d.units oCs-1\n",c);
      fprintf(tsalert.fp, "##\n");
      fprintf(tsalert.fp, "## COLUMN%1.1d.name  dS/dt\n", ++c);
      fprintf(tsalert.fp, "## COLUMN%1.1d.long_name Maximum rate of change of S\n",c);
      fprintf(tsalert.fp, "## COLUMN%1.1d.units psus-1\n",c);
      fprintf(tsalert.fp, "##\n");
      fprintf(tsalert.fp, "## COLUMN%1.1d.name  MechanicalEnergy\n", ++c);
      fprintf(tsalert.fp, "## COLUMN%1.1d.long_name Mechanical energy\n",c);
      fprintf(tsalert.fp, "## COLUMN%1.1d.units Jm-2\n",c);
      fprintf(tsalert.fp, "##\n");
      fprintf(tsalert.fp, "## COLUMN%1.1d.name  KineticEnergy\n", ++c);
      fprintf(tsalert.fp, "## COLUMN%1.1d.long_name Kinetic energy\n",c);
      fprintf(tsalert.fp, "## COLUMN%1.1d.units Jm-3\n",c);
      fprintf(tsalert.fp, "##\n");
      fprintf(tsalert.fp, "## COLUMN%1.1d.name  ExcessMass\n", ++c);
      fprintf(tsalert.fp, "## COLUMN%1.1d.long_name Excess mass\n",c);
      fprintf(tsalert.fp, "## COLUMN%1.1d.units m\n",c);
      fprintf(tsalert.fp, "##\n");
      for (nb = 0; nb < geom->nobc; nb++) {
	fprintf(tsalert.fp, "## COLUMN%1.1d.name  OBC%1.1d_EnergyFlux\n",
		  ++c, nb);
	fprintf(tsalert.fp, "## COLUMN%1.1d.long_name Boundary #%1.1d energy flux\n",
		c, nb);
	fprintf(tsalert.fp, "## COLUMN%1.1d.units m4s-3\n",c);
	fprintf(tsalert.fp, "##\n");
      }
      fprintf(tsalert.fp, "##\n");
    }
  }
}

/* END alerts_init()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate alert diagnostics and take action if         */
/* thresholds are volated.                                           */
/*-------------------------------------------------------------------*/
void alerts(master_t *master     /* Master data                      */
	    )
{
  FILE *afp, *fopen();         /* File handle for alerts output      */
  int ncflstep = -1;           /* No. hours to reduce CFL for        */
  int n;                       /* Counters                           */
  int dots = 0;                /* Print to ts file                   */
  double relax_rate = 600.0;   /* Relaxation rate for T/S            */
  int i, j, k, cg;
  char buf[MAXSTRLEN];
  geometry_t *geom = master->geom;

  if (!(master->alertf & NONE)) {
    if (master->alert_dt && master->t >= tsalert.tsout - DT_EPS) {
      dots = 1;
      tsalert.tsout += tsalert.tsdt;
      fflush(tsalert.fp);
    }
    /*---------------------------------------------------------------*/
    /* Calculate and print the alert diagnostics for this time-step  */
    if (strlen(master->opath))
      sprintf(buf,"%s%s.txt", master->opath, master->alertname);
    else
      sprintf(buf, "%s.txt", master->alertname);
    afp = fopen(buf, "w");
    if (afp != NULL) {
      fprintf(afp, "\nSimulation time = %f (days)\n\n",
	      master->t / 86400.0);

      /*-------------------------------------------------------------*/
      /* Print the maximum absolute sea level                        */
      fprintf(afp, "Maximum absolute sea level :\n");
      master->seta = max(master->seta, master->meta);
      cg = master->ceta;
      i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
      fprintf(afp, "%s %f at %d (%d %d %d)\n", "elevation : ",
	      master->meta, cg, i, j, k);

      /*-------------------------------------------------------------*/
      /* Print the maximum absolute velocity                         */
      fprintf(afp, "\nMaximum absolute velocity :\n");
      master->su1 = max(master->su1, master->mu1);
      cg = geom->e2ijk[master->cu1];
      i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
      fprintf(afp, "%s %f at %d (%d %d %d)\n", "u1 3D     : ",
	      master->mu1, cg, i, j, k);
      master->su1a = max(master->su1a, master->mu1a);
      cg = geom->e2ijk[master->cu1a];
      i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
      fprintf(afp, "%s %f at %d (%d %d %d)\n", "u1 2D     : ",
	      master->mu1a, cg, i, j, k);
      master->sw = max(master->sw, master->mw);
      cg = master->cw;
      i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
      fprintf(afp, "%s %f at %d (%d %d %d)\n", "w         : ",
	      master->mw, cg, i, j, k);

      /*-------------------------------------------------------------*/
      /* Print maximum divergences.                                  */
      /* 3D divergence = dw/dz                                       */
      master->sdw = max(master->sdw, master->mdw);
      cg = master->cdw;
      i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
      fprintf(afp, "%s %f at %d (%d %d %d)\n", "Div 3D    : ",
	      master->mdw, cg, i, j, k);
      /* 2D divergence = deta/dt                                     */
      master->sdeta = max(master->sdeta, master->mdeta);
      cg = master->cdeta;

      i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
      fprintf(afp, "%s %f at %d (%d %d %d)\n", "Div 2D    : ",
	      master->mdeta, cg, i, j, k);

      /*-------------------------------------------------------------*/
      /* Print the maximum tracer rates of change                    */
      master->ssmin = max(master->ssmin, master->msmin);
      cg = master->csmin;
      i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
      fprintf(afp, "%s %f at %d (%d %d %d)\n", "S min     : ",
	      master->msmin, cg, i, j, k);
      master->ssmax = max(master->ssmax, master->msmax);
      cg = master->csmax;
      i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
      fprintf(afp, "%s %f at %d (%d %d %d)\n", "S max     : ",
	      master->msmax, cg, i, j, k);
      master->sds = max(master->sds, master->mds);
      cg = master->cds;
      i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
      fprintf(afp, "%s %f at %d (%d %d %d)\n", "dS/dt     : ",
	      master->mds, cg, i, j, k);
      master->stmin = max(master->stmin, master->mtmin);
      cg = master->ctmin;
      i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
      fprintf(afp, "%s %f at %d (%d %d %d)\n", "T min     : ",
	      master->mtmin, cg, i, j, k);
      master->stmax = max(master->stmax, master->mtmax);
      cg = master->ctmax;
      i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
      fprintf(afp, "%s %f at %d (%d %d %d)\n", "T max     : ",
	      master->mtmax, cg, i, j, k);
      master->sdt = max(master->sdt, master->mdt);
      cg = master->cdt;
      i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
      fprintf(afp, "%s %f at %d (%d %d %d)\n", "dT/dt     : ",
	      master->mdt, cg, i, j, k);

      master->sshear = max(master->sshear, master->mshear);

      /* Take action on alert volations                              */
      if(master->alertf & ACTIVE && master->ts_f) {
	/* Rates of change exceed thresholds                         */
	/*    => relax harder to input data                          */
	int tno = master->tno;
	int sno = master->sno;
	if (master->mdt > master->dtmax) {
	  cg = master->cdt;
	  tr_relax_event_c(sched_get_even_by_name(schedule,
						  "tracer_relax:temp"),
			   relax_rate, master->dt, master->t, cg);
	  hd_warn("LEVEL5 alert at %f days : dT/dt exceeds threshold\n",
		  master->t);
	  master->nalert[6]++;
	}
	cg = master->ctmax;
	if (master->temp[cg] > master->maxtr[tno]) {
	  tr_relax_event_c(sched_get_even_by_name(schedule,
						  "tracer_relax:temp"),
			   relax_rate, master->dt, master->t, cg);
	  hd_warn("LEVEL5 alert at %f days : T exceeds maximum\n",
		  master->t);
	  master->nalert[6]++;
	}
	cg = master->ctmin;
	if (master->temp[cg] < master->mintr[tno]) {
	  tr_relax_event_c(sched_get_even_by_name(schedule,
						  "tracer_relax:temp"),
			   relax_rate, master->dt, master->t, cg);
	  hd_warn("LEVEL5 alert at %f days : T exceeds minimum\n",
		  master->t);
	  master->nalert[6]++;
	}
	if (master->mds > master->dsmax) {
	  cg = master->cds;
	  tr_relax_event_c(sched_get_even_by_name(schedule,
						  "tracer_relax:salt"),
			   relax_rate, master->dt, master->t, cg);
	  hd_warn("LEVEL5 alert at %f days : dS/dt exceeds threshold\n",
		  master->t);
	  master->nalert[6]++;
	}
	cg = master->csmax;
	if (master->sal[cg] > master->maxtr[sno]) {
	  tr_relax_event_c(sched_get_even_by_name(schedule,
						  "tracer_relax:salt"),
			   relax_rate, master->dt, master->t, cg);
	  hd_warn("LEVEL5 alert at %f days : S exceeds maximum\n",
		  master->t);
	  master->nalert[6]++;
	}
	cg = master->csmin;
	if (master->sal[cg] < master->mintr[sno]) {
	  tr_relax_event_c(sched_get_even_by_name(schedule,
						  "tracer_relax:salt"),
			   relax_rate, master->dt, master->t, cg);
	  hd_warn("LEVEL5 alert at %f days : S exceeds minimum\n",
		  master->t);
	  master->nalert[6]++;
	}
      }

      /* Print minimum CFL                                       */
      master->scfl = min(master->scfl, master->mcfl);
      cg = master->ccfl;
      i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
      fprintf(afp, "%s %5.2e at %d (%d %d %d)\n", "CFL       : ",
	      master->mcfl, cg, i, j, k);
      /* Take action on CFL alert volations. Model 3D time step  */
      /* is set to the CFL time-step if mcfl < grid_dt, and is   */
      /* kept at this value for ncflstep iterations after which  */
      /* it is re-set to ogrid_dt (original time-step).          */
      if(master->alertf & ACTIVE && master->cfl_f) {
	if (master->mcfl < master->grid_dt) {
	  reset_cfl_event(master, master->mcfl, CFL_ALERT);
	  master->ncfl = ncflstep * 3600 / master->grid_dt;
	} else {
	  if (master->ncfl > 0)
	    master->ncfl--;
	  else if(master->ncfl == 0) {
	    master->ncfl = -1;
	    master->grid_dt = master->ogrid_dt;
	    hd_warn("LEVEL4 alert at %f days : reset CLF timestep = %6.2f\n",
		    master->t/86400, master->grid_dt);
	    master->nalert[5]++;
	  }
	}
      }

      /* Print the maximum absolute tendencies                   */
      fprintf(afp, "\nMaximum absolute tendencies :\n");
      fprintf(afp, "u1 velocity\n");
      if (master->ma1 > 0.0) {
	master->sa1 = max(master->sa1, master->ma1);
	cg = geom->e2ijk[master->ca1];
	i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
	fprintf(afp, "%s %f at %d (%d %d %d)\n", "u1 adv    : ",
		master->ma1, cg, i, j, k);
      }
      if (master->mh1 > 0.0) {
	master->sh1 = max(master->sh1, master->mh1);
	cg = geom->e2ijk[master->ch1];
	i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
	fprintf(afp, "%s %f at %d (%d %d %d)\n", "u1 hdiff  : ",
		master->mh1, cg, i, j, k);
      }
      if (master->mv1 > 0.0) {
	master->sv1 = max(master->sv1, master->mv1);
	cg = geom->e2ijk[master->cv1];
	i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
	fprintf(afp, "%s %f at %d (%d %d %d)\n", "u1 vdiff  : ",
		master->mv1, cg, i, j, k);
      }
      if (master->mb1 > 0.0) {
	master->sb1 = max(master->sa1, master->mb1);
	cg = geom->e2ijk[master->cb1];
	i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
	fprintf(afp, "%s %f at %d (%d %d %d)\n", "u1 btp    : ",
		master->mb1, cg, i, j, k);
      }
      if (master->md1 > 0.0) {
	master->sd1 = max(master->sd1, master->md1);
	cg = geom->e2ijk[master->cd1];
	i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
	fprintf(afp, "%s %f at %d (%d %d %d)\n", "u1 bcp    : ",
		master->md1, cg, i, j, k);
      }
      if (master->mc1 > 0.0) {
	master->sc1 = max(master->sc1, master->mc1);
	cg = geom->e2ijk[master->cc1];
	i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
	fprintf(afp, "%s %f at %d (%d %d %d)\n", "u1 cor    : ",
		master->mc1, cg, i, j, k);
      }

      /* Print the mechanical energy, excess mass and energy      */
      /* flux through the open boundaries (Palma & Matano,        */
      /* 1998, section 2.3 or 2001 eqn 6).                        */
      fprintf(afp, "\nArea averaged energy :\n");
      fprintf(afp, "Mechanical energy  : %f (J/m2)\n",
	      master->me / geom->totarea);
      fprintf(afp, "Kinetic energy  : %f (J/m3)\n",
	      master->ke / geom->totarea) ;
      fprintf(afp, "Excess mass        : %f (m)\n", master->em / geom->totarea);
      if(dots) {
	fprintf(tsalert.fp, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f ",
		master->t/86400, master->meta, master->mu1, master->mu1a,
		master->mu2, master->mu2a, master->mw, master->mdeta,
		master->mdw, master->mshear, master->mcfl, master->mds, 
		master->mdt, master->me / geom->totarea, 
		master->ke / geom->totarea, master->em / geom->totarea);
      }

      fprintf(afp, "Energy flux (m3/s3) :\n");
      for (n = 0; n < geom->nobc; n++) {
	fprintf(afp, "    Boundary %d (%s) : %f\n", n,
		geom->open[n]->name, master->ef[n] / geom->open[n]->length);
	if(dots)
	  fprintf(tsalert.fp, "%f ", master->ef[n] / geom->open[n]->length);
      }
      if(dots) {
	fprintf(tsalert.fp, "\n");
	fflush(tsalert.fp);
      }

      /*-------------------------------------------------------------*/
      /* Print alert diagnostics for the simulation                  */
      fprintf(afp, "\n---------------------------------------\n");
      fprintf(afp, "Simulation maximum absolutes\n\n");
      fprintf(afp, "eta      : %f\n", master->seta);
      fprintf(afp, "u1 3D    : %f\n", master->su1);
      fprintf(afp, "u1 2D    : %f\n", master->su1a);
      fprintf(afp, "w        : %5.2e\n", master->sw);
      fprintf(afp, "div 3D   : %5.2e\n", master->sdw);
      fprintf(afp, "div 2D   : %5.2e\n", master->sdeta);
      fprintf(afp, "dS/dt    : %5.2e\n", master->sds);
      fprintf(afp, "dT/dt    : %5.2e\n", master->sdt);
      fprintf(afp, "S min    : %5.2e\n", master->ssmin);
      fprintf(afp, "S max    : %5.2e\n", master->ssmax);
      fprintf(afp, "T min    : %5.2e\n", master->stmin);
      fprintf(afp, "T max    : %5.2e\n", master->stmax);
      fprintf(afp, "CFL      : %5.2e\n", master->scfl);
      fprintf(afp, "shear    : %5.2e\n", master->sshear);
      fprintf(afp, "\n");
      if (master->sa1 > 0.0)
	fprintf(afp, "u1 adv   : %5.2e\n", master->sa1);
      if (master->sh1 > 0.0)
	fprintf(afp, "u1 hdif  : %5.2e\n", master->sh1);
      if (master->sv1 > 0.0)
	fprintf(afp, "u1 vdif  : %5.2e\n", master->sv1);
      if (master->sb1 > 0.0)
	fprintf(afp, "u1 btp   : %5.2e\n", master->sb1);
      if (master->sd1 > 0.0)
	fprintf(afp, "u1 bcp   : %5.2e\n", master->sd1);
      if (master->sc1 > 0.0)
	fprintf(afp, "u1 cor   : %5.2e\n", master->sc1);
      fprintf(afp, "\n");
      if (master->sa2 > 0.0)
	fprintf(afp, "u2 adv   : %5.2e\n", master->sa2);
      if (master->sh2 > 0.0)
	fprintf(afp, "u2 hdif  : %5.2e\n", master->sh2);
      if (master->sv2 > 0.0)
	fprintf(afp, "u2 vdif  : %5.2e\n", master->sv2);
      if (master->sb2 > 0.0)
	fprintf(afp, "u2 btp   : %5.2e\n", master->sb2);
      if (master->sd2 > 0.0)
	fprintf(afp, "u2 bcp   : %5.2e\n", master->sd2);
      if (master->sc2 > 0.0)
	fprintf(afp, "u2 cor   : %5.2e\n", master->sc2);
      fprintf(afp, "\n\nAlert summary\n");
      fprintf(afp, "LEVEL1 elevation alerts    : %d (%d%%)\n",
	      master->nalert[0], 100*master->nalert[0] / master->nstep2d);
      fprintf(afp, "LEVEL1 velocity 2D alerts  : %d (%d%%)\n",
	      master->nalert[1], 100*master->nalert[1] / master->nstep2d);
      fprintf(afp, "LEVEL1 velocity 3D alerts  : %d (%d%%)\n",
	      master->nalert[2], 100*master->nalert[2] / master->nstep);
      fprintf(afp, "LEVEL2 (divergence) alerts : %d (%d%%)\n",
	      master->nalert[3], 100*master->nalert[3] /
	      (master->nstep * master->iratio));
      fprintf(afp, "LEVEL3 (tendency) alerts   : %d (%d%%)\n",
	      master->nalert[4], 100*master->nalert[4] / master->nstep);
      fprintf(afp, "LEVEL4 (CFL) alerts        : %d (%d%%)\n",
	      master->nalert[5], 100*master->nalert[5] / master->nstep);
      fprintf(afp, "LEVEL5 (T/S) alerts        : %d (%d%%)\n",
	      master->nalert[6], 100*master->nalert[6] / master->nstep);
      fprintf(afp, "Horizontal mixing resets   : %d (%d%%)\n",
	      master->nalert[7], 100*master->nalert[7] / master->nstep);
      if (master->nstep != 0) {
	double mts = master->mclk / master->nstep;
	fprintf(afp, "Mean CPU / cell / iteration = %5.3e s\n",
		mts / master->geom->a3_t);
      }
      n = fclose(afp);
    }
    master->meta = master->mw = 0.0;
    master->mu1 = master->mu2 = 0.0;
    master->mu1a = master->mu2a = 0.0;
    master->mdeta = master->mdw = 0.0;
    master->mdt = master->mds = 0.0;
    master->msmin = master->mtmin = HUGE;
    master->msmax = master->mtmax = 0.0;
    master->ma1 = master->ma2 = 0.0;
    master->mh1 = master->mh2 = 0.0;
    master->mv1 = master->mv2 = 0.0;
    master->mb1 = master->mb2 = 0.0;
    master->md1 = master->md2 = 0.0;
    master->mc1 = master->mc2 = 0.0;
    master->me = master->ke = master->em = 0.0;
    master->mshear = 0.0;
    master->mcfl = HUGE;
    if (master->ef)
      memset(master->ef, 0, geom->nobc * sizeof(double));
    memset(master->nalert, 0, nna * sizeof(int));
  }
}

/* END alerts()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate alert diagnostics and take action if         */
/* thresholds are violated.                                          */
/*-------------------------------------------------------------------*/
void alerts_w(geometry_t *window,  /* Window geometry                */
	      int mode             /* Alert to activate              */
	      )
{
  window_t *windat = window->windat;    /* Window data               */
  win_priv_t *wincon = window->wincon;  /* Window constants          */
  int nb, cc, c, ee, e, j;     /* Counters                           */
  double *eta;                 /* Pointer to elevation               */
  double ta = 0.0;             /* Total area                         */
  double me = 0.0;             /* Mechanical energy                  */
  double ke = 0.0;             /* Kinetic energy                     */
  double em = 0.0;             /* Excess mass                        */
  double ef = 0.0;             /* OBC energy flux                    */
  double d1, d2, d3;           /* Dummies                            */
  double t;                    /* Time                               */
  double sf = 0.3;             /* Safety factor for eta maximums     */
  double relax_rate = 20.0;    /* Fixed relaxation rate for eta      */
  double rrmax = 20.0;         /* Maximum rate for linear relax_rate */
  double rrmin = 2.0;          /* Minimum rate for linear relax_rate */
  double etamax;               /* Maximum surface elevation          */
  int tide_r;                  /* Remove tide for maximum eta        */
  int rlx_diff = 1;            /* meta = eta - eta_rlx               */
  double meta, mvel;

  if (!(wincon->alertf & NONE)) {
    t = windat->t/86400;
    tide_r = wincon->tide_r;
    if (wincon->etadiff) {
      /* Standard ROAM implementation                                */
      etamax = wincon->etadiff;
      /* Alternative implementation; ramp the etadiff in to allow    */
      /* large flows during spinup.                                  */
      /*
      double dtr = wincon->rampend - wincon->rampstart;
      etamax = (windat->t - wincon->rampstart) * 
	(wincon->etadiff - wincon->etamax) / dtr + wincon->etamax;
      */
    } else
      etamax = wincon->etamax;

    /*---------------------------------------------------------------*/
    /* 2D mode elevation velocity                                    */
    if(mode & ETA_A) {

      eta = wincon->d7;
      relax_rate *= windat->dt2d;
      if (!wincon->csr_tide && tide_r == CSR_R) tide_r = 0;
      if (!(wincon->etarlx & ALERT)) rlx_diff = 0;
      if (tide_r) sf = (rlx_diff) ? 0.09 : 0.4;
      sf = (wincon->etadiff) ? 1.0 : sf;

      /* Get the maximum absolute elevation.                         */
      /* Remove the tide if required.                                */
      if (tide_r & CSR_R)
	remove_tide(window, wincon->d7);
      else if (tide_r & MEAN_R && windat->etam)
	memcpy(eta, windat->etam, window->szcS * sizeof(double));
      else {
	memcpy(eta, wincon->neweta, window->szcS * sizeof(double));
      }
      /* Remove the relaxation elevation if required                 */
      if (rlx_diff) {
	for (cc = 1; cc <= window->b2_t; cc++) {
	  c = window->w2_t[cc];
	  eta[c] -= windat->eta_rlx->val1[c];
	}
      }
      /* Find the maximum                                            */
      memset(windat->ceta, 0, window->szcS * sizeof(int));
      windat->meta = get_maxs(window, eta, 1, window->b2_t,
			     window->w2_t, windat->ceta, sf * etamax);

      /* Get the maximum 2D divergences = deta/dt                    */
      /*
      windat->mdeta = get_max(window, windat->detadt, 1, window->b2_t,
			      window->w2_t, &windat->cdeta);
      */
      windat->mdeta = get_max_div_2d(window, windat, wincon,
				     &windat->cdeta);

      /*-------------------------------------------------------------*/
      /* Take action on alert volations                              */
      if(wincon->alertf & ACTIVE) {

	/* 1. Elevation exceeds thresholds                           */
	/*    => relax hard to elevation                             */
	if (wincon->eta_f && windat->meta > sf * etamax) {
	  double d1 = sf * etamax, d2, etadif, rate;
	  if (wincon->etarlx & ALERT) {
	    for(cc = 1; (c = windat->ceta[cc]); cc++) {
	      /* Use a linear relax_rate where rate = rrmax.dt if    */
	      /* etam = sf.etamax and rate = rrmin.dt if etam >=     */
	      /* 2.sf.etamax.                                        */
	      d2 = min(windat->etam[c], 2.0 * d1);
	      etadif = fabs(eta[c]);
	      /*etadif = fabs(windat->etam[c] - windat->eta_rlx->val1[c]);*/
	      rate = windat->dt2d * max(((rrmin - rrmax) *
		     (etadif - d1) / d1 + rrmax), rrmin);
	      do_eta_relax(window, c, rate, tide_r);
	      windat->alert_a[c] = 1.0;
	      windat->alert_c[c] += 1.0;
	    }
	    windat->nalert[0]++;
	    hd_warn("LEVEL1 alert at %f days : eta exceeds threshold\n", t);
	  }
	}
	/* 2. Divergences exceeds thresholds                         */
	/*    => damp or smooth velocity field                       */
	if (wincon->div2d_f && windat->mdeta > wincon->detamax) {
	  int cl = windat->cdeta;
	  int j, e, cs = window->m2d[cl];
	  if (wincon->etarlx & ALERT) {

	    /* Relax hard to elevation                               */
	    do_eta_relax(window, cl, relax_rate, tide_r);

	    /* Increase horizontal viscosity locally                 */
	    wincon->alertf |= LEVEL2;

	    for (j = 1; j <= window->npe[cs]; j++) {
	      e = window->c2e[j][cl];
	      shapiro_smooth(window, windat->u1av, e, 0);
	    }

	    windat->nalert[3]++;
	    hd_warn("LEVEL2 alert at %f days : 2D divergence exceeds threshold\n", t);
	  }
	}
      }
    }

    /*---------------------------------------------------------------*/
    /* 2D mode velocity                                              */
    if(mode & VEL2D) {

      /* Get the maximum absolute 2D velocity                        */
      memset(windat->cu1a, 0, window->szeS * sizeof(int));
      windat->mu1a = get_maxs(window, windat->u1av, 1, window->b2_e1,
			     window->w2_e1, windat->cu1a, wincon->velmax2d);

      /* Get the mechanical energy, excess mass and energy flux      */
      /* through the open boundaries (Palma & Matano, 1998, section  */
      /* 2.3 or 2001 eqn 6).                                         */
      /* Note: Energy flux in m4s-3 is the flux in Wm-1/rho          */
      /* Wm-2 = kgs-3                                                */
      /* Jm-2 = kgs-2                                                */
      for (cc = 1; cc <= window->b2_t; cc++) {
	c = window->w2_t[cc];
	d2 = windat->eta[c] - window->botz[c];
	d3 = window->cellarea[c];
	ta += d3;
	d1 = 0.5 * (windat->uav[c] * windat->uav[c] +
		    windat->vav[c] * windat->vav[c]);
	ke += windat->dens[c] * d3 * d1;
	me += windat->dens[c] * d3 * 
	  (wincon->g * windat->eta[c] * windat->eta[c] + d2 * d1);
	em += windat->eta[c] * d3;
      }
      /* Note: the following must still be divided by total area (on */
      /* the master).                                                */
      windat->me = me;       /* Mechanical energy                    */
      windat->ke = ke;       /* Kinetic energy                       */
      windat->em = em;       /* Excess mass (mean sea level)         */

      /* OBC energy flux (note : may be positive or negative)        */
      for (nb = 0; nb < window->nobc; nb++) {
	open_bdrys_t *open = window->open[nb];
	double *cw;
	double *vel;
	int ce, *cells;
	double sgn;
	cw = window->h1au1;
	vel = windat->u1av;
	ce = open->no2_e1;
	cells = open->obc_e1;
	ef = d3 = d1 = 0.0;
	meta = mvel = 0.0;

	for (ee = 1; ee <= open->no2_e1; ee++) {
	  e = open->obc_e1[ee];
	  c = open->obc_e2[ee];
	  sgn = (open->ceni[ee]) ? 1.0 : -1.0;
	  /*if(window->s2j[c] > 18) continue;*/
	  d3 += cw[e];
	  ef += sgn * cw[e] * vel[e] * (windat->eta[c] - 
					window->botz[c] * wincon->Ds[c]) *
	    (wincon->g * windat->eta[c] +
	     0.5 * (windat->uav[c] * windat->uav[c] +
		    windat->vav[c] * windat->vav[c]));
	  meta += cw[e] * windat->eta[c];
	  mvel += vel[e] * cw[e] * (windat->eta[c] - window->botz[c] * wincon->Ds[c]);
	  d1 += cw[e] * (windat->eta[c] - window->botz[c] * wincon->Ds[c]);
	}
	windat->ef[nb] = ef;
	/*printf("%f %d %f %f\n",windat->days, nb, meta / d3, mvel / d1);*/
      }

      /*-------------------------------------------------------------*/
      /* Take action on alert volations                              */
      if(wincon->alertf & ACTIVE) {

	if (wincon->ic == 0) {
	  /* Reset the actual (non-cumulative) alert diagnostic for  */
	  /* this timestep. Note; alerts_w() is called first with    */
	  /* mode = VEL2D.                                           */
	  memset(windat->alert_a, 0, window->szcS * sizeof(double));

	  /* Convert the cumulative alert diagnostic from % to       */
	  /* actual number of alerts.                                */
	  if (windat->nstep) {
	    for(cc = 1; cc <= window->b2_t; cc++) {
	      c = window->w2_t[cc];
	      windat->alert_c[c] *= ((double)(windat->nstep - 1) / 100.0);
	    }
	  }
	}

	/* 1. Velocity exceeds thresholds                            */
	if (wincon->vel2d_f && windat->mu1a > wincon->velmax2d) {
	  /*    => cap the velocity to velmax                        */
	  for(ee = 1; (e = windat->cu1a[ee]); ee++) {
	    windat->u1av[e] = con_median(window, windat->u1av, wincon->velmax2d, 
					 e, ST_SQ3|ST_EDGE);
	    /*
	    windat->u1av[e] = (windat->u1av[e] > 0.0) ?
	                       wincon->velmax2d : -wincon->velmax2d;
	    */
	    /*shuman_smooth(window, windat->u1av, c);*/
	    c = window->e2c[e][0];
	    windat->alert_a[c] = 2.0;
	    windat->alert_c[c] += 1.0;
	  }
	  /*    => increase horizontal viscosity locally             */
	  wincon->alertf |= LEVEL1_UA;
	  hd_warn("LEVEL1 alert at %f days : u1av exceeds threshold\n", t);
	  windat->nalert[1]++;
	}
	if (wincon->shear_f) windat->mshear = 0.0;
	if (wincon->vel2d_f && wincon->shear_f) {
	  int ym1, cl = (windat->mu1a > wincon->velmax2d) ? ee : 1;
	  double shear;

	  vel_grad(window, windat, wincon, windat->u1av, windat->u2av,
		   wincon->d1, wincon->d2, GRAD_2D|GRAD_TAN);

	  for (ee = 1; ee <= window->v2_e1; ee++) {
	    e = window->w2_e1[ee];
	    shear = wincon->d1[e];
	    if (fabs(shear) > wincon->smax) {
	      windat->cu1a[cl] = e;
	      cl++;
	      wincon->alertf |= LEVEL1_UA;
	    }
	    if (fabs(shear) > windat->mshear)
	      windat->mshear = fabs(shear);
	  }
	}
      }
    }

    /*---------------------------------------------------------------*/
    /* 3D mode                                                       */
    if(mode & VEL3D) {

      /* Reset momentum omissions and alert level                    */
      wincon->u1_f = NONE;

      /* Get the maximum absolute velocity                           */
      memset(windat->cu1, 0, window->sze * sizeof(int));
      windat->mu1 = get_maxs(window, windat->u1, 1, window->b3_e1,
			     window->w3_e1, windat->cu1, wincon->velmax);

      /* Get the maximum absolute tendencies                         */
      if (windat->u1_adv) {
	windat->ma1 = get_max(window, windat->u1_adv, 1, window->b3_e1,
			      window->w3_e1, &windat->ca1);
      }
      if (windat->u1_hdif) {
	windat->mh1 = get_max(window, windat->u1_hdif, 1, window->b3_e1,
			      window->w3_e1, &windat->ch1);
      }
      if (windat->u1_vdif) {
	windat->mv1 = get_max(window, windat->u1_vdif, 1, window->b3_e1,
			      window->w3_e1, &windat->cv1);
      }
      if (windat->u1_btp) {
	windat->mb1 = get_max(window, windat->u1_btp, 1, window->b3_e1,
			      window->w3_e1, &windat->cb1);
      }
      if (windat->u1_bcp) {
	windat->md1 = get_max(window, windat->u1_bcp, 1, window->b3_e1,
			      window->w3_e1, &windat->cd1);
      }
      if (windat->u1_cor) {
	windat->mc1 = get_max(window, windat->u1_cor, 1, window->b3_e1,
			      window->w3_e1, &windat->cc1);
      }
      /* 3D divergence = dw/dz                                       */
      windat->mdw = 0.0;
      for (cc = 1; cc <= window->v2_t; cc++) {
	c = window->w2_t[cc];
	d1 = fabs((windat->wtop[c] - windat->w[c]) / wincon->dz[c]);
	if (d1 > windat->mdw) {
	  windat->mdw = d1;
	  windat->cdw = c;
	}
      }
      for (; cc <= window->v3_t; cc++) {
	c = window->w3_t[cc];
	d1 = fabs((windat->w[c] - windat->w[window->zm1[c]]) / wincon->dz[c]);
	if (d1 > windat->mdw) {
	  windat->mdw = d1;
	  windat->cdw = c;
	}
      }

      /*-------------------------------------------------------------*/
      /* Take action on alert volations                              */
      if(wincon->alertf & ACTIVE) {

	/* 1. Velocity or elevation exceeds thresholds               */
	if (wincon->vel3d_f && windat->mu1 > wincon->velmax) {

	  /*    => cap the velocity to velmax                        */
	  for(ee = 1; (e = windat->cu1[ee]); ee++) {
	    windat->u1[e] = con_median(window, windat->u1, wincon->velmax, 
				       e, ST_SQ3|ST_EDGE);
	    /*
	    windat->u1[e] = (windat->u1[e] > 0.0) ?
	                     wincon->velmax : -wincon->velmax;
	    */
	    /*shuman_smooth(window, windat->u1, c);*/
	    c = window->e2ijk[e];
	    windat->alert_a[window->m2d[c]] = 4.0;
	    windat->alert_c[window->m2d[c]] += 1.0;
	  }
	  /*    => increase horizontal viscosity locally             */
	  wincon->alertf |= LEVEL1_U;
	  hd_warn("LEVEL1 alert at %f days : u1 exceeds threshold\n", t);
	  windat->nalert[2]++;
	}

	if (wincon->vel3d_f && wincon->shear_f) {
	  int cl = (windat->mu1 > wincon->velmax) ? ee : 1;
	  double shear;

	  vel_grad(window, windat, wincon, windat->u1av, windat->u2av,
		   wincon->w5, wincon->w6, GRAD_2D|GRAD_TAN);

	  for (ee = 1; ee <= window->v3_e1; ee++) {
	    e = window->w3_e1[ee];
	    shear = wincon->w5[e];
	    if (fabs(shear) > wincon->smax) {
	      windat->cu1[cl] = c;
	      cl++;
	      wincon->alertf |= LEVEL1_U;
	    }
	    if (fabs(shear) > windat->mshear)
	      windat->mshear = fabs(shear);
	  }
	  hd_warn("LEVEL1 alert at %f days : u1 shear exceeds threshold\n", t);
	}

	/* 2. Tendencies exceeds threshold                           */
	/*    => omit the process                                    */
	if(wincon->tend_f) {
	  if (windat->ma1 > wincon->amax) {
	    wincon->u1_f |= ADVECT;
	    wincon->u1av_f |= ADVECT;
	    hd_warn("LEVEL3 alert at %f days : u1 adv tendency exceeds threshold\n", t);
	    windat->nalert[4]++;
	  }
	  if (windat->mh1 > wincon->hmax) {
	    wincon->u1_f |= HDIFF;
	    wincon->u1av_f |= HDIFF;
	    hd_warn("LEVEL3 alert at %f days : u1 hdif tendency exceeds threshold\n", t);
	    windat->nalert[4]++;
	  }
	  if (windat->mv1 > wincon->vmax) {
	    wincon->u1_f |= VDIFF;
	    wincon->u1av_f |= VDIFF;
	    hd_warn("LEVEL3 alert at %f days : u1 vdif tendency exceeds threshold\n", t);
	    windat->nalert[4]++;
	  }
	  if (windat->mb1 > wincon->btmax) {
	    wincon->u1_f |= PRESS_BT;
	    wincon->u1av_f |= PRESS_BT;
	    hd_warn("LEVEL3 alert at %f days : u1 btp tendency exceeds threshold\n", t);
	    windat->nalert[4]++;
	  }
	  if (windat->md1 > wincon->bcmax) {
	    wincon->u1_f |= PRESS_BC;
	    wincon->u1av_f |= PRESS_BC;
	    hd_warn("LEVEL3 alert at %f days : u1 bcp tendency exceeds threshold\n", t);
	    windat->nalert[4]++;
	  }
	  if (windat->mc1 > wincon->cmax) {
	    wincon->u1_f |= CORIOLIS;
	    wincon->u1av_f |= CORIOLIS;
	    hd_warn("LEVEL3 alert at %f days : u1 cor tendency exceeds threshold\n", t);
	    windat->nalert[4]++;
	  }
	}
	/* 3. Divergences exceeds thresholds                         */
	/*    => damp or smooth velocity field                       */
	if (wincon->div3d_f && windat->mdw > wincon->dwmax) {
	  /*    => increase horizontal viscosity locally             */
	  int cl = windat->cdw;
	  int j, e, cs = window->m2d[cl];
	  wincon->alertf |= LEVEL2w;

	  for (j = 1; j <= window->npe[cs]; j++) {
	    e = window->c2e[j][cl];
	    shapiro_smooth(window, windat->u1, e, 0);
	  }
	  hd_warn("LEVEL2 alert at %f days : 3D divergence exceeds threshold\n", t);
	  windat->nalert[3]++;
	}
      }
    }

    /*---------------------------------------------------------------*/
    /* Vertical velocity                                             */
    if(mode & WVEL) {
      /* Get the maximum absolute velocity                           */
      windat->mw = get_max(window, windat->w, 1, window->b3_t,
			   window->w3_t, &windat->cw);

      /*-------------------------------------------------------------*/
      /* Take action on alert volations                              */
      if(wincon->alertf & ACTIVE) {
	/* 1. Velocity or elevation exceeds thresholds               */
	if (wincon->wvel_f && windat->mw > wincon->wmax) {
	  /*    => increase horizontal viscosity locally             */
	  wincon->alertf |= LEVEL1_W;
	  hd_warn("LEVEL1 alert at %f days : w exceeds threshold\n", t);
	  windat->nalert[2]++;
	}
      }
      /*---------------------------------------------------------------*/
      /* Courant criterion                                             */
      windat->mcfl = calc_min_cfl_3d(window, windat, wincon, &windat->ccfl);
    }

    /*---------------------------------------------------------------*/
    /* Temperature and salinity                                      */
    if(mode & TRACERS) {
      /* Get the minimum and maximum T and S                         */
      windat->msmax = get_max(window, windat->sal, 1, window->b3_t,
			      window->w3_t, &windat->csmax);
      windat->mtmax = get_max(window, windat->temp, 1, window->b3_t,
			      window->w3_t, &windat->ctmax);
      windat->msmin = get_min(window, windat->sal, 1, window->b3_t,
			      window->w3_t, &windat->csmin);
      windat->mtmin = get_min(window, windat->temp, 1, window->b3_t,
			      window->w3_t, &windat->ctmin);
      /* Get the maximum rate of change of T                         */
      windat->mdt = get_ddt(window, windat->temp, windat->tempb,
			    &windat->cdt);
      /* Get the maximum rate of change of S                         */
      windat->mds = get_ddt(window, windat->sal, windat->salb,
			    &windat->cds);
      if(wincon->alertf & ACTIVE) {

	/* Convert the cumulative alert diagnostic from % to actual  */
	/* number of alerts. Note; alerts_w() is called last with    */
	/* mode = TRACERS.                                           */
	if (windat->nstep) {
	  for(cc = 1; cc <= window->b2_t; cc++) {
	    c = window->w2_t[cc];
	    windat->alert_c[c] *= (100.0 / (double)(windat->nstep));
	  }
	}
      }
    }
  }
}

/* END alerts_w()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to remove the tidal component from an elevation signal    */
/* using the csr tide model.                                         */
/*-------------------------------------------------------------------*/
void remove_tide(geometry_t *window,  /* Window geometry             */
		 double *eta          /* eta with tide removed       */
	      )
{
  window_t *windat = window->windat;    /* Window data               */
  win_priv_t *wincon = window->wincon;  /* Window constants          */

  if (wincon->tc.nt) {
    int c, cc;
    double gmt = windat->t - wincon->tz;
    double tide = 0.0;
    double ramp = (wincon->rampf & TIDALH) ? windat->rampval : 1.0;
    for(cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      /* Get an estimate of tidal height via csr tide model          */
      tide = ramp * csr_tide_eval(&wincon->tc, cc, gmt);
      eta[c] = wincon->neweta[c] - tide;
    }
  } else {
    memcpy(eta, wincon->neweta, window->szcS * sizeof(double));
  }
}

/* END remove_tide()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to implement a Shuman type smoothing on an array (see     */
/* Kowalik and Murty, p145, Eqn 3.135 over a star shaped stencil.    */
/*-------------------------------------------------------------------*/
void shuman_smooth(geometry_t *window,  /* Window geometry           */
		   double *a,     /* Array to smooth                 */
		   int cl         /* Location to smooth              */
		   )
{
  double v = 0.5;       /* Damping for 0 < v < 0.5                   */
  int *st, size = 7;    /* Inner stencil                             */
  int *sst, ssize = 5;  /* Outer stencil                             */
  int c, cc;
  double *b = window->wincon->w4;

  st = NULL;
  sst = stencil(window, cl, &ssize, ST_STAR, 0);
  for (cc = 0; cc < ssize; cc++) {
    c = sst[cc];
    st = stencil(window, c, &size, ST_SQ3, 0);

    b[c] = a[c] + 0.5 * v * (1.0 - v) *
    (a[st[1]] + a[st[2]] + a[st[3]] + a[st[4]] - 4.0 * a[st[0]]) +
    0.25 * v * v * (a[st[5]] + a[st[6]] + a[st[7]] + a[st[8]] -
		    4.0 * a[st[0]]);
 }
  for (cc = 0; cc < ssize; cc++) {
    c = sst[cc];
    a[c] = b[c];
  }
  i_free_1d(st);
  i_free_1d(sst);
}

/* END shuman_smooth()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to implement a Shuman type smoothing on an array (see     */
/* Kowalik and Murty, p145, Eqn 3.135.                               */
/*-------------------------------------------------------------------*/
double shuman(geometry_t *window,   /* Window geometry               */
	      double *a,            /* Array to smooth               */
	      int c                 /* Location to smooth            */
	      )
{
  double v = 0.5;       /* Damping for 0 < v < 0.5                   */
  int *st = NULL;       /* Inner stencil                             */
  int size = 7;         /* Inner stencil size                        */
  double ret;           /* Return value                              */

  st = stencil(window, c, &size, ST_SQ3, 0);
  ret = a[c] + 0.5 * v * (1.0 - v) *
    (a[st[1]] + a[st[2]] + a[st[3]] + a[st[4]] - 4.0 * a[st[0]]) +
    0.25 * v * v * (a[st[5]] + a[st[6]] + a[st[7]] + a[st[8]] -
		    4.0 * a[st[0]]);
  i_free_1d(st);
  return(ret);
}

/* END shuman()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to implement a Shuman type smoothing on an array (see     */
/* Kowalik and Murty, p145, Eqn 3.135.                               */
/*-------------------------------------------------------------------*/
void shuman_3d(geometry_t *window,  /* Window geometry             */
		 double *var,         /* Array to smooth             */
		 double v
		 )
{
  win_priv_t *wincon = window->wincon;
  int c, cc, cs, cn, j;
  double *a = wincon->w5, wgt, npe;

  memcpy(a, var, window->szc * sizeof(double));
  /* Note not a Shuman filter on unstructured grids                  */
  for (cc = 1; cc < window->b3_t; cc++) {
    c = window->w3_t[cc];
    cs = window->m2d[c];
    npe = (double)window->npe[cs];
    wgt = (npe + 1.0 - v) / (npe * (npe + 1.0));
    var[c] = v * a[c] / (npe + 1.0);
    for (j = 1; j <= window->npe[cs]; j++) {
      cn = window->c2c[j][c];
      var[c] += wgt * a[cn];
    }
  }
}

/* END shuman_3d()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to implement a Shapiro type smoothing on an array in the  */
/* x direction (see Kowalik and Murty, p146, Eqn 3.136).             */
/*-------------------------------------------------------------------*/
void shapiro_smooth(geometry_t *window,   /* Window geometry         */
		    double *a,            /* Array to smooth         */
		    int cl,               /* Location to smooth      */
		    int mode              /* 0=edge, >0 = cenre      */
		    )
{
  int *mp = (!mode) ? window->ep : window->c2c[mode];
  int *mm = (!mode) ? window->em : window->c2c[jo(mode, window->npe[window->m2d[cl]])];
  int m1 = mm[cl];
  int p1 = mp[cl];
  int m2 = mm[m1];
  int p2 = mp[p1];
  int m3 = mm[m2];
  int p3 = mp[p2];
  int m4 = mm[m3];
  int p4 = mp[p3];

  a[cl] = (186 * a[cl] + 56 * (a[m1] + a[p1]) -
	   28 * (a[m2] + a[p2]) + 8 * (a[m3] + a[p3]) -
	   (a[m4] + a[p4])) / 256.0;
}

/* END shapiro_smooth()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to implement a Shapiro type smoothing on an array. Taken  */
/* Shapiro, R. (1975) Linear Filtering. Mathematics of Computation,  */
/* 29, 1094-1097. Stencils are taken from Table 2. For 3D arrays.    */
/*-------------------------------------------------------------------*/
void shapiro(geometry_t *window,     /* Window geometry              */
	     double *a,              /* Array to smooth              */
	     double *buf,            /* Dummy buffer                 */
	     int *vec,               /* Cells to process             */
	     int nvec,               /* Size of vec                  */
	     int order,              /* Order of the filter          */
	     int mode,               /* 0=centres, 1=edges           */
	     int dir                 /* Smoothing direction          */
	     )
{
  int c, cc;
  int cp, cm, n;
  int *p1, *m1;
  double ret, sgn;

  order -= 1;
  if (mode)
    memcpy(buf, a, window->sze * sizeof(double));
  else
    memcpy(buf, a, window->szc * sizeof(double));
  if (dir & XDIR) {
    p1 = window->ep;
    m1 = window->em;
  } else if (dir & ZDIR) {
    p1 = window->zp1;
    m1 = window->zm1;
  }
  for (cc = 1; cc <= nvec; cc++) {
    c = vec[cc];
    a[c] = buf[c] * fsw[order][0];
    cp = cm = c;
    for (n = 1; n <= order; n++) {
      cp = p1[cp];
      cm = m1[cp];
      sgn = (pow(-1.0, (double)order) < 0.0) ? 1.0 : -1.0;
	a[c] += (sgn * fsw[order][n] * (buf[cm] + buf[cp]));
    }
    a[c] /= (pow(2.0, 2.0*(double)order+2.0));
  }
}

/* END shapiro()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to make an array containing the size x size locations of  */
/* a stencil centered on location c. 'size' must be odd.             */
/*-------------------------------------------------------------------*/
int *stencil(geometry_t *window,    /* Window geometry               */
             int cl,                /* Central sparse coordinate     */
	     int *size,             /* Stencil size                  */
	     int type,              /* Type of stencil               */
	     int edge               /* Edge number for edge stencils */
	     )
{
  int cc, c, c1s, j, m;
  int *st = NULL, *se = NULL;
  int *mask = window->sask;
  int sz = floor(*size / 2);
  int ss;
  int cs = window->m2d[cl];

  if (type & ST_STAR) {
    /* Star shaped stencil. Stencil filled from m=0; m < ss;           */
    ss = window->npem * sz + 1;
    st = i_alloc_1d(ss);
    *size = ss;

    /* Get the indices spreading outward from the centre             */
    ss = 0;
    st[ss++] = cl;
    for (j = 1; j <= window->npe[cs]; j++) {
      c = cl;
      for (m = 1; m <  sz; m++) {
	c = window->c2c[j][c];
	st[ss++] = c;
      }
    }
  } else if(type & ST_SQ3) {
    /* Stencil of first cells surrounding a cell centre (e.g.        */
    /* equivalent to a 9 point stencil on a square grid).            */
    int v, vs, i, ss, c1;
    ss = 1;
    memset(mask, 0, window->szcS * sizeof(int));
    mask[cs] = 1;
    /* Count the number of cells in the stencil                      */
    for (j = 1; j <= window->npe[cs]; j++) {
      v = window->c2v[j][cl];
      vs = window->m2dv[v];
      for (i = 1; i <= window->nvc[vs]; i++) {
	c1 = window->v2c[v][i];
	c1s = window->m2d[c1];
	if (type & ST_EDGE && window->wgst[c1]) continue;
	if (!mask[c1s]) {
	  ss++;
	  mask[c1s] = 1;
	}
      }
    }
    *size = ss;

    st = i_alloc_1d(ss);
    ss = 1;
    st[0] = cl;
    memset(mask, 0, window->szcS * sizeof(int));
    mask[cs] = 1;
    for (j = 1; j <= window->npe[cs]; j++) {
      v = window->c2v[j][cl];
      vs = window->m2dv[v];
      for (i = 1; i <= window->nvc[vs]; i++) {
	c1 = window->v2c[v][i];
	c1s = window->m2d[c1];
	if (type & ST_EDGE && window->wgst[c1]) continue;
	if (!mask[c1s]) {
	  st[ss++] = c1;
	  mask[c1s] = 1;
	}
      }
    }
    *size  = mergeSort(st, 0, *size-1) + 1;
  } else if(type & ST_SIZED) {
    /* Stencil surrounding a cell centre.                            */
    /* This stencil is defined by finding the vertices surrounding a */
    /* cell centre, then populating the stencil with unique cell     */
    /* centres associated with each vertex.                          */
    int i, j, mm, v, vs, c1;
    int *sta = i_alloc_1d(window->szcS);
    int sl, ps = 0;

    if (sz == 0) {
      *size = 1;
      st = i_alloc_1d(*size);
      st[0] = cl;
      return(st);
    }

    memset(mask, 0, window->szcS * sizeof(int));
    sta[0] = cl;
    ss = 1;
    mask[window->m2d[cl]] = 1;
    for (mm = 1; mm <= sz; mm++) {
      sl = ss;
      for (cc = ps; cc < sl; cc++) {
	c = sta[cc];
	cs = window->m2d[c];
	for (j = 1; j <= window->npe[cs]; j++) {
	  v = window->c2v[j][c];
	  vs = window->m2dv[v];
	  for (i = 1; i <= window->nvc[vs]; i++) {
	    c1 = window->v2c[v][i];
	    if (c1 && !mask[window->m2d[c1]] && window->cask[c1] & W_TRA) {
	      sta[ss++] = c1;
	      mask[window->m2d[c1]] = 1;
	    }
	  }
	}
      }
      ps = sl;
    }
    *size = ss;
    st = i_alloc_1d(ss);    
    for (cc = 0; cc < ss; cc++) {
      st[cc] = sta[cc];
    }
    i_free_1d(sta);
  } else if(type & ST_SIZEO) {
    /* Stencil surrounding a cell centre.                            */
    /* This stencil is defined by finding the vertices surrounding a */
    /* cell centre, then populating the stencil with unique cell     */
    /* centres associated with each vertex. Loops over the whole 2D  */
    /* grid, so it is not very efficient. Can be used for any        */
    /* stencil size.                                                 */
    int mm, v, vs, j, i, ss, c1;
    int *oask, *vask;

    if (sz == 0) {
      *size = 1;
      st = i_alloc_1d(*size);
      st[0] = cl;
      return(st);
    }
    oask = i_alloc_1d(window->szcS);
    vask = i_alloc_1d(window->szvS);
    ss = 1;
    /* Count the number of cells in the stencil                      */
    memset(mask, 0, window->szcS * sizeof(int));
    memset(oask, 0, window->szcS * sizeof(int));
    memset(vask, 0, window->szvS * sizeof(int));
    oask[window->m2d[cl]] = 1;
    for (mm = 1; mm <= sz; mm++) {
      memcpy(mask, oask, window->szcS * sizeof(int));
      for (cc = 1; cc <= window->b2_t; cc++) {
	c = window->w2_t[cc];
	if (mask[window->m2d[c]]) {
	  for (j = 1; j <= window->npe[c]; j++) {
	    v = window->c2v[j][c];
	    vs = window->m2dv[v];
	    if (!vask[vs]) {
	      vask[vs] = 1;
	      for (i = 1; i <= window->nvc[vs]; i++) {
		c1 = window->v2c[v][i];
		if (!oask[window->m2d[c1]]) {
		  oask[window->m2d[c1]] = 1;
		  ss++;
		}
	      }
	    }
	  }
	}
      }
    }
    *size = (*size == 1) ? 1 : ss;
    st = i_alloc_1d(ss);
    /* Allocate the stencil cells                                    */
    ss = 1;
    st[0] = cl;
    memset(mask, 0, window->szcS * sizeof(int));
    memset(oask, 0, window->szcS * sizeof(int));
    memset(vask, 0, window->szvS * sizeof(int));
    oask[window->m2d[cl]] = 1;
    for (mm = 1; mm <= sz; mm++) {
      memcpy(mask, oask, window->szcS * sizeof(int));
      for (cc = 1; cc <= window->b2_t; cc++) {
	c = window->w2_t[cc];
	if (mask[window->m2d[c]]) {
	  for (j = 1; j <= window->npe[c]; j++) {
	    v = window->c2v[j][c];
	    vs = window->m2dv[v];
	    if (!vask[vs]) {
	      vask[vs] = 1;
	      for (i = 1; i <= window->nvc[vs]; i++) {
		c1 = window->v2c[v][i];
		if (!oask[window->m2d[c1]]) {
		  oask[window->m2d[c1]] = 1;
		  st[ss++] = c1;
		}
	      }
	    }
	  }
	}
      }
    }
    i_free_1d(oask);
    i_free_1d(vask);
  }

  /*-----------------------------------------------------------------*/
  /* Extract edges if required                                       */
  if (type & ST_EDGE) {
    int e, es;
    /* Extract 'edge' from each of the centres in the stencil        */
    ss = 0;
    for (cc = 0; cc < *size; cc++) {
      c = st[cc];
      e = window->c2e[edge][c];
      if (e) {
	st[ss] = e;
	ss++;
      }
    }
    *size = ss;
    if (type & ST_ALLE || edge < 0) {
      /* Extract all edges around each centre in the stencil         */
      memset(mask, 0, window->szeS * sizeof(int));
      ss = 0;
      for (cc = 0; cc < *size; cc++) {
	c = st[cc];
	cs = window->m2d[c];
	for (j = 1; j <= window->npe[cs]; j++) {
	  e = window->c2e[j][c];
	  es = window->m2de[e];
	  if (e && !mask[es]) {
	    ss++;
	    mask[es] = 1;
	  }
	}
      }
      se = i_alloc_1d(ss);
      *size = ss;
      memset(mask, 0, window->szeS * sizeof(int));
      ss = 0;
      for (cc = 0; cc < *size; cc++) {
	c = st[cc];
	cs = window->m2d[c];
	for (j = 1; j <= window->npe[cs]; j++) {
	  e = window->c2e[j][c];
	  es = window->m2de[e];
	  if (e && !mask[es]) {
	    st[ss] = e;
	    ss++;
	    mask[es] = 1;
	  }
	}
      }
      i_free_1d(st);
      st = se;
    }
  }
  return(st);
}

/* END stencil()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Creates a map of cells surrounding a cell centre                  */
/*-------------------------------------------------------------------*/
void get_filter(geometry_t *geom) {
  int c, cc, c1, c2;
  int *st = NULL, sz, szc, szm, n;
  double d1, d2;
  df_filter_t *f;

  geom->filter = (df_filter_t **)malloc(sizeof(df_filter_t *) * geom->szc);

  for (cc = 1; cc <=geom->b3_t; cc++) {
    c = geom->w3_t[cc];
    f = (df_filter_t *)malloc(sizeof(df_filter_t));
    memset(f, 0, sizeof(df_filter_t));

    /* Get the size of 1 ring around the centre. Note: this may not  */
    /* equal to npe near solid boundaries.                           */
    sz = 3;
    st = stencil(geom, c, &sz, ST_SIZED, 0);
    f->map3 = sz;
    /* Get the stencil for 2 rings around the centre                 */
    sz = f->size = 5;
    st = stencil(geom, c, &sz, ST_SIZED, 0);
    f->map5 = sz;
    f->map = i_alloc_1d(f->map5);
    f->k = d_alloc_1d(f->map5);
    d1 = 0.0;
    for (n = 0; n < f->map5; n++) {
      f->map[n] = st[n];
      f->k[n] = 1.0 / (double)f->map5;
      d1 += f->k[n];
    }
    f->f = (double)f->map5 / (double)f->map3;
    geom->filter[c] = f;
  }
}

/* END get_filter()                                                  */
/*-------------------------------------------------------------------*/



int merge(int a[], int start, int new_mid, int mid, int end)
{
  int i,j,k;
  int temp[end-start+1];
 
  i = start; j = mid; k = 0;
 
  while(i <= new_mid && j <= end) {
    if(a[i]< a[j]) {
      temp[k++]= a[i++];
    }
    else if(a[i] > a[j]) {
      temp[k++] = a[j++];
    } else {
      temp[k++] = a[j++];
      i++;
    }
  }
  while(i <= new_mid) {
    temp[k++] = a[i++]; 
  }
  while(j <= end){
    temp[k++] = a[j++]; 
  }
  for(i=0; i<k; i++){
    a[i+start] = temp[i];
  }
  return k;
}

int mergeSort(int a[], int start, int end )
{
  int mid  = (start + end)/2;

  if(start < end) {
    int new_mid = mergeSort(a, start, mid);
    int new_end = mergeSort(a, mid+1, end);
    int k = merge(a, start, new_mid, mid+1, new_end);
    return start+k-1;
  }
  return end;
}

/*-------------------------------------------------------------------*/
/* Calculates the divergence of 3D fluxes (equal to dw/dz). This     */
/* code is replicated in set_flux_3d(), but needs to be recalculated */
/* here so that if divergences are too large then 3D velocities can  */
/* be smoothed spatially (in order to reduce divergence) before      */
/* being used in vel_w_update() to calculate the vertical velocity.  */
/* Using dwdz to get to 2D divergence means any smoothing lags the   */
/* computed wvel by one time-step (but is more computationally       */
/* efficient).                                                       */
/*-------------------------------------------------------------------*/
double get_max_div_3d(geometry_t *window,  /* Window geometry        */
		      window_t *windat,    /* Window data            */
		      win_priv_t *wincon,  /* Window constants       */
		      int *cl              /* Location of maximum    */
)
{
  double *u1flux = wincon->w4;
  double *u2flux = wincon->w5;
  double colflux, mdiv = 0.0;
  int ee, e, es, cc, c, cs, n;

  /*-----------------------------------------------------------------*/
  /* Get the fluxes at time t over the whole window                  */
  memset(u1flux, 0, window->sze * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Fluxes through the e1 faces                                     */
  for (ee = 1; ee <= window->a3_e1; ee++) {
    e = window->w3_e1[ee];
    es = window->m2de[e];
    u1flux[e] = windat->u1[e] * windat->dzu1[e] * window->h2au1[es] *
      wincon->mdx[es];
  }

  /*-----------------------------------------------------------------*/
  /* Get the divergence of velocity transport                        */
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    cs = window->m2d[c];

    colflux = 0.0;
    for (n = 1; n <= window->npe[cs]; n++) {
      e = window->c2e[n][c];
      colflux += window->eSc[n][c] * u1flux[e];
    }
    colflux += window->cellarea[cs] * (windat->w[window->zp1[c]] - windat->w[c]);
    if (fabs(colflux) > mdiv) {
      mdiv = fabs(colflux);
      *cl = c;
    }
  }
  return(mdiv);
}

/* END get_max_div_3d()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates the divergence of 2D fluxes (equal to detadt : ms-1).  */
/* This code is replicated in eta_step(), but needs to be re-        */
/* calculated here so that if divergences are too large then 2D      */
/* velocities can be smoothed spatially (in order to reduce          */
/* divergence) before being used in eta_step() to calculate the      */
/* elevation. Using detadt to get to 2D divergence means any         */
/* smoothing lags the computed elevation by one time-step.           */
/*-------------------------------------------------------------------*/
double get_max_div_2d(geometry_t *window,  /* Window geometry        */
		      window_t *windat,    /* Window data            */
		      win_priv_t *wincon,  /* Window constants       */
		      int *cl              /* Location of maximum    */
)
{
  double *colflux = wincon->d1;
  double *u1flux = wincon->d2;
  double mdiv = 0.0;
  int ee, e, cc, c, n;

  /*-----------------------------------------------------------------*/
  /* Get the fluxes at time t over the whole window                  */
  memset(u1flux, 0, window->szeS * sizeof(double));
  /* Calculate the flux at e1 wet and boundary cells                 */
  for (ee = 1; ee <= window->a2_e1; ee++) {
    e = window->w2_e1[ee];
    u1flux[e] = windat->u1av[e] * windat->depth_e1[e] * window->h1au1[e] *
      wincon->mdx[e];
  }
  /*-----------------------------------------------------------------*/
  /* Get the divergence of velocity transport                        */
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->w2_t[cc];

    colflux[c] = 0.0;
    for (n = 1; n <= window->npe[c]; n++) {
      e = window->c2e[n][c];
      colflux[c] += window->eSc[n][c] * u1flux[e];
    }
    colflux[c] /= window->cellarea[c];
    if (fabs(colflux[c]) > mdiv) {
      mdiv = fabs(colflux[c]);
      *cl = c;
    }
  }
  return(mdiv);
}

/* END get_max_div_2d()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to fine maximum rates of change                           */
/*-------------------------------------------------------------------*/
double get_ddt(geometry_t *window,         /* Window geometry        */
	       double *T,                  /* Tracer at time t       */
	       double *TB,                 /* Tracer at time t-1     */
	       int *cl                     /* Location of maximum    */
)
{
  window_t *windat = window->windat;
  double dTdt;
  double dtmax = 0.0;
  int cc, c;

  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    dTdt = fabs((T[c] - TB[c]) / windat->dt);
    if (dTdt > dtmax) {
      dtmax = dTdt;
      *cl = c;
    }
  }
  return(dtmax);
}

/* END get_ddt()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculate the minimum resolution for a given run time ratio and   */
/* execution speed.                                                  */
/*-------------------------------------------------------------------*/
void calc_cfl(geometry_t *window,   /* Window geometry               */
              window_t *windat,     /* Window data                   */
              win_priv_t *wincon    /* Window constants              */
  )
{
  double cfl2d, cfl3d;          /* CFL conditions                    */
  int lc, c, cc, zp1, j;        /* Sparse locations                  */
  double dh, u, v, w;           /* Grid spacing                      */
  double spd, spd2d, iws, cspd; /* Current / wave speeds             */
  double minval = 1e-10;        /* Minimum value for velocity        */
  double eps = 1e5;             /* Maximum CFL timestep              */
  int traf;

  traf = (wincon->trasc & (LAGRANGE|FFSL)) ? 1 : 0;

  /* Get the CFL conditions */
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = lc = window->w2_t[cc];

    /* Get the internal wave speed                                   */
    /* c=NH/(n.pi), n=1,2,3... (Gill, 1982, p159) for long waves     */
    /* constant N.                                                   */
    iws = int_wave_spd_cont_w(window, windat, wincon, c, cc);
    /* Use a 2 layer approximation for shallow water internal waves  */
    iws = int_wave_spd(window, windat, wincon, c, cc);

    /* Find the maximum velocity in the water column                 */
    windat->cour[c] = 0.0;
    windat->lips[c] = 0.0;
    windat->ahsb[c] = 0.0;
    cspd = 0;
    while (lc <= window->bot_t[cc]) {

      if (windat->eta[c] >= window->gridz[lc]) {
	/* Advective speed */
	spd = sqrt(windat->u[lc] * windat->u[lc] +
		   windat->v[lc] * windat->v[lc]);
	if (spd > cspd)
	  cspd = spd;
	windat->cour[c] = HUGE;
	windat->lips[c] = HUGE;
	for (j = 1; j < window->npe[c]; j++) {
	  int e = window->c2e[j][c];
	  int es = window->m2de[e];
	  u = fabs(windat->u1[e]);
	  windat->courn[c] = max(windat->courn[c], u * windat->dt / window->h2au1[es]);
	  /*windat->courn[c] = u * windat->dt / window->h2au1[es];*/
	  u = (u) ? window->h2au1[es] / u : HUGE;
	  windat->cour[c] = min(windat->cour[c], u);
	  u = fabs(windat->u[window->e2c[e][0]] - 
		   windat->u[window->e2c[e][1]]);
	  u = (u) ? window->h2au1[es] / u : HUGE;
	  windat->lips[c] = min(windat->lips[c], u);
	}
	w = fabs(0.5 * (windat->w[window->zp1[lc]] + windat->w[lc]));
	w = (w) ? wincon->dz[lc] / w : HUGE;
	windat->cour[c] = min(windat->cour[c], w);
	w = fabs(windat->w[window->zp1[lc]] - windat->w[lc]);
	w = (w) ? wincon->dz[lc] / w : HUGE;
	windat->lips[c] = min(windat->lips[c], w);
	dh = 0.5 * sqrt(window->cellarea[c]);
	dh = (wincon->u1kh[lc]) ? dh / (4.0 * wincon->u1kh[lc]) : 0.0;
	windat->ahsb[c] = max(windat->ahsb[c], dh);
      }
      lc = window->zm1[lc];
    }

    /* Get the grid spacing                                          */
    dh = sqrt((double)window->npe[c] / (2.0 * window->cellarea[c]));

    /* Total 3D speed                                                */
    spd = 2.0 * iws + cspd;
    /* Get the gravity wave + 2D current speed                       */
    spd2d = 2.0 * sqrt(-wincon->g * window->botz[c]) +
            sqrt(windat->uav[c] * windat->uav[c] +
                 windat->vav[c] * windat->vav[c]);

    /* Get the cfl conditions                                        */
    cfl2d = (spd2d > 0.0) ? 1.0 / (dh * spd2d) : eps;
    cfl3d = (spd > 0.0) ? 1.0 / (dh * spd) : eps;
    if (traf) {
      if (windat->cour[c] && windat->cour[c] < windat->mcfl2d)
	windat->mcfl2d = windat->cour[c];
      if (windat->lips[c] && windat->lips[c] < windat->mcfl3d) {
	windat->mcfl3d = windat->lips[c];
	windat->cflc = c;
      }
    } else {
      if (cfl2d < windat->mcfl2d)
	windat->mcfl2d = cfl2d;
      if (cfl3d < windat->mcfl3d) {
	windat->mcfl3d = cfl3d;
	windat->cflc = c;
      }
    }
    windat->cfl2d[c] = cfl2d;
    windat->cfl3d[c] = cfl3d;
  }

  /* Get the Courant number for vertical advection */
  if (wincon->cfl & WVEL) {
    for (cc = 1; cc <= wincon->vc; cc++) {
      c = wincon->s1[cc];
      zp1 = window->zp1[c];
      lc = window->m2d[c];
      spd = (cc <= wincon->vcs) ? 0.5 * (windat->wtop[lc] + windat->w[c]) :
        0.5 * (windat->w[c] + windat->w[zp1]);
      if (fabs(spd) < minval)
        spd = minval;
      cfl3d = fabs(wincon->dz[c] * wincon->Ds[lc] / spd);
      if (!traf && cfl3d < windat->mcfl3d) {
        windat->mcfl3d = cfl3d;
        windat->cflc = c;
      }
    }
  }
}

/* END calc_cfl()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculate the minimum resolution for a given run time ratio and   */
/* execution speed.                                                  */
/*-------------------------------------------------------------------*/
double calc_min_cfl_3d(geometry_t *window, /* Window geometry        */
		       window_t *windat,   /* Window data            */
		       win_priv_t *wincon, /* Window constants       */
		       int *cl             /* Location of maximum    */
  )
{
  int c, cc, cs, cb;            /* Centre coordinates                */
  int e, ee, es;                /* Edge coordinates                  */
  int j, zp1, zm1;              /* Sparse maps                       */
  double minval = 1e-10;        /* Minimum value for velocity        */
  double spd, iws, cspd;        /* Current / wave speeds             */
  double pi = 3.14159;          /* Value of pi                       */
  double *dh = wincon->d1;      /* Grid spacing                      */
  double *ws = wincon->d2;      /* 2 layer internal wave speed       */
  double cfl = HUGE;            /* Miniumum CFL                      */
  double N2, N;                 /* Buoyancy frequency                */
  double hs, vs;                /* Horizontal and vertical stability */
  double dz;                    /* Layer thickness                   */
  int *sur = wincon->i3;        /* Surface at beginning of timestep  */
  int *bot = wincon->i1;        /* Bottom coordinate                 */

  /* Get the grid spacing                                            */
  for (cc = 1; cc <= window->v2_t; cc++) {
    cs = window->w2_t[cc];

    dh[cs] = 0.0;
    for (j = 1; j <= window->npe[cs]; j++) {
      e = window->c2e[j][c];
      if (cs == window->e2c[e][0]) {
	if (window->h1acell[e])
	  dh[cs] += 1.0 / (geom->h1acell[e] * geom->h1acell[e]);
      }
    }
    dh[cs] = sqrt(dh[cs]);
    ws[cs] = int_wave_spd(window, windat, wincon, cs, cc);
  }

  /* Get the 3D CFL condition                                        */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = sur[cc];
    cb = bot[cc];
    cs = window->m2d[c];

    /* Loop for the surface position at the beginning of the time    */
    /* step to the bottom. wincon->i1 and i3 were set in set_dz().   */
    while (c != window->zm1[c]) {
      zp1 = window->zp1[c];
      zm1 = window->zm1[c];

      /* Get the Brunt Vaisala frequency                             */
      /* N2 = -(g/rho)(drho/dz) (s-1)                                */
      N2 = buoyancy_frequency2(window, windat, wincon, c);

      N = N2 > 0 ? sqrt(N2) : 0.0;
      /* Get the mode 1 internal wave speed                          */
      /* c=NH/(n.pi), n=1,2,3... (Gill, 1982, p159) for long waves   */
      /* constant N.                                                 */
      iws = N * (windat->eta[cs] - window->botz[cs]) / pi;
      /* Use a 2 layer approximation for shallow water internal waves*/
      iws = ws[cs];

      /* Get the current speed                                       */
      cspd = sqrt(windat->u[c] * windat->u[c] +
		  windat->v[c] * windat->v[c]);

      /* Get the CFL condition                                       */
      spd = max(minval, 2.0 * iws + cspd);
      hs = 1.0 / (spd * dh[cs]);

      /* Get the vertical Courant condition                          */
      spd = (c == sur[cc]) ? 0.5 * (windat->wtop[cs] + windat->w[c]) :
	0.5 * (windat->w[c] + windat->w[zp1]);
      spd = max(fabs(spd), minval);
      dz = max(wincon->hmin, wincon->dz[c] * wincon->Ds[cs]);
      vs = min(hs, fabs(dz / spd));
      if (vs < cfl) {
	cfl = vs;
	*cl = c;
      }
      c = window->zm1[c];
    }
  }
  return(cfl);
}

/* END calc_min_cfl_3d()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computs a normalized vertical profile of a tracer. Note: this     */
/* differs to the tracerstats version in that a no-gradient is set   */
/* above the surface, whereas all tracerstats routines only perform  */
/* operations on wet cells.                                          */
/*-------------------------------------------------------------------*/
void nor_vert_prof(geometry_t *window,       /* Window geometry       */
		   window_t *windat,         /* Window data           */
		   win_priv_t *wincon        /* Window constants      */
		   )
{
  int c, cc, cs;
  double trs;

  for (cc = 1; cc <= wincon->vcs2; cc++) {
    c = wincon->s2[cc];
    cs = c;
    trs = windat->tr_wc[wincon->nprof][cs];
    while (c != window->zm1[c]) {
      windat->nprof[c] = (trs) ? windat->tr_wc[wincon->nprof][c] / trs : 0.0;
      c = window->zm1[c];
    }
    trs = windat->nprof[cs];
    c = cs;
    while(c != window->zp1[c]) {
      c = window->zp1[c];
      windat->nprof[c] = trs;
    } 
  }
  if (wincon->nprof2d >= 0) {
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      cs = window->m2d[c];
      windat->nprof[c] *= windat->tr_wcS[wincon->nprof2d][cs];
    }
  }
}

/* END nor_vert_prof()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates diagnostic numbers and saves to tracers if required    */
/*-------------------------------------------------------------------*/
void diag_numbers(geometry_t *window,       /* Window geometry       */
		  window_t *windat,         /* Window data           */
		  win_priv_t *wincon        /* Window constants      */
		  )
{
  int c, cc, cs, cn;
  int e, e1, ee, es;
  int zm1;
  double depth;  /* Water depth                                      */
  double dzface; /* Cell thickness between cell centers              */
  double speed;  /* Current speed                                    */
  double N, N2;  /* Brunt Vaisala (buoyancy) frequency (s-1)         */
  double iws;    /* Internal wave phase speed (ms-1)                 */
  double Ri;     /* Gradient Richardson number                       */
  double Rf;     /* Flux Richardson number                           */
  double Re;     /* Reynolds number                                  */
  double Fr;     /* Froude number                                    */
  double Roe;    /* External Rossby radius (m)                       */
  double Roi;    /* Internal Rossby radius (m)                       */
  double du1;    /* e1 vertical shear                                */
  double du2;    /* e2 vertical shear                                */
  double me;     /* Mechanical energy (Jm-3)                         */
  double ke;     /* Kinetic energy (Jm-3)                         */
  double pi = 3.14159;          /* Value of pi */
  double co = 1493.0;           /* Constants for speed of sound */
  double ao =  3.0;
  double bo = -0.006;
  double go = -0.04;
  double d0 = 1.2;
  double eo = -0.01;
  double ho = 0.0164;
  double To = 10.0;
  double T1 = 18.0;
  double So = 35.0;
  double *P = wincon->w8;
  double d1, d2;

  /* Cell centered horizontal viscosity                              */
  get_u1vhc(window, windat, wincon);

  /* Tidal elevation from TPXO model                                 */
  if (windat->tpxotide) {
    for (cc = 1; cc <= window->b2_t; cc++) {
      double gmt = windat->t - wincon->tz;
      double ramp = (wincon->rampf & (TIDALH|TIDALC)) ? 	
	windat->rampval : 1.0;
      c = window->w2_t[cc];
      windat->tpxotide[c] = ramp * csr_tide_eval(&wincon->tc, c, gmt);
    }
  }
  if (wincon->numbers1 & (TPXOV|TPXOT)) {
    double gmt = windat->t - wincon->tz;
    double u1, u2, t1, t2;
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      u1 = u2 = d1 = 0.0;
      for (ee = 1; ee <= window->npe[c]; ee++) {
	e = window->c2e[ee][c];
	d2 = window->h1au1[e];
	u1 += d2 * csr_tide_eval(&wincon->tcu, e, gmt);
	u2 += d2 * csr_tide_eval(&wincon->tcv, e, gmt);
	t1 = u1 * window->botzu1[e];
	t2 = u2 * window->botzu1[e];
	d1 += d2;
      }
      if (windat->tpxovelu) windat->tpxovelu[c] = u1 / d1;
      if (windat->tpxovelv) windat->tpxovelv[c] = u1 / d1;
      if (windat->tpxotranu) windat->tpxotranu[c] = t1 / d1;
      if (windat->tpxotranv) windat->tpxotranv[c] = t2 / d1;
    }
  }
  if (wincon->numbers1 & TRAN2D) {
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      windat->uat[c] = windat->vat[c] = d1 = 0.0;
      for (ee = 1; ee <= window->npe[c]; ee++) {
	e = window->c2e[ee][c];
	d2 = 0.5 * window->h1au1[es] * window->h2au1[es];
	windat->uat[c] += (windat->u1av[e] * window->costhu1[e] + 
			   windat->u2av[e] * window->costhu2[e]) *
	  d2 * window->botzu1[e];
	windat->vat[c] += (windat->u1av[e] * window->sinthu1[e] + 
			   windat->u2av[e] * window->sinthu2[e]) *
	  d2 * window->botzu1[e];
	d1 += d2;
      }
      windat->uat[c] /= d1;
      windat->vat[c] /= d1;
    }
  }

  /* Get the cell centered wind speed and direction                    */
  if (windat->windcs && windat->windcd) {
    double *wu = wincon->d1;
    double *wv = wincon->d2;
    double windx, windy;

    vel_cen(window, windat, wincon, windat->wind1, NULL, wu, wv, NULL, NULL, 1);
    for (cc = 1; cc <= window->b2_t; ++cc) {
      c = window->w2_t[cc];
      windx = wu[c];
      windy = wv[c];
      stresswind(&windx, &windy, 10.0, 26.0, 0.00114, 0.00218);
      windat->windcs[c] = sqrt(windx * windx + windy * windy);
      windat->windcd[c] = 0.0;
      if (windat->windcs[c] > 0.0) {
	windat->windcd[c] = acos(windy / windat->windcs[c]);
	if (windx < 0)
	    windat->windcd[c] *= -1.0;
	windat->windcd[c] *= 180.0/M_PI;
      }
    }
  }

  /* Get the pressure (Bar) for sound calculations */
  if (windat->sound) {
    double cf = 1e-5;
    for (cc = 1; cc <= window->b2_t; cc++) {
      cs = c = window->w2_t[cc];
      P[c] = cf * windat->dens[c] * g * (windat->eta[c] - window->cellz[c]);
      while (c != window->zm1[c]) {
	c = window->zm1[c];
	P[c] -= cf * windat->dens[c] * g * window->cellz[c];
      }
    }
  }

  if (windat->speed_sq) {
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->bot_t[cc];
      cs = window->m2d[c];
      windat->speed_sq[cs] = sqrt(windat->u[c] * windat->u[c] +
				  windat->v[c] * windat->v[c]);
    }
  }

  if (windat->tau_bm) {
    memcpy(wincon->w8, windat->tau_be1, window->szeS * sizeof(double));
    vel_cen(window, windat, wincon, wincon->w8, NULL, windat->tau_be1, windat->tau_be2,
	    windat->tau_bm, NULL, 1);
  }

  if (windat->slope_x) {
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      windat->slope_x[c] = 0.0;
      for (ee = 1; ee <= window->npe[c]; ee++) {
	e = window->c2e[ee][c];
	cn = window->c2c[ee][c];
	windat->slope_x[c] = window->eSc[ee][c] *
	  (windat->eta[cn] - windat->eta[c]) /
	  window->h2au1[e];
      }
      windat->slope_x[c] /= (double)window->npe[c];
    }
  }

  if (windat->surfz) {
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      windat->surfz[c] = NaN;
    }
    for (cc = 1; cc <= wincon->vcs2; cc++) {
      c = wincon->s2[cc];
      cs = window->m2d[c];
      windat->surfz[cs] = window->s2k[c];
    }
  }
  if (windat->wetcell) {
    double bot, dd = DRY_FRAC * wincon->hmin;
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      windat->wetcell[c] = 100.0;
    }
    for (cc = 1; cc <= wincon->vcs2; cc++) {
      c = window->m2d[wincon->s2[cc]];
      double bot = dd + window->botz[c];
      windat->wetcell[c] = 100.0 * (1.0 - (windat->eta[c] - bot) / fabs(bot));
    }
  }

  if (windat->speed_2d || windat->tfront) {
    for (ee = 1; ee <= window->b3_e1; ee++) {
      e = window->w2_e1[ee];
	wincon->d1[e] = sqrt(windat->u1av[e] * windat->u1av[e] +
			     windat->u2av[e] * windat->u2av[e]);
    }
    /* Use natural neighbours to interpolate edge speed onto the     */
    /* cell centre.                                                  */
    /*
    interp_edge_us(window, wincon->d1, window->b2_t, window->w2_t, 
		   window->cellx, window->celly, windat->speed_2d, "nn_non_sibson");
    */
    /* Use area weighted average to interpolate edge speed onto the  */
    /* cell centre.                                                  */
    memset(windat->speed_2d, 0, window->szcS * sizeof(double));
    for (cc = 1; cc <= window->b2_t; cc++) {
      double sum = 0.0;
      c = window->w2_t[cc];
      for (ee = 1; ee <= window->npe[cs]; ee++) {
	e = window->c2e[ee][c];
	windat->speed_2d[c] += window->edgearea[e] * wincon->d1[e];
	sum += window->edgearea[e];
      }
      windat->speed_2d[c] /= sum;;
    }
  }

  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    cs = window->m2d[c];
    zm1 = window->zm1[c];

    /* Get the cell thickness between cell faces                     */
    dzface = 0.5 * (wincon->dz[zm1] + wincon->dz[c]);
    depth = windat->eta[cs] - window->botz[cs];
    speed = sqrt(windat->u[c] * windat->u[c] +
		 windat->v[c] * windat->v[c]);
    if (windat->speed_3d) windat->speed_3d[c] = speed;
    
    /* Get the Brunt Vaisala frequency                               */
    /* N2 = -(g/rho)(drho/dz) (s-1)                                  */
    N2 = buoyancy_frequency2(window, windat, wincon, c);
    N = N2 > 0 ? sqrt(N2) : 0.0;
    if (windat->brunt) windat->brunt[c] = N;

    /* Simpson-Hunter tidal front                                    */
    /* Hearn (1985) On the value of the mixing efficiency in the     */
    /* Simpson-Hunter h/u3 criterion. eps=3e-3                       */
    /* Pingree and Griffiths (1978) Tidal Fronts on the Shelf Seas   */
    /* Around the British Isles. JGR, 83 C9.                         */
    /* S=log10(H/Cd.U^3) = 1.5 [cgs]. S > 2 stratified & S <1 mixed. */
    if (windat->tfront && c == cs) {
      d1 = depth / (wincon->Cd[c] * pow(windat->speed_2d[cs], 3.0));
      windat->tfront[c] = log10(d1);
    }

    /* Get the mode 1 internal wave speed                            */
    /* c=NH/(n.pi), n=1,2,3... (Gill, 1982, p159) for long waves     */
    /* constant N.                                                   */
    iws = N * depth / pi;
    if (windat->int_wave) windat->int_wave[c] = iws;

    /* Vertical u gradient                                           */
    du1 = windat->u[c] - windat->u[zm1];

    /* Vertical u2 gradient                                          */
    du2 = windat->v[c] - windat->v[zm1];

    /* Vertical shear magnitude                                      */
    if (windat->shear_v)
      windat->shear_v[c] = sqrt((du1 / dzface) * (du1 / dzface) + 
				(du2 / dzface) * (du2 / dzface));

    /* Gradient Richardson number (Apel, 1987, p 239)                */
    /* Ri > 0 : stable                                               */
    /* Ri = 0 : neutral                                              */
    /* Ri < 0 : unstable                                             */
    Ri = (du1 / dzface) * (du1 / dzface) + (du2 / dzface) * (du2 / dzface);
    Ri = (Ri > 0.0) ? N2 / Ri : 0.0;
    if (windat->rich_gr) windat->rich_gr[c] = Ri;

    /* Flux Richardson number                                        */
    Rf = windat->Kz[c] * Ri / windat->Vz[c];
    if (windat->rich_fl) windat->rich_fl[c] = Rf;

    /* Reynolds number (vertical eddy viscosity, e.g. Pond and       */
    /* Pickard (1983) p 104).                                        */
    /* Re > 10e5 : turbulent                                         */
    /* Re < 2e3  : laminar                                           */
    Re = speed * wincon->dz[c] / windat->Vz[c];
    if (windat->reynolds) windat->reynolds[c] = Re;

    /* Interfacial Froude Number (Dyer (1997) p42)                   */
    /* Fr^2 = U^2 / c^2 (c=internal wave speed).                     */
    /* Also quoted as : Fr=U/NH, or Fr^2=1/Ri^2                      */
    /* Fr < 1 : subcritial                                           */
    /* Fr = 1 : critial (0.33 for continuous stratification, Dyer)   */
    /* Fr > 1 : supercritial                                         */
    Fr = (iws > 0.0) ? speed / iws : 0.0;
    if (windat->froude) windat->froude[c] = Fr;

    /* Sigma_t */
    if (windat->sigma_t) windat->sigma_t[c] = windat->dens_0[c] - 1000.0;

    /* Internal Rossby radius (m) (Gill, 1982, p 207, 159)           */
    Roi = iws / fabs(vertex_mean(window, windat->fv, cs));
    if (windat->rossby_in) windat->rossby_in[c] = Roi;

    /* External Rossby radius (m)                                    */
    Roe = sqrt(depth * wincon->g) / fabs(vertex_mean(window, windat->fv, cs));
    if (windat->rossby_ex) windat->rossby_ex[cs] = Roe;

    /* Mechanical energy (Jm-3) (Kowalik and Murty (1993) p6, Eqn.   */
    /* 1.24.                                                         */
    me = windat->dens[c] * (wincon->g * windat->eta[cs] +
			    0.5 * (windat->u[c] * windat->u[c] +
				   windat->v[c] * windat->v[c] +
				   windat->w[c] * windat->w[c]));
    if (windat->energy) windat->energy[c] = me;
    ke = windat->dens[c] * 0.5 * (windat->u[c] * windat->u[c] +
				  windat->v[c] * windat->v[c] +
				  windat->w[c] * windat->w[c]);
    if (windat->kenergy) windat->kenergy[c] = ke;

    /* Speed of sound (ms-1) (Apel, 1987, p 349)                     */
    if (windat->sound) {
      /*UR until the various speed of sound calculations are verified we
       * force the mackenzie calculation providing depth
       */
      double depthz = windat->eta[cs] - window->cellz[c];
      windat->sound[c] = mackenzie_speed_of_sound(windat->temp[c], windat->sal[c],
					depthz, 1, g);
   /*   * /
      windat->sound[c] = speed_of_sound(windat->temp[c], windat->sal[c],
					P[c], 0, g);
*/
      /*
      windat->sound[c] = co + ao * (windat->temp[c] - To) +
	bo * (windat->temp[c] - To) * (windat->temp[c] - To) +
	go * (windat->temp[c] - T1) * (windat->temp[c] - T1) +
	d0 * (windat->sal[c] - So) +
	eo * (windat->temp[c] - T1) * (windat->sal[c] - So) +
	ho * fabs(windat->eta[cs] - window->cellz[c]);
      */
    }
  }
  if (windat->sound)
    sound_channel(window, windat, wincon);
}

/* END diag_numbers()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computs the cell centered horizontal viscosity diagnostic         */
/*-------------------------------------------------------------------*/
void get_u1vhc(geometry_t *window,   /* Window geometry              */
	       window_t *windat,     /* Window data                  */
	       win_priv_t *wincon    /* Window constants             */
	       )
{
  int cc, c, cs;
  int ee, e, es;
  double d1, d2;

  if (windat->u1vhc) {
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      cs = window->m2d[c];
      windat->u1vhc[c] = 0.0;
      d1 = 0.0;
      for (ee = 1; ee <= window->npe[cs]; ee++) {
	e = window->c2e[ee][c];
	es = window->m2de[e];
	d2 = window->edgearea[es];
	if (wincon->diff_scale & NONE)  d2 = 1.0;
	if (wincon->diff_scale & LINEAR) d2 = sqrt(d2);
	if (wincon->diff_scale & CUBIC) d2 = d2 * sqrt(d2);
	if (wincon->diff_scale & SCALEBI)
	  windat->u1vhc[c] += (d2 * wincon->u1vh[e] / (0.125 * window->edgearea[es]));
	else
	  windat->u1vhc[c] += (d2 * wincon->u1vh[e]);
	/*windat->u1vhc[c] += (d2 * wincon->u2kh[e]);*/
	d1 += d2;
      }
      windat->u1vhc[c] /= d1;
    }
  }
}

/* END get_u1vhc()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computs the surface and bottom Ekman pumping                      */
/*-------------------------------------------------------------------*/
void ekman_pump_e1(geometry_t *window,   /* Window geometry          */
		   window_t *windat,     /* Window data              */
		   win_priv_t *wincon,   /* Window constants         */
		   double *taus,
		   double *taub
		   )
{
  int cc, c, e, j;
  double *dutdn = wincon->d1;
  double *dundt = wincon->d2;
  double *d1 = wincon->d3;

  /* Initialize                                                      */
  memset(windat->sep, 0, window->szcS * sizeof(double));
  memset(windat->bep, 0, window->szcS * sizeof(double));

  /* Gradients of surface stress                                     */
  vel_grad(window, windat, wincon, taus, NULL, d1, dutdn, GRAD_2D|GRAD_NOR);
  vel_grad(window, windat, wincon, taus, NULL, dundt, d1, GRAD_2D|GRAD_TAN);

  /* Compute wind stress curl and surface Ekman pumping              */
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->w2_t[cc];
    for (j = 1; j<= window->npe[c]; j++) {
      e = window->c2e[j][c];
      windat->sep[c] += (dutdn[e] - dundt[e]) * wincon->coriolis[c] *
	windat->dens[c];
    }
    windat->sep[c] /= (double)window->npe[c];
  }

  /* Gradients of bottom stress                                      */
  vel_grad(window, windat, wincon, taub, NULL, d1, dutdn, GRAD_2D|GRAD_NOR);
  vel_grad(window, windat, wincon, taub, NULL, dundt, d1, GRAD_2D|GRAD_TAN);

  /* Compute bottom stress curl and surface Ekman pumping            */
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->w2_t[cc];
    for (j = 1; j<= window->npe[c]; j++) {
      e = window->c2e[j][c];
      windat->bep[c] += (dutdn[e] - dundt[e]) * wincon->coriolis[c] *
	windat->dens[c];
    }
    windat->bep[c] /= (double)window->npe[c];
  }
}

/* END ekman_pump_e1()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computs the depth of sound channels                               */
/*-------------------------------------------------------------------*/
void sound_channel(geometry_t *window,  /* Window geometry           */
		   window_t *windat,    /* Window data               */
		   win_priv_t *wincon   /* Window constants          */
		   )
{
  int c, cc, cs;
  int zm1, zp1, zm2, zp2;
  double cz, db, df, face, grad, curv;
  double *tr;
  double *Fz = wincon->w4;
  double *depth = wincon->w8;
  double *dSdz = wincon->w9;
  double *ss = wincon->w5;
  double k[5];
  int *loc = wincon->i4;
  int n, npass = 3;

  tr = ss;
  if (windat->sound && windat->schan && windat->sonic) {
    memset(windat->schan, 0, window->szc * sizeof(double));
    memset(windat->sonic, 0, window->szcS * sizeof(double));
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      /*windat->schan[c] = NaN;*/
      windat->schan[c] = 0.0;
    }
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      /*windat->sonic[c] = NaN;*/
      windat->sonic[c] = 0.0;
      c = window->bot_t[cc];
      windat->sound[window->zm1[c]] = windat->sound[c];
    }
    /* Smooth the sound profile vertically                           */
    k[0] = k[4] = 0.1;
    k[1] = k[3] = 0.2;
    k[2] = 0.4;
    /*k[0] = k[1] = k[2] = k[3] = k[4] = 0.2;*/
    memcpy(Fz, windat->sound, window->szc * sizeof(double));
    for (n = 0; n < npass; n++) {
      for (cc = 1; cc <= window->b3_t; cc++) {
	c = window->w3_t[cc];
	zm1 = window->zm1[c];
	zp1 = window->zp1[c];
	zm2 = window->zm1[zm1];
	zp2 = window->zp1[zp1];
	ss[c] = k[0] * Fz[zp2] +
	k[1] * Fz[zp1] +
	k[2] * Fz[c] +
	k[3] * Fz[zm1] +
	k[4] * Fz[zm2];
      }
      for (cc = 1; cc <= window->b2_t; cc++) {
	c = window->bot_t[cc];
	ss[window->zm1[c]] = ss[c];
      }
      memcpy(Fz, ss, window->szc * sizeof(double));
    }
    memcpy(windat->sound, ss, window->szc * sizeof(double));

    memset(Fz, 0, window->szc * sizeof(double));
    /* Set the bottom boundary condition                             */
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->bot_t[cc];
      tr[window->zm1[c]] = tr[c];
    }

    /* Get the concentration values at the lower face                */
    for (cc = 1; cc <= wincon->vc; cc++) {
      c = wincon->s2[cc];
      zp1 = window->zp1[c];
      zm1 = window->zm1[c];
      zm2 = window->zm1[zm1];
      order4_do(Fz, tr, c, c, zp1, zm1, zm2);
    }

    /* Set the bottom boundary condition                             */
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->bot_t[cc];
      Fz[window->zm1[c]] = Fz[c];
    }

    /* Get the mean cell spacings and Courant number                 */
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      cs = window->m2d[c];
      zm1 = window->zm1[c];
      wincon->w5[c] = wincon->dz[zm1] / (wincon->dz[zm1] + wincon->dz[c]);
      wincon->w9[c] = 0.5 * (wincon->dz[c] + wincon->dz[zm1]);
      wincon->w8[c] = windat->dt / (wincon->w9[c] * wincon->Ds[cs]);
    }
    
    /* Get the derivative and depth                                  */
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      zp1 = window->zp1[c];
      dSdz[c] = (Fz[zp1] - Fz[c]) / wincon->dz[c];
      depth[c] = window->cellz[c];
    }
    /* Get the surface depth, and set a no gradient condition at the */
    /* surface.                                                      */
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      cs = window->m2d[c];
      loc[cs] = cs;
      depth[c] = 0.5 * (windat->eta[cs] + window->gridz[c]);
      dSdz[c] = dSdz[window->zm1[c]];
      c = window->bot_t[cc];
      dSdz[window->zm1[c]] = dSdz[c];
    }
    /* Get the sonic depth (maximum in sound, or dSdz changes from   */
    /* negative to positive) and place in the top layer.             */
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      zm1 = window->zm1[c];
      cs = window->m2d[c];
      if (dSdz[zm1]) {
	if (isnan(windat->sonic[cs]) && dSdz[c] < 0.0 && dSdz[zm1] > 0.0) {
	  windat->sonic[cs] = depth[c] - dSdz[c] * (depth[zm1] - depth[c]) /
	    (dSdz[zm1] - dSdz[c]);
	}
      }
    }
    /* Get the sound channel depth (minimum in sound, or dSdz        */
    /* changes from positive to negative) and place in layers below  */
    /* the top layer.                                                */
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      zm1 = window->zm1[c];
      cs = window->m2d[c];
      if (dSdz[zm1]) {
	if (dSdz[c] > 0.0 && dSdz[zm1] < 0.0) {
	  windat->schan[loc[cs]] = depth[c] - dSdz[c] * (depth[zm1] - depth[c]) /
	    (dSdz[zm1] - dSdz[c]);
	  loc[cs] = window->zm1[loc[cs]];
	}
      } else {
	if (dSdz[c]) {
	  windat->schan[loc[cs]] = depth[zm1];
	  loc[cs] = window->zm1[loc[cs]];
	}
      }
    }
  }
}

/* END sound_channel()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the Brunt Vaisala (buoyancy) frequency (s-1) */
/*-------------------------------------------------------------------*/
double buoyancy_frequency2(geometry_t *window, /* Window geometry    */
			window_t *windat,   /* Window data           */
			win_priv_t *wincon, /* Window constants      */
			int c               /* Sparse coordinate     */
)
{
  int zm1;
  double drho, N2;

  /* Density gradient                                                */
  zm1 = window->zm1[c];
  drho = windat->dens_0[c] - windat->dens_0[zm1];

  /* Brunt Vaisala frequency squared                                 */
  N2 = -(wincon->g / RHO_0) *
    (drho / (0.5 * (wincon->dz[zm1] + wincon->dz[c])));

  return(N2);
}

double buoyancy_frequency2_m(master_t *master,  /* Master data       */
			     geometry_t *geom,  /* Global geometry   */
			     double *dens,      /* Density           */
			     int c              /* Sparse coordinate */
)
{
  int zm1;
  double drho, N2, dz;

  /* Density gradient                                                */
  zm1 = geom->zm1[c];
  drho = dens[c] - dens[zm1];
  dz = fabs(geom->cellz[c] - geom->cellz[zm1]);

  /* Brunt Vaisala frequency squared                                 */
  N2 = -(master->g / RHO_0) * (drho / dz);

  return(N2);
}

/* END buoyancy_frequency2()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the internal wave speed for shallow water    */
/* internal waves.                                                   */
/* c=sqrt(g'h') where g' = (rho2-rho1)/rho2                          */
/*                    h' = h1.h2/(h1+h2)                             */
/*-------------------------------------------------------------------*/
double int_wave_spd(geometry_t *window, /* Window geometry           */
                    window_t *windat,   /* Window data               */
                    win_priv_t *wincon, /* Window constants          */
                    int c, int cc       /* Sparse coordinate         */
  )
{
  double rho1 = 0.0, rho2 = 0.0;  /* Density of top & bottom layers  */
  double h1, h2;                /* Depths of top and bottom layers   */
  double top_m;                 /* Depth of top of pycnocline        */
  double bot_m;                 /* Depth of bottom of pycnocline     */ 
  int kts, ktb;                 /* k index of pycnocline top/bottom  */
  double spd;                   /* Internal wave speed               */
  double d1, d2;                /* Dummies                           */
  int cs = window->m2d[c];
  int cb = window->bot_t[cc];
  int lc, nc;

  mld(window, windat, wincon, cs, cb, &top_m, &bot_m, &kts, &ktb);
  if (ktb > cb)
    ktb = cb;
  if (kts < cs)
    kts = cs;

  /* Mean density in the surface layer                               */
  nc = 0;
  lc = cs;
  while (lc <= kts) {
    rho1 += windat->dens[lc];
    lc = window->zm1[lc];
    nc++;
  }
  rho1 /= (double)nc;
  /* Mean density in the bottom layer                                */
  nc = 0;
  lc = ktb;
  while (lc <= cb) {
    rho2 += windat->dens[lc];
    lc = window->zm1[lc];
    nc++;
  }
  rho2 /= (double)nc;
  h1 = windat->eta[cs] - 0.5 * (top_m + bot_m);
  h2 = windat->eta[cs] - window->botz[cs] - h1;
  d1 = (rho2 - rho1) * h1 * h2;
  d2 = rho2 * (h1 + h2);
  spd = g * d1 / d2;
  return ((spd > 0) ? sqrt(spd) : 0.0);
}

/* END int_wave_spd()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the internal wave speed for shallow water    */
/* internal waves with continuous stratification (1st mode).         */
/* c=-drho/dz*gH/rho.pi                                              */
/* The advective speed is added to the wave speed.                   */
/*-------------------------------------------------------------------*/
double int_wave_spd_cont_w(geometry_t *window, /* Window geometry    */
			   window_t *windat,   /* Window data        */
			   win_priv_t *wincon, /* Window constants   */
			   int c, int cc       /* Sparse coordinate  */
  )
{
  double spd, d1;               /* Internal wave speed               */
  double drhodz;                /* Vertical density gradient         */
  double pi = 3.14159;          /* Value of pi                       */
  int cb = window->bot_t[cc];
  int zm1, c3 = c;

  spd = 0.0;
  while (c3 < cb) {
    zm1 = window->zm1[c3];
    drhodz = (windat->dens_0[c3] - windat->dens_0[zm1]) / wincon->dz[c3];
    d1 =
      drhodz * wincon->g * window->cellz[c3] / (windat->dens_0[c3] * pi);
    /* d1+=sqrt(windat->u1[c3]*windat->u1[c3]+windat->u2[c3]*windat->u2[c3]);
     */
    if (d1 > spd)
      spd = d1;
    c3 = zm1;
  }
  return ((spd > 0) ? spd : 0.0);
}

/* END int_wave_spd_cont_w()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the internal wave speed for shallow water    */
/* internal waves with continuous stratification (1st mode).         */
/* c=-drho/dz*gH/rho.pi (Gill, 1982, p159).                          */
/* Calculation performed on the master.                              */
/*-------------------------------------------------------------------*/
double int_wave_spd_cont_m(master_t *master,  /* Master data         */
                           int c, int cc      /* Sparse coordinate   */
  )
{
  geometry_t *geom = master->geom;
  double spd, d1;               /* Internal wave speed               */
  double drhodz;                /* Vertical density gradient         */
  double pi = 3.14159;          /* Value of pi                       */
  double g = 9.81;              /* Acceleration due to gravity       */
  double top, bot, dz;          /* Top, bottom, thickness of layers  */
  int zm1, c3, cb;              /* Sparse coordinates                */
  double r1, r2, h1, h2;        /* Mean densities, layer thicknesses */
  double thr = -0.01;           /* Density gradient threshold        */
  int ml = 0;                   /* Flag for pycnocline               */
  int mode = 0;                 /* Type of speed calculation         */

  spd = 0.0;
  c3 = c;
  cb = geom->bot_t[cc];
  /* Get the internal wave speed using the 2 layer approximation     */
  if (mode == 0) {
    h1 = r1 = r2 = h2 = 0.0;
    /* Find the pycnocline and mean densities in each layer          */
    while (c3 < cb) {
      zm1 = geom->zm1[c3];
      top = (c3 == geom->sur_t[cc]) ? 0.0 : geom->gridz[geom->zp1[c3]];
      bot =
        (c3 ==
         geom->bot_t[cc]) ? geom->botz[geom->m2d[c3]] : geom->gridz[c3];
      dz = top - bot;
      drhodz = (master->dens_0[c3] - master->dens_0[zm1]) / dz;
      if (drhodz < thr)
        ml = 1;
      if (ml) {
        h2 += dz;
        r2 += (master->dens_0[c3] * dz);
      } else {
        h1 += dz;
        r1 += (master->dens_0[c3] * dz);
      }
      c3 = zm1;
    }
    r1 /= h1;
    r2 /= h2;
    r1 = master->dens_0[c];
    r2 = master->dens_0[cb];
    spd = sqrt(g * (r2 - r1) * h1 * h2 / (r2 * (h1 + h2)));
  } else {
    while (c3 < cb) {
      zm1 = geom->zm1[c3];
      top = (c3 == geom->sur_t[cc]) ? 0.0 : geom->gridz[geom->zp1[c3]];
      bot =
        (c3 ==
         geom->bot_t[cc]) ? geom->botz[geom->m2d[c3]] : geom->gridz[c3];
      dz = top - bot;
      drhodz = (master->dens_0[c3] - master->dens_0[zm1]) / dz;
      d1 = drhodz * g * geom->cellz[c3] / (master->dens_0[c3] * pi);
      if (d1 > spd)
        spd = d1;
      c3 = zm1;
    }
  }
  return ((spd > 0) ? spd : 0.0);
}

/* END int_wave_spd_cont_m()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to compute global rms error                               */
/*-------------------------------------------------------------------*/
double get_global_rms(double *a1, double *a2, int vc, int *vec)
{
  int c, cc;
  double rms = 0.0;

  for (cc = 1; cc <= vc; cc++) {
    c = vec[cc];
    rms += ((a1[c] - a2[c]) * (a1[c] - a2[c]));
  }
  return(sqrt(rms / (double)vc));
}

/* END get_rms()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to compute global rms error                               */
/*-------------------------------------------------------------------*/
double get_profile_rms(double *a1, double *a2, int cs, int *zm1)
{
  int c = cs;
  double rms = 0.0;
  double n = 0.0;

  while (c != zm1[c]) {
    rms += ((a1[c] - a2[c]) * (a1[c] - a2[c]));
    c = zm1[c];
    n += 1.0;
 }
  return(sqrt(rms / n));
}

/* END get_rms()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set the timesteps to those calculated by the cfl       */
/* diagnostic after a specific time has elapsed.                     */
/*-------------------------------------------------------------------*/
void reset_cfl(master_t *master /* Master data                       */
  )
{
  int dt;
  double sf = 0.9;              /* Safety factor                     */
  double cfl2d;

  if (master->t < master->cfl_dt)
    return;

  dt = (int)(master->mcfl3d * sf) / 5;
  master->grid_dt = 5.0 * (double)dt;

  dt = (int)(master->mcfl2d * sf);
  while (dt > 1 && (dt % 2 != 0 && dt % 3 != 0 && dt % 5 != 0))
    dt--;
  cfl2d = (double)dt;
  if (cfl2d <= 1) {
    hd_warn("2D time-step is less than 1 (%f)\n", cfl2d);
    cfl2d = 1.0;
  }
  master->iratio = (int)(master->grid_dt / cfl2d);
  if (master->iratio % 2 != 0)
    master->iratio++;
  master->grid_dt = (double)master->iratio * cfl2d;
  fprintf(stderr,
          "Timesteps altered to cfl conditions %5.2f (3d) %5.2f (2d) at %8.2f days\n",
          master->grid_dt, cfl2d, master->t / 86400);
  master->cfl = PASSIVE;
}

/* END reset_cfl()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set the timesteps to those calculated by the cfl       */
/* diagnostic after a specific time has elapsed.                     */
/*-------------------------------------------------------------------*/
void reset_cfl_event(master_t *master,  /* Model grid data structure */
		     double mcfl,       /* Minimum CFL time          */
		     int mode           /* Warning print mode        */
  )
{
  int dt;                          /* New time step                  */
  double olddt = master->grid_dt;  /* Old grid time step             */
  double sf = 0.9;                 /* Safety factor                  */
  double max_dt2d;                 /* 2D time step                   */
  int i, j, k, c, cc;              /* Counters                       */
  geometry_t *geom = master->geom; /* Master geometry                */

  if (sf * mcfl < master->grid_dt) {
    max_dt2d = master->grid_dt / master->iratio;
    dt = (int)(mcfl * sf) / 2;
    master->grid_dt = max(2.0 * (double)dt, 2.0 * max_dt2d);
    master->iratio = (int)ceil(master->grid_dt / max_dt2d);
    /* Only even iratio's allowed with the leapfrog scheme so that   */
    /* fluxes summed on odd timesteps add integrally to the 3d step  */
    if (((master->iratio + 1) % 2 == 0))
      master->iratio++;
    /* Print warning if grid_dt was changed using the cfl diagnostic */
    if(mode & CFL_DIAG) {
      c = master->cflc;
      i = geom->s2i[c];
      j = geom->s2j[c];
      k = geom->s2k[c];
      fprintf(stderr,
	      "Timesteps altered to 3D cfl condition %5.2f (3d) %5.2f (2d) at (%3d %3d %2d) %8.2f days\n",
	      master->grid_dt, master->grid_dt / master->iratio, i, j, k,
	      master->t / 86400);
      if (master->t >= master->cfl_dt)
	master->cfl = PASSIVE;
    }
    /* Print warning if grid_dt changed using the alert diagnostic   */
    else if(mode & CFL_ALERT) {
      hd_warn("LEVEL4 alert at %f days : CLF timestep = %6.2f\n",
	      master->t/86400, master->grid_dt);
    }

    /* Reset the horizontal mixing for the new timestep              */
    for (cc = 1; cc <= geom->b2_e1; cc++) {
      c = geom->w2_e1[cc];
      master->u1kh[c] *= (master->grid_dt / olddt);
      master->u1vh[c] *= (master->grid_dt / olddt);
    }
    master->u1kh[0] = master->u2kh[0] = master->u1vh[0] = master->u2vh[0] =
    1;
  }
}

/* END reset_cfl_event()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate steric height relative to a specified level  */
/* of no motion, lnm. If this level is deeper than the bottom depth  */
/* at any location, a geopotential thickness of zero is specified.   */
/* Steric height is in units of metres. Note : steric height =       */
/* geopotential anomaly/acceleration due to gravity : Godfrey and    */
/* Ridgeway, JPO 15, 481-495 with units of m). See Pond & Pickard    */
/* (1983) p66 for geopotential anomaly.                              */
/* Geostropic currents relative to the lnm are calculated via :      */
/* v=g/f.d(steric)/dx                                                */
/*-------------------------------------------------------------------*/
void steric(geometry_t *window, /* Window geometry                   */
            window_t *windat,   /* Window data                       */
            win_priv_t *wincon  /* Window constants                  */
  )
{
  double lnm;             /* Level of no motion (m).                 */
  double pressu, pressl;  /* Pressure at upper and lower levels      */
  double dnu, dnl;        /* Densities at upper and lower depths     */
  double svu, svl;        /* Sigma-t at upper and lower depths       */
  double asvu, asvl;      /* Anomaly of specific volume at depth 1&2 */
  double geo;             /* Geopotential anomaly at each depth      */
  int cc, c, c2, cs;      /* sparse coordinates/counters             */
  int zm1;
  int cl;
  double *dzz;
  double d1;

  /* Set pointers and initialise                                     */
  lnm = wincon->lnm;
  dzz = wincon->w9;
  memset(windat->steric, 0, window->szcS * sizeof(double));
  memset(dzz, 0, window->szc * sizeof(double));

  /* Get the cell thickness between layers                           */
  /* Surface layer                                                   */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = wincon->s1[cc];
    c2 = window->m2d[c];
    dzz[c] = 0.5 * wincon->dz[c] * wincon->Ds[c2];
  }
  /* Subsurface layers                                               */
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];

    c2 = window->m2d[c];
    zm1 = window->zm1[c];
    dzz[zm1] = 0.5 * (wincon->dz[c] + wincon->dz[zm1]) * wincon->Ds[c2];
  }

  /* Find the nearest level less than lnm.                           */
  /* If lnm > bottom depth then set zero geopotential.               */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = cs = wincon->s1[cc];
    c2 = window->m2d[c];
    cl = ztoc_w(window, c, lnm);

    cl = c;
    geo = 0.0;
    while (cl != window->zm1[cl]) {
      geo += (wincon->dz[cl] * wincon->Ds[c2]);
      cl = window->zm1[cl];
      if (geo >= lnm)
        break;
    }
    cl = window->zp1[cl];

    if (window->botz[c2] > -lnm)
      continue;

    /* Get anomaly specific volume for the surface                   */
    pressu =
      windat->patm[c2] +
      wincon->g * windat->dens[c] * dzz[c] * wincon->Ds[c2];
    dnu = windat->dens[c];
    eos2(35.0, 0.0, pressu, &svu, &d1);
    asvu = 1.0 / dnu - 1.0 / svu;
    geo = pressu * asvu;

    /* Get geopotential anomalies for depth levels 2 to cl           */
    pressl = pressu;
    c = window->zm1[c];

    while (c <= cl) {
      pressl += wincon->g * windat->dens[c] * dzz[c] * wincon->Ds[c2];
      dnl = windat->dens[c];
      eos2(35.0, 0.0, pressl, &svl, &d1);
      asvl = 1.0 / dnl - 1.0 / svl;
      geo += (pressl - pressu) * (asvl + asvu) / 2.0;
      pressu = pressl;
      asvu = asvl;
      c = window->zm1[c];
    }

    /* Get geopotential anomaly for depth level cl                   */
    zm1 = window->zm1[c];
    /* pressl = pressure at layer face cl                            */
    pressl +=
      wincon->g * windat->dens[c] * 0.5 * wincon->dz[c] * wincon->Ds[c2];
    /* pressl = pressure at depth lnm                                */
    pressl += wincon->g * windat->dens[zm1] * (lnm + window->gridz[c]);
    dnl = windat->dens[zm1];
    eos2(35.0, 0.0, pressl, &svl, &d1);
    asvl = 1.0 / dnl - 1.0 / svl;
    geo += (pressl - pressu) * (asvl + asvu) / 2.0;

    /* Convert the integrated geopotential anomaly to steric height  */
    windat->steric[c2] = geo / wincon->g;
  }
}

/* END steric()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the vertical component of vorticity.         */
/* Planetary = f (s-1)                                               */
/* Relative : rv = du2/de1 - du1/de2 (s-1)                           */
/* Absolute : av = (f + rv) (s-1)                                    */
/* Potential : pv = (f+rv)/(H+eta) (m-1s-1)                          */
/* Note: conservation of potential vorticity states pv = constant    */
/* See Apel, J. R. (1987) Principles of Ocean Physics, p271-282.     */
/*-------------------------------------------------------------------*/
void vorticity(geometry_t *window,  /* Window geometry               */
               window_t *windat,    /* Window data                   */
               win_priv_t *wincon   /* Window constants              */
  )
{
  int c, cc, c2;                /* Cell coordinate/counter           */
  int e, ee, e1, e2, e3, e4;    /* Edge coordinates                  */
  int n, j, eoe;                /* Counters                          */
  double *rv;                   /* Relative vorticity                */
  double *prv;                  /* d(rv)/dt                          */
  double *ew1, *ns1;            /* EW and NS gradients               */
  double *ew2, *ns2;            /* EW and NS gradients               */
  double *cor_grad;             /* NS Coriolis gradient              */

  ew1 = wincon->d6;
  ns1 = wincon->d7;
  ew2 = wincon->d1;
  ns2 = wincon->d2;
  cor_grad = wincon->d3;
  prv = windat->rv_drvdt;

  /* Relative, planetary, absolute and potential vorticity are       */
  /* computed and saved in nonlin_coriolis_2d().                     */

  /* Get the vorticity tendencies                                    */
  if (wincon->vorticity & TENDENCY) {
    int c1, c2, zm1;
    double *adv = windat->rv_nonlin;
    double *beta = windat->rv_beta;
    double *stretch = windat->rv_strch;
    double *jebar = windat->rv_jebar;
    double *wsc = windat->rv_wsc;
    double *bsc = windat->rv_bsc;
    double *depth = wincon->d1;
    double *chi = wincon->d2;
    double f;                   /* Coriolis parameter                */
    double dHde1, dHde2;        /* Bottom gradients                  */
    double rho0 = 1024;         /* Standard reference density        */

    memset(chi, 0, window->szcS * sizeof(double));
    for (cc = 1; cc <= window->n2_t; cc++) {
      c = window->w2_t[cc];
      c2 = window->m2d[c];
      zm1 = window->zm1[c];
      /* Depth                                                       */
      depth[c2] = windat->eta[c2] - window->botz[c2];
      jebar[c2] = 1.0 / depth[c2];

      /* Anomaly of potential energy                                 */
      do {
        chi[c2] +=
          wincon->g * window->cellz[c] * windat->dens[c] * wincon->dz[c] /
          rho0;
        c = zm1;
        zm1 = window->zm1[c];
      } while (c != zm1);
    }

    /* Get the stress components divided by depth                    */
    for (ee = 1; ee <= window->n2_e1; ee++) {
      e = window->w2_e1[ee];
      c1 = window->e2c[e][0];
      c2 = window->e2c[e][1];
      bsc[e] = 2.0 * windat->rv_bsc[e] / (depth[c1] + depth[c2]);
      wsc[e] = 2.0 * windat->wind1[e] / (depth[c1] + depth[c2]);
    }

    /* Coriolis gradient                                             */
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      ew2[c] = vertex_mean(window, windat->fv, c);
    }
    tra_grad(window, windat, wincon, ew2, ew1, beta, GRAD_2D);

    /* Depth gradients                                               */
    tra_grad(window, windat, wincon, depth, ew1, ns1, GRAD_2D);

    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = window->m2d[wincon->s1[cc]];

      /* Cell centered Coriolis                                      */
      f = vertex_mean(window, windat->fv, c);

      /* 1. Rate of change of relative vorticity                     */
      prv[c] = (windat->rv[c] - prv[c]) / windat->dt;

      /* 2. Beta                                                     */
      beta[c] = 2.0 * beta[c] * windat->vav[c];

      /* 3. Stretch                                                  */
      stretch[c] = (windat->rv[c] + f) * 
	(windat->detadt[c] + windat->uav[c] * ew1[c] +
	                     windat->vav[c] * ns1[c]) / depth[c];
    }

    /* JEBAR gradient                                                */
    tra_grad(window, windat, wincon, jebar, ew1, ns1, GRAD_2D);
    /* Anomaly of potential energy gradient                          */
    tra_grad(window, windat, wincon, chi, ew2, ns2, GRAD_2D);
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = window->m2d[wincon->s1[cc]];
      /* 4. JEBAR                                                    */
      jebar[c] = ns1[c] * ew2[c] - ew1[c] * ns2[c];
    }

    /* Gradient of (bottom stress) / depth normal to edges           */
    vel_grad(window, windat, wincon, bsc, NULL, ew1, ns1, GRAD_2D|GRAD_NOR);
    /* Center the bottom stress gradient                             */
    vel_cen(window, windat, wincon, ew1, NULL, ew2, ns2, NULL, NULL, 1);
    /* 5. BSC                                                        */
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = window->m2d[wincon->s1[cc]];
      bsc[c] = (ns2[c] - ew2[c]) / rho0;
    }

    /* Gradient of (wind stress) / depth normal to edges             */
    vel_grad(window, windat, wincon, wsc, NULL, ew1, ns1, GRAD_2D|GRAD_NOR);
    /* Center the wind stress gradient                               */
    vel_cen(window, windat, wincon, ew1, NULL, ew2, ns2, NULL, NULL, 1);
    /* 6. WSC                                                        */
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = window->m2d[wincon->s1[cc]];
      wsc[c] = (ns2[c] - ew2[c]) / rho0;
      /* 7. Nonlinear advection and diffusion */
      adv[c] = stretch[c] + jebar[c] + wsc[c] + bsc[c] - prv[c] - beta[c];
    }
  }
}

/* END vorticity()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set up the totals output timeseries file.              */
/*-------------------------------------------------------------------*/
void init_totals(master_t *master,    /* Master data                 */
		 geometry_t **window, /* Window geometry             */
		 win_priv_t **wincon  /* Window constants            */
  )
{
  geometry_t *geom = master->geom;
  int n, m, nn, trn;

  /*---------------------------------------------------------------*/
  /* Initialize the output file                                    */
  if (master->totals && tstotal.fp == NULL) {
    if (geom->nobc)
      master->vf = d_alloc_1d(geom->nobc);
    memset(master->vf, 0, geom->nobc * sizeof(double));
    for (n = 1; n <= master->nwindows; n++) {
      window_t *windat = window[n]->windat;
      if (window[n]->nobc)
	windat->vf = d_alloc_1d(window[n]->nobc);
    }
    if (strlen(master->opath))
      sprintf(tstotal.pname, "%stotals.ts", master->opath);
    else
      sprintf(tstotal.pname, "totals.ts");
    if ((tstotal.fp = fopen(tstotal.pname, "w")) == NULL)
      hd_quit("init_totals: Can't open file %s\n", tstotal.pname);
    tstotal.tsdt = master->totals_dt;
    tstotal.master = master;
    tstotal.i = 1;
    tstotal.j = 1;
    tstotal.k = 1;
    tstotal.tsout = master->t;
    fprintf(tstotal.fp, "## COLUMNS %d\n", 3 + master->ntot + geom->nobc);
    fprintf(tstotal.fp, "##\n");
    fprintf(tstotal.fp, "## COLUMN1.name  Time\n");
    fprintf(tstotal.fp, "## COLUMN1.long_name  Time\n");
    fprintf(tstotal.fp,
	    "## COLUMN1.units  days since 1990-01-01 00:00:00 +10\n");
    fprintf(tstotal.fp, "## COLUMN1.missing_value -999\n");
    fprintf(tstotal.fp, "##\n");
    fprintf(tstotal.fp, "## COLUMN2.name  Total mass\n");
    fprintf(tstotal.fp,
	    "## COLUMN2.long_name  Total mass\n");
    fprintf(tstotal.fp, "## COLUMN2.units  kg\n");
    fprintf(tstotal.fp, "## COLUMN2.missing_value  0.000000\n");
    fprintf(tstotal.fp, "##\n");
    fprintf(tstotal.fp, "## COLUMN3.name  Total_volume\n");
    fprintf(tstotal.fp, "## COLUMN3.long_name  Total volume\n");
    fprintf(tstotal.fp, "## COLUMN3.units  m3\n");
    fprintf(tstotal.fp, "## COLUMN3.missing_value  0.000000\n");
    fprintf(tstotal.fp, "##\n");
    fprintf(tstotal.fp, "## COLUMN4.name  Total_heat\n");
    fprintf(tstotal.fp, "## COLUMN4.long_name  Total heat\n");
    fprintf(tstotal.fp, "## COLUMN4.units  deg C m3\n");
    fprintf(tstotal.fp, "## COLUMN4.missing_value  0.000000\n");
    fprintf(tstotal.fp, "##\n");
    fprintf(tstotal.fp, "## COLUMN5.name  Total_salt\n");
    fprintf(tstotal.fp, "## COLUMN5.long_name  Total salt\n");
    fprintf(tstotal.fp, "## COLUMN5.units    psu m3\n");
    fprintf(tstotal.fp, "## COLUMN5.missing_value  0.000000\n");
    fprintf(tstotal.fp, "##\n");
    n = 2;
    for (nn = 2; nn < master->ntot; nn++) { /* n=0:salt, n=1:temp    */
      if((trn = tracer_find_index(master->totname[nn], master->ntr,
				  master->trinfo_3d)) >= 0) {
	fprintf(tstotal.fp, "## COLUMN%d.name  Total_%s\n", n+4, 
		master->totname[n]);
	fprintf(tstotal.fp, "## COLUMN%d.long_name      Total %s\n", n+4,
		master->totname[n]);
	fprintf(tstotal.fp, "## COLUMN%d.units          (%s)(m3)\n", n+4, 
		master->trinfo_3d[trn].units);
	fprintf(tstotal.fp, "## COLUMN%d.missing_value  0.000000\n", n+4);
	fprintf(tstotal.fp, "##\n");
	n++;
      } else if((trn = tracer_find_index(master->totname[nn], master->ntrS,
				  master->trinfo_2d)) >= 0) {
	fprintf(tstotal.fp, "## COLUMN%d.name  Total_%s\n", n+4, 
		master->totname[n]);
	fprintf(tstotal.fp, "## COLUMN%d.long_name      Total %s\n", n+4,
		master->totname[n]);
	fprintf(tstotal.fp, "## COLUMN%d.units          %s\n", n+4, 
		master->trinfo_3d[trn].units);
	fprintf(tstotal.fp, "## COLUMN%d.missing_value  0.000000\n", n+4);
	fprintf(tstotal.fp, "##\n");
	n++;
      } else
	hd_quit("Can't recognize tracer %s for totals.\n", master->totname[n]);

    }
    n += 3;
    for (m = 0; m < geom->nobc; m++) {
      fprintf(tstotal.fp, "## COLUMN%1.1d.name  OBC[%1.1d]_volume_flux\n",
	      ++n, m);
      fprintf(tstotal.fp, "## COLUMN%1.1d.units m-3s-1\n",n);
      fprintf(tstotal.fp, "##\n");
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the tracers ids for all totals                              */
  if (master->ntot) {
    master->trtot = d_alloc_1d(master->ntot);
    for (n = 1; n <= master->nwindows; n++) {
      window_t *windat = window[n]->windat;
      windat->ntot = master->ntot;
      windat->totid = i_alloc_1d(master->ntot);
      windat->totid_sed = i_alloc_1d(master->ntot);
      windat->trtot = d_alloc_1d(master->ntot);
      for (m = 0; m < windat->ntot; m++) {
	if((trn = tracer_find_index(master->totname[m], master->ntr,
				    master->trinfo_3d)) >= 0)
	  windat->totid[m] = trn;
	if((trn = tracer_find_index(master->totname[m], master->ntrS,
				    master->trinfo_2d)) >= 0)
	  windat->totid[m] = -trn;
      }
      for (m = 0; m < windat->ntot; m++) {
	windat->totid_sed[m] = -1;
	if((trn = tracer_find_index(master->totname[m], master->nsed,
				    master->trinfo_sed)) >= 0)
	  windat->totid_sed[m] = trn;
      }
    }
  }
}

/* END init_totals()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates the total mass, heat and salt in a window              */
/*-------------------------------------------------------------------*/
void mass_diag(geometry_t *window,     /* Window geometry            */
	       window_t *windat,       /* Window data                */
	       win_priv_t *wincon      /* Window constants           */
  )
{
  int c, cc, cs, k, n, trn;
  int sn = windat->sno;
  int tn = windat->tno;
  int npor;

  /* Sum the water column mass                                       */
  windat->tvol = windat->tmass = windat->tarea = 0.0;
  memset(windat->trtot, 0, windat->ntot * sizeof(double));
  for (cc = 1; cc <= wincon->vc2; cc++) {
    c = wincon->s2[cc];
    cs = window->m2d[c];
    windat->tvol += window->cellarea[cs] * wincon->dz[c] * 
      wincon->Ds[cs];
    windat->tmass += windat->dens_0[c] * window->cellarea[cs] *
      wincon->dz[c] * wincon->Ds[cs];
    for (n = 0; n < windat->ntot; n++) {
      trn = windat->totid[n];
      if (trn >= 0)
	windat->trtot[n] += windat->tr_wc[trn][c] * window->cellarea[cs] *
	  wincon->dz[c] * wincon->Ds[cs];
    }
  }
  /* 2D tracers report the average.                                  */
  for (cc = 1; cc <= wincon->vcs2; cc++) {
    c = wincon->s2[cc];
    cs = window->m2d[c];
    windat->tarea += window->cellarea[cs];
    for (n = 0; n < windat->ntot; n++) {
      trn = windat->totid[n];
      if (trn < 0)
	windat->trtot[n] += (windat->tr_wcS[-trn][cs] * window->cellarea[cs]);
    }
  }
  for (n = 0; n < windat->ntot; n++) {
    trn = windat->totid[n];
    if (trn < 0)
      windat->trtot[n] /= ((double)window->nwindows * windat->tarea);
  }

  /* Add sediment mass if this tracer has a sediment component       */
  npor = tracer_find_index("porosity", wincon->ntr, wincon->trinfo_sed);
  for (n = 0; n < windat->ntot; n++) {
    if ((trn = windat->totid_sed[n]) >= 0) {
      for(cc = 1; cc <= window->b2_t; cc++) {
	c = window->w2_t[cc];
	for(k = 0; k <= window->sednz-1; k++) {
	  double porosity = (npor >= 0) ? windat->tr_sed[npor][k][c] : 0.4;
	  if (!wincon->trinfo_3d[trn].dissol) porosity = 1.0;
	  windat->trtot[n] += windat->tr_sed[trn][k][c] * 
	    window->cellarea[c]*
	    (window->gridz_sed[k+1][c] - window->gridz_sed[k][c]) * porosity;
	}
      }
    }
  }

  /* Compute boundary volume fluxes                                  */
  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    double *vel, v;
    double *hat;
    int sc = 1.0, *cv;
    int ee, e, es;

    vel = windat->u1;
    hat = window->h1au1;
    cv = open->obc_e1;

    windat->vf[n] = 0.0;
    for (ee = 1; ee <= open->no3_e1; ee++) {
      e = open->obc_e1[ee];
      es = window->m2de[e];
      c = open->obc_e2[ee];
      sc = -1.0 * (double)open->dir[ee];
      v = vel[e];
      windat->vf[n] += (sc * v * hat[es] * wincon->dz[c]);
    }
  }
}

/* END mass_diag()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Prints the total mass, heat and salt in the domain                */
/*-------------------------------------------------------------------*/
void print_total_mass(master_t *master) /* Model grid data structure */
{
  geometry_t *geom = master->geom;
  int n;

  if (master->t >= tstotal.tsout - DT_EPS) {
    fprintf(tstotal.fp, "%f %f %f ",
	    master->t/86400, master->tmass, master->tvol);
    for (n = 0; n < master->ntot; n++) {
      fprintf(tstotal.fp, "%f ", master->trtot[n]);
    }
    for (n = 0; n < geom->nobc; n++) {
      fprintf(tstotal.fp, "%f ", master->vf[n]);
      /*master->vf[n] = 0.0;*/
    }
    fprintf(tstotal.fp, "\n");
    fflush(tstotal.fp);
    tstotal.tsout += tstotal.tsdt;
  }
}

/* END print_total_mass()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set up the error norms output timeseries file.         */
/*-------------------------------------------------------------------*/
void init_error_norms(master_t *master)    /* Master data            */
{
  int n, m, nn, trn;

  /*---------------------------------------------------------------*/
  /* Initialize the output file                                    */
  if (master->errornorm & N_ERRN && tserror.fp == NULL) {
    if (strlen(master->opath))
      sprintf(tserror.pname, "%serrornorm.ts", master->opath);
    else
      sprintf(tserror.pname, "errornorm.ts");
    if ((tserror.fp = fopen(tserror.pname, "w")) == NULL)
      hd_quit("init_error_norms: Can't open file %s\n", tserror.pname);
    tserror.tsdt = master->enorm_dt;
    tserror.master = master;
    tserror.i = 1;
    tserror.j = 1;
    tserror.k = 1;
    tserror.tsout = master->t;
    fprintf(tserror.fp, "## COLUMNS %d\n", 4);
    fprintf(tserror.fp, "##\n");
    fprintf(tserror.fp, "## COLUMN1.name  Time\n");
    fprintf(tserror.fp, "## COLUMN1.long_name  Time\n");
    fprintf(tserror.fp,
	    "## COLUMN1.units  days since 1990-01-01 00:00:00 +10\n");
    fprintf(tserror.fp, "## COLUMN1.missing_value -999\n");
    fprintf(tserror.fp, "##\n");
    fprintf(tserror.fp, "## COLUMN2.name  cells\n");
    fprintf(tserror.fp, "## COLUMN2.long_name  Number of cells\n");
    fprintf(tserror.fp, "## COLUMN2.units  \n");
    fprintf(tserror.fp, "## COLUMN2.missing_value  0.000000\n");
    fprintf(tserror.fp, "##\n");
    fprintf(tserror.fp, "## COLUMN3.name  L1\n");
    fprintf(tserror.fp, "## COLUMN3.long_name  L1 error loss norm\n");
    fprintf(tserror.fp, "## COLUMN3.units  \n");
    fprintf(tserror.fp, "## COLUMN3.missing_value  0.000000\n");
    fprintf(tserror.fp, "##\n");
    fprintf(tserror.fp, "## COLUMN4.name  L2\n");
    fprintf(tserror.fp, "## COLUMN4.long_name  L2 error loss norm\n");
    fprintf(tserror.fp, "## COLUMN4.units  \n");
    fprintf(tserror.fp, "## COLUMN4.missing_value  0.000000\n");
    fprintf(tserror.fp, "##\n");
    fprintf(tserror.fp, "## COLUMN5.name  Linf\n");
    fprintf(tserror.fp, "## COLUMN5.long_name  L-infinity error loss norm\n");
    fprintf(tserror.fp, "## COLUMN5.units  \n");
    fprintf(tserror.fp, "## COLUMN5.missing_value  0.000000\n");
    fprintf(tserror.fp, "##\n");
  }
}

/* END init_error_norms()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Prints the error norms to file. These are relative error norms    */
/* computed as in Ringler et al (2010), Eq. 74 and 75.               */
/*-------------------------------------------------------------------*/
void print_error_norms(master_t *master) /* Global data              */
{
  geometry_t *geom = master->geom;
  double L1, L2, Linf;
  int n;

  if (master->t >= tserror.tsout - DT_EPS) {
    n = (master->errornorm & N_TRA2D) ? geom->b2_t : geom->b3_t;
    L1 = master->enorm[1] / master->enorm[2];
    L2 = sqrt(master->enorm[3] / master->enorm[0]) / 
      sqrt(master->enorm[4] / master->enorm[0]);
    Linf = master->enorm[5] / master->enorm[6];
    fprintf(tserror.fp, "%f %d %f %f %f\n",
	    master->t/86400, n, L1, L2, Linf);
    fflush(tserror.fp);
    tserror.tsout += tserror.tsdt;
    memset(master->enorm, 0, 9 * sizeof(double));
  }
}

/* END print_error_norms()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computs the error norm loss diagnostic on the windows             */
/*-------------------------------------------------------------------*/
void error_norms_w(geometry_t *window,   /* Window geometry          */
		   window_t *windat,     /* Window data              */
		   win_priv_t *wincon    /* Window constants         */
  )
{
  int c, cc, cs;
  double *tr1, *tr2;
  double vol, dif, *dz;
  int vc, *cells;

  if (wincon->errornorm & N_TRA2D) {
    vc = window->b2_t;
    cells = window->w2_t;
    dz = wincon->one;
    tr1 = windat->tr_wcS[wincon->normt1];
    tr2 = windat->tr_wcS[wincon->normt2];
  } else {
    vc = window->b3_t;
    cells = window->w3_t;
    dz = wincon->dz;
    if (wincon->errornorm & N_SURF) {
      vc = window->b2_t;
      cells = window->w2_t;
      dz = wincon->one;
    }
    tr1 = windat->tr_wc[wincon->normt1];
    tr2 = windat->tr_wc[wincon->normt2];
  }

  memset(wincon->enorm, 0, 9 * sizeof(double));
  for (cc = 1; cc <= vc; cc++) {
    c = cells[cc];
    cs = window->m2d[c];
    vol = window->cellarea[cs] * dz[c];
    dif = fabs(tr1[c] - tr2[c]);
    wincon->enorm[0] += vol;
    /* L1 norm                                                       */
    wincon->enorm[1] += dif * vol;
    wincon->enorm[2] += tr1[c] * vol;
    /* L2 norm                                                       */
    wincon->enorm[3] += dif * dif * vol;
    wincon->enorm[4] += tr1[c] * tr1[c] * vol;
    /* L infinity norm                                               */
    if (dif > wincon->enorm[5]) wincon->enorm[5] = dif;
    if (tr1[c] > wincon->enorm[6]) wincon->enorm[6] = tr1[c];
  }
}

/* END error_norms_w()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computs the error norm loss diagnostic on the master              */
/*-------------------------------------------------------------------*/
void error_norms_m(master_t *master,     /* Master data              */
		   geometry_t *window,   /* Window geometry          */
		   window_t *windat,     /* Window data              */
		   win_priv_t *wincon    /* Window constants         */
  )
{
  master->enorm[0] += wincon->enorm[0];
  master->enorm[1] += wincon->enorm[1];
  master->enorm[2] += wincon->enorm[2];
  master->enorm[3] += wincon->enorm[3];
  master->enorm[4] += wincon->enorm[4];
  if (wincon->enorm[5] > master->enorm[5]) master->enorm[5] = wincon->enorm[5];
  if (wincon->enorm[6] > master->enorm[6]) master->enorm[6] = wincon->enorm[6];
}


/* END error_norms_m()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Returns the maximum velocity threshold for constraining velocity  */
/*-------------------------------------------------------------------*/
double con_velmax(double vel, double velmax)
{
  return((vel > 0.0) ? velmax : -velmax);
}

/* END con_velmax()                                                  */
/*-------------------------------------------------------------------*/

void order (double *a, double *b)
{
  double val;
  if (*a > *b) {
    val = *a;
    *a = *b;
    *b = val;
  }
}

/*-------------------------------------------------------------------*/
/* Returns the median velocity in the vicinity of the cell for       */
/* constraining velocity.                                            */
/*-------------------------------------------------------------------*/
double con_median(geometry_t *window,     /* Window geometry         */
		  double *vel,            /* Velocity array          */
		  double velmax,          /* Velocity threshold      */
		  int c,                  /* Sparse coodinate        */
		  int mode                /* Filtering mode          */
)                  
{
  int i, j, m, nns, ns = 3;
  double *v, ret;
  int *st;
  int cl;

  /* Get the stencil                                                 */
  if (mode & ST_EDGE) {
    cl = window->e2ijk[c];
    j = window->e2e[c][0];
    st = stencil(window, cl, &ns, mode, j);
  } else {
    st = stencil(window, c, &ns, mode, 0);
  }

  m = ns / 2;
  v = d_alloc_1d(ns);
  nns = 0;
  for (i = 0; i < ns; i++) {
    if (!isnan(vel[st[i]])) {
      v[i] = vel[st[i]];
      nns++;
    }
  }
  m = nns / 2;

  for (i = 0; i < nns; i++)
    for (j = nns-1; i < j; --j)
      order(&v[j-1], &v[j]);

  ret = v[m];
  if (ret > velmax) ret = velmax;
  if (ret < -velmax) ret = -velmax;
  i_free_1d(st);
  d_free_1d(v);
  return(ret);
}

/* END con_median()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Returns the mean velocity in the vicinity of the cell for         */
/* constraining velocity.                                            */
/*-------------------------------------------------------------------*/
double con_mean(geometry_t *window,       /* Window geometry         */
		double *vel,              /* Velocity array          */
		double velmax,            /* Velocity threshold      */
		int c,                    /* Sparse coodinate        */
		int mode
		)
{
  int i, ns;
  double m, n, ret;
  double *v;
  int *st;

  /* Get the stencil                                                 */
  if (mode & ST_EDGE) {
    int cl = window->c2e[c][0];
    int j = window->e2e[c][0];
    st = stencil(window, cl, &ns, mode, j);
  } else {
    st = stencil(window, c, &ns, mode, 0);
  }

  v = d_alloc_1d(ns);
  for (i = 0; i < ns; i++) {
    v[i] = vel[st[i]];
  }

  m = n = 0.0;
  for (i = 0; i < ns; i++) {
    if (v[i] > -velmax && v[i] < velmax) {
      m += v[i];
      n += 1.0;
    }
  }
  ret = (n > 0.0) ? m / n : con_velmax(vel[c], velmax);
  d_free_1d(v);
  return(ret);
}

/* END con_mean()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to get the spatial percentiles of a snapshot of tracer    */
/* distribution.                                                     */
/*-------------------------------------------------------------------*/
void perc_diag(geometry_t *window, 
	       window_t *windat,   
	       win_priv_t *wincon 
	       )
{
  int cc, c, nc;
  int *ot = wincon->s4;
  double *p = wincon->w8;
  double *tr;
  int *vec, nv, sf;

  /* Set the cells to process                                        */
  if (wincon->percmsk[0]) {
    vec = wincon->s2;
    nv = wincon->vcs2;
    sf = 1;
  } else {
    vec = window->w3_t;
    nv = window->b3_t;
    sf = 0;
  }

  /* Initialize                                                      */
  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    windat->perc[c] = NaN;
  }
  nc = 0;
  tr = windat->tr_wc[wincon->trperc];

  /* Set the array to sort and sort                                  */
  for (cc = 1; cc <= nv; cc++) {
    c = vec[cc];
    if ((sf && wincon->percmsk[window->m2d[c]]) || (!sf && wincon->percmsk[c])) {
      nc++;
      p[nc] = tr[c];
      ot[nc] = c;
    }
  }
  quicksort(p, ot, 1, nc);

  /* Set the percentile vale at the grid location                    */
  windat->perc[ot[1]] = 0.0;
  windat->perc[ot[nc]] = 100.0;
  for (cc = 2; cc < nc; cc++) {
    c = (sf) ? window->m2d[ot[cc]] : ot[cc];
    windat->perc[c] = (double)cc * 100.0 / (double)(nc + 1);
  }
}

/* END perc_diag()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sorts an array into increasing order                              */
/*-------------------------------------------------------------------*/
void psorts(double *a, int *c, int n)
{
  int i, j;
  for (i = 1; i <= n; ++i)
    for (j = n; i < j; --j)
      porders(&a[j - 1], &a[j], &c[j - 1], &c[j]);
}

void porders(double *p, double *q, int *i, int *j)
{
  double t1;
  int t2;
  if (*p > *q) {
    t1 = *p;
    *p = *q;
    *q = t1;
    t2 = *i;
    *i = *j;
    *j = t2;
  }
}

/* END psorts()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sorts an array into increasing order using quicksort              */
/*-------------------------------------------------------------------*/
void quicksort(double *p, int *q, int l, int r)
{
  int j;

  if(l < r) {
    j = partition(p, q, l, r);
    quicksort(p, q, l, j-1);
    quicksort(p, q, j+1, r);
  }
}


int partition(double *p, int *q, int l, int r)
{
  int i, j, m;
  double pivot = p[l], t;
  i = l; j = r+1;
		
  while(1) {
     do ++i; while(p[i] <= pivot && i <= r );
     do --j; while(p[j] > pivot );
     if( i >= j ) break;
     t = p[i]; p[i] = p[j]; p[j] = t;
     m = q[i]; q[i] = q[j]; q[j] = m;
  }
  t = p[l]; p[l] = p[j]; p[j] = t;
  m = q[l]; q[l] = q[j]; q[j] = m;
  return j;
}

/* END quicksort()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routines to compute percentiles from a nominated output file      */
/*-------------------------------------------------------------------*/
static int compare_double(const void* v1, const void* v2)
{
double d1 = *(double *)v1;
double d2 = *(double *)v2;

    if (isnan(d1)) return 1;

    if (d1 > d2)
        return 1;
    else if (d1 < d2)
        return -1;
    return 0;
}

static void stats(double *data, int nn, double *results) {
  int i, n = 0;

  qsort(data, nn, sizeof(double), compare_double);

  /* Filter out NaN's */
  for (i=nn-1; i>=0; --i) {
    if (!isnan(data[i])) {
      n = i+1;
      break;
    }
  }

  for (i=0; i<Ns; ++i) {
    if (n != 0) {
      double r = percentiles[i] * (n - 1);
      int nr = (int) r;
      double fr = r - nr;
      results[i] = (nr >= n-1) ? data[n-1]
	: data[nr] * (1.0 - fr) + data[nr+1] * fr;
    } else {
      results[i] = NaN;
    }
  }
}

void calc_perc(FILE *fp) 
{
  FILE *dfp, *fopen();
  parameters_t *params;
  geometry_t *geom, *sgrid;
  dump_data_t *dumpdata;
  char pfile[MAXSTRLEN];
  char files[MAXNUMTSFILES][MAXSTRLEN];
  char *arg[MAXSTRLEN * MAXNUMARGS];
  char outfile[MAXSTRLEN];
  char infile[MAXSTRLEN];
  char vlist[MAXSTRLEN];
  char buf[MAXSTRLEN];
  char timeunit[MAXSTRLEN];
  char iunits[MAXSTRLEN];
  cstring *filenames;
  timeseries_t **tsfiles;
  datafile_t *df;
  size_t ndumps = 0;
  size_t start[4];
  size_t count[4];
  char **vars = NULL;
  int ntsfiles;
  int nvars = 0;
  int outfid;
  int i, j, k, kk, cc, c, cs, v, n, m, mm;
  int vids[Ns];
  int ndims;
  int *sdump;
  int *edump;
  int nargs;
  int spf;
  int first = 1;
  int mnc = 0;
  int *mncid;
  int r0, r1;
  int dsize = 0, asize;
  int varids[MAXNUMTSFILES];
  int vid;
  int ad;
  int inc;
  int *s2i, *s2j, *s2k, ***map;
  double *data, results[Ns];
  double stime, etime;
  double rfrac;
  double d1, d2, dtf, dt = 0.0;
  double mean, nm;
  double nn, ntot, clck, mclck;
  double elapsed = 0.0, start_time;
  long t;
  struct timeval tm1;
  int icdfid;
  dump_file_t *dfout;
  int oset, doset;
  int **k2c;
  int istart, kend;
  int nface2;

  params = params_alloc();
  set_default_param(params);

  /*----------------------------------------------------------------*/
  /* Build unstructured mappings based on the INPUT_FILE. The       */
  /* output file will be created with the same format.              */
  if (!prm_read_char(fp, "INPUT_FILE", params->idumpname)) return;
  icdfid = dump_open_us(params, params->idumpname, 1);
  meshstruct_us(params);
  sgrid = (geometry_t *)malloc(sizeof(geometry_t));
  memset(sgrid, 0, sizeof(geometry_t));
  build_sparse_grid_us(params, sgrid, params->nz, params->bathy,
		       params->layers);
  /*dumpdata_read_us(geom, params, master, dumpdata, icdfid, timeindex);*/
  dump_close(icdfid);
  geom = master->geom;
  dumpdata = master->dumpdata;
  oset = params->oset;
  doset = (dumpdata->start_index) ? 0 : 1;

  /*
  master = master_alloc();
  master->tsfile_caching = 0;
  master->lyear = 0;
  */
  master->means = params->means;
  start_time = time(NULL);
  time(&t);

  /*-----------------------------------------------------------------*/
  /* Build the reverse map for reading UGRID files                   */
  k2c = i_alloc_2d(geom->szcS, geom->nz);
  for (cc = 1; cc <= geom->b3_t; cc++) {
    c = geom->w3_t[cc];
    k = geom->s2k[c];
    cs = geom->m2d[c];
    k2c[k][cs-oset] = c;
  }

  /*-----------------------------------------------------------------*/
  /* Set the time unit                                               */
  prm_set_errfn(hd_silent_warn);
  prm_read_char(fp, "TIMEUNIT", params->timeunit);
  strcpy(master->timeunit, params->timeunit);
  if (!prm_read_char(fp, "OUTPUT_TIMEUNIT", params->output_tunit))
    strcpy(params->output_tunit, params->timeunit);
  if (prm_read_char(fp, "P_DT", buf))
    tm_scale_to_secs(buf, &dt);
  if (prm_read_char(fp, "PROJECTION", params->projection)) {
    strcpy(projection, params->projection);
    ts_set_default_proj_type(projection);
  }

  /*-----------------------------------------------------------------*/
  /* Read the input files                                            */
  prm_set_errfn(hd_warn);
  if (!prm_read_char(fp, "P_IFILE", pfile)) return;
  ntsfiles = parseline(pfile, (char **)files, MAXNUMTSFILES);
  filenames = (cstring *)malloc(sizeof(cstring) * ntsfiles);
  for (n = 0; n < ntsfiles; ++n)
    strcpy(filenames[n], ((char **)files)[n]);
  tsfiles = hd_ts_multifile_read(master, ntsfiles, filenames);
  if (tsfiles == NULL) {
    free(filenames);
    free(tsfiles);
    return;
  }
  hd_warn("percentiles: Input file %s read OK.\n", pfile);

  /*-----------------------------------------------------------------*/
  /* Read the output files                                           */
  if (!prm_read_char(fp, "P_OFILE", buf)) {
    free(filenames);
    free(tsfiles);
    hd_warn("percentiles: Output file name P_OFILE required.\n");
    return;
  }
  if (prm_read_char(fp, "OutputPath", pfile))
    sprintf(outfile,"%s%s", pfile, buf);
  else
    strcpy(outfile, buf);
  hd_warn("percentiles: Output file %s created OK.\n", outfile);

  /*-----------------------------------------------------------------*/
  /* Read the time limits for the percentile computations            */
  if (prm_read_char(fp, "P_STIME", buf)) {
    tm_scale_to_secs(buf, &stime);
    params->t = stime;
  } else {
    hd_warn("percentiles: Varible start time P_STIME required.\n");
    free(filenames);
    free(tsfiles);
    return;
  }
  if (prm_read_char(fp, "P_ETIME", buf)) {
    tm_scale_to_secs(buf, &etime);
  } else {
    hd_warn("percentiles: Varible start time P_ETIME required.\n");
    free(filenames);
    free(tsfiles);
    return;
  }

  /*-----------------------------------------------------------------*/
  /* Read the variables to get percentiles for                       */
  if (!prm_read_char(fp, "P_VARS", vlist)) {
    hd_warn("percentiles: Varible names P_VARS required.\n");
    free(filenames);
    free(tsfiles);
    return;
  }
  dumpdata->nvars = nvars = parseline(vlist, arg, MAXNUMARGS);
  dumpdata->vars = (char **)malloc(sizeof(char *)*nvars);
  j = 0;
  for (i=0; i<dumpdata->nvars; ++i) {
    if (hd_ts_multifile_get_index(ntsfiles, tsfiles,
				  filenames, arg[i], &n) <= 0)
      hd_warn("Can't find variable %s in files.\n", arg[i]);
    else {
      dumpdata->vars[j] = arg[i];
      j++;
    }
  }
  dumpdata->nvars = j;
  hd_warn("percentiles: variables %s read OK.\n", vlist);

  /*-----------------------------------------------------------------*/
  /* Get the dump numbers for each input file                        */
  spf = 0;
  params->ns2 = params->ns3 = params->nce1 = params->nce2 = 0;
  for (n = 0; n < ntsfiles; n++) {
    int nfiles;
    df_variable_t *var;
    df_multi_t *fd;
    spf = 0;
    df = tsfiles[n]->df;
    fd = (df_multi_t *)df->private_data;
    if (fd != NULL) { /* multi-netcdf-version */
      mnc = 1; 
      if (ntsfiles > 1) hd_quit("Can't process multiple 'multi-netcdf-version' files.\n");
    }
    if (n == 0) {
      nfiles = (mnc) ? fd->nfiles : ntsfiles;
      sdump = i_alloc_1d(nfiles);
      edump = i_alloc_1d(nfiles);
      memset(sdump, 0, nfiles * sizeof(int));
      memset(edump, 0, nfiles * sizeof(int));
    }
    /*---------------------------------------------------------------*/
    /* Read the input file dimensions                                */
    for (m = 0; m < df->nd; m++) {
      if (strcmp(df->dimensions[m].name, "ns2") == 0) {
	params->ns2 = df->dimensions[m].size;
	spf = 1;
      }
      if (strcmp(df->dimensions[m].name, "ns3") == 0) 
	params->ns3 = df->dimensions[m].size;
      if (strcmp(df->dimensions[m].name, "nMesh2_face") == 0) {
	params->ns2 = df->dimensions[m].size;
	spf = 2;
      }
      if (strcmp(df->dimensions[m].name, "Mesh2_layers") == 0) {
	params->nz = df->dimensions[m].size;
	spf = 2;
      }
      if (strcmp(df->dimensions[m].name, "i_centre") == 0) 
	params->nce1 = df->dimensions[m].size;
      if (strcmp(df->dimensions[m].name, "j_centre") == 0) 
	params->nce2 = df->dimensions[m].size;
      if (strcmp(df->dimensions[m].name, "k_centre") == 0) 
	params->nz = df->dimensions[m].size;
    }
    if (!spf) hd_quit("Input files are Cartesian : use SHOC for structured netCDF.\n");
    if (spf == 1) {
      hd_warn("Input files are sparse with size %d\n", params->ns3);
    }
    if (spf == 2) {
      hd_warn("Input and output files are UGRID with size %d x %d\n", params->ns2, params->nz);
    }

    /*---------------------------------------------------------------*/
    /* Get the dump numbers for each file                            */
    if (mnc) { /* Dump numbers for multi-netcdf-version files        */
      mncid = i_alloc_1d(fd->nfiles);
      for (m=0; m<fd->nfiles; ++m) {
	df_multi_file_t *f = &fd->files[m];
	/* Open the netcdf file                                      */
	if (nc_open(f->filename, NC_NOWRITE, &mncid[m]) == NC_NOERR) {
	  /* Copy the records to the df structure (adjusting for the */
	  /* timeunit and reset df->nrecords = f->nrecords.          */
	  memset(buf, 0, MAXSTRLEN);
	  nc_get_att_text(mncid[m], ncw_var_id(mncid[m], "t"), "units", buf);
	  ndumps = df->nrecords = f->nrecords;
	  for (j = 0; j < f->nrecords; j++) {
	    df->records[j] = f->records[j];
	    tm_change_time_units(buf, params->timeunit, &df->records[j], 1);
	  }
	  /* Find the record closest to t                            */
	  df_find_record(df, stime, &sdump[m], &r1, &rfrac);
	  df_find_record(df, etime, &edump[m], &r1, &rfrac);
	  if (edump[m] >= ndumps)
	    edump[m] = ndumps-1;
	  if (sdump[m] < 0)
	    sdump[m] = 0;  
	  dsize += edump[m]-sdump[m]+1;
	  /* Get the computation interval */
	  d1 = d2 = df->records[0];
	  if (edump[m] > 1) d2 = df->records[1];
	  dtf = d2 - d1;
	  inc = (dtf && dt > dtf) ? (int)(dt / dtf) : 1;
	} else
	  hd_quit("Can't open 'multi-netcdf file%d %s\n",m, f->filename);
      }
      df->ncid = mncid[0];
    } else {  /* Dump numbers for regular files                      */
      df_find_record(df, stime, &sdump[n], &r1, &rfrac);
      df_find_record(df, etime, &edump[n], &r1, &rfrac);
      ndumps = df->nrecords;
      if (edump[n] >= ndumps)
				edump[n] = ndumps-1;
      if (sdump[n] < 0)
				sdump[n] = 0;  
      dsize += edump[n]-sdump[n]+1;

      /* Get the computation interval */
      d1 = d2 = df->records[0];
      if (edump[n] > 1) d2 = df->records[1];
      dtf = d2 - d1;
      inc = (dtf && dt > dtf) ? (int)(dt / dtf) : 1;
    }
    if (inc > 1)
      hd_warn("File dt=%.0f, computation dt=%.0f : every %d data used.\n",
	      dtf, dt, inc);

    /*---------------------------------------------------------------*/
    /* Create the output file when the first input file is read      */
    if (first) {
      dfout = create_output(master, &outfid, outfile);
      first = 0;
    }
  }
  data = d_alloc_1d(dsize);

  /*-----------------------------------------------------------------*/
  /* Get the sparse maps and set up botz for sparse input files      */
  /*
  prm_read_int(fp, "LAYERFACES", &k);
  if (k-1 != params->nz) {
    if (spf == 1)
      hd_quit("Sparse input files must have 'k_centre = %d'\n",k-1);
    else
      hd_quit("UGRID input files must have 'Mesh2_layers = %d'\n",k-1);
  }
  */
  hd_warn("percentiles: Processing data.\n");
 


  /*-----------------------------------------------------------------*/
  /* Process the data                                                */ 
  nn = mclck = 0.0;
  ntot = (double)(params->nz * params->ns2 * dumpdata->nvars);
  for (v = 0; v < dumpdata->nvars; ++v) {
    /* Build variable percentile netCDF file id's */
    for (n = 0; n < Ns; ++n) {
      sprintf(buf, "%s_%3.3d", dumpdata->vars[v], (int)(percentiles[n]*100));
      vids[n] = ncw_var_id(outfid, buf);
    }
    mean = nm = 0.0;

    /* Get the number of dimensions for the variable. For sparse     */
    /* files variables are assumed to be 3D; for UGRID this may be   */
    /* 2D or 3D. Set layer end index and start and end indices of    */
    /* the input vector.                                             */
    oset = 0;
    if (mnc) {
      vid = ncw_var_id(mncid[0], dumpdata->vars[v]);
      nc_inq_varndims(mncid[0], vid, &ndims);
      if (spf == 2) nc_get_att_int(mncid[0], NC_GLOBAL, "start_index", &oset);
    } else {
      vid = ncw_var_id(df->ncid, dumpdata->vars[v]);
      nc_inq_varndims(df->ncid, vid, &ndims);
      if (spf == 2) nc_get_att_int(df->ncid, NC_GLOBAL, "start_index", &oset);
    }
    oset = (oset) ? 0 : 1;
    kend = 1;
    istart = 1;
    nface2 = params->ns3 + 1;
    if (spf == 2) {
      kend = (ndims == 3) ? params->nz : 1;
      istart = (oset) ? 0 : 1;
      /*nface2 = dumpdata->nface2 - dumpdata->start_index;*/
      nface2 = params->ns2;
    }

    /* Loop through the input file and extract time range for each   */
    /* cell.                                                         */
    for (k = 0; k < kend; ++k) {
      hd_warn("Percentiles: processing %s layer %d\n", dumpdata->vars[v], k);
      /* Unstructured input only (sparse or UGRID)                   */
      for (c = istart; c < nface2; c++) {
	kk = geom->nz - 1 - k;
	if (spf == 2 && !k2c[kk][c]) continue;
	kk = params->nz - 1 - k;
	nn += 1.0;
	gettimeofday(&tm1, NULL);
	clck = tm1.tv_sec + tm1.tv_usec * 1e-6;
	ad = 0;
	for (m = 0; m < ntsfiles; m++) {
	  df = tsfiles[m]->df;
	  if (mnc) { /* multi-netcdf-version                         */
	    df_multi_t *fd = (df_multi_t *)df->private_data;
	    for (mm=0; mm<fd->nfiles; ++mm) {
	      vid = ncw_var_id(mncid[mm], dumpdata->vars[v]);
	      start[0] = sdump[mm];
	      start[1] = c;
	      start[2] = 0;
	      count[2] = 0;
	      if (spf == 2 && ndims == 3) {
		start[1] = kk;
		start[2] = c;
		count[2] = 1;
	      }
	      count[0] = edump[mm]-sdump[mm]+1;
	      count[1] = 1;
	      nc_get_vara_double(mncid[mm], vid, start, count, &data[ad]);
	      if (!isnan(data[ad]) && data[ad] > 0.0) {
		mean += data[ad];
		nm += 1.0;
	      }
	      ad += count[0];
	    }
	  } else {
	    vid = ncw_var_id(df->ncid, dumpdata->vars[v]);
	    start[0] = sdump[m];
	    start[1] = c;
	    start[2] = 0;
	    count[2] = 0;
	    if (spf == 2 && ndims == 3) {
	      start[1] = kk;
	      start[2] = c;
	      count[2] = 1;
	    }
	    count[0] = edump[m]-sdump[m]+1;
	    count[1] = 1;
	    nc_get_vara_double(df->ncid, vid, start, count, &data[ad]);
	    j = edump[m]-sdump[m];
	    if (!isnan(data[j]) && data[j] > 0.0) {
	      mean += data[j];
	      nm += 1.0;
	    }
	    ad += count[0];
	  }

	  /* Resort the data into inc steps                          */
	  asize = 0;
	  for (m = 0; m < dsize; m++) {
	    if (m%inc == 0) {
	      data[asize] = data[m];
	      asize++;
	    }
	  }
	  /* Compute statistics                                      */
	  stats(data, asize, results);
	  /* Write data                                              */
	  if (spf == 1) {
	    start[0] = geom->s2k[c];
	    start[1] = geom->m2d[c];
	  } else {
	    if (ndims == 3) {
	      start[0] = kk;
	      start[1] = c;
	    } else {
	      start[0] = c;
	      start[1] = 0;
	    }
	  }
	  start[2] = 0;
	  count[0] = 1;
	  count[1] = 1;
	  count[2] = 1;
	  for (n = 0; n < Ns; ++n) {
	    nc_put_vara_double(outfid, vids[n],
			       start, count, &results[n]);
	  }	  
	  gettimeofday(&tm1, NULL);
	  clck = tm1.tv_sec + tm1.tv_usec * 1e-6 - clck;
	  mclck += clck;

	  dfp = fopen(diag_logfile, "w");
	  if (dfp != NULL) {
	    fprintf(dfp, "\n\nStart time :  %s\n", ctime(&t));
	    fprintf(dfp, "CPU time / iteration = %4.3f (%4.3f)\n", clck, mclck / nn);
	    elapsed = difftime(time(NULL), start_time) / 3600.0;
	    fprintf(dfp, "Elapsed time = %d day(s) %2.2d:%2.2d:%2.2d\n",
		    (int)(elapsed/24.0), (int)fmod(elapsed, 24.0), (int)fmod(elapsed*60, 60.0),
		    (int)fmod(elapsed*3600, 60.0));
	    elapsed = ((ntot - nn) * (mclck / nn)) / 3600.0;
	    fprintf(dfp, "Time remaining = %d day(s) %2.2d:%2.2d:%2.2d\n",
		    (int)(elapsed/24.0), (int)fmod(elapsed, 24.0), (int)fmod(elapsed*60, 60.0),
		    (int)fmod(elapsed*3600, 60.0));
	    fprintf(dfp, "(%d %d) of (%d %d) : %s(%d/%d) : %% complete = %5.2f\n", c, k,
	      params->ns2, params->nz, dumpdata->vars[v], v+1, nvars, 100.0 * nn / ntot);
	    fclose(dfp);
	  }
	}
      }
      hd_warn("  Domain mean of %s = %f\n", dumpdata->vars[v], mean / nm);
    }
  }

  df_ugrid_close(dumpdata, dfout);
  free((dump_file_t *) dfout);

  if (mnc) {
    for (m = 0; m < ntsfiles; m++) {
      df = tsfiles[m]->df;
      df_multi_t *fd = (df_multi_t *)df->private_data;
      for (mm=0; mm<fd->nfiles; ++mm) {
				nc_close(mncid[mm]);
      }
    }
    i_free_1d(mncid);
  }
  i_free_2d(k2c);
  d_free_1d(data);
  free(filenames);
  free(tsfiles);
  free(params);
  master_free(master);
  exit(0);
}

dump_file_t *create_output(master_t *master, int *fid, char *filename)
{
  geometry_t *geom = master->geom;
  parameters_t *params = master->params;
  char buf[MAXSTRLEN];
  int dims[4];
  size_t start[4];
  size_t count[4];
  dump_file_t *df = NULL;
  df_ugrid_data_t *data;

  df = (dump_file_t *)malloc(sizeof(dump_file_t));
  memset(df, 0, sizeof(dump_file_t));
  strcpy(df->name, filename);
  df->tout = 0.0;
  df->tinc = 1.0;
  df->tstop = df->tout + df->tinc;
  df->bpv = 8;

  df->nvars = 0;
  df->landfill = locate_landfill_function("default");
  df->ilower = 0;
  df->jlower = 0;
  /* Output files are always 3D, even though only one layer may be   */
  /* processed.                                                      */
  df->klower = 0;
  df->nce1 = dumpdata->nce1;
  df->nfe1 = dumpdata->nfe1;
  df->nce2 = dumpdata->nce2;
  df->nfe2 = dumpdata->nfe2;
  df->nz = params->nz;
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
  strcpy(df->tunit, dumpdata->output_tunit);
  strcpy(df->name, filename);

  df_ugrid_create(dumpdata, df);
  data = (df_ugrid_data_t *)df->private_data;
  /*
  dfout = df;
  printf("%x\n",dfout);
  return data->fid;
  */
  *fid = data->fid;
  return(df);
}

/* END percentiles()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to check that the vertical integral of 3D fluxes is equal */
/* to the 2D fluxes. Differences should be around 1e-8. Note: if     */
/* this routine is called after set_new_cells_ then the 1st and 3rd  */
/* terms will be large since dzu1 or dzu2 are reset using the        */
/* elevation at the forward timestep.                                */
/* The input coordinate, dc, must be the 2D (top) coordinate.        */
/*-------------------------------------------------------------------*/
void check_flux(geometry_t *window, /* Window geometry               */
                window_t *windat,   /* Window data                   */
                win_priv_t *wincon, /* Window constants              */
                int dw, int dc)
{
  int e, es;
  int j;
  int cs = window->m2d[dc];
  double d1, d2;

  if (window->wn == dw) {
    d1 = d2 = 0.0;
    for (j = 1; j <= window->npe[cs]; j++) {
      e = window->c2e[j][dc];
      es = window->m2de[e];
      d1 = windat->u1[e] * windat->dzu1[e] * wincon->mdx[es] * windat->dt *
	    window->h1au1[es];
      d2 = windat->u1flux3d[e];
      emstag(LDEBUG,"hd:vel3d:check_flux", "check_flux edge%d %e %e\n", j,
	     windat->u1flux[es] - d1,
	     windat->u1flux[es] - d2 * windat->dt);
      printf("check_flux cell%d edge%d %e %e\n", dc, j,
	     windat->u1flux[es] - d1,
	     windat->u1flux[es] - d2 * windat->dt);
    }
  }
}

/* END check_flux()                                                  */
/*-------------------------------------------------------------------*/

#define DEG2RAD(d) ((d)*M_PI/180.0)
#define RAD2DEG(r) ((r)*180.0/M_PI)

#define RADIUS 6370997.0
#define ECC 0.0
#define GEODESIC(x1, y1, x2, y2) (geod_inv_geod_fwd_sodanos(DEG2RAD(x1), DEG2RAD(y1),\
                                 DEG2RAD(x2), DEG2RAD(y2),\
                                 RADIUS, ECC))

/*-------------------------------------------------------------------*/
/* Compares bathymetry in a database with that on the grid           */
/*-------------------------------------------------------------------*/
void bathy_compare(master_t *master)
{
  FILE *op;
  geometry_t *geom = master->geom;
  int ncerr;
  int fid;
  size_t nce1, nce2;
  size_t start[4];
  size_t count[4];
  double **bathy;
  double *lat, *lon;
  double **gridx, **gridy;
  char bathyvar[MAXSTRLEN];
  char bathyname[MAXSTRLEN];
  double val, d1;
  double gm, g1, g2;
  int n, i, j, c, cc, cn, e;
  int mi, mj, ni, nj;
  int is, js, ie, je;
  int nz = geom->nz - 1;
  double latmx, latmn, lonmx, lonmn;
  int iinc = 1;
  int jinc = 1;
  char *fields[MAXSTRLEN * MAXNUMARGS];
  int dograd = 1;

  if (!strlen(master->bathystats)) return;

  strcpy(bathyvar, "height");
  i = parseline(master->bathystats, fields, MAXNUMARGS);
  strcpy(bathyname, fields[0]);
  if (i > 1) {
    iinc = atoi(fields[1]);
    jinc = atoi(fields[2]);
  }
  if (i > 3)
    dograd = is_true(fields[3]);

  /*-----------------------------------------------------------------*/
  /* Read the dimensions of the bathymetry grid                      */
  /* Open the dump file for reading                                  */
  if ((ncerr = nc_open(bathyname, NC_NOWRITE, &fid)) != NC_NOERR) {
    hd_warn("Can't find transport source grid file %s\n", master->bathystats);
    hd_quit((char *)nc_strerror(ncerr));
  }

  /* Get dimensions                                                  */
  nc_inq_dimlen(fid, ncw_dim_id(fid, "lat"), &nce2);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "lon"), &nce1);

  /* Read the vertical layer structure                               */
  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;
  count[0] = nce1;
  count[1] = 0;
  count[2] = 0;
  count[3] = 0;

  lon = d_alloc_1d(nce1);
  nc_get_vara_double(fid, ncw_var_id(fid, "lon"), start, count, lon);
  count[0] = nce2;
  lat = d_alloc_1d(nce2);
  nc_get_vara_double(fid, ncw_var_id(fid, "lat"), start, count, lat);
  count[0] = nce2;
  count[1] = nce1;
  bathy =  d_alloc_2d(nce1, nce2);
  nc_get_vara_double(fid, ncw_var_id(fid, "height"), start, count, bathy[0]);

  /* Get the model grid bounds                                       */
  latmx = -90.0;
  latmn = 90.0;
  lonmn = 360.0;
  lonmx = -360.0;
  for (cc = 1; cc < geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    if (geom->cellx[c] > lonmx) lonmx = geom->cellx[c];
    if (geom->cellx[c] < lonmn) lonmn = geom->cellx[c];
    if (geom->celly[c] > latmx) latmx = geom->celly[c];
    if (geom->celly[c] < latmn) latmn = geom->celly[c];
  }
  /* Truncate the bathy grid                                         */
  is=js=0;
  ie=nce1;
  je=nce2;
  for (i = 0; i < nce1; i++) {
    if (lon[i] > lonmn) {
      is = i;
      break;
    }
  }
  for (i = nce1-1; i >= 0; i--) {
    if (lon[i] < lonmx) {
      ie = i;
      break;
    }
  }
  for (j = 0; j < nce2; j++) {
    if (lat[j] > latmn) {
      js = j;
      break;
    }
  }
  for (j = nce2-1; j >= 0; j--) {
    if (lat[j] < latmx) {
      je = j;
      break;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the range of depth differences                              */
  memset(master->bathy_range_min, 0.0, geom->szcS);
  memset(master->bathy_range_max, 0.0, geom->szcS);
  for (j = js; j < je; j += jinc) {
    for (i = is; i < ie; i += iinc) {
      if (!(c = hd_grid_xytoij(master, lon[i], lat[j], &mi, &mj))) continue; 
      if (c <= 0 || c > geom->ewetS) continue;
      d1 = fabs(bathy[j][i]) - fabs(geom->botz[c]);
      /* Get the difference between model bathy and the deepest      */
      /* source bathy.                                               */
      if (d1 < 0 && fabs(d1) > master->bathy_range_min[c]) master->bathy_range_min[c] = fabs(d1);
      /* Get the difference between model bathy and the deepest      */
      /* source bathy.                                               */
      if (d1 > 0 && d1 > master->bathy_range_max[c]) master->bathy_range_max[c] = d1;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the range of depth gradients                                */
  dograd = 0;
  if (dograd) {
    double maxgrad, mingrad;
    for (c = 0; c < geom->sgnumS; c++) {
      master->bathy_grad_max[c] = 0.0;
      master->bathy_grad_min[c] = HUGE;
    }
    for (j = js+1; j < je-1; j += jinc) {
      for (i = is+1; i < ie-1; i += iinc) {
	if (!(c = hd_grid_xytoij(master, lon[i], lat[j], &mi, &mj))) continue; 
	/*if (mi < 1 || mi > geom->nce1-1 || mj < 1 || mj > geom->nce2-1) continue;*/
	maxgrad = 0.0;
	mingrad = HUGE;

	/* Get the model grid gradients                                */
	if (c <= 0 || c > geom->ewetS) continue;
	g1 = 0.0;
	for (j = 1; j <= geom->npe[c]; j++) {
	  cn = geom->c2c[j][c];
	  e = geom->c2e[j][c];
	  g1 += fabs((geom->botz[c] - geom->botz[cn]) / geom->h2au1[e]);
	}
	g1 /= (double)geom->npe[c];
	
	/* Get the bathy gradient across i grid faces                  */
	if (!(c = hd_grid_xytoij(master, lon[i-1], lat[j], &ni, &nj))) continue; 
	if (mi != ni) {
	  val = (bathy[j][i] - bathy[j][i-1]) / GEODESIC(lon[i], lat[j], lon[i - 1], lat[j]);
	  maxgrad = max(maxgrad, val);
	  mingrad = min(mingrad, val);
	}
	if (mj != nj) {
	  val = (bathy[j][i] - bathy[j-1][i]) / GEODESIC(lon[i], lat[j], lon[i], lat[j - 1]);
	  maxgrad = max(maxgrad, val);
	  mingrad = min(mingrad, val);
	}
	/* Get the bathy gradient across j grid faces                  */
	if (!(c = hd_grid_xytoij(master, lon[i], lat[j-1], &ni, &nj))) continue; 
	if (mi != ni) {
	  val = (bathy[j][i] - bathy[j][i-1]) / GEODESIC(lon[i], lat[j], lon[i - 1], lat[j]);
	  maxgrad = max(maxgrad, val);
	  mingrad = min(mingrad, val);
	}
	if (mj != nj) {
	  val = (bathy[j][i] - bathy[j-1][i]) / GEODESIC(lon[i], lat[j], lon[i], lat[j - 1]);
	  maxgrad = max(maxgrad, val);
	  mingrad = min(mingrad, val);
	}
	master->bathy_grad_max[c] = max(fabs(fabs(maxgrad) - g1),
					master->bathy_grad_max[c]);
	master->bathy_grad_min[c] = min(fabs(fabs(mingrad) - g1),
					master->bathy_grad_min[c]);
      }
    }
    for (c = 0; c < geom->sgnumS; c++) {
      if (master->bathy_grad_min[c] == HUGE) master->bathy_grad_min[c] = 0.0;
    }
  }
  dograd = 1;
  if (dograd) {
    memset(master->bathy_grad_min, 0.0, geom->szcS);
    memset(master->bathy_grad_max, 0.0, geom->szcS);
    /* Get the model grid gradient                                   */
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = geom->w2_t[cc];
      g1 = 0.0;
      for (j = 1; j <= geom->npe[c]; j++) {
	cn = geom->c2c[j][c];
	e = geom->c2e[j][c];
	g1 += fabs((geom->botz[c] - geom->botz[cn]) / geom->h2au1[e]);
      } 
      master->bathy_grad_min[c] = g1 / (double)geom->npe[c];
    }
    /* Get the bathymetry database gradient                          */
    for (j = js+1; j < je-1; j += jinc) {
      for (i = is+1; i < ie-1; i += iinc) {
	if (!(c = hd_grid_xytoij(master, lon[i], lat[j], &mi, &mj))) continue; 
	/*if (mi < 1 || mi > geom->nce1-1 || mj < 1 || mj > geom->nce2-1) continue;*/
	if (c <= 0 || c > geom->ewetS) continue;
	g1 = (bathy[j][i] - bathy[j][i-1]) / GEODESIC(lon[i], lat[j], lon[i - 1], lat[j]);
	g2 = (bathy[j][i] - bathy[j-1][i]) / GEODESIC(lon[i], lat[j], lon[i], lat[j - 1]);
	val = sqrt(g1 * g1 + g2 * g2);
	master->bathy_grad_max[c] = max(master->bathy_grad_max[c], val);
      }
    }
    /* Get the absolute value of gradient difference                 */
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = geom->w2_t[cc];
      master->bathy_grad_max[c] = fabs(master->bathy_grad_max[c] - master->bathy_grad_min[c]);
    }
  }

  /* Create a list of bathymetry differences for timeseries points   */
  if ((op = fopen("ts_bathy_stats.txt", "w")) != NULL) {
    fprintf(op, "       n                name                       obs_depth mod_depth obs_grad mod_grad\n");
    for (n = 0; n < nts; ++n) {
      master_t *master = tslist[n].master;
      geometry_t *geom = master->geom;
      double dx, dy, dist, mindist = HUGE;

      for (mj = 0; mj < nce2; mj++) {
	for (mi = 0; mi < nce1; mi++) {
	  dx = tslist[n].x - lon[mi];
	  dy = tslist[n].y - lat[mj];
	  dist = dx * dx + dy * dy;
	  if (bathy[mj][mi] < 0.0 && dist < mindist) {
	    i = mi;
	    j = mj;
	    mindist = dist;
	  }
	}
      }
      g1 = (bathy[j][i+1] - bathy[j][i-1]) / (2.0 * GEODESIC(lon[i + 1], lat[j], lon[i - 1], lat[j]));
      g2 = (bathy[j+1][i] - bathy[j-1][i]) / (2.0 * GEODESIC(lon[i], lat[j + 1], lon[i], lat[j - 1]));
      val = sqrt(g1 * g1 + g2 * g2);
      if (!(c = hd_grid_xytoij(master, tslist[n].x, tslist[n].y, &mi, &mj))) continue; 
      d1 = fabs(bathy[j][i]) - fabs(geom->botz[c]);
      g1 = 0.0;
      for (j = 1; j <= geom->npe[c]; j++) {
	cn = geom->c2c[j][c];
	e = geom->c2e[j][c];
	g1 += fabs((geom->botz[c] - geom->botz[cn]) / geom->h2au1[e]);
      } 
      gm = g1 / (double)geom->npe[c];
      fprintf(op, "%8d %40s %8.4f %8.4f %8.4f %8.4f\n",n, tslist[n].pname, bathy[j][i], geom->botz[c], val, gm);
    }
    fclose(op);
  }

  d_free_1d(lat);
  d_free_1d(lon);
  d_free_2d(bathy);
  nc_close(fid);
}

/* END bathy_compare()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Produces a map of a spatial decorrelation distance.               */
/*-------------------------------------------------------------------*/
void calc_decorr(geometry_t *window,   /* Processing window          */
		 double *a,
		 double *dex,
		 int sz,
		 int mode,
		 double scale
		 )
{
  int c, cc, c2, j;
  int *vec, nvec;
  
  nvec = (mode) == 2 ? window->b2_t : window->b3_t;
  vec = (mode) == 2 ? window->w2_t : window->w3_t;

  /* Set a no gradient over lateral ghost cells                      */
  c = (mode == 2) ? geom->nbptS : geom->nbpt;
  for (cc = 1; cc <= c; cc++) {
    a[window->bpt[cc]] = a[window->bin[cc]];
  }

  /* Compute mean decorrelation length scales in all directions      */
  for (cc = 1; cc < nvec; cc++) {
    c = vec[cc];
    c2 = window->m2d[c];
    dex[c] = 0.0;
    for (j = 1; j <= window->npe[c2]; j++) {
      dex[c] += decorr(window, a, sz, c, j, scale);
    }
    dex[c] /= (double)window->npe[c2];
  }
}

/* END calc_decorr()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to find the spatial decorrelation distance in e1 or e2    */
/* directions. Based on Romanou at al, (2006) Journal of Climate, 19 */
/* 3378-3393.                                                        */
/*-------------------------------------------------------------------*/
double decorr(geometry_t *window, /* Processing window               */
	      double *a,          /* Variable to decorrelate         */
	      int sz,             /* Dist within which to compute    */
	      int c,              /* Starting index for computation  */
	      int j,              /* Direction to compute            */
	      double scale        /* Scaling factor                  */
	      )
{
  int e = window->c2e[j][c];
  int N = ceil(((double)sz * scale) / window->h2au1[window->m2de[e]]);
  int n, m, cp, *map;
  double d1;
  double vm = a[c];
  double p[N+1], v[N+1], d[N+2];
  double eps = 1e-10;

  map = window->c2c[j];

  /* Get the mean                                                    */
  cp = c;
  v[1] = a[c];
  vm = a[c];
  for (n = 2; n <= N; n++) {
    cp = map[cp];
    vm += a[cp];
    v[n] = a[cp];
  }
  vm /= (double)N;

  /* Get the autocorrelated lags                                     */
  cp = c;
  for (m = 1; m <= N; m++) {
    e = window->m2de[window->c2e[j][cp]];
    d[m] = (m == 1) ? window->h2au1[e] : d[m-1] + window->h2au1[cp];
    p[m] = 0.0;
    for (n = 1; n <= N - m - 1; n++) {
      p[m] += (v[n] - vm) * (v[n+m] - vm);
    }
    p[m] /= (double)(N-m);
    cp = map[cp];
  }

  /* Find the zero crossing point                                    */
  d1 = p[1];
  for (m = 1; m <= N; m++) {
    if (fabs(p[m])< eps) p[m] = 0.0;
    p[m] = (d1) ? p[m]/d1 : 0.0;
  }
  for (m = 1; m < N; m++) {
    if ((p[m] > 0.0 && p[m+1] < 0.0) ||	(p[m] < 0.0 && p[m+1] > 0.0))
      return(d[m] / scale);
  }
  return(d[N] / scale);
}

/* END decorr()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computes monotinicity violations and saves as % of runtime        */
/*-------------------------------------------------------------------*/
void calc_monotonic(geometry_t *window,
		    window_t *windat,   
		    win_priv_t *wincon 
		    )
{
  int c, cc;
  int tn = wincon->monon;
  double d1;
  int don = 0;

  if (wincon->trasc & FFSL) {
    return;
  }

  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    if (windat->tr_wc[tn][c] < wincon->monomn ||
	windat->tr_wc[tn][c] > wincon->monomx) {
      d1 = 0.0;

      if (wincon->lw) {
	qweights *lw = &wincon->lw[c];
	double mx = lw->tmx;
	double mn = lw->tmn;
	if (windat->tr_wc[tn][c] < mn)
	  d1 = fabs(windat->tr_wc[tn][c] - mn);
	if (windat->tr_wc[tn][c] > mx)
	  d1 = fabs(windat->tr_wc[tn][c] - mx);
      } else {

	/*
	  d1 = windat->mono[c] * 100.0 * (double)windat->nstep;
	  windat->mono[c]  = 100.0 * (d1 + 1.0) / (double)(windat->nstep + 1);
	*/
	/* Count the violations                                        */
	/*windat->mono[c] += 1.0;*/
	/* Save maximum difference                                     */
      if (windat->tr_wc[tn][c] < wincon->monomn)
	d1 = fabs(windat->tr_wc[tn][c] - wincon->monomn);
      if (windat->tr_wc[tn][c] > wincon->monomx)
	d1 = fabs(windat->tr_wc[tn][c] - wincon->monomx);
      /*
      if (windat->tr_wc[tn][c] > wincon->monomx) d1 = windat->tr_wc[tn][c];
      if (windat->tr_wc[tn][c] < wincon->monomn) d1 = windat->tr_wc[tn][c];
      windat->mono[c] = d1;
      */
      }
      windat->mono[c] = max(windat->mono[c], d1);
    }
  }
}

/* END calc_monotonic()                                              */
/*-------------------------------------------------------------------*/
