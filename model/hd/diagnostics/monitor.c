/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/diagnostics/monitor.c
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
 *  $Id: monitor.c 5841 2018-06-28 06:51:55Z riz008 $
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

double percentiles[] = {0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40,
			0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85,
			0.90, 0.95, 1.00};

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

#define Ns (sizeof(percentiles) / sizeof(double))

void wbdrycustom(parameters_t *params, FILE * fp, int bnum,
                 bdry_details_t *data);
void reset_cfl(master_t *master);
void reset_cfl_event(master_t *master, double mcfl, int mode);
double get_max_div_2d(geometry_t *window, window_t *windat,
		      win_priv_t *wincon, int *cl);
double get_max_div_3d(geometry_t *window, window_t *windat,
		      win_priv_t *wincon, int *cl);
double get_div_2d(geometry_t *window, window_t *windat,
		  win_priv_t *wincon, int c);
double get_ddt(geometry_t *window, double *T, double *TB, int *cl);
double buoyancy_frequency2(geometry_t *window, window_t *windat,
			   win_priv_t *wincon, int c);
double calc_min_cfl_3d(geometry_t *window, window_t *windat,
		       win_priv_t *wincon, int *cl);
void remove_tide(geometry_t *window, double *eta);
void print_total_mass(master_t *master);
void order (double *a, double *b);
void psorts(double *a, int *c, int n);
void porders(double *p, double *q, int *i, int *j);
void sound_channel(geometry_t *window, window_t *windat, win_priv_t *wincon);
int read_mean_3d(master_t *master, double t, int sf, char *var, double *v, 
		 double ***vc);
void reset_means3d_m(master_t *master, int trm);
int read_mean_2d(master_t *master, double t, int sf, char *var, double *v, 
		 double **vc);
void reset_means2d_m(master_t *master, int trm);
void quicksort(double *p, int *q, int l, int r);
int partition(double *p, int *q, int l, int r);
static int create_output(parameters_t *params,   datafile_t *df,
			 char *filename, char** vars, int nvars, 
			 int *s2i, int *s2j, int *s2k);

/*-------------------------------------------------------------------*/
/* Routine to print simulation time details to file                  */
/*-------------------------------------------------------------------*/
void monitor(master_t *master,  /* Model grid data structure         */
	     geometry_t **window, /* Window geometry structure       */
             int mode           /* mode=0, initialise, mode=1 output */
  )
{
  static int prev_elaps;
  double d1 = 0.0, d2, elapsed = 0.0;        /* Dummies              */
  int n;
  FILE *dfp, *fopen();

  /* Only for the master master */
  if (master->mpi_rank)
    return;

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
    if (!(master->decf & NONE)) {
      int dim = (master->decf == DEC_ETA) ? 2 : 3;
      calc_decorr(geom, master->decv, master->decv1, master->decv2, 
		  master->decorr, dim, master->decs);
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
          /* Get the new window sizes */
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
              fprintf(dfp, "%4.5f ", geom->win_size[n]);
            }
            fprintf(dfp, "\nPredicted load balance = ");
            for (n = 1; n <= master->nwindows; n++) {
              size[n] /= d3;
              fprintf(dfp, "%4.5f ", size[n]);
            }
            fprintf(dfp, "\n");
          }
        }
/*UR        d2 =
          (schedule->stop_time - master->t) * d1 / (3600.0 * master->dt);
        if (d2 <= 24.0)
          fprintf(dfp, "Simulation end = %6.3f hours\n", d2);
        else {
          n = (int)(d2 / 24.0);
          fprintf(dfp, "Simulation end = %d (days) %5.2f (hours)\n",
                  n, d2 - (double)(24.0 * n));
        }
*/
        /* Estimated time to completion */
/*UR        d2 = (schedule->stop_time - master->t) * d1 / (3600.0 * master->dt);
        fprintf(dfp, "Time to completion = %d day(s) %2.2d:%2.2d:%2.2d\n",
              (int)(d2/24.0), (int)fmod(d2, 24.0), (int)fmod(d2*60, 60.0),
              (int)fmod(d2*3600, 60.0));
*/
        /* Elapsed time */
        elapsed = difftime(time(NULL), schedule->exec_start_time) / 3600.0;
        fprintf(dfp, "Elapsed time = %d day(s) %2.2d:%2.2d:%2.2d\n",
              (int)(elapsed/24.0), (int)fmod(elapsed, 24.0), (int)fmod(elapsed*60, 60.0),
              (int)fmod(elapsed*3600, 60.0));
        d1 = ((master->t - schedule->start_time)/3600) / elapsed;
        d2 = (schedule->stop_time - master->t)/ 3600.0;
	if (master->nstep > 1 && d1 < HUGE)
	  master->mrtr = (master->mrtr * (double)(master->nstep-1) + d1) / (double)master->nstep;
        fprintf(dfp, "Total time ratio = %.2f (%d:1)\n", d1, (int)master->mrtr);
	if (d1 > 0)
	  d2 = d2 /d1;
        fprintf(dfp, "Time to completion = %d day(s) %2.2d:%2.2d:%2.2d\n",
              (int)(d2/24.0), (int)fmod(d2, 24.0), (int)fmod(d2*60, 60.0),
              (int)fmod(d2*3600, 60.0));
 /*UR      if (d2 <= 24.0)
          fprintf(dfp, "Simulation end = %6.3f hours\n", d2);
        else {
          n = (int)(d2 / 24.0);
          fprintf(dfp, "Simulation end = %d (days) %5.2f (hours)\n",
                  n, d2 - (double)(24.0 * n));
        }
*/

        /* Percent complete */
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
	if (mode == 3)
	  fprintf(dfp, "** CRASHED **\n");
	else {
	  if (master->t == schedule->stop_time)
	    fprintf(dfp, "Run successful.\n");
	  else
	    fprintf(dfp, "Running...\n");
	}
        n = fclose(dfp);
	FLUSH_TIMING
      }
    }
    /*UR-ADDED runtime statistics log, this should only be
     * enabled during development
     */
    if(stat_log)
    {
      int elaps = (int)(elapsed*3600);
      if(master->t - master->dt - schedule->start_time < 1)
      {
        dfp = fopen(stat_logfile, "w");
        if (dfp != NULL)
            fprintf(dfp, "Run time (sec) \t Simulation time (sec) \t complete \n");

      }else if(elaps % STAT_INCREMENT == 0 && elaps > prev_elaps) /* only every 10 seconds */
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
    }/*end stat_log */
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
  int wn, cc, c, cs, cm;
  double max;

  if (wincon->dbc == 0) return;
  if (window->wn != wincon->dbw) return;
  if (wincon->dbgf & D_AFTER && windat->t < wincon->dbgtime) return;
  c = wincon->dbc;
  cs = window->m2d[c];
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
    fprintf(wincon->dbf, "cc = %d\n", cc);
    fprintf(wincon->dbf, "cs = %d\n", cs);
    fprintf(wincon->dbf, "cb = %d\n", window->bot_t[cc]);
    fprintf(wincon->dbf, "cg = %d\n", window->wsa[c]);
    fprintf(wincon->dbf, "xm1 = %d, xp1 = %d, ym1 = %d, yp1 = %d\n", 
	    window->xm1[c], window->xp1[c], window->ym1[c], window->yp1[c]);
    fprintf(wincon->dbf, "depth = %5.2f\n", window->botz[cs]);
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = wincon->s1[cc];
      if (cs == window->m2d[c]) break;
    }
    if (c == 0)c = cs;
    fprintf(wincon->dbf,"\nc       k   eta    u1av   u2av   u1     u2    temp  salt   dz    dzu1  dzu2  gridz\n");
    fprintf(wincon->dbf,"%-7d %-3d %-5.3f %-6.3f %-6.3f %-6.3f %-6.3f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n",
	    c, window->s2k[c], windat->eta[cs], windat->u1av[cs], windat->u2av[cs],
	    windat->u1[c], windat->u2[c], windat->temp[c], windat->sal[c],
	    wincon->dz[c], windat->dzu1[c], windat->dzu2[c], window->gridz[c]);    
    c = window->zm1[c];
    if (c != window->zm1[c]) {
      do {
	fprintf(wincon->dbf,"%-7d %-3d                     %-6.3f %-6.3f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n",
		c, window->s2k[c], windat->u1[c], windat->u2[c], windat->temp[c], 
		windat->sal[c], wincon->dz[c], windat->dzu1[c], windat->dzu2[c],
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
    fprintf(wincon->dbf, "Post-advection U1 = %f\n", windat->nu1[c]);
    fprintf(wincon->dbf, "  u1inter(3D)  = %f\n", wincon->u1inter[cs]);
  }
  if (var == D_U && mode == D_HDIFF) {
    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "u1 horizontal diffusion\n");
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "Post-diffusion U1 = %f\n", windat->nu1[c]);
    fprintf(wincon->dbf, "  u1inter (+2D) = %f\n", wincon->u1inter[cs]);
  }
  if (var == D_U && mode == D_PRESSURE) {
    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "u1 pressure\n");
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "Post-pressure term U1 = %f\n", windat->nu1[c]);
    fprintf(wincon->dbf, "  u1inter (+integral rho) = %f\n", wincon->u1inter[cs]);
  }
  if (var == D_U && mode == D_VZ) {
    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "u1 vertical mixing\n");
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "Post-vertical mixing U1 = %f\n", windat->nu1[c]);
    fprintf(wincon->dbf, "  u1inter (+wind stress) = %f\n", wincon->u1inter[cs]);
  }
  if (var == D_U && mode == D_POST) {
    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "u1 complete\n");
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "New U1 = %f\n", windat->nu1[c]);
    fprintf(wincon->dbf, "  u1inter = %f\n\n", wincon->u1inter[cs]);
  }

  if (var == D_V && mode == D_ADVECT) {
    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "u2 advection\n");
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "Post-advection U2 = %f\n", windat->nu2[c]);
    fprintf(wincon->dbf, "  u2inter (3D) = %f\n", wincon->u2inter[cs]);
  }
  if (var == D_V && mode == D_HDIFF) {
    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "u2 horizontal diffusion\n");
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "Post-diffusion U2 = %f\n", windat->nu2[c]);
    fprintf(wincon->dbf, "  u2inter (+2D) = %f\n", wincon->u2inter[cs]);
  }
  if (var == D_V && mode == D_PRESSURE) {
    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "u2 pressure\n");
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "Post-pressure term U2 = %f\n", windat->nu2[c]);
    fprintf(wincon->dbf, "  u2inter (+integral rho) = %f\n", wincon->u2inter[cs]);
  }
  if (var == D_V && mode == D_VZ) {
    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "u2 vertical mixing\n");
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "Post-vertical mixing U2 = %f\n", windat->nu2[c]);
    fprintf(wincon->dbf, "  u2inter (+wind stress) = %f\n", wincon->u2inter[cs]);
  }
  if (var == D_V && mode == D_POST) {
    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "u2 complete\n");
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "New U2 = %f\n", windat->nu2[c]);
    fprintf(wincon->dbf, "  u2inter = %f\n", wincon->u2inter[cs]);
  }

  if (var == D_UA && mode == D_ADVECT) {
    fprintf(wincon->dbf, "\n2D mode step %d\n", wincon->ic);
    fprintf(wincon->dbf, "---------------\n");
    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "u1av advection\n");
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "Post-advection U1AV = %f\n", windat->nu1av[cs]);
  }
  if (var == D_UA && mode == D_HDIFF) {
    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "u1av horizontal diffusion\n");
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "Post-diffusion U1AV = %f\n", windat->nu1av[cs]);
  }
  if (var == D_UA && mode == D_POST) {
    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "u1av update\n");
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "Barotropic pressure U1AV = %f\n", wincon->b1);
    fprintf(wincon->dbf, "Coriolis U1AV = %f\n", wincon->b2);
    fprintf(wincon->dbf, "Bottom friction U1AV = %f\n", wincon->b3);
    fprintf(wincon->dbf, "Wind stress/baroclinic pressure U1AV = %f\n", wincon->u1inter[cs]);
    fprintf(wincon->dbf, "U1AV + tendencies = %f\n", windat->nu1av[cs]);
  }
  if (var == D_UA && mode == D_BDRY) {
    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "u1av complete\n");
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "New U1AV = %f\n", windat->nu1av[cs]);
    max = 0;
    cm = 1;
    for (cc = 1; cc <= window->b2_e1; cc++) {
      c = window->w2_e1[cc];
      if (fabs(windat->nu1av[c]) > max) {
	max = fabs(windat->nu1av[c]);
	cm = c;
      }
    }
    fprintf(wincon->dbf, "Maximum |u1av| in window%d = %f @ (%d %d)\n\n",
	    window->wn, max, window->s2i[cm], window->s2j[cm]);
    fflush(wincon->dbf);
  }

  if (var == D_VA && mode == D_ADVECT) {
    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "u2av advection\n");
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "Post-advection U2AV = %f\n", windat->nu2av[cs]);
  }
  if (var == D_VA && mode == D_HDIFF) {
    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "u2av horizontal diffusion\n");
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "Post-diffusion U2AV = %f\n", windat->nu2av[cs]);
  }
  if (var == D_VA && mode == D_POST) {
    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "u2av update\n");
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "Barotropic pressure U2AV = %f\n", wincon->b1);
    fprintf(wincon->dbf, "Coriolis U2AV = %f\n", wincon->b2);
    fprintf(wincon->dbf, "Bottom friction U2AV = %f\n", wincon->b3);
    fprintf(wincon->dbf, "Wind stress/baroclinic pressure U2AV = %f\n", wincon->u2inter[cs]);
    fprintf(wincon->dbf, "U2AV + tendencies = %f\n", windat->nu2av[cs]);
  }
  if (var == D_VA && mode == D_BDRY) {
    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "u2av complete\n");
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "New U2AV = %f\n", windat->nu2av[cs]);
    max = 0;
    cm = 1;
    for (cc = 1; cc <= window->b2_e2; cc++) {
      c = window->w2_e2[cc];
      if (fabs(windat->nu2av[c]) > max) {
	max = fabs(windat->nu2av[c]);
	cm = c;
      }
    }
    fprintf(wincon->dbf, "Maximum |u2av| in window%d = %f @ (%d %d)\n\n",
	    window->wn, max, window->s2i[cm], window->s2j[cm]);
    fflush(wincon->dbf);
  }

  if (var == D_ETA && mode == D_POST) {
    double colflux = wincon->b1 - wincon->b2 + wincon->b3 - wincon->b4;
    double neweta = max(windat->etab[c] - colflux / window->cellarea[c],
			window->botz[c]);

    if (wincon->dbgf & D_STEP) {
      fprintf(wincon->dbf, "eta update\n");
      fflush(wincon->dbf);
      return;
    }
    fprintf(wincon->dbf, "          %f\n", windat->u2av[window->yp1[c]]);
    fprintf(wincon->dbf, "%f           %f\n", windat->u1av[c], 
	    windat->u1av[window->xp1[c]]);
    fprintf(wincon->dbf, "          %f\n", windat->u2av[c]);
    fprintf(wincon->dbf, "u1flux[xp1[c]] = %f\n", wincon->b1);
    fprintf(wincon->dbf, "u1flux[c] = %f\n", wincon->b2);
    fprintf(wincon->dbf, "u2flux[yp1[c]] = %f\n", wincon->b3);
    fprintf(wincon->dbf, "u2flux[c] = %f\n", wincon->b4);
    fprintf(wincon->dbf, "colflux = %f\n", colflux);
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
      memset(tendency, 0, window->sgsiz * sizeof(double));
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
  memcpy(ovel, vel, window->sgsiz * sizeof(double));
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
      windat->ceta = i_alloc_1d(window[nb]->sgsiz);
      windat->cu1 = i_alloc_1d(window[nb]->sgsiz);
      windat->cu2 = i_alloc_1d(window[nb]->sgsiz);
      windat->cu1a = i_alloc_1d(window[nb]->sgsizS);
      windat->cu2a = i_alloc_1d(window[nb]->sgsizS);
      windat->tempb = d_alloc_1d(window[nb]->sgsiz);
      windat->salb = d_alloc_1d(window[nb]->sgsiz);
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
      cg = master->cu1;
      i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
      fprintf(afp, "%s %f at %d (%d %d %d)\n", "u1 3D     : ",
	      master->mu1, cg, i, j, k);
      master->su1a = max(master->su1a, master->mu1a);
      cg = master->cu1a;
      i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
      fprintf(afp, "%s %f at %d (%d %d %d)\n", "u1 2D     : ",
	      master->mu1a, cg, i, j, k);
      master->su2 = max(master->su2, master->mu2);
      cg = master->cu2;
      i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
      fprintf(afp, "%s %f at %d (%d %d %d)\n", "u2 3D     : ",
	      master->mu2, cg, i, j, k);
      master->su2a = max(master->su2a, master->mu2a);
      cg = master->cu2a;
      i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
      fprintf(afp, "%s %f at %d (%d %d %d)\n", "u2 2D     : ",
	      master->mu2a, cg, i, j, k);
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
	cg = master->ca1;
	i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
	fprintf(afp, "%s %f at %d (%d %d %d)\n", "u1 adv    : ",
		master->ma1, cg, i, j, k);
      }
      if (master->mh1 > 0.0) {
	master->sh1 = max(master->sh1, master->mh1);
	cg = master->ch1;
	i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
	fprintf(afp, "%s %f at %d (%d %d %d)\n", "u1 hdiff  : ",
		master->mh1, cg, i, j, k);
      }
      if (master->mv1 > 0.0) {
	master->sv1 = max(master->sv1, master->mv1);
	cg = master->cv1;
	i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
	fprintf(afp, "%s %f at %d (%d %d %d)\n", "u1 vdiff  : ",
		master->mv1, cg, i, j, k);
      }
      if (master->mb1 > 0.0) {
	master->sb1 = max(master->sa1, master->mb1);
	cg = master->cb1;
	i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
	fprintf(afp, "%s %f at %d (%d %d %d)\n", "u1 btp    : ",
		master->mb1, cg, i, j, k);
      }
      if (master->md1 > 0.0) {
	master->sd1 = max(master->sd1, master->md1);
	cg = master->cd1;
	i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
	fprintf(afp, "%s %f at %d (%d %d %d)\n", "u1 bcp    : ",
		master->md1, cg, i, j, k);
      }
      if (master->mc1 > 0.0) {
	master->sc1 = max(master->sc1, master->mc1);
	cg = master->cc1;
	i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
	fprintf(afp, "%s %f at %d (%d %d %d)\n", "u1 cor    : ",
		master->mc1, cg, i, j, k);
      }
      fprintf(afp, "u2 velocity\n");
      if (master->ma2 > 0.0) {
	master->sa2 = max(master->sa2, master->ma2);
	cg = master->ca2;
	i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
	fprintf(afp, "%s %f at %d (%d %d %d)\n", "u2 adv    : ",
		master->ma2, cg, i, j, k);
      }
      if (master->mh2 > 0.0) {
	master->sh2 = max(master->sh2, master->mh2);
	cg = master->ch2;
	i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
	fprintf(afp, "%s %f at %d (%d %d %d)\n", "u2 hdiff  : ",
		master->mh2, cg, i, j, k);
      }
      if (master->mv2 > 0.0) {
	master->sv2 = max(master->sv2, master->mv2);
	cg = master->cv2;
	i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
	fprintf(afp, "%s %f at %d (%d %d %d)\n", "u2 vdiff  : ",
		master->mv2, cg, i, j, k);
      }
      if (master->mb2 > 0.0) {
	master->sb2 = max(master->sa2, master->mb2);
	cg = master->cb2;
	i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
	fprintf(afp, "%s %f at %d (%d %d %d)\n", "u2 btp    : ",
		master->mb2, cg, i, j, k);
      }
      if (master->md2 > 0.0) {
	master->sd2 = max(master->sd2, master->md2);
	cg = master->cd2;
	i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
	fprintf(afp, "%s %f at %d (%d %d %d)\n", "u2 bcp    : ",
		master->md2, cg, i, j, k);
      }
      if (master->mc2 > 0.0) {
	master->sc2 = max(master->sc2, master->mc2);
	cg = master->cc2;
	i = geom->s2i[cg]; j = geom->s2j[cg]; k = geom->s2k[cg];
	fprintf(afp, "%s %f at %d (%d %d %d)\n", "u2 cor    : ",
		master->mc2, cg, i, j, k);
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
      fprintf(afp, "u2 3D    : %f\n", master->su2);
      fprintf(afp, "u2 2D    : %f\n", master->su2a);
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
  int nb, cc, c;               /* Counters                           */
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
	memcpy(eta, windat->etam, window->sgsizS * sizeof(double));
      else {
	memcpy(eta, wincon->neweta, window->sgsizS * sizeof(double));
      }
      /* Remove the relaxation elevation if required                 */
      if (rlx_diff) {
	for (cc = 1; cc <= window->b2_t; cc++) {
	  c = window->w2_t[cc];
	  eta[c] -= windat->eta_rlx->val1[c];
	}
      }
      /* Find the maximum                                            */
      memset(windat->ceta, 0, window->sgsizS * sizeof(int));
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
	  if (wincon->etarlx & ALERT) {

	    /* Relax hard to elevation                               */
	    do_eta_relax(window, cl, relax_rate, tide_r);

	    /* Increase horizontal viscosity locally                 */
	    wincon->alertf |= LEVEL2;

	    /*
	    shapiro_smoothx(window, windat->u1av ,windat->cdeta);
	    shapiro_smoothy(window, windat->u2av, windat->cdeta);
	    shapiro_smoothx(window, windat->u1av ,window->xp1[windat->cdeta]);
	    shapiro_smoothy(window, windat->u2av, window->yp1[windat->cdeta]);

	    shuman_smooth(window, windat->u1av, windat->cdeta);
	    shuman_smooth(window, windat->u2av, windat->cdeta);
	    shuman_smooth(window, windat->u1av, window->xp1[windat->cdeta]);
	    shuman_smooth(window, windat->u2av, window->yp1[windat->cdeta]);
	    */
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
      memset(windat->cu1a, 0, window->sgsizS * sizeof(int));
      windat->mu1a = get_maxs(window, windat->u1av, 1, window->b2_e1,
			     window->w2_e1, windat->cu1a, wincon->velmax2d);
      memset(windat->cu2a, 0, window->sgsizS * sizeof(int));
      windat->mu2a = get_maxs(window, windat->u2av, 1, window->b2_e2,
			     window->w2_e2, windat->cu2a, wincon->velmax2d);

      /* Get the mechanical energy, excess mass and energy flux      */
      /* through the open boundaries (Palma & Matano, 1998, section  */
      /* 2.3 or 2001 eqn 6).                                         */
      /* Note: Energy flux in m4s-3 is the flux in Wm-1/rho          */
      /* Wm-2 = kgs-3                                                */
      /* Jm-2 = kgs-2                                                */
      for (cc = 1; cc <= window->b2_t; cc++) {
	c = window->w2_t[cc];
	/*if(window->s2i[c] < 18 || window->s2i[c] > 38 || window->s2j[c] > 18) continue;*/
	d2 = windat->eta[c] - window->botz[c];
	d3 = window->cellarea[c];
	ta += d3;
	d1 = 0.5 * (windat->u1av[c] * windat->u1av[c] +
		    windat->u2av[c] * windat->u2av[c]);
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
	if (open->type & U1BDRY) {
	  cw = window->h2au1;
	  vel = windat->u1av;
	  ce = open->no2_e1;
	  cells = open->obc_e1;
	} else {
	  cw = window->h1au2;
	  vel = windat->u2av;
	  ce = open->no2_e2;
	  cells = open->obc_e2;
	}
	ef = d3 = d1 = 0.0;
	meta = mvel = 0.0;
	sgn = (open->ocodex & L_EDGE || open->ocodey & B_EDGE) ? 1.0 : -1.0;

	for (cc = 1; cc <= ce; cc++) {
	  c = cells[cc];
	  /*if(window->s2j[c] > 18) continue;*/
	  d3 += cw[c];
	  ef += sgn * cw[c] * vel[c] * (windat->eta[c] - 
					window->botz[c] * wincon->Ds[c]) *
	    (wincon->g * windat->eta[c] +
	     0.5 * (windat->u1av[c] * windat->u1av[c] +
		    windat->u2av[c] * windat->u2av[c]));
	  meta += cw[c] * windat->eta[c];
	  mvel += vel[c] * cw[c] * (windat->eta[c] - window->botz[c] * wincon->Ds[c]);
	  d1 += cw[c] * (windat->eta[c] - window->botz[c] * wincon->Ds[c]);
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
	  memset(windat->alert_a, 0, window->sgsizS * sizeof(double));

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
	  for(cc = 1; (c = windat->cu1a[cc]); cc++) {
	    windat->u1av[c] = con_median(window, windat->u1av, wincon->velmax2d, c);
	    /*
	    windat->u1av[c] = (windat->u1av[c] > 0.0) ?
	                       wincon->velmax2d : -wincon->velmax2d;
	    */
	    /*shuman_smooth(window, windat->u1av, c);*/
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
	  int ym1, cl = (windat->mu1a > wincon->velmax2d) ? cc : 1;
	  double shear;
	  for (cc = 1; cc <= window->v2_e1; cc++) {
	    c = window->w2_e1[cc];
	    ym1 = window->ym1[c];
	    shear = (windat->u1av[c] - windat->u1av[ym1]) / window->h1au2[c];
	    if (fabs(shear) > wincon->smax) {
	      windat->cu1a[cl] = c;
	      cl++;
	      wincon->alertf |= LEVEL1_UA;
	    }
	    if (fabs(shear) > windat->mshear)
	      windat->mshear = fabs(shear);
	  }
	}

	if (wincon->vel2d_f && windat->mu2a > wincon->velmax2d) {
	  /*    => cap the velocity to velmax                        */
	  for(cc = 1; (c = windat->cu2a[cc]); cc++) {

	    windat->u2av[c] = con_median(window, windat->u2av, wincon->velmax2d, c);
	    /*
	    windat->u2av[c] = (windat->u2av[c] > 0.0) ?
	                       wincon->velmax2d : -wincon->velmax2d;
	    */
	    /*shuman_smooth(window, windat->u2av, c);*/
	    windat->alert_a[c] = 3.0;
	    windat->alert_c[c] += 1.0;
	  }
	  /*    => increase horizontal viscosity locally             */
	  wincon->alertf |= LEVEL1_VA;
	  hd_warn("LEVEL1 alert at %f days : u2av exceeds threshold\n", t);
	  windat->nalert[1]++;
	}
	if (wincon->vel2d_f && wincon->shear_f) {
	  int xm1, cl = (windat->mu2a > wincon->velmax2d) ? cc : 1;
	  double shear;
	  for (cc = 1; cc <= window->v2_e2; cc++) {
	    c = window->w2_e2[cc];
	    xm1 = window->xm1[c];
	    shear = (windat->u2av[c] - windat->u2av[xm1]) / window->h2au1[c];
	    if (fabs(shear) > wincon->smax) {
	      windat->cu2a[cl] = c;
	      cl++;
	      wincon->alertf |= LEVEL1_VA;
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
      wincon->u2_f = NONE;

      /* Get the maximum absolute velocity                           */
      memset(windat->cu1, 0, window->sgsiz * sizeof(int));
      windat->mu1 = get_maxs(window, windat->u1, 1, window->b3_e1,
			     window->w3_e1, windat->cu1, wincon->velmax);
      memset(windat->cu2, 0, window->sgsiz * sizeof(int));
      windat->mu2 = get_maxs(window, windat->u2, 1, window->b3_e2,
			     window->w3_e2, windat->cu2, wincon->velmax);
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
      if (windat->u2_adv) {
	windat->ma2 = get_max(window, windat->u2_adv, 1, window->b3_e2,
			      window->w3_e2, &windat->ca2);
      }
      if (windat->u2_hdif) {
	windat->mh2 = get_max(window, windat->u2_hdif, 1, window->b3_e2,
			      window->w3_e2, &windat->ch2);
      }
      if (windat->u2_vdif) {
	windat->mv2 = get_max(window, windat->u2_vdif, 1, window->b3_e2,
			      window->w3_e2, &windat->cv2);
      }
      if (windat->u2_btp) {
	windat->mb2 = get_max(window, windat->u2_btp, 1, window->b3_e2,
			      window->w3_e2, &windat->cb2);
      }
      if (windat->u2_bcp) {
	windat->md2 = get_max(window, windat->u2_bcp, 1, window->b3_e2,
			      window->w3_e2, &windat->cd2);
      }
      if (windat->u2_cor) {
	windat->mc2 = get_max(window, windat->u2_cor, 1, window->b3_e2,
			      window->w3_e2, &windat->cc2);
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
	  for(cc = 1; (c = windat->cu1[cc]); cc++) {
	    windat->u1[c] = con_median(window, windat->u1, wincon->velmax, c);
	    /*
	    windat->u1[c] = (windat->u1[c] > 0.0) ?
	                     wincon->velmax : -wincon->velmax;
	    */
	    /*shuman_smooth(window, windat->u1, c);*/
	    windat->alert_a[window->m2d[c]] = 4.0;
	    windat->alert_c[window->m2d[c]] += 1.0;
	  }
	  /*    => increase horizontal viscosity locally             */
	  wincon->alertf |= LEVEL1_U;
	  hd_warn("LEVEL1 alert at %f days : u1 exceeds threshold\n", t);
	  windat->nalert[2]++;
	}
	if (wincon->vel3d_f && wincon->shear_f) {
	  int ym1, cl = (windat->mu1 > wincon->velmax) ? cc : 1;
	  double shear;
	  for (cc = 1; cc <= window->v3_e1; cc++) {
	    c = window->w3_e1[cc];
	    ym1 = window->ym1[c];
	    shear = (windat->u1[c] - windat->u1[ym1]) /
	      window->h1au2[window->m2d[c]];
	    if (fabs(shear) > wincon->smax) {
	      windat->cu1[cl] = c;
	      cl++;
	      wincon->alertf |= LEVEL1_U;
	    }
	    if (fabs(shear) > windat->mshear)
	      windat->mshear = fabs(shear);
	  }
	}

	if (wincon->vel3d_f && windat->mu2 > wincon->velmax) {
	  /*    => cap the velocity to velmax                        */
	  for(cc = 1; (c = windat->cu2[cc]); cc++) {
	    windat->u2[c] = con_median(window, windat->u2, wincon->velmax, c);
	    /*
	    windat->u2[c] = (windat->u2[c] > 0.0) ?
	                     wincon->velmax : -wincon->velmax;
	    */
	    /*shuman_smooth(window, windat->u2, c);*/
	    windat->alert_a[window->m2d[c]] = 5.0;
	    windat->alert_c[window->m2d[c]] += 1.0;
	  }
	  /*    => increase horizontal viscosity locally             */
	  wincon->alertf |= LEVEL1_V;
	  hd_warn("LEVEL1 alert at %f days : u2 exceeds threshold\n", t);
	  windat->nalert[2]++;
	}

	if (wincon->vel3d_f && wincon->shear_f) {
	  int xm1, cl = (windat->mu2 > wincon->velmax) ? cc : 1;
	  double shear;
	  for (cc = 1; cc <= window->v3_e2; cc++) {
	    c = window->w3_e2[cc];
	    xm1 = window->xm1[c];
	    shear = (windat->u2[c] - windat->u2[xm1]) /
	      window->h2au1[window->m2d[c]];
	    if (fabs(shear) > wincon->smax) {
	      windat->cu2[cl] = c;
	      cl++;
	      wincon->alertf |= LEVEL1_V;
	    }
	    if (fabs(shear) > windat->mshear)
	      windat->mshear = fabs(shear);
	  }
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
	  if (windat->ma2 > wincon->amax) {
	    wincon->u2_f |= ADVECT;
	    wincon->u2av_f |= ADVECT;
	    hd_warn("LEVEL3 alert at %f days : u2 adv tendency exceeds threshold\n", t);
	    windat->nalert[4]++;
	  }
	  if (windat->mh2 > wincon->hmax) {
	    wincon->u2_f |= HDIFF;
	    wincon->u2av_f |= HDIFF;
	    hd_warn("LEVEL3 alert at %f days : u2 hdif tendency exceeds threshold\n", t);
	    windat->nalert[4]++;
	  }
	  if (windat->mv2 > wincon->vmax) {
	    wincon->u2_f |= VDIFF;
	    wincon->u2av_f |= VDIFF;
	    hd_warn("LEVEL3 alert at %f days : u2 vdif tendency exceeds threshold\n", t);
	    windat->nalert[4]++;
	  }
	  if (windat->mb2 > wincon->btmax) {
	    wincon->u2_f |= PRESS_BT;
	    wincon->u2av_f |= PRESS_BT;
	    hd_warn("LEVEL3 alert at %f days : u2 btp tendency exceeds threshold\n", t);
	    windat->nalert[4]++;
	  }
	  if (windat->md2 > wincon->bcmax) {
	    wincon->u2_f |= PRESS_BC;
	    wincon->u2av_f |= PRESS_BC;
	    hd_warn("LEVEL3 alert at %f days : u2 bcp tendency exceeds threshold\n", t);
	    windat->nalert[4]++;
	  }
	  if (windat->mc2 > wincon->cmax) {
	    wincon->u2_f |= CORIOLIS;
	    wincon->u2av_f |= CORIOLIS;
	    hd_warn("LEVEL3 alert at %f days : u2 cor tendency exceeds threshold\n", t);
	    windat->nalert[4]++;
	  }
	}
	/* 3. Divergences exceeds thresholds                         */
	/*    => damp or smooth velocity field                       */
	if (wincon->div3d_f && windat->mdw > wincon->dwmax) {
	  /*    => increase horizontal viscosity locally             */
	  wincon->alertf |= LEVEL2w;
	  /*
	  shapiro_smoothx(window, windat->u1, windat->cdw);
	  shapiro_smoothy(window, windat->u2, windat->cdw);
	  shapiro_smoothx(window, windat->u1, window->xp1[windat->cdw]);
	  shapiro_smoothy(window, windat->u2, window->yp1[windat->cdw]);

	  shuman_smooth(window, windat->u1, windat->cdw);
	  shuman_smooth(window, windat->u2, windat->cdw);
	  */
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
    memcpy(eta, wincon->neweta, window->sgsizS * sizeof(double));
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
  sst = stencil(window, cl, &ssize, 1);
  for (cc = 0; cc < ssize; cc++) {
    c = sst[cc];
    st = stencil(window, c, &size, 2);

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

  st = stencil(window, c, &size, 2);
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
  int c, cc;
  int xm1, xp1, ym1, yp1;
  double *a = wincon->w5;

  memcpy(a, var, window->sgsiz * sizeof(double));
  for (cc = 1; cc < window->b3_t; cc++) {
    c = window->w3_t[cc];
    xm1 = window->xm1[c];
    xp1 = window->xp1[c];
    ym1 = window->ym1[c];
    yp1 = window->yp1[c];
    var[c] = a[c] + 0.5 * v * (1.0 - v) *
      (a[xm1] + a[xp1] + a[ym1] + a[yp1] - 4.0 * a[c]) +
      0.25 * v * v * (a[window->xp1[ym1]] + a[window->xm1[ym1]] + 
		      a[window->xp1[yp1]] + a[window->xm1[yp1]] -
		      4.0 * a[c]);
  }
}

/* END shuman_3d()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to implement a Shapiro type smoothing on an array in the  */
/* x direction (see Kowalik and Murty, p146, Eqn 3.136).             */
/*-------------------------------------------------------------------*/
void shapiro_smoothx(geometry_t *window,  /* Window geometry         */
		     double *a,     /* Array to smooth               */
		     int cl         /* Location to smooth            */
		     )
{
  int xm1 = window->xm1[cl];
  int xp1 = window->xp1[cl];
  int xm2 = window->xm1[xm1];
  int xp2 = window->xp1[xp1];
  int xm3 = window->xm1[xm2];
  int xp3 = window->xp1[xp2];
  int xm4 = window->xm1[xm3];
  int xp4 = window->xp1[xp3];

  a[cl] = (186 * a[cl] + 56 * (a[xm1] + a[xp1]) -
	   28 * (a[xm2] + a[xp2]) + 8 * (a[xm3] + a[xp3]) -
	   (a[xm4] + a[xp4])) / 256.0;
}

double shapiro_smoothxr(geometry_t *window,  /* Window geometry         */
			double *a,     /* Array to smooth               */
			int cl         /* Location to smooth            */
		     )
{
  double ret;
  int xm1 = window->xm1[cl];
  int xp1 = window->xp1[cl];
  int xm2 = window->xm1[xm1];
  int xp2 = window->xp1[xp1];
  int xm3 = window->xm1[xm2];
  int xp3 = window->xp1[xp2];
  int xm4 = window->xm1[xm3];
  int xp4 = window->xp1[xp3];

  ret = (186 * a[cl] + 56 * (a[xm1] + a[xp1]) -
	   28 * (a[xm2] + a[xp2]) + 8 * (a[xm3] + a[xp3]) -
	   (a[xm4] + a[xp4])) / 256.0;

  return(ret);
}

/* END shapiro_smoothx()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to implement a Shapiro type smoothing on an array in the  */
/* y direction (see Kowalik and Murty, p146, Eqn 3.136).             */
/*-------------------------------------------------------------------*/
void shapiro_smoothy(geometry_t *window,  /* Window geometry         */
		     double *a,     /* Array to smooth               */
		     int cl         /* Location to smooth            */
		     )
{
  int ym1 = window->ym1[cl];
  int yp1 = window->yp1[cl];
  int ym2 = window->ym1[ym1];
  int yp2 = window->yp1[yp1];
  int ym3 = window->ym1[ym2];
  int yp3 = window->yp1[yp2];
  int ym4 = window->ym1[ym3];
  int yp4 = window->yp1[yp3];

  a[cl] = (186 * a[cl] + 56 * (a[ym1] + a[yp1]) -
	   28 * (a[ym2] + a[yp2]) + 8 * (a[ym3] + a[yp3]) -
	   (a[ym4] + a[yp4])) / 256.0;
}

double shapiro_smoothyr(geometry_t *window,  /* Window geometry         */
			double *a,     /* Array to smooth               */
			int cl         /* Location to smooth            */
			)
{
  double ret;
  int ym1 = window->ym1[cl];
  int yp1 = window->yp1[cl];
  int ym2 = window->ym1[ym1];
  int yp2 = window->yp1[yp1];
  int ym3 = window->ym1[ym2];
  int yp3 = window->yp1[yp2];
  int ym4 = window->ym1[ym3];
  int yp4 = window->yp1[yp3];

  ret = (186 * a[cl] + 56 * (a[ym1] + a[yp1]) -
	   28 * (a[ym2] + a[yp2]) + 8 * (a[ym3] + a[yp3]) -
	   (a[ym4] + a[yp4])) / 256.0;

  return(ret);
}

/* END shapiro_smoothy()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to implement a Shapiro type smoothing on an array. Taken  */
/* Sapiro, R. (1975) Linear Filtering. Mathematics of Computation,   */
/* 29, 1094-1097. Stencils are taken from Table 2. For 3D arrays.    */
/*-------------------------------------------------------------------*/
void shapiro(geometry_t *window,     /* Window geometry              */
	     double *a,              /* Array to smooth              */
	     double *buf,            /* Dummy buffer                 */
	     int *vec,               /* Cells to process             */
	     int nvec,               /* Size of vec                  */
	     int order,              /* Order of the filter          */
	     int dir                 /* Smoothing direction          */
	     )
{
  int c, cc;
  int cp, cm, n;
  int *p1, *m1;
  double ret, sgn;

  order -= 1;
  memcpy(buf, a, window->sgsiz * sizeof(double));
  if (dir & XDIR) {
    p1 = window->xp1;
    m1 = window->xm1;
  } else if (dir & YDIR) {
    p1 = window->yp1;
    m1 = window->ym1;
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
	     int type               /* Type of stencil               */
	     )
{
  int cc, c, m = 0;
  int *st = NULL;
  int sz = floor(*size / 2);
  int ss = type;

  /* Star shaped stencil. Stencil filled from m=0; m < ss;           */
  if (ss == 1) {
    int n, s, e, w, mm;
    for (cc = 1; cc <= sz; cc++)
      ss += 4 * cc;
    st = i_alloc_1d(ss);
    *size = ss;

    /* Get the furthest locations from the center                    */
    n = s = e = w = cl;
    for (cc = 1; cc <= sz; cc++) {
      n = window->yp1[n];
      s = window->ym1[s];
      e = window->xp1[e];
      w = window->xm1[w];
    }
    st[m++] = n; st[m++] = s; st[m++] = e; st[m++] = w;

    /* Get the remaining locations in order of decreasing distance  */
    for(mm = sz - 1; mm >= 1; mm--) {
      n = s = e = w = cl;
      for(cc = 1; cc <= mm; cc++) {
	n = window->yp1[n];
	s = window->ym1[s];
	e = window->xp1[e];
	w = window->xm1[w];
      }
      st[m++] = n; st[m++] = s; st[m++] = e; st[m++] = w;
      n = window->xm1[n];
      s = window->xp1[s];
      e = window->yp1[e];
      w = window->ym1[w];
      st[m++] = n; st[m++] = s; st[m++] = e; st[m++] = w;
      for(cc = 1; cc < mm; cc++) {
	n = window->ym1[window->xm1[n]];
	s = window->yp1[window->xp1[s]];
	e = window->yp1[window->xm1[e]];
	w = window->ym1[window->xp1[w]];
	st[m++] = n; st[m++] = s; st[m++] = e; st[m++] = w;
      }
    }
    st[m++] = cl;
  } else if(ss == 2) {
    /* 3 point square stencil ordered for Shuman filter              */
    *size = 9;
    st = i_alloc_1d(*size * *size);
    st[0] = cl;
    st[1] = window->xm1[cl]; st[2] = window->xp1[cl];
    st[3] = window->ym1[cl]; st[4] = window->yp1[cl];
    st[5] = window->xp1[window->ym1[cl]];
    st[6] = window->xm1[window->ym1[cl]];
    st[7] = window->xp1[window->yp1[cl]];
    st[8] = window->xm1[window->yp1[cl]];
  } else if(ss == 0) {
    /* Square stencil                                                */
    int mm, cs = cl;
    st = i_alloc_1d(*size * *size);
    /* Get the bottom left corner of the stencil                     */
    for (cc = 1; cc < sz; cc++)
      cs = window->xm1[window->ym1[cs]];

    /* Save the stencil clocations                                   */
    for (mm = 0; mm < *size; mm++) {
      c = cs;
      for (cc = 0; cc < *size; cc++) {
	st[m++] = c;
	c = window->xp1[c];
      }
      cs = window->yp1[cs];
    }
    *size *= *size;
  }

  return(st);
}

/* END stencil()                                                     */
/*-------------------------------------------------------------------*/


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
  int cc, c, cs;

  /*-----------------------------------------------------------------*/
  /* Get the fluxes at time t over the whole window                  */
  memset(u1flux, 0, window->sgsiz * sizeof(double));
  memset(u2flux, 0, window->sgsiz * sizeof(double));
  /*-----------------------------------------------------------------*/
  /* Fluxes through the e1 faces                                     */
  for (cc = 1; cc <= window->a3_e1; cc++) {
    c = window->w3_e1[cc];
    cs = window->m2d[c];
    u1flux[c] = windat->u1[c] * windat->dzu1[c] * window->h2au1[cs] *
      wincon->mdx[cs];
  }
  /* Fluxes through the e2 faces                                     */
  for (cc = 1; cc <= window->a3_e2; cc++) {
    c = window->w3_e2[cc];
    cs = window->m2d[c];
    u2flux[c] = windat->u2[c] * windat->dzu2[c] * window->h1au2[cs] *
      wincon->mdy[cs];
  }

  /*-----------------------------------------------------------------*/
  /* Get the divergence of velocity transport                        */
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    cs = window->m2d[c];

    colflux = u1flux[c] - u1flux[window->xp1[c]] +
      u2flux[c] - u2flux[window->yp1[c]] + window->cellarea[cs] *
      (windat->w[window->zp1[c]] - windat->w[c]);
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
  double *u2flux = wincon->d3;
  double mdiv = 0.0;
  int cc, c;

  /*-----------------------------------------------------------------*/
  /* Get the fluxes at time t over the whole window                  */
  memset(u1flux, 0, window->sgsizS * sizeof(double));
  memset(u2flux, 0, window->sgsizS * sizeof(double));
  /* Calculate the flux at e1 wet and boundary cells                 */
  for (cc = 1; cc <= window->a2_e1; cc++) {
    c = window->w2_e1[cc];
    u1flux[c] = windat->u1av[c] * windat->depth_e1[c] * window->h2au1[c] *
      wincon->mdx[c];
  }
  /* Calculate the flux at e2 wet and boundary cells                 */
  for (cc = 1; cc <= window->a2_e2; cc++) {
    c = window->w2_e2[cc];
    u2flux[c] = windat->u2av[c] * windat->depth_e2[c] * window->h1au2[c] *
      wincon->mdy[c];
  }
  /*-----------------------------------------------------------------*/
  /* Get the divergence of velocity transport                        */
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->w2_t[cc];
    colflux[c] = (u1flux[window->xp1[c]] - u1flux[c] +
		  u2flux[window->yp1[c]] - u2flux[c]) / window->cellarea[c];
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
/* Calculates the divergence of 2D fluxes (equal to detadt) at a     */
/* given location (units of ms-1).                                   */
/*-------------------------------------------------------------------*/
double get_div_2d(geometry_t *window,      /* Window geometry        */
		  window_t *windat,        /* Window data            */
		  win_priv_t *wincon,      /* Window constants       */
		  int c                    /* Sparse coordinate      */
)
{
  double u1flux1, u1flux2;
  double u2flux1, u2flux2;
  double colflux;
  int xp1 = window->xp1[c];
  int yp1 = window->yp1[c];

  /*-----------------------------------------------------------------*/
  /* Get the fluxes at time t                                        */
  u1flux1 = windat->u1av[c] * windat->depth_e1[c] * window->h2au1[c] *
    wincon->mdx[c];
  u1flux2 = windat->u1av[xp1] * windat->depth_e1[xp1] * window->h2au1[xp1] *
    wincon->mdx[xp1];
  u2flux1 = windat->u2av[c] * windat->depth_e2[c] * window->h1au2[c] *
    wincon->mdy[c] * windat->dt2d;
  u2flux2 = windat->u2av[yp1] * windat->depth_e2[yp1] * window->h1au2[yp1] *
    wincon->mdy[yp1] * windat->dt2d;

  /*-----------------------------------------------------------------*/
  /* Get the divergence of velocity transport                        */
  colflux = (u1flux2 - u1flux1 + u2flux2 - u2flux1) / window->cellarea[c];
  /*
  colflux = wincon->d2[window->xp1[c]] - wincon->d2[c] +
    wincon->d3[window->yp1[c]] - wincon->d3[c];
  */
  return(colflux);
}

/* END get_div_2d()                                                  */
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
  int lc, c, cc, zp1;           /* Sparse locations                  */
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
	spd = sqrt(windat->u1[lc] * windat->u1[lc] +
		   windat->u2[lc] * windat->u2[lc]);
	if (spd > cspd)
	  cspd = spd;
	u = fabs(0.5 * (windat->u1[window->xp1[lc]] + windat->u1[lc]));
	u = (u) ? window->h1acell[c] / u : HUGE;
	v = fabs(0.5 * (windat->u2[window->yp1[lc]] + windat->u2[lc]));
	v = (v) ? window->h2acell[c] / v : HUGE;
	w = fabs(0.5 * (windat->w[window->zp1[lc]] + windat->w[lc]));
	w = (w) ? wincon->dz[lc] / w : HUGE;
	windat->cour[c] = min(u, v);
	windat->cour[c] = min(windat->cour[c], w);
	u = fabs(windat->u1[window->xp1[lc]] - windat->u1[lc]);
	u = (u) ? window->h1acell[c] / u : HUGE;
	v = fabs(windat->u2[window->yp1[lc]] - windat->u2[lc]);
	v = (v) ? window->h2acell[c] / v : HUGE;
	w = fabs(windat->w[window->zp1[lc]] - windat->w[lc]);
	w = (w) ? wincon->dz[lc] / w : HUGE;
	windat->lips[c] = min(u, v);
	windat->lips[c] = min(windat->lips[c], w);
	dh = 1.0 / (1.0 / (window->h1acell[c] * window->h1acell[c]) +
		    1.0 / (window->h2acell[c] * window->h2acell[c]));
	dh = (wincon->u1kh[lc]) ? dh / (4.0 * wincon->u1kh[lc]) : 0.0;
	windat->ahsb[c] = max(windat->ahsb[c], dh);
      }
      lc = window->zm1[lc];
    }

    /* Get the grid spacing                                          */
    dh = sqrt(1.0 / (window->h1acell[c] * window->h1acell[c]) +
              1.0 / (window->h2acell[c] * window->h2acell[c]));

    /* Total 3D speed                                                */
    spd = 2.0 * iws + cspd;
    /* Get the gravity wave + 2D current speed                       */
    spd2d = 2.0 * sqrt(-wincon->g * window->botz[c]) +
            sqrt(windat->u1av[c] * windat->u1av[c] +
                 windat->u2av[c] * windat->u2av[c]);

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
  int c, cc, cs, cb;            /* Sparse coordinates                */
  int zp1, zm1;                 /* Sparse maps                       */
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
  for (cc = 1; cc <= window->b2_t; cc++) {
    cs = window->w2_t[cc];
    dh[cs] = sqrt(1.0 / (window->h1acell[cs] * window->h1acell[cs]) +
		  1.0 / (window->h2acell[cs] * window->h2acell[cs]));
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
      cspd = sqrt(windat->u1[c] * windat->u1[c] +
		  windat->u2[c] * windat->u2[c]);

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
/* Calculates diagnostic numbers and saves to tracers if required    */
/*-------------------------------------------------------------------*/
void diag_numbers(geometry_t *window,       /* Window geometry       */
		  window_t *windat,         /* Window data           */
		  win_priv_t *wincon        /* Window constants      */
		  )
{
  int c, cc, cs;
  int xp1, yp1;
  int zm1;
  double depth;  /* Water depth                                      */
  double dzface; /* Cell thickness between cell centers              */
  double speed;  /* Current speed                                    */
  double dir;    /* Current direction                                */
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
      xp1 = window->xp1[c];
      yp1 = window->yp1[c];
      du1 = 0.5 * (windat->u1[c] + windat->u1[xp1]);
      du2 = 0.5 * (windat->u2[c] + windat->u2[yp1]);
      windat->speed_sq[cs] = du1 * du1 + du2 * du2;
    }
  }

  if (windat->tau_bm) {
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      xp1 = window->xp1[c];
      yp1 = window->yp1[c];
      du1 = 0.5 * (windat->tau_be1[c] + windat->tau_be1[xp1]);
      du2 = 0.5 * (windat->tau_be2[c] + windat->tau_be2[yp1]);
      windat->tau_bm[c] = sqrt(du1 * du1 + du2 * du2);
    }
  }

  if (windat->slope_x && windat->slope_y) {
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      windat->slope_x[c] = (windat->eta[window->xp1[c]] - 
			    windat->eta[window->xm1[c]]) / 
	(window->h1au1[c] + window->h1au1[window->xp1[c]]);
      windat->slope_y[c] = (windat->eta[window->yp1[c]] - 
			    windat->eta[window->ym1[c]]) / 
	(window->h2au2[c] + window->h2au2[window->yp1[c]]);
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

  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    cs = window->m2d[c];
    xp1 = window->xp1[c];
    yp1 = window->yp1[c];
    zm1 = window->zm1[c];

    /* Get the cell thickness between cell faces                     */
    dzface = 0.5 * (wincon->dz[zm1] + wincon->dz[c]);
    depth = windat->eta[cs] - window->botz[cs];

    du1 = 0.5 * (windat->u1[c] + windat->u1[xp1]);
    du2 = 0.5 * (windat->u2[c] + windat->u2[yp1]);
    speed = sqrt(du1 * du1 + du2 * du2);
    dir =  fmod(atan2(du1, du2) * 180.0 / M_PI + 180.0, 360.0);
    if (windat->speed_3d) {
      windat->speed_3d[c] = speed;
      windat->speedd_3d[c] = dir;
    }
    if (windat->speed_2d && c == cs) {
      du1 = 0.5 * (windat->u1av[cs] + windat->u1av[xp1]);
      du2 = 0.5 * (windat->u2av[cs] + windat->u2av[yp1]);
      speed = sqrt(du1 * du1 + du2 * du2);
      dir =  fmod(atan2(du1, du2) * 180.0 / M_PI + 180.0, 360.0);
      windat->speed_2d[cs] = speed;
      windat->speedd_2d[c] = dir;
    }

    /* Get the Brunt Vaisala frequency                               */
    /* N2 = -(g/rho)(drho/dz) (s-1)                                  */
    N2 = buoyancy_frequency2(window, windat, wincon, c);
    N = N2 > 0 ? sqrt(N2) : 0.0;
    if (windat->brunt) windat->brunt[c] = N;

    /* Get the mode 1 internal wave speed                            */
    /* c=NH/(n.pi), n=1,2,3... (Gill, 1982, p159) for long waves     */
    /* constant N.                                                   */
    iws = N * depth / pi;
    if (windat->int_wave) windat->int_wave[c] = iws;

    /* Vertical u1 gradient                                          */
    du1 = 0.5 * (windat->u1[c] - windat->u1[zm1] +
		 windat->u1[xp1] - windat->u1[window->zm1[xp1]]);

    /* Vertical u2 gradient                                          */
    du2 = 0.5 * (windat->u2[c] - windat->u2[zm1] +
		 windat->u2[yp1] - windat->u2[window->zm1[yp1]]);

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
    Roi = iws /
      fabs(0.5 * (wincon->u1c5[cs] + wincon->u1c5[window->xp1[cs]]));
    if (windat->rossby_in) windat->rossby_in[c] = Roi;

    /* External Rossby radius (m)                                    */
    Roe = sqrt(depth * wincon->g) /
      fabs(0.5 * (wincon->u1c5[cs] + wincon->u1c5[window->xp1[cs]]));
    if (windat->rossby_ex) windat->rossby_ex[cs] = Roe;

    /* Mechanical energy (Jm-3) (Kowalik and Murty (1993) p6, Eqn.   */
    /* 1.24.                                                         */
    me = windat->dens[c] * (wincon->g * windat->eta[cs] +
			    0.5 * (windat->u1[c] * windat->u1[c] +
				   windat->u2[c] * windat->u2[c] +
				   windat->w[c] * windat->w[c]));
    if (windat->energy) windat->energy[c] = me;
    ke = windat->dens[c] * 0.5 * (windat->u1[c] * windat->u1[c] +
				  windat->u2[c] * windat->u2[c] +
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
/* Computs the surface and bottom Ekman pumping                      */
/*-------------------------------------------------------------------*/
void ekman_pump_e1(geometry_t *window,   /* Window geometry          */
		   window_t *windat,     /* Window data              */
		   win_priv_t *wincon,   /* Window constants         */
		   double *taus,
		   double *taub
		   )
{
  int cc;
  double dp, dm;
  double *tau = wincon->d3;
  double val, botdz, Cdu1;
  double *u1au2 = wincon->w7;

  /* Initialize                                                      */
  memset(windat->sep, 0, window->sgsizS * sizeof(double));
  memset(windat->bep, 0, window->sgsizS * sizeof(double));

  /* Compute the bottom stress at e2 locations                       */
  memset(tau, 0 , window->sgsizS * sizeof(double));
  for (cc = 1; cc <= window->a2_e2; cc++) {
    int c = window->w2_e2[cc];
    int cs = window->m2d[c];
    int cb = window->bot_e2[cc];

    /* Bottom flux */
    Cdu1 = 0.5 * (wincon->Cd[cs] + wincon->Cd[window->ym1[cs]]);
    val = sqrt(windat->u2[cb] * windat->u2[cb] + u1au2[cb] * u1au2[cb]);
    val = Cdu1 * max(wincon->uf, val);
    botdz = max(wincon->hmin, windat->dzu2[cb] * wincon->Ds[cs]);
    /* Quadratic bottom friction, truncated to ensure stability */
    if (val > botdz / windat->dt)
      val = botdz / windat->dt;
    tau[cs] = -val * u1au2[cb];
  }

  /* Compute e1 contribution to curl of the Ekman pumping            */
  for (cc = 1; cc <= window->b2_t; cc++) {
    int c = window->w2_t[cc];
    int xp1 = window->xp1[c];
    int yp1 = window->yp1[c];
    int ym1 = window->xm1[c];
    int xpyp1 = window->yp1[xp1];
    int xpym1 = window->ym1[xp1];

    dp = 0.25 * (taus[c] + taus[xp1] + taus[yp1] + taus[xpyp1]);
    dm = 0.25 * (taus[c] + taus[xp1] + taus[ym1] + taus[xpym1]);
    windat->sep[c] -= 0.5 * (dp - dm) / (window->h2acell[c] * 
					 wincon->coriolis[c] *
					 windat->dens[c]);

    windat->bep[c] -= 0.5 * (tau[yp1] - tau[c]) / (window->h2acell[c] * 
					 wincon->coriolis[c] *
					 windat->dens[c]);
  }
}

void ekman_pump_e2(geometry_t *window,   /* Window geometry          */
		   window_t *windat,     /* Window data              */
		   win_priv_t *wincon,   /* Window constants         */
		   double *taus,
		   double *taub
		   )
{
  int cc;
  double dp, dm;
  double *tau = wincon->d3;
  double val, botdz, Cdu2;
  double *u2au1 = wincon->w7;

  /* Compute the bottom stress at e1 locations                       */
  memset(tau, 0 , window->sgsizS * sizeof(double));
  for (cc = 1; cc <= window->a2_e1; cc++) {
    int c = window->w2_e1[cc];
    int cs = window->m2d[c];
    int cb = window->bot_e1[cc];

    /* Bottom flux */
    Cdu2 = 0.5 * (wincon->Cd[cs] + wincon->Cd[window->xm1[cs]]);
    val = sqrt(windat->u1[cb] * windat->u1[cb] + u2au1[cb] * u2au1[cb]);
    val = Cdu2 * max(wincon->uf, val);
    botdz = max(wincon->hmin, windat->dzu1[cb] * wincon->Ds[cs]);
    /* Quadratic bottom friction, truncated to ensure stability */
    if (val > botdz / windat->dt)
      val = botdz / windat->dt;
    tau[cs] = -val * u2au1[cb];
  }

  /* Compute e2 contribution to curl of the Ekman pumping            */
  for (cc = 1; cc <= window->b2_t; cc++) {
    int c = window->w2_t[cc];
    int xp1 = window->xp1[c];
    int yp1 = window->yp1[c];
    int ym1 = window->xm1[c];
    int xpyp1 = window->yp1[xp1];
    int xpym1 = window->ym1[xp1];

    dp = 0.25 * (taus[c] + taus[xp1] + taus[yp1] + taus[xpyp1]);
    dm = 0.25 * (taus[c] + taus[xp1] + taus[ym1] + taus[xpym1]);
    windat->sep[c] += 0.5 * (dp - dm) / (window->h1acell[c] * 
					 wincon->coriolis[c] *
					 windat->dens[c]);

    windat->bep[c] += 0.5 * (tau[xp1] - tau[c]) / (window->h1acell[c] * 
					 wincon->coriolis[c] *
					 windat->dens[c]);
  }
}

/* END ekman_pump()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computs the depth of sound channels                               */
/*-------------------------------------------------------------------*/
void sound_channel(geometry_t *window,  /* Processing window         */
                    window_t *windat,   /* Window data structure     */
                    win_priv_t *wincon  /* Window geometry / constants */
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
  int dosmooth = 0;
  int donan = 0;

  tr = ss;
  if (windat->sound && windat->schan && windat->sonic) {
    memset(windat->schan, 0, window->sgsiz * sizeof(double));
    memset(windat->sonic, 0, window->sgsizS * sizeof(double));
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      if (donan)
	windat->schan[c] = NaN;
      else
	windat->schan[c] = 0.0;
    }
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      if (donan)
	windat->sonic[c] = NaN;
      else
	windat->sonic[c] = 0.0;
      c = window->bot_t[cc];
      windat->sound[window->zm1[c]] = windat->sound[c];
    }
    /* Smooth the sound profile vertically */
    k[0] = k[4] = 0.1;
    k[1] = k[3] = 0.2;
    k[2] = 0.4;
    /*k[0] = k[1] = k[2] = k[3] = k[4] = 0.2;*/
    memcpy(Fz, windat->sound, window->sgsiz * sizeof(double));
    if (dosmooth) {
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
	memcpy(Fz, ss, window->sgsiz * sizeof(double));
      }
      memcpy(windat->sound, ss, window->sgsiz * sizeof(double));
    } else
      memcpy(ss, windat->sound, window->sgsiz * sizeof(double));

    memset(Fz, 0, window->sgsiz * sizeof(double));
    /* Set the bottom boundary condition */
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->bot_t[cc];
      tr[window->zm1[c]] = tr[c];
    }

    /* Get the concentration values at the lower face */
    order4_do(Fz, tr, 1, wincon->vc, wincon->s2, window->zm1, window->zp1);
    /* Set the bottom boundary condition */
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->bot_t[cc];
      Fz[window->zm1[c]] = Fz[c];
    }

    /* Get the mean cell spacings and Courant number */
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      cs = window->m2d[c];
      zm1 = window->zm1[c];
      wincon->w5[c] = wincon->dz[zm1] / (wincon->dz[zm1] + wincon->dz[c]);
      wincon->w9[c] = 0.5 * (wincon->dz[c] + wincon->dz[zm1]);
      wincon->w8[c] = windat->dt / (wincon->w9[c] * wincon->Ds[cs]);
    }
    
    /* Get the derivative and depth */
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      zp1 = window->zp1[c];
      dSdz[c] = (Fz[zp1] - Fz[c]) / wincon->dz[c];
      depth[c] = window->cellz[c];
    }
    /* Get the surface depth, and set a no gradient condition at the surface */
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      cs = window->m2d[c];
      loc[cs] = cs;
      depth[c] = 0.5 * (windat->eta[cs] + window->gridz[c]);
      dSdz[c] = dSdz[window->zm1[c]];
      c = window->bot_t[cc];
      dSdz[window->zm1[c]] = dSdz[c];
    }
    /* Get the sonic depth (maximum in sound, or dSdz changes from negative  */
    /* to positive) and place in the top layer.                              */
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      zm1 = window->zm1[c];
      cs = window->m2d[c];
      if (dSdz[zm1]) {
        if (donan) {
          if (isnan(windat->sonic[cs]) && dSdz[c] < 0.0 && dSdz[zm1] > 0.0) {
             windat->sonic[cs] = depth[c] - dSdz[c] * (depth[zm1] - depth[c]) /
               (dSdz[zm1] - dSdz[c]);
          }
        } else {
          if (windat->sonic[cs] == 0 && dSdz[c] < 0.0 && dSdz[zm1] > 0.0) {
             windat->sonic[cs] = depth[c] - dSdz[c] * (depth[zm1] - depth[c]) /
               (dSdz[zm1] - dSdz[c]);
          }
        }
      }
    }
    /* Get the sound channel depth (minimum in sound, or dSdz changes from   */
    /* positive to negative) and place in layers below the top layer.        */
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
double buoyancy_frequency2(geometry_t *window, /* Window geometry       */
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
  double rho1 = 0.0, rho2 = 0.0;  /* Density of the top and bottom layers */
  double h1, h2;                /* Depths of the top and bottom layers */
  double top_m;                 /* Depth of top of pycnocline */
  double bot_m;                 /* Depth of bottom of pycnocline */
  int kts, ktb;                 /* k index of pycnocline top/bottom */
  double spd;                   /* Internal wave speed */
  double d1, d2;                /* Dummies */
  int cs = window->m2d[c];
  int cb = window->bot_t[cc];
  int lc, nc;

  mld(window, windat, wincon, cs, cb, &top_m, &bot_m, &kts, &ktb);
  if (ktb > cb)
    ktb = cb;
  if (kts < cs)
    kts = cs;

  /* Mean density in the surface layer */
  nc = 0;
  lc = cs;
  while (lc <= kts) {
    rho1 += windat->dens[lc];
    lc = window->zm1[lc];
    nc++;
  }
  rho1 /= (double)nc;
  /* Mean density in the bottom layer */
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
  double spd, d1;               /* Internal wave speed */
  double drhodz;                /* Vertical density gradient */
  double pi = 3.14159;          /* Value of pi */
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
double int_wave_spd_cont_m(master_t *master,  /* Master data structure */
                           int c, int cc  /* Sparse coordinate */
  )
{
  geometry_t *geom = master->geom;
  double spd, d1;               /* Internal wave speed */
  double drhodz;                /* Vertical density gradient */
  double pi = 3.14159;          /* Value of pi */
  double g = 9.81;              /* Acceleration due to gravity */
  double top, bot, dz;          /* Top, bottom, thickness of layers */
  int zm1, c3, cb;              /* Sparse coordinates */
  double r1, r2, h1, h2;        /* Mean densities, layer thicknesses */
  double thr = -0.01;           /* Density gradient threshold */
  int ml = 0;                   /* Flag for pycnocline */
  int mode = 0;                 /* Type of speed calculation */

  spd = 0.0;
  c3 = c;
  cb = geom->bot_t[cc];
  /* Get the internal wave speed using the 2 layer approximation */
  if (mode == 0) {
    h1 = r1 = r2 = h2 = 0.0;
    /* Find the pycnocline and mean densities in each layer */
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
void reset_cfl(master_t *master /* Model grid data structure */
  )
{
  int dt;
  double sf = 0.9;              /* Safety factor */
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
    for (cc = 1; cc <= geom->b2_e2; cc++) {
      c = geom->w2_e2[cc];
      master->u2kh[c] *= (master->grid_dt / olddt);
      master->u2vh[c] *= (master->grid_dt / olddt);
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
void steric(geometry_t *window, /* Processing window */
            window_t *windat,   /* Window data structure */
            win_priv_t *wincon  /* Window geometry / constants */
  )
{
  double lnm;                   /* Level of no motion (m).  */
  double pressu, pressl;        /* Pressure at upper and lower levels */
  double dnu, dnl;              /* Densities at upper and lower depths */
  double svu, svl;              /* Sigma-t at upper and lower depths */
  double asvu, asvl;            /* Anomaly of specific volume at depth 1&2
                                 */
  double geo;                   /* Geopotential anomaly at each depth
                                   level */
  int cc, c, c2, cs;            /* sparse coordinates/counters */
  int zm1;
  int cl;
  double *dzz;
  double d1;

  /* Set pointers and initialise */
  lnm = wincon->lnm;
  dzz = wincon->w9;
  memset(windat->steric, 0, window->sgsizS * sizeof(double));
  memset(dzz, 0, window->sgsiz * sizeof(double));

  /* Get the cell thickness between layers */
  /* Surface layer */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = wincon->s1[cc];
    c2 = window->m2d[c];
    dzz[c] = 0.5 * wincon->dz[c] * wincon->Ds[c2];
  }
  /* Subsurface layers */
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];

    c2 = window->m2d[c];
    zm1 = window->zm1[c];
    dzz[zm1] = 0.5 * (wincon->dz[c] + wincon->dz[zm1]) * wincon->Ds[c2];
  }

  /* Find the nearest level less than lnm.  */
  /* If lnm > bottom depth then set zero geopotential.  */
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

    /* Get anomaly specific volume for the surface */
    pressu =
      windat->patm[c2] +
      wincon->g * windat->dens[c] * dzz[c] * wincon->Ds[c2];
    dnu = windat->dens[c];
    eos2(35.0, 0.0, pressu, &svu, &d1);
    asvu = 1.0 / dnu - 1.0 / svu;
    geo = pressu * asvu;

    /* Get geopotential anomalies for depth levels 2 to cl */
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

    /* Get geopotential anomaly for depth level cl */
    zm1 = window->zm1[c];
    /* pressl = pressure at layer face cl */
    pressl +=
      wincon->g * windat->dens[c] * 0.5 * wincon->dz[c] * wincon->Ds[c2];
    /* pressl = pressure at depth lnm */
    pressl += wincon->g * windat->dens[zm1] * (lnm + window->gridz[c]);
    dnl = windat->dens[zm1];
    eos2(35.0, 0.0, pressl, &svl, &d1);
    asvl = 1.0 / dnl - 1.0 / svl;
    geo += (pressl - pressu) * (asvl + asvu) / 2.0;

    /* Convert the integrated geopotential anomaly to steric height */
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
void vorticity(geometry_t *window,  /* Processing window */
               window_t *windat,  /* Window data structure */
               win_priv_t *wincon /* Window geometry / constants */
  )
{
  int c, cc, c2;                /* Sparse coordinate/counter */
  double *rv;                   /* Relative vorticity */
  double *prv;                  /* d(rv)/dt */
  int xp1, yp1;                 /* Sparse location at x+1,x-1 */
  int xm1, ym1;                 /* Sparse location at y+1,y-1 */
  int xpym1, xmyp1;             /* Sparse location at (x+1,y-1) */
  double *u1au2, *u2au1;

  u1au2 = wincon->d1;
  u2au1 = wincon->d2;
  rv = wincon->d3;
  prv = windat->rv_drvdt;

  /* Get the velocity values at the correct faces */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = window->m2d[wincon->s1[cc]];
    ym1 = window->ym1[c];
    xp1 = window->xp1[c];
    xpym1 = window->xpym1[c];
    u1au2[c] = 0.25 * (windat->u1av[ym1] + windat->u1av[xpym1] +
                       windat->u1av[c] + windat->u1av[xp1]);

    xm1 = window->xm1[c];
    yp1 = window->yp1[c];
    xmyp1 = window->xmyp1[c];
    u2au1[c] = 0.25 * (windat->u2av[xm1] + windat->u2av[xmyp1] +
                       windat->u2av[c] + windat->u2av[yp1]);
  }
  /* Calculate relative vorticity */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = window->m2d[wincon->s1[cc]];
    xp1 = window->xp1[c];
    yp1 = window->yp1[c];
    rv[c] = (u2au1[xp1] - u2au1[c]) / window->h1acell[c] -
      (u1au2[yp1] - u1au2[c]) / window->h2acell[c];
  }
  /* Save the relative vorticity */
  if (wincon->vorticity & (RELATIVE | TENDENCY)) {
    /* Save the relative vorticity at the previous timestep */
    if (wincon->vorticity & TENDENCY)
      memcpy(prv, windat->rv, window->sgsizS * sizeof(double));
    memcpy(windat->rv, rv, window->sgsizS * sizeof(double));
  }
  /* Save the absolute vorticity */
  if (wincon->vorticity & ABSOLUTE) {
    memset(windat->av, 0, window->sgsizS * sizeof(double));
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = window->m2d[wincon->s1[cc]];
      windat->av[c] =
        rv[c] + 0.5 * (wincon->u1c5[c] + wincon->u1c5[window->xp1[c]]);
    }
  }

  /* Save the potential vorticity */
  if (wincon->vorticity & POTENTIAL) {
    memset(windat->pv, 0, window->sgsizS * sizeof(double));
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = window->m2d[wincon->s1[cc]];
      windat->pv[c] = (rv[c] + 0.5 * (wincon->u1c5[c] +
                                      wincon->u1c5[window->xp1[c]])) /
        (windat->eta[c] - window->botz[c]);
    }
  }

  /* Get the vorticity tendencies */
  if (wincon->vorticity & TENDENCY) {
    int zm1;
    double *adv = windat->rv_nonlin;
    double *beta = windat->rv_beta;
    double *stretch = windat->rv_strch;
    double *jebar = windat->rv_jebar;
    double *wsc = windat->rv_wsc;
    double *bsc = windat->rv_bsc;
    double *depth = wincon->d1;
    double *chi = wincon->d2;
    double *tsx = windat->wind1;
    double *tsy = windat->wind2;
    double *tbx = windat->rv_bsc; /* Saved in vdiff_u1() */
    double *tby = windat->rv_wsc; /* Saved in vdiff_u2() */
    double f;                   /* Coriolis parameter */
    double dHde1, dHde2;        /* Bottom gradients */
    double rho0 = 1024;         /* Standard reference density */
    double de1, de2;            /* 2xgrid spacings */
    double Dpe1, Dme1;          /* Face centered D at c, xp1 */
    double Dpe2, Dme2;          /* Face centered D at c, yp1 */

    memset(chi, 0, window->sgsizS * sizeof(double));
    for (cc = 1; cc <= window->n2_t; cc++) {
      c = window->w2_t[cc];
      c2 = window->m2d[c];
      zm1 = window->zm1[c];
      /* Depth */
      depth[c2] = windat->eta[c2] - window->botz[c2];
      /* Anomaly of potential energy */
      do {
        chi[c2] +=
          wincon->g * window->cellz[c] * windat->dens[c] * wincon->dz[c] /
          rho0;
        c = zm1;
        zm1 = window->zm1[c];
      } while (c != zm1);
    }

    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = window->m2d[wincon->s1[cc]];
      xp1 = window->xp1[c];
      xm1 = window->xm1[c];
      yp1 = window->yp1[c];
      ym1 = window->ym1[c];

      /* Cell centered Coriolis */
      f = 0.5 * (wincon->u1c5[c] + wincon->u1c5[xp1]);

      /* Bottom bathymetry gradients */
      de1 = window->h1au1[c] + window->h1au1[xp1];
      de2 = window->h2au2[c] + window->h2au2[yp1];
      dHde1 = (depth[xp1] - depth[xm1]) / de1;
      dHde2 = (depth[yp1] - depth[ym1]) / de2;

      /* 1. Rate of change of relative vorticity */
      prv[c] = (windat->rv[c] - prv[c]) / windat->dt;

      /* 2. Beta */
      beta[c] = (wincon->u2c5[yp1] - wincon->u2c5[c]) *
        (windat->u2av[c] + windat->u2av[yp1]) / window->h2acell[c];

      /* 3. Stretch */
      stretch[c] = (windat->rv[c] + f) * (windat->detadt[c] +
                                          windat->u1av[c] * dHde1 +
                                          windat->u2av[c] * dHde2) /
        depth[c];

      /* 4. JEBAR */
      dHde1 = (1.0 / depth[xp1] - 1.0 / depth[xm1]) / de1;
      dHde2 = (1.0 / depth[yp1] - 1.0 / depth[ym1]) / de2;
      jebar[c] = dHde2 * (chi[xp1] - chi[xm1]) / de1 -
        dHde1 * (chi[yp1] - chi[ym1]) / de2;

      /* 5. BSC */
      Dpe1 = 2.0 / (depth[c] + depth[xp1]);
      Dme1 = 2.0 / (depth[c] + depth[xm1]);
      Dpe2 = 2.0 / (depth[c] + depth[yp1]);
      Dme2 = 2.0 / (depth[c] + depth[ym1]);
      bsc[c] = ((tby[yp1] * Dpe2 - tby[c] * Dme2) / window->h2acell[c] -
                (tbx[xp1] * Dpe1 -
                 tbx[c] * Dme1) / window->h2acell[c]) / rho0;

      /* 6. WSC */
      wsc[c] = ((tsy[yp1] * Dpe2 - tsy[c] * Dme2) / window->h2acell[c] -
                (tsx[xp1] * Dpe1 -
                 tsx[c] * Dme1) / window->h2acell[c]) / rho0;

      /* 7. Nonlinear advection and diffusion */
      adv[c] = stretch[c] + jebar[c] + wsc[c] + bsc[c] - prv[c] - beta[c];
    }
  }
}

/* END vorticity()                                                   */
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
  double ns;                    /* Time increment                    */ 
  double ttime = 24;            /* Predominant tidal period (hours)  */
  double time;                  /* Cumulative time counter           */
  double *detadt = wincon->d4;  /* Mean forward eta difference       */

  if (wincon->means & NONE || !wincon->means_dt)
    return;
  if (mode & ALL) return;                  /* Legacy code            */

  /*-----------------------------------------------------------------*/
  /* Get the mean wind if required */
  if (mode & WIND && wincon->means & WIND) {
    int cc, c;
    double ns = windat->dtf;
    if (windat->w1m) {
      for (cc = 1; cc <= window->b2_e1; cc++) {
        c = window->w2_e1[cc];
	windat->w1m[c] = (windat->w1m[c] * windat->meanc[c] + 
			  windat->wind1[c] * ns) / (windat->meanc[c] + ns);
      }
    }
    if (windat->w2m) {
      for (cc = 1; cc <= window->b2_e2; cc++) {
        c = window->w2_e2[cc];
	windat->w2m[c] = (windat->w2m[c] * windat->meanc[c] + 
			  windat->wind2[c] * ns) / (windat->meanc[c] + ns);
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the mean T/S if required */
  if (mode & TS && wincon->means & (TS|MTRA3D|MTRA2D)) {
    int cc, c, cs;
    double ns = windat->dtf;
    if (windat->tempm && windat->saltm) {
      if (wincon->means & MMM) {
	for (cc = 1; cc <= window->b3_t; cc++) {
	  c = window->w3_t[cc];
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
	  windat->tempm[c] = (windat->tram[c] * windat->meanc[cs] + 
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
      memcpy(windat->odeta, detadt, window->sgsizS * sizeof(double));
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
/* Re-initializes means on the master                                */
/*-------------------------------------------------------------------*/
void reset_means_m(master_t *master)
{
  geometry_t *geom = master->geom;
  int mon;
  double itc;
  int n;

  if (master->means & NONE) return;

  if (master->t >= master->means_next) {
    master->means |= RESET;  /* Set the RESET flag for transfers     */
    itc =  master->meanc[1]; /* Save the time counter                */
    /* Re-initialize the time counter.                               */
    memset(master->meanc, 0, geom->sgsizS * sizeof(double));
    /* Set the time counter flag; this is set to 0 when meanc is     */
    /* re-initialized in reset_means_3d_m() or reset_means_2d_m()    */
    /* for SEASONAL or MONTHLY means below.                          */
    master->meanc[0] = 1.0;  
    /* Re-initialize 3D mean tracers                                 */
    for (n = 0; n < master->ntm_3d; n++)
      reset_means3d_m(master, master->tm_3d[n]);
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
    }
    else
      master->means_next = master->t + master->means_dt;
  }
}

/* END reset_means_m()                                               */
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
    if (read_mean_3d(master, t, sf, vname, master->tr_wc[trm], 
		     dumpdata->tr_wc[trd])) {
      if (master->meanc[0]) {
	for (cc = 1; cc <= geom->enonS; cc++)
	  master->meanc[cc] = master->meancs[mon];
	master->meanc[0] = 0.0;
      }
    }
  } else {
    memset(master->tr_wc[trm], 0, geom->sgsiz * sizeof(double));
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
  size_t nce1, nce2, nz, ndumps;
  int ti;
  int n, m, nn;
  size_t start[4];
  size_t count[4];
  int fid;
  double tvals;
  char timeunits[MAXSTRLEN];
  char incname[MAXSTRLEN];
  int ncerr;
  
  strcpy(incname, "SEASONAL");
  if (sf == MONTHLY)
    strcpy(incname, "MONTHLY");

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
	  nc_inq_dimlen(fid, ncw_dim_id(fid, "k_centre"), &nz);
	  nc_inq_dimlen(fid, ncw_dim_id(fid, "j_centre"), &nce2);
	  nc_inq_dimlen(fid, ncw_dim_id(fid, "i_centre"), &nce1);
	  if (nce1 != geom->nce1 || nce2 != geom->nce2 || nz != geom->nz)
	    continue;
	  nc_inq_dimlen(fid, ncw_dim_id(fid, "record"), &ndumps);
	  /* Get time units                                          */
	  memset(timeunits, 0, MAXSTRLEN);
	  nc_get_att_text(fid, ncw_var_id(fid, "t"), "units", timeunits);
	  /* Find the dump number corresponding to t                 */
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
	  /* Dump file doesn't contain the time; zero and return     */
	  if (ti < 0) {
	    memset(v, 0, geom->sgsiz * sizeof(double));
	    return(0);
	  }
	  /* Read in the data and map to the sparse variable         */
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
	}
      }
    }
  }
  memset(v, 0, geom->sgsiz * sizeof(double));
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
    if (read_mean_2d(master, t, sf, vname, master->tr_wcS[trm], 
		     dumpdata->tr_wcS[trd])) {
      if (master->meanc[0]) {
	for (cc = 1; cc <= geom->enonS; cc++)
	  master->meanc[cc] = master->meancs[mon];
	master->meanc[0] = 0.0;
      }
    }
  } else
    memset(master->tr_wcS[trm], 0, geom->sgsizS * sizeof(double));
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
  size_t nce1, nce2, ndumps;
  int ti;
  int n, m, nn;
  size_t start[4];
  size_t count[4];
  int fid;
  double tvals;
  char timeunits[MAXSTRLEN];
  char incname[MAXSTRLEN];
  int ncerr;

  strcpy(incname, "SEASONAL");
  if (sf == MONTHLY)
    strcpy(incname, "MONTHLY");

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
	  nc_inq_dimlen(fid, ncw_dim_id(fid, "j_centre"), &nce2);
	  nc_inq_dimlen(fid, ncw_dim_id(fid, "i_centre"), &nce1);
	  if (nce1 != geom->nce1 || nce2 != geom->nce2)
	    continue;
	  nc_inq_dimlen(fid, ncw_dim_id(fid, "record"), &ndumps);
	  /* Get time units                                          */
	  memset(timeunits, 0, MAXSTRLEN);
	  nc_get_att_text(fid, ncw_var_id(fid, "t"), "units", timeunits);
	  /* Find the dump number corresponding to t                 */
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
	  /* Dump file doesn't contain the time; zero and return     */
	  if (ti < 0) {
	    memset(v, 0, geom->sgsizS * sizeof(double));
	    return(0);
	  }
	  /* Read in the data and map to the sparse variable         */
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
	}
      }
    }
  }
  memset(v, 0, geom->sgsizS * sizeof(double));
  return(0);
}

/* END read_mean_2d_m()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set up the totals output timeseries file.              */
/*-------------------------------------------------------------------*/
void init_totals(master_t *master,  /* Master data structure */
		 geometry_t **window, /* Processing window */
		 win_priv_t **wincon  /* Window geometry / constants */
  )
{
  geometry_t *geom = master->geom;
  int n, m, nn, trn;

  /*---------------------------------------------------------------*/
  /* Initialize the output file */
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
    for (nn = 2; nn < master->ntot; nn++) { /* n=0:salt, n=1:temp */
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

  /*---------------------------------------------------------------*/
  /* Get the tracers ids for all totals                            */
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
    if (open->ocodex & (L_EDGE | R_EDGE)) {
      vel = windat->u1;
      hat = window->h2au1;
      cv = open->obc_e1;
    } else {
      vel = windat->u2;
      hat = window->h1au2;
      cv = open->obc_e2;
    }
    if (open->ocodex & R_EDGE || open->ocodey & F_EDGE) sc = -1.0;
    windat->vf[n] = 0.0;
    for (cc = 1; cc <= open->no3_t; cc++) {
      c = open->obc_t[cc];
      v = vel[cv[cc]];	
      cs = window->m2d[cv[cc]];
      windat->vf[n] += (sc * v * hat[cs] * wincon->dz[c]);
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
/* Returns the maximum velocity threshold for constraining velocity  */
/*-------------------------------------------------------------------*/
double con_velmax(double vel, double velmax)
{
  return((vel > 0.0) ? velmax : -velmax);
}

/* END con_velmax()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Returns the median velocity in the vicinity of the cell for       */
/* constraining velocity.                                            */
/*-------------------------------------------------------------------*/
double con_median(geometry_t *window,     /* Window geometry         */
		  double *vel,            /* Velocity array          */
		  double velmax,          /* Velocity threshold      */
		  int c)                  /* Sparse coodinate        */
{
  int i, j;
  double v[9];
  v[0] = vel[window->xm1[c]];
  v[1] = vel[window->xp1[c]];
  v[2] = vel[window->ym1[c]];
  v[3] = vel[window->yp1[c]];
  v[4] = vel[window->xm1[window->ym1[c]]];
  v[5] = vel[window->xp1[window->ym1[c]]];
  v[6] = vel[window->xm1[window->yp1[c]]];
  v[7] = vel[window->xp1[window->yp1[c]]];
  v[8] = vel[c];

  for (i = 0; i < 8; i++)
    for (j = 8; i < j; --j)
      order(&v[j-1], &v[j]);
  if (v[4] > velmax) v[4] = velmax;
  if (v[4] < -velmax) v[4] = -velmax;
  return(v[4]);
}

double con_med(geometry_t *window,     /* Window geometry         */
	       double *vel,            /* Velocity array          */
	       int c)                  /* Sparse coodinate        */
{
  int i, j;
  double v[9];
  v[0] = vel[window->xm1[c]];
  v[1] = vel[window->xp1[c]];
  v[2] = vel[window->ym1[c]];
  v[3] = vel[window->yp1[c]];
  v[4] = vel[window->xm1[window->ym1[c]]];
  v[5] = vel[window->xp1[window->ym1[c]]];
  v[6] = vel[window->xm1[window->yp1[c]]];
  v[7] = vel[window->xp1[window->yp1[c]]];
  v[8] = vel[c];

  for (i = 0; i < 8; i++)
    for (j = 8; i < j; --j)
      order(&v[j-1], &v[j]);
  return(v[4]);
}

void order (double *a, double *b)
{
  double val;
  if (*a > *b) {
    val = *a;
    *a = *b;
    *b = val;
  }
}

/* END con_median()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Returns the mean velocity in the vicinity of the cell for         */
/* constraining velocity.                                            */
/*-------------------------------------------------------------------*/
double con_mean(geometry_t *window,     /* Window geometry         */
		  double *vel,            /* Velocity array          */
		  double velmax,          /* Velocity threshold      */
		  int c)                  /* Sparse coodinate        */
{
  int i;
  double v[8], m, n, ret;
  v[0] = vel[window->xm1[c]];
  v[1] = vel[window->xp1[c]];
  v[2] = vel[window->ym1[c]];
  v[3] = vel[window->yp1[c]];
  v[4] = vel[window->xm1[window->ym1[c]]];
  v[5] = vel[window->xp1[window->ym1[c]]];
  v[6] = vel[window->xm1[window->yp1[c]]];
  v[7] = vel[window->xp1[window->yp1[c]]];

  m = n = 0.0;
  for (i = 0; i < 8; i++) {
    if (v[i] > -velmax && v[i] < velmax) {
      m += v[i];
      n += 1.0;
    }
  }
  ret = (n > 0.0) ? m / n : con_velmax(vel[c], velmax);
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
  master_t *master;
  char pfile[MAXSTRLEN];
  char files[MAXNUMTSFILES][MAXSTRLEN];
  char *arg[MAXSTRLEN * MAXNUMARGS];
  char outfile[MAXSTRLEN];
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
  int i, j, k, kk, c, v, n, m, mm;
  int vids[Ns];
  int ndims;
  int *sdump;
  int *edump;
  int nargs;
  int sparse, spf;
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
  double **botz, *bots;
  double d1, d2, dtf, dt = 0.0;
  double mean, nm;
  double nn, ntot, clck, mclck;
  double elapsed = 0.0, start_time;
  long t;
  struct timeval tm1;
  int *layers, nl;

  params = params_alloc();
  master = master_alloc();
  master->tsfile_caching = 0;
  master->lyear = 0;
  start_time = time(NULL);
  time(&t);

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
  nvars = parseline(vlist, arg, MAXNUMARGS);
  vars = (char **)malloc(sizeof(char *)*nvars);
  j = 0;
  for (i=0; i<nvars; ++i) {
    if (hd_ts_multifile_get_index(ntsfiles, tsfiles,
				  filenames, arg[i], &n) <= 0)
      hd_warn("Can't find variable %s in files.\n", arg[i]);
    else {
      vars[j] = arg[i];
      j++;
    }
  }
  nvars = j;
  hd_warn("percentiles: variables %s read OK.\n", vlist);

  /*-----------------------------------------------------------------*/
  /* Get the dump numbers for each input file                        */
  sparse = 0;
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
    for (m=0; m<df->nd;m++) {
      if (strcmp(df->dimensions[m].name, "ns2") == 0) {
	params->ns2 = df->dimensions[m].size;
	sparse = spf = 1;
      }
      if (strcmp(df->dimensions[m].name, "ns3") == 0) 
	params->ns3 = df->dimensions[m].size;
      if (strcmp(df->dimensions[m].name, "i_centre") == 0) 
	params->nce1 = df->dimensions[m].size;
      if (strcmp(df->dimensions[m].name, "j_centre") == 0) 
	params->nce2 = df->dimensions[m].size;
      if (strcmp(df->dimensions[m].name, "k_centre") == 0) 
	params->nz = df->dimensions[m].size;
    }
    if (sparse && !spf) hd_quit("Can't mix sparse and Cartesian input files.\n");
    if (sparse) {
      hd_warn("Input files are sparse with size %d\n", params->ns3);
      hd_warn("Output files are Cartesian with size %dx%dx%d\n", params->nce1,
	      params->nce2, params->nz);
    }
    else
      hd_warn("Input/output files are Cartesian with size %dx%dx%d\n", params->nce1,
	      params->nce2, params->nz);

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
      if (sparse) {
	s2i = i_alloc_1d(params->ns3);
	s2j = i_alloc_1d(params->ns3);
	s2k = i_alloc_1d(params->ns3);
	outfid = create_output(params, df, outfile, vars, nvars, 
				 s2i, s2j, s2k);
      } else {
	outfid = create_output(params, df, outfile, vars, nvars, 
				 NULL, NULL, NULL);
      }
      if (sparse) {
	bots = d_alloc_1d(params->ns2);
	botz = d_alloc_2d(params->nce1, params->nce2);
	start[0] = 0;
	start[1] = 0;
	start[2] = 0;
	start[3] = 0;
	count[0] = params->ns2;
	count[1] = 0;
	count[2] = 0;
	count[3] = 0;
	nc_get_vara_double(df->ncid, ncw_var_id(df->ncid, "botz"), 
			   start, count, bots);
      } else {
	botz = d_alloc_2d(params->nce1, params->nce2);
	start[0] = 0;
	start[1] = 0;
	start[2] = 0;
	start[3] = 0;
	count[0] = params->nce2;
	count[1] = params->nce1;
	count[2] = 0;
	count[3] = 0;
	nc_get_vara_double(df->ncid, ncw_var_id(df->ncid, "botz"), 
			   start, count, botz[0]);
      }
      first = 0;
    }
  }
  data = d_alloc_1d(dsize);

  /*-----------------------------------------------------------------*/
  /* Get the sparse maps and set up botz for sparse input files      */
  if (sparse) {
    prm_read_int(fp, "LAYERFACES", &k);
    if (k-1 != params->nz)
      hd_quit("Sparse input files must have 'k_centre = %d'\n",k-1);
    map = i_alloc_3d(params->nce1, params->nce2, params->nz);
    for (j=0; j<params->nce2; ++j) {
      for (i=0; i<params->nce1; ++i) {
				botz[j][i] = LANDCELL;
				for (k=0; k<params->nz; ++k)
					map[k][j][i] = 0;
      }
    }
    for (c=0; c<params->ns2; ++c) {
      i = s2i[c];
      j = s2j[c];      
      botz[j][i] = bots[c];
    }
    d_free_1d(bots);
    for (c=0; c<params->ns3; ++c) {
      i = s2i[c];
      j = s2j[c];
      k = s2k[c];
      if (botz[j][i] == LANDCELL || botz[j][i] == NOTVALID ||
					isnan(botz[j][i]))
				map[k][j][i] = 0;
      else
				map[k][j][i] = c;
    }
  }
  vid = ncw_var_id(outfid, "botz");
  count[0] = params->nce2;
  count[1] = params->nce1;
  count[2] = 0;
  count[3] = 0;
  nc_put_vara_double(outfid,vid,start,count,botz[0]);
  hd_warn("percentiles: Processing data.\n");
 
  /*-----------------------------------------------------------------*/
  /* Read the layers to process                                      */
  if (prm_read_char(fp, "P_LAYERS", buf)) {
    nl = parseline(buf, arg, MAXNUMARGS);
    layers = i_alloc_1d(nl);
    if (nl == 1 && strcmp(buf, "surf") == 0) {
      layers[0] = params->nz - 1;
    } else {
      for (i=0; i<nl; ++i)
	layers[i] = atoi(arg[i]);
    }
  } else {
    nl = params->nz;
    layers = i_alloc_1d(nl);
    for (i=0; i<nl; ++i)
      layers[i] = i;
  }

  /*-----------------------------------------------------------------*/
  /* Process the data                                                */ 
  nn = mclck = 0.0;
  ntot = (double)(nl * params->nce1 * params->nce2 * nvars);
  for (v=0; v<nvars; ++v) {
    /* Build variable percentile netCDF file id's */
    for (n=0; n<Ns; ++n) {
      sprintf(buf, "%s_%3.3d", vars[v], (int)(percentiles[n]*100));
      vids[n] = ncw_var_id(outfid, buf);
    }
    mean = nm = 0.0;
    for (kk=0; kk<nl; ++kk) {
      k = layers[kk];
      hd_warn("Percentiles: processing %s layer %d\n", vars[v], k);
      for (j=0; j<params->nce2; ++j) {
	for (i=0; i<params->nce1; ++i) {
	  nn += 1.0;
	  gettimeofday(&tm1, NULL);
	  clck = tm1.tv_sec + tm1.tv_usec * 1e-6;
	  c = (sparse) ? map[k][j][i] : 1;
	  /* Non-valid cells get a NaN                             */
	  if (c && (botz[j][i] == LANDCELL || botz[j][i] == NOTVALID ||
		    isnan(botz[j][i]))) {
	    start[0] = k;
	    start[1] = j;
	    start[2] = i;
	    count[0] = 1;
	    count[1] = 1;
	    count[2] = 1;
	    for (n=0; n<Ns; ++n) {
	      results[n] = NaN;
	      nc_put_vara_double(outfid, vids[n],
				 start, count, &results[n]);
	    }
	    continue;
	  }
          /*---------------------------------------------------------*/
	  /* Read data                                               */
	  if (sparse) {   /* Sparse file input                       */
	    ad = 0;
	    for (m = 0; m < ntsfiles; m++) {
	      df = tsfiles[m]->df;
	      if (mnc) { /* multi-netcdf-version */
		df_multi_t *fd = (df_multi_t *)df->private_data;
		for (mm=0; mm<fd->nfiles; ++mm) {
		  vid = ncw_var_id(mncid[mm], vars[v]);
		  start[0] = sdump[mm];
		  start[1] = c;
		  count[0] = edump[mm]-sdump[mm]+1;
		  count[1] = 1;
		  nc_get_vara_double(mncid[mm], vid, start, count, &data[ad]);
		  ad += count[0];
		}
	      } else { 
		vid = ncw_var_id(df->ncid, vars[v]);
		start[0] = sdump[m];
		start[1] = c;
		count[0] = edump[m]-sdump[m]+1;
		count[1] = 1;
		nc_get_vara_double(df->ncid, vid, start, count, &data[ad]);
		ad += count[0];
	      }
	    }
	  } else {   /* Cartesian file input                         */
	    ad = 0;
	    for (m = 0; m < ntsfiles; m++) {
	      df = tsfiles[m]->df;
	      if (mnc) { /* multi-netcdf-version */
		df_multi_t *fd = (df_multi_t *)df->private_data;
		for (mm=0; mm<fd->nfiles; ++mm) {
		  vid = ncw_var_id(mncid[mm], vars[v]);
		  nc_inq_varndims(mncid[mm], vid, &ndims);
		  start[0] = sdump[mm];
		  count[0] = edump[mm]-sdump[mm]+1;
		  if (ndims == 3) {
		    start[1] = j;
		    start[2] = i;
		    start[3] = 0;
		    count[1] = 1;
		    count[2] = 1;
		    count[3] = 0;
		    nc_get_vara_double(mncid[mm], vid, start, count, &data[ad]);
		    ad += count[0];
		  } else {
		    start[1] = k;
		    start[2] = j;
		    start[3] = i;
		    count[1] = 1;
		    count[2] = 1;
		    count[3] = 1;
		    nc_get_vara_double(mncid[mm], vid, start, count, &data[ad]);
		    if (!isnan(data[ad]) && data[ad] > 0.0) {
		      mean += data[ad];
		      nm += 1.0;
		    }
		    ad += count[0];
		  }
		}
	      } else { 
		vid = ncw_var_id(df->ncid, vars[v]);
		nc_inq_varndims(df->ncid, vid, &ndims);
		start[0] = sdump[m];
		count[0] = edump[m]-sdump[m]+1;
		if (ndims == 3) {
		  start[1] = j;
		  start[2] = i;
		  count[1] = 1;
		  count[2] = 1;
		  nc_get_vara_double(df->ncid, vid, start, count, &data[ad]);
		  ad += count[0];
		} else {
		  start[1] = k;
		  start[2] = j;
		  start[3] = i;
		  count[1] = 1;
		  count[2] = 1;
		  count[3] = 1;
		  nc_get_vara_double(df->ncid, vid, start, count, &data[ad]);
		  ad += count[0];
		}
	      }
	    }
	  }
	  /* Resort the data into inc steps */
	  asize = 0;
	  for (m = 0; m < dsize; m++) {
	    if (m%inc == 0) {
	      data[asize] = data[m];
	      asize++;
	    }
	  }
	  /* Compute statistics */
	  stats(data, asize, results);
	  /* Write data */
	  start[0] = k;
	  start[1] = j;
	  start[2] = i;
	  count[0] = 1;
	  count[1] = 1;
	  count[2] = 1;
	  for (n=0; n<Ns; ++n) {
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
	    fprintf(dfp, "(%d %d %d) of (%d %d %d) : %s(%d/%d) : %% complete = %5.2f\n", i, j, k, 
		    params->nce1, params->nce2, params->nz, vars[v], v+1, nvars, 100.0 * nn / ntot);
	    fclose(dfp);
	  }
	}
      }
      hd_warn("  Domain mean of %s = %f\n", vars[v], mean / nm);
    }
  }
  
  nc_close(outfid);
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
  if (sparse) {
    i_free_1d(s2i);
    i_free_1d(s2j);
    i_free_1d(s2k);
  }
  d_free_2d(botz);
  d_free_1d(data);
  free(filenames);
  free(tsfiles);
  free(params);
  master_free(master);
  exit(0);
}


static int create_output(parameters_t *params,   datafile_t *df,
			 char *filename, char** vars, int nvars, 
			 int *s2i, int *s2j, int *s2k)
{
  int cdfid = -1;
  int ingeog;
  int recdimid, e1id, e2id, e3id, xid, yid, zid, tid, vid, ns3id, i, j;
  int s2iid, s2jid, s2kid, kgridid;
  char buf[MAXSTRLEN];
  int dims[4];
  size_t start[4];
  size_t count[4];
  double **cell, *cellz, *gridz;
  int has_proj = (strlen(projection) > 0);
  int is_geog = has_proj && (strcasecmp(projection, GEOGRAPHIC_TAG) == 0);

  /* Now that all the data has been parameters have been read.
   * Create the netCDF file. */
  if (nc_create(filename,NC_NOCLOBBER, &cdfid) != NC_NOERR)
    hd_quit("Couldn't create output dump file %s\n", filename);

  /* define dimensions */
  nc_def_dim(cdfid,"time",NC_UNLIMITED,&recdimid);
  nc_def_dim(cdfid,"i", params->nce1,&e1id);
  nc_def_dim(cdfid,"j", params->nce2,&e2id);
  nc_def_dim(cdfid,"k", params->nz,&e3id);
  /*nc_def_dim(cdfid,"k_grid", params->nz+1,&kgridid);*/
  
  dims[0] = e2id;
  dims[1] = e1id;

  /*  if (is_geog) {
    nc_def_var(cdfid,"longitude",NC_DOUBLE,2,dims,&xid);
    } else {*/
    nc_def_var(cdfid,"x",NC_DOUBLE,2,dims,&xid);
    if (is_geog) {
      write_text_att(cdfid,xid,"long_name", "Longitude");
      write_text_att(cdfid,xid,"coordinate_type", "longitude");
      write_text_att(cdfid,xid,"units", "degrees_east");
    } else {
      write_text_att(cdfid,xid,"long_name", "X coordinate");
      write_text_att(cdfid,xid,"coordinate_type", "X");
      write_text_att(cdfid,xid,"units", "metres");
    }
  
    /*if (is_geog) {
    nc_def_var(cdfid,"latitude",NC_DOUBLE,2,dims,&yid);
    } else {*/
    nc_def_var(cdfid,"y",NC_DOUBLE,2,dims,&yid);
    if (is_geog) {
      write_text_att(cdfid,yid,"long_name", "Latitude");
      write_text_att(cdfid,yid,"coordinate_type", "latitude");
      write_text_att(cdfid,yid,"units", "degrees_north");
    } else {
      write_text_att(cdfid,yid,"long_name", "Y coordinate");
      write_text_att(cdfid,yid,"coordinate_type", "Y");
      write_text_att(cdfid,yid,"units", "metres");
    }
  
  dims[0] = e3id;
  nc_def_var(cdfid,"zc",NC_DOUBLE,1,dims,&zid);
  write_text_att(cdfid,zid,"long_name", "Z coordinate at grid layer centre");
  write_text_att(cdfid,zid,"units", "metres");
  write_text_att(cdfid,zid,"coordinate_type", "Z");
  /*
  dims[0] = kgridid;
  nc_def_var(cdfid,"z_grid",NC_DOUBLE,1,dims,&zid);
  write_text_att(cdfid,zid,"long_name", "Z coordinate at grid layer faces");
  write_text_att(cdfid,zid,"units", "metres");
  write_text_att(cdfid,zid,"coordinate_type", "Z");
  */
  /* time dependent variables */
  dims[0] = recdimid;
  nc_def_var(cdfid,"t",NC_DOUBLE,1,dims,&tid);
  write_text_att(cdfid,tid,"units", params->timeunit);
  write_text_att(cdfid,tid,"long_name", "");
  write_text_att(cdfid,tid,"coordinate_type", "time");
  
  dims[0] = e2id;
  dims[1] = e1id;
  nc_def_var(cdfid,"botz",NC_DOUBLE,2,dims,&vid);
  write_text_att(cdfid,vid,"long_name", "Z coordinate at sea-bed at cell centre");
  write_text_att(cdfid,vid,"units", "metres");
  /*if (is_geog)
    write_text_att(cdfid,vid,"coordinate", "longitude, latitude");
    else*/
    write_text_att(cdfid,vid,"coordinate", "x, y");

  for (i=0; i<nvars; ++i) {
    for (j=0; j<Ns; ++j) {
      sprintf(buf, "%s_%3.3d", vars[i], (int)(percentiles[j]*100));
      dims[0] = e3id;
      dims[1] = e2id;
      dims[2] = e1id;
      nc_def_var(cdfid,buf,NC_DOUBLE,3,dims,&vid);
      write_text_att(cdfid,vid,"units", "");
      sprintf(buf, "%s %g percentile", vars[i], percentiles[j]);
      write_text_att(cdfid,vid,"long_name", buf);
      /*if (is_geog)
	write_text_att(cdfid,vid,"coordinate", "longitude, latitude, z");
	else*/
	write_text_att(cdfid,vid,"coordinates", "x, y, z");
    }
  }
  
  nc_enddef(cdfid);

  /* Read position information */
  cellz = d_alloc_1d(params->nz);
  gridz = d_alloc_1d(params->nz+1);
  cell = d_alloc_2d(params->nce1, params->nce2);
  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;
  count[0] = params->nce2;
  count[1] = params->nce1;
  count[2] = 0;
  count[3] = 0;

  nc_get_vara_double(df->ncid, ncw_var_id(df->ncid, "x_centre"), 
		     start, count, cell[0]);
  nc_put_vara_double(cdfid,xid,start,count,cell[0]);
  nc_get_vara_double(df->ncid, ncw_var_id(df->ncid, "y_centre"), 
		     start, count, cell[0]);
  nc_put_vara_double(cdfid,yid,start,count,cell[0]);

  count[0] = params->nz;
  count[1] = 0;
  nc_get_vara_double(df->ncid, ncw_var_id(df->ncid, "z_centre"), 
		     start, count, cellz);
  nc_put_vara_double(cdfid,zid,start,count,cellz);
  /*
  count[0] = params->nz+1;
  count[1] = 0;
  nc_get_vara_double(df->ncid, ncw_var_id(df->ncid, "z_grid"), 
		     start, count, gridz);
  nc_put_vara_double(cdfid,zid,start,count,gridz);
  */
  count[0] = params->ns3;
  if (s2i != NULL)
    nc_get_vara_int(df->ncid, ncw_var_id(df->ncid, "s2i"), 
		       start, count, s2i);
  if (s2j != NULL)
    nc_get_vara_int(df->ncid, ncw_var_id(df->ncid, "s2j"), 
		       start, count, s2j);
  if (s2k != NULL)
    nc_get_vara_int(df->ncid, ncw_var_id(df->ncid, "s2k"), 
		       start, count, s2k);
  
  /* Output the start time */
  start[0] = 0;
  count[0] = 1;
  nc_put_vara_double(cdfid,tid,start,count,&params->t);

  d_free_1d(cellz);
  d_free_1d(gridz);
  d_free_2d(cell);

  return cdfid;
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
/* Nore: 2D fluxes in zoomed windows boundaries should be u1flux_z   */
/*-------------------------------------------------------------------*/
void check_flux(geometry_t *window, /* Window geometry               */
                window_t *windat,   /* Window data                   */
                win_priv_t *wincon, /* Window constants              */
                int dw, int dc)
{
  int cc, c, cs;
  double d1, d2, d3, d4;
  if (window->wn == dw) {
    d1 = d2 = d3 = d4 = 0;
    for (cc = 1; cc <= window->enon; cc++) {
      c = cc;
      cs = window->m2d[c];
      if (cs == dc) {
        d1 +=
          windat->u1[c] * windat->dzu1[c] * wincon->mdx[cs] * windat->dt *
          window->h2au1[cs];
        d2 += windat->u1flux3d[c];
        d3 +=
          windat->u2[c] * windat->dzu2[c] * wincon->mdy[cs] * windat->dt *
          window->h1au2[cs];
        d4 += windat->u2flux3d[c];
      }
    }
    emstag(LDEBUG,"hd:vel3d:check_flux", "check_flux %e %e %e %e\n",
            windat->u1flux[window->m2d[dc]] - d1,
            windat->u1flux[window->m2d[dc]] - d2 * windat->dt,
            windat->u2flux[window->m2d[dc]] - d3,
            windat->u2flux[window->m2d[dc]] - d4 * windat->dt);
    printf("check_flux (c) %e %e %e %e\n",
            windat->u1flux[window->m2d[dc]] - d1,
            windat->u1flux[window->m2d[dc]] - d2 * windat->dt,
            windat->u2flux[window->m2d[dc]] - d3,
            windat->u2flux[window->m2d[dc]] - d4 * windat->dt);

    d1 = d2 = d3 = d4 = 0;
    for (cc = 1; cc <= window->enon; cc++) {
      c = cc;
      cs = window->m2d[c];
      if (cs == window->xp1[dc]) {
        d1 +=
          windat->u1[c] * windat->dzu1[c] * wincon->mdx[cs] * windat->dt *
          window->h2au1[cs];
        d2 += windat->u1flux3d[c];
	printf("check %d %d %f %f %f %f %f %10.12f\n",window->s2k[c],c,d2,windat->u1flux3d[c],windat->u1flux[window->m2d[c]]/windat->dt,windat->u1flux[window->m2d[c]],windat->dzu1[c],windat->eta[cs]);
      }
      if (cs == window->yp1[dc]) {
        d3 +=
          windat->u2[c] * windat->dzu2[c] * wincon->mdy[cs] * windat->dt *
          window->h1au2[cs];
        d4 += windat->u2flux3d[c];
      }
    }
    printf("check_flux (xp1, yp1) %e %e %e %e\n",
            windat->u1flux[window->m2d[window->xp1[dc]]] - d1,
            windat->u1flux[window->m2d[window->xp1[dc]]] - d2 * windat->dt,
            windat->u2flux[window->m2d[window->yp1[dc]]] - d3,
            windat->u2flux[window->m2d[window->yp1[dc]]] - d4 * windat->dt);

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
  xytoij_tree_t *xyij_tree = master->xyij_tree;
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
  int n, i, j, c, cc;
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

  /* Read the vertical layer structure */
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

  /* Get the model grid bounds */
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
  /* Truncate the bathy grid */
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
  memset(master->bathy_range_min, 0.0, geom->sgsizS);
  memset(master->bathy_range_max, 0.0, geom->sgsizS);
  for (j = js; j < je; j += jinc) {
    for (i = is; i < ie; i += iinc) {
      if (grid_xytoij(xyij_tree, lon[i], lat[j], &mi, &mj) == 0) continue;
      c = geom->map[nz][mj][mi];
      if (c <= 0 || c > geom->ewetS) continue;
      d1 = fabs(bathy[j][i]) - fabs(geom->botz[c]);
      /* Get the difference between model bathy and the deepest source bathy */
      if (d1 < 0 && fabs(d1) > master->bathy_range_min[c]) master->bathy_range_min[c] = fabs(d1);
      /* Get the difference between model bathy and the deepest source bathy */
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
	if (grid_xytoij(xyij_tree, lon[i], lat[j], &mi, &mj) == 0) continue;
	if (mi < 1 || mi > geom->nce1-1 || mj < 1 || mj > geom->nce2-1) continue;
	maxgrad = 0.0;
	mingrad = HUGE;

	/* Get the model grid gradients                                */
	c = geom->map[nz][mj][mi];
	if (c <= 0 || c > geom->ewetS) continue;
	g1 = fabs((geom->botz[c] - geom->botz[geom->xm1[c]]) / geom->h1au1[c]);
	g2 = fabs((geom->botz[c] - geom->botz[geom->ym1[c]]) / geom->h2au2[c]);
	
	/* Get the bathy gradient across i grid faces                  */
	if (grid_xytoij(xyij_tree, lon[i-1], lat[j], &ni, &nj) == 0) continue;
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
	if (grid_xytoij(xyij_tree, lon[i], lat[j-1], &ni, &nj) == 0) continue;
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
	master->bathy_grad_max[c] = max(fabs(fabs(maxgrad) - max(g1, g2)),
					master->bathy_grad_max[c]);
	master->bathy_grad_min[c] = min(fabs(fabs(mingrad) - min(g1, g2)),
					master->bathy_grad_min[c]);
      }
    }
    for (c = 0; c < geom->sgnumS; c++) {
      if (master->bathy_grad_min[c] == HUGE) master->bathy_grad_min[c] = 0.0;
    }
  }
  dograd = 1;
  if (dograd) {
    memset(master->bathy_grad_min, 0.0, geom->sgsizS);
    memset(master->bathy_grad_max, 0.0, geom->sgsizS);
    /* Get the model grid gradient                                   */
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = geom->w2_t[cc];
      g1 = fabs((geom->botz[c] - geom->botz[geom->xm1[c]]) / geom->h1au1[c]);
      g2 = fabs((geom->botz[c] - geom->botz[geom->ym1[c]]) / geom->h2au2[c]);
      master->bathy_grad_min[c] = sqrt(g1 * g1 + g2 * g2);
    }
    /* Get the bathymetry database gradient                          */
    for (j = js+1; j < je-1; j += jinc) {
      for (i = is+1; i < ie-1; i += iinc) {
	if (grid_xytoij(xyij_tree, lon[i], lat[j], &mi, &mj) == 0) continue;
	if (mi < 1 || mi > geom->nce1-1 || mj < 1 || mj > geom->nce2-1) continue;
	c = geom->map[nz][mj][mi];
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
      if (grid_xytoij(xyij_tree, tslist[n].x, tslist[n].y, &mi, &mj) == 0) continue;
      c = geom->map[nz][mj][mi];
      d1 = fabs(bathy[j][i]) - fabs(geom->botz[c]);
      g1 = fabs((geom->botz[geom->xp1[c]] - geom->botz[geom->xm1[c]]) / 
		(geom->h1au1[c] + geom->h1au1[geom->xp1[c]]));
      g2 = fabs((geom->botz[geom->yp1[c]] - geom->botz[geom->ym1[c]]) / 
		(geom->h1au1[c] + geom->h1au1[geom->yp1[c]]));
      gm = sqrt(g1 * g1 + g2 * g2);
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
		 double *dey,
		 int sz,
		 int mode,
		 double scale
		 )
{
  int c, cc;
  int *vec, nvec;
  
  nvec = (mode) == 2 ? window->b2_t : window->b3_t;
  vec = (mode) == 2 ? window->w2_t : window->w3_t;

  /* Set a no gradient over lateral ghost cells                      */
  c = (mode == 2) ? geom->nbptS : geom->nbpt;
  for (cc = 1; cc <= c; cc++) {
    a[window->bpt[cc]] = a[window->bin[cc]];
  }

  /* Compute e1 and e2 decorrelation length scales                   */
  for (cc = 1; cc < nvec; cc++) {
    c = vec[cc];
    dex[c] = decorr(window, a, sz, c, 0, scale);
    dey[c] = decorr(window, a, sz, c, 1, scale);
  }
}

/* END calc_decorr()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to find the spatial decorrelation distance in e1 or e2    */
/* directions. Based on Romanou at al, (2006) Journal of Climate, 19 */
/* 3378-3393.                                                        */
/*-------------------------------------------------------------------*/
double decorr(geometry_t *window,      /* Processing window          */
	      double *a,
	      int sz,
	      int c,
	      int mode,
	      double scale
	      )
{
  double *dist = (mode) ? window->h2acell : window->h1acell;
  int N2 = ceil(((double)sz * scale) / dist[window->m2d[c]]) / 2;
  int N = 2 * N2 + 1;
  int n, m, cp, cm, *map, *rmap;
  double d1;
  double vm = a[c];
  double p[N+1], v[N+1], d[N+2];
  double eps = 1e-10;

  if (mode) {
    map = window->yp1;
    rmap = window->ym1;
  } else {
    map = window->xp1;
    rmap = window->xm1;
  }

  /* Get the mean */
  cp = cm = c;
  v[N2+1] = a[c];
  vm = a[c];
  for (n = 1; n <= N2; n++) {
    cp = map[cp];
    cm = rmap[cm];
    vm += (a[cp] + a[cm]);
    v[N2+1-n] = a[cm];
    v[N2+1+n] = a[cp];
  }
  vm /= (double)N;

  /* Get the autocorrelated lags */
  for (m = 1; m <= N; m++) {
    d[m] = (m == 1) ? dist[window->m2d[map[cm]]] / 2.0 : d[m-1] + dist[window->m2d[cm]];
    p[m] = 0.0;
    for (n = 1; n <= N - m - 1; n++) {
      p[m] += (v[n] - vm) * (v[n+m] - vm);
    }
    p[m] /= (double)(N-m);
    cm = map[cm];
  }

  /* Find the zero crossing point */
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
