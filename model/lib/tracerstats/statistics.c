/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/tracerstats/statistics.c
 *  
 *  Description:
 *  Tracer statistics routines
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: statistics.c 6519 2020-04-09 12:46:12Z her127 $
 *
 */

#include "tracerstats.h"

/* defined in tracerstats.c */
void append_ts(char* filename, double* vals, char** format, int nvals);
int in_section(void* model, sectioncoord_t* data, int col);
/* returns ij, ij must be valid length 3 int array */
int* i_get_windowijk(void* model, int col, int* ij);

/* defined in this file */
static void sum_section_tracer_u1(trs_t *trs, sectioncoord_t *section,
				                             int index, int n);
static void sum_section_tracer_u2(trs_t *trs, sectioncoord_t *section,
				                             int index, int n);
static void sum_section_tracer_w(trs_t *trs, sectioncoord_t *section,
				                             int index, int n);

#define SIGNUM(x) ( (x < 0.0) ? -1 : 1 )

/*
 * End Local functions
 --------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates the mean of 3d tracer n+1 and stores in tracer n       */
/*-------------------------------------------------------------------*/
void tracer_mean_3d(trs_t *trs, stat3d_work_arrays *w_arr,
		    double ***tr_vals, int n, 
		    int botk, int topk, int **map, void* data) {
  int k;

  /* Get the tracer value                                            */
  double **tr = tr_vals[n + 1];

  /*reset the child flag in a new step */
  if(trs->new_step) {
    w_arr->w2[n] += trs->dt;
  }

  /* Re-initialise at the averaging interval                         */
  if((w_arr->w2[n] > w_arr->w1[n])  && trs->new_step) {
    w_arr->w2[n] = trs->dt;
  }

  /* Get the mean                                                    */
  for(k = botk; k <= topk; k++) {
    *tr_vals[n][k] = (*tr_vals[n][k] * (w_arr->w2[n] - trs->dt) +
		      *tr[k] * trs->dt) / w_arr->w2[n];
  }
}
/* END tracer_mean_3d()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates the run mean of 3d tracer n+1 and stores in tracer n   */
/*-------------------------------------------------------------------*/
void tracer_run_mean_3d(trs_t *trs, stat3d_work_arrays *w_arr,
			double ***tr_vals, int n, 
			int botk, int topk, int **map, void* data) {
  int k;

  /* Get the tracer value                                            */
  double **tr = tr_vals[n + 1];

  /* Get the run mean, averaging over T */
  for(k = botk; k <= topk; k++) {
    *tr_vals[n][k] = (*tr_vals[n][k] * (w_arr->w1[n] - trs->dt) +
		      *tr[k] * trs->dt) / w_arr->w1[n];
  }
}
/* END tracer_run_mean_3d()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates the rmse of 3d tracer n+1 & n+2 and stores in tracer n */
/*-------------------------------------------------------------------*/
void tracer_rmse_3d(trs_t *trs, stat3d_work_arrays *w_arr,
		    double ***tr_vals, int n, 
		    int botk, int topk, int **map, void* data) {
  int k;
  double rmse;

  /* Get the tracer value                                            */
  double **tr1 = tr_vals[n + 1];
  double **tr2 = tr_vals[n + 2];

  /* Re-initialise at the averaging interval                         */
  if((w_arr->w2[n] > w_arr->w1[n])  && trs->new_step) {
    w_arr->w2[n] = trs->dt;
  }

  /* Get the mean                                                    */
  if (w_arr->w2[n]) {
    for(k = botk; k <= topk; k++) {
    /*
    rmse = sqrt((*tr1[k] - *tr2[k]) * (*tr1[k] - *tr2[k]));
    *tr_vals[n][k] = (*tr_vals[n][k] * (w_arr->w2[n] - trs->dt) +
		      rmse * trs->dt) / w_arr->w2[n];
    */
      double mse = *tr_vals[n][k] * *tr_vals[n][k];
      rmse = (*tr1[k] - *tr2[k]) * (*tr1[k] - *tr2[k]);
      *tr_vals[n][k] = sqrt((mse * w_arr->w2[n] + rmse * trs->dt) /
			    (w_arr->w2[n] + trs->dt));

    }
  } else {
    for(k = botk; k <= topk; k++)
      *tr_vals[n][k] = sqrt((*tr1[k] - *tr2[k]) * (*tr1[k] - *tr2[k]));
  }


  /*reset the child flag in a new step */
  if(trs->new_step) {
    w_arr->w2[n] += trs->dt;
  }
}
/* END tracer_rmse_3d()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates the max of 3d tracer n+1 and stores in tracer n        */
/*-------------------------------------------------------------------*/
/*
 * w2 is the running time within the interval
 * w4 is a flag to keep track of the very first value in an interval
 * w1 contains the value of the time interval
 */
void tracer_max_3d(trs_t *trs, stat3d_work_arrays *w_arr, double ***tr_vals, int n, int botk, int topk, int **map, void* data) {
  int k;

  // Setup/Clear flags.  We only do this at the beginning of a new
  // time step otherwise things get messed up as this function is
  // called many times, one for each 2d cell, in fact
  if (trs->new_step) {
    // Check if we're starting a new period or its the very beginning
    if ( (w_arr->w2[n] > w_arr->w1[n]) || (w_arr->w2[n] == 0) ) {
      w_arr->w4[n] = 1;        // set initial flag
      w_arr->w2[n] = trs->dt;  // reset dt count
    }
    
    /* Clear the initial flag, if this is the second time step */
    if (w_arr->w2[n] > trs->dt) {
      w_arr->w4[n] = 0;
    }
    
    /* Move time forward */
    w_arr->w2[n] += trs->dt;
  }

  /* This is the beginning of the interval */
  if (w_arr->w4[n]) {
    /* This is the very first value i.e. its max is itself */
    for(k = botk; k <= topk; k++) {
      *tr_vals[n][k] = *tr_vals[n+1][k];
    }
  } else {
    /* Get the max */
    for(k = botk; k <= topk; k++) {
      *tr_vals[n][k] = max(*tr_vals[n][k], *tr_vals[n+1][k]);
    }
  }
}

/* END tracer_max_3d()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates the min of 3d tracer n+1 and stores in tracer n        */
/*-------------------------------------------------------------------*/
/*
 * w2 is the running time within the interval
 * w4 is a flag to keep track of the very first value in an interval
 * w1 contains the value of the time interval
 *
 * Note: This is an exact replica of the max function above
 */
void tracer_min_3d(trs_t *trs, stat3d_work_arrays *w_arr, double ***tr_vals, int n, int botk, int topk, int **map, void* data) {
  int k;

  // Setup/Clear flags.  We only do this at the beginning of a new
  // time step otherwise things get messed up as this function is
  // called many times, one for each 2d cell, in fact
  if (trs->new_step) {
    // Check if we're starting a new period or its the very beginning
    if ( (w_arr->w2[n] > w_arr->w1[n]) || (w_arr->w2[n] == 0) ) {
      w_arr->w4[n] = 1;        // set initial flag
      w_arr->w2[n] = trs->dt;  // reset dt count
    }
    
    /* Clear the initial flag, if this is the second time step */
    if (w_arr->w2[n] > trs->dt) {
      w_arr->w4[n] = 0;
    }
    
    /* Move time forward */
    w_arr->w2[n] += trs->dt;
  }

  /* This is the beginning of the interval */
  if (w_arr->w4[n]) {
    /* This is the very first value i.e. its min is itself */
    for(k = botk; k <= topk; k++) {
      *tr_vals[n][k] = *tr_vals[n+1][k];
    }
  } else {
    /* Get the min */
    for(k = botk; k <= topk; k++) {
      *tr_vals[n][k] = min(*tr_vals[n][k], *tr_vals[n+1][k]);
    }
  }
}

/* END tracer_min_3d()                                              */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Calculates the variance of 3d tracer map_3d[n]+1 and stores in    */
/* tracer n using the mean located in map_3d[n].                     */
/*-------------------------------------------------------------------*/
void tracer_var_3d(trs_t *trs, stat3d_work_arrays *w_arr, double ***tr_vals, int n, int botk, int topk, int **map, void* data) 
{
  int k;
  int trn = map[n][0];
  double **tr;
  double **var;
  double **mean;
  double sum_t;
  double sum, sum2;
  double nstep;
  
  /* Get the tracer value                                            */
  tr    = tr_vals[trn + 1];
  mean  = tr_vals[trn];
  var   = tr_vals[n];
  nstep = (double)trs->nstep[n];

  /* Re-initialise at the averaging interval                         */
  if(w_arr->w2[trn] == trs->dt) {
    for(k = botk; k <= topk; k++)
      *tr_vals[n][k] = 0.0;
    trs->nstep[n] = 0.0;
  }
  
  /* Get the tracer variance                                         */
  /* var = (sum(x.x) - sum(x).sum(x)/N)/(N-1)                        */
  /*
   * FR: Var = (1/N-1) {sum over N(x-sqr) - N*mean(x)} [1]
   * At the current step we have :
   *    mean including the current x
   *    current x
   *    previous variance
   *
   * The algorithm below expands out [1] in terms of the above
   * variables and computes the variance for the current step.
   *
   * Note: Apparently there are more robust (numercially stable) ways
   * of calculating the variance as outlined in The Art of Computer
   * programming (and Numerical recipies but this one is harder to
   * implement as a running scheme) by Donald Knuth but I haven't had
   * the chance to follow up on this yet.
   * -Farhan 08/09
   */
  if(nstep > 1.0) {
    for(k = botk; k <= topk; k++) {
      /* Get the sum of tracer trn+1 at the current nstep            */
      sum_t = *mean[k] * nstep;
      /* Get the sum of tracer trn+1 up to the previous nstep        */
      sum = sum_t - *tr[k];
      /* Get the sum of trn+1 squared up to the previous nstep       */
      sum = sum * sum / (nstep - 1.0);
      sum2 = *var[k] * (nstep - 2.0) + sum;
      /* Increment the sum of trn+1 squared (i.e. get the sum of     */
      /* trn+1 squared up to the current nstep).                     */
      sum2 += *tr[k] * *tr[k];
      /* Get the variance                                            */
      sum =  (sum_t * sum_t) / nstep;
      *var[k] = (sum2 - sum) / (nstep - 1.0);
    }
  }
  else {
    for(k = botk; k <= topk; k++) {
      *var[k] = 0.0;
    }
  }
}
/* END tracer_var_3d()                                               */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Calculates the standard deviation of 3d tracer map_3d[n]+1 and    */
/* stores in tracer n using the mean located in map_3d[n].           */
/* Calculated using the variance function.                           */
/*-------------------------------------------------------------------*/
void tracer_stdev_3d(trs_t *trs, stat3d_work_arrays *w_arr, double ***tr_vals, int n, int botk, int topk, int **map, void* data) 
{
  int k;
  double **tr = tr_vals[n];
  
  for(k = botk; k <= topk; k++) {
    *tr[k] = *tr[k] * *tr[k];
  }
  tracer_var_3d(trs, w_arr, tr_vals, n, botk, topk, map, data);
  for(k = botk; k <= topk; k++) {
    *tr[k] = (*tr[k] > 0.0) ? sqrt(*tr[k]) : 0.0;
  }
}
/* END tracer_stdev_3d()                                              */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Calculates the sum of 3d tracer n+1 and n+2 and stores in  */
/* tracer n.                                                         */
/*-------------------------------------------------------------------*/
void tracer_sum_3d(trs_t *trs, stat3d_work_arrays *w_arr, double ***tr_vals, int n, int botk, int topk, int **map, void* data) 
{
  int i,k;
  /* double **tr1, **tr2; */
  int nsumtr = *((int*)data);
  
  /* sum the tracer                                              */
  for(k = botk; k <= topk; k++) {
    *tr_vals[n][k] = 0.0;
    for(i = 1;i <= nsumtr;i++) {
      *tr_vals[n][k] += *tr_vals[n+i][k];
    }
  }
}
/* END tracer_sum_3d()                                              */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Calculates the exposure of tracer n+1 using the threshold in n+2, */
/* and stores in tracer n, with the time of exposure in n+3.         */
/*-------------------------------------------------------------------*/
/* This formulation re-initialises in blocks of dt                   */
void tracer_exposureo(trs_t *trs, int* nn, void* data, int c) {
  int k,n;
  double **tr, **tr1, **tr2;
  double fact;

  n = *nn;
  fact = 1.0 / trs->w_wc.w6[n];

  /* Get the tracer value                                            */
  tr = trs->tr_wc[n + 1];
  tr1 = trs->tr_wc[n + 2];
  tr2 = trs->tr_wc[n + 3];

  /* Reset the child flag in a new step. Note a new step is the      */
  /* start of a new time-step (tracerstats are called multiple times */
  /* in a timestep, here we want the first.                          */
  if(trs->new_step ) {
    trs->w_wc.w4[n] = 0;
  }

  /* Reinitialise the exposeure and exposure time if exposure dt has */
  /* elapsed for cells c=2 to sgsiz.                                 */
  if( trs->w_wc.w4[n] ) {
    for(k = trs->botk_wc; k <= trs->topk_wc; k++) {
      if (*tr2[k] <= 0.0) {
	*trs->tr_wc[n][k] = 0.0;
      }
      *tr2[k] = 0.0;
    }
    trs->w_wc.w2[n] = trs->time;
  }

  /* If the exposure dt has elapsed, then check if no exposure has   */
  /* occurred (tr2 = 0). If, so then reinitialise the exposeure and  */
  /* exposure time. This is done for the first point in the grid,    */
  /* (c=1), so set the flag to indicate exposure dt has elapsed      */
  /* (w4[n]=1), then check for no exposure c=2,sgsiz above.          */
  if (trs->time >= trs->w_wc.w1[n] + trs->w_wc.w2[n] && trs->new_step) {
    trs->w_wc.w4[n] = 1;
    for(k = trs->botk_wc; k <= trs->topk_wc; k++) {
      if (*tr2[k] <= 0.0) {
	*trs->tr_wc[n][k] = 0.0;
      }
      *tr2[k] = 0.0;
    }
  }

  /* Get the exposure */
  for(k = trs->botk_wc; k <= trs->topk_wc; k++) {
    if (trs->w_wc.w5[n] & EX_TR) {
      if ((trs->w_wc.w5[n] & EX_GT && *tr[k] > *tr1[k]) ||
	  (trs->w_wc.w5[n] & EX_LT && *tr[k] < *tr1[k])) {
	*trs->tr_wc[n][k] += fabs(*tr[k] - *tr1[k]) * trs->dt * fact;
	*tr2[k] += (fact * trs->dt);
      }
    }
    if (trs->w_wc.w5[n] & EX_VA) {
      if ((trs->w_wc.w5[n] & EX_GT && *tr[k] > trs->w_wc.w3[n]) ||
	  (trs->w_wc.w5[n] & EX_LT && *tr[k] < trs->w_wc.w3[n])) {
	/**trs->tr_wc[n][k] += fabs(*tr[k] - trs->w_wc.w3[n]) * trs->dt * fact;*/
	*trs->tr_wc[n][k] += fabs(*tr[k]) * trs->dt;
	*tr2[k] += (fact * trs->dt);
      }
    }
  }
}

/*-------------------------------------------------------------------*/
/* Calculates the exposure of tracer n+1 using the threshold in n+2, */
/* and stores in tracer n, with the time of exposure in n+3.         */
/*-------------------------------------------------------------------*/
void tracer_reeftemp(trs_t *trs, int* nn, void* data, int c) {
  int k,n;
  double **tr, **tr1, **tr2;
  double fact;
  double tm;
  int ys, mos, ds, hs, mis, ss;  /* Start year, month, day           */
  int startm = 12;               /* Start month for ReefTemp         */
  int endm = 2;                  /* End month for ReefTemp           */

  n = *nn;
  fact = 1.0 / trs->w_wc.w6[n];

  /* Get the tracer value                                            */
  tr = trs->tr_wc[n + 1];
  tr1 = trs->tr_wc[n + 2];

  /* If the calandar date is before 01/12 or after 01/03 then the    */
  /* ReefTemp exposure is zero.                                      */
  tm = trs->time / 86400.0 + trs->w_wc.w2[n];
  tm_to_julsecs(tm, &ys, &mos, &ds, &hs, &mis, &ss);
  if (mos >= endm + 1 && mos <= startm - 1) {
    for(k = trs->botk_wc; k <= trs->topk_wc; k++) {
      *trs->tr_wc[n][k] = 0.0;
    }
    return;
  }

  /* Get the exposure */
  for(k = trs->botk_wc; k <= trs->topk_wc; k++) {
    if (trs->w_wc.w5[n] & EX_TR) {
      if (*tr[k] > *tr1[k])
	*trs->tr_wc[n][k] += fabs(*tr[k] - *tr1[k]) * trs->dt * fact;
    }
    if (trs->w_wc.w5[n] & EX_VA) {
      if (*tr[k] > trs->w_wc.w3[n])
	*trs->tr_wc[n][k] += fabs(*tr[k] - trs->w_wc.w3[n]) * trs->dt * fact;
    }
  }
}
/* END tracer_reeftemp()                                             */
/*-------------------------------------------------------------------*/


/* This formulation re-initialises dt after the threshold is not     */
/* violated.                                                         */
void tracer_exposure(trs_t *trs, int* nn, void* data, int c) {
  int k,n;
  double **tr, **tr1, **tr2;
  double fact;

  n = *nn;
  fact = 1.0 / trs->w_wc.w6[n];

  /* Get the tracer value                                            */
  tr = trs->tr_wc[n + 1];
  tr1 = trs->tr_wc[n + 2];
  tr2 = trs->tr_wc[n + 3];
  /* If the exposure dt has elapsed, then check if no exposure has   */
  /* occurred (tr2 > dt). If, so then reinitialise the exposeure and */
  /* exposure time.                                                  */
  for(k = trs->botk_wc; k <= trs->topk_wc; k++) {
    if (*tr2[k] > trs->w_wc.w1[n] * fact) {
      *trs->tr_wc[n][k] = 0.0;
      *tr2[k] = 0.0;
    }
  }

  /* Get the exposure */
  for(k = trs->botk_wc; k <= trs->topk_wc; k++) {
    if (trs->w_wc.w5[n] & EX_TR) {
      if ((trs->w_wc.w5[n] & EX_GT && *tr[k] > *tr1[k]) ||
	  (trs->w_wc.w5[n] & EX_LT && *tr[k] < *tr1[k])) {
	*trs->tr_wc[n][k] += fabs(*tr[k] - *tr1[k]) * trs->dt * fact;
	*tr2[k] = 0.0;
      } else
	*tr2[k] += (fact * trs->dt);
    }
    if (trs->w_wc.w5[n] & EX_VA) {
      if ((trs->w_wc.w5[n] & EX_GT && *tr[k] > trs->w_wc.w3[n]) ||
	  (trs->w_wc.w5[n] & EX_LT && *tr[k] < trs->w_wc.w3[n])) {
	*trs->tr_wc[n][k] += fabs(*tr[k] - trs->w_wc.w3[n]) * trs->dt * fact;
	*tr2[k] = 0.0;
      } else
	*tr2[k] += (fact * trs->dt);
    }    
  }
}
/* END tracer_exposure()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates the Degree Heating Week of tracer n+1 using the        */
/* threshold in n+2, and stores in tracer n. DHW is a measure of     */
/* coral bleaching.                                                  */
/*-------------------------------------------------------------------*/
void tracer_dhw(trs_t *trs, int* nn, void* data, int c) {
  int k,n;
  double **tr, **tr1;
  double init_scale = 12.0; /* Re-initialisation scale               */
  double dhw, hotspot;

  n = *nn;

  /* Get the tracer value                                            */
  tr = trs->tr_wc[n + 1];
  tr1 = trs->tr_wc[n + 2];

  /* reset the child flag in a new step                              */
  if(trs->new_step ) {
    trs->w_wc.w4[n] = 0;
    trs->w_wc.w2[n] += trs->dt;
  }

  /* Re-initialise at the averaging interval                         */
  if( trs->w_wc.w4[n] ) {
    /* reset the tracers for the other columns than the first one */
    for(k = trs->botk_wc; k <= trs->topk_wc; k++) {
      *trs->tr_wc[n][k] = 0.0;
    }
  }

  /* Re-initialise at the averaging interval.                        */
  if((trs->w_wc.w2[n] > init_scale * trs->w_wc.w1[n]  && trs->new_step)) {
    trs->w_wc.w2[n] = trs->dt;
    trs->w_wc.w4[n] = 1;
    /* reset the tracers for this column */
    for(k = trs->botk_wc; k <= trs->topk_wc; k++) {
      *trs->tr_wc[n][k] = 0.0;
    }
  }

  /* Get the mean                                                    */
  if(trs->w_wc.w2[n]) {
    for(k = trs->botk_wc; k <= trs->topk_wc; k++) {
      hotspot = *tr[k] - *tr1[k];
      dhw = (hotspot > 1.0) ? hotspot : 0.0;
      *trs->tr_wc[n][k] = (*trs->tr_wc[n][k] * (trs->w_wc.w2[n] - trs->dt) +
			   dhw * trs->dt) / trs->w_wc.w2[n];
    }
  } else {
    for(k = trs->botk_wc; k <= trs->topk_wc; k++) {
      hotspot = *tr[k] - *tr1[k];
      dhw = (hotspot > 1.0) ? hotspot : 0.0;
      *trs->tr_wc[n][k] = dhw;
    }
  }
}
/* END tracer_dhw()                                                  */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Calculates the time integral at intervals of dt of the 3d tracer  */
/* n+1 and stores in tracer n.                                       */
/*-------------------------------------------------------------------*/
void tracer_sum_dt_3d(trs_t *trs, int* nn, void* data, int c) {
  int k,n;
  double **tr;

  n = *nn;

  /* Get the tracer value                                            */
  tr = trs->tr_wc[n + 1];
  /*reset the child flag in a new step */
  if(trs->new_step ) {
    trs->w_wc.w4[n] = 0;
    trs->w_wc.w2[n] += trs->dt;
  }

  /* Re-initialise at the averaging interval                         */
  if( trs->w_wc.w4[n] ) {
    /* reset the tracers for the other columns than the first one */
    for(k = trs->botk_wc; k <= trs->topk_wc; k++) {
      *trs->tr_wc[n][k] += *tr[k];
    }    
  }

  /* Re-initialise at the averaging interval                         */
  if((trs->w_wc.w2[n] > trs->w_wc.w1[n]  && trs->new_step)) {
    trs->w_wc.w2[n] = trs->dt;
    trs->w_wc.w4[n] = 1;
    /* reset the tracers for this column */
    for(k = trs->botk_wc; k <= trs->topk_wc; k++) {
      *trs->tr_wc[n][k] += *tr[k];
    }
  }
}
/* END tracer_integrate_3d()                                         */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Calculates the mean of 3d tracer n+1 and stores in tracer n       */
/*-------------------------------------------------------------------*/
void tracer_meanflux_3d(trs_t *trs, double *flux, void* data, int* nn, int c) {
  int k,n;
  double **tr;
  double f;
  n = *nn;

/* temp * /
int doit = 0;
int testij[2][2]= {{18,70},{49,15}};

/ * end temp*/

  /* Get the tracer value                                            */
  tr = trs->tr_wc[n + 1];
/*  tr = ((double**)data)[n+1]; */
/*
if(i_is_ijk(trs->model,c,testij[1],0))
  doit = 1;
else
	doit = 0;
*/
  /*reset the child flag in a new step */
	if(trs->new_step )
	{
		/* emstag(LDEBUG,"tmp","resetting at new timestep %.4f (%.4f) - %u ",trs->time,trs->w_wc.2[n],c); */
		/* reset the re-initialise flag and add the new dt at the first call*/
		trs->w_wc.w4[n] = 0;
		trs->w_wc.w2[n] += trs->dt;
	}

  /* Re-initialise at the averaging interval                         */
  if( trs->w_wc.w4[n] )
  {
/*  	if(doit)
  		emstag(LDEBUG,"tmp","resetting follow on columns %f - %u ",trs->time,c);
*/
  	/* reset the tracers for the other columns than the first one    */
    for(k = trs->botk_wc; k <= trs->topk_wc; k++)
    {
      *trs->tr_wc[n][k] = 0.0;
    }
  }



  /* Re-initialise the averaging interval for the first column       */
  if((trs->w_wc.w2[n] > trs->w_wc.w1[n]  && trs->new_step))
  {
/*  	if(doit)
			emstag(LDEBUG,"tmp","resetting first column %.4f (%.4f) - %u ",trs->time,trs->w_wc.2[n],c);
*/
    trs->w_wc.w2[n] = trs->dt;
    trs->w_wc.w4[n] = 1;

    /* reset the tracers for this column */
    for(k = trs->botk_wc; k <= trs->topk_wc; k++)
    {
      *trs->tr_wc[n][k] = 0.0;
    }
  }

  if((trs->w_wc.w2[n] - trs->dt) && !trs->w_wc.w4[n]) {

    for(k = trs->botk_wc; k <= trs->topk_wc; k++) {
      f = *tr[k] * flux[k];
      *trs->tr_wc[n][k] = ((*trs->tr_wc[n][k] * (trs->w_wc.w2[n] - trs->dt)) + (f * trs->dt)) / trs->w_wc.w2[n];
    }
/*			if(doit)
    		emstag(LDEBUG,"tmp","calculated at %.4f (%.4f) columns %u - %f * %f result %f %f  ",trs->time,trs->w_wc.2[n],c,*tr[trs->topk_wc],flux[trs->topk_wc],*trs->tr_wc[n][trs->topk_wc], (*trs->tr_wc[n][trs->topk_wc] * trs->w_wc.2[n]));
*/
  }
  else {

    for(k = trs->botk_wc; k <= trs->topk_wc; k++) {
      *trs->tr_wc[n][k] = *tr[k] * flux[k];
    }
/*    if(doit)
   		emstag(LDEBUG,"tmp","calculated first at %.4f columns %u - %f * %f result %f %f ",trs->time,c,*trs->tr_wc[n+1][trs->topk_wc],flux[trs->topk_wc],*trs->tr_wc[n][trs->topk_wc], (*trs->tr_wc[n][trs->topk_wc] * trs->w_wc.2[n]));
*/
  }
}
/* END tracer_meanflux_3d()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates the covariance of 3d tracer map_3d[n][0]+1 and         */
/* map_3d[n][1]+1 and stores in tracer n using the means found in    */
/* map_3d[n][0] and map_3d[n][1].                                    */
/*-------------------------------------------------------------------*/
void tracer_cov_3d(trs_t *trs, int* nn, void* data, int c) {
  int k,n = *nn;
  double **tr1, **tr2;
  double **mean1, **mean2;
  int tm1 = trs->map_3d[n][0];
  int tm2 = trs->map_3d[n][1];
  double sum_t1, sum_t2;
  double sum1, sum2;
  double sum, cov;
  double nstep ;


  /* Get the tracer values                                           */
  tr1 = trs->tr_wc[tm1 + 1];
  tr2 = trs->tr_wc[tm2 + 1];
  mean1 = trs->tr_wc[tm1];
  mean2 = trs->tr_wc[tm2];

  nstep = (double)trs->nstep[n];

  /* Re-initialise at the averaging interval                         */
  if(trs->w_wc.w2[tm1] == trs->dt) {
    for(k = trs->botk_wc; k <= trs->topk_wc; k++)
      *trs->tr_wc[n][k] = 0.0;
    trs->nstep[n] = 0.0;
  }

  /* Get the covariance                                              */
  /* cov = sum(x.y) - sum(x).sum(y)/N                                */
  if(nstep > 1.0) {
    for(k = trs->botk_wc; k <= trs->topk_wc; k++) {
      /* Get the sum of tracer tm1+1 at the current nstep            */
      sum_t1 = *mean1[k] * nstep;
      sum_t2 = *mean2[k] * nstep;
      /* Get the sum of tm1+1 * sum of tm2+1 up to previous nstep    */
      sum1 = sum_t1 - *tr1[k];
      sum2 = sum_t2 - *tr2[k];
      /* Get the sum of tm1+1 * tm2+1 up to the previous nstep       */
      sum = sum1 * sum2 / (nstep - 1.0);
      cov = *trs->tr_wc[n][k] * (nstep - 2.0) + sum;
      /* Increment the sum of tm1+1 * tm2+1 (i.e. get the sum of     */
      /* tm1+1 * tm2+1 up to the current nstep).                     */
      cov += *tr1[k] * *tr2[k];
      /* Get the covariance                                          */
      sum =  sum_t1 * sum_t2 / nstep;
      *trs->tr_wc[n][k] = (cov - sum) / (nstep - 1.0);
    }
  }
  else {
    for(k = trs->botk_wc; k <= trs->topk_wc; k++) {
      *trs->tr_wc[n][k] = 0.0;
    }
  }
}
/* END tracer_cov_3d()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates the correlation coefficient of 3d tracers              */
/* map_3d[n][0]+1 and map_3d[n][1]+1 and stores in tracer n using    */
/* the means found in  map_3d[n][0] and map_3d[n][1] and standard    */
/* deviations found in map_3d[n][2] and map_3d[n][3].                */
/*-------------------------------------------------------------------*/
void tracer_corr_3d_pre(trs_t *trs, int* nn, void* dat, int c) {
  int k,n= *nn;
  double **std1, **std2;

  /* Get the tracer values                                           */
  std1 = trs->tr_wc[trs->map_3d[n][2]];
  std2 = trs->tr_wc[trs->map_3d[n][3]];

  /* Get the covariance by multiplying by standard deviation         */
  /* corr = cov(stdev1.stdev2)                                       */
  /*      = (sum(x.y)-sum(x).sum(y)/N)/((N-1).stdev1.stdev2)         */
  for(k = trs->botk_wc; k <= trs->topk_wc; k++) {
      *trs->tr_wc[n][k] *= (*std1[k] * *std1[k]);
  }
}


void tracer_corr_3d_post(trs_t *trs, int* nn, void* data, int c) {
  int k,n = *nn;
  double **std1, **std2;
  double stdev;


  /* Get the tracer values                                           */
  std1 = trs->tr_wc[trs->map_3d[n][2]];
  std2 = trs->tr_wc[trs->map_3d[n][3]];

  /* Get the correlation coefficient                                 */
  /* corr = cov(stdev1.stdev2)                                       */
  /*      = (sum(x.y)-sum(x).sum(y)/N)/((N-1).stdev1.stdev2)         */
  for(k = trs->botk_wc; k <= trs->topk_wc; k++) {
    stdev = *std1[k] * *std2[k];
    *trs->tr_wc[n][k] = (stdev) ? *trs->tr_wc[n][k] / stdev : 0.0;
  }
}

/* END tracer_corr_3d()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates the difference of 3d tracer n+1 and n+2 and stores in  */
/* tracer n.                                                         */
/*-------------------------------------------------------------------*/
void tracer_diff_3d(trs_t *trs, int* nn, void* data, int c) {
  int k,n = *nn;
  double **tr1, **tr2;

  /* Get the tracer values                                           */
  tr1 = trs->tr_wc[n + 1];
  tr2 = trs->tr_wc[n + 2];

  /* Get the difference                                              */
  for(k = trs->botk_wc; k <= trs->topk_wc; k++) {
    *trs->tr_wc[n][k] = *tr1[k] - *tr2[k];
  }
}
/* END tracer_diff_3d()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates the flux of a 3D tracer in the e1 direction            */
/*-------------------------------------------------------------------*/
void tracer_flux3d_e1(trs_t *trs, int* nn, void* data, int c) {
  int k,n;
  double **tr;
  n = *nn;

  /* Get the tracer values                                           */
  tr = trs->tr_wc[n + 1];

  /* Get the difference                                              */
  for(k = trs->botk_wc; k <= trs->topk_wc; k++) {
    *trs->tr_wc[n][k] = *tr[k] * trs->u1flux3d[k];
  }
}

/* END tracer_flux_3d_e1()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates the flux of a 3D tracer in the e2 direction            */
/*-------------------------------------------------------------------*/
void tracer_flux3d_e2(trs_t *trs, int* nn, void* data, int c) {
  int k,n;
  double **tr;
  n = *nn;
  /* Get the tracer values                                           */
  tr = trs->tr_wc[n + 1];

  /* Get the difference                                              */
  for(k = trs->botk_wc; k <= trs->topk_wc; k++) {
    *trs->tr_wc[n][k] = *tr[k] * trs->u2flux3d[k];
  }
}

/* END tracer_flux_3d_e2()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates the flux of a 3D tracer in the vertical direction      */
/*-------------------------------------------------------------------*/
void tracer_flux3d_w(trs_t *trs, int* nn, void* data, int c) {
  int k,n;
  double **tr;
  n = *nn;

  /* Get the tracer values                                           */
  tr = trs->tr_wc[n + 1];

  /* Get the difference                                              */
  for(k = trs->botk_wc; k <= trs->topk_wc; k++) {
    *trs->tr_wc[n][k] = *tr[k] * trs->w[k] * trs->area_w;
  }
}

/* END tracer_flux_3d_w()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates the normalized vertical profile of a tracer            */
/*-------------------------------------------------------------------*/
void tracer_vprof(trs_t *trs, int* nn, void* data, int c) {
  int k,n;
  double **tr;
  double trt, trz;
  n = *nn;

  /* Get the tracer values                                           */
  tr = trs->tr_wc[n + 1];
  trt = *tr[trs->topk_wc];

  /* Get the difference                                              */
  for(k = trs->botk_wc; k <= trs->topk_wc; k++) {
    trz = (trt) ? *tr[k] / trt : 0.0;
    *trs->tr_wc[n][k] = trz;
  }
}

/* END tracer_flux_3d_w()                                            */
/*-------------------------------------------------------------------*/
 
/*----------------------------------------------------------------------*/
/* Calculates the section fluxes of a column given a set of coordinates */
/* in a given direction, SECTIONH for vertical, SECTIONE1 through the   */
/* face e1, SECTIONE2 through the face e2                               */
/*----------------------------------------------------------------------*/
void tracer_section(void* model, trs_t *trs, int* nn, void* data, int c)
{
  int n;
  sectioncoord_t* section = (sectioncoord_t*)data;
  int dir;

  int topk,botk,index;
  n = *nn;

  if(model == NULL)
  {
    emstag(LFATAL,"tracerstats:statistics:tracer_section"," Model null - exiting!" );
    exit(0);
  }

  /* this must stay here -do not move below next lines */
  section->time = trs->time;
  if(section->time < section->startt)
    return;

  index = in_section(model,section,c);
  if(index < 0)
    return;

  topk = section->data[2][index];
  botk = section->data[3][index];
  dir  = section->data[4][index];

  /* 
   * Get the tracer values, note that n refers to the tracer we are
   * interested in and NOT to the target - there is no necessarily a
   * target tracer
   *
   * dir takes precedence here
   */
  if (dir == SECTIONW || section->dir == SECTIONW) {
    for(n = 0;n < section->ntrs;n++) {
      sum_section_tracer_w(trs, section, index, n);
    }
    
  } else if (dir == SECTIONE1 || section->dir == SECTIONE1) {
    for(n = 0;n < section->ntrs;n++) {
      sum_section_tracer_u1(trs, section, index, n);
    }
    
  } else if (dir == SECTIONE2 || section->dir == SECTIONE2) {
    for(n = 0;n < section->ntrs;n++) {
      sum_section_tracer_u2(trs, section, index, n);
    }
  } else {
    emstag(LERROR,
	   "tracerstats:statistics:tracer_section",
	   " Direction is not specified for either the section as a whole OR on a per cell basis!");
    exit(1);
  }
}

/* END tracer_section()                                            */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Calculates the mean of 2d tracer n+1 and stores in tracer n       */
/*-------------------------------------------------------------------*/
void tracer_mean_2d(trs_t *trs, int* nn, void* data) {
  double *tr;
  int n=*nn;
  /* Get the tracer value                                            */
  tr = trs->tr_in[n + 1];

  /* Re-initialise at the averaging interval                         */
  if(trs->w2S[n] > trs->w1S[n] && trs->new_step) {
    trs->w2S[n] = trs->dt;
  }

  /* Get the mean                                                    */
  if(trs->w2S[n]) {
    *trs->tr_in[n] = (*trs->tr_in[n] * trs->w2S[n] +
         *tr * trs->dt) /
                     (trs->w2S[n] + trs->dt);
  }
  else {
    *trs->tr_in[n] = *tr;
  }
  if (trs->new_step)
    trs->w2S[n] += trs->dt;
}

/* END tracer_mean_2d()                                              */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Calculates the run mean of 2d tracer n+1 and stores in tracer n   */
/*-------------------------------------------------------------------*/
void tracer_run_mean_2d(trs_t *trs, int* nn, void* data) {
  double *tr;
  int n=*nn;
  /* Get the tracer value                                            */
  tr = trs->tr_in[n + 1];

  /* Get the run mean                                                */
  *trs->tr_in[n] = (*trs->tr_in[n] * (trs->w1S[n] -  trs->dt) +
		    *tr * trs->dt) / trs->w1S[n] ;
}

/* END tracer_run_mean_2d()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates the rmse of 2d tracer n+1 & n+2 and stores in tracer n */
/*-------------------------------------------------------------------*/
void tracer_rmse_2d(trs_t *trs, int* nn, void* data) {
  double *tr1, *tr2;
  double rmse;
  int n=*nn;
  /* Get the tracer value                                            */
  tr1 = trs->tr_in[n + 1];
  tr2 = trs->tr_in[n + 2];

  /* Re-initialise at the averaging interval                         */
  if(trs->w2S[n] > trs->w1S[n] && trs->new_step) {
    trs->w2S[n] = trs->dt;
  }

  /* Get the mean                                                    */
  if(trs->w2S[n]) {
    /*
    rmse = sqrt((*tr1 - *tr2) * (*tr1 - *tr2));
    *trs->tr_in[n] = (*trs->tr_in[n] * trs->w2S[n] +
         rmse * trs->dt) /
                     (trs->w2S[n] + trs->dt);
    */
    double mse = *trs->tr_in[n] * *trs->tr_in[n];
    rmse = (*tr1 - *tr2) * (*tr1 - *tr2);
    *trs->tr_in[n] = sqrt((mse * trs->w2S[n] + rmse * trs->dt) /
			  (trs->w2S[n] + trs->dt));
  }
  else {
    *trs->tr_in[n] = sqrt((*tr1 - *tr2) * (*tr1 - *tr2));
  }
  if (trs->new_step)
    trs->w2S[n] += trs->dt;
}

/* END tracer_rmse_2d()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates the max of 2d tracer n+1 and stores in tracer n        */
/*                                                                   */
/* Note: This is an exact copy of the mean code above                */
/*-------------------------------------------------------------------*/
void tracer_max_2d(trs_t *trs, int* nn, void* data) {
  double *tr;
  int n=*nn;

  /* Get the tracer value                                            */
  tr = trs->tr_in[n + 1];

  /* Re-initialise at the averaging interval                         */
  if(trs->w2S[n] > trs->w1S[n] && trs->new_step) {
    trs->w2S[n] = trs->dt;
  }

  /* Get the max */
  if(trs->w2S[n]) {
    *trs->tr_in[n] = max(*trs->tr_in[n], *tr);
  } else {
    *trs->tr_in[n] = *tr;
  }

  if (trs->new_step)
    trs->w2S[n] += trs->dt;
}
/* END tracer_max_2d()                                               */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Calculates the min of 2d tracer n+1 and stores in tracer n        */
/*                                                                   */
/* Note: This is an exact copy of the max code above                 */
/*-------------------------------------------------------------------*/
void tracer_min_2d(trs_t *trs, int* nn, void* data) {
  double *tr;
  int n=*nn;

  /* Get the tracer value                                            */
  tr = trs->tr_in[n + 1];

  /* Re-initialise at the averaging interval                         */
  if(trs->w2S[n] > trs->w1S[n] && trs->new_step) {
    trs->w2S[n] = trs->dt;
  }

  /* Get the min */
  if(trs->w2S[n]) {
    *trs->tr_in[n] = min(*trs->tr_in[n], *tr);
  } else {
    *trs->tr_in[n] = *tr;
  }

  if (trs->new_step)
    trs->w2S[n] += trs->dt;
}
/* END tracer_min_2d()                                               */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Calculates the variance of 2d tracer map_2d[n]+1 and stores in    */
/* tracer n using the mean located in map_2d[n].                     */
/*-------------------------------------------------------------------*/
void tracer_var_2d(trs_t *trs, int* nn, void* data) {
  int n=*nn;
  int trn = trs->map_2d[n][0];
  double *tr;
  double sum_t;
  double sum, sum2;
  double nstep ;

  /* Get the tracer value                                            */
  tr = trs->tr_in[trn + 1];

  nstep = (double)trs->nstep[trs->ntr+n];
  /* Re-initialise at the averaging interval                         */
  if(trs->w2S[trn] == trs->dt) {
    *trs->tr_in[n] = 0.0;
    trs->nstep[trs->ntr+n] = 0.0;
  }

  /* Get the tracer standard deviation                               */
  if(nstep > 1.0) {
    /* Get the sum of tracer trn+1 at the current nstep              */
    sum_t = *trs->tr_in[trn] * nstep;
    /* Get the sum of tracer trn+1 up to the previous nstep          */
    sum = sum_t - *trs->tr_in[trn + 1];
    sum = sum * sum / (nstep - 1.0);
    /* Get the sum of trn+1 squared up to the previous nstep         */
    sum2 = *trs->tr_in[n] * (nstep - 2.0) + sum;
    /* Increment the sum of trn+1 squared (i.e. get the sum of trn+1 */
    /* squared up to the current nstep).                             */
    sum2 += *tr * *tr;
    /* Get the variance                                              */
    sum =  (sum_t * sum_t) / nstep;
    *trs->tr_in[n] = (sum2 - sum) / (nstep - 1.0);
  }
  else {
    *trs->tr_in[n] = 0.0;
  }
}
/* END tracer_var_2d()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates the standard deviation of 2d tracer map_2d[n]+1 and    */
/* stores in tracer n using the mean located in map_2d[n].           */
/* Calculated using the variance function.                           */
/*-------------------------------------------------------------------*/
void tracer_stdev_2d(trs_t *trs, int* nn, void* data) {
  int n=*nn;
  double *tr = trs->tr_in[n];


  *tr = *tr * *tr;
  tracer_var_2d(trs, nn, data);
  *tr = (*tr > 0.0) ? sqrt(*tr) : 0.0;
}
/* END tracer_stdev_2d()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates vertical maximum of 3D tracer n3 and stores in 2D      */
/* tracer n.                                                         */
/*-------------------------------------------------------------------*/
void vertical_max(trs_t *trs, int* nn, void* data) {
  int k,n = *nn;
  double **tr, max;
  int n3 = *((int*)data);

  /* Get the tracer value                                            */
  tr = trs->tr_wc[n3];

  /* Initialise                                                      */
  *trs->tr_in[n] = -1e10;

  /* Integrate the water column                                      */
  for(k = trs->botk_wc; k <= trs->topk_wc; k++) {
    *trs->tr_in[n] = max(*tr[k],*trs->tr_in[n]);
  }
}
/* END vertical_max()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates vertical integral of 3D tracer n3 and stores in 2D     */
/* tracer n.                                                         */
/*-------------------------------------------------------------------*/
void vertical_integrals(trs_t *trs, int* nn, void* data) {
  int k,n = *nn;
  double **tr;
  int n3 = *((int*)data);

  /* Get the tracer value                                            */
  tr = trs->tr_wc[n3];

  /* Initialise                                                      */
  *trs->tr_in[n] = 0.0;

  /* Integrate the water column                                      */
  for(k = trs->botk_wc; k <= trs->topk_wc; k++) {
    *trs->tr_in[n] += *tr[k] * trs->dz_wc[k];
  }
}
/* END vertical_integrals()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates vertical mean of 3D tracer n3 and stores in 2D     */
/* tracer n.                                                         */
/*-------------------------------------------------------------------*/
void vertical_mean(trs_t *trs, int* nn, void* data) {
  int k,n3,n= *nn;
  double **tr;
  n3 = *((int*)data);
  /* Get the tracer value                                            */
  tr = trs->tr_wc[n3];

  /* Initialise                                                      */
  *trs->tr_in[n] = 0.0;

  /* Integrate the water column                                      */
  for(k = trs->botk_wc; k <= trs->topk_wc; k++) {
    *trs->tr_in[n] += *tr[k] * trs->dz_wc[k];
  }
  *trs->tr_in[n] = (trs->depth_wc) ? *trs->tr_in[n] / trs->depth_wc : 0.0;

}
/* END vertical_mean()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Extracs a copy of a layer of 3D tracer n3 and stores in 2D        */
/* tracer n.                                                         */
/*-------------------------------------------------------------------*/
void copy_layer(trs_t *trs, int* nn, void* data) {
  int k,n3,n= *nn;
  double **tr;
  int layer = (int)trs->w1S[n];

  n3 = *((int*)data);
  /* Get the tracer value                                            */
  tr = trs->tr_wc[n3];

  /* Initialise                                                      */
  *trs->tr_in[n] = 0.0;

  if (layer == -1)
    *trs->tr_in[n] = *tr[trs->topk_wc];
  else if (layer == -2)
    *trs->tr_in[n] = *tr[trs->botk_wc];
  else 
    *trs->tr_in[n] = *tr[layer];

}
/* END copy_layer()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates the difference of 2d tracer n+1 and n+2 and stores in  */
/* tracer n.                                                         */
/*-------------------------------------------------------------------*/
void tracer_diff_2d(trs_t *trs, int* nn, void* data) {
  int n=*nn;
  double *tr1, *tr2;

  /* Get the tracer values                                           */
  tr1 = trs->tr_in[n + 1];
  tr2 = trs->tr_in[n + 2];

  /* Get the difference                                              */
  *trs->tr_in[n] = *tr1 - *tr2;
}
/* END tracer_diff_2d()                                              */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Calculates the sum of 2d tracer n+1 and n+2 and stores in  */
/* tracer n.                                                         */
/*-------------------------------------------------------------------*/
void tracer_sum_2d(trs_t *trs, int* nn, void* data) {
  int i,n=*nn;
/*  double *tr1, *tr2; */
  int nsumtr = *((int*)data);
  
	*trs->tr_in[n] = 0.0;
	
  /* sum the tracer                                              */
	for(i = 1;i <= nsumtr;i++)
	{	
  	*trs->tr_in[n] += *trs->tr_in[n+i];
	}



  /* Get the tracer values                                           */
/*  tr1 = trs->tr_in[n + 1];
  tr2 = trs->tr_in[n + 2];
*/
  /* sum the tracer                                          */
/*
  *trs->tr_in[n] = *tr1 + *tr2;
  */
}
/* END tracer_sum_2d()                                              */
/*-------------------------------------------------------------------*/



/*-------------------------------------------------------------------*/
/* Calculates the covariance of 2d tracer map_2d[n][0]+1 and         */
/* map_2d[n][1]+1 and stores in tracer n using the means found in    */
/* map_2d[n][0] and map_2d[n][1].                                    */
/*-------------------------------------------------------------------*/
void tracer_cov_2d(trs_t *trs, int* nn, void* data) {
  double *tr1, *tr2;
  int n=*nn;
  int tm1 = trs->map_2d[n][0];
  int tm2 = trs->map_2d[n][1];
  double sum_t1, sum_t2;
  double sum1, sum2;
  double sum, cov;
  double nstep ;

  /* Get the tracer values                                           */
  tr1 = trs->tr_in[tm1 + 1];
  tr2 = trs->tr_in[tm2 + 1];
  nstep = (double)trs->nstep[trs->ntr+n];

  /* Re-initialise at the averaging interval                         */
  if(trs->w2S[tm1] == trs->dt) {
    *trs->tr_in[n] = 0.0;
    trs->nstep[trs->ntr+n] = 0.0;
  }

  /* Get the covariance                                              */
  if(nstep > 1.0) {
    /* Get the sum of tracer tm1+1 at the current nstep              */
    sum_t1 = *trs->tr_in[tm1] * nstep;
    sum_t2 = *trs->tr_in[tm2] * nstep;
    /* Get the sum of tm1+1 * sum of tm2+1 up to the previous nstep  */
    sum1 = sum_t1 - *tr1;
    sum2 = sum_t2 - *tr2;
    sum = sum1 * sum2 / (nstep - 1.0);
    /* Get the sum of tm1+1 * tm2+1 up to the previous nstep         */
    cov = *trs->tr_in[n] * (nstep - 2.0) + sum;
    /* Increment the sum of tm1+1 * tm2+1 (i.e. get the sum of tm1+1 */
    /* * tm2+1 up to the current nstep).                             */
    cov += *tr1 * *tr2;
    /* Get the covariance                                            */
    sum =  sum_t1 * sum_t2 / nstep;
    *trs->tr_in[n] = (cov - sum) / (nstep - 1.0);
  }
  else {
    *trs->tr_in[n] = 0.0;
  }

}
/* END tracer_cov_2d()                                               */
/*-------------------------------------------------------------------*/



/*-------------------------------------------------------------------*/
/* Calculates the correlation coefficient of 2d tracers              */
/* map_2d[n][0]+1 and map_2d[n][1]+1 and stores in tracer n using    */
/* the means found in  map_2d[n][0] and map_2d[n][1] and standard    */
/* deviations found in map_2d[n][2] and map_2d[n][3].                */
/*-------------------------------------------------------------------*/
void tracer_corr_2d_pre(trs_t *trs, int* nn, void* data) {
  double *std1, *std2;
  int n=*nn;

  /* Get the tracer values                                           */
  std1 = trs->tr_in[trs->map_2d[n][2]];
  std2 = trs->tr_in[trs->map_2d[n][3]];

  /* Get the covariance by multiplying by standard deviation         */
  /* corr = cov(stdev1.stdev2)                                       */
  /*      = (sum(x.y)-sum(x).sum(y)/N)/((N-1).stdev1.stdev2)         */
  *trs->tr_in[n] *= (*std1 * *std1);
}

void tracer_corr_2d_post(trs_t *trs, int* nn, void* data) {
  double *std1, *std2;
  double stdev;
  int n=*nn;

  /* Get the tracer values                                           */
  std1 = trs->tr_in[trs->map_2d[n][2]];
  std2 = trs->tr_in[trs->map_2d[n][3]];

  /* Get the correlation coefficient                                 */
  /* corr = cov(stdev1.stdev2)                                       */
  /*      = (sum(x.y)-sum(x).sum(y)/N)/((N-1).stdev1.stdev2)         */
  stdev = *std1 * *std2;
  *trs->tr_in[n] = (stdev) ? *trs->tr_in[n] / stdev : 0.0;
}

/* END tracer_corr_2d()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Copys the value of 3d tracer n+1 and stores in tracer n       */
/*-------------------------------------------------------------------*/
/*
 * void tracer_copy_3d(trs_t *trs, int n, void* data, int c) {
  int k;
  double **tr;

  / * Get the tracer value                                            * /

  tr = trs->tr_wc[n + 1];

  for(k = trs->botk_wc; k <= trs->topk_wc; k++)
      *trs->tr_wc[n][k] = *trs->tr_wc[n+1][k];


}*/
/* END tracer_copy_3d()                                              */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Calculates the vertical difference of cell content of a column
 *  given two ranges of depth/layers                                 */
/*-------------------------------------------------------------------*/
void tracer_vertical_strat_vol(trs_t *trs, int* nn, void* data)
{
  int n;
  int k;
  double **tr;
  double tt,bb;
  vdiff_range_t* vdiff=(vdiff_range_t*)data;
  n = *nn;
  tt = 0.0;
  bb = 0.0;
  /* Get the 3d tracer value                                            */
  tr = trs->tr_wc[vdiff->ntr3d];

  /* Initialise                                                      */
  *trs->tr_in[n] = 0.0;


  /* Integrate the water column                                      */
  for(k = 0; k < vdiff->nt; k++) {
    if(trs->topk_wc >= vdiff->toprange[k] && vdiff->toprange[k] >= trs->botk_wc )/*&& tr[vdiff->toprange[k]] != NULL)*/
    {
        tt +=  *tr[vdiff->toprange[k]] * trs->dz_wc[vdiff->toprange[k]] * trs->area_w;
    }else if(vdiff->strict)
    {
      /*emstag(LTRACE,"tracerstats:statistics:tracer_vertical_strat"," top incomplete skipping %d  - %f",vdiff->toprange[k],tt);*/
      tt = 0.0;
      bb = -1;
      break;
    }
  }

  if(bb == -1)
    bb = 0.0;
  else
  {
    for(k = 0; k < vdiff->nb; k++) {
      if(trs->topk_wc >= vdiff->botrange[k] && vdiff->botrange[k] >= trs->botk_wc )
      {
          bb +=  *tr[vdiff->botrange[k]] * trs->dz_wc[vdiff->botrange[k]] * trs->area_w;
      }else if(vdiff->strict)
      {
        /*emstag(LTRACE,"tracerstats:statistics:tracer_vertical_strat"," bottom incomplete skipping %d - %f",vdiff->botrange[k],bb);*/
        bb = 0.0;
        tt = 0.0;
        break;
      }
    }
  }

  *(trs->tr_in[n]) = tt - bb;
}


/*-------------------------------------------------------------------*/
/* Calculates the vertical difference of values in a cell of a
 * column given two ranges of depth/layers                           */
/*-------------------------------------------------------------------*/
void tracer_vertical_strat(trs_t *trs, int* nn, void* data)
{
  int n;
  int k;
  double **tr;
  double tt,bb;
  vdiff_range_t* vdiff=(vdiff_range_t*)data;
  n = *nn;
  tt = 0.0;
  bb = 0.0;
  /* Get the 3d tracer value                                            */
  tr = trs->tr_wc[vdiff->ntr3d];

  /* Initialise                                                      */
  *trs->tr_in[n] = 0.0;


  /* Integrate the water column                                      */
  for(k = 0; k < vdiff->nt; k++) {
    if(trs->topk_wc >= vdiff->toprange[k] && vdiff->toprange[k] >= trs->botk_wc )/*&& tr[vdiff->toprange[k]] != NULL)*/
    {
        tt +=  *tr[vdiff->toprange[k]];
    }else if(vdiff->strict)
    {
      /*emstag(LTRACE,"tracerstats:statistics:tracer_vertical_strat"," top incomplete skipping %d  - %f",vdiff->toprange[k],tt);*/
      tt = 0.0;
      bb = -1;
      break;
    }
  }

  if(bb == -1)
    bb = 0.0;
  else
  {
    for(k = 0; k < vdiff->nb; k++) {
      if(trs->topk_wc >= vdiff->botrange[k] && vdiff->botrange[k] >= trs->botk_wc )
      {
          bb +=  *tr[vdiff->botrange[k]];
      }else if(vdiff->strict)
      {
        /*emstag(LTRACE,"tracerstats:statistics:tracer_vertical_strat"," bottom incomplete skipping %d - %f",vdiff->botrange[k],bb);*/
        bb = 0.0;
        tt = 0.0;
        break;
      }
    }
  }

  *(trs->tr_in[n]) = tt - bb;
}

/*-------------------------------------------------------------------*/
/* Calculates the vertical ratio of integrated values of a cell
 * within a column
 * given two ranges of depth/layers                                  */
/*-------------------------------------------------------------------*/
void tracer_vertical_strat_voldiv(trs_t *trs, int* nn, void* data)
{
  int n;
  int k;
  double **tr;
  double tt,bb;
  vdiff_range_t* vdiff=(vdiff_range_t*)data;
  n = *nn;
  tt = 0.0;
  bb = 0.0;
  /* Get the 3d tracer value                                            */
  tr = trs->tr_wc[vdiff->ntr3d];

  /* Initialise                                                      */
  *trs->tr_in[n] = 0.0;


  /* Integrate the water column                                      */
  for(k = 0; k < vdiff->nt; k++) {
    if(trs->topk_wc >= vdiff->toprange[k] && vdiff->toprange[k] >= trs->botk_wc )/*&& tr[vdiff->toprange[k]] != NULL)*/
    {
        tt +=  *tr[vdiff->toprange[k]] * trs->dz_wc[vdiff->toprange[k]]* trs->area_w;
    }else if(vdiff->strict)
    {
      /*emstag(LTRACE,"tracerstats:statistics:tracer_vertical_strat_voldiff"," top incomplete skipping %d  - %f",vdiff->toprange[k],tt);*/
      tt = 0.0;
      bb = -1;
      break;
    }
  }

  if(bb == -1)
    bb = 0.0;
  else
  {
    for(k = 0; k < vdiff->nb; k++) {
      if(trs->topk_wc >= vdiff->botrange[k] && vdiff->botrange[k] >= trs->botk_wc )
      {
          bb +=  *tr[vdiff->botrange[k]] * trs->dz_wc[vdiff->botrange[k]]* trs->area_w;
      }else if(vdiff->strict)
      {
        /*emstag(LTRACE,"tracerstats:statistics:tracer_vertical_strat_voldiff"," bottom incomplete skipping %d - %f",vdiff->botrange[k],bb);*/
        bb = 0.0;
        tt = 0.0;
        break;
      }
    }
  }

  *(trs->tr_in[n]) = tt / bb;
}


/*-------------------------------------------------------------------*/
/* Calculates the vertical ratio of values of a cell within a column
 * given two                                                         */
/* ranges of depth/layers                                            */
/*-------------------------------------------------------------------*/
void tracer_vertical_strat_div(trs_t *trs, int* nn, void* data)
{
  int n;
  int k;
  double **tr;
  double tt,bb;
  vdiff_range_t* vdiff=(vdiff_range_t*)data;
  n = *nn;
  tt = 0.0;
  bb = 0.0;
  /* Get the 3d tracer value                                            */
  tr = trs->tr_wc[vdiff->ntr3d];

  /* Initialise                                                      */
  *trs->tr_in[n] = 0.0;


  /* Integrate the water column                                      */
  for(k = 0; k < vdiff->nt; k++) {
    if(trs->topk_wc >= vdiff->toprange[k] && vdiff->toprange[k] >= trs->botk_wc )/*&& tr[vdiff->toprange[k]] != NULL)*/
    {
        tt +=  *tr[vdiff->toprange[k]];
    }else if(vdiff->strict)
    {
      /*emstag(LTRACE,"tracerstats:statistics:tracer_vertical_strat_diff"," top incomplete skipping %d  - %f",vdiff->toprange[k],tt);*/
      tt = 0.0;
      bb = -1;
      break;
    }
  }

  if(bb == -1)
    bb = 0.0;
  else
  {
    for(k = 0; k < vdiff->nb; k++) {
      if(trs->topk_wc >= vdiff->botrange[k] && vdiff->botrange[k] >= trs->botk_wc )
      {
          bb +=  *tr[vdiff->botrange[k]];
      }else if(vdiff->strict)
      {
        /*emstag(LTRACE,"tracerstats:statistics:tracer_vertical_strat_diff"," bottom incomplete skipping %d - %f",vdiff->botrange[k],bb);*/
        bb = 0.0;
        tt = 0.0;
        break;
      }
    }
  }

  *(trs->tr_in[n]) = tt / bb;
}

void tracer_vertical_sedflux1(trs_t *trs, int* nn, void* data)
{
  int i, n;

  vdiff_range_t* vdiff=(vdiff_range_t*)data;
  n = *nn;
  /* calculate the mass at this point  in the bottom cell */

  *trs->tr_in[n] = 0.0;

  for(i = 0;i < trs->sednz; i++)
  {
    *trs->tr_in[n] += *trs->tr_sed[vdiff->ntr3d][i] * trs->dz_sed[i] * trs->area_w;
  }
}


void tracer_vertical_sedflux2(trs_t *trs, int* nn, void* data)
{
  int n,i;
  double tt = 0.0;
  vdiff_range_t* vdiff=(vdiff_range_t*)data;
  n = *nn;

  /* last reset + flux dt > current time */
  if(vdiff->botrange[1] + vdiff->botrange[0] < trs->time )
  {
    *trs->tr_in[n] = 0.0;
    vdiff->botrange[1] = trs->time;
    vdiff->botrange[2] = 0;
  }
  /* calculate the mass at this point  in the bottom cell */

  for(i = 0;i < trs->sednz; i++)
  {
    tt += *trs->tr_sed[vdiff->ntr3d][i] * trs->dz_sed[i] * trs->area_w;
  }


  /* *trs->tr_wc[vdiff->ntr3d][trs->botk_wc] * trs->dz_wc[trs->botk_wc] * trs->area;
  */
  if(vdiff->mode == VDIFFSEDFLUX3)
  {

    *trs->tr_in[n] = -1 * ((*trs->tr_in[n] * vdiff->botrange[1] +
                     (tt - *trs->tr_in[vdiff->nb])) / (vdiff->botrange[2] + trs->dt));
    vdiff->botrange[2] += trs->dt;
  }
  else
    *trs->tr_in[n] += -1 *((tt - *trs->tr_in[vdiff->nb]));

}


/*
 * Local helper functions
 * These sum up all the cells that make up the section
 *
 * There are 2 modes:
 *    1) The original mode integrates over time and space as:
 *        = conc(units/m3)  x area(m2) x normal vel (m/s) * dt (s) [1]
 *                            ^^^^^^^^^^^^^^^^^^^^^^^^^^^
 *                          this is defined as the vel flux
 *          ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 *      this is defined as the concentration flux of the tracer
 *        = result will be in 'units'
 *
 *    2) The second mode is the average flux through the face as
 *       follows:
 *         A   = sum[i]a(i) - total area of the face, i:=cell number
 *         Nt := number of time steps
 *              Note: in the code this is sometimes written as:
 *              Nt * dt / dt
 *         c := tracer concentration as a function of (i,t)
 *         c_flux := c * vel_flux;
 *
 *         Mean concentration on the face
 *         C  = sum[i]sum[t]a(i)c(i,t) / (A.Nt)
 *         units: units of tracer concentration, := unts
 *
 *         Mean concentration flux through the face
 *         CF = sum[i]sum[t]a(i)c_flux / (A.Nt)
 *              Note: this is implemented as :
 *              [1] / (A.Nt.dt)
 *         units: unts/m2/s
 *
 *         Covariance of (C,Usgn) 
 *         CUsgn = sum[i]sum[t]a(i)sgn(u(i))(c(i,t) - C)
 *                 Note: Since we don't have C ahead of time, we
 *                       calculate and store the expanded terms and do
 *                       the mult x C of the second and take away from
 *                       the first in section_scatter
 *         units: unts^2
 *    
 *     And then there are some velocity calculations for 2) as well:
 *       U    = sum[i]sum[t]a(i)u(i,t)      / (A.Nt)
 *       Uabs = sum[i]sum[t]a(i)abs(u(i,t)) / (A.Nt)
 *       Usq  = sum[i]sum[t]a(i)u(i,t)^2    / (A.Nt)
 *
 */
static void sum_section_tracer_u1(trs_t *trs, sectioncoord_t *section, 
				  int index, int n)
{
  int k,kv;
  int topk = section->data[2][index];
  int botk = section->data[3][index];
  int ntr  = section->trsn[n];

  if(topk < botk)
    kv = trs->topk_wc;
  else
    kv = topk;
  if(botk < trs->botk_wc)
    botk = trs->botk_wc;

  for(k = botk; k <= kv && k <= trs->topk_wc ; k++) {
    // This is the integrated conc flux
    section->results[n] += *trs->tr_wc[ntr][k] * trs->u1flux3d[k] * trs->dt;

    if (section->mode == SECTION_AVERAGE) {
      // averaged concentrations
      section->results_c[n] += *trs->tr_wc[ntr][k] * trs->area_e1[k];

      // velocities
      section->u    +=      trs->u1[k]     * trs->area_e1[k];
      section->uabs += fabs(trs->u1[k])    * trs->area_e1[k];
      section->usq  += pow( trs->u1[k], 2) * trs->area_e1[k];

      // Tracer squared
      section->results_sq[n] += pow(*trs->tr_wc[ntr][k], 2) * trs->area_e1[k];

      // CUsgn part a
      section->results_cusgn_a[n] += trs->area_e1[k] * SIGNUM(trs->u1[k]) *
	                                                 *trs->tr_wc[ntr][k];
      // CUsgn part b
      section->results_cusgn_b[n] += trs->area_e1[k] * SIGNUM(trs->u1[k]);

      // Sum the area
      section->Tot_area += trs->area_e1[k];
    }
  }
}

static void sum_section_tracer_u2(trs_t *trs, sectioncoord_t *section, 
				  int index, int n)
{
  int k,kv;
  int topk = section->data[2][index];
  int botk = section->data[3][index];
  int ntr  = section->trsn[n];

  if(topk < botk)
    kv = trs->topk_wc;
  else
    kv = topk;
  if(botk < trs->botk_wc)
    botk = trs->botk_wc;

  for(k = botk; k <= kv && k <= trs->topk_wc ; k++) {
    // This is the integrated conc flux
    section->results[n] += *trs->tr_wc[ntr][k] * trs->u2flux3d[k] * trs->dt;

    if (section->mode == SECTION_AVERAGE) {
      // averaged concentration
      section->results_c[n]  += *trs->tr_wc[ntr][k] * trs->area_e2[k];
      // averaged concentration squared
      section->results_sq[n] += pow(*trs->tr_wc[ntr][k],2) * trs->area_e2[k];
      // velocities
      section->u    +=      trs->u2[k]     * trs->area_e2[k];
      section->uabs += fabs(trs->u2[k])    * trs->area_e2[k];
      section->usq  += pow( trs->u2[k], 2) * trs->area_e2[k];

      // CUsgn part a
      section->results_cusgn_a[n] += trs->area_e2[k] * SIGNUM(trs->u2[k]) *
	                                                 *trs->tr_wc[ntr][k];
      // CUsgn part b
      section->results_cusgn_b[n] += trs->area_e2[k] * SIGNUM(trs->u2[k]);

      // Sum the area
      section->Tot_area += trs->area_e2[k];
    }
  }
}


static void sum_section_tracer_w(trs_t *trs, sectioncoord_t *section, 
				 int index, int n)
{
  int kv;
  int topk = section->data[2][index];
  int botk = section->data[3][index];
  int ntr  = section->trsn[n];

  if(topk < 0) {
    /* if we have negative values, deduct from top layer */
    kv = trs->topk_wc + topk;
    /* 
     * if we are to low take the second bottom layer and get the flux
     * between the last two water layers 
     */
    if(kv <= trs->botk_wc)
      kv = trs->botk_wc + 1;
  } else
    kv = topk;
  
  // This is the integrated conc flux
  // We don't have wflux3d so we need to multiply by the area * w
  section->results[n] += *trs->tr_wc[ntr][kv] * trs->w[kv] * trs->area_w 
                                                                  * trs->dt;

  // Keep track of velocities, if reqd
  if (section->mode == SECTION_AVERAGE) {
    // averaged concentrations
    section->results_c[n] += *trs->tr_wc[ntr][kv] * trs->area_w;

    // velocities
    section->u    +=      trs->w[kv]     * trs->area_w;
    section->uabs += fabs(trs->w[kv])    * trs->area_w;
    section->usq  += pow( trs->w[kv], 2) * trs->area_w;
    // Tracer squared
    section->results_sq[n] += pow(*trs->tr_wc[ntr][kv],2) * trs->area_w;
    // CUsgn part a
    section->results_cusgn_a[n] += trs->area_w * SIGNUM(trs->w[kv]) *
                                                      *trs->tr_wc[ntr][kv];
    // CUsgn part b
    section->results_cusgn_b[n] += trs->area_w * SIGNUM(trs->w[kv]);
    
    // Sum the area
    section->Tot_area += trs->area_w;
  }
}

// EOF
