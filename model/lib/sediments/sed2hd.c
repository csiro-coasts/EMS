/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/sediments/sed2hd.c
 *  
 *  Description:
 *  Routine to handle data flux from sediment to hydromodule
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: sed2hd.c 7572 2024-05-30 03:34:39Z riz008 $
 *
 */

#ifdef __cplusplus
#include "stdafx.h"
extern "C" {
#endif

#include "sediments.h"
#if defined(HAVE_SEDIMENT_MODULE)


int sinterface_put_ustrcw(void* hmodel, int c, double ustrcw);
static void sed2hd_internal(sediment_t *sediment, sed_column_t *sm, int c);
void sinterface_putgridz_sed(void* hmodel, int c, double *gridz_sed);
void sinterface_putcellz_sed(void* hmodel, int c, double *cellz_sed);
void sinterface_putgridz_wc(void* hmodel, int c, double *gridz_wc);
void sinterface_putcellz_wc(void* hmodel, int c, double *cellz_wc);
void sinterface_putdz_wc(void* hmodel, int c, double *dz_wc) ;
void sinterface_puttopz_wc(void* hmodel, int c, double topz_wc);
void sinterface_putbotz_wc(void* hmodel, int c, double botz_wc);
int sinterface_gettopk_wc(void* hmodel, int c) ;
int sinterface_getbotk_wc(void* hmodel, int c) ;
void sinterface_putdiagtracer_2d(void* hmodel, int c,
double hripple,int n_hripple,  double lripple, int n_lripple,
double ustrcw_skin, int n_ustrcw_skin,
double depth_sedi, int n_depth_sedi, double dzactive,  int n_dzactive,
double erdepflux_total, int n_erdepflux_total, int n_erdepflux_total_ac,
double erdepflux_oxygen, int n_erdepflux_oxygen, int n_erdepflux_oxygen_ac);

void cohsedcontent(sediment_t *sediment, sed_column_t *sm);
void diagnostics(sediment_t *sediment, sed_column_t *sm, const char *comment);
void mass_balance(sediment_t *sediment, sed_column_t *sm);

double *sinterface_getpointerBtracer(void* hmodel, int n, int c);


/* Update hydro using sediment data */
void sed2hd(sediment_t *sediment, sed_column_t *sm, int c)
{
  void *hmodel = sediment->hmodel;
  int k, n, m;
  int tk,bk;
  sed_params_t *param = sediment->msparam;
  int col_index = sm->col_number-1;
  /* Copy internal spatial sediment variables */
  sed2hd_internal(sediment, sm, col_index);
  /* wc */
  sm->sed_start=0;
  k = sinterface_put_ustrcw(hmodel, c, sm->ustrcw_skin);
  /* put wc tracers */
  tk =sinterface_gettopk_wc( hmodel, c) ;
  bk =sinterface_getbotk_wc( hmodel, c) ;
  cohsedcontent(sediment, sm);
 
  for(n=0; n < param->ntr; n++) {
     sed_tracer_t *tracer = &sediment->mstracers[n];
     if (tracer->diagn) {
       if (strcmp(tracer->name, "tss") == 0 ) {
	   for(k=sm->botk_wc;k<=sm->topk_wc;k++)
	     sm->tr_wc[n][k]=sm->tss_wc[k];
	   for(k=sm->botk_sed;k<=sm->topk_sed;k++)
	     sm->tr_sed[n][k]=sm->tss_sed[k];
	 }
       else if (strcmp(tracer->name, "svel_floc") == 0) {
         for(k=sm->botk_wc;k<=sm->topk_wc;k++)
	   sm->tr_wc[n][k]=sm->svel_floc[k];
         for(k=sm->botk_sed;k<=sm->topk_sed;k++)
	   sm->tr_sed[n][k]=0.5*(sm->svel_consolid[k]+sm->svel_consolid[k+1]);
       }
       else if (strcmp(tracer->name, "porosity") == 0) {
       for(k=sm->botk_wc;k<=sm->topk_wc;k++)
         sm->tr_wc[n][k]=sm->por_wc[k];
       for(k=sm->botk_sed;k<=sm->topk_sed;k++)
         sm->tr_sed[n][k]=sm->por_sed[k];
       }
       else if (strcmp(tracer->name, "cohsed") == 0) {
       for(k=sm->botk_wc;k<=sm->topk_wc;k++)
         sm->tr_wc[n][k]=sm->coh_wc[k];
       for(k=sm->botk_sed;k<=sm->topk_sed;k++)
         sm->tr_sed[n][k]=sm->coh_sed[k];
       }
     }
  }

  /* fill up empty layers with surface value */
  if(sm->topk_wc < param->nz-2){
    for(n=0; n < param->ntr; n++) {
    sed_tracer_t *tracer = &sediment->mstracers[n];
    for(k=sm->topk_wc+1;k<param->nz-1;k++)
      sm->tr_wc[n][k] = sm->tr_wc[n][sm->topk_wc];
    }
  }

  // move 1d arrays to 3d
  if(tk > bk)
  {
     for(n=0; n < param->ntr; n++) {
       sed_tracer_t *tr = &sediment->mstracers[n];
        /* do not update water column T,S in hd wc,
     if mecosed is decoupled from hd*/
       if(strcmp("temp",tr->name) != 0 && strcmp("salt",tr->name) != 0) {
          for(k=sm->botk_wc;k<=sm->topk_wc;k++)
            *sm->ptr_wc[n][k] = tr->u_scale*sm->tr_wc[n][k];
        }
        else
          continue;
      }
   }
   else  if (tk == bk) {
    /* If there was one water colum layer in hydromodule, then
       merge two layers, created by mecosed into one and move
       data to hydromodule. Note that the mixing procedure used
       below is correct, as meco assumes that
       porosity in water column = 1 (geomorph=0) */
      for(n=0; n < param->ntr; n++) {
        sed_tracer_t *tr = &sediment->mstracers[n];
        if(strcmp("temp",tr->name) != 0 && strcmp("salt",tr->name) != 0) {
          *sm->ptr_wc[n][tk] = tr->u_scale*0.5*(sm->tr_wc[n][0]+sm->tr_wc[n][1]);
        }
        else
          continue;
      }
    } else if (tk<bk) {
       sedtag(LWARN,"sed:sed2hd:sed2hd"," k_top < k_bot in wc \n");
    }
      /*NMYfix*/
    if (param->geomorph) {
      /*
      sinterface_putgridz_wc(hmodel, c, sm->gridz_wc);
      sinterface_putcellz_wc(hmodel, c, sm->cellz_wc);
      */
      sinterface_putdz_wc(hmodel, c, sm->dz_wc);
      sinterface_puttopz_wc(hmodel, c, sm->topz_wc);
      sinterface_putbotz_wc(hmodel, c, sm->botz_wc);
    }
   
    /* sed bed */
    /* Diagnostics */
    /* Fluxes */
    sm->erdepflux_total = 0.;
    for(n=0; n < param->ntr; n++) {
      sed_tracer_t *tr = &sediment->mstracers[n];
      //	if (!tr->prmspatial) {
	  for(k=sm->botk_sed;k<=sm->topk_sed;k++)
	    *sm->ptr_sed[n][k] = tr->u_scale*sm->tr_sed[n][k];
	  if (!tr->dissol && !tr->diagn) {
	    sm->erdepflux_total += sm->erdepflux[n];
	  }
	  //	}
          if(strcmp("Oxygen",tr->name) == 0 )
              sm->erdepflux_oxygen = sm->erdepflux[n];
    }

    sinterface_putgridz_sed( hmodel, c, sm->gridz_sed);
    sinterface_putcellz_sed( hmodel, c, sm->cellz_sed);
    sinterface_putdiagtracer_2d(hmodel, c,
	   sm->hripples, param->n_hripple,
	   sm->lripples, param->n_lripple,
	   sm->ustrcw_skin, param->n_ustrcw_skin,
	   sm->depth_sedi, param->n_depth_sedi,
	   sm->dzactive,  param->n_dzactive,
	   sm->erdepflux_total, param->n_erdepflux_total,
	   param->n_erdepflux_total_ac,
           sm->erdepflux_oxygen, param->n_erdepflux_oxygen, param->n_erdepflux_oxygen_ac);

//NMY 2018
// material fluxes across water and sediments
 for(n=0; n < param->ntrB; n++) {
   double *point = NULL;
     if (param->fluxsedimap_inst[n] > 0) {
         point = sinterface_getpointerBtracer(hmodel, n, c); //get pointer to 2D diag tracer
         m = param->fluxsedimap_inst[n]; // get the number of the corresponding 3D tracer (i.e. erdepflux[m])
        *point = sm->erdepflux[m];  // update value of the 2D diagn tracer
     } else {
         if (param->fluxsedimap_ac[n] > 0) {
             point = sinterface_getpointerBtracer(hmodel, n, c); //get pointer to 2D diag tracer
             m = param->fluxsedimap_ac[n]; // get number of the corresponding 3D tracer
            *point += sm->erdepflux[m]*param->dt;  // update value of the 2D diagn tracer
         }
    }
  }


  if (param->verbose_sed == 2)
    mass_balance(sediment, sm);
  /*
    sed_limits(sediment, c);
  */
}


/************************************************************/

/* Check that the sediment data lies within valid limits */
/* MH 07/2012: included instability handling. Moved the limit checking from */
/* sed2hd() into its own routine. */
void sed_limits(sediment_t *sediment, sed_column_t *sm, int c)
{
  void *hmodel = sediment->hmodel;
  int k, n;
  int tk,bk;
  sed_params_t *param = sediment->msparam;
  int col_index = sm->col_number-1;
  char etext[MAXSTRLEN];

  if (param->verbose_sed == 1) {
    /* some checks */
    for(n=0; n < param->ntr; n++) {
      sed_tracer_t *tr = &sediment->mstracers[n];
      for(k=sm->botk_wc;k<=sm->topk_wc;k++) {
	//	if (!tr->prmspatial)
	  if (tr->diagn < 1 && sm->tr_wc[n][k] < 0.) {
	    sedtag(LWARN,"sed:sed2hd:sed2hd"," Negative tracer output in wc n=%d ", n);
	    sprintf(etext, "sed:sed2hd:sed_limits: k=%d topk=%d botk=%d topz=%f botz=%f col_numb=%d nstep=%d\n",
		    k, sm->topk_wc, sm->botk_wc,  sm->topz_wc, sm->botz_wc, sm->col_number, param->nstep);
	    i_set_error(hmodel, sm->col_number, LFATAL, etext);
	    diagnostics(sediment,sm," sed:sed2hd:sed2hd: Negative tracer output in wc"); 
	  }
        if (tr->diagn < 1 && finite(sm->tr_wc[n][k])==0 ) {
          sedtag(LWARN,"sed:sed2hd:sed2hd"," NaN or INF tracer output in wc n=%d ", n);
	    sprintf(etext, "sed:sed2hd:sed_limits: k=%d topk=%d botk=%d topz=%f botz=%f col_numb=%d nstep=%d\n",
		    k, sm->topk_wc, sm->botk_wc,  sm->topz_wc, sm->botz_wc, sm->col_number, param->nstep);
	    i_set_error(hmodel, sm->col_number, LFATAL, etext);
	    diagnostics(sediment, sm, " sed:sed2hd:sed2hd:  NaN or INF tracer output in wc"); 
	}
      } /* end for k=sm->botk_wc */
      for(k=sm->botk_sed;k<=sm->topk_sed;k++) {
	//	if (!tr->prmspatial)
	  if (tr->diagn < 1 && sm->tr_sed[n][k] < 0.) {
	    sedtag(LWARN,"sed:sed2hd:sed2hd"," Negative tracer output in sed n=%d \n", n);
	    sprintf(etext, "sed:sed2hd:sed_limits: k=%d topk=%d botk=%d topz=%f botz=%f col_numb=%d nstep=%d tr_sed=%e \n",
		   k, sm->topk_sed, sm->botk_sed,  sm->topz_sed,
		   sm->botz_sed, sm->col_number, param->nstep, sm->tr_sed[n][k]);
	    i_set_error(hmodel, sm->col_number, LFATAL, etext);
	    diagnostics(sediment, sm, " sed:sed2hd:sed2hd: Negative tracer output in sed"); 
	  }
	if (tr->diagn < 1 && finite(sm->tr_sed[n][k])==0 ) {
          sedtag(LWARN,"sed:sed2hd:sed2hd"," NaN or INF tracer output in sed n=%d \n", n);
	  sprintf(etext, "sed:sed2hd:sed_limits: k=%d topk=%d botk=%d topz=%f botz=%f col_numb=%d nstep=%d tr_sed=%e \n",
		 k, sm->topk_sed, sm->botk_sed,  sm->topz_sed,
		 sm->botz_sed, sm->col_number, param->nstep, sm->tr_sed[n][k]);
	  i_set_error(hmodel, sm->col_number, LFATAL, etext);
	  diagnostics(sediment, sm, " sed:sed2hd:sed2hd: NaN or INF tracer output in sed"); 
        }
      } /*end for k=sm->botk_sed */
    }
  }/* end if verbose*/
}
  /* END MH */

static void sed2hd_internal(sediment_t *sediment, sed_column_t *sm, int c)
{
  sed_params_t *param = sediment->msparam;
  sed_spatial_t *spatial = sediment->spatial;

  int n;
  spatial->hripples[c] = sm->hripples;
  spatial->lripples[c] = sm->lripples;
  for (n = 0; n < param->ntr; n++) {
      spatial->erdeprate[n][c] = sm->erdepflux[n];
      spatial->cbnm1[n][c] = sm->cbnm1[n];
      spatial->cbnm2[n][c] = sm->cbnm2[n];
      spatial->cbnm3[n][c] = sm->cbnm3[n];
      spatial->cbnm4[n][c] = sm->cbnm4[n];
  }

}
#endif


#ifdef __cplusplus
}
#endif
