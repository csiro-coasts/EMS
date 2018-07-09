/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/sediments/hd2sed.c
 *  
 *  Description:
 *  Routine to handle data flow from hydromodule to sediment
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: hd2sed.c 5848 2018-06-29 05:01:15Z riz008 $
 *
 */

#ifdef __cplusplus
#include "stdafx.h"
extern "C" {
#endif

#include "sediments.h"
#if defined(HAVE_SEDIMENT_MODULE)

extern void update_porosity_wc(sediment_t *sediment, sed_column_t *sm);
extern void update_porosity_sed(sediment_t *sediment, sed_column_t *sm);
static void hd2sed_internal(sediment_t *sediment, sed_column_t *sm, int c);

void reef_scale_depth(sediment_t *sediment, sed_column_t *sm);

double sinterface_getmodeltime(void* hmodel) ;
double sinterface_getcellarea(void* hmodel, int c);
double sinterface_getwind1(void* hmodel, int c) ;
double sinterface_getwind2(void* hmodel, int c) ;
double sinterface_getwindspeed(void* hmodel, int c);
double sinterface_gethripples(void* hmodel, int c) ;
double sinterface_getlripples(void* hmodel, int c) ;
double sinterface_getwave_ub(void* hmodel, int c) ;
double sinterface_getwave_period(void* hmodel, int c);
double sinterface_getwave_dir(void* hmodel, int c) ;
//double sinterface_gettheta(void* hmodel, int c);
int sinterface_gettopk_wc(void* hmodel, int c) ;
int sinterface_getbotk_wc(void* hmodel, int c) ;
double sinterface_getbotz_wc(void* hmodel, int c);
double sinterface_gettopz_wc(void* hmodel, int c) ;

int sinterface_getcelli(void* hmodel, int c);
int sinterface_getcellj(void* hmodel, int c);

void sinterface_getdz_wc(void* hmodel, int c, double *dz_wc) ;
void sinterface_getgridz_sed(void* hmodel, int c,double *gridz_sed );
void sinterface_getu1_wc(void* hmodel, int c, double *u1_wc);
void sinterface_getu2_wc(void* hmodel, int c, double *u2_wc);
void sinterface_getkz_wc(void* hmodel, int c, double *Kz_wc) ;
void si_gettracer_wc(void* hmodel, int c,
         int n, int m, double ***ptr_wc);
void si_gettracer_sed(void* hmodel, int c,
          int n, int m, double ***ptr_sed);
void diagnostics(sediment_t *sediment, sed_column_t *sm, const char *comment);
void mass_balance(sediment_t *sediment, sed_column_t *sm);

 double sinterface_getvalueofBtracer(void* hmodel, int n, int c);
double sinterface_get_srf_flux(void* hmodel, char *name, int c);

/* Fill sed_column_t structure */
void hd2sed(sediment_t *sediment, sed_column_t *sm, int c)
{
  int kone, k, n,m;
  void *hmodel = sediment->hmodel;
  sed_params_t *param = sediment->msparam;
  int col_index;
  FILE *flog;
  double a, rab1,rab2,ship_prob;
  char etext[MAXSTRLEN];

  /*UR-CHANGED , u_scale; */
  /* update column index */
  col_index=sm->col_number-1;

  sm->hd_col_index = c;
  sm->sed_start=1;

  /* Copy internal spatial sediment variables */
  hd2sed_internal(sediment, sm, col_index);
  // param->t=sinterface_getmodeltime(hmodel);
  sm->i=sinterface_getcelli(hmodel, c);
  sm->j=sinterface_getcellj(hmodel, c);
  sm->area = sinterface_getcellarea(hmodel, c);
  /* sm->z0_skin= const specified during initialisation ; */
  sm->wind1= sinterface_getwind1( hmodel, c);
  sm->wind2= sinterface_getwind2( hmodel, c) ;
  sm->windspeed= sinterface_getwindspeed( hmodel, c) ;
  sm->wave_ub = sinterface_getwave_ub( hmodel, c) ;
  sm->wave_period = sinterface_getwave_period( hmodel, c) ;
  sm->wave_dir = sinterface_getwave_dir( hmodel, c) ;
  //sm->theta = sinterface_gettheta( hmodel, c) ;
  sm->css_dep = param->css_dep;
  for(n=0; n < param->ntr; n++) {
      sm->erdepflux_ac[n] = 0.;
      sm->erdepflux[n] = 0.;
  }

  /* wc */
  sm->topk_wc =sinterface_gettopk_wc( hmodel, c) ;
  sm->botk_wc =sinterface_getbotk_wc( hmodel, c) ;
  sm->botz_wc =sinterface_getbotz_wc( hmodel, c);
  sm->topz_wc = sinterface_gettopz_wc( hmodel, c);
  sinterface_getdz_wc(hmodel, c, sm->dz_wc);

  sinterface_getu1_wc( hmodel, c, sm->u1_wc) ;
  sinterface_getu2_wc( hmodel, c, sm->u2_wc) ;
  sinterface_getkz_wc(hmodel, c, sm->Kz_wc);

  /*
  // test 2011
  if(sm->col_number == 1) {
    flog = fopen("uvtest.txt","a");
    fprintf(flog,"%f %f %f %f \n ", param->t,sm->u1_wc[sm->botk_wc],sm->u2_wc[sm->botk_wc],sm->Kz_wc[sm->botk_wc+1]);
    fclose(flog);
  }
  */

  sm->depth_wc = -sm->botz_wc + sm->topz_wc;

  /* get wc tracers */
  for(n=0; n < param->ntr; n++) {
      sed_tracer_t *tr = &sediment->mstracers[n];
      si_gettracer_wc(hmodel, c, tr->n, tr->n_hd_wc, sm->ptr_wc) ;
      // NMY 2013 Nov
      sm->tr_srf_flux[n] = sinterface_get_srf_flux(hmodel, tr->name, c);

  }

  for(n=0; n < param->ntr; n++) {
      sed_tracer_t *tr = &sediment->mstracers[n];
      sm->tr_srf_flux[n] /= tr->u_scale;
      for(k=sm->botk_wc;k<=sm->topk_wc;k++) {
	sm->tr_wc[n][k] = *sm->ptr_wc[n][k]/tr->u_scale;
        // NM: June 2017: tmp fix to get rid of negative concentrations
        // produced in EMS 
        if(sm->tr_wc[n][k] < 0) sm->tr_wc[n][k] = 0.;
      }     
  }

  /* Calculate sm->gridz_wc as it differs from window->gridz.
     (for example water surface and top of the window->gridz
     interface does not necessarily coincide, while in
     mecosed the grid starts from the surface) */
  if (sm->topk_wc > sm->botk_wc) {
    sm->gridz_wc[sm->topk_wc+1] = sm->topz_wc ;
    for(k=sm->topk_wc;k>=sm->botk_wc;k--)
      sm->gridz_wc[k] = sm->gridz_wc[k+1] - sm->dz_wc[k];
      //sm->gridz_wc[sm->botk_wc] = sm->botz_wc;
    for(k=sm->topk_wc;k>=sm->botk_wc;k--)
      sm->gridvel_wc[k] = 0.;
  }
  else
  { /* If there is one layer in water column split it on two layers.
       Note that only mecosed uses these two layers and
       when data are calculateted the layers are again merged
       into one layer and transferred to hydromodule */
    if (sm->topk_wc == sm->botk_wc) {
      kone=sm->topk_wc;
      sm->botk_wc = 0;
      sm->topk_wc = 1;
      sm->gridz_wc[2] = sm->topz_wc ;
      sm->gridz_wc[0] = sm->botz_wc;
      sm->gridz_wc[1] = 0.5*(sm->gridz_wc[0]+sm->gridz_wc[2]);
      for(k=0;k<2;k++) {
          sm->dz_wc[k] = sm->gridz_wc[k+1] - sm->gridz_wc[k] ;
          sm->u1_wc[k] = sm->u1_wc[kone];
          sm->u2_wc[k] = sm->u2_wc[kone];
          sm->Kz_wc[k] =   sm->Kz_wc[kone];
          sm->gridvel_wc[k] = 0.;
          for(n=0; n < param->ntr; n++)
	    sm->tr_wc[n][k] =sm->tr_wc[n][kone];
      }

      /* switch off srf fluxes for one-layer shallow water
// this is done in   sinterface_get_srf_flux routine
      if(sm->depth_wc < 0.01) {
	for(n=0; n < param->ntr; n++) {
	  sed_tracer_t *tr = &sediment->mstracers[n];
	  tr->srf_flux = 0;
	}
      }
      */

    }
  }


  for(k=sm->botk_wc;k<=sm->topk_wc;k++){
    //sm->dz_wc[k] = sm->gridz_wc[k+1] - sm->gridz_wc[k] ;
    sm->dzold_wc[k] = sm->dz_wc[k];
    sm->cellz_wc[k] = 0.5*(sm->gridz_wc[k]+ sm->gridz_wc[k+1]);
  }


  /* 
  {// 1d test ********************* 2010 May
    double zr=1;
    double ur=0;
    double tideT=0.5;
    double uA=1;
    ur = uA*cos(param->nstep*2.*PI/(tideT*24) +PI/2);
    ur += exp(-(param->nstep-5*24)*(param->nstep-5*24)/(2*24*2*24))+
          0.5*exp(-(param->nstep-10*24)*(param->nstep-10*24)/(1*24*1*24));
    for(k=sm->botk_wc;k<=sm->topk_wc;k++) {
      sm->u1_wc[k] = ur * log(-sm->gridz_wc[k]/sm->z0_skin)/
	log(zr/sm->z0_skin);
      sm->Kz_wc[k] =   0.001;
    }
    //  fprintf(stderr,"ur=%f u_top=%f \n",ur, sm->u1_wc[sm->topk_wc]); 
    //  fprintf(stderr,"z_bot=%f z_top=%f \n",sm->gridz_wc[sm->botk_wc], sm->gridz_wc[sm->topk_wc]); 
  }
*/

  /* sedibed */
  sm->botk_sed = 0;
  sm->topk_sed = param->sednz - 2;
    sinterface_getgridz_sed(hmodel, c, sm->gridz_sed);
  for(k=0;k<param->sednz;k++) {
    sm->dissol_kz[k] = (param->bi_dissol_kz_spv ? param->bi_dissol_kz_spv[col_index] : param->bi_dissol_kz);
    sm->partic_kz[k] = (param->bt_partic_kz_spv ? param->bt_partic_kz_spv[col_index] : param->bt_partic_kz);
    sm->css[k] = param->css * (param->css_scale_spv ? param->css_scale_spv[col_index] : param->css_scale);
    sm->gridvel_sed[k] = 0.;
  }
  sm->dissol_kz_i = (param->bi_dissol_kz_i_spv ? param->bi_dissol_kz_i_spv[col_index] : param->bi_dissol_kz_i);
  sm->partic_kz_i = (param->bt_partic_kz_i_spv ? param->bt_partic_kz_i_spv[col_index] : param->bt_partic_kz_i);


  /* Recalculate */
  sm->botz_sed = sm->gridz_sed[sm->botk_sed];
  sm->topz_sed = sm->gridz_sed[sm->topk_sed+1];
  sm->depth_sedi = - sm->botz_sed + sm->topz_sed;
  for(k=0;k<param->sednz-1;k++) {
    sm->dz_sed[k] = sm->gridz_sed[k+1] - sm->gridz_sed[k];
    sm->dzold_sed[k] = sm->dz_sed[k];
    sm->cellz_sed[k] = 0.5*(sm->gridz_sed[k]+ sm->gridz_sed[k+1]);
  }


  /*get sedibed tracers */
  for(n=0; n < param->ntr; n++) {
    sed_tracer_t *tr = &sediment->mstracers[n];
    si_gettracer_sed(hmodel, c, tr->n, tr->n_hd_sed, sm->ptr_sed) ;
  }
  for(n=0; n < param->ntr; n++) {
    sed_tracer_t *tr = &sediment->mstracers[n];
    for(k=sm->botk_sed;k<=sm->topk_sed;k++) {
     	sm->tr_sed[n][k] = *sm->ptr_sed[n][k]/tr->u_scale;
      	//NMY August 2017
	if (sm->tr_sed[n][k] < 0) sm->tr_sed[n][k] = 0.;
    }

  }

  update_porosity_wc(sediment, sm);
  update_porosity_sed(sediment, sm);
  
  //2010 css varying with sediment depth
  if (param->cssmode == 4) {
    // FR 06-2015: Note we are no longer setting param->css here
    double css = param->css_er_val[sm->topk_sed];
    for(k=sm->topk_sed;k>=sm->botk_sed;k--)
      if(sm->gridz_sed[sm->topk_sed+1] < param->css_er_depth[k])
	css=param->css_er_val[k];
    for(k=0;k<param->sednz;k++)
      sm->css[k] = css * (param->css_scale_spv ? param->css_scale_spv[col_index] : param->css_scale);
  }
  

  /* 2012 nov11 subregions: handle spatially varying parameters */
   if(param->nprmS)
     for (m=0;m<param->nprmS;m++){
       k=param->prmindexS[m]; 
       n=param->trindexS[m];    
       *param->prmpointS[k][col_index] = sinterface_getvalueofBtracer(hmodel, n, c);
       //  fprintf(stderr,"hd2sed:spatial prm %s; tarcer number %d, val %lf \n",
       //  param->prmname[k],n, *param->prmpoint[k] );
     }
  // end nov11

   // FR 06-2015: Moved from sed_init:alloc_sed_spatial
  {
    sed_spatial_t *spatial = sediment->spatial;
    double physriph = param->physriph;
    double physripl = param->physripl;
    if (param->physriph_spv != NULL)
      physriph = param->physriph_spv[col_index];
    if (param->physripl_spv != NULL)
      physripl = param->physripl_spv[col_index];

    spatial->hripples[col_index] = max(physriph, param->bioriph);
    spatial->lripples[col_index] = max(physripl, param->bioripl);
  }

 //2012 hardsub
   sm->hardsub_scale = 1.0;
   if(param->hardsub_numb > -1) {
     n = param->hardsub_numb;
     sm->hardsub_scale = 0.01 * fabs(100.0 - sinterface_getvalueofBtracer(hmodel, n, c) );

     //     if(sm->hardsub_scale > 1) 
     // fprintf(stderr, "hardsub_scale = %d, i=%d, j=%d \n",sm->hardsub_scale, sm->i, sm->j);
     //   if(sm->hardsub_scale < 0.1)sm->hardsub_scale = 0.1; //2012 tmp


   }

   //2015 ships
   if(param->ship_N) {
     ship_prob = ( param->ship_T * param->ship_N ) / (24 * param->ship_C);
     for(n=0;n<param->ship_C;n++) {
       if(param->ship_i[n]==sm->i && param->ship_j[n]==sm->j) {
         rab1= (1.*rand())/(RAND_MAX+1.);
	 if(rab1 <  ship_prob) {
	   for(k=sm->botk_wc;k<=sm->topk_wc+1;k++)
	     sm->Kz_wc[k] += param->ship_Kz;
	 }
       }
     }
   }


  /* Print section */
 if (param->verbose_sed == 2) {
   mass_balance(sediment, sm);
   if (sm->col_number==1 && param->nstep<2)
     diagnostics(sediment,sm,"Diagnostics print out for column 1 and step 1; verbose_sed=2");  
 }


  if (param->verbose_sed == 1) {
    if (sm->col_number==1 && param->nstep<2) {
      diagnostics(sediment,sm,"Diagnostics print out for column 1 and step 1; verbose_sed=1");  
    }
  }

/* some checks */
    if (sm->topk_wc < sm->botk_wc) {
      /* MH 07/2012: included instability handling */
      i_set_error(hmodel, sm->col_number, LFATAL, "sed:hd2sed: Error topk_wc < botk_wc \n");
      /*sedtag(LFATAL,"sed:hd2sed:hd2sed", "Error topk_wc < botk_wc \n");*/
      diagnostics(sediment,sm,"if(sm->topk_wc < sm->botk_wc)");  
      /*exit(1);*/
      /* END MH */
    }
    if (sm->topz_wc < sm->botz_wc) {
      /* MH 07/2012: included instability handling */
      i_set_error(hmodel, sm->col_number, LFATAL, "sed:hd2sed: topz_wc < botz_wc \n");
      /*sedtag(LFATAL,"sed:hd2sed:hd2sed", " topz_wc < botz_wc \n");*/
      diagnostics(sediment,sm, "if (sm->topz_wc < sm->botz_wc)");  
      /*exit(1);*/
      /* END MH */
    }
    if (fabs(sm->gridz_wc[sm->topk_wc+1] - sm->topz_wc) > 0.00001) {
      /* MH 07/2012: included instability handling */
      i_set_error(hmodel, sm->col_number, LFATAL, "sed:hd2sed: topz_wc =/= gridz_wc(topk_wc+1) \n");
      /*sedtag(LFATAL,"sed:hd2sed:hd2sed", " topz_wc =/= gridz_wc(topk_wc+1) \n");*/
      diagnostics(sediment,sm,"if (fabs(sm->gridz_wc[sm->topk_wc+1] - sm->topz_wc) > 0.00001)");  
      /*exit(1);*/
      /* END MH */
    }
    if (fabs(sm->gridz_wc[sm->botk_wc] - sm->botz_wc)> 0.00001) {
      /* MH 07/2012: included instability handling */
      i_set_error(hmodel, sm->col_number, LFATAL, "sed:hd2sed: Error botz_wc =/= gridz_wc(botz_wc) \n");
      /*sedtag(LFATAL,"sed:hd2sed:hd2sed", "Error botz_wc =/= gridz_wc(botz_wc) \n");*/
      diagnostics(sediment,sm," if (fabs(sm->gridz_wc[sm->botk_wc] - sm->botz_wc)> 0.00001)");
      /*exit(1);*/
      /* END MH */
    }

    a=0.;
    for(k=sm->botk_wc;k<=sm->topk_wc;k++)
      a=a+sm->dz_wc[k];

    if (fabs(a - sm->depth_wc) > 0.00001) {
      sedtag(LWARN,"sed:hd2sed:hd2sed", "Sum(dz[k]) =/= sm->depth_wc  col_numb= %d nstep=%d, c-index=%d, i=%d, j=%d \n",
	     sm->col_number, param->nstep, c, sm->i,sm->j);

      /* MH 07/2012: included instability handling */
      sprintf(etext, "sed:hd2sed: Sum(dz[k]) = %f, sm->depth_wc=%f \n", a, sm->depth_wc);
      i_set_error(hmodel, sm->col_number, LFATAL, etext);

      /*
      sprintf(etext, "sed:hd2sed:  botz = %f, topz = %f, sm->depth_wc=%f \n",  sm->botz_wc, sm->topz_wc, sm->depth_wc);
      i_set_error(hmodel, sm->col_number, LFATAL, etext);
      for(k=sm->botk_wc;k<=sm->topk_wc;k++) {
	sprintf(etext, "sed:hd2sed: k = %d, dz[k]=%f \n", k, sm->dz_wc[k]); 
	i_set_error(hmodel, sm->col_number, LFATAL, etext);
      }
      */

 
      /*sedtag(LFATAL,"sed:hd2sed:hd2sed", "Sum(dz[k]) = %f, sm->depth_wc=%f \n",
	a, sm->depth_wc);*/
      diagnostics(sediment,sm,"Sum(dz[k]) =/= sm->depth_wc");  
      /*exit(1);*/
      /* END MH */
    }

    for(k=sm->botk_wc;k<=sm->topk_wc;k++)
    {
      if ( sm->dz_wc[k] < 0.0) {
	/* MH 07/2012: included instability handling */
	sprintf(etext, "sed:hd2sed: dz_wc < 0, k=%d  \n", k);
	i_set_error(hmodel, sm->col_number, LFATAL, etext);
        /*sedtag(LFATAL,"sed:hd2sed:hd2sed", " dz_wc < 0, k=%d  \n", k);*/
	diagnostics(sediment,sm,"  if ( sm->dz_wc[k] < 0.0)");  
        /*exit(1);*/
	/* END MH */
      }
    }

    for(k=sm->botk_sed;k<=sm->topk_sed;k++) {
      if(sm->dz_sed[k] <= 0) {
	/*fprintf(stderr, "dz_sed<=0"); */
	/* MH 07/2012: included instability handling */
	i_set_error(hmodel, sm->col_number, LFATAL, "sed:hd2sed: dz_sed<=0");
	diagnostics(sediment,sm," if(sm->dz_sed[k] <= 0)");
	/*exit(1);*/
	/* END MH */
      }
      if(sm->dz_sed[k] > 1) {
	//	sm->css_dep = 0; // do not deposit
	//	fprintf(stderr, "dz_sed>100: t=%d i=%d j=%d k =%d hd_col_index=%d col_number=%d dz_sed=%f \n", param->nstep, sm->i, sm->j, k, sm->hd_col_index, sm->col_number, sm->dz_sed[k]);
	// diagnostics(sediment, "  if(sm->dz_sed[k] > 100)");  
	// exit(1);
      }

      if(sm->por_sed[k] < 0) {
	/* MH 07/2012: included instability handling */
	i_set_error(hmodel, sm->col_number, LFATAL, "sed:hd2sed: por_sed<0"); 
	/*fprintf(stderr, "por_sed<0"); */
	diagnostics(sediment,sm,"  if(sm->por_sed[k] < 0)");  
	/*exit(1);*/
	/* END MH */
      }
      if(sm->por_sed[k] >1) {
	/* MH 07/2012: included instability handling */
	i_set_error(hmodel, sm->col_number, LFATAL, "sed:hd2sed: por_sed>1");
	/*fprintf(stderr, "por_sed>1"); */
	diagnostics(sediment,sm," if(sm->por_sed[k] >1)");  
	/*exit(1);*/
	/* END MH */
      }
      if(sm->partic_kz[k] <= 0) {
	/* MH 07/2012: included instability handling */
	i_set_error(hmodel, sm->col_number, LFATAL, "sed:hd2sed: kz_sed<=0");
	/*fprintf(stderr, "kz_sed<=0"); */
	diagnostics(sediment,sm," if(sm->partic_kz[k] <= 0) ");  
	/*exit(1);*/
	/* END MH */
      }
    }

    for(n=0; n < param->ntr; n++) {
      sed_tracer_t *tr = &sediment->mstracers[n];
      for(k=sm->botk_wc;k<=sm->topk_wc;k++) {
	//	if (!tr->prmspatial)
        if (tr->diagn < 1 && sm->tr_wc[n][k] < 0.) {
          sedtag(LWARN,"sed:hd2sed:hd2sed",
          " Negative tracer input in wc n=%d \n", n);
	  /* MH 07/2012: included instability handling */
	  sprintf(etext, "sed:hd2sed: k=%d tiopk=%d botk=%d topz=%f botz=%f col_numb=%d nstep=%d \n",
		  k, sm->topk_wc, sm->botk_wc,  sm->topz_wc, sm->botz_wc, sm->col_number, param->nstep);
	  i_set_error(hmodel, sm->col_number, LFATAL, etext);
          /*sedtag(LFATAL,"sed:hd2sed:hd2sed",
          "k=%d tiopk=%d botk=%d topz=%f botz=%f col_numb=%d nstep=%d \n",
          k, sm->topk_wc, sm->botk_wc,  sm->topz_wc, sm->botz_wc, sm->col_number, param->nstep);*/
	  diagnostics(sediment,sm,"  if (tr->diagn < 1 && sm->tr_wc[n][k] < 0.)");  
          /*exit(1);*/
	  /* END MH */
        }
	if (tr->diagn < 1 && finite(sm->tr_wc[n][k])==0 ) {
          sedtag(LWARN,"sed:hd2sed:hd2sed",
          " NaN or INF tracer input into mecosed wc-layer n=%d \n", n);
	  /* MH 07/2012: included instability handling */
	  sprintf(etext, "sed:hd2sed: k=%d tiopk=%d botk=%d topz=%f botz=%f col_numb=%d nstep=%d \n",
		  k, sm->topk_wc, sm->botk_wc,  sm->topz_wc, sm->botz_wc, sm->col_number, param->nstep);
	  i_set_error(hmodel, sm->col_number, LFATAL, etext);
          /*sedtag(LFATAL,"sed:hd2sed:hd2sed",
          "k=%d tiopk=%d botk=%d topz=%f botz=%f col_numb=%d nstep=%d \n",
          k, sm->topk_wc, sm->botk_wc,  sm->topz_wc, sm->botz_wc, sm->col_number, param->nstep);*/
	  diagnostics(sediment,sm,"if (tr->diagn < 1 && finite(sm->tr_wc[n][k])==0 ) ");  
          /*exit(1);*/
	  /* END MH */
        }
      } /* end for k=sm->botk_wc */

      for(k=sm->botk_sed;k<=sm->topk_sed;k++) {
	//	if (!tr->prmspatial)
        if (tr->diagn < 1 && sm->tr_sed[n][k] < 0.) {
          sedtag(LWARN,"sed:hd2sed:hd2sed",
          "Negative tracer input in sed n=%d \n", n);
	  /* MH 07/2012: included instability handling */
	  sprintf(etext, "sed:hd2sed: k=%d tiopk=%d botk=%d topz=%f botz=%f col_numb= %d nstep=%d\n",
		  k, sm->topk_sed, sm->botk_sed,  sm->topz_sed, sm->botz_sed, sm->col_number, param->nstep);
	  i_set_error(hmodel, sm->col_number, LFATAL, etext);
          /*sedtag(LFATAL,"sed:hd2sed:hd2sed",
          "k=%d tiopk=%d botk=%d topz=%f botz=%f col_numb= %d nstep=%d\n",
          k, sm->topk_sed, sm->botk_sed,  sm->topz_sed,
          sm->botz_sed, sm->col_number, param->nstep);*/
	  diagnostics(sediment,sm," if (tr->diagn < 1 && sm->tr_sed[n][k] < 0.) ");  
          /*exit(1);*/
	  /* END MH */
        }
	if (tr->diagn < 1 && finite(sm->tr_sed[n][k])==0 ) {
          sedtag(LWARN,"sed:hd2sed:hd2sed",
          "NaN or INF tracer input into mecosed: sed-layer n=%d \n", n);
	  /* MH 07/2012: included instability handling */
	  sprintf(etext, "sed:hd2sed: k=%d tiopk=%d botk=%d topz=%f botz=%f col_numb= %d nstep=%d\n",
		  k, sm->topk_sed, sm->botk_sed,  sm->topz_sed, sm->botz_sed, sm->col_number, param->nstep);
	  i_set_error(hmodel, sm->col_number, LFATAL, etext);
          /*sedtag(LFATAL,"sed:hd2sed:hd2sed",
          "k=%d tiopk=%d botk=%d topz=%f botz=%f col_numb= %d nstep=%d\n",
          k, sm->topk_sed, sm->botk_sed,  sm->topz_sed,
          sm->botz_sed, sm->col_number, param->nstep);*/
	  diagnostics(sediment,sm, "if (tr->diagn < 1 && finite(sm->tr_sed[n][k])==0 ) ");  
          /*exit(1);*/
	  /* END MH */
        }
      } /* end for k=sm->botk_sed */
    }
}

/***************************************************************/
//  spatial->cbnm1[n][c] holds concentration of a tracer n in the 
// top sediment layer integrated over a number of previous time steps 
// to get rid of hf oscillations when eroding sediments

// NMY Feb 2015 - hripples amd lripples are redefined in bbl.c !
static void hd2sed_internal(sediment_t *sediment, sed_column_t *sm, int c)
{
  int n;
  sed_params_t *param = sediment->msparam;
  sed_spatial_t *spatial = sediment->spatial;

  sm->theta=spatial->theta[c];
  sm->hripples = spatial->hripples[c];
  sm->lripples = spatial->lripples[c];
  for (n = 0; n < param->ntr; n++) {
      sm->cbnm1[n] = spatial->cbnm1[n][c];
      sm->cbnm2[n] = spatial->cbnm2[n][c];
      sm->cbnm3[n] = spatial->cbnm3[n][c];
      sm->cbnm4[n] = spatial->cbnm4[n][c];
  }
}

/***************************************************************/

void diagnostics(sediment_t *sediment, sed_column_t *sm, const char *comment)
{
  int kone, k, n;
  sed_params_t *param = sediment->msparam;
  int col_index;
  FILE *flog;
  double a;
 /* Print section */

#ifdef HAVE_OMP
#pragma omp critical
  {
#endif
    if (param->verbose_sed) {
      flog = fopen("sedlog.txt","a");

      fprintf(flog, " \n hd2sed \n \n");
      fprintf(flog, " %s \n",comment);
      fprintf(flog, "sm->i=%d, sm->j=%d, sm->hd_col_index=%d \n"
        ,sm->i,sm->j,sm->hd_col_index);
      fprintf(flog, "sm->hd_col_index=%d, sm->col_number=%d, param->nstep=%d \n"
	      ,sm->hd_col_index,sm->col_number,param->nstep);
      fprintf(flog, "sm->botk_wc=%d, sm->topk_wc=%d \n"
        ,sm->botk_wc,sm->topk_wc);
      fprintf(flog, "sm->botz_wc=%f, sm->topz_wc=%f \n \n"
        ,sm->botz_wc,sm->topz_wc);
      for(k=sm->botk_wc;k<=sm->topk_wc+1;k++)
        fprintf(flog, "k=%d, sm->gridz_wc=%f,  sm->Kz=%f \n",
            k, sm->gridz_wc[k], sm->Kz_wc[k]);
        fprintf(flog, "\n");
      for(k=sm->botk_wc;k<=sm->topk_wc;k++)
        fprintf(flog, "k=%d, sm->dz_wc=%f, sm->cellz_wc=%f \n",
            k, sm->dz_wc[k], sm->cellz_wc[k]);
      fprintf(flog, "\n");
      for(k=sm->botk_wc;k<=sm->topk_wc;k++)
        fprintf(flog, "k=%d, sm->u1_wc=%f, sm->u2_wc=%f, \n",
            k ,sm->u1_wc[k],sm->u2_wc[k]);
      fprintf(flog, "\n");

      for(n=0; n < param->ntr; n++) {
        for(k=sm->botk_wc;k<=sm->topk_wc;k++) {
          fprintf(flog, "n=%d, k=%d, sm->tr_wc[n][k] = %f \n",
          n, k, sm->tr_wc[n][k]);
        }
        fprintf(flog, "\n");
      }

      fprintf(flog, " \n ");

      fprintf(flog, "sm->botk_sed=%d, sm->topk_sed=%d \n"
        ,sm->botk_sed,sm->topk_sed);
      fprintf(flog, "sm->botz_sed=%f, sm->topz_sed=%f \n"
        ,sm->botz_sed,sm->topz_sed);
      fprintf(flog, "\n");
      for(k=sm->botk_sed;k<=sm->topk_sed+1;k++)
        fprintf(flog, "k=%d, sm->gridz_sed=%f, sm->dissol_kz=%e, sm->partic_kz=%e\n", k, sm->gridz_sed[k], sm->dissol_kz[k], sm->partic_kz[k]);
      fprintf(flog, "\n");
      for(k=sm->botk_sed;k<=sm->topk_sed;k++)
        fprintf(flog, "k=%d, sm->dz_sed=%f, sm->cellz_sed=%f, sm->por_sed=%f \n", k, sm->dz_sed[k], sm->cellz_sed[k] ,sm->por_sed[k]);
      fprintf(flog, "\n");

      for(n=0; n < param->ntr; n++) {
        for(k=sm->botk_sed;k<=sm->topk_sed;k++) {
          fprintf(flog, "n=%d, k=%d, sm->tr_sed[n][k] = %f \n",
          n, k, sm->tr_sed[n][k]);
        }
        fprintf(flog, "\n");
      }
      fprintf(flog, "\n");
      fprintf(flog, "sm->depth_wc=%f, sm->depth_sedi=%f \n"
        ,sm->depth_wc,sm->depth_sedi);

      fclose(flog);
    }
#ifdef HAVE_OMP
  }
#endif

} //end diagnostics

/*****************/
void mass_balance(sediment_t *sediment, sed_column_t *sm)
{
  int k, n;
  sed_params_t *param = sediment->msparam;
  FILE *flog;
 /* Print section */

      // ********* 2010 sediment mass balance check before and after 
      // the sediment cycle
#ifdef HAVE_OMP
#pragma omp critical
  {
#endif
      // initialise
      if(sm->col_number == 1 && sm->sed_start == 1) {
	for(n=0; n < param->ntr; n++) {
	  sed_tracer_t *tr = &sediment->mstracers[n];
	  tr->mass_start_wc = 0;
	  tr->mass_start_sed = 0;
	  tr->mass_end_wc = 0;
	  tr->mass_end_sed = 0;
	}
      }
      // mass at the beiginning of the sed cycle
      if(sm->sed_start == 1) {
	for(n=0; n < param->ntr; n++) {
	  sed_tracer_t *tr = &sediment->mstracers[n];
	  for(k=sm->botk_wc;k<=sm->topk_wc;k++) {
	    if(tr->dissol)
	      tr->mass_start_wc += sm->tr_wc[n][k]*sm->dz_wc[k]*sm->por_wc[k]*sm->area/1e6;
	    else
	      tr->mass_start_wc += sm->tr_wc[n][k]*sm->dz_wc[k]*sm->area/1e6;
	  }
	  for(k=sm->botk_sed;k<=sm->topk_sed;k++) {
	    if(tr->dissol)
	      tr->mass_start_sed += sm->tr_sed[n][k]*sm->dz_sed[k]*sm->por_sed[k]*sm->area/1e6;
	    else
	      tr->mass_start_sed += sm->tr_sed[n][k]*sm->dz_sed[k]*sm->area/1e6;
	  }
	}
      }
      // mass after the sed cycle
      if(sm->sed_start == 0) {
	for(n=0; n < param->ntr; n++) {
	  sed_tracer_t *tr = &sediment->mstracers[n];
	  for(k=sm->botk_wc;k<=sm->topk_wc;k++) {
	    if(tr->dissol)
	      tr->mass_end_wc += sm->tr_wc[n][k]*sm->dz_wc[k]*sm->por_wc[k]*sm->area/1e6;
	    else
	      tr->mass_end_wc += sm->tr_wc[n][k]*sm->dz_wc[k]*sm->area/1e6;
	  }
	  for(k=sm->botk_sed;k<=sm->topk_sed;k++) {
	    if(tr->dissol)
	      tr->mass_end_sed += sm->tr_sed[n][k]*sm->dz_sed[k]*sm->por_sed[k]*sm->area/1e6;
	    else
	      tr->mass_end_sed += sm->tr_sed[n][k]*sm->dz_sed[k]*sm->area/1e6;
	  }
	}
      }
      // print the mass balance
      if(sm->col_number == sediment->msparam->ncol && sm->sed_start == 0) {
	flog = fopen("sed_mass_start.txt","a");
	fprintf(flog,"%f ", param->t);
	for(n=0; n < param->ntr; n++) {
	  sed_tracer_t *tr = &sediment->mstracers[n];
	  fprintf(flog, "%f ", tr->mass_start_wc +  tr->mass_start_sed);
	}
	fprintf(flog, "\n");
	fclose(flog);

	/**/
	flog = fopen("sed_mass_end.txt","a");
	fprintf(flog,"%f ", param->t);
	for(n=0; n < param->ntr; n++) {
	  sed_tracer_t *tr = &sediment->mstracers[n];
	  fprintf(flog, "%f ", tr->mass_end_wc +  tr->mass_end_sed);
	}
	fprintf(flog, "\n");
	fclose(flog);
	/**/

	flog = fopen("sed_mass_diff.txt","a");
	fprintf(flog,"%f ", param->t);
	for(n=0; n < param->ntr; n++) {
	  sed_tracer_t *tr = &sediment->mstracers[n];
	  fprintf(flog, "%e ", tr->mass_end_wc +  tr->mass_end_sed - tr->mass_start_wc -  tr->mass_start_sed);
	}
	fprintf(flog, "\n");
	fclose(flog);
      }

      //*************
#ifdef HAVE_OMP
  }
#endif

} //end mass_balance


/************************/
// scales ripple height and benthic diff coefficients with depth
// implemented in bbl.c to scale roughness and in vtransp.c to scale kz and kz_i
// invoked with REEF_SCALE_DEPTH key
void reef_scale_depth(sediment_t *sediment, sed_column_t *sm)
{
  sed_params_t *param = sediment->msparam;
  int k;
  double rab1,rab2, mindep,maxdep;

  mindep=20.;
  maxdep=param->reef_scale_depth;

  if(maxdep < 1e-11) 
    return;

 /* eReefs scaling with water depth NMY 2013 */
 if (sm->depth_wc < mindep) { 
   rab1=1;
   rab2=1;
 }
 else
   if (sm->depth_wc>=mindep && sm->depth_wc<maxdep ) {
     rab1=1-(sm->depth_wc - mindep)/(maxdep-mindep);
     rab2=1+3*(sm->depth_wc - mindep)/(maxdep-mindep);
   }
   else {
     rab1=0.;
     rab2=4.;
   }

  for(k=0;k<param->sednz;k++) {
       sm->dissol_kz[k] = sm->dissol_kz[k]*rab1+1e-25;
       sm->partic_kz[k] = sm->partic_kz[k]*rab1+1e-25;
       sm->css[k] = sm->css[k]*rab2;
  }
 sm->dissol_kz_i = sm->dissol_kz_i*(rab1+0.1);
 sm->partic_kz_i = sm->partic_kz_i*(rab1+0.1);
 // sm->hripples = param->physriph * (rab1+0.1); 

 /*  
  if(sm->i == 361 && sm->j==51)
    fprintf(stderr,"hd2sed depth_wc= %lf, css=%lf, kz_i=%e, hripples=%lf \n", sm->depth_wc,sm->css[0], sm->partic_kz_i,sm->hripples);
 */

}
 /* end eReef scaling */

#endif



#ifdef __cplusplus
}
#endif

