/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/inputs/readparam_r.c
 *  
 *  Description:
 *  Routine to overwrite ROAM model parameters
 *  from a file.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: readparam_r.c 6757 2021-04-07 01:00:44Z her127 $
 *
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hd.h"

#define LAYERMIN   0.5

/*-------------------------------------------------------------------*/
/* Routine to overwrite input parameters to ROAM values              */
/*-------------------------------------------------------------------*/
void auto_params_roam_pre1(FILE * fp, parameters_t *params)
{
  int n;                        /* Counters */
  open_bdrys_t *open;           /* Pointer to boundary structure */
  char keyword[MAXSTRLEN];      /* Input dummy */
  char buf[MAXSTRLEN];          /* Input dummy */

  /* Parameter header (optional) */
  sprintf(keyword, "PARAMETERHEADER");
  if (!(prm_read_char(fp, keyword, params->parameterheader))) {
    strcpy(params->parameterheader, "ROAMv1 grid");
  }
  strcpy(parameterheader, params->parameterheader);

  /*-----------------------------------------------------------------*/
  /* ROAM forcing data                                               */
  /* OFAM_DATA : used for obc forcing and IC's                       */
  /* INIT_DATA : if present used for IC instead of OFAM_DATA         */
  /* TEMP_DATA : if present used for T IC instead of INIT_DATA       */
  /* SALT_DATA : if present used for S IC instead of INIT_DATA       */
  /* ETA_DATA : if present used for eta IC instead of INIT_DATA      */
  char *fields[MAXSTRLEN * MAXNUMARGS];
  int nf;
  sprintf(params->rdata, "%c", '\0');
  sprintf(params->odata, "%c", '\0');
  sprintf(params->idata, "%c", '\0');
  sprintf(params->tdata, "%c", '\0');
  sprintf(params->sdata, "%c", '\0');
  sprintf(params->edata, "%c", '\0');
  sprintf(params->vdata, "%c", '\0');
  prm_read_char(fp, "RAMS_DATA", params->rdata);
  prm_read_char(fp, "OFAM_DATA", params->odata);
  strcpy(buf, params->odata);
  nf = parseline(buf, fields, MAXNUMARGS);
  if (!(prm_read_char(fp, "INIT_DATA", params->idata))) {
    if (nf >= 1)
      strcpy(params->idata, fields[0]);
  }
  if(!(prm_read_char(fp, "TEMP_DATA", params->tdata)))
    strcpy(params->tdata, params->idata);
  if(!(prm_read_char(fp, "SALT_DATA", params->sdata)))
    strcpy(params->sdata, params->idata);
  if(!(prm_read_char(fp, "ETA_DATA", params->edata)))
    strcpy(params->edata, params->idata);
  if(!(prm_read_char(fp, "VELOCITY_DATA", params->vdata)))
    strcpy(params->vdata, params->idata);
  if(!strlen(params->tdata) || !strlen(params->sdata) ||
     !strlen(params->edata))
    hd_quit("auto_params: Insufficient initialisation data specified\n");
  if(!strlen(params->odata))
    sprintf(params->odata,"%s %s %s",params->edata, params->tdata,
	    params->sdata);
  params->rampf = WIND|TIDALH|INV_BARO;
  params->save_force = OTEMP|OSALT|OETA;
  params->save_force = NONE;
  params->save_force = OTEMP|OSALT;
  params->save_force = OTEMP|OSALT|ROAM;
  params->save_force = OTEMP|OSALT|OVELU|OVELV|ROAM;
  params->save_force = OVELU|OVELV|ROAM;
  params->atr += 2;    /* ovelu and ovelv */

  params->rsalt = params->rtemp = 1;
  params->atr += 2;    /* rtemp and rsalt */

  params->robust = 6;
  params->speed = 5;
  params->compatible |= V1670;

  /* Mixing scheme (optional) */
  sprintf(keyword, "MIXING_SCHEME");
  if (!prm_read_char(fp, keyword, params->mixsc)) {
    strcpy(params->mixsc, "mellor_yamada_2_0");
    params->vz0 = 1e-5;
    params->kz0 = 1e-5;
    params->zs = 1.0;
    params->atr -= 2; /* Remove default k-e auto tracers */
  }
  /* Surface height */
  if (strlen(params->edata))
    strcpy(params->eta_init, params->edata);

  /* Wind and pressure */
  if (strlen(params->rdata)) {
    if (strlen(params->wind) == 0) {
      strcpy(params->wind, params->rdata);
      params->wind_dt = 600.0;
    }
    if (strlen(params->patm) == 0) {
      strcpy(params->patm, params->rdata);
      params->patm_dt = 600.0;
    }
  }

  /* Heatflux */
  sprintf(keyword, "HEATFLUX");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "COMP_HEAT_NONE") == 0) {
      params->heatflux = COMP_HEAT_NONE;
      sprintf(params->hf, "%c", '\0');
      prm_set_errfn(hd_quit);
      sprintf(keyword, "HEATFLUX_FILE");
      prm_read_char(fp, keyword, params->hf);
      sprintf(keyword, "HEATFLUX_DT");
      prm_get_time_in_secs(fp, keyword, &params->hf_dt);
      params->ntrS += 4+2; // 4 x hf comps + precip and evap
    } else {
      params->smooth_VzKz = 1;
      /* Set defaults */
      /*params->albedo = 0.0;*/
      strcpy(params->swr_attn,"0.2");
      strcpy(params->swr_tran,"0.42");
      params->water_type = TYPE_II;
      params->hf_ramp = -43200.0;
      /* Override as required */
      read_swr(params, fp, 0);
      sprintf(keyword, "ALBEDO");
      prm_read_double(fp, keyword, &params->albedo);
      
      if (strcmp(buf, "COMP_HEAT") == 0 && strlen(params->rdata)) {
	params->heatflux = COMP_HEAT;
	strcpy(params->hf, params->rdata);
	params->hf_dt = 600.0;
	strcpy(params->swr, params->rdata);
	params->swr_dt = 600.0;
	strcpy(params->rh, params->rdata);
	params->rh_dt = 600.0;
	params->ntrS += 5;
      }
      if ((strcmp(buf, "ADVANCED") == 0 || strcmp(buf, "BULK") == 0)
	  && strlen(params->rdata)) {
	params->heatflux = ADVANCED;
	n = 3;
	if (params->water_type != NONE) n -= 2;
	if (prm_skip_to_end_of_key(fp, "SWR_ATTENUATION")) n -= 1;
	if (prm_skip_to_end_of_key(fp, "SWR_TRANSMISSION")) n -= 1;
	if (prm_skip_to_end_of_key(fp, "SWR_BOT_ABSORB")) n -= 1;
	read_swr(params, fp, 0);
	params->ntrS += n;
	read_hf_bulk(params, fp);
	strcpy(params->airtemp, params->rdata);
	params->airtemp_dt = 600.0;
	strcpy(params->wetb, params->rdata);
	params->wetb_dt = 600.0;
	strcpy(params->cloud, params->rdata);
	params->cloud_dt = 600.0;
	/*
	  strcpy(params->swr, params->rdata);
	  params->swr_dt = 600.0;
	*/
	params->ntrS += 5;
      }
    }
  }

  /*
   * Add dummy string to force writing of transport file
   */
  if (params->runmode & PRE_MARVL)
    sprintf(params->trans_data, "PRE_MARVL");

  /*
   * Read precipitation and evaporation, if specified
   */
  sprintf(params->precip, "%c", '\0');
  sprintf(keyword, "PRECIPITATION");
  if (prm_read_char(fp, keyword, params->precip)) {
    sprintf(keyword, "PRECIPITATION_INPUT_DT");
    prm_get_time_in_secs(fp, keyword, &params->precip_dt);
  }
  sprintf(params->evap, "%c", '\0');
  sprintf(keyword, "EVAPORATION");
  if (prm_read_char(fp, keyword, params->evap)) {
    sprintf(keyword, "EVAPORATION_INPUT_DT");
    prm_get_time_in_secs(fp, keyword, &params->evap_dt);
  }

  /* ROAM speed (optional) */
  sprintf(keyword, "SPEED");
  prm_read_int(fp, keyword, &params->speed);


  /* ROAM robustness (optional) */
  sprintf(keyword, "ROBUST");
  prm_read_int(fp, keyword, &params->robust);
  if (params->robust == 1) 
    sprintf(params->robusttext, "Smagorinsky viscosity/diffusivity = 0.1");
  if (params->robust == 2) 
    sprintf(params->robusttext, "Smagorinsky viscosity/diffusivity = 0.2");
  if (params->robust == 3) 
    sprintf(params->robusttext, "Smagorinsky viscosity/diffusivity = 0.3");
  if (params->robust == 4) 
    sprintf(params->robusttext, "Smagorinsky viscosity/diffusivity = 0.4");
  if (params->robust == 5) 
    sprintf(params->robusttext, "Smagorinsky viscosity/diffusivity = 0.5");
  if (params->robust == 6) 
    sprintf(params->robusttext, "Smagorinsky diffusivity = 0.1");
  if (params->robust == 7) 
    sprintf(params->robusttext, "Smagorinsky diffusivity = 0.1, start from rest, time-step *= 0.875");
  if (params->robust == 8) 
    sprintf(params->robusttext, "Constant horizontal viscosity/diffusivity, start from rest, time-step *= 0.75");
  if (params->robust == 9) 
    sprintf(params->robusttext, "Constant horizontal viscosity/diffusivity, start from rest, time-step *= 0.625");
  if (params->robust == 10) 
    sprintf(params->robusttext, "Constant horizontal viscosity/diffusivity, start from rest, time-step *= 0.5");

  if ((params->robust >= 1 && params->robust <= 7) && 
      params->smagorinsky == 0.0)
    params->atr++;
  if (params->robust < 1)
    params->robust = 1;
  if (params->robust > 10)
    params->robust = 10;
  if (params->robust == 1) {
      params->smagorinsky = 0.1;
  } else if (params->robust == 2) {
    params->smagorinsky = 0.2;
    params->smag_smooth = 1;
  } else if (params->robust == 3) {
    params->smagorinsky = 0.3;
    params->smag_smooth = 2;
  } else if (params->robust == 4) {
    params->smagorinsky = 0.4;
    params->smag_smooth = 3;
  } else if (params->robust == 5) {
    params->smagorinsky = 0.5;
    params->smag_smooth = 4;
  } else if (params->robust == 6 || params->robust == 7) {
    params->smagorinsky = 0.1;
    params->smag_smooth = 1;
  }  else
    params->smagorinsky = 0.0;
  params->sue1 = params->kue1 = 1.0;
  params->bsue1 = params->bkue1 = 0.0;
  sprintf(params->smag, "%f", params->smagorinsky);

  /* Velocity */
  if (params->robust <= 6 && strlen(params->vdata))
    strcpy(params->vel_init, params->vdata);

  /* Thresholds */
  params->velmax = 3.87;
  params->velmax2d = 2.55;
  if (!params->etadiff)
    params->etadiff = 0.60;
  if (contains_token(params->alert, "PASSIVE") == NULL &&
      contains_token(params->alert, "ACTIVE") == NULL) {
    strcpy(params->alert, "ACTIVE");
    params->ntrS += 4;
  }

  /* Means (optional) */
  if (!(params->means & ETA_M)) {
    if (params->means == NONE)
      params->means = ETA_M;
    else
      params->means |= ETA_M;
    params->ntrS += 1;
    sprintf(params->means_dt, "%f days", get_run_length(fp) + 1);
  }
  params->tide_r = MEAN_R;

  /* Bathymetry checks (optional) */
  sprintf(keyword, "MAXGRAD");
  if (!(prm_read_double(fp, keyword, &params->maxgrad)))
    params->maxgrad = 0.05;
  sprintf(keyword, "SMOOTHING");
  if (!(prm_read_int(fp, keyword, &params->smooth)))
    params->smooth = 1;

  /* Horizontal mixing */
  params->u1vh = params->u2vh = params->u1kh = params->u2kh = -1.0;
  if (params->robust == 6 || params->robust == 7) params->u1vh = params->u2vh = 1.0;

  /* Surface relaxation (optional) */
  if (strlen(params->edata))
    strcpy(params->etarlxn, params->edata);
  if (prm_read_char(fp, "eta_relaxation_input_dt", buf) > 0 &&
      prm_read_char(fp, "eta_relaxation_time_constant", keyword) > 0) {
    tm_scale_to_secs(buf, &params->etarlxdt);
    tm_scale_to_secs(keyword, &params->etarlxtc);
    strcpy(params->etarlxtcs, keyword);
    params->etarlx = RELAX;
  } else {
    params->etarlxdt = 3600.0;
    params->etarlxtc = 0.0;
    params->etarlx = ALERT;
  }
  params->ntrS++;

  /* netCDF dumpfile name */
  if (params->runmode & RE_ROAM) {
    sprintf(keyword, "RESTART_FILE");
    prm_read_char(fp, keyword, params->idumpname);
  }

  sprintf(codeheader, "SHOC ROAMv1 : ROBUST %d", params->robust);
}

/* END auto_params_roam_pre1()                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to overwrite input parameters to ROAM values.             */
/* Developed Nov 2013 to include updated ROBUST parameterisation.    */
/* ROBUST=1: OFAM currents + Smagorinsky = 0.1                       */
/* ROBUST=2: OFAM currents + Smagorinsky = 0.2                       */
/* ROBUST=3: Geostrophic currents + Smagorinsky = 0.1                */
/* ROBUST=4: Rest + hard T/S ramp relaxation + Smagorinsky = 0.1     */
/* ROBUST=5: Start from rest + Smagorinsky = 0.1                     */
/* ROBUST=6: Start from rest + Smagorinsky = 0.2                     */
/* ROBUST=7: OFAM currents + constant horizontal viscosity           */
/* ROBUST=8: Rest + hard T/S ramp relaxation + constant viscosity    */
/* ROBUST=9: Start from rest + constant horizontal viscosity         */
/* ROBUST=10: As for ROBUST=9 with reduced time-step                 */
/*-------------------------------------------------------------------*/
void auto_params_roam_pre2(FILE * fp, parameters_t *params)
{
  int n;                        /* Counters */
  open_bdrys_t *open;           /* Pointer to boundary structure */
  char keyword[MAXSTRLEN];      /* Input dummy */
  char buf[MAXSTRLEN];          /* Input dummy */

  /* Parameter header (optional) */
  sprintf(keyword, "PARAMETERHEADER");
  if (!(prm_read_char(fp, keyword, params->parameterheader))) {
    strcpy(params->parameterheader, "ROAMv2 grid");
  }
  strcpy(parameterheader, params->parameterheader);

  /*-----------------------------------------------------------------*/
  /* ROAM forcing data                                               */
  /* OFAM_DATA : used for obc forcing and IC's                       */
  /* INIT_DATA : if present used for IC instead of OFAM_DATA         */
  /* TEMP_DATA : if present used for T IC instead of INIT_DATA       */
  /* SALT_DATA : if present used for S IC instead of INIT_DATA       */
  /* ETA_DATA : if present used for eta IC instead of INIT_DATA      */
  char *fields[MAXSTRLEN * MAXNUMARGS];
  int nf;
  sprintf(params->rdata, "%c", '\0');
  sprintf(params->odata, "%c", '\0');
  sprintf(params->idata, "%c", '\0');
  sprintf(params->tdata, "%c", '\0');
  sprintf(params->sdata, "%c", '\0');
  sprintf(params->edata, "%c", '\0');
  sprintf(params->vdata, "%c", '\0');
  prm_read_char(fp, "RAMS_DATA", params->rdata);
  prm_read_char(fp, "OFAM_DATA", params->odata);
  strcpy(buf, params->odata);
  nf = parseline(buf, fields, MAXNUMARGS);
  if (!(prm_read_char(fp, "INIT_DATA", params->idata))) {
    if (nf >= 1)
      strcpy(params->idata, fields[0]);
  }
  if(!(prm_read_char(fp, "TEMP_DATA", params->tdata)))
    strcpy(params->tdata, params->idata);
  if(!(prm_read_char(fp, "SALT_DATA", params->sdata)))
    strcpy(params->sdata, params->idata);
  if(!(prm_read_char(fp, "ETA_DATA", params->edata)))
    strcpy(params->edata, params->idata);
  if(!(prm_read_char(fp, "VELOCITY_DATA", params->vdata)))
    strcpy(params->vdata, params->idata);
  if(!strlen(params->tdata) || !strlen(params->sdata) ||
     !strlen(params->edata))
    hd_quit("auto_params: Insufficient initialisation data specified\n");
  if(!strlen(params->odata))
    sprintf(params->odata,"%s %s %s",params->edata, params->tdata,
	    params->sdata);
  params->rampf = WIND|TIDALH|INV_BARO;
  params->save_force = OTEMP|OSALT|OETA;
  params->save_force = NONE;
  params->save_force = OTEMP|OSALT;
  params->save_force = OTEMP|OSALT|ROAM;
  params->save_force = OTEMP|OSALT|OVELU|OVELV|ROAM;
  params->save_force = OVELU|OVELV|ROAM;
  params->atr += 2;    /* ovelu and ovelv */

  params->rsalt = params->rtemp = 1;
  params->atr += 2;    /* rtemp and rsalt */

  params->robust = 7;
  params->speed = 5;
  params->compatible |= V1670;

  /* Mixing scheme (optional) */
  sprintf(keyword, "MIXING_SCHEME");
  if (!prm_read_char(fp, keyword, params->mixsc)) {
    strcpy(params->mixsc, "mellor_yamada_2_0");
    params->vz0 = 1e-5;
    params->kz0 = 1e-5;
    params->zs = 1.0;
    params->atr -= 2; /* Remove default k-e auto tracers */
  }
  /* Surface height */
  if (strlen(params->edata))
    strcpy(params->eta_init, params->edata);

  /* Wind and pressure */
  if (strlen(params->rdata)) {
    if (strlen(params->wind) == 0) {
      strcpy(params->wind, params->rdata);
      params->wind_dt = 600.0;
    }
    if (strlen(params->patm) == 0) {
      strcpy(params->patm, params->rdata);
      params->patm_dt = 600.0;
    }
  }

  /* Heatflux */
  sprintf(keyword, "HEATFLUX");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "COMP_HEAT_NONE") == 0) {
      params->heatflux = COMP_HEAT_NONE;
      sprintf(params->hf, "%c", '\0');
      prm_set_errfn(hd_quit);
      sprintf(keyword, "HEATFLUX_FILE");
      prm_read_char(fp, keyword, params->hf);
      sprintf(keyword, "HEATFLUX_DT");
      prm_get_time_in_secs(fp, keyword, &params->hf_dt);
      params->ntrS += 4+2; // 4 x hf comps + precip and evap
    } else {
      params->smooth_VzKz = 1;
      /* Set defaults */
      /*params->albedo = 0.0;*/
      strcpy(params->swr_attn,"0.2");
      strcpy(params->swr_tran,"0.42");
      params->water_type = TYPE_II;
      params->hf_ramp = -43200.0;
      /* Override as required */
      read_swr(params, fp, 0);
      sprintf(keyword, "ALBEDO");
      prm_read_double(fp, keyword, &params->albedo);
      
      if (strcmp(buf, "COMP_HEAT") == 0 && strlen(params->rdata)) {
	params->heatflux = COMP_HEAT;
	strcpy(params->hf, params->rdata);
	params->hf_dt = 600.0;
	strcpy(params->swr, params->rdata);
	params->swr_dt = 600.0;
	strcpy(params->rh, params->rdata);
	params->rh_dt = 600.0;
	params->ntrS += 5;
      }
      if ((strcmp(buf, "ADVANCED") == 0 || strcmp(buf, "BULK") == 0)
	  && strlen(params->rdata)) {
	params->heatflux = ADVANCED;
	n = 3;
	if (params->water_type != NONE) n -= 2;
	if (prm_skip_to_end_of_key(fp, "SWR_ATTENUATION")) n -= 1;
	if (prm_skip_to_end_of_key(fp, "SWR_TRANSMISSION")) n -= 1;
	if (prm_skip_to_end_of_key(fp, "SWR_BOT_ABSORB")) n -= 1;
	read_swr(params, fp, 0);
	params->ntrS += n;
	read_hf_bulk(params, fp);
	strcpy(params->airtemp, params->rdata);
	params->airtemp_dt = 600.0;
	strcpy(params->wetb, params->rdata);
	params->wetb_dt = 600.0;
	strcpy(params->cloud, params->rdata);
	params->cloud_dt = 600.0;
	/*
	  strcpy(params->swr, params->rdata);
	  params->swr_dt = 600.0;
	*/
	params->ntrS += 5;
      }
    }
  }

  /*
   * Add dummy string to force writing of transport file
   */
  if (params->runmode & PRE_MARVL)
    sprintf(params->trans_data, "PRE_MARVL");

  /*
   * Read precipitation and evaporation, if specified
   */
  sprintf(params->precip, "%c", '\0');
  sprintf(keyword, "PRECIPITATION");
  if (prm_read_char(fp, keyword, params->precip)) {
    sprintf(keyword, "PRECIPITATION_INPUT_DT");
    prm_get_time_in_secs(fp, keyword, &params->precip_dt);
  }
  sprintf(params->evap, "%c", '\0');
  sprintf(keyword, "EVAPORATION");
  if (prm_read_char(fp, keyword, params->evap)) {
    sprintf(keyword, "EVAPORATION_INPUT_DT");
    prm_get_time_in_secs(fp, keyword, &params->evap_dt);
  }

  /* ROAM speed (optional) */
  sprintf(keyword, "SPEED");
  prm_read_int(fp, keyword, &params->speed);

  /* ROAM robustness (optional) */
  sprintf(keyword, "ROBUST");
  prm_read_int(fp, keyword, &params->robust);
  if (params->robust == 0) 
    sprintf(params->robusttext, "Manually optimised configuration.\n");
  if (params->robust == 1) 
    sprintf(params->robusttext, "OFAM current initialisation + Smagorinsky = 0.1");
  if (params->robust == 2) 
    sprintf(params->robusttext, "OFAM current initialisation + Smagorinsky = 0.2");
  if (params->robust == 3) 
    sprintf(params->robusttext, "Geostrophic current initialisation + Smagorinsky = 0.1");
  if (params->robust == 4) 
    sprintf(params->robusttext, "Start from rest with hard T/S ramp relaxation + Smagorinsky = 0.1");
  if (params->robust == 5) 
    sprintf(params->robusttext, "Start from rest + Smagorinsky = 0.1");
  if (params->robust == 6) 
    sprintf(params->robusttext, "Start from rest + Smagorinsky = 0.2");
  if (params->robust == 7) 
    sprintf(params->robusttext, "OFAM current initialisation + constant horizontal viscosity");
  if (params->robust == 8) 
    sprintf(params->robusttext, "GStart from rest with hard T/S ramp relaxation + constant horizontal viscosity");
  if (params->robust == 9) 
    sprintf(params->robusttext, "Start from rest + constant horizontal viscosity");
  if (params->robust == 10) 
    sprintf(params->robusttext, "Start from rest + constant horizontal viscosity + 50%% reduced time-step");

  if (params->robust < 0)
    params->robust = 1;
  if (params->robust > 10)
    params->robust = 10;
  if (params->robust == 0 ||
      params->robust == 1 ||
      params->robust == 3 ||
      params->robust == 4 ||
      params->robust == 5) {
      params->smagorinsky = 0.1;
      params->smag_smooth = 2;
      params->atr++;
  } else if (params->robust == 2 ||
	     params->robust == 6) {
    params->smagorinsky = 0.2;
    params->smag_smooth = 4;
    params->atr++;
  } else
    params->smagorinsky = 0.0;

  params->sue1 = params->kue1 = 1.0;
  params->bsue1 = params->bkue1 = 0.0;
  sprintf(params->smag, "%f", params->smagorinsky);

  /* Velocity */
  if (strlen(params->vdata)) {
    if (params->robust > 0 && params->robust <= 2 || params->robust == 7)
      strcpy(params->vel_init, params->vdata);
    else if (params->robust == 3) {
      strcpy(params->vel_init, "GEOSTROPHIC");
      params->rampf |= CUSTOM;
    } else 
      sprintf(params->vel_init, "%c", '\0');
  }

  /* Thresholds */
  params->velmax = 3.87;
  params->velmax2d = 2.55;
  if (!params->etadiff)
    params->etadiff = 0.60;
  if (params->robust == 0) {
    params->speed = 10;
    strcpy(params->alert, "NONE");
    params->trasc = QUICKEST;
    params->ultimate = 1;
    strcpy(params->mixsc, "k-w");
    params->atr += 2;
  } else {
    if (contains_token(params->alert, "PASSIVE") == NULL &&
	contains_token(params->alert, "ACTIVE") == NULL) {
      strcpy(params->alert, "ACTIVE");
      params->ntrS += 4;
    }

    /* Means (optional) */
    if (!(params->means & ETA_M)) {
      if (params->means == NONE)
	params->means = ETA_M;
      else
	params->means |= ETA_M;
      params->ntrS += 1;
      sprintf(params->means_dt, "%f days", get_run_length(fp) + 1);
    }
    params->tide_r = MEAN_R;
  }

  /* Bathymetry checks (optional) */
  sprintf(keyword, "MAXGRAD");
  if (!(prm_read_double(fp, keyword, &params->maxgrad)))
    params->maxgrad = 0.05;
  sprintf(keyword, "SMOOTHING");
  if (!(prm_read_int(fp, keyword, &params->smooth)))
    params->smooth = 1;

  /* Horizontal mixing */
  params->u1vh = params->u2vh = params->u1kh = params->u2kh = -1.0;
  if (params->robust >= 7) params->u1vh = params->u2vh = 1.0;

  /* Surface relaxation (optional) */
  if (params->robust) {
    if (strlen(params->edata))
      strcpy(params->etarlxn, params->edata);
    if (prm_read_char(fp, "eta_relaxation_input_dt", buf) > 0 &&
	prm_read_char(fp, "eta_relaxation_time_constant", keyword) > 0) {
      tm_scale_to_secs(buf, &params->etarlxdt);
      tm_scale_to_secs(keyword, &params->etarlxtc);
      strcpy(params->etarlxtcs, keyword);
      params->etarlx = RELAX;
    } else {
      params->etarlxdt = 3600.0;
      params->etarlxtc = 0.0;
      params->etarlx = ALERT;
    }
    params->ntrS++;
    if (params->robust == 4) {
      strcpy(params->etarlxtcs, "temporal -4 1 day 0 20 days");
      params->etarlx |= RELAX;
    }
  } else
    params->compatible &= ~V1670;

  /* netCDF dumpfile name */
  if (params->runmode & RE_ROAM) {
    sprintf(keyword, "RESTART_FILE");
    prm_read_char(fp, keyword, params->idumpname);
  }
  sprintf(codeheader, "SHOC ROAMv2 : ROBUST %d", params->robust);
}

/* END auto_params_roam_pre2()                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to overwrite input parameters to ROAM values.             */
/* Developed Nov 2013 to include updated ROBUST parameterisation.    */
/* ROBUST=1: Initialisation with OFAM currents                       */
/* ROBUST=2: Initialisation with geostrophic currents                */
/* ROBUST=3: Rest + hard T/S ramp relaxation                         */
/* ROBUST=4: Start from rest                                         */
/* ROBUST=5: As for ROBUST=4 with reduced time-step                  */
/*-------------------------------------------------------------------*/
void auto_params_roam_pre3(FILE * fp, parameters_t *params)
{
  int n;                        /* Counters */
  open_bdrys_t *open;           /* Pointer to boundary structure */
  char keyword[MAXSTRLEN];      /* Input dummy */
  char buf[MAXSTRLEN];          /* Input dummy */

  /* Parameter header (optional) */
  sprintf(keyword, "PARAMETERHEADER");
  if (!(prm_read_char(fp, keyword, params->parameterheader))) {
    strcpy(params->parameterheader, "ROAMv3 grid");
  }
  strcpy(parameterheader, params->parameterheader);

  /*-----------------------------------------------------------------*/
  /* ROAM forcing data                                               */
  /* OFAM_DATA : used for obc forcing and IC's                       */
  /* INIT_DATA : if present used for IC instead of OFAM_DATA         */
  /* TEMP_DATA : if present used for T IC instead of INIT_DATA       */
  /* SALT_DATA : if present used for S IC instead of INIT_DATA       */
  /* ETA_DATA : if present used for eta IC instead of INIT_DATA      */
  char *fields[MAXSTRLEN * MAXNUMARGS];
  int nf;
  sprintf(params->rdata, "%c", '\0');
  sprintf(params->odata, "%c", '\0');
  sprintf(params->idata, "%c", '\0');
  sprintf(params->tdata, "%c", '\0');
  sprintf(params->sdata, "%c", '\0');
  sprintf(params->edata, "%c", '\0');
  sprintf(params->vdata, "%c", '\0');
  prm_read_char(fp, "RAMS_DATA", params->rdata);
  prm_read_char(fp, "OFAM_DATA", params->odata);
  strcpy(buf, params->odata);
  nf = parseline(buf, fields, MAXNUMARGS);
  if (!(prm_read_char(fp, "INIT_DATA", params->idata))) {
    if (nf >= 1)
      strcpy(params->idata, fields[0]);
  }
  if(!(prm_read_char(fp, "TEMP_DATA", params->tdata)))
    strcpy(params->tdata, params->idata);
  if(!(prm_read_char(fp, "SALT_DATA", params->sdata)))
    strcpy(params->sdata, params->idata);
  if(!(prm_read_char(fp, "ETA_DATA", params->edata)))
    strcpy(params->edata, params->idata);
  if(!(prm_read_char(fp, "VELOCITY_DATA", params->vdata)))
    strcpy(params->vdata, params->idata);
  if(!strlen(params->tdata) || !strlen(params->sdata) ||
     !strlen(params->edata))
    hd_quit("auto_params: Insufficient initialisation data specified\n");
  if(!strlen(params->odata))
    sprintf(params->odata,"%s %s %s",params->edata, params->tdata,
	    params->sdata);
  params->rampf = WIND|TIDALH|TIDALC|INV_BARO;
  params->save_force = OTEMP|OSALT|OETA;
  params->save_force = NONE;
  params->save_force = OTEMP|OSALT;
  params->save_force = OTEMP|OSALT|ROAM;
  params->save_force = OTEMP|OSALT|OVELU|OVELV|ROAM;
  params->save_force = OVELU|OVELV|ROAM;
  params->atr += 2;    /* ovelu and ovelv */

  params->rsalt = params->rtemp = 1;
  params->atr += 2;    /* rtemp and rsalt */

  params->robust = 4;
  params->speed = 5;
  /*params->compatible |= V1670;*/

  /* Mixing scheme (optional) */
  sprintf(keyword, "MIXING_SCHEME");
  if (!prm_read_char(fp, keyword, params->mixsc)) {
    strcpy(params->mixsc, "mellor_yamada_2_0");
    params->vz0 = 1e-5;
    params->kz0 = 1e-5;
    params->zs = 1.0;
    params->atr -= 2; /* Remove default k-e auto tracers */
  }
  /* Surface height */
  if (strlen(params->edata))
    strcpy(params->eta_init, params->edata);

  /* Wind and pressure */
  if (strlen(params->rdata)) {
    if (strlen(params->wind) == 0) {
      strcpy(params->wind, params->rdata);
      params->wind_dt = 600.0;
    }
    if (strlen(params->patm) == 0) {
      strcpy(params->patm, params->rdata);
      params->patm_dt = 600.0;
    }
  }

  /* Heatflux */
  sprintf(keyword, "HEATFLUX");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "COMP_HEAT_NONE") == 0) {
      params->heatflux = COMP_HEAT_NONE;
      sprintf(params->hf, "%c", '\0');
      prm_set_errfn(hd_quit);
      sprintf(keyword, "HEATFLUX_FILE");
      prm_read_char(fp, keyword, params->hf);
      sprintf(keyword, "HEATFLUX_DT");
      prm_get_time_in_secs(fp, keyword, &params->hf_dt);
      params->ntrS += 4+2; // 4 x hf comps + precip and evap
    } else {
      if (strcmp(buf, "NONE") != 0) {
	params->smooth_VzKz = 1;
	/* Set defaults */
	/*params->albedo = 0.0;*/
	if (params->swr_type & NONE) {
	  strcpy(params->swr_attn,"0.2");
	  strcpy(params->swr_tran,"0.42");
	  params->water_type = TYPE_II;
	}
	params->hf_ramp = -43200.0;
	/* Override as required 
	read_swr(params, fp, 0);*/
	sprintf(keyword, "ALBEDO");
	prm_read_double(fp, keyword, &params->albedo);
      }
      if (strcmp(buf, "COMP_HEAT") == 0 && strlen(params->rdata)) {
	params->heatflux = COMP_HEAT;
	strcpy(params->hf, params->rdata);
	params->hf_dt = 600.0;
	strcpy(params->swr, params->rdata);
	params->swr_dt = 600.0;
	strcpy(params->rh, params->rdata);
	params->rh_dt = 600.0;
	params->ntrS += 5;
      }
      if ((strcmp(buf, "ADVANCED") == 0 || strcmp(buf, "BULK") == 0)
	  && strlen(params->rdata)) {
	params->heatflux = ADVANCED;
	/*
	n = 3;
	if (params->water_type != NONE) n -= 2;
	if (prm_skip_to_end_of_key(fp, "SWR_ATTENUATION")) n -= 1;
	if (prm_skip_to_end_of_key(fp, "SWR_TRANSMISSION")) n -= 1;
	if (prm_skip_to_end_of_key(fp, "SWR_BOT_ABSORB")) n -= 1;
	read_swr(params, fp, 0);
	params->ntrS += n;
	*/
	read_hf_bulk(params, fp);
	strcpy(params->airtemp, params->rdata);
	params->airtemp_dt = 600.0;
	strcpy(params->wetb, params->rdata);
	params->wetb_dt = 600.0;
	strcpy(params->cloud, params->rdata);
	params->cloud_dt = 600.0;
	/*
	  strcpy(params->swr, params->rdata);
	  params->swr_dt = 600.0;
	*/
	params->ntrS += 5;
      }

    }
  }

  /*
   * Add dummy string to force writing of transport file
   */
  if (params->runmode & PRE_MARVL)
    sprintf(params->trans_data, "PRE_MARVL");

  /*
   * Read precipitation and evaporation, if specified
   */
  sprintf(params->precip, "%c", '\0');
  sprintf(keyword, "PRECIPITATION");
  if (prm_read_char(fp, keyword, params->precip)) {
    sprintf(keyword, "PRECIPITATION_INPUT_DT");
    prm_get_time_in_secs(fp, keyword, &params->precip_dt);
  }
  sprintf(params->evap, "%c", '\0');
  sprintf(keyword, "EVAPORATION");
  if (prm_read_char(fp, keyword, params->evap)) {
    sprintf(keyword, "EVAPORATION_INPUT_DT");
    prm_get_time_in_secs(fp, keyword, &params->evap_dt);
  }

  /* ROAM speed (optional) */
  sprintf(keyword, "SPEED");
  prm_read_int(fp, keyword, &params->speed);

  /* ROAM robustness (optional) */
  sprintf(keyword, "ROBUST");
  prm_read_int(fp, keyword, &params->robust);
  if (params->robust == 0)
    sprintf(params->robusttext, "Manually optimised configuration.\n");
  if (params->robust == 1) 
    sprintf(params->robusttext, "OFAM current initialisation");
  if (params->robust == 2) 
    sprintf(params->robusttext, "Geostrophic current initialisation");
  if (params->robust == 3) 
    sprintf(params->robusttext, "Start from rest with hard T/S ramp relaxation");
  if (params->robust == 4) 
    sprintf(params->robusttext, "Start from rest");
  if (params->robust == 5) 
    sprintf(params->robusttext, "Start from rest + 50%% reduced time-step");

  if (params->robust < 0)
    params->robust = 1;
  if (params->robust > 5)
    params->robust = 5;
  if (!strlen(params->smag)) {
    strcpy(params->smag, "100.0 0.1");
    params->sue1 = 0.0;
    params->bsue1 = 100.0;
    params->kue1 = 0.1;
    params->bkue1 = 0.0;
    params->smagorinsky = 1.0;
    params->atr ++;
  }
  params->bsue1 /= 100.0;
  params->bkue1 /= 100.0;

  /* Velocity */
  if (strlen(params->vdata)) {
    if (params->robust > 0 && params->robust == 1)
      strcpy(params->vel_init, params->vdata);
    else if (params->robust == 2) {
      strcpy(params->vel_init, "GEOSTROPHIC");
      params->rampf |= CUSTOM;
    } else 
      sprintf(params->vel_init, "%c", '\0');
  }

  /* Thresholds */
  params->velmax = 3.87;
  params->velmax2d = 2.55;
  if (!params->etadiff)
    params->etadiff = 0.60;
  if (params->robust == 0) {
    params->speed = 10;
    strcpy(params->alert, "NONE");
    params->trasc = QUICKEST;
    params->ultimate = 1;
    strcpy(params->mixsc, "k-w");
    params->atr += 2;
  } else {
    if (contains_token(params->alert, "PASSIVE") == NULL &&
	contains_token(params->alert, "ACTIVE") == NULL) {
      strcpy(params->alert, "ACTIVE");
      params->ntrS += 4;
    }

    /* Means (optional) */
    if (!(params->means & ETA_M)) {
      if (params->means == NONE)
	params->means = ETA_M;
      else
	params->means |= ETA_M;
      params->ntrS += 1;
      sprintf(params->means_dt, "%f days", get_run_length(fp) + 1);
    }
    params->tide_r = MEAN_R;
  }

  /* Bathymetry checks (optional) */
  sprintf(keyword, "MAXGRAD");
  if (!(prm_read_double(fp, keyword, &params->maxgrad)))
    params->maxgrad = 0.05;
  sprintf(keyword, "SMOOTHING");
  if (!(prm_read_int(fp, keyword, &params->smooth)))
    params->smooth = 1;

  /* Horizontal mixing */
  params->u1vh = params->u2vh = params->u1kh = params->u2kh = -1.0;

  /* Surface relaxation (optional) */
  if (params->robust) {
    if (strlen(params->edata))
      strcpy(params->etarlxn, params->edata);
    if (prm_read_char(fp, "eta_relaxation_input_dt", buf) > 0 &&
	prm_read_char(fp, "eta_relaxation_time_constant", keyword) > 0) {
      tm_scale_to_secs(buf, &params->etarlxdt);
      tm_scale_to_secs(keyword, &params->etarlxtc);
      strcpy(params->etarlxtcs, keyword);
      params->etarlx = RELAX;
    } else {
      params->etarlxdt = 3600.0;
      params->etarlxtc = 0.0;
      params->etarlx = ALERT;
    }
    params->ntrS++;
    if (params->robust == 3) {
      strcpy(params->etarlxtcs, "temporal -4 1 day 0 20 days");
      params->etarlx |= RELAX;
    }
  } else
    params->compatible &= ~V1670;

  /* netCDF dumpfile name */
  if (params->runmode & RE_ROAM) {
    sprintf(keyword, "RESTART_FILE");
    prm_read_char(fp, keyword, params->idumpname);
  }
  sprintf(codeheader, "SHOC ROAMv3 : ROBUST %d", params->robust);
}

/* END auto_params_roam_pre3()                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to overwrite input parameters to ROAM values.             */
/* Developed Nov 2013 to include updated ROBUST parameterisation.    */
/* ROBUST=1: Initialisation with OFAM currents, dual relaxation      */
/* ROBUST=2: Initialisation with OFAM currents, single relaxation    */
/* ROBUST=3: Rest + hard T/S ramp relaxation, dual relaxation        */
/* ROBUST=4: Rest + hard T/S ramp relaxation, single relaxation      */
/* ROBUST=5: Start from rest, dual relaxation                        */
/* ROBUST=6: Start from rest, single relaxation                      */
/*-------------------------------------------------------------------*/
void auto_params_roam_pre4(FILE * fp, parameters_t *params)
{
  int n;                        /* Counters */
  open_bdrys_t *open;           /* Pointer to boundary structure */
  char keyword[MAXSTRLEN];      /* Input dummy */
  char buf[MAXSTRLEN];          /* Input dummy */

  /* Parameter header (optional) */
  sprintf(keyword, "PARAMETERHEADER");
  if (!(prm_read_char(fp, keyword, params->parameterheader))) {
    strcpy(params->parameterheader, "ROAMv3 grid");
  }
  strcpy(parameterheader, params->parameterheader);

  /*-----------------------------------------------------------------*/
  /* ROAM forcing data                                               */
  /* OFAM_DATA : used for obc forcing and IC's                       */
  /* INIT_DATA : if present used for IC instead of OFAM_DATA         */
  /* TEMP_DATA : if present used for T IC instead of INIT_DATA       */
  /* SALT_DATA : if present used for S IC instead of INIT_DATA       */
  /* ETA_DATA : if present used for eta IC instead of INIT_DATA      */
  char *fields[MAXSTRLEN * MAXNUMARGS];
  int nf;
  sprintf(params->rdata, "%c", '\0');
  sprintf(params->odata, "%c", '\0');
  sprintf(params->idata, "%c", '\0');
  sprintf(params->tdata, "%c", '\0');
  sprintf(params->sdata, "%c", '\0');
  sprintf(params->edata, "%c", '\0');
  sprintf(params->vdata, "%c", '\0');
  prm_read_char(fp, "RAMS_DATA", params->rdata);
  prm_read_char(fp, "OFAM_DATA", params->odata);
  strcpy(buf, params->odata);
  nf = parseline(buf, fields, MAXNUMARGS);
  if (!(prm_read_char(fp, "INIT_DATA", params->idata))) {
    if (nf >= 1)
      strcpy(params->idata, fields[0]);
  }
  if(!(prm_read_char(fp, "TEMP_DATA", params->tdata)))
    strcpy(params->tdata, params->idata);
  if(!(prm_read_char(fp, "SALT_DATA", params->sdata)))
    strcpy(params->sdata, params->idata);
  if(!(prm_read_char(fp, "ETA_DATA", params->edata)))
    strcpy(params->edata, params->idata);
  if(!(prm_read_char(fp, "VELOCITY_DATA", params->vdata)))
    strcpy(params->vdata, params->idata);
  if(!strlen(params->tdata) || !strlen(params->sdata) ||
     !strlen(params->edata))
    hd_quit("auto_params: Insufficient initialisation data specified\n");
  if(!strlen(params->odata))
    sprintf(params->odata,"%s %s %s",params->edata, params->tdata,
	    params->sdata);
  params->rampf = WIND|TIDALH|TIDALC|INV_BARO;
  params->save_force = OTEMP|OSALT|OETA;
  params->save_force = NONE;
  params->save_force = OTEMP|OSALT;
  params->save_force = OTEMP|OSALT|ROAM;
  params->save_force = OTEMP|OSALT|OVELU|OVELV|ROAM;
  params->save_force = OVELU|OVELV|ROAM;
  params->atr += 2;    /* ovelu and ovelv */

  params->rsalt = params->rtemp = 1;
  params->atr += 2;    /* rtemp and rsalt */

  params->robust = 5;
  params->speed = 5;
  params->compatible |= V1670;

  /* Mixing scheme (optional) */
  sprintf(keyword, "MIXING_SCHEME");
  if (!prm_read_char(fp, keyword, params->mixsc)) {
    strcpy(params->mixsc, "mellor_yamada_2_0");
    params->vz0 = 1e-5;
    params->kz0 = 1e-5;
    params->zs = 1.0;
    params->atr -= 2; /* Remove default k-e auto tracers */
  }
  /* Surface height */
  if (strlen(params->edata))
    strcpy(params->eta_init, params->edata);

  /* Wind and pressure */
  if (strlen(params->rdata)) {
    if (strlen(params->wind) == 0) {
      strcpy(params->wind, params->rdata);
      params->wind_dt = 600.0;
    }
    if (strlen(params->patm) == 0) {
      strcpy(params->patm, params->rdata);
      params->patm_dt = 600.0;
    }
  }

  /* Heatflux */
  sprintf(keyword, "HEATFLUX");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "COMP_HEAT_NONE") == 0) {
      params->heatflux = COMP_HEAT_NONE;
      sprintf(params->hf, "%c", '\0');
      prm_set_errfn(hd_quit);
      sprintf(keyword, "HEATFLUX_FILE");
      prm_read_char(fp, keyword, params->hf);
      sprintf(keyword, "HEATFLUX_DT");
      prm_get_time_in_secs(fp, keyword, &params->hf_dt);
      params->ntrS += 4+2; // 4 x hf comps + precip and evap
    } else {
      params->smooth_VzKz = 1;
      /* Set defaults */
      /*params->albedo = 0.0;*/
      if (params->swr_type & NONE) {
	strcpy(params->swr_attn,"0.2");
	strcpy(params->swr_tran,"0.42");
	params->water_type = TYPE_II;
      }
      params->hf_ramp = -43200.0;
      /* Override as required 
	 read_swr(params, fp, 0); */
      sprintf(keyword, "ALBEDO");
      prm_read_double(fp, keyword, &params->albedo);
      
      if (strcmp(buf, "COMP_HEAT") == 0 && strlen(params->rdata)) {
	params->heatflux = COMP_HEAT;
	strcpy(params->hf, params->rdata);
	params->hf_dt = 600.0;
	strcpy(params->swr, params->rdata);
	params->swr_dt = 600.0;
	strcpy(params->rh, params->rdata);
	params->rh_dt = 600.0;
	params->ntrS += 5;
      }
      if ((strcmp(buf, "ADVANCED") == 0 || strcmp(buf, "BULK") == 0)
	  && strlen(params->rdata)) {
	params->heatflux = ADVANCED;
	/*
	n = 3;
	if (params->water_type != NONE) n -= 2;
	if (prm_skip_to_end_of_key(fp, "SWR_ATTENUATION")) n -= 1;
	if (prm_skip_to_end_of_key(fp, "SWR_TRANSMISSION")) n -= 1;
	if (prm_skip_to_end_of_key(fp, "SWR_BOT_ABSORB")) n -= 1;
	read_swr(params, fp, 0);
	params->ntrS += n;
	*/
	read_hf_bulk(params, fp);
	strcpy(params->airtemp, params->rdata);
	params->airtemp_dt = 600.0;
	strcpy(params->wetb, params->rdata);
	params->wetb_dt = 600.0;
	strcpy(params->cloud, params->rdata);
	params->cloud_dt = 600.0;
	/*
	  strcpy(params->swr, params->rdata);
	  params->swr_dt = 600.0;
	*/
	params->ntrS += 5;
      }
    }
  }

  /*
   * Add dummy string to force writing of transport file
   */
  if (params->runmode & PRE_MARVL)
    sprintf(params->trans_data, "PRE_MARVL");

  /*
   * Read precipitation and evaporation, if specified
   */
  sprintf(params->precip, "%c", '\0');
  sprintf(keyword, "PRECIPITATION");
  if (prm_read_char(fp, keyword, params->precip)) {
    sprintf(keyword, "PRECIPITATION_INPUT_DT");
    prm_get_time_in_secs(fp, keyword, &params->precip_dt);
  }
  sprintf(params->evap, "%c", '\0');
  sprintf(keyword, "EVAPORATION");
  if (prm_read_char(fp, keyword, params->evap)) {
    sprintf(keyword, "EVAPORATION_INPUT_DT");
    prm_get_time_in_secs(fp, keyword, &params->evap_dt);
  }

  /* ROAM speed (optional) */
  sprintf(keyword, "SPEED");
  prm_read_int(fp, keyword, &params->speed);

  /* ROAM robustness (optional) */
  sprintf(keyword, "ROBUST");
  prm_read_int(fp, keyword, &params->robust);
  if (params->robust == 0) 
    sprintf(params->robusttext, "Manually optimised configuration.\n");
  if (params->robust == 1) 
    sprintf(params->robusttext, "OFAM current initialisation, dual flux adjustment");
  if (params->robust == 2) 
    sprintf(params->robusttext, "OFAM current initialisation\n");
  if (params->robust == 3) 
    sprintf(params->robusttext, "Start from rest with hard T/S ramp relaxation, dual flux adjustment");
  if (params->robust == 4) 
    sprintf(params->robusttext, "Start from rest with hard T/S ramp relaxation\n");
  if (params->robust == 5) 
    sprintf(params->robusttext, "Start from rest, dual flux adjustment");
  if (params->robust == 6) 
    sprintf(params->robusttext, "Start from rest\n");

  if (params->robust < 0)
    params->robust = 1;
  /* MH test
  if (params->robust > 5)
    params->robust = 5;
  */
  if (!strlen(params->smag)) {
    strcpy(params->smag, "100.0 100.0 0.1 0.1");
    params->sue1 = params->sue2 = 0.0;
    params->bsue1 = params->bsue2 =100.0;
    params->kue1 = params->kue2 = 0.1;
    params->bkue1 = params->bkue2 = 0.0;
    params->smagorinsky = 1.0;
    params->atr ++;
  }
  params->bsue1 /= 100.0;
  params->bsue2 /= 100.0;
  params->bkue1 /= 100.0;
  params->bkue2 /= 100.0;

  /* Velocity */
  if (strlen(params->vdata)) {
    if (params->robust == 1 || params->robust == 2)
      strcpy(params->vel_init, params->vdata);
    else 
      sprintf(params->vel_init, "%c", '\0');
  }

  /* Thresholds */
  params->velmax = 3.87;
  params->velmax2d = 2.55;
  if (!params->etadiff)
    params->etadiff = 0.60;
  if (params->robust == 0) {
    params->speed = 10;
    strcpy(params->alert, "NONE");
    params->trasc = QUICKEST;
    params->ultimate = 1;
    strcpy(params->mixsc, "k-w");
    params->atr += 2;
  } else {
    if (contains_token(params->alert, "PASSIVE") == NULL &&
	contains_token(params->alert, "ACTIVE") == NULL) {
      strcpy(params->alert, "ACTIVE");
      params->ntrS += 4;
    }

    /* Means (optional) */
    if (!(params->means & ETA_M)) {
      if (params->means == NONE)
	params->means = ETA_M;
      else
	params->means |= ETA_M;
      params->ntrS += 1;
      sprintf(params->means_dt, "%f days", get_run_length(fp) + 1);
    }
    params->tide_r = MEAN_R;
  }

  /* Bathymetry checks (optional) */
  sprintf(keyword, "MAXGRAD");
  if (!(prm_read_double(fp, keyword, &params->maxgrad)))
    params->maxgrad = 0.05;
  sprintf(keyword, "SMOOTHING");
  if (!(prm_read_int(fp, keyword, &params->smooth)))
    params->smooth = 1;

  /* Horizontal mixing */
  params->u1vh = params->u2vh = params->u1kh = params->u2kh = -1.0;

  /* Surface relaxation (optional) */
  if (params->robust == 3 || params->robust == 4) {
    strcpy(params->etarlxn, params->edata);
    params->etarlxdt = 3600.0;
    params->etarlxtc = 0.0;
    params->etarlx = ALERT;
    params->ntrS++;
    strcpy(params->etarlxtcs, "temporal -4 1 day 0 20 days");
    params->etarlx |= RELAX;
  } else
    params->compatible &= ~V1670;

  /* netCDF dumpfile name */
  if (params->runmode & RE_ROAM) {
    sprintf(keyword, "RESTART_FILE");
    prm_read_char(fp, keyword, params->idumpname);
  }
  sprintf(codeheader, "SHOC ROAMv4 : ROBUST %d", params->robust);
}

/* END auto_params_roam_pre4()                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to overwrite input parameters to ROAM values              */
/*-------------------------------------------------------------------*/
void auto_params_recom_pre1(FILE * fp, parameters_t *params)
{
  int n;                        /* Counters */
  open_bdrys_t *open;           /* Pointer to boundary structure */
  char keyword[MAXSTRLEN];      /* Input dummy */
  char buf[MAXSTRLEN];          /* Input dummy */

  params->runmode |= DUMP;
  sprintf(params->sourcefile, "%c", '\0');

  /* Parameter header (optional) */
  sprintf(keyword, "PARAMETERHEADER");
  if (!(prm_read_char(fp, keyword, params->parameterheader))) {
    strcpy(params->parameterheader, "RECOM grid");
  }
  strcpy(parameterheader, params->parameterheader);

  /*-----------------------------------------------------------------*/
  /* ROAM forcing data                                               */
  /* OFAM_DATA : used for obc forcing and IC's                       */
  /* INIT_DATA : if present used for IC instead of OFAM_DATA         */
  /* TEMP_DATA : if present used for T IC instead of INIT_DATA       */
  /* SALT_DATA : if present used for S IC instead of INIT_DATA       */
  /* ETA_DATA : if present used for eta IC instead of INIT_DATA      */
  char *fields[MAXSTRLEN * MAXNUMARGS];
  int nf;
  sprintf(params->rdata, "%c", '\0');
  sprintf(params->odata, "%c", '\0');
  sprintf(params->idata, "%c", '\0');
  sprintf(params->tdata, "%c", '\0');
  sprintf(params->sdata, "%c", '\0');
  sprintf(params->edata, "%c", '\0');
  sprintf(params->vdata, "%c", '\0');
  prm_read_char(fp, "MET_DATA", params->rdata);
  prm_read_char(fp, "OCEAN_DATA", params->odata);
  if (prm_read_char(fp, "SEDIMENT_DATA", buf))
    sprintf(params->odata, "%s %s", params->odata, buf);
  if (prm_read_char(fp, "ECOLOGY_DATA", buf))
    sprintf(params->odata, "%s %s", params->odata, buf);
  strcpy(buf, params->odata);
  nf = parseline(buf, fields, MAXNUMARGS);
  if (!(prm_read_char(fp, "INIT_DATA", params->idata))) {
    if (nf >= 1)
      strcpy(params->idata, fields[0]);
  }
  prm_read_char(fp, "TRACER_DATA", params->tracerdata);
  if(!(prm_read_char(fp, "TEMP_DATA", params->tdata)))
    strcpy(params->tdata, params->idata);
  if(!(prm_read_char(fp, "SALT_DATA", params->sdata)))
    strcpy(params->sdata, params->idata);
  if(!(prm_read_char(fp, "ETA_DATA", params->edata)))
    strcpy(params->edata, params->idata);
  if(!(prm_read_char(fp, "VELOCITY_DATA", params->vdata)))
    strcpy(params->vdata, params->idata);
  prm_read_char(fp, "BOUNDARY_DATA", params->vdata);
  if(!strlen(params->tdata) || !strlen(params->sdata) ||
     !strlen(params->edata))
    hd_quit("auto_params: Insufficient initialisation data specified\n");
  if(!strlen(params->odata))
    sprintf(params->odata,"%s %s %s",params->edata, params->tdata,
	    params->sdata);
  params->rampf = WIND|TIDALH|INV_BARO;
  params->save_force = NONE;
  params->rsalt = params->rtemp = 0;
  params->robust = 1;
  params->speed = 10;

  /* Mixing scheme (optional) */
  sprintf(keyword, "MIXING_SCHEME");
  if (!prm_read_char(fp, keyword, params->mixsc)) {
    strcpy(params->mixsc, "k-e");
    params->vz0 = 1e-6;
    params->kz0 = 1e-6;
    params->zs = 0.2;
  }
  /* Surface height */
  if (strlen(params->edata))
    strcpy(params->eta_init, params->edata);

  /* Wind and pressure */
  if (strlen(params->rdata)) {
    if (strlen(params->wind) == 0) {
      strcpy(params->wind, params->rdata);
      params->wind_dt = 600.0;
    }
    if (strlen(params->patm) == 0) {
      strcpy(params->patm, params->rdata);
      params->patm_dt = 600.0;
    }
  }

  /* Heatflux */
  params->heatflux = ADVANCED;
  params->smooth_VzKz = 1;
  /* Set defaults */
  strcpy(params->swr_attn,"0.2");
  strcpy(params->swr_tran,"0.42");
  params->water_type = TYPE_II;
  n = 3;
  if (params->water_type != NONE) n -= 2;
  if (prm_skip_to_end_of_key(fp, "SWR_ATTENUATION")) n -= 1;
  if (prm_skip_to_end_of_key(fp, "SWR_TRANSMISSION")) n -= 1;
  if (prm_skip_to_end_of_key(fp, "SWR_BOT_ABSORB")) n -= 1;
  /*params->albedo = 0.0;*/
  read_swr(params, fp, 0);
  params->ntrS += n;
  read_hf_bulk(params, fp);
  strcpy(params->airtemp, params->rdata);
  params->airtemp_dt = 600.0;
  strcpy(params->wetb, params->rdata);
  params->wetb_dt = 600.0;
  strcpy(params->cloud, params->rdata);
  params->cloud_dt = 600.0;
  params->ntrS += 5;

  /* ROAM speed (optional) */
  if ((params->robust >= 1 && params->robust <= 7) && 
      params->smagorinsky == 0.0)
    params->atr++;
  if (params->robust < 1)
    params->robust = 1;
  if (params->robust > 10)
    params->robust = 10;
  if (params->robust == 1)
      params->smagorinsky = 0.1;
  else if (params->robust == 2)
    params->smagorinsky = 0.2;
  else if (params->robust == 3)
    params->smagorinsky = 0.3;
  else if (params->robust == 4)
    params->smagorinsky = 0.4;
  else if (params->robust == 5)
    params->smagorinsky = 0.5;
  else if (params->robust == 6 || params->robust == 7)
    params->smagorinsky = 0.1;
  else
    params->smagorinsky = 0.0;
  params->sue1 = params->kue1 = 1.0;
  params->bsue1 = params->bkue1 = 0.0;

  /* Bathymetry checks (optional) */
  params->maxgrad = 0.05;
  params->smooth = 1.0;

  /* Horizontal mixing */
  params->u1vh = params->u2vh = params->u1kh = params->u2kh = -1.0;
  if (params->robust == 6 || params->robust == 7) params->u1vh = params->u2vh = 1.0;

  /* Transport */
  params->trout = 1;
  strcpy(params->trans_data, "trans.mnc");
  prm_read_char(fp, "TRANS_DATA", params->trans_data);
  /* Force flux-form */
  params->tmode = SP_FFSL;
  params->atr += 2; /* u1vmean and u2vmean */

  /*
   * Waves, Sediments and Ecology already read in auto_params
   */

}

/* END auto_params_recom_pre1()                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to overwrite input parameters to ROAM values              */
/* Developed May 2015 to include updated ROBUST parameterisation.    */
/* ROBUST=1: Standard parameterisation, Smagorinsky, no alerts       */
/* ROBUST=2: Standard parameterisation, const. viscosity, no alerts  */
/* ROBUST=3: Standard parameterisation, Smagorinsky, active alerts   */
/* ROBUST=4: Standard parameterisation, const. viscosity, alerts     */
/* ROBUST=5: ROAM parameteristion, Smagorinsky = 0.1, rest start     */
/* ROBUST=6: ROAM parameteristion, Smagorinsky = 0.1, OFAM start     */
/* ROBUST=7: Rest + hard T/S ramp relaxation + Smagorinsky = 0.1     */
/* ROBUST=8: Rest + hard T/S ramp relaxation + constant viscosity    */
/* ROBUST=9: OFAM currents + constant horizontal viscosity           */
/* ROBUST=10: Start from rest + constant horizontal viscosity        */
/*-------------------------------------------------------------------*/
void auto_params_recom_pre2(FILE * fp, parameters_t *params)
{
  int n;                        /* Counters */
  open_bdrys_t *open;           /* Pointer to boundary structure */
  char keyword[MAXSTRLEN];      /* Input dummy */
  char buf[MAXSTRLEN];          /* Input dummy */

  params->runmode |= DUMP;
  sprintf(params->sourcefile, "%c", '\0');

  /* Parameter header (optional) */
  sprintf(keyword, "PARAMETERHEADER");
  if (!(prm_read_char(fp, keyword, params->parameterheader))) {
    strcpy(params->parameterheader, "RECOM grid");
  }
  strcpy(parameterheader, params->parameterheader);

  /*-----------------------------------------------------------------*/
  /* ROAM forcing data                                               */
  /* OFAM_DATA : used for obc forcing and IC's                       */
  /* INIT_DATA : if present used for IC instead of OFAM_DATA         */
  /* TEMP_DATA : if present used for T IC instead of INIT_DATA       */
  /* SALT_DATA : if present used for S IC instead of INIT_DATA       */
  /* ETA_DATA : if present used for eta IC instead of INIT_DATA      */
  char *fields[MAXSTRLEN * MAXNUMARGS];
  int nf;
  sprintf(params->rdata, "%c", '\0');
  sprintf(params->odata, "%c", '\0');
  sprintf(params->idata, "%c", '\0');
  sprintf(params->tdata, "%c", '\0');
  sprintf(params->sdata, "%c", '\0');
  sprintf(params->edata, "%c", '\0');
  sprintf(params->vdata, "%c", '\0');
  prm_read_char(fp, "MET_DATA", params->rdata);
  prm_read_char(fp, "OCEAN_DATA", params->odata);
  if (prm_read_char(fp, "SEDIMENT_DATA", buf))
    sprintf(params->odata, "%s %s", params->odata, buf);
  if (prm_read_char(fp, "ECOLOGY_DATA", buf))
    sprintf(params->odata, "%s %s", params->odata, buf);
  strcpy(buf, params->odata);
  nf = parseline(buf, fields, MAXNUMARGS);
  if (!(prm_read_char(fp, "INIT_DATA", params->idata))) {
    if (nf >= 1)
      strcpy(params->idata, fields[0]);
  }

  read_trfilter(params, fp);

  prm_read_char(fp, "TRACER_DATA", params->tracerdata);
  if(!(prm_read_char(fp, "TEMP_DATA", params->tdata)))
    strcpy(params->tdata, params->idata);
  if(!(prm_read_char(fp, "SALT_DATA", params->sdata)))
    strcpy(params->sdata, params->idata);
  if(!(prm_read_char(fp, "ETA_DATA", params->edata)))
    strcpy(params->edata, params->idata);
  if(!(prm_read_char(fp, "VELOCITY_DATA", params->vdata)))
    strcpy(params->vdata, params->idata);
  prm_read_char(fp, "BOUNDARY_DATA", params->vdata);
  if(!strlen(params->tdata) || !strlen(params->sdata) ||
     !strlen(params->edata))
    hd_quit("auto_params: Insufficient initialisation data specified\n");
  if(!strlen(params->odata))
    sprintf(params->odata,"%s %s %s",params->edata, params->tdata,
	    params->sdata);
  params->rampf = WIND|TIDALH|INV_BARO;
  params->save_force = NONE;
  params->rsalt = params->rtemp = 0;
  params->robust = 1;
  params->speed = 10;

  /* Surface height */
  if (strlen(params->edata))
    strcpy(params->eta_init, params->edata);

  /* Wind and pressure */
  if (strlen(params->rdata)) {
    if (strlen(params->wind) == 0) {
      strcpy(params->wind, params->rdata);
      params->wind_dt = 600.0;
    }
    if (strlen(params->patm) == 0) {
      strcpy(params->patm, params->rdata);
      params->patm_dt = 600.0;
    }
  }

  /* Read transport timestep */
  prm_read_char(fp, "TRANS_DT", params->trans_dt);

  /* Heatflux */
  params->heatflux = ADVANCED;
  params->smooth_VzKz = 1;
  /* Set defaults */
  strcpy(params->swr_attn,"0.2");
  strcpy(params->swr_tran,"0.42");
  params->water_type = TYPE_II;
  n = 3;
  if (params->water_type != NONE) n -= 2;
  if (prm_skip_to_end_of_key(fp, "SWR_ATTENUATION")) n -= 1;
  if (prm_skip_to_end_of_key(fp, "SWR_TRANSMISSION")) n -= 1;
  if (prm_skip_to_end_of_key(fp, "SWR_BOT_ABSORB")) n -= 1;
  /*params->albedo = 0.0;*/
  read_swr(params, fp, 0);
  params->ntrS += n;
  read_hf_bulk(params, fp);
  strcpy(params->airtemp, params->rdata);
  params->airtemp_dt = 600.0;
  strcpy(params->wetb, params->rdata);
  params->wetb_dt = 600.0;
  strcpy(params->cloud, params->rdata);
  params->cloud_dt = 600.0;
  params->ntrS += 5;

  /* ROAM speed (optional) */
  sprintf(keyword, "SPEED");
  prm_read_int(fp, keyword, &params->speed);

  /* ROAM robustness (optional) */
  sprintf(keyword, "ROBUST");
  prm_read_int(fp, keyword, &params->robust);
  if (params->robust == 1) 
    sprintf(params->robusttext, "Standard configuration, no constraint, Smagorinsky, start from rest. No MAX_GRAD\n");
  if (params->robust == 2) 
    sprintf(params->robusttext, "Standard configuration, no constraint, constant viscosity, start from rest.\n");
  if (params->robust == 3) 
    sprintf(params->robusttext, "Standard configuration, constraint invoked, Smagorinsky, start from rest.\n");
  if (params->robust == 4) 
    sprintf(params->robusttext, "Standard configuration, constraint invoked, constant viscosity, start from rest.\n");
  if (params->robust == 5) 
    sprintf(params->robusttext, "ROAM configuration, Smagorinsky, start from rest.\n");
  if (params->robust == 6) 
    sprintf(params->robusttext, "ROAM configuration, Smagorinsky, start with non-zeo currents.\n");
  if (params->robust == 7) 
    sprintf(params->robusttext, "ROAM configuration, Smagorinsky, start from rest with hard T/S relaxation.\n");
  if (params->robust == 8) 
    sprintf(params->robusttext, "ROAM configuration, constant viscosity, start from rest with hard T/S relaxation.\n");
  if (params->robust == 9) 
    sprintf(params->robusttext, "ROAM configuration, constant viscosity, start with non-zeo currents.\n");
  if (params->robust == 10) 
    sprintf(params->robusttext, "ROAM configuration, constant viscosity, start from rest.\n");

  if (params->robust < 0)
    params->robust = 1;
  if (params->robust > 10)
    params->robust = 10;
  if (params->robust >= 3)
    params->compatible |= V1670;
  if (params->robust == 1 ||
      params->robust == 3 ||
      params->robust == 4 ||
      params->robust == 5 ||
      params->robust == 6 ||
      params->robust == 7 ||
      params->robust == 8) {
      params->smagorinsky = 0.1;
      params->atr++;
  } else
    params->smagorinsky = 0.0;
  params->sue1 = params->kue1 = 1.0;
  params->bsue1 = params->bkue1 = 0.0;

  /* Velocity */
  if (strlen(params->vdata)) {
    if (params->robust == 6 || params->robust == 9)
      strcpy(params->vel_init, params->vdata);
    else 
      sprintf(params->vel_init, "%c", '\0');
  }

  /* Advection scheme */
  if (params->robust <= 2) {
    params->trasc = QUICKEST;
    params->ultimate = 1;
  }

  /* Mixing scheme (optional) */
  sprintf(keyword, "MIXING_SCHEME");
  if (!prm_read_char(fp, keyword, params->mixsc)) {
    if (params->robust <= 4) {
      strcpy(params->mixsc, "k-e");
      params->vz0 = 1e-6;
      params->kz0 = 1e-6;
      params->zs = 0.2;
    } else {
      strcpy(params->mixsc, "mellor_yamada_2_0");
      params->vz0 = 1e-5;
      params->kz0 = 1e-5;
      params->zs = 1.0;
      params->atr -= 2; /* Remove default k-e auto tracers */
    }
  }

  /* Thresholds */
  if (params->robust >= 3) {
    params->velmax = 3.87;
    params->velmax2d = 2.55;
    if (!params->etadiff)
      params->etadiff = 0.60;
    if (contains_token(params->alert, "PASSIVE") == NULL &&
	contains_token(params->alert, "ACTIVE") == NULL) {
      strcpy(params->alert, "ACTIVE");
      params->ntrS += 4;
    }

    /* Means (optional) */
    /*
    if (!(params->means & ETA_M)) {
      if (params->means == NONE)
        params->means = ETA_M;
      else
        params->means |= ETA_M;
      params->ntrS += 1;
      sprintf(params->means_dt, "%f days", get_run_length(fp) + 1);
    }
    */
    params->tide_r = NONE;
    params->shear_f = 0;
    params->eta_f = 0;
  }

  /* Bathymetry checks (optional) */
  params->smooth  = 1.0;
  params->maxgrad = 0.05;
  if (params->robust == 1)
    params->maxgrad = 0.0;

  /* Horizontal mixing  (note: reset in autoset_roam()) */
  params->u1vh = params->u2vh = params->u1kh = params->u2kh = -1.0;
  if (params->robust == 1 || params->robust == 4 || params->robust >= 8) params->u1vh = params->u2vh = 1.0;

  /* Surface relaxation (optional) */
  if (params->robust >= 3) {
    if (strlen(params->edata))
      strcpy(params->etarlxn, params->edata);
    if (prm_read_char(fp, "eta_relaxation_input_dt", buf) > 0 &&
        prm_read_char(fp, "eta_relaxation_time_constant", keyword) > 0) {
      tm_scale_to_secs(buf, &params->etarlxdt);
      tm_scale_to_secs(keyword, &params->etarlxtc);
      strcpy(params->etarlxtcs, keyword);
      params->etarlx = RELAX;
    } else {
      params->etarlxdt = 3600.0;
      params->etarlxtc = 0.0;
      params->etarlx = ALERT;
    }
    params->ntrS++;
    /*
    if (params->robust == 7 || params->robust == 8) {
      sprintf(params->etarlxtcs, "temporal %f 1 day %f 20 days",params->t/86400, params->t/86400 + 4);
      params->etarlx |= RELAX;
    }
    */
  }


  /* Transport */
  params->trout = 1;
  strcpy(params->trans_data, "trans.mnc");
  prm_read_char(fp, "TRANS_DATA", params->trans_data);
  /* Force flux-form */
  params->tmode = SP_FFSL;
  params->atr += 2; /* u1vmean and u2vmean */

  /*
   * Waves, Sediments and Ecology already read in auto_params
   */
}

/* END auto_params_recom_pre2()                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to overwrite input parameters to ROAM values              */
/*-------------------------------------------------------------------*/
void auto_params_roam_post1(FILE * fp, parameters_t *params)
{
  int n, tn;                    /* Counters */
  open_bdrys_t *open;           /* Pointer to boundary structure */
  char keyword[MAXSTRLEN];      /* Input dummy */
  char buf[MAXSTRLEN];          /* Input dummy */

  /* Tracer relaxation */
  for (n = 0; n < params->ntr; n++) {
    int tm = params->trinfo_3d[n].m;
    tracer_info_t *tracer = &params->trinfo_3d[n];
    if (strcmp(tracer->name, "temp") == 0 && strlen(params->tdata))
      strcpy(buf, params->tdata);
    else if (strcmp(tracer->name, "salt") == 0 && strlen(params->sdata))
      strcpy(buf, params->sdata);
    else
      continue;
    strcpy(tracer->data, buf);
    if (strlen(params->trrlxn[n]) == 0) {
      strcpy(params->trrlxn[n], buf);
      strcpy(tracer->relax_file, buf);
      params->trrlxdt[n] = 3600.0;
      strcpy(tracer->relax_dt, "1 hour");
      strcpy(params->trrlxtc[n], "20 days");
      strcpy(tracer->r_rate, "20 days");
    }
  }

  /* Tracer reset */
  for (n = 0; n < params->ntr; n++) {
    int tm = params->trinfo_3d[n].m;
    if (strlen(params->trinfo_3d[n].reset_file))
      strcpy(params->trrest[n], params->trinfo_3d[n].reset_file);
    if (strlen(params->trinfo_3d[n].reset_dt))
      tm_scale_to_secs(params->trinfo_3d[n].reset_dt, &params->trrestdt[n]);
  }

  /* Bathymetry limits */
  if (params->bmin < 1.0)
    params->bmin = 1.0;
  if (params->bmax > 20.0 && params->bmin < 2.0)
    params->bmin = 2.0;
  if (params->bmax > 200.0 && params->bmin < 4.0)
    params->bmin = 4.0;

  /* Vertical grid */
  for (n = params->nz - 1; n >= 0; n--) {
    if (params->layers[n] > -LAYERMIN) {
      params->layers[n] = params->layers[params->nz];
      params->nz -= 1;
    }
    if (params->layers[n] >= -LAYERMIN)
      break;
  }
  params->hmin = min(0.1, 0.07 * (params->layers[params->nz] -
				  params->layers[params->nz - 1]));

  /*-----------------------------------------------------------------*/
  /* Open boundaries */
  for (n = 0; n < params->nobc; n++) {
    open = params->open[n];
    if (strlen(params->odata))
      strcpy(open->tsfn, params->odata);
    if(prm_read_char(fp, "TIDE_CSR_ORTHOWEIGHTS", params->orthoweights)) {
      if (prm_read_char(fp, "TIDE_CSR_CON_DIR", params->nodal_dir)) 
	open->bcond_ele = FILEIN|TIDALH;

    } else
      open->bcond_ele = FILEIN;
    open->sponge_zone_h = 8;
    open->stagger = OUTFACE;

    for (tn = 0; tn < params->ntr; tn++) {
      tracer_info_t *tracer = &params->trinfo_3d[tn];
      if ((strcmp(tracer->name, "salt") == 0) ||
	  (strcmp(tracer->name, "temp") == 0))
        open->bcond_tra[tn] = UPSTRM | FILEIN;
      else
        open->bcond_tra[tn] = NOGRAD;
    }
  }
}

/* END auto_params_roam_post1()                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to overwrite input parameters to ROAM values              */
/*-------------------------------------------------------------------*/
void auto_params_roam_post2(FILE * fp, parameters_t *params)
{
  int n, tn;                    /* Counters */
  open_bdrys_t *open;           /* Pointer to boundary structure */
  char keyword[MAXSTRLEN];      /* Input dummy */
  char buf[MAXSTRLEN];          /* Input dummy */


  /* Changes to _pre fromulation */
  params->rampf |= ETA_RELAX;

  /* Tracer relaxation */
  for (n = 0; n < params->ntr; n++) {
    int tm = params->trinfo_3d[n].m;
    tracer_info_t *tracer = &params->trinfo_3d[n];
    if (strcmp(tracer->name, "temp") == 0 && strlen(params->tdata))
      strcpy(buf, params->tdata);
    else if (strcmp(tracer->name, "salt") == 0 && strlen(params->sdata))
      strcpy(buf, params->sdata);
    else
      continue;
    strcpy(tracer->data, buf);
    if (strlen(params->trrlxn[n]) == 0) {
      strcpy(params->trrlxn[n], buf);
      strcpy(tracer->relax_file, buf);
      params->trrlxdt[n] = 3600.0;
      strcpy(tracer->relax_dt, "1 hour");
      strcpy(params->trrlxtc[n], "20 days");
      strcpy(tracer->r_rate, "20 days");
    }
  }

  /* Tracer reset */
  for (n = 0; n < params->ntr; n++) {
    int tm = params->trinfo_3d[n].m;
    if (strlen(params->trinfo_3d[n].reset_file))
      strcpy(params->trrest[n], params->trinfo_3d[n].reset_file);
    if (strlen(params->trinfo_3d[n].reset_dt))
      tm_scale_to_secs(params->trinfo_3d[n].reset_dt, &params->trrestdt[n]);
  }

  /* Bathymetry limits */
  if (params->bmin < 1.0)
    params->bmin = 1.0;
  if (params->bmax > 20.0 && params->bmin < 2.0)
    params->bmin = 2.0;
  if (params->bmax > 200.0 && params->bmin < 4.0)
    params->bmin = 4.0;

  /* Vertical grid */
  for (n = params->nz - 1; n >= 0; n--) {
    if (params->layers[n] > -LAYERMIN) {
      params->layers[n] = params->layers[params->nz];
      params->nz -= 1;
    }
    if (params->layers[n] >= -LAYERMIN)
      break;
  }
  params->hmin = min(0.1, 0.07 * (params->layers[params->nz] -
				  params->layers[params->nz - 1]));

  /*-----------------------------------------------------------------*/
  /* Open boundaries */
  for (n = 0; n < params->nobc; n++) {
    open = params->open[n];
    open->bcond_tan = RAYMND;
    open->bcond_tan2d = CLAMPD;
    if (strlen(params->odata))
      strcpy(open->tsfn, params->odata);
    if(prm_read_char(fp, "TIDE_CSR_ORTHOWEIGHTS", params->orthoweights)) {
      if (prm_read_char(fp, "TIDE_CSR_CON_DIR", params->nodal_dir))
	open->bcond_ele = FILEIN|TIDALH;
      /*open->bcond_ele = FILEIN|TIDALH|RAYMND;*/

    } else
      open->bcond_ele = FILEIN;
    /*
      open->bcond_ele = FILEIN|RAYMND;
    open->relax_timei = 72.0 * 60.0;
    open->relax_time = 2.0 * 86400.0;
    open->relax_ele = 8; 
    open->rele_b = 10;
    open->rele_i = 100;
    */
    open->sponge_zone_h = 8;
    if(strlen(params->patm)) open->inverse_barometer = 1;
    open->stagger = OUTFACE;

    for (tn = 0; tn < params->ntr; tn++) {
      tracer_info_t *tracer = &params->trinfo_3d[tn];
      if ((strcmp(tracer->name, "salt") == 0) ||
	  (strcmp(tracer->name, "temp") == 0))
        open->bcond_tra[tn] = UPSTRM | FILEIN;
      else
        open->bcond_tra[tn] = NOGRAD;
    }
  }
}

/* END auto_params_roam_post2()                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to overwrite input parameters to ROAM values              */
/*-------------------------------------------------------------------*/
void auto_params_roam_post3(FILE * fp, parameters_t *params)
{
  int n, tn;                    /* Counters */
  open_bdrys_t *open;           /* Pointer to boundary structure */
  char keyword[MAXSTRLEN];      /* Input dummy */
  char buf[MAXSTRLEN];          /* Input dummy */
  int nvf = 1;                  /* 0 : bcond_nor2d = FLATHR|CLAMPD  */
                                /* 1 : bcond_nor2d = FLATHR|CUSTOM  */

  /* Changes to _pre fromulation */
  if (strlen(params->vel_init))
    params->rampf |= ETA_RELAX;
  else
    params->rampf |= (ETA_RELAX|CUSTOM);

  /* Tracer relaxation */
  for (n = 0; n < params->ntr; n++) {
    int tm = params->trinfo_3d[n].m;
    tracer_info_t *tracer = &params->trinfo_3d[n];
    if (strcmp(tracer->name, "temp") == 0 && strlen(params->tdata))
      strcpy(buf, params->tdata);
    else if (strcmp(tracer->name, "salt") == 0 && strlen(params->sdata))
      strcpy(buf, params->sdata);
    else
      continue;
    strcpy(tracer->data, buf);
    if (strlen(params->trrlxn[n]) == 0) {
      strcpy(params->trrlxn[n], buf);
      strcpy(tracer->relax_file, buf);
      params->trrlxdt[n] = 3600.0;
      strcpy(tracer->relax_dt, "1 hour");
      strcpy(params->trrlxtc[n], "20 days");
      strcpy(tracer->r_rate, "20 days");
    }
  }

  /* Tracer reset */
  for (n = 0; n < params->ntr; n++) {
    int tm = params->trinfo_3d[n].m;
    if (strlen(params->trinfo_3d[n].reset_file))
      strcpy(params->trrest[n], params->trinfo_3d[n].reset_file);
    if (strlen(params->trinfo_3d[n].reset_dt))
      tm_scale_to_secs(params->trinfo_3d[n].reset_dt, &params->trrestdt[n]);
  }

  /* Bathymetry limits */
  if (params->bmin < 1.0)
    params->bmin = 1.0;
  if (params->bmax > 20.0 && params->bmin < 2.0)
    params->bmin = 2.0;
  if (params->bmax > 200.0 && params->bmin < 4.0)
    params->bmin = 4.0;

  /* Vertical grid */
  for (n = params->nz - 1; n >= 0; n--) {
    if (params->layers[n] > -LAYERMIN) {
      params->layers[n] = params->layers[params->nz];
      params->nz -= 1;
    }
    if (params->layers[n] >= -LAYERMIN)
      break;
  }
  params->hmin = min(0.1, 0.07 * (params->layers[params->nz] -
				  params->layers[params->nz - 1]));

  /*-----------------------------------------------------------------*/
  /* Open boundaries */
  for (n = 0; n < params->nobc; n++) {
    open = params->open[n];
    if (nvf == 0)
      open->bcond_nor2d = FLATHR|CLAMPD;
    else {
      open->bcond_nor2d = FLATHR|CUSTOM;
      if (open->type & U1BDRY) 
	sprintf(open->cusname_u1av, "uv_to_u1av %s", params->vdata);
      /*sprintf(open->cusname_u1av, "hdstd_to_u1av %s", params->vdata);*/
      else if (open->type & U2BDRY)
	sprintf(open->cusname_u2av, "uv_to_u2av %s", params->vdata);
      /*sprintf(open->cusname_u2av, "hdstd_to_u2av %s", params->vdata);*/
    }
    open->bcond_tan = GRAVTY;
    open->bcond_tan2d = GRAVTY;
    if (strlen(params->odata))
      strcpy(open->tsfn, params->odata);
    if(prm_read_char(fp, "TIDE_CSR_ORTHOWEIGHTS", params->orthoweights)) {
      if (prm_read_char(fp, "TIDE_CSR_CON_DIR", params->nodal_dir))
	open->bcond_ele = FLATHE|FILEIN|TIDALH|GRAVTY;

    } else
      open->bcond_ele = FLATHE|FILEIN|GRAVTY;

    open->sponge_zone_h = 8;
    open->stagger = INFACE;
    open->bathycon = 1;

    for (tn = 0; tn < params->ntr; tn++) {
      tracer_info_t *tracer = &params->trinfo_3d[tn];
      if ((strcmp(tracer->name, "salt") == 0) ||
	  (strcmp(tracer->name, "temp") == 0))
        open->bcond_tra[tn] = UPSTRM | FILEIN;
      else
        open->bcond_tra[tn] = NOGRAD;
    }
  }
}

/* END auto_params_roam_post3()                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to overwrite input parameters to ROAM values              */
/*-------------------------------------------------------------------*/
void auto_params_roam_post4(FILE * fp, parameters_t *params)
{
  int n, tn;                    /* Counters */
  open_bdrys_t *open;           /* Pointer to boundary structure */
  char keyword[MAXSTRLEN];      /* Input dummy */
  char buf[MAXSTRLEN];          /* Input dummy */
  int nvf = 0;                  /* 0 : bcond_nor = CUSTOM            */
                                /* 1 : bcond_nor2d = CUSTOM          */

  /* Changes to _pre fromulation */
  if (strlen(params->vel_init))
    params->rampf |= ETA_RELAX;
  else
    params->rampf |= (ETA_RELAX|CUSTOM);

  /* Tracer relaxation */
  for (n = 0; n < params->ntr; n++) {
    int tm = params->trinfo_3d[n].m;
    tracer_info_t *tracer = &params->trinfo_3d[n];
    if (strcmp(tracer->name, "temp") == 0 && strlen(params->tdata))
      strcpy(buf, params->tdata);
    else if (strcmp(tracer->name, "salt") == 0 && strlen(params->sdata))
      strcpy(buf, params->sdata);
    else
      continue;
    strcpy(tracer->data, buf);
    if (strlen(params->trrlxn[n]) == 0) {
      strcpy(params->trrlxn[n], buf);
      strcpy(tracer->relax_file, buf);
      params->trrlxdt[n] = 3600.0;
      strcpy(tracer->relax_dt, "1 hour");
      strcpy(params->trrlxtc[n], "20 days");
      strcpy(tracer->r_rate, "20 days");
    }
  }

  /* Tracer reset */
  for (n = 0; n < params->ntr; n++) {
    int tm = params->trinfo_3d[n].m;
    if (strlen(params->trinfo_3d[n].reset_file))
      strcpy(params->trrest[n], params->trinfo_3d[n].reset_file);
    if (strlen(params->trinfo_3d[n].reset_dt))
      tm_scale_to_secs(params->trinfo_3d[n].reset_dt, &params->trrestdt[n]);
  }

  /* Bathymetry limits */
  if (params->bmin < 1.0)
    params->bmin = 1.0;
  if (params->bmax > 20.0 && params->bmin < 2.0)
    params->bmin = 2.0;
  if (params->bmax > 200.0 && params->bmin < 4.0)
    params->bmin = 4.0;

  /* Vertical grid */
  for (n = params->nz - 1; n >= 0; n--) {
    if (params->layers[n] > -LAYERMIN) {
      params->layers[n] = params->layers[params->nz];
      params->nz -= 1;
    }
    if (params->layers[n] >= -LAYERMIN)
      break;
  }
  params->hmin = min(0.1, 0.07 * (params->layers[params->nz] -
				  params->layers[params->nz - 1]));
 
  /*-----------------------------------------------------------------*/
  /* Open boundaries */
  for (n = 0; n < params->nobc; n++) {
    open = params->open[n];

    if (strlen(open->bflow)) continue;

    if (nvf == 0) {
      open->bcond_nor = CUSTOM;
      if (open->type & U1BDRY) 
	sprintf(open->cusname_u1, "uv_to_u1 %s", params->vdata);
      else if (open->type & U2BDRY)
	sprintf(open->cusname_u2, "uv_to_u2 %s", params->vdata);
      open->bcond_nor2d = VERTIN;
      
      open->bcond_tan = CUSTOM;
      if (open->type & U1BDRY) 
	sprintf(open->cusname_u2, "uv_to_u2 %s", params->vdata);
      else if (open->type & U2BDRY)
	sprintf(open->cusname_u1, "uv_to_u1 %s", params->vdata);
      open->bcond_tan2d = VERTIN;
    } else {
      open->bcond_nor = NOGRAD;
      open->bcond_nor2d = CUSTOM;
      if (open->type & U1BDRY) 
	sprintf(open->cusname_u1av, "uv_to_u1av %s", params->vdata);
      else if (open->type & U2BDRY)
	sprintf(open->cusname_u2av, "uv_to_u2av %s", params->vdata);

      open->bcond_tan = NOGRAD;
      open->bcond_tan2d = CUSTOM;
      if (open->type & U1BDRY) 
	sprintf(open->cusname_u2av, "uv_to_u2av %s", params->vdata);
      else if (open->type & U2BDRY)
	sprintf(open->cusname_u1av, "uv_to_u1av %s", params->vdata);
    }
    if (strlen(params->odata))
      strcpy(open->tsfn, params->odata);
    if(prm_read_char(fp, "TIDE_CSR_ORTHOWEIGHTS", params->orthoweights)) {
      if (prm_read_char(fp, "TIDE_CSR_CON_DIR", params->nodal_dir))
	open->bcond_ele = NOTHIN|TIDALH;
    } else
      open->bcond_ele = NOTHIN;
    /*
    open->relax_zone_nor = 8;
    open->relax_zone_tan = 8;
    */
    open->adjust_flux = -1.2;  /* 10 times dt2d */
    open->sponge_zone_h = 8;
    /*open->sponge_f = 5;*/
    if(strlen(params->patm)) open->inverse_barometer = 1;
    open->stagger = OUTFACE;

    for (tn = 0; tn < params->ntr; tn++) {
      tracer_info_t *tracer = &params->trinfo_3d[tn];
      if ((strcmp(tracer->name, "salt") == 0) ||
	  (strcmp(tracer->name, "temp") == 0))
        open->bcond_tra[tn] = UPSTRM | FILEIN;
      else
        open->bcond_tra[tn] = NOGRAD;
    }
  }
}

/* END auto_params_roam_post4()                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to overwrite input parameters to ROAM values.             */
/* Included Nov 2013. Uses sponge factor = 5, and FILEIN for         */
/* boundary eta.                                                     */
/*-------------------------------------------------------------------*/
void auto_params_roam_post5(FILE * fp, parameters_t *params)
{
  int n, tn;                    /* Counters */
  open_bdrys_t *open;           /* Pointer to boundary structure */
  char keyword[MAXSTRLEN];      /* Input dummy */
  char buf[MAXSTRLEN];          /* Input dummy */
  int nvf = 0;                  /* 0 : bcond_nor = CUSTOM            */
                                /* 1 : bcond_nor2d = CUSTOM          */

  /* Changes to _pre fromulation */
  if (!strlen(params->vel_init))
    params->rampf |= CUSTOM;

  /* Tracer relaxation */
  for (n = 0; n < params->ntr; n++) {
    int tm = params->trinfo_3d[n].m;
    tracer_info_t *tracer = &params->trinfo_3d[n];

    if (strcmp(tracer->name, "temp") == 0) {
      if (strlen(params->tdata)) strcpy(buf, params->tdata);
      strcpy(tracer->i_rule, "nn_sibson");
      if (params->sharp_pyc) {
	strcpy(tracer->tag, "hipass_vert:1");
	tracer->flag |= V_HP;
	tracer->scale = 1.0;
      }
    } else if (strcmp(tracer->name, "salt") == 0) {
      if (strlen(params->sdata)) strcpy(buf, params->sdata);
      strcpy(tracer->i_rule, "nn_sibson");
      if (params->sharp_pyc) {
	strcpy(tracer->tag, "hipass_vert:1");
	tracer->flag |= V_HP;
	tracer->scale = 1.0;
      }
    } else
      continue;
    strcpy(tracer->data, buf);
    if (strlen(params->trrlxn[n]) == 0) {
      strcpy(params->trrlxn[n], buf);
      strcpy(tracer->relax_file, buf);
      params->trrlxdt[n] = 3600.0;
      strcpy(tracer->relax_dt, "1 hour");
      strcpy(params->trrlxtc[n], "20 days");
      strcpy(tracer->r_rate, "20 days");
      if (params->robust == 4 || params->robust == 8) {
	strcpy(params->trrlxtc[n], "temporal -4 1 hour 0 20 days");
	strcpy(tracer->r_rate, "temporal -4 1 hour 0 20 days");
      }
    }
  }

  /* Tracer reset */
  for (n = 0; n < params->ntr; n++) {
    int tm = params->trinfo_3d[n].m;
    if (strlen(params->trinfo_3d[n].reset_file))
      strcpy(params->trrest[n], params->trinfo_3d[n].reset_file);
    if (strlen(params->trinfo_3d[n].reset_dt))
      tm_scale_to_secs(params->trinfo_3d[n].reset_dt, &params->trrestdt[n]);
  }

  /* Bathymetry limits */
  if (params->bmin < 1.0)
    params->bmin = 1.0;
  if (params->bmax > 20.0 && params->bmin < 2.0)
    params->bmin = 2.0;
  if (params->bmax > 200.0 && params->bmin < 4.0)
    params->bmin = 4.0;

  /* Vertical grid */
  for (n = params->nz - 1; n >= 0; n--) {
    if (params->layers[n] > -LAYERMIN) {
      params->layers[n] = params->layers[params->nz];
      params->nz -= 1;
    }
    if (params->layers[n] >= -LAYERMIN)
      break;
  }
  params->hmin = min(0.1, 0.07 * (params->layers[params->nz] -
				  params->layers[params->nz - 1]));
 
  /*-----------------------------------------------------------------*/
  /* Open boundaries */
  for (n = 0; n < params->nobc; n++) {
    open = params->open[n];

    if (strlen(open->bflow)) continue;

    if (nvf == 0) {
      /*
      sprintf(open->bstd[0], "NEST1WAY %s %s %s %s %s\n",
	      params->odata, params->vdata, params->vdata);
      */
      open->bcond_nor = CUSTOM;
      if (open->type & U1BDRY) 
	sprintf(open->cusname_u1, "uv_to_u1 %s", params->vdata);
      else if (open->type & U2BDRY)
	sprintf(open->cusname_u2, "uv_to_u2 %s", params->vdata);
      open->bcond_nor2d = VERTIN;
      
      open->bcond_tan = CUSTOM;
      if (open->type & U1BDRY) 
	sprintf(open->cusname_u2, "uv_to_u2 %s", params->vdata);
      else if (open->type & U2BDRY)
	sprintf(open->cusname_u1, "uv_to_u1 %s", params->vdata);
      open->bcond_tan2d = VERTIN;
    } else {
      open->bcond_nor = NOGRAD;
      open->bcond_nor2d = CUSTOM;
      if (open->type & U1BDRY) 
	sprintf(open->cusname_u1av, "uv_to_u1av %s", params->vdata);
      else if (open->type & U2BDRY)
	sprintf(open->cusname_u2av, "uv_to_u2av %s", params->vdata);

      open->bcond_tan = NOGRAD;
      open->bcond_tan2d = CUSTOM;
      if (open->type & U1BDRY) 
	sprintf(open->cusname_u2av, "uv_to_u2av %s", params->vdata);
      else if (open->type & U2BDRY)
	sprintf(open->cusname_u1av, "uv_to_u1av %s", params->vdata);
    }
    if (strlen(params->odata))
      strcpy(open->tsfn, params->odata);
    if(prm_read_char(fp, "TIDE_CSR_ORTHOWEIGHTS", params->orthoweights)) {
      if (prm_read_char(fp, "TIDE_CSR_CON_DIR", params->nodal_dir))
	open->bcond_ele = NOTHIN|TIDALH|FILEIN;
    } else
      open->bcond_ele = NOTHIN;
    /*
    open->relax_zone_nor = 8;
    open->relax_zone_tan = 8;
    */
    if (params->robust <= 1) {
      open->adjust_flux = -1.2;  /* 10 times dt2d */
      open->adjust_flux_s = 1.5;
    } else if (params->robust > 1 && params->robust <= 3) {
      open->adjust_flux = 1.8;
    } else {
      open->adjust_flux = -1.2;  /* 10 times dt2d */
    }
    /*open->adjust_flux = -1.2;*/  /* 10 times dt2d */
    if(params->robust < 3 || params->robust >= 4) {
      open->sponge_zone_h = 40e3;
      open->sponge_f = 5;
    }
    if(strlen(params->patm)) open->inverse_barometer = 1;
    open->stagger = OUTFACE;
    open->options = OP_OBCCM;

    for (tn = 0; tn < params->ntr; tn++) {
      tracer_info_t *tracer = &params->trinfo_3d[tn];
      if ((strcmp(tracer->name, "salt") == 0) ||
	  (strcmp(tracer->name, "temp") == 0))
        open->bcond_tra[tn] = TRCONC | FILEIN;
      else
        open->bcond_tra[tn] = NOGRAD;
    }
    open->bgz = laux;
  }
}

/* END auto_params_roam_post5()                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to overwrite input parameters to ROAM values.             */
/* Included Nov 2019. Uses TPXO tide.                                */
/*-------------------------------------------------------------------*/
void auto_params_roam_post6(FILE * fp, parameters_t *params)
{
  int n, tn;                    /* Counters */
  open_bdrys_t *open;           /* Pointer to boundary structure */
  char keyword[MAXSTRLEN];      /* Input dummy */
  char buf[MAXSTRLEN];          /* Input dummy */
  int nvf = 0;                  /* 0 : bcond_nor = CUSTOM            */
                                /* 1 : bcond_nor2d = CUSTOM          */

  /* Changes to _pre fromulation */
  if (!strlen(params->vel_init))
    params->rampf |= CUSTOM;

  /* Tracer relaxation */
  for (n = 0; n < params->ntr; n++) {
    int tm = params->trinfo_3d[n].m;
    tracer_info_t *tracer = &params->trinfo_3d[n];

    if (strcmp(tracer->name, "temp") == 0) {
      if (strlen(params->tdata)) strcpy(buf, params->tdata);
      strcpy(tracer->i_rule, "nn_sibson");
      if (params->sharp_pyc) {
	strcpy(tracer->tag, "hipass_vert:1");
	tracer->flag |= V_HP;
	tracer->scale = 1.0;
      }
    } else if (strcmp(tracer->name, "salt") == 0) {
      if (strlen(params->sdata)) strcpy(buf, params->sdata);
      strcpy(tracer->i_rule, "nn_sibson");
      if (params->sharp_pyc) {
	strcpy(tracer->tag, "hipass_vert:1");
	tracer->flag |= V_HP;
	tracer->scale = 1.0;
      }
    } else
      continue;
    strcpy(tracer->data, buf);
    if (strlen(params->trrlxn[n]) == 0) {
      strcpy(params->trrlxn[n], buf);
      strcpy(tracer->relax_file, buf);
      params->trrlxdt[n] = 3600.0;
      strcpy(tracer->relax_dt, "1 hour");
      strcpy(params->trrlxtc[n], "20 days");
      strcpy(tracer->r_rate, "20 days");
      if (params->robust == 4 || params->robust == 8) {
	strcpy(params->trrlxtc[n], "temporal -4 1 hour 0 20 days");
	strcpy(tracer->r_rate, "temporal -4 1 hour 0 20 days");
      }
    }
  }

  /* Tracer reset */
  for (n = 0; n < params->ntr; n++) {
    int tm = params->trinfo_3d[n].m;
    if (strlen(params->trinfo_3d[n].reset_file))
      strcpy(params->trrest[n], params->trinfo_3d[n].reset_file);
    if (strlen(params->trinfo_3d[n].reset_dt))
      tm_scale_to_secs(params->trinfo_3d[n].reset_dt, &params->trrestdt[n]);
  }

  /* Bathymetry limits */
  if (params->bmin < 1.0)
    params->bmin = 1.0;
  if (params->bmax > 20.0 && params->bmin < 2.0)
    params->bmin = 2.0;
  if (params->bmax > 200.0 && params->bmin < 4.0)
    params->bmin = 4.0;

  /* Vertical grid */
  for (n = params->nz - 1; n >= 0; n--) {
    if (params->layers[n] > -LAYERMIN) {
      params->layers[n] = params->layers[params->nz];
      params->nz -= 1;
    }
    if (params->layers[n] >= -LAYERMIN)
      break;
  }
  params->hmin = min(0.1, 0.07 * (params->layers[params->nz] -
				  params->layers[params->nz - 1]));
 
  /*-----------------------------------------------------------------*/
  /* Open boundaries */
  for (n = 0; n < params->nobc; n++) {
    open = params->open[n];

    if (strlen(open->bflow)) continue;

    if (nvf == 0) {
      /*
      sprintf(open->bstd[0], "NEST1WAY %s %s %s %s %s\n",
	      params->odata, params->vdata, params->vdata);
      */
      open->bcond_nor = CUSTOM;
      if (open->type & U1BDRY) 
	sprintf(open->cusname_u1, "uv_to_u1 %s", params->vdata);
      else if (open->type & U2BDRY)
	sprintf(open->cusname_u2, "uv_to_u2 %s", params->vdata);
      open->bcond_nor2d = VERTIN;
      
      open->bcond_tan = CUSTOM;
      if (open->type & U1BDRY) 
	sprintf(open->cusname_u2, "uv_to_u2 %s", params->vdata);
      else if (open->type & U2BDRY)
	sprintf(open->cusname_u1, "uv_to_u1 %s", params->vdata);
      open->bcond_tan2d = VERTIN;
    } else {
      open->bcond_nor = NOGRAD;
      open->bcond_nor2d = CUSTOM;
      if (open->type & U1BDRY) 
	sprintf(open->cusname_u1av, "uv_to_u1av %s", params->vdata);
      else if (open->type & U2BDRY)
	sprintf(open->cusname_u2av, "uv_to_u2av %s", params->vdata);

      open->bcond_tan = NOGRAD;
      open->bcond_tan2d = CUSTOM;
      if (open->type & U1BDRY) 
	sprintf(open->cusname_u2av, "uv_to_u2av %s", params->vdata);
      else if (open->type & U2BDRY)
	sprintf(open->cusname_u1av, "uv_to_u1av %s", params->vdata);
    }
    if (strlen(params->odata))
      strcpy(open->tsfn, params->odata);
    if(prm_read_char(fp, "TIDE_CONSTITUENTS", params->tide_con_file)) {
      if (prm_read_char(fp, "TIDE_CSR_CON_DIR", params->nodal_dir)) {
	open->bcond_ele = NOTHIN|TIDALC|FILEIN;
	strcpy(open->tide_con, "M2 S2 N2 K2 K1 O1 P1 Q1");
      }
    } else
      open->bcond_ele = NOTHIN;
    /*
    open->relax_zone_nor = 8;
    open->relax_zone_tan = 8;
    */
    if (params->robust <= 1) {
      open->adjust_flux = -1.2;  /* 10 times dt2d */
      open->adjust_flux_s = 1.5;
    } else if (params->robust > 1 && params->robust <= 3) {
      open->adjust_flux = 1.8;
    } else {
      open->adjust_flux = -1.2;  /* 10 times dt2d */
    }
    /*open->adjust_flux = -1.2;*/  /* 10 times dt2d */
    if(params->robust < 3 || params->robust > 4) {
      open->sponge_zone_h = 8;
      open->sponge_f = 5;
    }
    if(strlen(params->patm)) open->inverse_barometer = 1;
    open->stagger = OUTFACE;
    open->options = OP_OBCCM;

    for (tn = 0; tn < params->ntr; tn++) {
      tracer_info_t *tracer = &params->trinfo_3d[tn];
      if ((strcmp(tracer->name, "salt") == 0) ||
	  (strcmp(tracer->name, "temp") == 0))
        open->bcond_tra[tn] = TRCONC | FILEIN;
      else
        open->bcond_tra[tn] = NOGRAD;
    }
    open->bgz = laux;
  }
}

/* END auto_params_roam_post6()                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to overwrite input parameters to ROAM values.             */
/* Included Sep 2020. Uses new ROBUST and OBCs.                      */
/*-------------------------------------------------------------------*/
void auto_params_roam_post7(FILE * fp, parameters_t *params)
{
  geometry_t *geom = master->geom;
  int n, tn;                    /* Counters */
  open_bdrys_t *open;           /* Pointer to boundary structure */
  char keyword[MAXSTRLEN];      /* Input dummy */
  char buf[MAXSTRLEN];          /* Input dummy */
  int nvf = 0;                  /* 0 : bcond_nor = CUSTOM            */

                                /* 1 : bcond_nor2d = CUSTOM          */

  /* Changes to _pre fromulation */
  if (!strlen(params->vel_init))
    params->rampf |= CUSTOM;

  /* Tracer relaxation */
  for (n = 0; n < params->ntr; n++) {
    int tm = params->trinfo_3d[n].m;
    tracer_info_t *tracer = &params->trinfo_3d[n];

    if (strcmp(tracer->name, "temp") == 0) {
      if (strlen(params->tdata)) strcpy(buf, params->tdata);
      if (params->sharp_pyc) {
	strcpy(tracer->tag, "hipass_vert:1");
	tracer->flag |= V_HP;
	tracer->scale = 1.0;
      }
    } else if (strcmp(tracer->name, "salt") == 0) {
      if (strlen(params->sdata)) strcpy(buf, params->sdata);
      if (params->sharp_pyc) {
	strcpy(tracer->tag, "hipass_vert:1");
	tracer->flag |= V_HP;
	tracer->scale = 1.0;
      }
    } else
      continue;
    strcpy(tracer->data, buf);
    if (strlen(params->trrlxn[n]) == 0) {
      strcpy(params->trrlxn[n], buf);
      strcpy(tracer->relax_file, buf);
      params->trrlxdt[n] = 3600.0;
      strcpy(tracer->relax_dt, "1 hour");
      strcpy(params->trrlxtc[n], "20 days");
      strcpy(tracer->r_rate, "20 days");
      if (params->robust == 3 || params->robust == 4) {
	strcpy(params->trrlxtc[n], "temporal -4 1 hour 0 20 days");
	strcpy(tracer->r_rate, "temporal -4 1 hour 0 20 days");
      }
    }
  }

  /* Tracer reset */
  for (n = 0; n < params->ntr; n++) {
    int tm = params->trinfo_3d[n].m;
    if (strlen(params->trinfo_3d[n].reset_file))
      strcpy(params->trrest[n], params->trinfo_3d[n].reset_file);
    if (strlen(params->trinfo_3d[n].reset_dt))
      tm_scale_to_secs(params->trinfo_3d[n].reset_dt, &params->trrestdt[n]);
  }

  /* Bathymetry limits */
  if (params->bmin < 1.0)
    params->bmin = 1.0;
  if (params->bmax > 20.0 && params->bmin < 2.0)
    params->bmin = 2.0;
  if (params->bmax > 200.0 && params->bmin < 4.0)
    params->bmin = 4.0;

  /* Vertical grid */
  for (n = params->nz - 1; n >= 0; n--) {
    if (params->layers[n] > -LAYERMIN) {
      params->layers[n] = params->layers[params->nz];
      params->nz -= 1;
    }
    if (params->layers[n] >= -LAYERMIN)
      break;
  }
  params->hmin = min(0.1, 0.07 * (params->layers[params->nz] -
				  params->layers[params->nz - 1]));
 
  /*-----------------------------------------------------------------*/
  /* Open boundaries */
  for (n = 0; n < params->nobc; n++) {
    open = params->open[n];

    if (strlen(open->bflow)) continue;

    if (nvf == 0) {
      open->bcond_nor = CUSTOM;
      if (open->type & U1BDRY) 
	sprintf(open->cusname_u1, "uv_to_u1 %s", params->vdata);
      else if (open->type & U2BDRY)
	sprintf(open->cusname_u2, "uv_to_u2 %s", params->vdata);
      open->bcond_nor2d = VERTIN|TIDALC;

      open->bcond_tan = CUSTOM;
      if (open->type & U1BDRY) 
	sprintf(open->cusname_u2, "uv_to_u2 %s", params->vdata);
      else if (open->type & U2BDRY)
	sprintf(open->cusname_u1, "uv_to_u1 %s", params->vdata);
      open->bcond_tan2d = VERTIN|TIDALC;
    } else {
      open->bcond_nor = NOGRAD;
      open->bcond_nor2d = CUSTOM|TIDALC;
      if (open->type & U1BDRY) 
	sprintf(open->cusname_u1av, "uv_to_u1av %s", params->vdata);
      else if (open->type & U2BDRY)
	sprintf(open->cusname_u2av, "uv_to_u2av %s", params->vdata);

      open->bcond_tan = NOGRAD;
      open->bcond_tan2d = CUSTOM|TIDALC;
      if (open->type & U1BDRY) 
	sprintf(open->cusname_u2av, "uv_to_u2av %s", params->vdata);
      else if (open->type & U2BDRY)
	sprintf(open->cusname_u1av, "uv_to_u1av %s", params->vdata);
    }
    if (strlen(params->odata))
      strcpy(open->tsfn, params->odata);
    if(prm_read_char(fp, "TIDE_CONSTITUENTS", params->tide_con_file)) {
      if (prm_read_char(fp, "TIDE_CSR_CON_DIR", params->nodal_dir)) {
	open->bcond_ele = NOTHIN|TIDALC|FILEIN;
	strcpy(open->tide_con, "M2 S2 N2 K2 K1 O1 P1 Q1");
      }
    } else
      open->bcond_ele = NOTHIN;
    /*
    open->relax_zone_nor = 8;
    open->relax_zone_tan = 8;
    */

    open->options |= OP_OBCCM;          /* Corner means              */
    params->tidef |= TD_TRAN;           /* Trasports tide file       */
    if (params->robust == 1 || params->robust == 3 || params->robust == 5) {
     open->adjust_flux_s = -10.0;       /* Tidal relaxation          */
     open->adjust_flux = 1.2;           /* Constant relaxation       */
     /*open->adjust_flux = -1.2;*/        /* Default relaxation        */
     open->options |= (OP_FAT|OP_FAS);  /* Scale to minimum OBC CFL  */
    } else {
      open->adjust_flux = -1.2;         /* Default relaxation        */
      open->options |= OP_FAS;          /* Scale to minimum OBC CFL  */
    }

    /*
    if(params->robust < 3) {
      open->sponge_zone_h = 50000;
      open->sponge_f = 5;
    }
    */
    if(strlen(params->patm)) open->inverse_barometer = 1;
    open->stagger = OUTFACE;

    for (tn = 0; tn < params->ntr; tn++) {
      tracer_info_t *tracer = &params->trinfo_3d[tn];
      if ((strcmp(tracer->name, "salt") == 0) ||
	  (strcmp(tracer->name, "temp") == 0))
        open->bcond_tra[tn] = UPSTRM | FILEIN;
      else
        open->bcond_tra[tn] = NOGRAD;
    }
  }
}

/* END auto_params_roam_post7()                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to overwrite input parameters to RECOM values             */
/*-------------------------------------------------------------------*/
void auto_params_recom_post1(FILE * fp, parameters_t *params)
{
  int n, m;                     /* Counters */
  open_bdrys_t *open;           /* Pointer to boundary structure */
  char keyword[MAXSTRLEN];      /* Input dummy */
  char buf[MAXSTRLEN];          /* Input dummy */

  /* Changes to _pre fromulation */
  params->rampf |= CUSTOM;

  /* Tracer reset */
  for (n = 0; n < params->ntr; n++) {
    int tm = params->trinfo_3d[n].m;
    if (strlen(params->trinfo_3d[n].reset_file))
      strcpy(params->trrest[n], params->trinfo_3d[n].reset_file);
    if (strlen(params->trinfo_3d[n].reset_dt))
      tm_scale_to_secs(params->trinfo_3d[n].reset_dt, &params->trrestdt[n]);
  }

  /* Vertical grid */
  for (n = params->nz - 1; n >= 0; n--) {
    if (params->layers[n] > -LAYERMIN) {
      params->layers[n] = params->layers[params->nz];
      params->nz -= 1;
    }
    if (params->layers[n] >= -LAYERMIN)
      break;
  }
  params->hmin = min(0.1, 0.07 * (params->layers[params->nz] -
				  params->layers[params->nz - 1]));
 
  /*-----------------------------------------------------------------*/
  /* Open boundaries */
  for (n = 0; n < params->nobc; n++) {
    open = params->open[n];

    if (strlen(open->bflow)) {
      sprintf(buf, "RIVER %s %s",open->bflow, params->odata);
      std_bdry_init(open, buf, params->grid_dt, 
		    params->grid_dt / (double)params->iratio,
		    params->bdrypath, params->trinfo_3d);
    } else {
      char vnor[MAXSTRLEN], vtan[MAXSTRLEN];
      char *fields[MAXSTRLEN * MAXNUMARGS];
      strcpy(buf, params->vdata);
      m = parseline(buf, fields, MAXNUMARGS);
      if (m == 2) {
	strcpy(vnor, fields[0]);
	strcpy(vtan, fields[1]);
      } else if (m == 1) {
	strcpy(vnor, fields[0]);
	strcpy(vtan, fields[0]);
      }
      if (m >= 3)
	sprintf(buf, "NEST1WAY|STD %s",params->vdata);
      else
	sprintf(buf, "NEST1WAY|STD %s %s %s",params->odata, vnor, vtan);
      open->adjust_flux = -1.2;
      std_bdry_init(open, buf, params->grid_dt, 
		    params->grid_dt / (double)params->iratio,
		    params->bdrypath, params->trinfo_3d);
      if(prm_read_char(fp, "TIDE_CSR_ORTHOWEIGHTS", params->orthoweights)) {
	if (prm_read_char(fp, "TIDE_CSR_CON_DIR", params->nodal_dir))
	  open->bcond_ele |= TIDALH;
      }
    }
    open->bgz = laux;
    /*open->sponge_zone_h = 8;*/
    /*open->sponge_f = 5;*/
  }
}

/* END auto_params_recom_post1()                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to overwrite input parameters to RECOM values             */
/*-------------------------------------------------------------------*/
void auto_params_recom_post2(FILE * fp, parameters_t *params)
{
  int n, m;                     /* Counters */
  open_bdrys_t *open;           /* Pointer to boundary structure */
  char keyword[MAXSTRLEN];      /* Input dummy */
  char buf[MAXSTRLEN];          /* Input dummy */

  /* Changes to _pre fromulation */
  if (!strlen(params->vel_init))
    params->rampf |= CUSTOM;

  /* Tracer relaxation */
  for (n = 0; n < params->ntr; n++) {
    int tm = params->trinfo_3d[n].m;
    tracer_info_t *tracer = &params->trinfo_3d[n];
    if (strcmp(tracer->name, "temp") == 0 && strlen(params->tdata))
      strcpy(buf, params->tdata);
    else if (strcmp(tracer->name, "salt") == 0 && strlen(params->sdata))
      strcpy(buf, params->sdata);
    else
      continue;
    if (params->robust == 7 || params->robust == 8) {
      if (strlen(params->trrlxn[n]) == 0) {
	strcpy(params->trrlxn[n], buf);
	strcpy(tracer->relax_file, buf);
	params->trrlxdt[n] = 3600.0;
	strcpy(tracer->relax_dt, "1 hour");
	strcpy(params->trrlxtc[n], "20 days");
	strcpy(tracer->r_rate, "20 days");
	sprintf(params->trrlxtc[n], "temporal %f 1 day %f 20 days",params->t/86400, params->t/86400 + 4);
	sprintf(tracer->r_rate, "temporal %f 1 day %f 20 days",params->t/86400, params->t/86400 + 4);
      }
    }
  }

  /* Tracer reset */
  for (n = 0; n < params->ntr; n++) {
    int tm = params->trinfo_3d[n].m;
    if (strlen(params->trinfo_3d[n].reset_file))
      strcpy(params->trrest[n], params->trinfo_3d[n].reset_file);
    if (strlen(params->trinfo_3d[n].reset_dt))
      tm_scale_to_secs(params->trinfo_3d[n].reset_dt, &params->trrestdt[n]);
  }

  /* Bathymetry limits */
 if (params->robust >= 6) {
   if (params->bmin < 1.0)
     params->bmin = 1.0;
   if (params->bmax > 20.0 && params->bmin < 2.0)
     params->bmin = 2.0;
   if (params->bmax > 200.0 && params->bmin < 4.0)
     params->bmin = 4.0;
 }

  /* Vertical grid */
  if (params->robust <= 6) {
    double *nlayers = d_alloc_1d(params->nz + 3);
    for (n = 0; n <= params->nz; n++)
      nlayers[n] = params->layers[n];
    nlayers[n] = -params->layers[n - 2];
    nlayers[n + 1] = -params->layers[n - 3];
    params->nz += 2;
    d_free_1d(params->layers);
    params->layers = d_alloc_1d(params->nz + 1);
    for (n = 0; n <= params->nz; n++)
      params->layers[n] = nlayers[n];
    d_free_1d(nlayers);
  } else {
    for (n = params->nz - 1; n >= 0; n--) {
      if (params->layers[n] > -LAYERMIN) {
	params->layers[n] = params->layers[params->nz];
      params->nz -= 1;
      }
      if (params->layers[n] >= -LAYERMIN)
	break;
    }
  }
  params->hmin = min(0.1, 0.07 * (params->layers[params->nz] -
				  params->layers[params->nz - 1]));

 
  /*-----------------------------------------------------------------*/
  /* Open boundaries */
  for (n = 0; n < params->nobc; n++) {
    int riv = 0;
    open = params->open[n];

    if (strlen(open->bflow)) {
      sprintf(buf, "RIVER %s %s",open->bflow, params->odata);
      std_bdry_init(open, buf, params->grid_dt, 
		    params->grid_dt / (double)params->iratio,
		    params->bdrypath, params->trinfo_3d);
      riv = 1;
    } else {
      char vnor[MAXSTRLEN], vtan[MAXSTRLEN];
      char *fields[MAXSTRLEN * MAXNUMARGS];
      strcpy(buf, params->vdata);
      m = parseline(buf, fields, MAXNUMARGS);
      if (m == 2) {
	strcpy(vnor, fields[0]);
	strcpy(vtan, fields[1]);
      } else if (m == 1) {
	strcpy(vnor, fields[0]);
	strcpy(vtan, fields[0]);
      }
      if (m >= 3)
	sprintf(buf, "NEST1WAY|STD %s",params->vdata);
      else
	sprintf(buf, "NEST1WAY|STD %s %s %s",params->odata, vnor, vtan);
      open->adjust_flux = -1.2;
      std_bdry_init(open, buf, params->grid_dt, 
		    params->grid_dt / (double)params->iratio,
		    params->bdrypath, params->trinfo_3d);
      if(prm_read_char(fp, "TIDE_CSR_ORTHOWEIGHTS", params->orthoweights)) {
	if (prm_read_char(fp, "TIDE_CSR_CON_DIR", params->nodal_dir))
	  open->bcond_ele |= TIDALH;
      }
      if (params->robust >= 7) {
	for (m = 0; m < params->ntr; m++) {
	  if (open->bcond_tra[m] == (TRCONC|FILEIN)) {
	    open->bcond_tra[m] = (UPSTRM|FILEIN);
	  }
	}
      }
    }
    open->bgz = laux;
    if(!riv && params->robust >= 3) {
      open->sponge_zone_h = 8;
      open->sponge_f = 5;
    }
  }
}

/* END auto_params_recom_post2()                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to overwrite automated parameters for ROAM                */
/*-------------------------------------------------------------------*/
void autoset_roam(parameters_t *params, master_t *master, geometry_t **window)
{
  geometry_t *geom = master->geom;

  int n, cc, c, e, lc;          /* Sparse counters */
  int i, j, k;                  /* Cartesian counters */
  double cif = 2.0;             /* Scaling factor for internal wave speed */
  double sf = 0.4;              /* Safety factor for cfl calculations */
  double hf = 0.05;             /* Factor for horizontal diffusion */
  double sfr = 1.0;             /* Additional ROAM safety factor */
  double eps = 1e5;             /* Large value for cfl calculations */
  double cfl2d, cfl3d;          /* CFL timesteps */
  double hmax;                  /* Maximum horizontal diffusion */
  double lat, lon;              /* Latitude and longitude for Coriolis */
  double etas;                  /* Eta scaling */
  double d1, d2, d3;            /* Dummies */
  double btws, btwsm = 0.0;
  double bcws, bcwsm = 0.0;

  /*-----------------------------------------------------------------*/
  /* Re-calculate the time steps based on conservative ROAM safety   */
  /* factors.                                                        */
  cfl2d = cfl3d = 1e10;
  /* Linear from sf *= 1 for robust = 1, sf *= 0.5 for robust = 10   */
  /* sf *= (1.0 - ((double)(params->robust) - 1.0) / 18.0);*/
  /* Linear from sf *= 1 for robust = 6, sf *= 0.5 for robust = 10   */
  if (params->roammode <= A_ROAM_R1) {
    if (params->robust >= 6)
      sf *= (1.0 - 0.5 * ((double)(params->robust) - 6.0) / 4.0);
  } else if (params->roammode == A_ROAM_R2) {
    if (params->robust == 10)
      sf *= 0.5;
  } else if (params->roammode & (A_ROAM_R3|A_ROAM_R4|A_ROAM_R5)) {
    if (params->robust == 5)
      sf *= 0.5;
  }

  for (cc = 1; cc <= geom->v2_t; cc++) {
    c = lc = geom->w2_t[cc];
    
    if (params->us_type & US_IJ) {
      d1 = 0.0;
      for (j = 1; j <= geom->npe[c]; j++) {
	e = geom->c2e[j][c];
	if (c == geom->e2c[e][0]) {
	  if (geom->h1acell[e])
	    d1 += 1.0 / (geom->h1acell[e] * geom->h1acell[e]);
	}
      }
      d1 = sqrt(d1);
    } else
      d1 = sqrt((double)geom->npe[c] / (2.0 * geom->cellarea[c]));

    bcws = int_wave_spd_cont_m(master, c, cc);
    if(bcws == 0.0) bcws = 1.0;
	
    bcwsm = (bcws > bcwsm) ? bcws : bcwsm;
    btws = sqrt(-params->g * geom->botz[c]);
    btwsm = (btws > btwsm) ? btws : btwsm;
    d2 = 2.0 * cif * bcws + params->velmax;
    d3 = 2.0 * btws + params->velmax2d;
    d3 = (d3 > 0.0 && d1 > 0.0) ? 1.0 / (d1 * d3) : eps;
    if (d3 < cfl2d)
      cfl2d = d3;
    d2 = (d2 > 0.0 && d1 > 0.0) ? 1.0 / (d1 * d2) : eps;
    if (d2 < cfl3d)
      cfl3d = d2;
  }

  if (DEBUG("init_m")) {
    dlog("autoset_roam : auto dt","Maximum advective velocity expected = %3.1f", params->velmax);
    dlog("autoset_roam : auto dt","Internal wave speed scaling factor= %3.1f",cif);
    dlog("autoset_roam : auto dt","Safety factor for cfl calculations = %3.1f",sf);
    dlog("autoset_roam : auto dt","Maximum barotropic wave speed = %3.1f",btwsm);
    dlog("autoset_roam : auto dt","Maximum baroclinic wave speed = %3.1f",bcwsm);
    dlog("autoset_roam : auto dt","cfl3d = %5.1f, cfl2d = %3.1f",cfl3d,cfl2d);
  }

  if (params->bmax > 3000)
    sfr = 0.35;
  else if (params->bmax > 1000 && params->bmax <= 3000)
    sfr = 0.5;
  sf *= sfr;
  if(params->speed > 1)
    sf = (0.9 - sf) * ((double)params->speed - 1.0) / 9.0 + sf;
  if (sf * cfl3d > 300)
    d1 = 30;
  else
    d1 = 5;

  c = (int)cfl3d *sf / d1;
  params->grid_dt = d1 * (double)c;
  c = (int)cfl2d *sf;
  if (c > 1.0) {
    while (c > 1 && (c % 2 != 0 && c % 3 != 0 && c % 5 != 0))
      c--;
    if (c)
      cfl2d = (double)c;
    else
      cfl2d = floor(cfl2d);
  } else {
    c = (int)(cfl2d *10) * sf;
    cfl2d = (double)c / 10.0;
    hd_warn("floor(2D time-step) <= 1 (%f)\n", cfl2d);
  }

  params->iratio = max(1, (int)(params->grid_dt / cfl2d));
  if (params->iratio % 2 != 0)
    params->iratio++;
  params->grid_dt = (double)params->iratio * cfl2d;
  if ( (params->roammode & (A_RECOM_R1|A_RECOM_R2)) && (params->grid_dt < 1) )
    hd_quit("RECOM 3D time step smaller than 1 sec : %f\n", params->grid_dt);
  
  if (DEBUG("init_m")) {
      dlog("autoset_roam : auto dt","Scaled cfl3d = %5.1f",params->grid_dt);
      dlog("autoset_roam : auto dt","Scaled cfl2d = %5.1f",cfl2d);
  }

  /*-----------------------------------------------------------------*/
  /* Horizontal mixing coefficients */
  n = 0;
  d1 = 1e-10;
  hmax = 1e10;
  for (cc = 1; cc <= geom->v2_t; cc++) {
    c = geom->w2_t[cc];

    if (params->us_type & US_IJ) {
      d2 = 0.0;
      for (j = 1; j <= geom->npe[c]; j++) {
	e = geom->c2e[j][c];
	if (c == geom->e2c[e][0]) {
	  /*d1 += geom->h1acell[e];*/
	  if (geom->h1acell[e])
	    d2 += 1.0 / (geom->h1acell[e] * geom->h1acell[e]);
	}
      }
    } else {
      d2 = 1.0 / geom->cellarea[c];
    }
    d1 += geom->cellarea[c];
    d3 = (d2) ? 1.0 / (d2 * 4.0 * params->grid_dt) : 1e10;
    if (d3 < hmax)
      hmax = d3;
    n += 1;
  }
  if (hmax > 5.0) {
    c = (int)hmax / 5;
    hmax = 5.0 * (double)c;
  }
  if (n) {
    int u1khf = 0, u1vhf = 0;
    u1khf = 1; u1vhf = 1;
    d1 /= (double)n;
    d1 = 0.01 * d1 / params->grid_dt;
    if (d1 > 5.0) {
      c = (int)d1 / 5;
      if (u1khf)
	params->u1kh = 5.0 * (double)c;
      if (u1vhf)
	params->u1vh = 5.0 * (double)c;
    } else {
      if (u1khf)
	params->u1kh = d1;
      if (u1vhf)
	params->u1vh = d1;
    }
    /* Set limits */
    if (u1khf) {
      if (params->u1kh > hmax)
        params->u1kh = hmax;
      if (hf * hmax > 5.0 && params->u1kh < hf * hmax) {
        c = (int)(hf * hmax) / 5;
        params->u1kh = 5.0 * (double)c;
      }
    }
    if (u1vhf) {
      /* Linear from u1vh *= 0.5 for robust = 6, u1vh *= 1 for       */
      /* robust = 10                                                 */
      /*
      if (params->robust >= 6)
	      params->u1vh *= (((double)params->robust - 10.0) / 8.0 + 1.0);
      */
      if (params->u1vh > hmax)
        params->u1vh = hmax;
      if (hf * hmax > 5.0 && params->u1vh < hf * hmax) {
        c = (int)(hf * hmax) / 5;
        params->u1vh = 5.0 * (double)c;
      }
    }
  } else
    hd_warn("Horizontal mixing is set to zero\n");
  if(params->smagorinsky > 0.0) {
    if (params->roammode == A_ROAM_R1 && params->robust != 6 && params->robust != 7) {
      params->u1vh *= -1.0;
    }
    if (params->roammode == A_ROAM_R2 && params->robust <=6) {
      params->u1vh *= -1.0;
    }
    if (params->roammode & (A_ROAM_R3|A_ROAM_R4|A_ROAM_R5)) {
      params->u1vh *= -1.0;
    }
    if (params->roammode == A_RECOM_R2 && params->robust != 4) {
      params->u1vh *= -1.0;
    }
    params->u1kh *= -1.0;
  }
  if (params->roammode & (A_ROAM_R3|A_ROAM_R4|A_ROAM_R5)) {
    if (!(params->diff_scale & AUTO)) {
      params->bsue1 = floor(params->bsue1 * fabs(params->u1vh));
      params->bkue1 = floor(params->bkue1 * fabs(params->u1kh));
    } else {
      /* AUTO mixing requires a scaling percentage, but we've        */
      /* already converted u1vh to a fraction in ROAM. Convert back. */
      params->bsue1 *= 100.0;
      params->bkue1 *= 100.0;
    }
    sprintf(params->smag, "%f %f", fabs(params->bsue1) + params->sue1,
	    fabs(params->bkue1) + params->kue1);
  }

  /* Scale the elevation initial condition                           */
  etas = 0.0;
  if (params->roammode & A_ROAM_R5) {
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = geom->w2_t[cc];
      etas += master->eta[c];
    }
    etas /= (double)geom->b2_t;
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = geom->w2_t[cc];
      master->eta[c] -= etas;
    }
  }

  /* Set time-step dependent parameters */
  if (params->roammode & (A_ROAM_R1|A_ROAM_R2|A_ROAM_R3|A_ROAM_R4|A_ROAM_R5|A_RECOM_R2)) {
    int n, nn, c1;
    d1 = params->grid_dt / (double)params->iratio;

    for (n = 0; n < geom->nobc; n++) {
      double btws, f, fg, h;
      double cfl2d = HUGE;
      double mws = HUGE;
      double miws = 0.0;
      double spd = 2.0;
      open_bdrys_t *open = geom->open[n];
      fa_info_t *fas;
      int checkf = 0;

      /* Custom adjustment timescale
      flux_adjust_init(open);
      fas = open->fas;
      fas->tc0 = 1.0 / 10.0;
      fas->tc1 = 1.0 / 1.2;
      fas->dv0 = fabs(open->mindep);
      fas->dv1 = fabs(open->maxdep);
      fas->slope = (fas->tc0 - fas->tc1) / (fas->dv0 - fas->dv1);
      flux_adjust_init(params->open[n]);
      flux_adjust_copy(open, params->open[n]);
      */

      d1 = params->grid_dt / (double)params->iratio;
      params->open[n]->adjust_flux *= d1;
      open->adjust_flux *= d1;
      params->open[n]->adjust_flux_s *= d1;
      open->adjust_flux_s *= d1;

      /* Boundary eta scaling */
      if (params->roammode & A_ROAM_R5)
	sprintf(params->open[n]->scale_e, "-%f", etas);

      /* Distribute to windows */
      for (nn = 1; nn <= params->nwindows; nn++) {
	c1 = geom->owc[n][nn];
	if (c1 != -1) {
	  window[nn]->open[c1]->adjust_flux = open->adjust_flux;
	  window[nn]->open[c1]->adjust_flux_s = open->adjust_flux_s;
	  /*
	  flux_adjust_init(window[nn]->open[c1]);
	  flux_adjust_copy(open, window[nn]->open[c1]);
	  */
	}
      }
    }
  }
}

/* END autoset_roam()                                                */
/*-------------------------------------------------------------------*/
