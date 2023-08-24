/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/inputs/readparam_t.c
 *  
 *  Description:
 *  Routine to read model parameters for '-t' option
 *  from a file.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: readparam_t.c 7276 2022-12-14 05:39:52Z her127 $
 *
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hd.h"

void check_wave_vars(parameters_t *params);

/*-------------------------------------------------------------------*/
/* Routine to set input parameters for transport option.             */
/*-------------------------------------------------------------------*/

parameters_t *params_read_t(fp)
FILE *fp;
{
  parameters_t *params;
  char buf[MAXSTRLEN];
  char buf1[MAXSTRLEN];
  char keyword[MAXSTRLEN];
  int m, n;
  int cc;

  /* Allocate memory for the parameter data structure                */
  params = params_alloc();

  /* Set defaults                                                    */
  params->runmode = TRANS;
  set_default_param(params);

  /* Read in the name of the input file                              */
  params->prmfd = fp;
  prm_set_errfn(hd_quit);
  if (forced_restart) {
    if (prm_read_char(fp, "restart_name", buf))
      strcpy(params->idumpname, buf);
    else
      strcpy(params->idumpname, "restart.nc");
    params->t = get_restart_time(params->idumpname, schedule->units);
    sprintf(params->start_time, "%g days", params->t / 86400);
  } else {
    prm_read_char(fp, "INPUT_FILE", params->idumpname);
    /* Read the start and stop time (used only for diagnostics)        */
    sprintf(keyword, "START_TIME");
    prm_read_char(fp, keyword, params->start_time);
    prm_get_time_in_secs(fp, "START_TIME", &params->t);
  }
  sprintf(keyword, "PROJECTION");
  if (prm_read_char(fp, keyword, params->projection)) {
    strcpy(projection, params->projection);
    ts_set_default_proj_type(projection);
  }

  /* Time                                                            */
  sprintf(keyword, "DT");
  prm_get_time_in_secs(fp, keyword, &params->grid_dt);
  params->iratio = 2;
  sprintf(keyword, "TRATIO");
  prm_read_double(fp, keyword, &params->tratio);
  sprintf(keyword, "OUTPUT_TIMEUNIT");
  if (prm_read_char(fp, keyword, params->output_tunit)) {
    char *p = strstr(params->output_tunit, "since");
    if (p == NULL)
      sprintf(params->output_tunit, "%s %s",
              strdup(params->output_tunit),
              strstr(params->timeunit, "since"));
  } else
    strcpy(params->output_tunit, params->timeunit);

  read_tmode_params(fp, params);

  sprintf(keyword, "STOP_TIME");
  prm_read_char(fp, keyword, params->stop_time);
  sprintf(keyword, "TIMEUNIT");
  prm_read_char(fp, keyword, params->timeunit);
  sprintf(keyword, "OUTPUT_TIMEUNIT");
  if (prm_read_char(fp, keyword, params->output_tunit)) {
    char *p = strstr(params->output_tunit, "since");
    if (p == NULL)
      sprintf(params->output_tunit, "%s %s",
              strdup(params->output_tunit),
              strstr(params->timeunit, "since"));
  } else
    strcpy(params->output_tunit, params->timeunit);

  /* Grid size, layer structure, time unit and length unit are read  */
  /* in dump_open_us().                                              */
  /* Bathymetry is set in create_mesh().                             */

  /* Mandatory parameters                                            */
  prm_set_errfn(quit);
  sprintf(keyword, "HMIN");
  prm_read_double(fp, keyword, &params->hmin);

  /* Grid description */
  prm_set_errfn(hd_silent_warn);
  strcpy(params->prmname, prmname);
  sprintf(keyword, "DESCRIPTION");
  prm_read_char(fp, keyword, params->grid_desc);
  sprintf(keyword, "NAME");
  prm_read_char(fp, keyword, params->grid_name);
  sprintf(keyword, "CODEHEADER");
  prm_read_char(fp, keyword, params->codeheader);
  strcpy(codeheader, params->codeheader);
  sprintf(keyword, "PARAMETERHEADER");
  prm_read_char(fp, keyword, params->parameterheader);
  strcpy(parameterheader, params->parameterheader);
  sprintf(keyword, "SEQUENCE");
  prm_read_char(fp, keyword, params->sequence);
  /* ID number and revision */
  sprintf(keyword, "ID_NUMBER");
  prm_read_char(fp, keyword, params->runnoc);
  params->runno = atof(params->runnoc);
  sprintf(keyword, "REVISION");
  prm_read_char(fp, keyword, params->rev);

  /* Optional parameters */
  prm_read_char(fp, "HISTORY", buf);
  if (contains_token(buf, "NONE") != NULL) {
    params->history = NONE;
  } else {
    params->history = 0;
    if (contains_token(buf, "LOG") != NULL)
      params->history |= HST_LOG;
    if (contains_token(buf, "DIFF") != NULL)
      params->history |= HST_DIF;
    if (contains_token(buf, "RESET") != NULL)
      params->history |= HST_RESET;
  }

  /* Bottom roughness */
  z0_init(params);

  /* Compatibility */
  read_compatible(params, fp);

  /* Advection scheme */
  params->trasc = LAGRANGE;
  sprintf(keyword, "TRA_SCHEME");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "ORIGINAL") == 0)
      params->trasc = ORIGINAL;
    if (strcmp(buf, "ORDER1") == 0)
      params->trasc = ORDER1;
    if (strcmp(buf, "ORDER2") == 0)
      params->trasc = ORDER2;
    if (strcmp(buf, "ORDER4") == 0)
      params->trasc = ORDER4;
    if (strcmp(buf, "QUICKEST") == 0)
      params->trasc = QUICKEST;
    if (strcmp(buf, "QUICKEST|US") == 0)
      params->trasc = QUICKEST|HIORDER;
    if (strcmp(buf, "VANLEER") == 0)
      params->trasc = VANLEER;
    if (strcmp(buf, "LAGRANGE") == 0)
      params->trasc = LAGRANGE;
    if (strcmp(buf, "FFSL") == 0)
      params->trasc = FFSL;
    if (strcmp(buf, "NONE") == 0)
      params->trasc = NONE;
    if (strcmp(buf, "LAGRANGE|VANLEER") == 0) {
      params->trasc = LAGRANGE|VANLEER;
      params->trsplit = 1;
    }
  }
  if (params->trasc & LAGRANGE)
    params->ntr += 1;         /* Volume error                        */
  sprintf(keyword, "ULTIMATE");
  if (prm_read_char(fp, keyword, buf))
    params->ultimate = is_true(buf);
  /* Runge-Kutta stages                                              */
  sprintf(keyword, "RUNGE-KUTTA");
  prm_read_int(fp, keyword, &params->rkstage);
  sprintf(keyword, "ORDER_SL");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "LINEAR") == 0)
      params->osl = L_LINEAR;
    else if (strcmp(buf, "NN_SIBSON") == 0)
      params->osl = L_SIB;
    else if (strcmp(buf, "NN_NON_SIBSON") == 0)
      params->osl = L_NONSIB;
    else if (strcmp(buf, "CUBIC") == 0)
      params->osl = L_CUBIC;
    else if (strcmp(buf, "QUADRATIC") == 0)
      params->osl = L_LSQUAD;
    else if (strcmp(buf, "LINEARLSQ") == 0)
      params->osl = L_LSLIN;
    else if (strcmp(buf, "BILINEAR") == 0)
      params->osl = L_BILIN;
    else if (strcmp(buf, "BAYCENTRIC") == 0)
      params->osl = L_BAYLIN;
  }
  if (params->trasc & FFSL && params->osl & (L_BILIN|L_BAYLIN))
    hd_quit("Cannot use FFSL advection with bilinear or baycentric interpolation, due to streamline tracking from edges.\n");

  /* Sub-stepping method */
  sprintf(keyword, "STABILITY");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "NONE") == 0)
      params->stab = NONE;
    if (strcmp(buf, "ORIGINAL") == 0)
      params->stab = ORIGINAL;
    if (strcmp(buf, "SUB-STEP") == 0)
      params->stab = SUB_STEP;
    if (strcmp(buf, "SUB-STEP-NOSURF") == 0)
      params->stab = SUB_STEP_NOSURF;
    if (strcmp(buf, "SUB-STEP-TRACER") == 0)
      params->stab = SUB_STEP_TRACER;
  }
  params->thin_merge = 0;
  sprintf(keyword, "MERGE_THIN");
  if (prm_read_char(fp, keyword, buf))
    params->thin_merge = is_true(buf);

  params->means = NONE;

  /* Windows */
  read_window_info(params, fp);

  sprintf(keyword, "PARRAY_INTPL_BOTZ");
  if(prm_read_char(fp, keyword, buf))
  	params->parray_inter_botz = is_true(buf);

  /* Mixed layers */
  sprintf(keyword, "MIX_LAYER");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "NONE") == 0)
      params->mixlayer = NONE;
    if (strcmp(buf, "DENS_MIX") == 0)
      params->mixlayer = DENS_MIX;
    if (strcmp(buf, "TKE_MIX") == 0)
      params->mixlayer = TKE_MIX;
    if (strcmp(buf, "TEMP_MIX") == 0)
      params->mixlayer = TEMP_MIX;
    if (!(params->mixlayer & NONE))
      params->ntrS++;
  }

  /* Vertical mixing scheme */
  sprintf(params->mixsc, "%c", '\0');
  sprintf(params->s_func, "%c", '\0');
  strcpy(params->mixsc, "constant");
  params->min_tke = 0.0;
  params->min_diss = 0.0;
  params->vz0 = 0.0;
  params->kz0 = 0.0;
  prm_set_errfn(hd_quit);
  sprintf(keyword, "MIXING_SCHEME");
  if (prm_read_char(fp, keyword, params->mixsc))
    params->do_closure = 1;    
  if (strcmp(params->mixsc, "constant") == 0 &&
      params->mixlayer == TKE_MIX) {
    params->mixlayer = NONE;
    hd_warn("MIX_LAYER = TKE_MIX not functional with MIXING_SCHEME = constant.\n");
  }
  sprintf(keyword, "KZ0");
  prm_read_double(fp, keyword, &params->kz0);
  sprintf(keyword, "VZ0");
  prm_read_double(fp, keyword, &params->vz0);
  prm_set_errfn(hd_silent_warn);
  sprintf(keyword, "KZ_ALPHA");
  prm_read_double(fp, keyword, &params->kz_alpha);
  sprintf(keyword, "VZ_ALPHA");
  prm_read_double(fp, keyword, &params->vz_alpha);
  sprintf(keyword, "ZS");
  if (!(prm_read_double(fp, keyword, &params->zs))) {
    if (strcmp(params->mixsc, "mellor_yamada_2_0") == 0 ||
        strcmp(params->mixsc, "mellor_yamada_2_5") == 0 ||
        strcmp(params->mixsc, "k-e") == 0 ||
        strcmp(params->mixsc, "k-w") == 0 ||
        strcmp(params->mixsc, "W88") == 0 ||
        strcmp(params->mixsc, "mellor_yamada_2_0_estuarine") == 0)
      hd_quit("%s requires surface length scale, ZS.\n", params->mixsc);
  }
  sprintf(keyword, "STABILITY_FUNC");
  prm_read_char(fp, keyword, params->s_func);
  if (prm_read_char(fp, "SMOOTH_VzKz", buf))
    params->smooth_VzKz = is_true(buf);
  sprintf(keyword, "LMIN");
  prm_read_double(fp, keyword, &params->Lmin);
  sprintf(keyword, "E");
  prm_read_double(fp, keyword, &params->eparam);
  sprintf(keyword, "WAVE_ALPHA");
  prm_read_double(fp, keyword, &params->wave_alpha);
  sprintf(keyword, "WAVE_B1");
  prm_read_double(fp, keyword, &params->wave_b1);
  sprintf(keyword, "WAVE_HEIGHT_FACT");
  prm_read_double(fp, keyword, &params->wave_hf);
  sprintf(keyword, "MIN_TKE");
  prm_read_double(fp, keyword, &params->min_tke);
  sprintf(keyword, "MIN_DISS");
  prm_read_double(fp, keyword, &params->min_diss);
  if (strcmp(params->mixsc, "k-e") == 0)
    params->ntr += 2;
  if (strcmp(params->mixsc, "k-w") == 0)
    params->ntr += 2;
  if (strcmp(params->mixsc, "W88") == 0)
    params->ntr += 2;
  if (strcmp(params->mixsc, "mellor_yamada_2_5") == 0)
    params->ntr += 4;

  /* Flushing tracer */
  sprintf(keyword, "FLUSHING_TR");
  if (prm_read_char(fp, keyword, buf)) {
    params->trflsh = is_true(buf);
    params->ntr += 1;
  }
  sprintf(keyword, "AGE_TR");
  if (prm_read_char(fp, keyword, params->trage)) {
    params->ntr += 1;
  }
  /* Steric height */
  sprintf(keyword, "STERIC_HEIGHT");
  prm_read_double(fp, keyword, &params->lnm);
  if (params->lnm != 0.0)
    params->ntrS++;
  read_profile(params, fp);
  read_debug(params, fp);
  /* Totals diagnostics */
  read_totals(params, fp);
  /* Regions */
  read_region_info(params, fp);
  /* AVHRR SST */
  if (prm_read_char(fp, "AVHRR", params->avhrr_path)) {
    create_avhrr_list(params);
    params->avhrr = 1;
    params->ntrS++;
  }

  /* GHRSST SST */
  if (prm_read_char(fp, "GHRSST", params->ghrsst_path)) {
    prm_read_char(fp, "GHRSST_OPTIONS", params->ghrsst_opt);
    prm_read_char(fp, "GHRSST_INTERP", params->ghrsst_irule);
    create_ghrsst_list(params);
    params->ntrS+=2;
  }
  /* Error norm diagnostic */
  if (prm_read_char(fp, "ERROR_NORM", params->errornorm)) {
    if (prm_read_char(fp, "ERROR_NORM_DT", buf))
      tm_scale_to_secs(buf, &params->enorm_dt);
    else
      params->enorm_dt = 3600.0;
  }
  /* Diagnistic numbers */
  params->ntr += numbers_init(params);
  params->ntr += import_init(params, fp);

  read_decorr(params, fp, 0);
  read_monotone(params, fp, 0);

  /* Auto point source                                               */
  if (prm_read_char(fp, "pss", buf)) {
    params->numbers |= PASS;
    params->ntr += 1;
  }

  /* Means */
  read_means(params, fp, 0);
  /* Alerts */
  sprintf(keyword, "ALERT");
  if (prm_read_char(fp, keyword, params->alert)) {
    if (contains_token(params->alert, "NONE") != NULL)
      memset(params->alert, 0, sizeof(params->alert));
    if (contains_token(params->alert, "ACTIVE") != NULL)
      params->ntrS += 4;
    sprintf(keyword, "ALERT_CODE");
    if (prm_read_char(fp, keyword, params->alertc)) {
      char *fields[MAXSTRLEN * MAXNUMARGS];
      n = parseline(params->alertc, fields, MAXNUMARGS);
      if (n == 11) {
	params->eta_f = atoi(fields[0]);
	params->vel2d_f = atoi(fields[1]);
	params->vel3d_f = atoi(fields[2]);
	params->wvel_f = atoi(fields[3]);
	params->tend_f = atoi(fields[4]);
	params->div2d_f = atoi(fields[5]);
	params->div3d_f = atoi(fields[6]);
	params->cfl_f = atoi(fields[7]);
	params->ts_f = atoi(fields[8]);
	params->shear_f = atoi(fields[9]);
	params->hdiff_f = atoi(fields[10]);
      }
    }
  }
  sprintf(keyword, "ALERT_DT");
  prm_read_char(fp, keyword, params->alert_dt);
  params->vorticity = 0;
  /* Vorticity */
  sprintf(keyword, "VORTICITY");
  if (prm_read_char(fp, keyword, buf)) {
    if (contains_token(buf, "ABSOLUTE") != NULL) {
      params->vorticity |= ABSOLUTE;
      params->ntrS++;
    }
    if (contains_token(buf, "RELATIVE") != NULL) {
      params->vorticity |= RELATIVE;
      params->ntrS++;
    }
    if (contains_token(buf, "POTENTIAL") != NULL) {
      params->vorticity |= POTENTIAL;
      params->ntrS++;
    }
    if (contains_token(buf, "TENDENCY") != NULL) {
      params->vorticity |= TENDENCY;
      params->ntrS += 7;
      if (!(params->vorticity & RELATIVE)) {
        params->vorticity |= RELATIVE;
        params->ntrS++;
      }
    }
  } else
    params->vorticity = NONE;

  sprintf(keyword, "CFL");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "NONE") == 0)
      params->cfl = NONE;
    if (strcmp(buf, "PASSIVE") == 0)
      params->cfl = PASSIVE;
    if (strcmp(buf, "PASSIVE|WVEL") == 0)
      params->cfl = PASSIVE | WVEL;
    if (strcmp(buf, "ACTIVE3D") == 0)
      params->cfl = ACTIVE3D;
    if (strcmp(buf, "ACTIVE3D|WVEL") == 0)
      params->cfl = ACTIVE3D | WVEL;
    if (strcmp(buf, "ACTIVE") == 0)
      params->cfl = ACTIVE;
    if (strcmp(buf, "ACTIVE|WVEL") == 0)
      params->cfl = ACTIVE | WVEL;
    if (params->cfl & (ACTIVE | ACTIVE3D))
      prm_set_errfn(hd_quit);
    sprintf(keyword, "CFL_DT");
    prm_read_char(fp, keyword, params->cfl_dt);
    prm_set_errfn(hd_silent_warn);
  }
  if (!(params->cfl & NONE))
    params->ntrS += 6;

  /* Waves */
#if defined(HAVE_WAVE_MODULE)
  params->do_wave = NONE;
  sprintf(keyword, "DO_WAVES");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "NONE") == 0)
      params->do_wave = NONE;
    if (strcmp(buf, "FILE") == 0)
      params->do_wave = W_FILE;
    if (strcmp(buf, "COMP") == 0)
      params->do_wave = W_COMP;
    if (strcmp(buf, "SWAN") == 0)
      params->do_wave = (W_SWAN|W_SWANM);
    if (strcmp(buf, "SWAN_W") == 0)
      params->do_wave = (W_SWAN|W_SWANW);
    prm_get_time_in_secs(fp, "WAVES_DT", &params->wavedt);
  }
  read_waves(params, fp, 0);
#endif

  /* Library error handling */
  sprintf(keyword, "LIB_ERROR_FCN");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "QUIT") == 0)
      params->gint_errfcn = LFATAL;
    if (strcmp(buf, "WARN") == 0)
      params->gint_errfcn = LWARN;
  }

  /* Tracer perceltiles */
  sprintf(keyword, "CALC_PERCS");
  if (prm_read_char(fp, keyword, params->trperc)) {
    if (strcmp(params->trperc, "NONE") != 0) {
      params->ntr += 1;
      sprintf(keyword, "PERC_REGION");
      prm_read_char(fp, keyword, params->trpercr);
    }
  } else 
    sprintf(params->trperc, "NONE");

  /* Particle tracking */
  if(prm_read_char(fp, "PT_InputFile", params->ptinname)) {
    params->do_pt = 1;
    params->ntr++;
    /* Auto particle source                                          */
    if (prm_read_char(fp, "particles", params->particles))
      params->do_pt = 2;
    params->do_lag = params->do_pt;
  }

  /* Fixed constants */
  params->g = 9.81;
  params->spec_heat = 3990;
  air_dens = params->air_dens = 1.225;
  params->ambpress = 100800;
  params->u1kh = 0.0;
  params->u2kh = 0.0;
  read_hdiff(params, fp, 0);
  /*
  if (params->smagorinsky > 0.0)
    params->ntr++;
  */
  sprintf(keyword, "ETAMAX");
  prm_read_double(fp, keyword, &params->etamax);
  sprintf(keyword, "VELMAX");
  prm_read_double(fp, keyword, &params->velmax);
  sprintf(keyword, "VELMAX_2D");
  if(!prm_read_double(fp, keyword, &params->velmax2d))
    params->velmax2d = params->velmax;

  /* MOM grid conversion */
  if (prm_read_char(fp, "MOM_CONVERT", params->momfile))
    params->domom |= (CDF|RDGRID);

  /* ROMS grid conversion */
  if (prm_read_char(fp, "ROMS_CONVERT", params->romsfile))
    params->doroms |= (CDF|RDGRID);

  /* Get the Z to Sigma levels scaling factor for ROMS */
  params->roms_z2s = 1.0;
  prm_read_double(fp, "ROMS_Z2S_FACTOR", &params->roms_z2s);

  /* Wind */
  sprintf(params->wind, "%c", '\0');
  params->wind_dt = params->storm_dt = params->wind_scale = 0.0;
  prm_set_errfn(hd_silent_warn);
  params->wind_type = SPEED;
  if (prm_read_char(fp, "WIND_TS", params->wind)) {
    prm_set_errfn(hd_quit);
    prm_get_time_in_secs(fp, "WIND_INPUT_DT", &params->wind_dt);
    prm_read_double(fp, "WIND_SPEED_SCALE", &params->wind_scale);
    prm_set_errfn(hd_silent_warn);
    if (prm_read_char(fp, "WIND_STRESS_FCTN", buf)) {
      prm_set_errfn(hd_quit);
      if(strcmp(buf, "L&P") == 0) {
	params->stress_fn = LARGEPOND;
	prm_read_double(fp, "WIND_STRESS_REFH", &params->dlv0);
      } else if(strcmp(buf, "B") == 0) {
	params->stress_fn = BUNKER;
	prm_read_double(fp, "WIND_STRESS_REFH", &params->dlv0);
      } else if(strcmp(buf, "K/W") == 0) {
	params->stress_fn = KITIAG;
	prm_read_double(fp, "WIND_STRESS_REFH", &params->dlv0);
      } else {
	prm_set_errfn(hd_quit);
	params->stress_fn = ORIGINAL;
	prm_read_double(fp, "DRAG_LAW_V0", &params->dlv0);
	prm_read_double(fp, "DRAG_LAW_V1", &params->dlv1);
	prm_read_double(fp, "DRAG_LAW_CD0", &params->dlc0);
	prm_read_double(fp, "DRAG_LAW_CD1", &params->dlc1);
      }
    } else {
      prm_set_errfn(hd_quit);
      params->stress_fn = ORIGINAL;
      prm_read_double(fp, "DRAG_LAW_V0", &params->dlv0);
      prm_read_double(fp, "DRAG_LAW_V1", &params->dlv1);
      prm_read_double(fp, "DRAG_LAW_CD0", &params->dlc0);
      prm_read_double(fp, "DRAG_LAW_CD1", &params->dlc1);
    }
    if (prm_read_char(fp, "WIND_TYPE", buf)) {
      if(strcmp(buf, "STRESS") == 0)
	params->wind_type = STRESS;
      else
	params->wind_type = SPEED;
    }
  }
  if (prm_read_int(fp, "NSTORM", &params->nstorm)) {
    prm_set_errfn(hd_quit);
    prm_get_time_in_secs(fp, "STORM_INPUT_DT", &params->storm_dt);
    if (!params->wind_scale) {
      prm_read_double(fp, "WIND_SPEED_SCALE", &params->wind_scale);
      prm_read_double(fp, "DRAG_LAW_V0", &params->dlv0);
      prm_read_double(fp, "DRAG_LAW_V1", &params->dlv1);
      prm_read_double(fp, "DRAG_LAW_CD0", &params->dlc0);
      prm_read_double(fp, "DRAG_LAW_CD1", &params->dlc1);
    }
  }

  /* Atmospherics */
  prm_set_errfn(hd_silent_warn);
  sprintf(params->patm, "%c", '\0');
  sprintf(keyword, "PRESSURE");
  prm_read_char(fp, keyword, params->patm);
  sprintf(keyword, "PRESSURE_INPUT_DT");
  prm_get_time_in_secs(fp, keyword, &params->patm_dt);
  sprintf(params->precip, "%c", '\0');
  sprintf(keyword, "PRECIPITATION");
  prm_read_char(fp, keyword, params->precip);
  sprintf(keyword, "PRECIPITATION_INPUT_DT");
  prm_get_time_in_secs(fp, keyword, &params->precip_dt);
  sprintf(params->evap, "%c", '\0');
  sprintf(keyword, "EVAPORATION");
  prm_read_char(fp, keyword, params->evap);
  sprintf(keyword, "EVAPORATION_INPUT_DT");
  prm_get_time_in_secs(fp, keyword, &params->evap_dt);
  sprintf(params->airtemp, "%c", '\0');
  sprintf(keyword, "AIRTEMP");
  prm_read_char(fp, keyword, params->airtemp);
  sprintf(keyword, "AIRTEMP_INPUT_DT");
  prm_get_time_in_secs(fp, keyword, &params->airtemp_dt);
  sprintf(params->rh, "%c", '\0');
  sprintf(keyword, "HUMIDITY");
  prm_read_char(fp, keyword, params->rh);
  sprintf(keyword, "HUMIDITY_INPUT_DT");
  prm_get_time_in_secs(fp, keyword, &params->rh_dt);
  sprintf(params->cloud, "%c", '\0');
  sprintf(keyword, "CLOUD");
  prm_read_char(fp, keyword, params->cloud);
  sprintf(keyword, "CLOUD_INPUT_DT");
  prm_get_time_in_secs(fp, keyword, &params->cloud_dt);
  sprintf(params->swr, "%c", '\0');
  params->albedo = -1;
  sprintf(keyword, "RADIATION");
  prm_read_char(fp, keyword, params->swr);
  sprintf(keyword, "RADIATION_INPUT_DT");
  prm_get_time_in_secs(fp, keyword, &params->swr_dt);
  prm_read_double(fp, "ALBEDO", &params->albedo);
  sprintf(params->light, "%c", '\0');
  params->albedo_l = -1;
  sprintf(keyword, "LIGHT");
  prm_read_char(fp, keyword, params->light);
  sprintf(keyword, "LIGHT_INPUT_DT");
  prm_get_time_in_secs(fp, keyword, &params->light_dt);
  prm_read_double(fp, "ALBEDO_LIGHT", &params->albedo_l);
  sprintf(params->wetb, "%c", '\0');
  sprintf(keyword, "WET_BULB");
  prm_read_char(fp, keyword, params->wetb);
  sprintf(keyword, "WET_BULB_INPUT_DT");
  prm_get_time_in_secs(fp, keyword, &params->wetb_dt);
  sprintf(keyword, "DEW_POINT");
  prm_read_char(fp, keyword, params->wetb);
  sprintf(keyword, "DEW_POINT_INPUT_DT");
  prm_get_time_in_secs(fp, keyword, &params->wetb_dt);

  /* Heatflux */
  prm_set_errfn(hd_silent_warn);
  params->heatflux = NONE;
  sprintf(params->hftemp, "%c", '\0');
  sprintf(keyword, "HEATFLUX");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "NONE") == 0)
      params->heatflux = NONE;
    if (strcmp(buf, "ORIGINAL") == 0) {
      hd_warn
        ("ORIGINAL HEATFLUX formulation no longer supported : using ADVANCED formulation\n");
      params->heatflux = ADVANCED;
    }
    if ((strcmp(buf, "NET_HEAT") == 0) || (strcmp(buf, "COMP_HEAT") == 0)) {
    if (strcmp(buf, "COMP_HEAT") == 0)
      params->heatflux = COMP_HEAT;
    if (strcmp(buf, "NET_HEAT") == 0)
      params->heatflux = NET_HEAT;

      /* Read the swr parameters                                     */
    read_swr(params, fp, 0);
    }
    if (strcmp(buf, "COMP_HEAT_NONE") == 0) {
      params->heatflux = COMP_HEAT_NONE;
      sprintf(params->hf, "%c", '\0');
      prm_set_errfn(hd_quit);
      sprintf(keyword, "HEATFLUX_FILE");
      prm_read_char(fp, keyword, params->hf);
      sprintf(keyword, "HEATFLUX_DT");
      prm_get_time_in_secs(fp, keyword, &params->hf_dt);
      params->ntrS += 4+2; // 4 x hf comps + precip and evap
    }

    if (strcmp(buf, "ADVANCED") == 0 || params->heatflux & ADVANCED) {
      params->heatflux = ADVANCED;

      /* Read the swr parameters                                     */
      read_swr(params, fp, 0);

      sprintf(keyword, "HEATFLUX_REFH");
      if (!(prm_read_double(fp, keyword, &params->zref)))
        params->zref = 10.0;
      /* Note : the input codes 0 - 4 are included for backwards */
      /* compatibility.                                          */
      if (prm_read_char(fp, "BULK_SCHEME", buf)) {
	if((strcmp(buf, "L&P") == 0) || (strcmp(buf, "1") == 0))
	  params->bulkf = LARGEPOND;
	else if((strcmp(buf, "B") == 0) || (strcmp(buf, "4") == 0))
	  params->bulkf = BUNKER;
	else if((strcmp(buf, "K/W") == 0) || (strcmp(buf, "3") == 0))
	  params->bulkf = KITIAG;
	else if((strcmp(buf, "Ko") == 0) || (strcmp(buf, "0") == 0))
	  params->bulkf = KONDO;
	else if((strcmp(buf, "M") == 0) || (strcmp(buf, "2") == 0))
	  params->bulkf = MASAG;
      } else
        params->bulkf = KONDO;
      params->hfadf = 0;
      sprintf(keyword, "HEATFLUX_ADVECT");
      if (prm_read_char(fp, keyword, buf1))
        params->hfadf = is_true(buf1);
      sprintf(keyword, "HEATFLUX_TEMP");
      prm_read_char(fp, keyword, params->hftemp);
      sprintf(keyword, "HEATFLUX_TEMP_DT");
      prm_get_time_in_secs(fp, keyword, &params->hftemp_dt);
      sprintf(keyword, "HEATFLUX_TC");
      prm_get_time_in_secs(fp, keyword, &params->hftc);
      params->ntrS += 5;
    }
    if (strcmp(buf, "SURF_RELAX") == 0) {
      params->heatflux = SURF_RELAX;
      sprintf(keyword, "HEATFLUX_TEMP");
      if (!(prm_read_char(fp, keyword, params->hftemp)))
	hd_warn
	  ("params_read() : SURF_RELAX heatflux requires HEATFLUX_TEMP file.\n");
      sprintf(keyword, "HEATFLUX_TEMP_DT");
      if (!(prm_get_time_in_secs(fp, keyword, &params->hftemp_dt)))
	hd_warn
	  ("params_read() : SURF_RELAX heatflux requires HEATFLUX_TEMP_DT constant.\n");
      sprintf(keyword, "HEATFLUX_TC");
      if (!(prm_get_time_in_secs(fp, keyword, &params->hftc)))
	hd_warn
	  ("params_read() : SURF_RELAX heatflux requires HEATFLUX_TC constant.\n");
    }
    if (strcmp(buf, "AVHRR") == 0) {
      params->heatflux = AVHRR;
      if (!strlen(params->avhrr_path)) {
	hd_warn
	  ("params_read() : AVHRR heatflux requires an AVHRR diagnostic path name.\n");
	params->heatflux = NONE;
      }
      sprintf(keyword, "HEATFLUX_TC");
      if (!(prm_get_time_in_secs(fp, keyword, &params->hftc)))
	hd_warn
	  ("params_read() : AVHRR heatflux requires HEATFLUX_TC constant.\n");
    }
    if (strcmp(buf, "INVERSE") == 0) {
      params->heatflux = INVERSE;
      sprintf(keyword, "HEATFLUX_TEMP");
      prm_read_char(fp, keyword, params->hftemp);
      sprintf(keyword, "HEATFLUX_TEMP_DT");
      prm_get_time_in_secs(fp, keyword, &params->hftemp_dt);
      sprintf(keyword, "HEATFLUX_TC");
      prm_get_time_in_secs(fp, keyword, &params->hftc);
      if(params->mixlayer & NONE) {
	hd_quit_and_dump
	  ("params_read() : INVERSE heatflux requires MIX_LAYER to be set.\n");
      }
      params->ntrS += 1;
    }
    if (strcmp(buf, "NET_HEAT") == 0) {
      params->heatflux = NET_HEAT;
      sprintf(params->hf, "%c", '\0');
      prm_set_errfn(hd_quit);
      sprintf(keyword, "HEATFLUX_FILE");
      prm_read_char(fp, keyword, params->hf);
      sprintf(keyword, "HEATFLUX_DT");
      prm_get_time_in_secs(fp, keyword, &params->hf_dt);
      prm_set_errfn(hd_silent_warn);
      sprintf(keyword, "HEATFLUX_TC");
      prm_get_time_in_secs(fp, keyword, &params->hftc);
      params->ntrS += 2;
    }
    if (strcmp(buf, "COMP_HEAT") == 0) {
      params->heatflux = COMP_HEAT;
      sprintf(params->hf, "%c", '\0');
      prm_set_errfn(hd_quit);
      sprintf(keyword, "HEATFLUX_FILE");
      prm_read_char(fp, keyword, params->hf);
      sprintf(keyword, "HEATFLUX_DT");
      prm_get_time_in_secs(fp, keyword, &params->hf_dt);
      prm_set_errfn(hd_silent_warn);
      sprintf(keyword, "HEATFLUX_TC");
      prm_get_time_in_secs(fp, keyword, &params->hftc);
      params->ntrS += 5;
    }
    sprintf(keyword, "HEATFLUX_RAMP");
    if (prm_read_char(fp, keyword, buf))
      tm_scale_to_secs(buf, &params->hf_ramp);
    else
      params->hf_ramp = params->t;
  }

  /* Saltflux */
  prm_set_errfn(hd_silent_warn);
  params->saltflux = NONE;
  sprintf(keyword, "SALTFLUX");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "NONE") == 0)
      params->saltflux = NONE;
    if (strcmp(buf, "ORIGINAL") == 0) {
      params->saltflux = ORIGINAL;
      params->ntrS += 1;
    }
    if (strcmp(buf, "ADVANCED") == 0) {
      params->saltflux = ADVANCED;
      params->ntrS += 1;
    }
    if (strcmp(buf, "BULK") == 0) {
      params->saltflux = BULK;
      params->ntrS += 1;
    }
  }

  prm_set_errfn(hd_silent_warn);

  /* Sediment geometry */
  params->sednz = 0;
  /* Note : sediment layers are read in assuming the first layer is */
  /* closest to the water column. This is reversed when copying to */
  /* the master so that the sediment origin is the deepest layer.  */
  if (prm_read_darray(fp, "LAYERFACES_SED", &params->gridz_sed,
                      &params->sednz)) {
    if (--params->sednz < 1)
      hd_quit("Number of sediment layers must be 2 or more\n");
  } else {
    double *dz_sed = NULL;
    if (prm_read_int(fp, "NSEDLAYERS", &params->sednz)) {
      dz_sed = d_alloc_1d(params->sednz);
      if (prm_read_darray(fp, "NSEDLAYERS", &dz_sed, &params->sednz)) {
        if (params->sednz) {
          params->gridz_sed = d_alloc_1d(params->sednz + 1);

          params->gridz_sed[params->sednz] = 0.0;
          for (m = params->sednz - 1; m >= 0; m--) {
            params->gridz_sed[m] = params->gridz_sed[m + 1] -
              dz_sed[params->sednz - m - 1];
          }
        }
      }
      if(dz_sed)
        d_free_1d(dz_sed);
    }
  }

  sprintf(keyword, "RENDER_ECOSED_NAME");
  if (prm_read_char(fp, keyword, params->rendername)) {
    sprintf(keyword, "RENDER_ECOSED_DESC");
    prm_read_char(fp, keyword, params->renderdesc);
  }
  sprintf(keyword, "RENDER_ECOSED_PATH");
  prm_read_char(fp, keyword, params->renderpath);
  sprintf(keyword, "RENDER_ECOSED_REMOVE");
  prm_read_char(fp, keyword, params->renderrem);
  sprintf(keyword, "ECOSED_CONFIG");
  prm_read_char(fp, keyword, params->ecosedconfig);

#if defined(HAVE_SEDIMENT_MODULE)
  read_sediments(params, fp, &params->ntr);
#endif

#if defined(HAVE_ECOLOGY_MODULE)
  read_ecology(params, fp, &params->ntr);
#endif

  /*
    sprintf(keyword, "LIB_ERROR_FCN");
    if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "QUIT") == 0)
    params->gint_errfcn = LFATAL;
    if (strcmp(buf, "WARN") == 0)
    params->gint_errfcn = LWARN;
    }
  */

  read_trfilter(params, fp);

  /* 3D Tracer constants and variables.  */
  /* Tracers relating to certain diagnostics etc. are automatically */
  /* generated if the relevant flag is set. If these tracers are */
  /* manually defined in the input parameter file, then the manual */
  /* definition over-rides the auto generation, and the defined */
  /* tracer_info is used.  */
  /* Always assume salt and temp are auto tracers */
  params->ntr += 2;
  params->atr = params->ntr;
  tracer_setup(params, fp);
  create_tracer_3d(params);
  for (n = 0; n < params->atr; n++) {
    params->trinfo_3d[n].m = -1;
    params->trinfo_3d[n].n = n;
  }

  /* Get the tracers which undergo diffusion (horizontal and */
  /* vertical)  and store in tdif_h[] and tdif_v[]. */
  params->ntdif_h = params->ntdif_v = 0;
  for (n = 0; n < params->ntr; n++) {
    tracer_info_t *tracer = &params->trinfo_3d[n];
    if (!tracer->diagn && tracer->type&WATER && tracer->diffuse) {
      params->ntdif_h++;
#if defined(HAVE_SEDIMENT_MODULE)
      if (params->do_sed) {
        if (strcmp(tracer->name, "salt") == 0 ||
            strcmp(tracer->name, "temp") == 0)
          params->ntdif_v++;
      } else
        params->ntdif_v++;
#else
      params->ntdif_v++;
#endif
    }
  }
  if (params->ntdif_h)
    params->tdif_h = i_alloc_1d(params->ntdif_h);
  m = 0;
  for (n = 0; n < params->ntr; n++) {
    tracer_info_t *tracer = &params->trinfo_3d[n];
    if (!tracer->diagn && tracer->type&WATER && tracer->diffuse) {
      params->tdif_h[m] = n;
      m++;
    }
  }
  if (params->ntdif_v)
    params->tdif_v = i_alloc_1d(params->ntdif_v);
  m = 0;
  for (n = 0; n < params->ntr; n++) {
    tracer_info_t *tracer = &params->trinfo_3d[n];
    if (!tracer->diagn && tracer->type&WATER && tracer->diffuse) {
#if defined(HAVE_SEDIMENT_MODULE)
      if (params->do_sed) {
        if (strcmp(tracer->name, "salt") == 0 ||
            strcmp(tracer->name, "temp") == 0) {
          params->tdif_v[m] = n;
          m++;
        }
      } else {
        params->tdif_v[m] = n;
        m++;
      }
#else
      params->tdif_v[m] = n;
      m++;
#endif
    }
  }

  /* Sediment tracer constants and variables */
  prm_set_errfn(hd_silent_warn);
  tracer_read(fp, NULL, SEDIM, hd_quit, hd_warn, hd_silent_warn,
        &params->nsed, &params->trinfo_sed);
  if (params->nsed && !params->sednz)
    hd_quit("Sediment tracers requires NSEDLAYERS > 0\n");

  /* 2D Tracer constants and variables */
  prm_set_errfn(hd_silent_warn);
  params->atrS = params->ntrS;
  tracer_read(fp, NULL, INTER, hd_quit, hd_warn, hd_silent_warn,
              &params->ntrS, &params->trinfo_2d);

  /* Timeseries file caching */
  prm_set_errfn(hd_silent_warn);
  sprintf(keyword, "CACHE_TSFILES");
  if (prm_read_char(fp, keyword, buf) > 0)
    params->tsfile_caching = is_true(buf);
  else
    params->tsfile_caching = 1;
  emstag(LDEBUG,"hd:readparam:params_read","Setting ts_file_cachinf to: %s",(params->tsfile_caching?"true":"false"));

  /* Explicit mappings */
  read_explicit_maps(params, fp);

  /* Open boudaries */
  prm_set_errfn(hd_quit);
  get_bdry_params(params, fp);
  /* Allocate memory for the open boundary data structure */
  params->open =
    (open_bdrys_t **)malloc(sizeof(open_bdrys_t *) * params->nobc);

  for (n = 0; n < params->nobc; n++) {

    /* Allocate memory for the boundary structures */
    params->open[n] = OBC_alloc();

    /* Read the open boudary conditions */
    params->open[n]->ntr = params->ntr;
    params->open[n]->atr = params->atr;
    get_OBC_conds(params, params->open[n], fp, n, params->trinfo_3d);
    get_obc_list(params->open[n], fp, n, "BOUNDARY");
    prm_set_errfn(hd_quit);
  }

  /* Read the csr tide model paths if required */
  for (n = 0; n < params->nobc; n++) {
    if(params->open[n]->bcond_ele & TIDALH) {
      prm_set_errfn(hd_quit);
      sprintf(keyword, "TIDE_CSR_ORTHOWEIGHTS");
      prm_read_char(fp, keyword, params->orthoweights);
      sprintf(keyword, "TIDE_CSR_CON_DIR");
      prm_read_char(fp, keyword, params->nodal_dir);
      prm_set_errfn(hd_silent_warn);
      break;
    }
    if(params->open[n]->bcond_ele & TIDALC) {
      prm_set_errfn(hd_quit);
      sprintf(keyword, "TIDE_CSR_CON_DIR");
      prm_read_char(fp, keyword, params->nodal_dir);
      sprintf(keyword, "TIDE_CONSTITUENTS");
      prm_read_char(fp, keyword, params->tide_con_file);
      break;
    }
  }

  /* Determine if velocity data is to be read in */
  for (n = 0; n < params->nobc; n++) {
    if(params->open[n]->bcond_nor & (FILEIN|CUSTOM))
      params->tmode |= DO_OBC;
    if(params->open[n]->bcond_tan & (FILEIN|CUSTOM))
      params->tmode |= DO_OBC;
  }

  /* Process exclusion */
  read_exclude_points(params, fp);

  /* Output files to write to setup.txt */
  get_output_path(params, fp);
  sprintf(keyword, "OutputFiles");
  if (prm_read_char(fp, keyword, buf)) {
    if (fopen(buf, "r") == NULL) {
      params->ndf = atoi(buf);
      params->d_name = (char **)malloc(params->ndf * sizeof(char *));
      for (n = 0; n < params->ndf; n++) {
	params->d_name[n] = (char *)malloc(sizeof(char)*MAXSTRLEN);
	sprintf(keyword, "file%1.1d.name",n);
	prm_read_char(fp, keyword, params->d_name[n]);
	sprintf(keyword, "file%1.1d.filetype",n);
	prm_read_char(fp, keyword, buf);
	if (strcmp(buf, "mom") == 0) params->domom |= RDGRID;
	if (strcmp(buf, "roms") == 0) params->doroms |= RDGRID;
      }
    }
  }
  sprintf(keyword, "TRANS_OUTPUT");
  if (prm_read_char(fp, keyword, buf))
    params->trout = is_true(buf);
  if (!prm_read_char(fp, "OutputTransport", params->trkey)) {
    if (params->trout) {
      strcpy(buf, params->prmname);
      stripend(buf);
      strcpy(params->trkey, buf);
    }
  }

  /* Whether we read do threaded I/O on transport file resets */
  params->thIO = 0;
  if (prm_read_char(fp, "SCHED_MODE", buf)) {
    if (strcasecmp(buf, "pthreads") == 0) {
#ifdef HDF5_THREADSAFE
      params->thIO = 1;
#else
      hd_quit("SCHED_MODE error: This executable is not built with the correct HDF5 threadsafe library version via NETCDF4");
#endif
    }
  }

  if (DEBUG("init_m"))
    dlog("init_m", "Input parameters read OK\n");
  return (params);
}

/* END params_read_t()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads in basic transport mode parameters                          */
/*-------------------------------------------------------------------*/
void read_tmode_params(FILE *fp, parameters_t *params)
{
  char buf[MAXSTRLEN];
  char keyword[MAXSTRLEN];

  strcpy(params->sourcefile, params->idumpname);
  prm_read_char(fp, "TRANS_DATA", params->trans_data);
  prm_read_char(fp, "TRANS_VARS", params->trvars);
  check_wave_vars(params);
  params->tmode = SP_EXACT;
  sprintf(keyword, "TRANS_MODE");
  if (prm_read_char(fp, keyword, buf)) {
    if (contains_token(buf, "XYZ_INTERP") != NULL)
      params->tmode = XYZ_TINT;
    if (contains_token(buf, "SP_INTERP") != NULL)
      params->tmode = SP_TINT;
    if (contains_token(buf, "SP_EXACT") != NULL)
      params->tmode = SP_EXACT;
    if (contains_token(buf, "SP_EXACT|INEXACT") != NULL)
      params->tmode = SP_EXACT|INEXACT;
    if (contains_token(buf, "SP_INTERP|INEXACT") != NULL)
      params->tmode = SP_TINT|INEXACT;
    if (contains_token(buf, "STREAMLINE") != NULL)
      params->tmode = SP_ORIGIN|SP_EXACT;
    if (contains_token(buf, "NONE") != NULL)
      params->tmode = NONE;
    if (contains_token(buf, "SP_DUMP") != NULL)
      params->tmode = SP_DUMP;
    if (contains_token(buf, "SP_CHECK") != NULL)
      params->tmode = (SP_EXACT|SP_CHECK);
    if (contains_token(buf, "TR_CHECK") != NULL)
      params->tmode = (SP_EXACT|TR_CHECK);
    if (contains_token(buf, "GLOBAL") != NULL)
      params->tmode = (GLOBAL|XYZ_TINT);
    if (contains_token(buf, "SIMPLE") != NULL)
      params->tmode = (SP_SIMPLE|XYZ_TINT);
    if (contains_token(buf, "SIMPLEU") != NULL)
      params->tmode = (SP_SIMPLEU|XYZ_TINT);
    if (contains_token(buf, "SP_FFSL") != NULL) {
      params->tmode = (SP_FFSL|SP_EXACT);
      if (params->runmode & TRANS) params->tmode |= SP_U1VM;
      params->ntrS += 1;    
    }
    if (contains_token(buf, "SP_FFSLS") != NULL) {
      params->tmode = (SP_FFSL|SP_EXACT|SP_STRUCT);
      if (params->runmode & TRANS) params->tmode |= SP_U1VM;
      params->ntrS += 1;    
    }
    if (contains_token(buf, "SP_FFSLU") != NULL) {
      params->tmode = (SP_FFSL|SP_EXACT|SP_UGRID);
      if (params->runmode & TRANS) params->tmode |= SP_U1VM;
      params->ntrS += 1;
    }
  }

  if (params->runmode & TRANS)
    params->fillf = MONOTONIC;
  if (params->tmode & NONE) params->fillf = NONE;
  sprintf(keyword, "FILL_METHOD");
  if (prm_read_char(fp, keyword, buf)) {
    if (contains_token(buf, "NONE") != NULL) {
      params->fillf = NONE;
    } else {
      params->fillf = 0;
      if (contains_token(buf, "GLOBAL") != NULL)
	params->fillf |= GLOBAL;
      if (contains_token(buf, "MONOTONIC") != NULL)
	params->fillf |= MONOTONIC;
      if (contains_token(buf, "CLIP") != NULL)
	params->fillf |= CLIP;
      if (contains_token(buf, "OBC_ADJUST") != NULL)
	params->fillf |= OBC_ADJUST;
      if (contains_token(buf, "WEIGHTED") != NULL)
	params->fillf |= WEIGHTED;
      if (contains_token(buf, "DIAGNOSE") != NULL)
	params->fillf |= DIAGNOSE;
      if (contains_token(buf, "DIAGNOSE_BGC") != NULL)
	params->fillf |= DIAGNOSE_BGC;
      if (contains_token(buf, "MONO+GLOB") != NULL)
	params->fillf |= (MONGLOB|MONOTONIC);
    }
  }
  if (params->fillf & (WEIGHTED|MONOTONIC))
    params->ntrS += 1;
  
  params->pssinput = WIMPLICIT;
  sprintf(keyword, "PSS_INPUT");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "EXPLICIT") == 0)
      params->pssinput = EXPLICIT;
    if (strcmp(buf, "IMPLICIT") == 0)
      params->pssinput = WIMPLICIT;
  }
  params->conserve = NONE;
  sprintf(keyword, "CONSERVATION");
  if (prm_read_char(fp, keyword, buf)) {
    if (contains_token(buf, "NONE") != NULL) {
      params->conserve = NONE;
    } else {
      params->conserve = 0;
      if (contains_token(buf, "RE_INITIALIZE") != NULL)
	params->conserve |= REINIT;
      else if (contains_token(buf, "NO_INITIALIZE") != NULL)
	params->conserve |= NOINIT;
      else
	params->conserve |= REINIT;
      if (contains_token(buf, "W") != NULL)
	params->conserve |= CONS_W;
      if (contains_token(buf, "WSTAB") != NULL)
	params->conserve |= (CONS_W|CONS_WS);
      if (contains_token(buf, "ETA") != NULL)
	params->conserve |= CONS_ETA;
      if (contains_token(buf, "W") == NULL && contains_token(buf, "ETA") == NULL)
	params->conserve |= (CONS_W|CONS_ETA);
      if (contains_token(buf, "NO_PSS_FLOW") != NULL)
	params->conserve |= CONS_NOF;
      if (contains_token(buf, "NOGRAD") != NULL)
	params->conserve |= CONS_NGR;
      if (contains_token(buf, "MERGED") != NULL)
	params->conserve |= CONS_MRG;
    }
    if (params->tratio < 1.0) params->conserve |= CONS_SUB;
  }

  sprintf(keyword, "NO_LEAP_YEARS");
  if (prm_read_char(fp, keyword, buf))
    params->lyear = is_true(buf);

#ifdef HAVE_OMP
  sprintf(keyword, "TRANS_OMP_NUM_THREADS");
  params->trans_num_omp = 1; /* default to 1 */
  prm_read_int(fp, keyword, &params->trans_num_omp);
#endif

}

/* END read_tmode_params()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to check if wave variables are in the TRANS_VARS, and if  */
/* so set flags for autotracers and sediment code.                   */
/*-------------------------------------------------------------------*/
void check_wave_vars(parameters_t *params)
{
  char trvars[MAXSTRLEN];

  strcpy(trvars, params->trvars);
  if (strlen(trvars)) {
    int nfiles, tt, wf = 0;
    char files[MAXNUMTSFILES][MAXSTRLEN], trfile[MAXSTRLEN];
    nfiles = parseline(trvars, (char **)files, MAXNUMTSFILES);
    if(nfiles) {
      for (tt = 0; tt < nfiles; ++tt) {
	strcpy(trfile, ((char **)files)[tt]);
	if (strcmp(trfile, "wave_amp") == 0) wf++;
	if (strcmp(trfile, "wave_dir") == 0) wf++;
	if (strcmp(trfile, "wave_period") == 0) wf++;
	if (strcmp(trfile, "wave_ub") == 0) wf++;
      }
      if (wf == 4) {
	params->do_wave = W_FILE;
	params->webf_dt = 0.0;
	params->waves |= ORBITAL;
	params->ntrS += 5;
	hd_warn("Wave variables read from TRANS_DATA\n");
      }
    }
  }
}

/* END check_wave_vars()                                             */
/*-------------------------------------------------------------------*/
    
