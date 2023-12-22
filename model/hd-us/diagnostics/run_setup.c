/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/diagnostics/run_setup.c
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
 *  $Id: run_setup.c 7458 2023-12-13 03:50:02Z her127 $
 *
 */

#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "hd.h"
#include "eqn_parser.h"

void wbdrycustom(FILE * fp, parameters_t *params, int bnum,
                 bdry_details_t *data);
void write_grid_specs(FILE *op, parameters_t *params);
/* defined in readparam.c */
char *btname(int m);
char *heatfluxname(int m);

/*-------------------------------------------------------------------*/
/* Routine to print a summary of a simulation setup                  */
/*-------------------------------------------------------------------*/
void write_run_setup(hd_data_t *hd_data)
{
  master_t *master = hd_data->master;
  geometry_t *geom = hd_data->geom;
  parameters_t *params = hd_data->params;
  int n, nn, ntr;
  int c, cc;
  long t;
  double d1, d2, d3, d4;
  double maxvel = 1.0;
  double int_wave_speed = 2.0;
  char bname[MAXSTRLEN], buf[MAXSTRLEN];
  struct timeval tm1;
  FILE *fp;

  if (!setup_log)
    return;

  if(!(fp = fopen(setup_logfile, "w")))
     hd_warn("Can't open setup file %s\n", setup_logfile);
  fprintf(fp, "COMPAS Simulation Summary\n");
  fprintf(fp, "EMS Version : %s\n", version);
  fprintf(fp, "Input file = %s\n", params->prmname);
  fprintf(fp, "Executable file = %s\n",executable );
  getcwd(bname, MAXSTRLEN);
  fprintf(fp, "Working directory = %s\n",bname);
  fprintf(fp, "%s\n", params->parameterheader);
  if (params->runno) fprintf(fp, "Identifier # %s\n", params->runnoc);
  if (strlen(params->runcode)) {
    char *tok;
    fprintf(fp, "Run code %s\n", params->runcode);
    strcpy(bname, params->runcode);
    tok = strtok(bname, "|");
    if (tok != NULL) fprintf(fp, "  Grid name : %s\n", tok);
    strcpy(bname, params->runcode);
    get_idcodec(bname, "G", buf);
    fprintf(fp, "  Grid ID : %s\n", buf);
    get_idcodec(bname, "H", buf);
    fprintf(fp, "  Hydrodynamic model ID : %s\n", buf);
    get_idcodec(bname, "S", buf);
    fprintf(fp, "  Sediment model ID : %s\n", buf);
    get_idcodec(bname, "B", buf);
    fprintf(fp, "  Biogeochemical model ID : %s\n", buf);
  }
  if (strlen(params->rev)) fprintf(fp, "Parameter file revision %s\n", params->rev);
  if (strlen(params->trl)) fprintf(fp, "Technology_Readiness_Level %s\n", params->trl);
  if (strcmp(params->trl, "TR1") == 0)
    fprintf(fp, "  Basic principles observed and reported.\n");
  if (strcmp(params->trl, "TR2") == 0)
    fprintf(fp, "  Technology concept formulated.\n");
  if (strcmp(params->trl, "TR3") == 0)
    fprintf(fp, "  Experimental proof of concept.\n");
  if (strcmp(params->trl, "TR4") == 0)
    fprintf(fp, "  Technology validated in laboratory environment.\n");
  if (strcmp(params->trl, "TR5") == 0)
    fprintf(fp, "  Technology validated in relevant environment.\n");
  if (strcmp(params->trl, "TR6") == 0)
    fprintf(fp, "  Technology demonstrated in relevant environment (pilot model).\n");
  if (strcmp(params->trl, "TR7") == 0)
    fprintf(fp, "  System prototype demonstration in operational environment (prototype model).\n");
  if (strcmp(params->trl, "TR8") == 0)
    fprintf(fp, "  System complete and qualified (calibrated model).\n");
  if (strcmp(params->trl, "TR9") == 0)
    fprintf(fp, "  Actual system proven in operational environment (operational model).\n");
  if (strlen(params->sequence)) fprintf(fp, "Run # %s\n", params->sequence);
  if (strlen(params->rev)) fprintf(fp, "Parameter file revision : %s\n",
				   params->rev);
  fprintf(fp, "Grid description : %s\n\n", params->grid_desc);

  if (strlen(params->notes))
    fprintf(fp, "Version notes   :  %s\n", params->notes);
  time(&t);
  fprintf(fp, "Simulation start time :  %s\n", ctime(&t));

  gettimeofday(&tm1, NULL);
  params->initt = tm1.tv_sec + tm1.tv_usec * 1e-6 - params->initt;
  if (params->initt < 3600.0)
    fprintf(fp, "Time taken for initialisation = %5.2f sec\n\n", params->initt);
  else
    fprintf(fp, "Time taken for initialisation = %5.2f hr\n\n", params->initt / 3600.0);

  /*-----------------------------------------------------------------*/
  /* Transport mode                                                  */
  if (params->tmode & SP_DUMP) 
    fprintf(fp, "\nTRANSPORT dump mode: Transport file created for forcing dumps.\n");
  if (params->runmode & TRANS ) {
    if (params->tmode & SP_FFSL) 
      fprintf(fp, "\nFFSL TRANSPORT: Data input from exact transport files.\n");
    else if (params->tmode & SP_EXACT)
      fprintf(fp, "\nTRANSPORT MODE: Data input from exact transport files.\n");
    else if (params->tmode & SP_TINT)
      fprintf(fp, "\nTRANSPORT MODE: Data input temporally interpolated from transport files.\n");
    else if (params->tmode & XYZ_TINT) {
      fprintf(fp, "\nTRANSPORT MODE: Data input interpolated from geographic files.\n");
      fprintf(fp, "                (velocities assumed to be easterly (u) and northerly (v))\n");
    } else if (params->tmode & NONE)
      fprintf(fp, "\nTRANSPORT MODE: Data IO only performed.\n");
    else if (params->tmode & SP_ORIGIN)
      fprintf(fp, "\nTRANSPORT MODE: Streamline origin input mode.\n");
    else if (params->tmode & SP_CHECK)
      fprintf(fp, "\nTRANSPORT MODE: Input data checking only performed.\n");
    if (!(params->tmode & SP_STRUCT) && !(params->tmode & SP_UGRID))
      fprintf(fp, "  Transport data in UGRID3 format.\n");
    if (params->tmode & SP_STRUCT)
      fprintf(fp, "  Transport data in 'sparse' format.\n");
    if (params->tmode & SP_UGRID)
      fprintf(fp, "  Transport data in UGRID layered tologogy format.\n");
    if (params->tmode & TR_CHECK)
      fprintf(fp, "  Transport data checking performed.\n");
    if (params->fillf & GLOBAL) 
      fprintf(fp, "Global filling invoked - non monotonic\n");
    if (params->fillf & WEIGHTED)
      fprintf(fp, "Weighted additive global filling invoked\n");
    if (params->fillf & OBC_ADJUST)
      fprintf(fp, "Open boundary flux adjustment performed\n");
    if (params->fillf & DIAGNOSE)
      fprintf(fp, "  Filling diagnostics written to file 'trans.ts'\n");
    if (params->fillf & DIAGNOSE_BGC)
      fprintf(fp, "  BGC filling diagnostics written to file 'trans.ts'\n");
    if (params->conserve & CONS_W && !(params->conserve & CONS_WS))
      fprintf(fp, "Vertical velocity computed inversely\n");
    if (params->conserve & CONS_W && params->conserve & CONS_WS)
      fprintf(fp, "Vertical velocity computed inversely if Lipschitz condition not violated\n");
    if (params->conserve & CONS_ETA)
      fprintf(fp, "Sea level computed inversely\n");
    if (params->conserve & CONS_MRG)
      fprintf(fp, "Surface layer horizontal fluxes merged with layer below\n");
    if (params->conserve & NOINIT)
      fprintf(fp, "  No sea level re-initialisation from file every timestep\n");
    if (params->conserve & REINIT)
      fprintf(fp, "  Sea level re-initialised from file each timestep\n");
    if (master->tmode & DO_DZ)
      fprintf(fp, "  Cell thicknesses read from transport file.\n");
    if (params->lyear)
      fprintf(fp, "Leap years unaccounted for in input files.\n");
    fprintf(fp, "\n");

#ifdef HAVE_OMP
    fprintf(fp, "Number of OpenMP threads = %d\n", params->trans_num_omp);
#endif
  }

  /*-----------------------------------------------------------------*/
  /* Parameters and flags                                            */
  if (strlen(params->regulate))
    fprintf(fp, "Run regulation file = %s\n", params->regulate);
  if (params->mode2d)
    fprintf(fp, "Operating in 2D mode\n");
  else
    fprintf(fp, "Operating in 3D mode\n");
  fprintf(fp, "3D time step = %5.3f\n", params->grid_dt);
  fprintf(fp, "2D time step = %5.3f\n", params->grid_dt / params->iratio);
  if (params->tratio)
    fprintf(fp, "Tracer time step = %5.3f\n",
            params->grid_dt * params->tratio);
  d2 = d4 = 1e10;
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    d1 = 1.0 / (4.0 * geom->hacell[1][c] * geom->hacell[1][c]) +
         1.0 / (4.0 * geom->hacell[2][c] * geom->hacell[2][c]);
    d3 = (1.0 / sqrt(d1)) / (maxvel + 2.0 * int_wave_speed);
    d1 =
      (1.0 / sqrt(d1)) / (maxvel +
                          2.0 * sqrt(g *
                                     fabs(master->topz[c] -
                                          geom->botz[c]) * master->Ds[c]));
    if (d1 < d2)
      d2 = d1;
    if (d3 < d4)
      d4 = d3;
  }
  /*
     fprintf(fp,"Maximum allowable 3D time step = %4.1f\n",d4);
     fprintf(fp,"Maximum allowable 2D time step = %4.1f\n",d2); */
  if (params->nwindows > 1) {
    fprintf(fp, "Number of windows = %d\n",params->nwindows);
    if (params->win_type & STRIPE_E1)
      fprintf(fp, "  Striped in the e1 direction\n");
    if (params->win_type & STRIPE_E2)
      fprintf(fp, "  Striped in the e2 direction\n");
    if (params->win_type & BLOCK_E1) {
      if (params->win_block)
	fprintf(fp, "  Rectangularly blocked in the e1 direction (%d)\n", params->win_block);
      else
	fprintf(fp, "  Blocked in the e1 direction\n");
    }
    if (params->win_type & BLOCK_E2) {
      if (params->win_block)
	fprintf(fp, "  Rectangularly blocked in the e2 direction (%d)\n", params->win_block);
      else
	fprintf(fp, "  Blocked in the e2 direction\n");
    }
    if (params->runmode & (DUMP|MANUAL)) {
      if (params->map_type & GEOM_DUMP && strlen(params->geom_file))
	fprintf(fp, "  Geometry map dumped to file %s\n", params->geom_file);
      if (params->map_type & WIN_DUMP && strlen(params->wind_file))
	fprintf(fp, "  Window map dumped to file %s\n", params->wind_file);
    }
    if (params->runmode & MANUAL) {
      if (params->map_type & GEOM_READ && strlen(params->geom_file))
	fprintf(fp, "  Geometry map read from file %s\n", params->geom_file);
      if (params->map_type & WIN_READ && strlen(params->win_file))
	fprintf(fp, "  Window map read from file %s\n", params->win_file);
    }
    if(params->win_size) {
	fprintf(fp, "  Window sizes = ");
      for (n = 0; n < params->nwindows; n++)
	fprintf(fp, "%2.1f ",params->win_size[n]);
      fprintf(fp,"\n");
    }
    fprintf(fp, "\n");
  }
#if defined(HAVE_OMP)
  if (strcasecmp(params->dp_mode,"openmp") == 0)
    fprintf(fp, "Distributed processing : OpenMP\n");
#endif
#if defined(HAVE_PTHREADS)
  if (strcasecmp(params->dp_mode,"pthreads") == 0)
    fprintf(fp, "Distributed processing : pthreads\n");
#endif
#if defined(HAVE_ECOLOGY_MODULE)
  if (params->do_eco)
    fprintf(fp, "Ecology invoked\n");
#endif
#if defined(HAVE_SEDIMENT_MODULE)
  if (params->do_sed)
    fprintf(fp, "Sediment transport invoked\n");
#endif
  if (params->compatible & V1246)
    fprintf(fp, "PRE-V1246 compatibility: global boundary cells include R_EDGE & F_EDGE OUTSIDE cells.\n");
  if (params->compatible & V1283) {
    fprintf(fp, "PRE-V1283 compatibility: Numerous bugfixes for multiple windows not included.\n");
    fprintf(fp, "   refer to Revision History Nov 16 2009, v1283-1331.\n");
  }
  if (params->compatible & V1562) {
    fprintf(fp, "PRE-V1562 compatibility: swr input.\n");
    fprintf(fp, "  swr added explicitly to the water column.\n");
  }
  if (params->compatible & V1598) {
    fprintf(fp, "PRE-V1598 compatibility: momentum advection.\n");
    if (!(master->runmode & TRANS)) {
      fprintf(fp, "  u1 and u2 = 0 above free surface for horizontal fluxes.\n");
      fprintf(fp, "  wtop uses 2D detadt & low order approximations.\n");
    } else
      fprintf(fp, "  Ghost cells & associated mappings not created for sediment.\n");
  }
  if (params->compatible & V1652)
    fprintf(fp, "PRE-V1652 compatibility: isotropic grid refinement maps used.\n");
  if (params->compatible & V1670)
    fprintf(fp, "PRE-V1670 compatibility: OBC flux adjustment eta read from relaxation file.\n");
  if (params->compatible & V1957)
    fprintf(fp, "PRE-V1957 compatibility: No ghost cells identified for open boundaries.\n");
  if (params->compatible & V4201)
    fprintf(fp, "PRE-V4201 compatibility: Define all netcdf files as netCDF classic\n");
  if (params->compatible & V5342)
    fprintf(fp, "PRE-V5342 compatibility: Turbulence closure quantities vertically diffused in both closure and vertical diffusion schemes.\n");
  if (params->compatible & V6257)
    fprintf(fp, "PRE-V6257 compatibility: Momentum tendencies added sequentially to velocity.\n");
  if (params->compatible & V6898)
    fprintf(fp, "PRE-V6898 compatibility: Set a surface no-gradient in FFSL scheme.\n");
  if (params->compatible & V7367)
    fprintf(fp, "PRE-V7367 compatibility: Using original transport scheduling.\n");

  if (params->stab & NONE)
    fprintf(fp, "No stability compensation\n");
  if (params->stab & SUB_STEP)
    fprintf(fp, "Sub-stepping stability compensation\n");
  if (params->stab & SUB_STEP_NOSURF)
    fprintf(fp,
            "Sub-stepping stability compensation; excluding surface layer\n");
  if (params->stab & SUB_STEP_TRACER)
    fprintf(fp, "Sub-stepping stability compensation; tracers only\n");
  if (params->thin_merge)
    fprintf(fp, "Thin layer adjustment implemented (HMIN=%4.2f)\n", params->hmin);
  if (!(params->fatal & NONE)) {
    if (params->fatal & ETA_A)
      fprintf(fp, "Exit on fatal eta instabilities when |eta| > %6.2f\n", 
	      params->etamax);
    if (params->fatal & VEL2D)
      fprintf(fp, "Exit on fatal vel2d instabilities when |vel| > %6.2f\n", 
	      params->velmax2d);
    if (params->fatal & VEL3D)
      fprintf(fp, "Exit on fatal vel3d instabilities when |vel| > %6.2f\n", 
	      params->velmax);
    if (params->fatal & (TS|NANF))
      fprintf(fp, "Exit on fatal T/S instabilities when T or S = NaN\n"); 
    if (params->fatal & NANF)
      fprintf(fp, "Exit when above variables = NaN\n");
  }
  if (master->dbc)
    fprintf(fp, "Debugging info for location (%d,%d,%d) written to file debug.txt\n",
	    params->dbi, params->dbj, params->dbk);

  /*-----------------------------------------------------------------*/
  /* Diagnostics                                                     */
  fprintf(fp, "\n");
  if (params->mixlayer & DENS_MIX)
    fprintf(fp,
            "Mixed layer depth diagnostic invoked (using density gradient)\n");
  if (params->mixlayer & TKE_MIX)
    fprintf(fp, "Mixed layer depth diagnostic invoked (using tke)\n");
  if (params->cfl & PASSIVE)
    fprintf(fp, "CFL time-step diagnostic invoked\n");
  if (params->cfl & ACTIVE)
    fprintf(fp, "Active CFL time-step diagnostic invoked : set at %s\n",
            params->cfl_dt);
  if (params->cfl & ACTIVE3D)
    fprintf(fp, "Active 3D CFL time-step diagnostic invoked until %s\n",
            params->cfl_dt);
  if (params->cfl & WVEL)
    fprintf(fp,
            "  Vertical velocity Courant violations included in CFL\n");
  if (strcmp(params->trflux, "NONE") != 0)
    fprintf(fp, "Tracer flux '%s' diagnostic invoked\n", params->trflux);
  if (strcmp(params->trperc, "NONE") != 0)
    fprintf(fp, "Tracer percentile '%s' diagnostic invoked\n", params->trperc);
  if (params->trflsh)
    fprintf(fp, "Flushing diagnostic invoked for tracer 'flush'\n");
  if (params->means & TRANSPORT)
    fprintf(fp, "Mean quantities used for inline transport\n");
  if (params->means & VEL3D)
    fprintf(fp, "Mean 3D velocity diagnostic invoked\n");
  if (params->means & VEL2D)
    fprintf(fp, "Mean 2D velocity diagnostic invoked\n");
  if (params->means & ETA_M)
    fprintf(fp, "Mean elevation diagnostic invoked\n");
  if (params->means & TS)
    fprintf(fp, "Mean T/S diagnostic invoked\n");
  if (params->means & WIND)
    fprintf(fp, "Mean wind diagnostic invoked\n");
  if (params->means & TENDENCY)
    fprintf(fp, "Mean momentum tendency diagnostic invoked\n");
  if (params->means & (MTRA2D|MTRA3D))
    fprintf(fp, "Mean tracer %s diagnostic invoked\n", params->means_tra);
  if (params->means & FLUX) {
    if (strcmp(params->trflux, "NONE") != 0)
    	fprintf(fp, "Mean flux diagnostic invoked for tracer '%s'\n",
		params->trflux);
    else /*UR-ADDED 12/2005*/
      hd_quit
        ("write_run_setup: Mean Flux is selected but no tracer is provided (calc_flux #) exiting... \n");

  }
  if (!(params->means & NONE))
    if (params->means & TRANSPORT)
      fprintf(fp, "  Averaging period = %3.0f timesteps\n", params->tratio);
    else
      fprintf(fp, "  Averaging period = %s\n", params->means_dt);
  if (params->lnm != 0.0)
    fprintf(fp, "Steric height diagnostic invoked : lnm = -%6.2f\n",
            params->lnm);
  if (!(params->decf & NONE)) {
    if (params->decf & DEC_ETA) 
      fprintf(fp, "Decorrelation length scale diagnosed on eta with sample size %d %s.\n", 
	      params->decorr, params->decs);
    if (params->decf & DEC_U1) 
      fprintf(fp, "Decorrelation length scale diagnosed on u1 velocity with sample size %d %s.\n", 
	      params->decorr, params->decs);
    if (params->decf & DEC_U2) 
      fprintf(fp, "Decorrelation length scale diagnosed on u2 velocity with sample size %d %s.\n", 
	      params->decorr, params->decs);
    if (params->decf & DEC_TRA) 
      fprintf(fp, "Decorrelation length scale diagnosed on tracer %s with sample size %d %s.\n", 
	      params->decv, params->decorr, params->decs);
  }
  if (params->vorticity & ABSOLUTE)
    fprintf(fp, "Absolute vorticity diagnostic invoked\n");
  if (params->vorticity & RELATIVE)
    fprintf(fp, "Relative vorticity diagnostic invoked\n");
  if (params->vorticity & POTENTIAL)
    fprintf(fp, "Potential vorticity diagnostic invoked\n");
  if (!(params->u1_f & NONE)) {
    fprintf(fp, "u1 momentum terms omitted :\n");
    if (params->u1_f & ADVECT)
      fprintf(fp, "  ADVECT\n");
    if (params->u1_f & HDIFF)
      fprintf(fp, "  HORIZONTAL DIFFUSION\n");
    if (params->u1_f & VDIFF)
      fprintf(fp, "  VERTICAL DIFFUSION\n");
    if (params->u1_f & CORIOLIS)
      fprintf(fp, "  CORIOLIS\n");
    if (params->u1_f & PRESS_BT)
      fprintf(fp, "  BAROTROPIC PRESSURE\n");
    if (params->u1_f & PRESS_BC)
      fprintf(fp, "  BAROCLINIC PRESSURE\n");
  }
  if (!(params->u1av_f & NONE)) {
    fprintf(fp, "u1av momentum terms omitted :\n");
    if (params->u1av_f & ADVECT)
      fprintf(fp, "  ADVECT\n");
    if (params->u1av_f & HDIFF)
      fprintf(fp, "  HORIZONTAL DIFFUSION\n");
    if (params->u1av_f & VDIFF)
      fprintf(fp, "  VERTICAL DIFFUSION\n");
    if (params->u1av_f & CORIOLIS)
      fprintf(fp, "  CORIOLIS\n");
    if (params->u1av_f & PRESS_BT)
      fprintf(fp, "  BAROTROPIC PRESSURE\n");
    if (params->u1av_f & PRESS_BC)
      fprintf(fp, "  BAROCLINIC PRESSURE\n");
  }
  if (master->alertf & PASSIVE)
    fprintf(fp, "Passive alert tracking written to file %s.txt\n",
	    master->alertname);
  if (master->alertf & ACTIVE)
  {
    fprintf(fp, "Active alert tracking written to file %s.txt\n",
	    master->alertname);
    if (master->eta_f) {
      if (master->etadiff)
	fprintf(fp, "  Implement alert action on elevation : eta difference = %5.2f\n",
		master->etadiff);
	else
	  fprintf(fp, "  Implement alert action on elevation : etamax = %5.2f\n",
		  master->etamax);
      if(master->etarlx & ALERT) {
	fprintf(fp, "    Elevation relaxation performed.\n");
	fprintf(fp, "      Input data %s read in every %s\n", params->etarlxn,
		otime(params->etarlxdt, bname));
      }
    }
    if (master->vel2d_f)
      fprintf(fp, "  Implement alert action on 2d velocity : velmax = %5.2f\n",
	      master->velmax2d);
    if (master->vel3d_f)
      fprintf(fp, "  Implement alert action on 3d velocity : velmax = %5.2f\n",
	      master->velmax);
    if (master->wvel_f)
      fprintf(fp, "  Implement alert action on vertical velocity : wmax = %5.2f\n",
	      master->wmax);
    if (master->div2d_f)
      fprintf(fp, "  Implement alert action on 2d divergence : detamax = %5.3e\n",
	      master->detamax);
    if (master->div2d_f)
      fprintf(fp, "  Implement alert action on 3d divergence : dwmax = %5.3e\n",
	      master->dwmax);
    if (master->tend_f) {
      fprintf(fp, "  Implement alert action on tendencies : maximums\n");
      fprintf(fp, "    advection = %5.3e\n", master->amax);
      fprintf(fp, "    horizontal diffusion = %5.3e\n", master->hmax);
      fprintf(fp, "    vertical diffusion = %5.3e\n", master->vmax);
      fprintf(fp, "    barotropic pressure = %5.3e\n", master->btmax);
      fprintf(fp, "    baroclinic pressure = %5.3e\n", master->bcmax);
      fprintf(fp, "    Coriolis = %5.3e\n", master->cmax);
    }
    if (master->cfl_f)
      fprintf(fp, "  Implement alert action CFL condition\n");
  }
  if (master->alert_dt)
    fprintf(fp, "Alert tracking time series written to file %s.ts\n",
	    master->alertname);
  if (params->runmode & ROAM) {

    if (params->roammode == A_ROAM_CPD1)   /* Standard ROAM */
      fprintf(fp, "Standard ROAM configutration : Clamped OBCs.\n");
    if (params->roammode == A_ROAM_CPD2)
      fprintf(fp, "Standard ROAM configutration : Raymond & Kuo OBCs.\n");
    if (params->roammode == A_ROAM_FLA)
      fprintf(fp, "Standard ROAM configutration : Flather OBCs.\n");
    if (params->roammode == A_ROAM_R1)
      fprintf(fp, "ROAM configutration : DRI OBCs, robust suite #1.\n");
    if (params->roammode == A_ROAM_R2)
      fprintf(fp, "ROAM configutration : DRI OBCs, robust suite #2.\n");
    if (params->roammode == A_RECOM_R1)
      fprintf(fp, "Standard RECOM configutration.\n");
    if (params->roammode == A_RECOM_R2)
      fprintf(fp, "RECOM configutration with robust levels.\n");

    if (params->robust > 1) {
      fprintf(fp, "Robustness level %d implemented\n", params->robust);
      fprintf(fp, "  %s\n", params->robusttext);
    }
    if (params->speed > 1)
      fprintf(fp, "ROAM speed level %d implemented\n", params->speed);
  }
  if (params->porusplate) {
    if (sscanf(params->reef_frac, "%lf", &d1) == 1)
      fprintf(fp, "Porus plate reef parameterisation included using blocking = %f\n", atof(params->reef_frac));
    else 
      fprintf(fp, "Porus plate reef parameterisation included using blocking files: %s\n", master->reef_frac);
  }

  /*-----------------------------------------------------------------*/
  /* Grid structure                                                  */
  fprintf(fp, "\n");
  fprintf(fp, "Grid dimension : %d x %d x %d\n", params->nce1,
          params->nce2, params->nz);
  fprintf(fp, "Vertical structure\n");
  if (params->sigma) {
    fprintf(fp, "  Vertical coordinate system = sigma\n");
    for (n = 0; n < params->nz; n++)
      fprintf(fp, "%4.3f ", params->layers[n]);
    fprintf(fp, "%4.3f ", 0.0);
  } else {
    fprintf(fp, "  Vertical coordinate system = 'z'\n");
    for (n = 0; n < params->nz; n++)
      fprintf(fp, "%3.1f ", params->layers[n]);
  }
  fprintf(fp, "\n\n");
  if (params->sednz) {
    fprintf(fp, "  Sediment layers\n");
    for (n = 0; n <= params->sednz; n++)
      fprintf(fp, "%5.3f ", params->gridz_sed[n]);
  }
  if (params->runmode & (AUTO | DUMP)) {
    if (params->bathyfill & B_REMOVE)
      fprintf(fp, "Land values removed from the mesh.\n");
    if (params->bathyfill & B_MINIMUM)
      fprintf(fp, "Land values in the mesh replaced with the minimum value.\n");
    if (params->bathyfill & B_AVERAGE)
      fprintf(fp, "Land values in the mesh replaced with the mean surrounding values.\n");
  }
  fprintf(fp, "\n\n");

  /*-----------------------------------------------------------------*/
  /* Ramp                                                            */
  fprintf(fp, "Ramp start = %f days\n", params->rampstart / 86400.0);
  fprintf(fp, "Ramp end = %f days\n", params->rampend / 86400.0);
  fprintf(fp, "Ramp variables = ");
  if(params->rampf & WIND) fprintf(fp, "WIND ");
  if(params->rampf & FILEIN) fprintf(fp, "FILEIN ");
  if(params->rampf & CUSTOM) fprintf(fp, "CUSTOM ");
  if(params->rampf & TIDALH) fprintf(fp, "TIDALH ");
  if(params->rampf & TIDALC) fprintf(fp, "TIDALC ");
  if(params->rampf & TIDEBC) fprintf(fp, "TIDEBC ");
  if(params->rampf & FLATHR) fprintf(fp, "FLATHR ");
  if(params->rampf & INV_BARO) fprintf(fp, "INV_BARO ");
  if(params->rampf & ETA_RELAX) fprintf(fp, "ETA_RELAX ");
  if(params->rampf & FLUX_ADJUST) fprintf(fp, "FLUX_ADJUST ");
  fprintf(fp,"\n\n");

  /*-----------------------------------------------------------------*/
  /* Advection schemes                                               */
  if (params->momsc & LAGRANGE) {
    int osl = (params->osl) ? params->osl : 1;
    fprintf(fp, "Momentum advection using %d order semi-Lagrange scheme.\n", osl);
  } else {
    if (params->momsc & ORDER1)
      fprintf(fp, "1st order momentum advection scheme.\n");
    if (params->momsc & ORDER2)
      fprintf(fp, "2nd order momentum advection scheme.\n");
    if (params->momsc & VANLEER) {
      fprintf(fp, "Van Leer momentum advection scheme.\n");
    }
    if (params->momsc & RINGLER)
      fprintf(fp, "Vector invariant (with nonlinear Coriolis) momentum advection scheme.\n");
    if (params->momsc & PV_ENEUT)
      fprintf(fp, "  Vorticity computed using energy neutral formulation.\n");
    if (params->momsc & PV_ENSCO)
      fprintf(fp, "  Vorticity computed using enstrophy conserving formulation.\n");
    if (params->momsc & PV_ENSDS)
      fprintf(fp, "  Vorticity computed using enstrophy dissipating formulation.\n");
    if (params->momsc & PV_APVM)
      fprintf(fp, "  Vorticity computed using Anticipated Potential Vorticity Method (APVM).\n");
    if (params->momsc & PV_LUST)
      fprintf(fp, "  Vorticity computed using Linear-Upwind Stabilized Transport (LUST).\n");
    if (params->momsc & PV_CLUST)
      fprintf(fp, "  Vorticity computed using Continuous Linear-Upwind Stabilized Transport (CLUST).\n");
    if (params->kinetic & K_GASS)
      fprintf(fp, "  Kinetic energy gradient uses the formulation of Skamarock et al (2012).\n");
    if (params->kinetic & K_YU)
      fprintf(fp, "  Kinetic energy gradient uses the formulation of Yu et al (2020) with weighting %5.2f.\n", master->kfact);
    if (params->kinetic & ORDER1)
      fprintf(fp, "  Kinetic energy gradient computed using 1st order upwind scheme.\n");
    if (params->kinetic & ORDER2)
      fprintf(fp, "  Kinetic energy gradient computed using 2nd order centered scheme.\n");
    if (params->kinetic & ORDER2_UW)
      fprintf(fp, "  Kinetic energy gradient computed using 2nd order upwind scheme.\n");
    if (params->kinetic & ORDER4)
      fprintf(fp, "  Kinetic energy gradient computed using 4th order centered scheme.\n");
  }
  if (params->momsc & WIMPLICIT)
    fprintf(fp, "  Implicit vertical momentum advection.\n");
  if (params->momsc & ADVECT_FORM)
    fprintf(fp, "  Momentum advection scheme cast horizontally in the advection form.\n");
  if (params->momsc & WTOP_O2)
    fprintf(fp, "  wtop uses 2D detadt & low order approximations.\n");
  if (params->momsc & WTOP_O4)
    fprintf(fp, "  wtop uses 3D detadt & high order approximations.\n");
  if (params->momsc & ZERO_DRYK)
    fprintf(fp, "  u1 and u2 = 0 above free surface for horizontal fluxes.\n");
  if (params->momsc & SHAPIRO || params->filter & ADVECT)
    fprintf(fp, "  Horizontal 1st order Shapiro filtering used on advective tendency.\n");

  if (params->trasc == ORDER1)
    fprintf(fp, "1st order upwind tracer advection scheme.\n");
  if (params->trasc == ORDER2)
    fprintf(fp, "2nd order tracer advection scheme.\n");
  if (params->trasc == ORDER4)
    fprintf(fp, "4th order tracer advection scheme.\n");
  if (params->trasc == QUICKEST)
    fprintf(fp, "QUICKEST (flux form : variable grid) tracer advection scheme.\n");
  if (params->trasc == (QUICKEST|HIORDER))
    fprintf(fp, "QUICKEST (flux form : unstructured grid) tracer advection scheme.\n");
  if (params->trasc == VANLEER)
    fprintf(fp, "Van Leer tracer advection scheme (structured).\n");
  if (params->trasc == (VANLEER|HIORDER))
    fprintf(fp, "Van Leer tracer advection scheme (unstructured).\n");
  if (params->trasc == ORDER3US)
    fprintf(fp, "3rd order unstructured tracer advection scheme.\n");
  if (params->trasc == ORDER4US)
    fprintf(fp, "4th order unstructured tracer advection scheme.\n");
  if (params->trasc & FCT) {
    fprintf(fp, "Flux-corrected transport tracer advection scheme.\n");
    if (params->trasc& ORDER2)
      fprintf(fp, "  High order scheme = 2nd order.\n");
    if (params->trasc& ORDER3US)
      fprintf(fp, "  High order scheme = 3rd order unstructured.\n");
    if (params->trasc& ORDER4US)
      fprintf(fp, "  High order scheme = 4th order unstructured.\n");
  }
  fprintf(fp, "Runge-Kutta time integration = %d stages.\n", params->rkstage);
  if (params->trasc == FFSL) {
    win_priv_t *wincon = hd_data->wincon[1];
    fprintf(fp, "Flux-form semi-lagrangian tracer advection scheme.\n");
    fprintf(fp, "Velocity interpolation uses %s scheme.\n", wincon->momsr);
    if (params->fillf & CLIP) 
      fprintf(fp, "Clipping invoked to ensure monotonicity.\n");
    if (params->fillf & MONOTONIC)
      fprintf(fp, "Monotonic multiplicative global filling invoked.\n");
  }
  if (params->trasc == LAGRANGE) {
    if (params->osl & L_LINEAR)
      fprintf(fp, "Linear semi-Lagrangian tracer advection scheme.\n");
    else if (params->osl & L_SIB)
      fprintf(fp, "Sibson natural neighbours semi-Lagrangian tracer advection scheme.\n");
    else if (params->osl & L_NONSIB)
      fprintf(fp, "Non-Sibson natural neighbours semi-Lagrangian tracer advection scheme.\n");
    else if (params->osl & L_CUBIC)
      fprintf(fp, "Cubic semi-Lagrangian tracer advection scheme.\n");
    else if (params->osl & L_LSQUAD)
      fprintf(fp, "Quadratic least squares semi-Lagrangian tracer advection scheme.\n");
    else if (params->osl & L_LSLIN)
      fprintf(fp, "Linear least squares semi-Lagrangian tracer advection scheme.\n");
    else if (params->osl & L_BAYLIN)
      fprintf(fp, "Baycentric linear semi-Lagrangian tracer advection scheme.\n");
    else if (params->osl & L_BILIN)
      fprintf(fp, "Bi-linear semi-Lagrangian tracer advection scheme (for quad meshes only).\n");
    else if (params->osl == 0)
      fprintf(fp, "Semi-Lagrangian tracer advection scheme.\n");
    else if (params->osl == 1)
      fprintf(fp, "Tri-linear semi-Lagrangian tracer advection scheme.\n");
    else if (params->osl == 2)
      fprintf(fp, "Tri-quadratic semi-Lagrangian tracer advection scheme.\n");
    else if (params->osl == 3)
      fprintf(fp, "Tri-cubic semi-Lagrangian tracer advection scheme.\n");
    else if (params->osl == 4)
      fprintf(fp, "Tri-quartic semi-Lagrangian tracer advection scheme.\n");
  }
  if (params->trasc == (LAGRANGE|VANLEER))
    fprintf(fp, "Split Semi-Lagrangian/Van Leer tracer advection scheme.\n");
  if (params->trasc & FFSL && params->trasc & VANLEER)
    fprintf(fp, "Split FFSL/Van Leer tracer advection scheme.\n");
  if (params->ultimate) {
    if (params->trasc == FFSL)
      fprintf(fp, "Universal flux limiter invoked.\n");
    else
      fprintf(fp, "Ultimate filter invoked.\n");

  }
  fprintf(fp, "\n");

  /*-----------------------------------------------------------------*/
  /* Horizontal mixing                                               */
  if (params->slipprm == 1.0)
    fprintf(fp, "Free slip condition.\n");
  if (params->slipprm == 0.0)
    fprintf(fp, "Half slip condition.\n");
  if (params->slipprm == -1.0)
    fprintf(fp, "No slip condition.\n");

  if (params->diff_scale & KH_REG) {
    fprintf(fp, "Regionalized horizontal diffusivity\n");
    fprintf(fp, "  %s\n", params->u1khc);
  }
  if (params->u1kh < 0.0 && params->smagorinsky != 0.0) {
    double smag = (params->smagorinsky == 1 && params->kue1 != 1) ? params->kue1 : params->smagorinsky;
    fprintf(fp,
            "Smagorinsky horizontal diffusion; constant = %5.3f\n", smag);
    if (params->diff_scale & KH_REG)
      fprintf(fp,"  regionalized base rate applied\n");
    else if (params->bkue1) {
      double kh = 0.0;
      for (cc = 1; cc <= geom->b2_t; cc++)
	kh += (master->basek[geom->w2_t[cc]] / (double)geom->b2_t);
      fprintf(fp,"  mean base rate applied = %5.3f\n", kh);
    }
  } else {
    if (params->diff_scale & LINEAR)
      fprintf(fp, "Horizontal diffusion = %5.3f : linearly scaled\n",
            params->u1kh);
    else if (params->diff_scale & NONLIN)
      fprintf(fp, "Horizontal diffusion = %5.3f : non-linearly scaled\n",
            params->u1kh);
    else if (params->diff_scale & AREAL)
      fprintf(fp, "Horizontal diffusion = %5.3f : areal scaling\n",
            params->u1kh);
    else if (params->diff_scale & NONE && !(params->diff_scale & AUTO))
      fprintf(fp, "Horizontal diffusion = %5.3f : un-scaled\n",
            params->u1kh);
  }
  if ((params->u1kh >= 0.0 || params->bkue1 >= 0.0 || params->diff_scale & KH_REG) && 
      params->diff_scale & AUTO) {
    double kh = -1.0;
    if (params->u1kh >= 0.0) kh = floor(params->u1kh);
    if (params->bkue1 > 0.0) kh = floor(params->bkue1);
    if (kh >= 0.0)
      fprintf(fp,"Grid-optimized horizontal diffusion (%5.1f%%)\n", kh);
    else
      fprintf(fp,"Grid-optimized horizontal diffusion\n");
  }
  if (params->diff_scale & VH_REG) {
    fprintf(fp, "Regionalized horizontal viscosity\n");
    fprintf(fp, "  %s\n", params->u1vhc);
  }
  if (params->u1vh < 0.0 && params->smagorinsky != 0.0) {
    double smag = (params->smagorinsky == 1 && params->sue1 != 1) ? params->sue1 : params->smagorinsky;
    fprintf(fp,
            "Smagorinsky horizontal viscosity; constant = %5.3f\n", smag);
    if (params->diff_scale & VH_REG)
      fprintf(fp,"  regionalized base rate applied\n");
    else if (params->bsue1) {
      double vh = 0.0;
      for (cc = 1; cc <= geom->b2_t; cc++)
	vh += (master->basev[geom->w2_t[cc]] / (double)geom->b2_t);
      fprintf(fp,"  mean base rate applied = %5.3f\n", vh);
    }
  } else {
    double ah = params->u1vh;
    if (params->visc_method & US_BIHARMONIC) {
      double ahl;
      if (params->diff_scale & SCALEBI)
	ahl = ah * (0.125 * master->hmean1 * master->hmean1);
      else
	ahl = ah / (0.125 * master->hmean1 * master->hmean1);
      if (params->diff_scale & LINEAR)
	fprintf(fp, "Horizontal viscosity = %5.3f (m2s-1) ~ %5.3f (m4s-1) : linearly scaled\n", ah, ahl);
      else if (params->diff_scale & NONLIN)
	fprintf(fp, "Horizontal viscosity = %5.3f (m2s-1) ~ %5.3f (m4s-1) : non-linearly scaled\n", ah, ahl);
      else if (params->diff_scale & AREAL)
	fprintf(fp, "Horizontal viscosity = %5.3f (m2s-1) ~ %5.3f (m4s-1) : areal scaling\n", ah, ahl);
      else if (params->diff_scale & CUBIC)
	fprintf(fp, "Horizontal viscosity = %5.3f (m2s-1) ~ %5.3f (m4s-1) : cubic scaling\n", ah, ahl);
      else if (params->diff_scale & NONE && !(params->diff_scale & AUTO))
	fprintf(fp, "Horizontal viscosity = %5.3f (m2s-1) ~ %5.3f (m4s-1) : un-scaled\n", ah, ahl);
      if (params->diff_scale & SCALE2D)
	fprintf(fp, "2D horizontal viscosity scaled by IRATIO = %d\n", params->iratio);
    } else {
      if (params->diff_scale & LINEAR)
	fprintf(fp, "Horizontal viscosity = %5.3f (m2s-1) : linearly scaled\n", ah);
      else if (params->diff_scale & NONLIN)
	fprintf(fp, "Horizontal viscosity = %5.3f (m2s-1) : non-linearly scaled\n", ah);
      else if (params->diff_scale & AREAL)
	fprintf(fp, "Horizontal viscosity = %5.3f (m2s-1) : areal scaling\n", ah);
      else if (params->diff_scale & CUBIC)
	fprintf(fp, "Horizontal viscosity = %5.3f (m2s-1) : cubic scaling\n", ah);
      else if (params->diff_scale & NONE && !(params->diff_scale & AUTO))
	fprintf(fp, "Horizontal viscosity = %5.3f (m2s-1) : un-scaled\n", ah);
      if (params->diff_scale & SCALE2D)
	fprintf(fp, "2D horizontal viscosity scaled by IRATIO = %d\n", params->iratio);
    }
  }
  if ((params->u1vh >= 0.0 || params->bsue1 >= 0.0 || params->diff_scale & VH_REG) && 
      params->diff_scale & AUTO) {
    double vh = -1.0;
    if (params->u1vh >= 0.0) vh = floor(params->u1vh);
    if (params->bsue1 > 0.0) vh = floor(params->bsue1);
    if (vh >= 0.0)
      fprintf(fp,"Grid-optimized horizontal viscosity (%5.1f%%)\n", vh);
    else
      fprintf(fp, "Grid-optimized horizontal viscosity\n");
  }
  if (params->visc_method & PRE794)
    fprintf(fp, "Horizontal viscosity using Laplacian scheme (full form - pre v794)\n");
  else if (params->visc_method & LAPLACIAN)
    fprintf(fp, "Horizontal viscosity using Laplacian scheme (full form)\n");
  else if (params->visc_method & SIMPLE)
    fprintf(fp, "Horizontal viscosity using Laplacian scheme (simple form)\n");
  else if (params->visc_fact && params->visc_method & US_BIHARMONIC)
    fprintf(fp, "Horizontal viscosity using unstructured f * Laplacian scheme + (1-f) * biharmonic scheme: f = %f\n",
	    master->visc_fact);
  else if (params->visc_method & US_LAPLACIAN)
    fprintf(fp, "Horizontal viscosity using unstructured Laplacian scheme\n");
  else if (params->visc_method & US_BIHARMONIC)
    fprintf(fp, "Horizontal viscosity using unstructured biharmonic scheme\n");
  d1 = 0.01 * (master->hmean1 * master->hmean1) / master->grid_dt;
  fprintf(fp, "Mean horizontal edge length = %8.2f m (VH ~ %5.0f)\n", master->hmean1, d1);
  d1 = 0.01 * (master->hmean2 * master->hmean2) / master->grid_dt;
  fprintf(fp, "Mean horizontal distance between centres = %8.2f m (VH ~ %5.0f)\n", master->hmean2, d1);
  fprintf(fp, "Minimum horizontal distance between centres = %8.2f m\n", master->minres);
  fprintf(fp, "Maximum horizontal distance between centres = %8.2f m\n", master->maxres);
  d1 = 0.01 * master->amean / master->grid_dt;
  fprintf(fp, "Mean cell area = %5.4e m^2 (VH ~ %5.0f)\n", master->amean, d1);
  d1 = 0.01 * master->edmean / master->grid_dt;
  fprintf(fp, "Mean edge area = %5.4e m^2 (VH ~ %5.0f)\n", master->edmean, d1);
  if (params->smag_smooth) {
    fprintf(fp, "Smagorinsky clipping and smoothing (%d passes)\n", params->smag_smooth);
    fprintf(fp, "Viscosity limited to %5.2f\n", -master->u1vh0);
    fprintf(fp, "Diffusivity limited to %5.2f\n", -master->u1kh0);
  }
  fprintf(fp, "\n");

  /*-----------------------------------------------------------------*/
  /* Vertical mixing scheme                                          */
  fprintf(fp, "Vertical mixing scheme : %s\n", params->mixsc);
  if (strcmp("constant", params->mixsc) == 0) {
    fprintf(fp, "  Vertical diffusivity = %6.3e\n", params->kz0);
    fprintf(fp, "  Vertical viscosity = %6.3e\n", params->vz0);
  } else if (strcmp("csanady", params->mixsc) == 0) {
    fprintf(fp, "  Background diffusivity = %6.3e\n", params->kz0);
    fprintf(fp, "  Alpha diffusivity parameter = %6.3f\n",
            params->kz_alpha);
    fprintf(fp, "  Background viscosity = %6.3e\n", params->vz0);
    fprintf(fp, "  Alpha viscosity parameter = %6.3f\n", params->vz_alpha);
  } else if (strcmp("mellor_yamada_2_0", params->mixsc) == 0) {
    fprintf(fp, "  Surface roughness length scale = %6.3f\n", params->zs);
    fprintf(fp, "  Background diffusivity = %6.3e\n", params->kz0);
    fprintf(fp, "  Background viscosity = %6.3e\n", params->vz0);
  } else if (strcmp("mellor_yamada_2_0_estuarine", params->mixsc) == 0) {
    fprintf(fp, "  Surface roughness length scale = %6.3f\n", params->zs);
    fprintf(fp, "  Stratified layer minimum length scale = %6.3f\n",
            params->Lmin);
    fprintf(fp, "  Stratified layer stability parameter = %6.3f\n",
            params->eparam);
    fprintf(fp, "  Background diffusivity = %6.3e\n", params->kz0);
    fprintf(fp, "  Background viscosity = %6.3e\n", params->vz0);

  } else if (strcmp("k-e", params->mixsc) == 0) {
    fprintf(fp, "  Background diffusivity = %6.3e\n", params->kz0);
    fprintf(fp, "  Background viscosity = %6.3e\n", params->vz0);
    fprintf(fp, "  Surface roughness length scale = %6.3f\n", params->zs);
    if (params->waves & VERTMIX) {
      if (params->waves & JONES) {
	fprintf(fp, "Wave enhanced mixing included using Jones & Monosmith, 2008, JGR, 113.\n"); 
	if (params->wave_alpha)
	  fprintf(fp, "  Wave parameter (alpha) = %4.1e\n", params->wave_alpha);
	if (params->wave_hf)
	  fprintf(fp, "  Wave height scaling = %4.1e\n", params->wave_hf);
      }
      if (params->waves & WOM) {
	fprintf(fp, "Wave enhanced mixing included using Wave Orbital Method.\n"); 
	fprintf(fp, "Babanin & Haus, (2009), JPO, 39, 2675-2679.\n");
	fprintf(fp, "  Wave parameter b1 = %4.1e\n", params->wave_b1);
      }
      if (params->waves & BVM) {
	fprintf(fp, "Wave enhanced mixing included using BV monochromatic wave method.\n"); 
	fprintf(fp, "Qiao et al, (2009), Ocean Dynamics, 60, 1339-1355.\n");
	fprintf(fp, "  Wave parameter alpha = %4.1e\n", params->wave_alpha);
      }
    }
  } else if (strcmp("k-w", params->mixsc) == 0) {
    fprintf(fp, "  Background diffusivity = %6.3e\n", params->kz0);
    fprintf(fp, "  Background viscosity = %6.3e\n", params->vz0);
    if (params->waves & VERTMIX) {
      if (params->wave_alpha)
	fprintf(fp, "  Wave parameter (alpha) = %4.1e\n", params->wave_alpha);
      if (params->wave_hf)
	fprintf(fp, "  Wave height scaling = %4.1e\n", params->wave_hf);
    }
  } else if (strcmp("W88", params->mixsc) == 0) {
    fprintf(fp, "  Background diffusivity = %6.3e\n", params->kz0);
    fprintf(fp, "  Background viscosity = %6.3e\n", params->vz0);
    if (params->waves & VERTMIX) {
      if (params->wave_alpha)
	fprintf(fp, "  Wave parameter (alpha) = %4.1e\n", params->wave_alpha);
      if (params->wave_hf)
	fprintf(fp, "  Wave height scaling = %4.1e\n", params->wave_hf);
    }
  } else if (strcmp("mellor_yamada_2_5", params->mixsc) == 0) {
    fprintf(fp, "  Background diffusivity = %6.3e\n", params->kz0);
    fprintf(fp, "  Background viscosity = %6.3e\n", params->vz0);
  } else if (strcmp("harcourt", params->mixsc) == 0) {
    fprintf(fp, "  Background diffusivity = %6.3e\n", params->kz0);
    fprintf(fp, "  Background viscosity = %6.3e\n", params->vz0);
  }
  if (params->s_func != NULL) {
    if (strcmp(params->s_func, "CANUTO_A") == 0)
      fprintf(fp, "  Canuto et. al. (2001) A stability functions used.\n");
      master->s_func = s_canutoA;
    if (strcmp(params->s_func, "CANUTO_B") == 0)
      fprintf(fp, "  Canuto et. al. (2001) B stability functions used.\n");
    if (strcmp(params->s_func, "KANTHA&CLAYSON") == 0)
      fprintf(fp, "  Kantha & Clayson (1994) stability functions used.\n");
    if (strcmp(params->s_func, "GALPERIN") == 0)
      fprintf(fp, "  Galperin et. al. (1988) stability functions used.\n");
    if (strcmp(params->s_func, "POM") == 0)
      fprintf(fp, "  Mellor (1992) stability functions used.\n");
    if (strcmp(params->s_func, "MUNK&ANDERSON") == 0)
      fprintf(fp, "  Munk and Anderson (1948) stability functions used.\n");
    if (strcmp(params->s_func, "EIFLER&SCHRIMPF") == 0)
      fprintf(fp, "  Eifler and Schrimpf (1992) stability functions used.\n");
    if (strcmp(params->s_func, "SCHUMANN&GERZ") == 0)
      fprintf(fp, "  Schumann and Gerz (1995) stability functions used.\n");
  }
  fprintf(fp, "Bottom roughness length scale = %f\n", params->z0);
  fprintf(fp, "Mean bottom drag coefficient = %f\n", master->quad_bfc);
  fprintf(fp, "\n");

  /*-----------------------------------------------------------------*/
  /* Tracers                                                         */
  if (params->ntr) {
    fprintf(fp, "Number of 3D tracers = %d\n", params->ntr);
    for (nn = 0; nn < params->ntr; nn++) {
      tracer_info_t *tr = &master->trinfo_3d[nn];
      char *adv = (tr->advect) ? "advect" : "";
      char *dif = (ANY0(nn, master->tdif_h, master->ntdif_h)) ? "h_diffuse" : "";
      char *vdif = (ANY0(nn, master->tdif_v, master->ntdif_v)) ? "v_diffuse" : "";
      /*char *dif = (tr->diffuse) ? "diffuse" : "";*/
      char *sed = (tracer_find_index(tr->name, master->nsed, master->trinfo_sed) >= 0) ? "[sediment]" : "";
      if (ANY0(nn, master->tm_3d, master->ntm_3d))
        fprintf(fp, "  Tracer #%d : %s [%4.2e : %4.2e] : mean tracer\n",
                nn, tr->name,
                tr->valid_range_wc[0], tr->valid_range_wc[1]);
      else
        fprintf(fp, "  Tracer #%d : %s [%4.2e : %4.2e] [%s %s %s] %s\n",
                nn, tr->name,
                tr->valid_range_wc[0], tr->valid_range_wc[1], adv, dif, vdif, sed);
    }
    fprintf(fp, "\n");
  }
  if (params->ntrS) {
    fprintf(fp, "Number of 2D tracers = %d\n", params->ntrS);
    for (nn = 0; nn < params->ntrS; nn++) {
      tracer_info_t *tr = &master->trinfo_2d[nn];
      if (ANY0(nn, master->tm_2d, master->ntm_2d))
        fprintf(fp, "  Tracer #%d : %s [%4.2e : %4.2e] : mean tracer\n",
                nn, tr->name,
                tr->valid_range_wc[0], tr->valid_range_wc[1]);
      else
	fprintf(fp, "  Tracer #%d : %s [%4.2e : %4.2e]\n",
		nn, tr->name,
		tr->valid_range_wc[0], tr->valid_range_wc[1]);
    }
    fprintf(fp, "\n");
  }

  /*-----------------------------------------------------------------*/
  /* Open boundaries                                                 */
  if (geom->nobc) {
    fprintf(fp, "Number of open boundaries = %d\n", geom->nobc);
    for (n = 0; n < geom->nobc; n++) {
      open_bdrys_t *open = geom->open[n];
      fprintf(fp, "Boundary #%d : %s\n", n, open->name);
      bcname(open->bcond_nor, bname);
      fprintf(fp, "  Normal velocity = %s\n", bname);
      if (open->bcond_nor2d & (CUSTOM|FILEIN) ||
	  open->bcond_nor2d != open->bcond_nor) {
        bcname(open->bcond_nor2d, bname);
        fprintf(fp, "  2D normal velocity = %s\n", bname);
      }
      if (open->bcond_nor & (CUSTOM | FILEIN)) {
        if (open->type & U1BDRY)
          wbdrycustom(fp, params, n, &open->datau1);
	/*
        if (open->type & U2BDRY)
          wbdrycustom(fp, params, n, &open->datau2);
	*/
      }
      if (open->bcond_nor2d & (CUSTOM | FILEIN)) {
        if (open->type & U1BDRY)
          wbdrycustom(fp, params, n, &open->datau1av);
	/*
        if (open->type & U2BDRY)
          wbdrycustom(fp, params, n, &open->datau2av);
	*/
      }

      bcname(open->bcond_tan, bname);
      fprintf(fp, "  Tangential velocity = %s\n", bname);
      if (open->bcond_tan2d & (CUSTOM|FILEIN) ||
	  open->bcond_tan2d != open->bcond_tan) {
        bcname(open->bcond_tan2d, bname);
        fprintf(fp, "  2D tangential velocity = %s\n", bname);
      }
      if (open->bcond_tan & (CUSTOM | FILEIN)) {
        if (open->type & U1BDRY)
          wbdrycustom(fp, params, n, &open->datau2);
	/*
        if (open->type & U2BDRY)
          wbdrycustom(fp, params, n, &open->datau1);
	*/
      }
      if (open->bcond_tan2d & (CUSTOM | FILEIN)) {
        if (open->type & U1BDRY)
          wbdrycustom(fp, params, n, &open->datau2av);
	/*
        if (open->type & U2BDRY)
          wbdrycustom(fp, params, n, &open->datau1av);
	*/
      }
      bcname(open->bcond_ele, bname);
      fprintf(fp, "  Elevation = %s\n", bname);
      if (open->bcond_ele & (CUSTOM | FILEIN))
        wbdrycustom(fp, params, n, &open->etadata);

      bcname(open->bcond_w, bname);
      if (!(open->bcond_w & NOTHIN))
	fprintf(fp, "  Vertical velocity = %s\n", bname);
      
      for (nn = 0; nn < params->ntr; nn++) {
	scale_details_t *scale = &open->sdata_t[nn];
        bcname(open->bcond_tra[nn], bname);
        fprintf(fp, "  Tracer #%d (%s) = %s\n", nn,
                master->trinfo_3d[nn].name, bname);
        if (open->bcond_tra[nn] & (CUSTOM | FILEIN))
          wbdrycustom(fp, params, n, &open->bdata_t[nn]);
	if (strlen(scale->name))
	  fprintf(fp, "    Scaling: %s\n", scale->name);
      }

      if (params->do_wave & (W_SWAN|W_SWANM)) {
	bcname(open->bcond_wav, bname);
	fprintf(fp, "  Wave variables = %s\n", bname);
      }

      if (open->relax_time) {
	if (open->relax_time != open->relax_timei) {
	  fprintf(fp, "    Out-going relaxation constant = %6.3f (hours)\n",
		  open->relax_time / 3600.0);
	  fprintf(fp, "    In-coming relaxation constant = %6.3f (hours)\n",
		  open->relax_timei / 3600.0);
	} else
	  fprintf(fp, "    Relaxation constant = %6.3f (hours)\n",
		  open->relax_time / 3600.0);
      }
      if (open->sponge_zone_h) {
        c = open->obc_e1[1];
        fprintf(fp, "    Horizontal sponge zone : e1=%6.3f\n",
                master->u1vh[c]);
      }
      if (open->inverse_barometer)
        fprintf(fp, "    Inverse barometer compensation included\n");
      else
        fprintf(fp, "    No inverse barometer compensation\n");
      if (open->sponge_zone) {
        fprintf(fp, "    Vertical sponge zone\n");
      }
      if (open->relax_zone_nor) {
        fprintf(fp, "    Normal velocity flow relaxation zone = %d\n",
		open->relax_zone_nor);
      }
      if (open->relax_zone_tan) {
        fprintf(fp, "    Tangential velocity flow relaxation zone = %d\n",
		open->relax_zone_tan);
      }
      if (open->relax_ele) {
        fprintf(fp, "    Relaxing elevation over %d cells\n", open->relax_ele);
	fprintf(fp, "      Input data %s read in every %s\n", 
		params->etarlxn, otime(params->etarlxdt, bname));
        fprintf(fp, "      Relaxation timescale on the boundary = %s\n",
		otime(open->rele_b, bname));
        fprintf(fp, "      Relaxation timescale in the interior = %s\n",
		otime(open->rele_i, bname));
      }
      if (strlen(open->i_rule))
        fprintf(fp, "    OBC data interpolated using %s scheme\n", open->i_rule);
      if (open->file_dt)
	fprintf(fp, "    Boundary data read from file every %s\n", otime(open->file_dt, bname));
      if (open->adjust_flux) {
	if (open->adjust_flux > 0.0)
	  fprintf(fp, "    Adjusting 2D boundary flux on a timescale of %s\n", otime(open->adjust_flux, bname));
	else
	  fprintf(fp, "    Adjusting 2D boundary flux using default timescale\n");
	if (open->adjust_flux_s > 0.0)
	  fprintf(fp, "    Adjusting tidal boundary flux on a timescale of %s\n", otime(open->adjust_flux_s, bname));
	if (open->adjust_flux_s < 0.0)
	  fprintf(fp, "    Adjusting tidal boundary flux using default timescale\n");
	if (params->compatible & V1670)
	  fprintf(fp, "      Input data %s read in every %s\n", 
		  params->etarlxn, otime(params->etarlxdt, bname));
	else {
	  bcname(open->bcond_ele, bname);
	  if (open->file_dt)
	    fprintf(fp, "      eta read from boundary data using %s every %s\n",bname, otime(open->file_dt, bname));
	  else
	    fprintf(fp, "      eta read from boundary data using %s\n", bname);
	}
      }
      if (open->linear_zone_nor) {
        fprintf(fp, "    Normal velocity is linear for %d interior boundary cells\n",
		open->linear_zone_nor);
      }
      if (open->stagger & INFACE) {
        fprintf(fp, "    Boundary stagger uses interior face for normal velocity\n");
      }
      if (open->linear_zone_tan) {
        fprintf(fp, "    Tangential velocity is linear for %d interior boundary cells\n",
		open->linear_zone_tan);
      }
      if (strlen(open->tide_con))
        fprintf(fp, "    Custom tidal constituents used: %s\n", open->tide_con);
      if (open->spf) {
        fprintf(fp, "    Phase speed smoothing using factor %3.1f applied to elevation\n", open->spf);
      }
      for (nn = 0; nn < params->ntr; nn++) {
        bcname(open->bcond_tra[nn], bname);
	if (open->relax_zone_tra[nn]) {
	  fprintf(fp, "    Tracer %s flow relaxation zone = %d\n",
		  bname,open->relax_zone_tra[nn]);
	}
      }
      if (open->relax_zone_ele) {
        fprintf(fp, "    Elevation flow relaxation zone = %d\n",
		open->relax_zone_ele);
      }
      if (!(open->upmeth & INTERIOR)) {
	if (open->upmeth & FACE)
	  fprintf(fp, "    Cell face velocity used in UPSTRM OBC.\n");
	if (open->upmeth & CENTER)
	  fprintf(fp, "    Cell center velocity used in UPSTRM OBC.\n");
	if (open->upmeth & ADAPTIVE)
	  fprintf(fp, "    Adaptively selected velocity used in UPSTRM OBC.\n");
      }
      for (nn = 0; nn < open->ntsfiles; ++nn)
        fprintf(fp, "    Boundary data file #%d : %s\n", nn,
                open->filenames[nn]);
      fprintf(fp, "  Boundary ghost zone = %d cells\n", open->bgz);
      if (open->options & OP_UPSTRM)
	fprintf(fp, "    Upstream advection applied to TRCONC/TRCONF OBC ghost cells.\n");
      if (open->options & OP_GEOSTR)
	fprintf(fp, "    Geostrophic flow computed for RIVER OBCs.\n");
      if (open->options & OP_DYNAHC)
	fprintf(fp, "    RIVER pycnocline depth computed dynamically.\n");
      if (open->options & OP_NOHDIF)
	fprintf(fp, "    No horizontal diffusion applied in boundary cell.\n");
      if (open->options & OP_YANKOVSKY)
	fprintf(fp, "    River OBCs use Yankovsky (2000).\n");
      if (open->options & OP_NOSALT)
	fprintf(fp, "    Inflow salinity is not modified using a salt mass balance.\n");
      if (open->options & OP_FDEPTH)
	fprintf(fp, "    Inflow velocity computed from flow assuming full depth distribution.\n");
      if (open->options & OP_IFRESH)
	fprintf(fp, "    Density of inflow assumed to be 1000.\n");
      if (open->options & OP_TRUNCL)
	fprintf(fp, "    Inflow profile is truncated to next deepest layer.\n");
      if (open->options & OP_MULTF)
	fprintf(fp, "    Inflow profile multiplicatively scaled to flow.\n");
      if (open->options & OP_ETAFIL)
	fprintf(fp, "    Boundary sea level median filtered at startup.\n");
      if (open->options & OP_OBCCM)
	fprintf(fp, "    Boundary overlap with other OBCs averaged for eta, T/S.\n");
      if (open->options & OP_ISPNG)
	fprintf(fp, "    Boundary sponges applied in both e1 and e2 directions.\n");
      if (open->options & OP_PARFLOW)
	fprintf(fp, "    River flow delivered with unmodified parabolic profile.\n");
      if (open->options & OP_MACREADY)
	fprintf(fp, "    Inflow salinity is modified using MacCready & Geyer, 2010, Annu. Rev. Mar. Sci.\n");
      if (open->options & OP_OWRITE)
	fprintf(fp, "    Open boundary location overwritten with external data using TRCONC.\n");
      if (open->options & OP_TILED) {
	fprintf(fp, "    Open boundary data exchange configured for 2-way tiled nesting.\n");
	if (master->obcf & DF_BARO)
	  fprintf(fp, "      Coupled at the barotropic level.\n");
      }
      fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    write_bdry(params, geom);
  }

  /*-----------------------------------------------------------------*/
  /* Source/sinks                                                    */
  if (master->npss) {
    fprintf(fp,"Number of point sourcesinks %d\n", master->npss);
    for (n = 0; n < master->npss; n++) {
      pss_t *p = &master->pss[n];

      if (p->v_offset == 1)
	strcpy(bname,"bottom");
      else if (p->v_offset == 0)
	strcpy(bname,"msl");
      else if (p->v_offset == 2)
	strcpy(bname,"surface");
      if (p->watertsid >= 0)
	if (p->vc == 1)
	  fprintf(fp,"%d %s : flow, %s referenced, loc=(%5.2f, %5.2f) at (%d %d)\n",n, p->name, bname, 
		  p->x, p->y, geom->s2i[p->e1[0]], geom->s2j[p->e1[0]]);
	else
	  fprintf(fp,"%d %s : flow, %s referenced, plane input\n",n, p->name, bname);
      else
	if (p->vc == 1)
	  fprintf(fp,"%d %s : flux, %s referenced, loc=(%5.2f, %5.2f) at (%d %d)\n",n, p->name, bname, 
		  p->x, p->y, geom->s2i[p->e1[0]], geom->s2j[p->e1[0]]);
	else
	  fprintf(fp,"%d %s : flux, %s referenced, plane input\n",n, p->name, bname);

      if (p->flag & PSS_AW)
	fprintf(fp,"  Area weighted input.\n");
      if (p->flag & PSS_VW)
	fprintf(fp,"  Volume weighted input.\n");

      if(p->vc > 1) {
	double zhigh_i;
	double zlow_i;
	double pdz;
	int cs, cb, c2, zp1;
	
	d1 = d2 = d3 = 0.0;
	for (cc = 0; cc < p->vc; cc++) {
	  c = p->e1[cc];
	  cs = p->e2[cc];
	  cb = p->e3[cc];
	  zhigh_i = p->zhigh;
	  zlow_i  = p->zlow;
	  ref_depth_m(master, p, &zlow_i, &zhigh_i, cs);
	  c2 = geom->m2d[c];
	  zp1 = geom->zp1[c];	    
	  if (c == 0) continue;
	  
	  zlow_i = max(zlow_i, geom->botz[c2]);
	  zhigh_i = min(zhigh_i, master->eta[c2]);
	  pdz = zhigh_i - zlow_i;

	  if (pdz > 0.0) {
	    c = c2;
	    while (c != geom->zm1[c]) {
	      double ctop = (c == cs) ? master->eta[c2] : geom->gridz[zp1];
	      double cbot = (c == cb) ? geom->botz[c2] : geom->gridz[c];
	      double zlow;
	      double zhigh;
	      double frac;
	      if (zlow_i <= cbot)
		zlow = cbot;
	      else if (zlow_i >= ctop)
		zlow = ctop;
	      else
		zlow = zlow_i;
	      if (zhigh_i <= cbot)
		zhigh = cbot;
	      else if (zhigh_i >= ctop)
		zhigh = ctop;
	      else
		zhigh = zhigh_i;
		
	      frac = (zhigh - zlow);
	      if (frac) {
		d1 += geom->cellarea[c2];
	      }
	      zp1 = c;
	      c = geom->zm1[c];
	    }
	  }
	}
	fprintf(fp,"  horizontal planar area = %5.2f km^2\n", d1/1e6);
      }
    }
    fprintf(fp,"\n");
  }

  /*-----------------------------------------------------------------*/
  /* Initial conditions                                              */
  if (params->runmode & AUTO) {
    if (strlen(params->eta_init))
      fprintf(fp, "Sea level initialisation from file %s\n", params->eta_init);
    if (strlen(params->vel_init)) {
      if (strcmp(params->vel_init, "GEOSTROPIC") == 0)
	fprintf(fp, "Geostrophic velocity initialisation\n");
      else
	fprintf(fp, "Velocity initialisation from file %s\n", params->vel_init);
    }
    fprintf(fp,"\n");
  }

  /*-----------------------------------------------------------------*/
  /* Forcing                                                         */
  fprintf(fp, "Wind forcing from file %s\n", params->wind);
  fprintf(fp, "  Wind speed scale = %5.2f\n", params->wind_scale);
  if (params->wind_type & SPEED)
    fprintf(fp, "  Wind speed input\n");
  else
    fprintf(fp, "  Wind stress input\n");
  if (params->stress_fn & ORIGINAL) {
    fprintf(fp, "  Wind speed threshold #1 = %5.2f\n", params->dlv0);
    fprintf(fp, "  Wind speed threshold #2 = %5.2f\n", params->dlv1);
    fprintf(fp, "  Surface drag coefficient #1 = %8.5f\n", params->dlc0);
    fprintf(fp, "  Surface drag coefficient #2 = %8.5f\n", params->dlc1);
  } else if (params->stress_fn == LARGEPOND) {
    fprintf(fp, "  Bulk scheme = Large and Pond (1982)\n");
    fprintf(fp, "  Reference height for wind = %5.2f\n", params->dlv0);
  } else if (params->stress_fn == BUNKER) {
    fprintf(fp, "  Bulk scheme = Bunker (1976)\n");
    fprintf(fp, "  Reference height for wind = %5.2f\n", params->dlv0);
  } else if (params->stress_fn == KITIAG) {
    fprintf(fp, "  Bulk scheme = Kitiagorodskii et al (1973)\n");
    fprintf(fp, "  Reference height for wind = %5.2f\n", params->dlv0);
  } else if (params->stress_fn == KONDO) {
    fprintf(fp, "  Bulk scheme = Kondo (1975)\n");
    fprintf(fp, "  Reference height for wind = %5.2f\n", params->dlv0);
    if (params->neutral || master->sh_f == NONE)
      fprintf(fp, "  Drag computed for neutral conditions.\n");
  }
  fprintf(fp, "\n");

  if (!(params->heatflux & NONE)) {
    fprintf(fp, "Heat flux calculated\n");
    if (params->heatflux & GHRSST)
      fprintf(fp, "  Relaxation to surface GHRSST temperature.\n");
    if (params->heatflux & (SURF_RELAX|AVHRR)) {
      fprintf(fp, "  Relaxation to surface temperature.\n");
      if (master->hftemp)
	fprintf(fp, "  Surface temperature file : %s\n", params->hftemp);
      else
	fprintf(fp, "  SST derived from AVHRR diagnostic.\n");
    } else {
      if (params->heatflux & (NET_HEAT|COMP_HEAT|COMP_HEAT_MOM)) {
        fprintf(fp, "Net heat flux file input\n");
	if (params->heatflux & COMP_HEAT)
	  fprintf(fp, "  Net heat flux assembled from components\n");
	if (params->heatflux & COMP_HEAT_MOM) {
	  fprintf(fp, "  Net heat flux assembled from MOM4 output components.\n");
	  fprintf(fp, "    swr = daily mean\n");
	  fprintf(fp, "    lhf = evaporation in kg/m2/s\n");
	}
        fprintf(fp, "  Net heat flux file = %s \n", params->hf);
        if (master->swr_attn) {
          fprintf(fp, "  Depth dependent short wave radiation implemented.\n");
	  if (params->water_type == TYPE_I)
	    fprintf(fp, "  Water type = I\n");
	  if (params->water_type == TYPE_IA)
	    fprintf(fp, "  Water type = IA\n");
	  if (params->water_type == TYPE_IB)
	    fprintf(fp, "  Water type = IB\n");
	  if (params->water_type == TYPE_II)
	    fprintf(fp, "  Water type = II\n");
	  if (params->water_type == TYPE_III)
	    fprintf(fp, "  Water type = III\n");
          fprintf(fp, "  Short wave attenuation = %4.2f\n", master->swr_attn[1]);
          if(master->swr_attn1) {
            fprintf(fp, "  Deep short wave attenuation = %4.2f\n",
              master->swr_attn1[1]);
            fprintf(fp, "  Short wave shallow - deep fraction = %4.2f\n",
              master->swr_tran[1]);
          } else {
            fprintf(fp, "  Short wave transmission = %4.2f\n",
              master->swr_tran[1]);
          }
        }
        if (params->albedo > 0)
          fprintf(fp, "  Short wave radiation albedo = %4.2f\n", params->albedo);
        if (master->swr_babs)
          fprintf(fp, "  Short wave radiation bottom absorption parameter = %4.2f\n", master->swr_babs[1]);
      }
      if (params->heatflux & ADVANCED) {
        fprintf(fp, "Heat flux calculated : bulk formulation\n");
        if (params->bulkf == KONDO)
          fprintf(fp, "  Bulk scheme = Kondo (1975)\n");
        if (params->bulkf == LARGEPOND)
          fprintf(fp, "  Bulk scheme = Large and Pond (1982)\n");
        if (params->bulkf == MASAG)
          fprintf(fp, "  Bulk scheme = Masagutov (1981)\n");
        if (params->bulkf == KITIAG)
          fprintf(fp, "  Bulk scheme = Kitiagorodskii et al (1973)\n");
        if (params->bulkf == BUNKER)
          fprintf(fp, "  Bulk scheme = Bunker (1976)\n");
        fprintf(fp,
                "  Reference height for air temperature/humidity = %5.2f\n",
                params->zref);
        if (master->swr_attn) {
          fprintf(fp, "  Depth dependent short wave radiation implemented.\n");
	  if (params->water_type == TYPE_I)
	    fprintf(fp, "  Water type = I\n");
	  if (params->water_type == TYPE_IA)
	    fprintf(fp, "  Water type = IA\n");
	  if (params->water_type == TYPE_IB)
	    fprintf(fp, "  Water type = IB\n");
	  if (params->water_type == TYPE_II)
	    fprintf(fp, "  Water type = II\n");
	  if (params->water_type == TYPE_III)
	    fprintf(fp, "  Water type = III\n");
          fprintf(fp, "  Short wave attenuation = %4.2f\n", master->swr_attn[1]);
	  if(master->swr_attn1) {
	    fprintf(fp, "  Deep short wave attenuation = %4.2f\n",
		    master->swr_attn1[1]);
	    fprintf(fp, "  Short wave shallow - deep fraction = %4.2f\n",
		    master->swr_tran[1]);
	  } else {
	    fprintf(fp, "  Short wave transmission = %4.2f\n",
		    master->swr_tran[1]);
	  }
        }
        if (master->swr_babs)
          fprintf(fp, "  Short wave radiation bottom absorption parameter = %4.2f\n", master->swr_babs[1]);
	if (strlen(params->swr_regions)) {
	  double atn0, atns, atni, trn0, trns, trni;
	  /* Same algorithm as in swr_params_event()                 */
	  get_swr_ensemble(params->swr_ens, &atn0, &atns, &atni,
			   &trn0, &trns, &trni);
	  fprintf(fp, "  Short wave radiation parameter estimation:\n");
	  fprintf(fp, "    Regions = %s\n", params->swr_regions);
	  fprintf(fp, "    Update dt = %5.1f hours\n", params->swreg_dt / 3600.0);
	  fprintf(fp, "    Data = %s\n", params->swr_data);
	  fprintf(fp, "    attn ensemble = [%3.2f:%3.2f:%3.2f]\n", 
		  atn0, atni, atns);
	  fprintf(fp, "    tran ensemble = [%3.2f:%3.2f:%3.2f]\n", 
		  trn0, trni, trns);
	}
      }
      fprintf(fp, "\n");
      if (master->airtemp)
        fprintf(fp, "  Air temperature file : %s\n", params->airtemp);
      if (master->wetb)
        fprintf(fp, "  Wet bulb temperature file : %s\n", params->wetb);
      if (strlen(params->swr))
        fprintf(fp, "  Short wave radiation file : %s\n", params->swr);
      if (master->cloud)
        fprintf(fp, "  Cloud cover file : %s\n", params->cloud);
      if (master->patm)
        fprintf(fp, "  Atmospheric pressure file : %s\n", params->patm);
      if (master->rh)
        fprintf(fp, "  Relative humidity file : %s\n", params->rh);
      if (master->precip)
        fprintf(fp, "  Precipitation file : %s\n", params->precip);
      if (master->evap)
        fprintf(fp, "  Evaporation file : %s\n", params->evap);
    }
    if (params->hf_ramp != params->t)
      fprintf(fp, "  Heat flux applied after %5.1f days\n", params->hf_ramp / 86400.0);
    fprintf(fp, "\n");
  } else
    fprintf(fp, "No heat flux specified\n\n");

  if (!(params->saltflux & NONE)) {
    fprintf(fp, "Salt flux calculated\n");
    if (params->saltflux & BULK) {
      fprintf(fp, "  Evaporation estimated using bulk latent heat\n");
    }
    fprintf(fp,"\n");
   } else
    fprintf(fp, "No salt flux specified\n\n");

  if(master->etarlx & (RELAX|ETA_TPXO)) {
    relax_info_t *rlx = master->eta_rlx;
    fprintf(fp, "Elevation relaxation performed.\n");
    fprintf(fp, "  Input data %s read in every %s\n", params->etarlxn,
	    otime(params->etarlxdt, bname));
    if (rlx->tctype & RLX_CONS)
      fprintf(fp, "  Relaxation time constant = %s\n",
	      otime(params->etarlxtc, bname));
    if (rlx->tctype & RLX_FILE)
      fprintf(fp, "  Relaxation time constant read from file %s\n", params->etarlxtcs);
    if (rlx->tctype & RLX_ADPT) {
      if (rlx->tctype & RLX_LINR)
	fprintf(fp, "  Linear adaptive time constant : rate = (|deta| - %4.2f) * %4.2f + %4.2f\n", rlx->dv0, rlx->slope / 86400.0, rlx->tc0 / 86400.0);
      if (rlx->tctype & RLX_EXP)
	fprintf(fp, "  Exponential adaptive time constant : rate = exp(%4.2f/|deta|)\n", rlx->slope);
      if (rlx->tctype & RLX_TIME)
	fprintf(fp, "  Temporal linear time constant : rate = (t - %4.2f) * %4.2f + %4.2f\n", rlx->dv0, rlx->slope / 86400.0, rlx->tc0 / 86400.0);
      if (rlx->tctype & RLX_DEP)
	fprintf(fp, "  Depth scaled linear time constant : rate = (depth - %4.2f) * %4.2f + %4.2f\n", rlx->dv0, rlx->slope / 86400.0, rlx->tc0 / 86400.0);
      if (rlx->tctype & RLX_EDEP) {
	fprintf(fp, "  Exponential depth scaled relaxation time constant\n");
	fprintf(fp, "    rate = (r0 - r1)exp(-d/d0) + (r1 - r0.exp(-d1/d0))\n");
	fprintf(fp, "    r0 = %5.1f at d = 0, r1 = %5.1f at d = %5.1f\n", rlx->tc0 / 86400.0, rlx->tc1 / 86400.0, rlx->dv1); 
	fprintf(fp, "    Constant d0 = %5.1f\n", rlx->dv0);
      }
      if (rlx->tctype & RLX_CDEP) {
	fprintf(fp, "  Cosine depth scaled relaxation time constant\n");
	fprintf(fp, "    rate = 0.5*((r0-r1)*cos(d*pi/(d1-d0)-d0*PI/(d1-d0))+(t0+t1))\n");
	fprintf(fp, "    r0 = %5.1f, r1 = %5.1f\n", rlx->tc0 / 86400.0, rlx->tc1 / 86400.0); 
	fprintf(fp, "    d0 = %5.1f, d1 = %5.1f\n", rlx->dv0, rlx->dv1); 
      }
    }
  }
  if (params->tide_r & MEAN_R && master->etam)
    fprintf(fp, "  Tide removed by relaxing to mean elevation.\n");
  if (params->tide_r & CSR_R)
    fprintf(fp, "  Tide removed using CSR tide model.\n");
  fprintf(fp, "\n");

  if(master->velrlx & RELAX) {
    relax_info_t *rlx = master->vel_rlx;
    fprintf(fp, "Velocity relaxation performed.\n");
    fprintf(fp, "  Input data %s read in every %s\n", params->velrlxn,
	    otime(params->velrlxdt, bname));
    if (rlx->tctype & RLX_CONS)
      fprintf(fp, "  Relaxation time constant = %s\n", params->velrlxtcs);
    if (rlx->tctype & RLX_FILE)
      fprintf(fp, "  Relaxation time constant read from file %s\n", params->velrlxtcs);
    if (rlx->tctype & RLX_ADPT) {
      if (rlx->tctype & RLX_LINR)
	fprintf(fp, "  Linear adaptive time constant : rate = (|dvel| - %4.2f) * %4.2f + %4.2f\n", rlx->dv0, rlx->slope / 86400.0, rlx->tc0 / 86400.0);
      if (rlx->tctype & RLX_EXP)
	fprintf(fp, "  Exponential adaptive time constant : rate = exp(%4.2f/|dvel|)\n", rlx->slope);
    }
  }
  if (params->tide_r & MEAN_R && master->u1m && master->u2m &&
      master->u1am && master->u2am)
    fprintf(fp, "  Tide removed by relaxing to mean velocities.\n");
  if (params->tide_r & CSR_R)
    fprintf(fp, "  Tide removed using CSR tide model.\n");
  fprintf(fp, "\n");

  if (strlen(params->webf)) {
    if (params->webf_dt <= 0.0)
      fprintf(fp,"Wind waves invoked using the wave library \n");
    else 
      webf_write_setup(fp, sched_get_even_by_name(schedule, "forcings:webf"));
  }
#if defined(HAVE_WAVE_MODULE)
  /* Write wave library info */
  if (!(params->do_wave & NONE)) {
    win_priv_t *wincon = hd_data->wincon[1];
    wave_run_setup(fp, wincon->wave);
  }
#endif
  if (params->waves & BOT_STR)
    fprintf(fp, "Wave enhanced bottom friction invoked\n");
  if (params->waves & TAN_RAD)
    fprintf(fp, "Tangential radiation stresses (waves) included in 2D mode\n");
  if (params->waves & WAVE_FOR)
    fprintf(fp, "Wave-induced forces included in 2D mode\n");
  if (params->waves & VERTMIX)
    fprintf(fp, "Wave enhanced vertical mixing invoked\n");
  if (params->waves & STOKES_DRIFT) {
    fprintf(fp, "Wave stokes drift velocity invoked\n");
    if (master->tau_w1 && master->tau_diss1 && master->tau_w2 && master->tau_diss2)
      fprintf(fp, "  Wind stress modified by wave to ocean stress and wave supported wind stress.\n");
    if (params->waves & NEARSHORE)
      fprintf(fp, "Nearshore (surf zone) wave processes invoked\n");
  }
  if (params->waves & STOKES_MIX && strcmp("harcourt", params->mixsc) == 0)
    fprintf(fp, "Wave enhanced vertical mixing (Langmuir) invoked\n");
  fprintf(fp,"\n");
  if (params->tidep) {
    fprintf(fp, "Tidal body force included.\n");
    fprintf(fp, "  Tidal body force constant constant = %5.2f\n\n", params->eqt_beta);
    fprintf(fp, "  Tidal self-attraction / loading (SAL) constant = %5.2f\n", params->eqt_alpha);
  }

  if (params->ghrsst_path) {
    fprintf(fp, "Remote sensing SST import activated.\n");
    fprintf(fp, " Path for SST: %s\n", params->ghrsst_path);
    if (strlen(params->ghrsst_opt))
      fprintf(fp, " SST options = %s\n", params->ghrsst_opt);
    if (strlen(params->ghrsst_irule))
      fprintf(fp, " SST interpolation rule = %s\n", params->ghrsst_irule);
    fprintf(fp, "\n");
  }

  /*-----------------------------------------------------------------*/
  /* Output files                                                    */
  if(master->dumpdata->ndf) {
    fprintf(fp,"Number of output file dumps = %d\n", master->dumpdata->ndf);
    if (strlen(params->opath))fprintf(fp,"Output path = %s\n", master->opath);
    for(n = 0; n < master->dumpdata->ndf; n++) {
      fprintf(fp,"Output file #%d : %s",n,master->dumpdata->dumplist[n].name);
      if (master->dumpdata->dumplist[n].compress)
	fprintf(fp," (compressed)");
      if (master->dumpdata->dumplist[n].filter)
	fprintf(fp," (filter: %s)",  master->dumpdata->dumplist[n].filter->name);
      fprintf(fp,"\n");
      /*
      if (params->d_filetype != NULL) {
	if (strcmp(params->d_filetype[n], "memory") == 0) {
	  fprintf(fp,"  2-way nesting output: tinc=%s, sync=%s\n",params->d_tinc[n], params->d_sync[n]);
	}
      }
      */
      if (master->dumpdata->dumplist[n].type != NULL) {
	if (strcmp(master->dumpdata->dumplist[n].type, "memory") == 0) {
	  char sinc[MAXSTRLEN];
	  fprintf(fp,"  2-way nesting output: tinc=%s, sync=%s\n",
		  otime(master->dumpdata->dumplist[n].tinc, bname), 
		  otime(master->dumpdata->dumplist[n].sinc, sinc));
	}
      }
    }
    fprintf(fp,"\n");
  }

  /* Time series observation comparison                              */
  for (n = 0; n < nts; ++n) {
    char buf[MAXSTRLEN];
    if (tslist[n].ndata) {
      fprintf(fp,"Time series observation comparison done using file:\n");
      for (nn = 0; nn < tslist[n].ndata; nn++)
	fprintf(fp,"  %s\n", tslist[n].tsdata[nn]->name);
      if (tslist[n].metric & TS_CLOS)
	strcpy(buf, "CLOSEST");
      if (tslist[n].metric & TS_DIFF)
	strcpy(buf, "DIFF");
      if (tslist[n].metric & TS_MEAN)
	strcpy(buf, "MEAN");
      if (tslist[n].metric & TS_RMSE)
	strcpy(buf, "RMSE");
      if (tslist[n].metric & TS_CAT)
	strcpy(buf, "CATEGORICAL");
      if (tslist[n].metric & TS_HK)
	strcpy(buf, "TRUE_SKILL");
      if (tslist[n].metric & TS_TS)
	strcpy(buf, "CRITICAL_SUCCESS");
      if (tslist[n].metric & TS_H)
	strcpy(buf, "HIT_RATE");
      if (tslist[n].metric & TS_F)
	strcpy(buf, "FALSE_ALARM_RATE");
      if (tslist[n].metric & TS_PRED)
	strcpy(buf, "PREDICTION_RATE");
      for (nn = 0; nn < tslist[n].dnvars; nn++) {
	if (tslist[n].thresh == NULL)
	  fprintf(fp, "Variable %s, metric %s\n", tslist[n].dvars[nn], buf);
	else if (tslist[n].thresh[nn][0] == 0)
	  fprintf(fp, "Variable %s, metric %s, threshhold = 3D tracer %s\n", 
		  tslist[n].dvars[nn], buf, master->trinfo_3d[(int)tslist[n].thresh[nn][1]].name);
	else if (tslist[n].thresh[nn][0] == 1)
	  fprintf(fp, "Variable %s, metric %s, threshhold = 2D tracer %s\n", 
		  tslist[n].dvars[nn], buf, master->trinfo_2d[(int)tslist[n].thresh[nn][1]].name);
	else if (tslist[n].thresh[nn][0] == 2)
	  fprintf(fp, "Variable %s, metric %s, threshhold = %f\n", 
		  tslist[n].dvars[nn], buf, tslist[n].thresh[nn][1]);
      }
    }
    fprintf(fp,"\n");
  }

  /*-----------------------------------------------------------------*/
  /* Regions                                                         */
  if (geom->nregions) {
    fprintf(fp,"Number of regions = %d\n", geom->nregions);
    fprintf(fp,"Region interval = %s\n", params->region_dt);
    fprintf(fp,"Region variables : %s\n", params->region_vars);
    if (master->region_mode & RG_BDRY && master->region_obcz)
      fprintf(fp,"  Region OBCs located %d cells into the interior\n", master->region_obcz);
    if (master->region_mode & RG_AREA && master->region_obcz)
      fprintf(fp,"  OBC zones created %d cells wide\n", master->region_obcz);
    for (n = 0; n < geom->nregions; n++) {
      int ncells;
      region_t *region = geom->region[n];
      fprintf(fp,"  -------------------------\n");
      fprintf(fp,"  region%d size = %d cells\n", n, region->nvec);
      fprintf(fp,"  region fluxes\n");
      for (nn = 0; nn < region->nboundaries; nn++) {
	ncells = (region->nbz[region->bmap[nn]] ? region->nbz[region->bmap[nn]] :
		  region->nbe1[region->bmap[nn]]);
	fprintf(fp,"    %d %s (%d cells)\n", nn, region->fluxname[nn], ncells);
      }
    }
    fprintf(fp,"\n");
  }

  /*-----------------------------------------------------------------*/
  /* Process exclusion                                               */
  if (params->prex) {
    fprintf(fp, "Process exclusion at points :\n");
    for (cc = 1; cc <= params->prex; cc++) {
      fprintf(fp, "  (%d %d) : ",params->prxi[cc], params->prxj[cc]);
      if (params->prxf[cc] & EX_BGC) fprintf(fp, "BGC ");
      if (params->prxf[cc] & EX_SED) fprintf(fp, "SEDIMENTS ");
      if (params->prxf[cc] & EX_TRAN) fprintf(fp, "TRANSPORT ");
      if (params->prxf[cc] & EX_WAVE) fprintf(fp, "WAVES ");
      if (params->prxf[cc] & EX_TRST) fprintf(fp, "TRACER_STAT ");
      fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
  }

  /*-----------------------------------------------------------------*/
  /* Tidal energy extraction                                         */
  if (master->nturb) {
    FILE *op;
    int cs;
    op = fopen("turb.site", "w");
    fprintf(fp, "Tidal extraction using %d turbines\n", master->nturb);
    for (n = 0; n < master->nturb; n++) {
      fprintf(fp, "Turbine%d @ %d[%f %f], %f m, Cext=%f", n, master->turb[n],
	      params->turbv[0][n], params->turbv[1][n], 
	      params->turbv[2][n], master->cturb[n]);
      /*fprintf(op, "%f %f t%d", params->turbv[0][n], params->turbv[1][n], n);*/
      c = master->turb[n];
      cs = geom->m2d[c];
      fprintf(op, "%f %f t%d", geom->cellx[cs], geom->celly[cs], n);
    }
    fprintf(fp, "\n");
    fclose(op);
  }

  /*-----------------------------------------------------------------*/
  /* Grid size                                                       */
  if (params->us_type & US_IJ) {
    fprintf(fp,"\nNumber of 2D wet cells = %d (%2.0f%% of total grid)\n",
	    geom->b2_t, 100.0*geom->b2_t/(geom->nce1 * geom->nce2));
    fprintf(fp,"Number of 3D wet cells = %d (%2.0f%% of total grid)\n",
	    geom->b3_t, 100.0*geom->b3_t/(geom->nz * geom->nce1 * geom->nce2));
  } else {
    fprintf(fp,"\nNumber of 2D wet cells = %d\n", geom->b2_t);
    fprintf(fp,"Number of 3D wet cells = %d\n", geom->b3_t);
  }
  fprintf(fp,"Percent of cells beneath the sea bed = %d%%\n", 100 * geom->bdry / 
	  (geom->b2_t * geom->nz));
  d1 = d2 = 0.0;
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    d1 += geom->cellarea[c];
    d2 += geom->cellarea[c] * (master->eta[c] - geom->botz[c]);
  }
  fprintf(fp,"Domain surface area = %e m^2\n", d1);
  fprintf(fp,"Domain initial volume = %e m^3\n", d2);
  fprintf(fp,"  (Use -debug init_m on the compas command line to get full stats).\n");

  /* 3D Tracer relaxation */
  if (master->nrlx)
    fprintf(fp, "\nRelaxing %d 3D tracers\n", master->nrlx);
  else
    fprintf(fp, "\nNO 3D tracer relaxtions\n");

  for (n=0; n<master->nrlx; n++) {
    int tr_rlx = master->relax[n];
    fprintf(fp, "  Relax #%2d (%s)\n", n, master->trinfo_3d[tr_rlx].name);
    fprintf(fp, "   File/Tracer = %s\n", master->trinfo_3d[tr_rlx].relax_file);
    fprintf(fp, "   DT     = %s\n", master->trinfo_3d[tr_rlx].relax_dt);
    fprintf(fp, "   TCONST = %s\n", master->trinfo_3d[tr_rlx].r_rate);
  }
  fprintf(fp, "\n");

  /*-----------------------------------------------------------------*/
  /* Diagnostic files                                                */
  fprintf(fp, "Model runtime statistics written to file diag.txt\n");
  sprintf(buf, "%scrash.site", master->opath);
  fprintf(fp, "File containing location of model instabilities = %s\n", buf);
  if (crash_restart)
    fprintf(fp, "CRASH RECOVERY mode invoked: log file = %s\n", params->crashname);
  if (!(params->dbgf & NONE))
    fprintf(fp, "Debugging file written to: debug.txt\n");
  if (!(params->history & NONE))
    fprintf(fp, "History log file written to: %s\n", params->histname);
  if (!(params->history & HST_MASTER))
    fprintf(fp, "Master history log file written to: %s\n", params->histnamem);
  if (params->history & HST_NOTES) {
    char pname[MAXSTRLEN];
    strcpy(buf, params->histname);
    n = strlen(buf);
    for (nn = 0; nn < n-4; nn++) pname[nn] = buf[nn];
    pname[nn] = '\0';
    strcat(pname, "notes");
    fprintf(fp, "Historical notes summary written to: %s\n", pname);
  }
  if (params->meshinfo) {
    mesh_ofile_name(params, buf);
    strcat(buf, "_meshfiles.txt");
    fprintf(fp, "Mesh information files written when using -g option.\n");
    fprintf(fp, "    Mesh summary in file: %s\n", buf);
  }
  fprintf(fp, "\n");

  /******/
  /* DA */
  /******/
#ifdef HAVE_DA
  fprintf(fp, "\n");
  if (master->da == NONE) {
    fprintf(fp,"Data Assimilation is NOT invoked\n");
  } else {
    fprintf(fp,"Data Assimilation is invoked\n");
    dassim_run_setup(fp, sched_get_even_by_name(schedule, "dassim"));
  }
#endif
  fprintf(fp, "\n");

  /* Sediments */
#if defined(HAVE_SEDIMENT_MODULE)
  if (params->do_sed) {
    win_priv_t *wincon = hd_data->wincon[1];
    fprintf(fp, "############################################################################\n");
    fprintf(fp, "# Sediment specification\n\n");
    ntr = tracer_write(params, fp);
    ntr = sediment_autotracer_write(master, fp, ntr);
    trans_write_sed(params, wincon->sediment, fp);
    fprintf(fp, "\n");
  }
#endif

#if defined(HAVE_ECOLOGY_MODULE)
  if (params->do_eco) {
    win_priv_t *wincon = hd_data->wincon[1];
    fprintf(fp, "############################################################################\n");
    fprintf(fp, "# Ecology specification\n\n");
    trans_write_eco(params->eco_vars, params->eco_defs, params->ecodt, wincon->e, fp);
    ntr = ecology_autotracer_write(master, fp, ntr);
    fprintf(fp, "\n");
  }
#endif

  fclose(fp);
}

void wbdrycustom(FILE * fp, parameters_t *params, int bnum,
                 bdry_details_t *data)
{
  if (data->explct) {
    if (data->custom_m == NULL)
      fprintf(fp, "    Boundary fill value = %g\n", data->fill_value);
    else {
      fprintf(fp, "    Custom %s specification : %s\n", data->name,
	      data->custom_tag);
      if (strcmp(data->custom_tag, "use_eqn") == 0) {
	void **ptr = (void **)data->custdata;
	fprintf(fp, "\t %s \n", EqnDisplayStr(ptr[0]));
      } else if (strcmp(data->custom_tag, "u1flowbdry") == 0 ||
		 strcmp(data->custom_tag, "u2flowbdry") == 0) {
	flow_data_t *d = data->custdata;
	fprintf(fp, "\t Flow profile depth = %f \n", d->hc);
	if (d->rlen)
	  fprintf(fp, "\t River length = %f \n", d->rlen);
      } else
	fprintf(fp, "      Forcing file : %s\n", data->args[0]);
    }
  }
}

/* END print_diag()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Writes the geographic locations of the open boundaries to file    */
/*-------------------------------------------------------------------*/
void write_bdry(parameters_t *params,  /* Input parameter data */
		geometry_t *geom)
{
  FILE *fp;
  int cc, c, n, nn, nt, tn, nf, m = 0;
  char key[MAXSTRLEN];
  char suf[MAXSTRLEN];
  char tinc[MAXSTRLEN];
  char type[MAXSTRLEN];
  char fill[MAXSTRLEN];
  char filter[MAXSTRLEN];
  char sync[MAXSTRLEN];
  int mpk = 0;
  int nel = 0, nts = 0, rf = 0;
  int bcond = 0;

  if(strlen(params->bdry_file) == 0) return;

  if(!(fp = fopen(params->bdry_file, "w")))
    hd_warn("Can't open boundary setup file %s\n", params->bdry_file);

  /*
   * Special case for RECOM 
   * 
   * Just write out the cell centres for SWAN forcing
   */
  if (params->roammode & (A_RECOM_R1|A_RECOM_R2)) {
    fprintf(fp, "# Cell centre locations of all OPEN boundary points, excluding RIVERs\n");
    fprintf(fp, "#   X         Y\n");
    for (n = 0; n < geom->nobc; n++) {
      open_bdrys_t *open = geom->open[n];
      /* Ignore RIVER boundaries */
      if (strlen(open->bflow)) continue;
      for (cc = 1; cc <= open->no2_t; cc++) {
	c = open->obc_t[cc];
	fprintf(fp, "%f %f\n", geom->cellx[c], geom->celly[c]);
      }
    }
    fprintf(fp,"\n");
    fclose(fp);
    return;
  }

  strcpy(suf, "nc");
  strcpy(type, "parray");
  strcpy(tinc, "1 hour");
  sprintf(fill, "%c", '\0');
  sprintf(filter, "%c", '\0');
  sprintf(sync, "%c", '\0');

  key[0] = params->bdry_file[0];
  key[1] = params->bdry_file[1];
  key[2] = '\0';
  if (strcmp(key, "lr") == 0) rf = 1;
  if (strcmp(key, "hr") == 0) rf = 2;

  /* Determine standard boundary types */
  /* Flather */
  bcond = NEST_FLA|NEST_BARO;
  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    if (!open->sbcond) continue;
    if (!(open->sbcond & NEST_BARO)) bcond &= ~NEST_BARO;
    if (!(open->sbcond & NEST_FLA)) bcond = 0;
  }
  /* Radiation */
  if (!bcond) {
    bcond = NEST_RAD|NEST_BARO;
    for (n = 0; n < geom->nobc; n++) {
      open_bdrys_t *open = geom->open[n];
      if (!open->sbcond) continue;
      if (!(open->sbcond & NEST_BARO)) bcond &= ~NEST_BARO;
      if (!(open->sbcond & NEST_RAD)) bcond = 0;
      
    }
  }
  /* Clamped */
  if (!bcond) {
    bcond = NEST_CPD|NEST_BARO;
    for (n = 0; n < geom->nobc; n++) {
      open_bdrys_t *open = geom->open[n];
      if (!open->sbcond) continue;
      if (!(open->sbcond & NEST_BARO)) bcond &= ~NEST_BARO;
      if (!(open->sbcond & NEST_CPD)) bcond = 0;
      
    }
  }
  /* 1 or 2 way nesting */
  if (!bcond) {
    bcond = NEST2WAY|NEST1WAY|NEST_BARO;
    for (n = 0; n < geom->nobc; n++) {
      open_bdrys_t *open = geom->open[n];
      if (!open->sbcond) continue;
      if (!(open->sbcond & NEST_BARO)) bcond &= ~NEST_BARO;
      if (!(open->sbcond & (NEST2WAY|NEST1WAY))) bcond = 0;
    }
  }
  if (!bcond) bcond = NEST1WAY;

  if (endswith(params->bdry_file, ".mpk")) {
    mpk = 1;
    for (n = 0; n < strlen(params->bdry_file) - 4; n++)
      key[n] = params->bdry_file[n];
    key[n] = '\0';
    strcpy(suf, "mpk");
    strcpy(type, "memory");
    sprintf(tinc, "%5.1f seconds", params->grid_dt);
    strcpy(fill, "cascade_search");
    strcpy(filter, "copy");
    strcpy(sync, "0 seconds");
  } else if (endswith(params->prmname, ".prm")) {
    for (n = 0; n < strlen(params->prmname) - 4; n++)
      key[n] = params->prmname[n];
    key[n] = '\0';
  } else
    strcpy(key, params->prmname);

  /* Count the number of .ets points */
  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    for (cc = 1; cc <= open->no2_t; cc++) {
      nel++;
      nts++;
      nf = 0;
      for (tn = 0; tn < params->ntr; tn++)
	if (open->bcond_tra[tn] & TRCONC) nf = 1;
      if (nf) {
	for (nn = 0; nn < open->bgz; nn++) {
	  nts++;
	}
      }
    }
  }

  if (mpk) {
    if (bcond & (NEST1WAY|NEST2WAY|NEST_CPD)) {
      if (bcond & NEST_BARO)
	fprintf(fp, "OutputFiles          %d\n\n", geom->nobc * 4 + 2);
      else
	fprintf(fp, "OutputFiles          %d\n\n", geom->nobc * 2 + 1);
    } else if (bcond & NEST_FLA) {
      if (bcond & NEST_BARO)
	fprintf(fp, "OutputFiles          %d\n\n", geom->nobc + 2);
      else
	fprintf(fp, "OutputFiles          %d\n\n", geom->nobc + 1);
    } else if (bcond & NEST_RAD) {
      if (bcond & NEST_BARO)
	fprintf(fp, "OutputFiles          2\n\n");
      else
	fprintf(fp, "OutputFiles          1\n\n");
    }
  }

  /* Sea level in separate file for memory output                    */
  if (mpk && bcond & NEST_BARO) {
    fprintf(fp,"file%d.name           %s_eta.%s\n",m,key,suf);
    fprintf(fp,"file%d.filetype       %s\n", m, type);
    fprintf(fp,"file%d.tstart         %f days\n",m, 
	    schedule->start_time/86400);
    fprintf(fp,"file%d.tinc           %s\n",m, tinc);
    fprintf(fp,"file%d.tstop          %f days\n",m,
	    schedule->stop_time/86400);
    if (rf == 2) {
      fprintf(fp,"#file%d.sync_dt        %s\n",m, sync);
      fprintf(fp,"#file%d.filter         weighted5\n",m);
    } else {
      fprintf(fp,"file%d.sync_dt        %s\n",m, sync);
      fprintf(fp,"file%d.filter         %s\n",m, filter);
    }
    fprintf(fp,"file%d.fill_rule      %s\n",m, fill);
    fprintf(fp,"file%d.bytespervalue  4\n",m);
    fprintf(fp,"file%d.vars           eta\n",m);
    fprintf(fp,"file%d.points         %d\n", m,nel);
    for (n = 0; n < geom->nobc; n++) {
      open_bdrys_t *open = geom->open[n];
      for (cc = 1; cc <= open->no2_t; cc++) {
	c = open->obc_t[cc];
	fprintf(fp, "%f %f\n", geom->cellx[c], geom->celly[c]);
      }
    }
    fprintf(fp,"\n");
    m++;
  }

  /* Temperature and salinity (and sea level)                        */
  if (mpk && bcond & NEST_BARO)
    fprintf(fp,"file%d.name           %s_ts.%s\n",m,key,suf);
  else
    fprintf(fp,"file%d.name           %s_ets.%s\n",m,key,suf);
  fprintf(fp,"file%d.filetype       %s\n", m, type);
  fprintf(fp,"file%d.tstart         %f days\n",m, 
	  schedule->start_time/86400);
  fprintf(fp,"file%d.tinc           %s\n",m, tinc);
  fprintf(fp,"file%d.tstop          %f days\n",m,
	  schedule->stop_time/86400);
  if (mpk) {
    if (rf == 2) {
      fprintf(fp,"#file%d.sync_dt        %s\n",m, sync);
      fprintf(fp,"#file%d.filter         weighted5\n",m);
    } else {
      fprintf(fp,"file%d.sync_dt        %s\n",m, sync);
      fprintf(fp,"file%d.filter         %s\n",m, filter);
    }
    fprintf(fp,"file%d.fill_rule      %s\n",m, fill);
  }
  fprintf(fp,"file%d.bytespervalue  4\n",m);
  if (mpk && bcond & NEST_BARO)
    fprintf(fp,"file%d.vars           temp salt\n",m);
  else
    fprintf(fp,"file%d.vars           eta temp salt\n",m);
  fprintf(fp,"file%d.points         %d\n", m,nts);

  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    for (cc = 1; cc <= open->no2_t; cc++) {
      int e, j;
      c = open->obc_t[cc];
      fprintf(fp, "%f %f\n", geom->cellx[c], geom->celly[c]);
      nf = 0;
      for (tn = 0; tn < params->ntr; tn++)
	if (open->bcond_tra[tn] & TRCONC) nf = 1;
      if (nf) {
	for (j = 1; j <= geom->npe[c]; j++) {
	  if ((e = open->bec[j][c])) {
	    for (nn = 0; nn < open->bgz; nn++) {
	      c = open->omape[e][c];
	      fprintf(fp, "%f %f\n", geom->cellx[c], geom->celly[c]);
	    }
	  }
	}
      }
    }
  }
  fprintf(fp,"\n");
  m++;
  if (bcond & NEST_RAD && bcond & NEST_BARO) bcond &= ~NEST_BARO;

  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    /* .ets dumped for individual boundaries
    fprintf(fp,"file%d.name           %s_ets.nc\n",m,open->name);
    fprintf(fp,"file%d.filetype       parray\n", m);
    fprintf(fp,"file%d.tstart         %f days\n",m, 
	    schedule->start_time/86400);
    fprintf(fp,"file%d.tinc           1 hour\n",m);
    fprintf(fp,"file%d.tstop          %f days\n",m,
	    schedule->stop_time/86400);
    fprintf(fp,"file%d.bytespervalue  4\n",m);
    fprintf(fp,"file%d.vars           eta temp salt\n",m);
    fprintf(fp,"file%d.points         %d\n", m,open->no2_t);
    for (cc = 1; cc <= open->no2_t; cc++) {
      c = open->obc_t[cc];
      fprintf(fp, "%f %f\n", geom->cellx[c], geom->celly[c]);
    }
    fprintf(fp,"\n");
    m++;
    */

    /* 3D normal and tangential velocity                             */
    if (bcond & (NEST1WAY|NEST2WAY|NEST_CPD)) {
      fprintf(fp,"file%d.name           %s_%s_uv_nor.%s\n",m,key,open->name,suf);
      fprintf(fp,"file%d.filetype       %s\n", m, type);
      fprintf(fp,"file%d.tstart         %f days\n",m, 
	      schedule->start_time/86400);
      fprintf(fp,"file%d.tinc           %s\n",m, tinc);
      fprintf(fp,"file%d.tstop          %f days\n",m,
	      schedule->stop_time/86400);
      if (mpk) {
	if (rf == 2) {
	  fprintf(fp,"#file%d.sync_dt        %s\n",m, sync);
	  fprintf(fp,"file%d.filter         average3\n",m);
	} else {
	  fprintf(fp,"file%d.sync_dt        %s\n",m, sync);
	  fprintf(fp,"file%d.filter         %s\n",m, filter);
	}
	fprintf(fp,"file%d.fill_rule      %s\n",m, fill);
      }
      fprintf(fp,"file%d.bytespervalue  4\n",m);
      fprintf(fp,"file%d.vars           u v\n",m);

      fprintf(fp,"file%d.points         %d\n", m,open->no2_e1);
      for (cc = 1; cc <= open->no2_e1; cc++) {
	c = open->obc_e1[cc];
	/* Require generalization */
	if (open->stagger & INFACET) c = open->nmap[c];
	fprintf(fp, "%f %f\n", geom->u1x[c], geom->u1y[c]);
      }
      fprintf(fp,"\n");
      m++;

      nt = open->to2_e1 - open->no3_e1;
      if (nt) {
	fprintf(fp,"file%d.name           %s_%s_uv_tan.%s\n",m,key,open->name,suf);
	fprintf(fp,"file%d.filetype       %s\n", m, type);
	fprintf(fp,"file%d.tstart         %f days\n",m, 
		schedule->start_time/86400);
	fprintf(fp,"file%d.tinc           %s\n",m, tinc);
	fprintf(fp,"file%d.tstop          %f days\n",m,
		schedule->stop_time/86400);
	if (mpk) {
	  if (rf == 2) {
	    fprintf(fp,"#file%d.sync_dt        %s\n",m, sync);
	    fprintf(fp,"file%d.filter         weighted5\n",m);
	  } else {
	    fprintf(fp,"file%d.sync_dt        %s\n",m, sync);
	    fprintf(fp,"file%d.filter         %s\n",m, filter);
	  }
	  fprintf(fp,"file%d.fill_rule      %s\n",m, fill);
	}
	fprintf(fp,"file%d.bytespervalue  4\n",m);
	fprintf(fp,"file%d.vars           u v\n",m);

	fprintf(fp,"file%d.points         %d\n", m,nt);
	for (cc = open->no3_e1 + 1; cc <= open->to2_e1; cc++) {
	  c = open->obc_e1[cc];
	  fprintf(fp, "%f %f\n", geom->u1x[c], geom->u1y[c]);
	}
	fprintf(fp,"\n");
	m++;
      }
    }

    /* 2D normal and tangential velocity                             */
    if (mpk && bcond & (NEST_FLA|NEST_BARO)) {
      fprintf(fp,"file%d.name           %s_%s_uvav_nor.%s\n",m,key,open->name,suf);
      fprintf(fp,"file%d.filetype       %s\n", m, type);
      fprintf(fp,"file%d.tstart         %f days\n",m, 
	      schedule->start_time/86400);
      fprintf(fp,"file%d.tinc           %s\n",m, tinc);
      fprintf(fp,"file%d.tstop          %f days\n",m,
	      schedule->stop_time/86400);
      if (rf == 2) {
	fprintf(fp,"#file%d.sync_dt        %s\n",m, sync);
	fprintf(fp,"file%d.filter         average3\n",m);
      } else {
	fprintf(fp,"file%d.sync_dt        %s\n",m, sync);
	fprintf(fp,"file%d.filter         %s\n",m, filter);
      }
      fprintf(fp,"file%d.fill_rule      %s\n",m, fill);
      fprintf(fp,"file%d.bytespervalue  4\n",m);
      fprintf(fp,"file%d.vars           uav vav\n",m);

      fprintf(fp,"file%d.points         %d\n", m,open->no2_e1);
      for (cc = 1; cc <= open->no2_e1; cc++) {
	c = open->obc_e1[cc];
	/* Require generalization */
	if (open->stagger & INFACET) c = open->nmap[c];
	fprintf(fp, "%f %f\n", geom->u1x[c], geom->u1y[c]);
      }
      fprintf(fp,"\n");
      m++;
      nt = open->to2_e1 - open->no3_e1;
      if (nt && bcond & NEST_BARO) {
	fprintf(fp,"file%d.name           %s_%s_uvav_tan.%s\n",m,key,open->name,suf);
	fprintf(fp,"file%d.filetype       %s\n", m, type);
	fprintf(fp,"file%d.tstart         %f days\n",m, 
		schedule->start_time/86400);
	fprintf(fp,"file%d.tinc           %s\n",m, tinc);
	fprintf(fp,"file%d.tstop          %f days\n",m,
		schedule->stop_time/86400);
	if (rf == 2) {
	  fprintf(fp,"#file%d.sync_dt        %s\n",m, sync);
	  fprintf(fp,"file%d.filter         weighted5\n",m);
	} else {
	  fprintf(fp,"file%d.sync_dt        %s\n",m, sync);
	  fprintf(fp,"file%d.filter         %s\n",m, filter);
	}
	fprintf(fp,"file%d.fill_rule      %s\n",m, fill);
	fprintf(fp,"file%d.bytespervalue  4\n",m);
	fprintf(fp,"file%d.vars           uav vav\n",m);

	fprintf(fp,"file%d.points         %d\n", m,nt);
	for (cc = open->no3_e1 + 1; cc <= open->to2_e1; cc++) {
	  c = open->obc_e1[cc];
	  fprintf(fp, "%f %f\n", geom->u1x[c], geom->u1y[c]);
	}
	fprintf(fp,"\n");
	m++;
      }
    }
  }
  fclose(fp);
}

/* END write_bdry()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Writes the sparse - Cartesian window map to file                  */
/*-------------------------------------------------------------------*/
void write_window_map(geometry_t **windows, /* Window geometery      */
		      parameters_t *params  /* Input parameters data */
		      )
{
  FILE *fp = NULL;       /* File handles for printing data           */
  int n;                 /* Window counters                          */
  int c, cc, ck, c1;     /* Sparse location counters                 */
  int i, j, k;           /* (i,j,k) locations of sparse coordinates  */
  char buf[10];          /* Buffer for cell status                   */
  int nwindows;          /* Number of windows to make                */
  long t;                /* Time variable                            */

  if (windows_log) {
    fp = fopen(window_map_logfile, "w");
    if (fp == NULL)
      hd_quit("window_build: Can't write to file %s\n",
              window_map_logfile);

    nwindows = params->nwindows;
    fprintf(fp, "\nInput file = %s\n\n", params->prmname);
    time(&t);
    fprintf(fp, "Number of windows = %d\n", nwindows);
    fprintf(fp, "Written at :  %s\n", ctime(&t));

    for (n = 1; n <= nwindows; n++) {
      fprintf(fp, "\nWINDOW #%d sparse locations\n\n", n);
      fprintf(fp, "2D cells to process vectors\n");
      fprintf(fp, "   Local     w2_t      w2_e1    w2_e2\n\n");
      ck = max(window[n]->n2_t, window[n]->n2_e1);
      ck = max(ck, window[n]->n2_e2);
      for (c = 1; c <= ck; c++) {
        fprintf(fp, "  %4d    ", c);
        if (c <= window[n]->v2_t)
          fprintf(fp, " %4d     ", window[n]->w2_t[c]);
        else if (c > window[n]->v2_t && c <= window[n]->b2_t)
          fprintf(fp, " %4d (b) ", window[n]->w2_t[c]);
        else if (c > window[n]->b2_t && c <= window[n]->a2_t)
          fprintf(fp, " %4d (a) ", window[n]->w2_t[c]);
        else if (c > window[n]->a2_t && c <= window[n]->n2_t)
          fprintf(fp, " %4d (g) ", window[n]->w2_t[c]);
        else
          fprintf(fp, "          ");
        if (c <= window[n]->v2_e1)
          fprintf(fp, " %4d     ", window[n]->w2_e1[c]);
        else if (c > window[n]->v2_e1 && c <= window[n]->b2_e1)
          fprintf(fp, " %4d (b) ", window[n]->w2_e1[c]);
        else if (c > window[n]->b2_e1 && c <= window[n]->a2_e1)
          fprintf(fp, " %4d (a) ", window[n]->w2_e1[c]);
        else if (c > window[n]->a2_e1 && c <= window[n]->x2_e1)
          fprintf(fp, " %4d (x) ", window[n]->w2_e1[c]);
        else if (c > window[n]->x2_e1 && c <= window[n]->n2_e1)
          fprintf(fp, " %4d (g) ", window[n]->w2_e1[c]);
        else
          fprintf(fp, "          ");
        if (c <= window[n]->v2_e2)
          fprintf(fp, " %4d     \n", window[n]->w2_e2[c]);
        else if (c > window[n]->v2_e2 && c <= window[n]->b2_e2)
          fprintf(fp, " %4d (b) \n", window[n]->w2_e2[c]);
        else if (c > window[n]->b2_e2 && c <= window[n]->a2_e2)
          fprintf(fp, " %4d (a) \n", window[n]->w2_e2[c]);
        else if (c > window[n]->a2_e2 && c <= window[n]->x2_e2)
          fprintf(fp, " %4d (x) \n", window[n]->w2_e2[c]);
        else if (c > window[n]->x2_e2 && c <= window[n]->n2_e2)
          fprintf(fp, " %4d (g) \n", window[n]->w2_e2[c]);
        else
          fprintf(fp, "          \n");
      }

      fprintf(fp, "\n2D surface and bottom vectors\n");
      fprintf(fp,
              "   Local     sur_t  bot_t    sur_e1  bot_e1   sur_e2  bot_e2\n\n");
      ck = max(window[n]->v2_t, window[n]->v2_e1);
      ck = max(ck, window[n]->v2_e2);
      for (c = 1; c <= ck; c++) {
        fprintf(fp, "  %4d    ", c);
        if (c <= window[n]->v2_t)
          fprintf(fp, " %4d     %4d  ", window[n]->nsur_t[c],
                  window[n]->bot_t[c]);
        else
          fprintf(fp, "                ");
        if (c <= window[n]->v2_e1)
          fprintf(fp, " %4d      %4d   ", window[n]->sur_e1[c],
                  window[n]->bot_e1[c]);
        else
          fprintf(fp, "                  ");
        if (c <= window[n]->v2_e2) {
          fprintf(fp, " %4d      %4d \n", window[n]->sur_e2[c],
                  window[n]->bot_e2[c]);
        } else
          fprintf(fp, "                \n");
      }

      fprintf(fp, "\nSparse coordinates\n");
      fprintf(fp,
              "   Local    Global       (i,j,k)     status  window\n\n");
      for (c = 1; c <= window[n]->enon; c++) {
        cc = window[n]->wsa[c];
        i = geom->s2i[cc];
        j = geom->s2j[cc];
        k = geom->s2k[cc];
        c1 = window[n]->wsa[window[n]->zp1[c]];
        strcpy(buf, "WET");
        if (geom->fm[cc].wn != n)
          strcpy(buf, "AUX");
        if (geom->fm[c1].wn == n && c == window[n]->zm1[c])
          strcpy(buf, "SED");
        if (i == NOTVALID && j == NOTVALID && k == NOTVALID)
          strcpy(buf, "GHOST");
        fprintf(fp, "  %4d      %5d     %3d,%3d,%3d   %6s   %3d\n", c, cc,
                i, j, k, buf, n);
      }
    }
    fclose(fp);
  }
}
/* END write_window_map()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Writes a transport file from an input file configuration          */
/*-------------------------------------------------------------------*/
void trans_write(hd_data_t *hd_data)
{
  parameters_t *params = hd_data->params;
  dump_data_t *dumpdata = hd_data->dumpdata;
  master_t *master = hd_data->master;
  geometry_t *geom = master->geom;
  FILE *op, *fopen();
  int n, tn, i, j;
  double d1, dt = HUGE;
  char tag[MAXSTRLEN];
  char key[MAXSTRLEN];
  char bname[MAXSTRLEN];
  char trdata[MAXSTRLEN];
  char TRANS_DT[MAXSTRLEN];
  char fvars[MAXSTRLEN];

  if (params->runmode & TRANS) return;
  if(!strlen(params->trans_data))
    if (!(params->tmode & SP_DUMP))
      return;

  strcpy(trdata, params->trans_data);

  if (strlen(params->trans_dt))
    sprintf(TRANS_DT, "%s", params->trans_dt);
  else
    sprintf(TRANS_DT, "1 hour");
  if (params->tmode & SP_DUMP) {
    if (params->wind) dt = min(dt, params->wind_dt);
    if (params->patm) dt = min(dt, params->patm_dt);
    for (i = 0; i < geom->nobc; i++) {
      open_bdrys_t *open = geom->open[i];
      if (open->file_dt) dt = min(dt, open->file_dt);
    }
    strcpy(TRANS_DT, otime(dt, key));
  }

  strcpy(tag, params->idumpname);
  stripend(tag);
  strcat(tag,".tran");

  if ((op = fopen(tag, "w")) == NULL) {
    hd_warn("trans_write: Can't open file %s\n", tag);
    return;
  }

  fprintf(op, "# COMPAS transport file\n");
  fprintf(op, "# Created from %s\n", params->prmname);
  fprintf(op, "# Use with INPUT_FILE %s\n", params->idumpname);
  fprintf(op, "# EMS Version : %s\n", version);
  fprintf(op, "CODEHEADER           %s\n", params->codeheader);
  fprintf(op, "PARAMETERHEADER      %s\n", params->parameterheader);
  fprintf(op, "DESCRIPTION          %s\n", params->grid_desc);
  fprintf(op, "NAME                 %s\n", params->grid_name);
  fprintf(op, "TIMEUNIT             %s\n", params->timeunit);
  fprintf(op, "OUTPUT_TIMEUNIT      %s\n", params->output_tunit);
  fprintf(op, "LENUNIT              %s\n", params->lenunit);
  if (strlen(params->projection))
    fprintf(op, "PROJECTION           %s\n", params->projection);
  fprintf(op, "START_TIME           %s\n", params->start_time);
  fprintf(op, "STOP_TIME            %s\n\n", params->stop_time);

  strcpy(tag, params->idumpname);
  stripend(tag);
  strcat(tag,".nc");
  fprintf(op, "INPUT_FILE           %s\n\n", tag);

  /* Force global for PRE_MARVL */
  if (params->runmode & PRE_MARVL) {
    /* Velocity data and eta should've been specified in the auto file */
    fprintf(op, "TRANS_DATA  %s %s\n", params->vdata, params->edata);
    fprintf(op, "TRANS_MODE  GLOBAL\n");
    fprintf(op, "FILL_METHOD NO_FILL\n");
    /* Z to sigma levels scaling factor */
    fprintf(op, "\nROMS_Z2S_FACTOR  %.2f\n", params->roms_z2s);
  } else if (params->tmode & SP_DUMP) {
    fprintf(op, "TRANS_MODE           SP_DUMP\n");
  } else {
    if(params->tmode & SP_FFSL)
      fprintf(op, "TRANS_DATA            %s(u1=umean)(w=wmean)(Kz=Kzmean)(u1vm=u1vmean)\n", trdata);
    else
      fprintf(op, "TRANS_DATA           %s(u1=umean)(w=wmean)(Kz=Kzmean)\n", trdata);
    /*
    if(!(params->tmode & SP_FFSL) && strlen(params->sourcefile))
      fprintf(op, "SOURCE_GRID          %s\n", params->sourcefile);
    */
    if(strlen(params->trvars))
      fprintf(op, "TRANS_VARS           %s\n", params->trvars);

    if(params->tmode == SP_EXACT)
      fprintf(op, "TRANS_MODE           SP_EXACT\n");
    else if(params->tmode & SP_TINT)
      fprintf(op, "TRANS_MODE           SP_INTERP\n");
    else if(params->tmode & XYZ_TINT)
      fprintf(op, "TRANS_MODE           XYZ_INTERP\n");
    else if(params->tmode & SP_FFSL) {
      fprintf(op, "TRANS_MODE            SP_FFSL\n");
      fprintf(op, "FILL_METHOD           NONE\n");
      fprintf(op, "TRA_SCHEME            FFSL\n");
      fprintf(op, "STABILITY             SUB-STEP-NOSURF\n");
      fprintf(op, "MERGE_THIN            YES\n");
      fprintf(op, "CONSERVATION          ETA W\n");
      if (params->roammode & (A_RECOM_R1|A_RECOM_R2))
	fprintf(op, "TRANS_OMP_NUM_THREADS %d\n", params->nwindows);
      else 
	fprintf(op, "TRANS_OMP_NUM_THREADS 1\n");
    } else
      fprintf(op, "TRANS_MODE           SP_EXACT\n");
  }

  /*
   * FR : Commenting out the following code as it was causing
   *      overmixing in the Fitzroy model so assuming it will also be
   *      an issue in RECOM
   if(params->tmode & SP_FFSL) {
   if ((params->u1kh < 0.0 || params->u2kh < 0.0) && params->smagorinsky != 0.0) {
   fprintf(op, "\n# Horizontal mixing\n");
   fprintf(op, "DIFF_SCALE           LINEAR\n");
   fprintf(op, "U1KH                 %5.3f\n", params->u1kh);
   fprintf(op, "U2KH                 %5.3f\n", params->u2kh);
   fprintf(op, "SMAGORINSKY          %5.3f\n\n", params->smagorinsky);
   }
   }
  */

  if (!(params->tmode & SP_DUMP)) {
    if(params->fillf & GLOBAL)
      fprintf(op, "FILL_METHOD          GLOBAL");
    else if(params->fillf & MONOTONIC)
      fprintf(op, "FILL_METHOD          MONOTONIC");
    if(params->fillf & OBC_ADJUST)
      fprintf(op, " OBC_ADJUST");
    if(params->fillf & DIAGNOSE)
      fprintf(op, " DIAGNOSE");
    if(params->fillf & DIAGNOSE_BGC)
      fprintf(op, " DIAGNOSE_BGC");
    fprintf(op, "\n\n");
  }

  if (params->runmode & PRE_MARVL)
    write_grid_specs(op, params);

  fprintf(op, "DT                   %s\n", TRANS_DT);
  if (!(params->tmode & SP_DUMP)) {
    /* Turn on CFL diagnostics for RECOM */
    if (params->roammode & (A_RECOM_R1|A_RECOM_R2))
      fprintf(op, "CFL                  PASSIVE\n");
    fprintf(op, "HMIN                 %-6.4f\n", params->hmin);
    fprintf(op, "Z0                   %-6.4f\n", params->z0);
    fprintf(op, "NUMBERS              RESOLUTION CELL_INDEX\n");
  }
  fprintf(op, "\n");

  /* Add wind and pressure for PRE_MARVL and WIND for RECOM */
  if (params->runmode & PRE_MARVL || 
      params->tmode & SP_DUMP ||
      params->roammode & (A_RECOM_R1|A_RECOM_R2)) {
    if (strlen(params->wind)) {
      fprintf(op, "WIND_TS              %s\n", params->wind);
      if (params->roammode & (A_RECOM_R1|A_RECOM_R2)) {
	// No faster than the transport dt
	fprintf(op, "WIND_INPUT_DT        %s\n", TRANS_DT);
      } else
	fprintf(op, "WIND_INPUT_DT        %s\n", otime(params->wind_dt, tag));
      fprintf(op, "WIND_SPEED_SCALE     %-6.1f\n", params->wind_scale);
      if (params->wind_type == SPEED) {
	fprintf(op, "DRAG_LAW_V0          %-6.1f\n", params->dlv0);
	fprintf(op, "DRAG_LAW_V1          %-6.1f\n", params->dlv1);
	fprintf(op, "DRAG_LAW_CD0         %-8.5f\n", params->dlc0);
	fprintf(op, "DRAG_LAW_CD1         %-8.5f\n", params->dlc1);
      } else
	fprintf(op, "WIND_TYPE            STRESS\n");
      fprintf(op, "\n");
      strcat(fvars, " wind1");
    }

    if (params->tmode & SP_DUMP || params->runmode & PRE_MARVL) {
      if (strlen(params->hf)) {
	fprintf(op, "HEATFLUX             %s\n", heatfluxname(params->heatflux));
	fprintf(op, "HEATFLUX_FILE        %s\n", params->hf);
	fprintf(op, "HEATFLUX_DT          %s\n", otime(params->hf_dt, tag));
	fprintf(op, "\n");
      }
      if (strlen(params->precip)) {
	fprintf(op, "PRECIPITATION          %s\n", params->precip);
	fprintf(op, "PRECIPITATION_INPUT_DT %s\n", otime(params->precip_dt, tag));
	fprintf(op, "\n");
      strcat(fvars, " precipitation");
      }
      if (strlen(params->evap)) {
	fprintf(op, "EVAPORATION          %s\n", params->evap);
	fprintf(op, "EVAPORATION_INPUT_DT %s\n", otime(params->evap_dt, tag));
	fprintf(op, "\n");
      }
      if (strlen(params->patm)) {
	fprintf(op, "PRESSURE             %s\n", params->patm);
	fprintf(op, "PRESSURE_INPUT_DT    %s\n", otime(params->patm_dt, tag));
	fprintf(op, "\n");
	strcat(fvars, " patm");
      }
      if (strlen(params->airtemp)) {
	fprintf(op, "AIRTEMP              %s\n", params->airtemp);
	fprintf(op, "AIRTEMP_INPUT_DT     %s\n", otime(params->airtemp_dt, tag));
	fprintf(op, "\n");
      strcat(fvars, " air_temp");
      }
      if (strlen(params->wetb)) {
	if (master->sh_f & WETBULB) {
	  fprintf(op, "WET_BULB           %s\n", params->wetb);
	  fprintf(op, "WET_BULB_INPUT_DT  %s\n", otime(params->wetb_dt, tag));
	  fprintf(op, "\n");
	  strcat(fvars, " wet_bulb");
	}
	if (master->sh_f & DEWPOINT) {
	  fprintf(op, "DEW_POINT          %s\n", params->wetb);
	  fprintf(op, "DEW_POINT_INPUT_DT %s\n", otime(params->wetb_dt, tag));
	  fprintf(op, "\n");
	  strcat(fvars, " dew_point");
	}
      }
      if (strlen(params->cloud)) {
	fprintf(op, "CLOUD                %s\n", params->cloud);
	fprintf(op, "CLOUD_INPUT_DT       %s\n", otime(params->cloud_dt, tag));
	fprintf(op, "\n");
	strcat(fvars, " cloud");
      }
    }
  }

  if (params->tmode & SP_DUMP) {
    if (strlen(params->opath))
      fprintf(op,"OutputPath             %s\n", params->opath);
    n = 0;
    /* Count the non-river boundaries                                */
    for (i = 0; i < geom->nobc; i++) {
      open_bdrys_t *open = geom->open[i];
      if (strcmp(open->cusname_u1, "u1flowbdry") != 0) n++;
    }
    fprintf(op,"OutputFiles          %d\n\n", n + 1);
    fprintf(op,"file%1.1d.name           force.nc\n",n);
    fprintf(op,"file%1.1d.filetype       ugrid\n",n);
    fprintf(op,"file%1.1d.tinc           %s\n",n,otime(params->wind_dt, tag));
    fprintf(op,"file%1.1d.bytespervalue  4\n",n);
    fprintf(op,"file%1.1d.vars          %s\n",n++, fvars);
    fprintf(op,"\n");
    n = 0;
    for (i = 0; i < geom->nobc; i++) {
      open_bdrys_t *open = geom->open[i];
      if (strcmp(open->cusname_u1, "u1flowbdry") == 0) continue;
      fprintf(op,"file%1.1d.name           bdry_%s.nc\n",n,open->name);
      fprintf(op,"file%1.1d.filetype       sparse\n",n);
      fprintf(op,"file%1.1d.tinc           %s\n",n,TRANS_DT);
      fprintf(op,"file%1.1d.bytespervalue  4\n",n);
      fprintf(op,"file%1.1d.vars           eta temp salt u1\n",n);
      fprintf(op,"file%1.1d.points         %s\n",n++, open->name);
      fprintf(op,"\n");
    }
  } else if(params->ndf && params->runmode & (AUTO|DUMP)) {
    if (strlen(params->opath))
      fprintf(op,"OutputPath             %s\n", params->opath);
    fprintf(op,"# Output files\n");
    if (params->runmode & PRE_MARVL) {
      /* Special handling for PRE_MARVL */
      int actual_ndf, fcnt = 0;
      /* 
       * PRE_MARVL adds dummy file types for MOM and ROMS - remove these 
       */
      actual_ndf = params->ndf;
      for(n = 0; n < params->ndf; n++)
	if (startswith(params->d_filetype[n], "mom_") ||
	    startswith(params->d_filetype[n], "roms_"))
	  actual_ndf--;
      
      fprintf(op,"OutputFiles            %d\n\n", actual_ndf);
      for(n = 0; n < params->ndf; n++) {
	sprintf(tag, "%s", params->d_name[n]);
	if (startswith(params->d_filetype[n], "mom_") ||
	    startswith(params->d_filetype[n], "roms_"))
	  continue;
	fprintf(op,"file%1.1d.name           %s\n",fcnt,tag);
	fprintf(op,"file%1.1d.filetype       %s\n",fcnt,params->d_filetype[n]);
	if (strlen(params->d_tstart[n]))
	  fprintf(op,"file%1.1d.tstart         %s\n",fcnt,params->d_tstart[n]);
	fprintf(op,"file%1.1d.tinc           %s\n",fcnt,params->d_tinc[n]);
	if (strlen(params->d_tstop[n]))
	  fprintf(op,"file%1.1d.tstop          %s\n",fcnt,params->d_tstop[n]);
	fprintf(op,"file%1.1d.bytespervalue  4\n", fcnt);
	fprintf(op,"file%1.1d.vars           %s\n",fcnt,params->d_vars[n]);
	fprintf(op,"\n");
	fcnt++;
      }
    } else {
      fprintf(op,"OutputFiles            %d\n\n", params->ndf);
      for(n = 0; n < params->ndf; n++) {
	sprintf(tag, "tran_%s", params->d_name[n]);
	fprintf(op,"file%1.1d.name           %s\n",n,tag);
	fprintf(op,"file%1.1d.filetype       %s\n",n,params->d_filetype[n]);
	if (strlen(params->d_tstart[n]))
	  fprintf(op,"file%1.1d.tstart         %s\n",n,params->d_tstart[n]);
	fprintf(op,"file%1.1d.tinc           %s\n",n,params->d_tinc[n]);
	if (strlen(params->d_tstop[n]))
	  fprintf(op,"file%1.1d.tstop          %s\n",n,params->d_tstop[n]);
	fprintf(op,"file%1.1d.bytespervalue  4\n",n);
	fprintf(op,"file%1.1d.vars           %s\n",n,params->d_vars[n]);
	fprintf(op,"\n");
      }
    }
  } else {
    fprintf(op, "# Output files\n");
    fprintf(op, "OutputFiles 1\n\n");
    fprintf(op, "file0.name           tran_out.nc\n");
    fprintf(op, "file0.filetype       ugrid\n");
    fprintf(op, "file0.tstart         %s\n", params->start_time);
    fprintf(op, "file0.tinc           1 day\n");
    fprintf(op, "file0.tstop          %s\n", params->stop_time);
    fprintf(op, "file0.bytespervalue  4\n");
    fprintf(op, "file0.vars           ALL\n\n");
  }

  tracer_write(params, op);

  if (params->roammode & (A_ROAM_R4|A_ROAM_R5)) {
    fprintf(op, "\n# Point source specification for auto tracer 'passive'\n");
    fprintf(op, "#pss   %-10.4f %-10.4f -1 0\n",
	    geom->cellx[geom->n2_t / 2], geom->celly[geom->n2_t / 2]);

    fprintf(op, "\n# Particle specification\n");
    fprintf(op, "#particles   %-10.4f %-10.4f -1 0\n",
	    geom->cellx[geom->n2_t / 2], geom->celly[geom->n2_t / 2]);
    fprintf(op, "#PT_InputFile   pt.nc\n\n");
  }

  if (!(params->tmode & SP_DUMP)) {
    fprintf(op, "# Time series\n");
    prm_set_errfn(hd_silent_warn);
    if (prm_read_int(params->prmfd, "TSPOINTS", &tn)) {
      fprintf(op, "TSPOINTS             %d\n\n", tn);
      for (n = 0; n < tn; n++) {
	sprintf(key, "TS%1d.name", n);
	prm_read_char(params->prmfd, key, tag);
	fprintf(op, "TS%1d.name              %s\n", n, tag);
	sprintf(key, "TS%1d.location", n);
	prm_read_char(params->prmfd, key, tag);
	fprintf(op, "TS%1d.location          %s\n", n, tag);
	sprintf(key, "TS%1d.dt", n);
	prm_read_char(params->prmfd, key, tag);
	fprintf(op, "TS%1d.dt                %s\n", n, tag);
	sprintf(key, "TS%1d.reference", n);
	if (prm_read_char(params->prmfd, key, tag))
	  fprintf(op, "TS%1d.reference         %s\n", n, tag);
	else
	  fprintf(op, "TS%1d.reference         msl\n", n);
	fprintf(op, "\n");
      }
    } else {
      fprintf(op, "TSPOINTS             1\n\n");
      fprintf(op, "TS0.name             loc1.ts\n");
      if (params->us_type & US_IJ) {
	fprintf(op, "TS0.location         %-10.4f %-10.4f 0\n",
		dumpdata->cellx[params->nce2 / 2][params->nce1 / 2],
		dumpdata->celly[params->nce2 / 2][params->nce1 / 2]);
      } else {
	fprintf(op, "TS0.location         %-10.4f %-10.4f 0\n",
		geom->cellx[geom->n2_t / 2], geom->celly[geom->n2_t / 2]);
      }
      fprintf(op, "TS0.dt               1 hour\n");
      fprintf(op, "TS0.reference        msl\n\n");
    }
  }


  if (params->roammode & (A_RECOM_R1|A_RECOM_R2)) {
    /* Hard-code specific SWAN configuration */
    fprintf(op, "WAVE_VARS              swan/out_swan.nc\n");
    fprintf(op, "WAVE_VARS_INPUT_DT     1 hour\n");
    fprintf(op, "WAVE_VARS_INTERP_TYPE  nn_non_sibson\n");
  } else
    if (strlen(params->webf)) {
      fprintf(op, "# Activate Waves\n");
      if (params->waves & (WAVE_FOR|TAN_RAD|BOT_STR|VERTMIX)) {
	fprintf(op, "DO_WAVES                 YES\n"); 
	fprintf(op, "WAVES_DT                 1 hour\n\n");
      }
      fprintf(op, "WAVE_VARS                %s\n", params->webf); 
      fprintf(op, "WAVE_VARS_INPUT_DT       1 hour\n");
      if (strlen(params->webf_interp))
	fprintf(op, "WAVE_VARS_INTERP_TYPE    %s\n", params->webf_interp);
      fprintf(op, "\n");
    }

#if defined(HAVE_SEDIMENT_MODULE)
  if (params->do_sed) {
    win_priv_t *wincon = hd_data->wincon[1];
    trans_write_sed(params, wincon->sediment, op);
  }
#endif

#if defined(HAVE_ECOLOGY_MODULE)
  if (params->do_eco) {
    win_priv_t *wincon = hd_data->wincon[1];
    /* Set up multithreading for RECOM */
    if (params->roammode & (A_RECOM_R1|A_RECOM_R2))
      eco_set_omp_num_threads(wincon->e, params->nwindows);
    trans_write_eco(params->eco_vars, params->eco_defs, params->ecodt,wincon->e, op);
    /* RECOM */
    if (params->roammode & (A_RECOM_R1|A_RECOM_R2)) {
      fprintf(op, "# Load LIGHT field from TRANSPORT file\n");
      fprintf(op, "LIGHT  file\n");
    }
  }
#endif

  fprintf(op, "\n# Open boundaries\n");
  fprintf(op, "NBOUNDARIES           %d\n\n", params->nobc);

  for (n = 0; n < params->nobc; n++) {
    open_bdrys_t *open = params->open[n];
    fprintf(op, "BOUNDARY%1.1d.NAME          %s\n", n, open->name);
    fprintf(op, "BOUNDARY%1.1d.TYPE          %s\n", n, btname(open->type));
    /* Forcing dump transport                                        */
    if (params->tmode & SP_DUMP) {
      if (strcmp(open->cusname_u1, "u1flowbdry") != 0) {
	fprintf(op, "BOUNDARY%1.1d.BCOND_ELE     FILEIN\n", n);
	if (strlen(open->cusname_u1)) {
	  fprintf(op, "BOUNDARY%1.1d.BCOND_NOR     CUSTOM\n", n);
	  fprintf(op, "BOUNDARY%1.1d.CUSTOM.u1     %s\n", n, open->cusname_u1);
	} else
	  fprintf(op, "BOUNDARY%1.1d.BCOND_NOR     NOTHIN\n", n);
	if (strlen(open->cusname_u2)) {
	  fprintf(op, "BOUNDARY%1.1d.BCOND_TAN     CUSTOM\n", n);
	  fprintf(op, "BOUNDARY%1.1d.CUSTOM.u2     %s\n", n, open->cusname_u2);
	} else
	  fprintf(op, "BOUNDARY%1.1d.BCOND_TAN     NOTHIN\n", n);
	fprintf(op, "BOUNDARY%1.1d.OPTIONS       OBCWRITE\n", n);
	/*
	  if (open->file_dt)
	  fprintf(op, "BOUNDARY%1.1d.FILEIN_DT     %s\n", n, otime(open->file_dt, bname));
	*/
	for (tn = 0; tn < params->ntr; tn++) {
	  bcname(open->bcond_tra[tn], bname);
	  /* Specify using tracer name */
	  if (!(open->bcond_tra[tn] & NOGRAD))
	    fprintf(op, "BOUNDARY%1.1d.BCOND_%s    %s\n", n, 
		    master->trinfo_3d[tn].name, bname);
	  if (open->bcond_tra[tn] & CUSTOM)
	    fprintf(op, "BOUNDARY%1.1d.CUSTOM.%s   %s\n", n, master->trinfo_3d[tn].name, open->cusname_t[tn]);
	}
	if (strlen(open->tsfn))
	  fprintf(op, "BOUNDARY%1.1d.DATA          %s\n", n, open->tsfn);

      } else {
	fprintf(op, "BOUNDARY%1.1d.BCOND_ELE     NOTHIN\n", n);
	fprintf(op, "BOUNDARY%1.1d.BCOND_NOR     NOTHIN\n", n);
	fprintf(op, "BOUNDARY%1.1d.BCOND_TAN     NOTHIN\n", n);
	fprintf(op, "BOUNDARY%1.1d.BCOND_TRA_ALL NOGRAD\n", n);
      }
    } else {
      /* Other transport files                                       */
      fprintf(op, "BOUNDARY%1.1d.BCOND_ELE     NOTHIN\n", n);
      fprintf(op, "BOUNDARY%1.1d.BCOND_NOR     NOTHIN\n", n);
      fprintf(op, "BOUNDARY%1.1d.BCOND_TAN     NOTHIN\n", n);
      fprintf(op, "BOUNDARY%1.1d.BCOND_TRA_ALL NOGRAD\n", n);
      for (tn = 0; tn < params->ntr; tn++) {
	bcname(open->bcond_tra[tn], bname);
	/* Leave salt/temp out for RECOM */
	if (params->roammode >= A_RECOM_R1|A_RECOM_R2)
	  if (strcmp(params->trinfo_3d->name, "salt") == 0 ||
	      strcmp(params->trinfo_3d->name, "temp") == 0 )
	    continue;
	/* Specify using tracer name */
	if (!(open->bcond_tra[tn] & NOGRAD))
	  fprintf(op, "BOUNDARY%1.1d.BCOND_%s    %s\n", n, 
		  master->trinfo_3d[tn].name, bname);
	if (open->bcond_tra[tn] & CUSTOM)
	  fprintf(op, "BOUNDARY%1.1d.CUSTOM.%s   %s\n", n, master->trinfo_3d[tn].name, open->cusname_t[tn]);
      }
      bcname(open->bcond_Vz, bname);
      fprintf(op, "BOUNDARY%1.1d.BCOND_VZ      %s\n", n, bname);
      bcname(open->bcond_Kz, bname);
      fprintf(op, "BOUNDARY%1.1d.BCOND_KZ      %s\n", n, bname);
      if (params->roammode & (A_RECOM_R1|A_RECOM_R2)) {
	/* Override RIVER boundary forcing */
	if (strlen(open->bflow) && strlen(params->rivldir))
	  fprintf(op, "BOUNDARY%1.1d.DATA          %s/%s.ts\n", n, params->rivldir, open->name);
	else if (strlen(open->tsfn))
	  fprintf(op, "BOUNDARY%1.1d.DATA          %s\n", n, open->tsfn);
      } else
	if (strlen(open->tsfn))
	  fprintf(op, "BOUNDARY%1.1d.DATA          %s\n", n, open->tsfn);
    }
    /* Boundary locations                                            */
    if (params->us_type & US_IJ) {
      fprintf(op, "BOUNDARY%1.1d.POINTS        %d\n", n, open->npts);
      for (tn = 0; tn < open->npts; tn++)
	fprintf(op, "%d %d\n", open->iloc[tn], open->jloc[tn]);
    } else {
      mesh_t *mesh = params->mesh;
      int npts = mesh->npts[n];
      if (npts == 2 && mesh->loc[n][1] == mesh->loc[n][2]) npts = 1;
      fprintf(op, "BOUNDARY%d.UPOINTS     %d\n", n, npts);
      for (tn = 1; tn <= npts; tn++) {
	fprintf(op, "%d (%lf,%lf)-(%lf,%lf)\n", mesh->loc[n][tn], 
		mesh->xloc[mesh->obc[n][tn][0]], mesh->yloc[mesh->obc[n][tn][0]],
		mesh->xloc[mesh->obc[n][tn][1]], mesh->yloc[mesh->obc[n][tn][1]]);
      }
    }
    fprintf(op, "\n");
  }
  fclose(op);
}

/* END trans_write()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Maintains the history log                                         */
/*-------------------------------------------------------------------*/
void history_log(master_t *master, int mode)
{
  parameters_t *params = master->params;
  FILE *dp;
  char buf[MAXSTRLEN], buf1[MAXSTRLEN], key[MAXSTRLEN], pname[MAXSTRLEN];
  int i, j, n, nn;
  long t;
  int maxdiff = 100;

  if (params->history & NONE) return;
  params->hrun = 1;

  /* Get the name of the parameter file                              */
  if (endswith(params->prmname, ".tran")) {
    j = 5;
    strcpy(buf, params->prmname);
  } else if (endswith(params->prmname, ".prm")) {
    j = 4;
    strcpy(buf, params->prmname);
  } else
    return;
  n = strlen(buf);
  for (i = 0; i < n-j; i++)
    pname[i] = buf[i];
  pname[i] = '\0';

  /* Open or create the history log                                  */
  if (mode == HST_PRE) {

    if (params->history & HST_LOG) {
      sprintf(buf, "%s.hist", pname);
      params->hstfd = NULL;
      strcpy(params->histname, buf);

      if (params->history & HST_RESET) {
	sprintf(key, "rm %s", buf); 
	if (params->history & HST_DIF) params->history &= ~HST_DIF;
	system(key);
      }
      if ((params->hstfd = fopen(params->histname, "r+")) == NULL) {
	params->hstfd = fopen(params->histname, "a");
	hd_warn("Creating history log file %s\n", buf);
	fprintf(params->hstfd, "\nHISTORY LOG for parameter file %s\n", params->prmname);
	if (j == 4) fprintf(params->hstfd, "Hydrodynamic model\n");
	if (j == 5) fprintf(params->hstfd, "Transport model\n");
	params->history |= HST_FST;
      } else {
	rewind(params->hstfd);
	sprintf(key, "Run%d:", params->hrun);
	while (prm_read_char(params->hstfd, key, buf)) {
	  params->hrun += 1;
	  sprintf(key, "Run%d:", params->hrun);
	}
	fclose(params->hstfd);
	params->hstfd = fopen(params->histname, "a");
      }
      fprintf(params->hstfd, "\n---------------------------\n");
      time(&t);
      fprintf(params->hstfd, "Run%d:              %s", params->hrun, ctime(&t));
      fprintf(params->hstfd, "EMS Version:       %s\n", version);
      fprintf(params->hstfd, "Executable file:   %s\n",executable );
      getcwd(buf, MAXSTRLEN);
      fprintf(params->hstfd, "Working directory: %s\n",buf);
      fprintf(params->hstfd, "ID_CODE:           %s\n", params->runcode);
      fprintf(params->hstfd, "Input file:        %s\n", params->idumpname);
      if (strlen(params->opath)) fprintf(params->hstfd, "Output path:       %s\n", params->opath);
      fprintf(params->hstfd, "Parameter header:  %s\n", params->parameterheader);
      if (strlen(params->notes))
      fprintf(params->hstfd, "Version notes:     %s\n", params->notes);
      fprintf(params->hstfd, "Start time:        %s\n", params->start_time);
      fprintf(params->hstfd, "Stop time:         %s\n", params->stop_time);
      if (forced_restart) fprintf(params->hstfd, "Restart run.\n");
      if (crash_restart) fprintf(params->hstfd, "Crash recovery run.\n");
      fflush(params->hstfd);

      if (params->history & HST_MASTER) {
	if ((dp = fopen(params->histnamem, "a+")) != NULL) {
	  fprintf(dp, "\n---------------------------\n");
	  time(&t);
	  if (j == 4) fprintf(dp, "Hydrodynamic model\n");
	  if (j == 5) fprintf(dp, "Transport model\n");
	  fprintf(dp, "Run%d:              %s", params->hrun, ctime(&t));
	  fprintf(dp, "EMS Version:       %s\n", version);
	  fprintf(dp, "Executable file:   %s\n",executable);
	  getcwd(buf, MAXSTRLEN);
	  fprintf(dp, "Working directory: %s\n",buf);
	  fprintf(dp, "Parameter file:    %s\n", params->prmname);
	  fprintf(dp, "ID_CODE:           %s\n", params->runcode);
	  fprintf(dp, "Input file:        %s\n", params->idumpname);
	  if (strlen(params->opath)) fprintf(dp, "Output path:       %s\n", params->opath);
	  fprintf(dp, "Parameter header:  %s\n", params->parameterheader);
	  if (strlen(params->notes))
	    fprintf(dp, "Version notes:     %s\n", params->notes);
	  fprintf(dp, "Start time:        %s\n", params->start_time);
	  fprintf(dp, "Stop time:         %s\n", params->stop_time);
	  
	  fclose(dp);
	}
      }

      if (params->history & HST_NOTES) {
	sprintf(buf, "%s.notes", pname);
	if (params->hrun == 1)
	  dp = fopen(buf, "w");
	else
	  dp = fopen(buf, "a+");
	if (!(prm_skip_to_end_of_key(dp, "Version notes summary"))) {
	  fprintf(dp, "\nVERSION NOTES SUMMARY for parameter file %s\n", params->prmname);
	  if (j == 4) fprintf(dp, "Hydrodynamic model\n");
	  if (j == 5) fprintf(dp, "Transport model\n");
	  fprintf(dp, "\n---------------------------\n");
	  fprintf(dp, "Run %d: ID_CODE: %s\n", params->hrun, params->runcode);
	  if (strlen(params->notes))
	    fprintf(dp, "  notes: %s\n", params->notes);
	  else
	    fprintf(dp, "  No NOTES specified for this run.\n");
	} else {
	  sprintf(key, "Run %d:", params->hrun-1);
	  prm_read_char(dp, key, buf);
	  if (prm_read_char(dp, "notes:", buf)) {
	    if (strlen(params->notes)) {
	      if (strcmp(buf, params->notes) != 0) {
		fprintf(dp, "Run %d: ID_CODE: %s\n", params->hrun, params->runcode);
		fprintf(dp, "  notes: %s\n", params->notes);
	      }
	    }
	  } else if (prm_skip_to_end_of_key(dp, "No notes.")) {
	    if (strlen(params->notes)) {
	      fprintf(dp, "Run %d: ID_CODE: %s\n", params->hrun, params->runcode);
	      fprintf(dp, "  notes: %s\n", params->notes);
	    }
	  }
	}
	fclose(dp);
      }
    }
    if (params->history & HST_DIF) {
      sprintf(key, "%s.txt", pname);
      if ((dp = fopen(key, "r")) == NULL) {
	sprintf(buf, "cp setup.txt %s", key); 
	if (system(buf)) return;
      } else {
	fclose(dp);
	if (params->history & HST_FST) {
	  sprintf(buf, "rm %s", key); 
	  system(buf);
	}
      }
    }
  }

  if (mode == HST_POST && params->history & HST_DIF) {
    if (!(params->history & HST_FST)) {
      char first[MAXSTRLEN];
      char pdiff[MAXSTRLEN];
      fprintf(params->hstfd, "Difference summary\n");
      sprintf(key, "%s.txt", pname);
      sprintf(buf, "diff setup.txt %s > setup.diff", key);
      system(buf);
      dp = fopen("setup.diff", "r");
      if (prm_read_char(dp, "<", first)) {
	fprintf(params->hstfd, "  New:    %s\n", first);
	if (prm_read_char(dp, ">", buf))
	  fprintf(params->hstfd, "  Old:    %s\n", buf);
	strcpy(pdiff, first);
      }
      j = 1;
      while(j) {
	if (prm_read_char(dp, "<", buf)) {
	  if (strcmp(buf, pdiff) == 0) {
	    fprintf(params->hstfd, "  recursive difference error: exiting....\n");
	    j = 0;
	    continue;
	  }
	  if (strcmp(buf, first) != 0) {
	    strcpy(pdiff, buf);
	    fprintf(params->hstfd, "  New:    %s\n",buf);
	    if (prm_read_char(dp, ">", buf))
	      fprintf(params->hstfd, "  Old:    %s\n",buf);
	    j++;
	  } else
	    j = 0;
	}
	if (j > maxdiff) {
	  fprintf(params->hstfd, "  truncating differences....\n");
	  j = 0;
	}
      }
      fflush(params->hstfd);
      fclose(dp);
      sprintf(buf, "cp setup.txt %s", key); 
      system(buf);
      sprintf(buf, "rm setup.diff");
      system(buf);
    } else
      params->history &= ~HST_FST;
  }

  if (mode == HST_OK && master->t == schedule->stop_time) {
    time(&t);
    if (params->hstfd != NULL) {
      fprintf(params->hstfd, "Run successful at %s\n", ctime(&t));
      fclose(params->hstfd);
    }
  }

  if (mode == HST_NOK) {
    time(&t);
    if (params->hstfd != NULL) {
      fprintf(params->hstfd, "Crashed %4.3f days: %s\n", master->days, ctime(&t));
      fclose(params->hstfd);
    }
  }
}

/* END history_log()                                                 */
/*-------------------------------------------------------------------*/
