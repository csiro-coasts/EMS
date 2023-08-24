/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/light_spectral_wc.c
 *  
 *  Description:
 *  
 *  Water column spectrally-resolved optical model. Model calculates:
 *  
 *  - scalar irradiance and absorption by each model autotroph in precalc (output time - ECOLOGY_DT).
 *  - simulated fluoresecence and turbidity in postcalc (output time)
 *  - the component of light returned to the surface from each layer for calculating simulated 
 *    satellite products in postcalc (output time).
 *
 *  CAVEATS: 1. variables last calculated in precalc (PAR etc.) are output at the output time but are the value from 
 *              output time - ECOLOGY_DT from calculation at variables calculated in postcalc are at the correct time.
 *           2. PAR_z is the downwelling irradiance on the z-grid, but when it is written to file, its dimension is
 *              given as cell z_centre. Infact its location is at the interface above.  
 *
 *  Options: light_spectral_wc(G|H|C|M|T,G|H,X|M)
 *  
 *  First argument:  G - GBR waters, B2p0 configuration.
 *                   H - GBR waters, B3p0 configuration (uses cdom_gbr tracer if available).
 *                   C - Chilean waters.
 *                   M - Macquarie Harbour.
 *                   T - Tasmanian shelf water (uses pale, amber and dark CDOM tracers).
 * 
 *  Second argument: G - Gaussian approx. of chl-a and non chl-a specific absorption.
 *                   H - HPLC-determined of pigment-specific absorption.
 *  
 *  Third argument:  X - no Moonlight calculation.
 *                   M - Moonlight from ROLO model.
 * 
 *  To output remote-sensing relectance at wavelenth XXX requires a 2D tracer named: R_XXX 
 *
 *  Optical model described in:
 *  
 *  Baird, M. E., N. Cherukuru, E. Jones, N. Margvelashvili, M. Mongin, K. Oubelkheir, P. J. Ralph, F. Rizwi, 
 *  B. J. Robson, T. Schroeder, J. Skerratt, A. D. L. Steven and K. A. Wild-Allen (2016) Remote-sensing 
 *  reflectance and true colour produced by a coupled hydrodynamic, optical, sediment, biogeochemical 
 *  model of the Great Barrier Reef, Australia: comparison with satellite data. Env. Model. Software 78: 79-96.
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: light_spectral_wc.c 7351 2023-04-30 03:42:14Z bai155 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ems.h>
#include "ecology_internal.h"
#include "ecofunct.h"
#include "utils.h"
#include "eprocess.h"
#include "column.h"
#include "cell.h"
#include "light_spectral_wc.h"
#include "constants.h"

double ginterface_calc_zenith(void *model, double t, int b);
double ginterface_getlighttop(void *model, int b);

typedef struct {

  char domain;
  char pig;
  char moon;

  /*
   * parameters
   */
  double m_l;
  double m_s;
  double m_MPB;
  double m_PhyD;
  double m_Tricho;
  double rad_l;
  double rad_s;
  double rad_MPB;
  double rad_PhyD;
  double rad_Tricho;
  double xan2chl_l;
  double xan2chl_s;
  double xan2chl_MPB;
  double xan2chl_PhyD;
  double xan2chl_Tricho;
  double bphy;
  double NtoCHL;
  int num_waves;
  double *wave;

  double *yC_s;
  double *yC_l;
  double *yC_Tricho;
  double *yC_MPB;
  double *yC_D;

  double *yC_s_rsr;
  double *yC_l_rsr;
  double *yC_Tricho_rsr;
  double *yC_MPB_rsr;
  double *yC_D_rsr;

  /*
   * tracers
   */
  int Kd_i;
  int PhyL_N_i;
  int PhyS_N_i;
  int PhyL_I_i;
  int PhyS_I_i;
  int PhyS_Chl_i;
  int PhyL_Chl_i;
  int DOR_C_i;
  int EFI_i;
  int CarbSand_i;
  int DetBL_N_i;
  int DetPL_N_i;
  int DetR_C_i;

  int Mud_mineral_i;
  int Mud_carbonate_i;
  int FineSed_i;
  int Mud_i;
  int Dust_i;

  int cdom_pale_i;
  int cdom_amber_i;
  int cdom_dark_i;
  int cdom_gbr_i;

  double Scdom_pale;
  double Scdom_amber;
  double Scdom_dark;

  double a440cdom_pale;
  double a440cdom_amber;
  double a440cdom_dark;

  int MPB_N_i;
  int MPB_Chl_i;
  int MPB_I_i;

  int PhyD_N_i;
  int PhyD_Chl_i;
  int PhyD_I_i;

  int Tricho_N_i;
  int Tricho_Chl_i;
  int Tricho_I_i;

  int at_440_i;
  int bt_550_i;
  int Kd_490_i;
  int bb_700_i;
  
  int PAR_i;
  int PAR_z_i;

  int Ed_440_i;
  int Ed_490_i;
  int Ed_550_i;
  int Ed_670_i;
  
  int if_glider;

  int ap_670_i;
  int Turbidity_i;
  
  int wPAR_bot;
  int w440;
  int w470;
  int w490;
  int w550;
  int w590;
  int w670;
  int w488;
  int w678;
  int w700;
  int wPAR_top;

  int w709;
  int w665;

  int Fluorescence_i;
  int BLP_i;

  int temp_i;
  int salt_i;

  int K_heat_i;

  int SMA_N_i;
  /*
   * Local variables
   */

  double MAleafden;

  /*
   * common variables
   */
  int cv_lighttop_s_i;
  int cv_zenith_i;

  int cv_lunar_zenith_i;
  int cv_lunar_phase_i;
  int cv_moonlight_i;
  int cv_moon_fulldisk_i;

  int z_secchi_i;
  int E_secchi_i;

  int KI_s_i;
  int KI_l_i;
  int KI_MPB_i;
  int KI_PhyD_i;
  int KI_Tricho_i;
  int KI_SMA_i;

  int yCfac_s_i;
  int yCfac_l_i;
  int yCfac_MPB_i;
  int yCfac_PhyD_i;
  int yCfac_Tricho_i;

  /*
   * Reflectances
   */
  int w_bot_i;
  int u_surf_i;

  double bbcell1;
  double bbcell2;
  double S_CDOM;
  double bbw_ratio;
  double bbnap_ratio;
  double SWRscale;
  double FLtoChl;

  /*
   * output variables
   */
  int zenith_i;
  
} workspace;


/* ginterface routines */
void ginterface_moonvars(void *hmodel, int c,
			 double *mlon, double *mlat,
			 double *dist_earth_sun, double *dist_moon_earth,
			 double *lunar_angle, double *sun_angle,double *moon_phase,double *lunar_dec);
double ginterface_get_cloud(void *hmodel, int c);

void light_spectral_wc_init(eprocess* p)
{
  ecology* e = p->ecology;
  stringtable* tracers = e->tracers;
  workspace* ws = malloc(sizeof(workspace));
  int large;
  char buf[MAXSTRLEN];
  int i,count_Li, count_Kd;

  char* prm = NULL;
  
  if (p->prms->n){
    prm = p->prms->se[0]->s;
    ws->domain = prm[0];
  }else{
    ws->domain = 'G';   // default is GBR waters, version 2p0.
  }

  ws->Mud_carbonate_i = -1;
  ws->Mud_mineral_i = -1;
  ws->FineSed_i = -1;
  ws->Dust_i = -1;
  ws->Mud_i = -1;
  ws->cdom_gbr_i = -1;
  
  switch (ws->domain){
  case 'G':
    eco_write_setup(e,"\nUsing BGC2p0 GBR optical properties \n");
    break;
  case 'H' :
    ws->Mud_carbonate_i = e->find_index(e->tracers, "Mud-carbonate", e);
    ws->Mud_mineral_i = e->find_index(e->tracers, "Mud-mineral", e);
    ws->FineSed_i = e->find_index(e->tracers, "FineSed", e);
    ws->Dust_i = e->find_index(e->tracers, "Dust", e);
    eco_write_setup(e,"\nUsing BGC3p0 mineral-based optical properties \n");
    ws->cdom_gbr_i = e->try_index(e->tracers, "cdom_gbr", e);
    if (ws->cdom_gbr_i > -1)
      eco_write_setup(e,"\nUsing cdom_gbr river tracer to determine CDOM absorption \n");
    break;
  case 'C' :
    eco_write_setup(e,"\nUsing Chilean optical properties \n");
    break;
  case 'T' :
    eco_write_setup(e,"\nUsing Tasmanian coastal waters optical properties: requires cdom_pale, cdom_amber, cdom_dark tracers \n");
    ws->cdom_pale_i = e->find_index(e->tracers, "cdom_pale", e);
    ws->cdom_amber_i = e->find_index(e->tracers, "cdom_amber", e);
    ws->cdom_dark_i = e->find_index(e->tracers, "cdom_dark", e);
    break;
  case 'M' :
    eco_write_setup(e,"\nUsing Macquarie Harbour optical properties: requires cdom_dark \n");
    ws->cdom_dark_i = e->find_index(e->tracers, "cdom_dark", e);
    break;
  default: 
    e->quitfn("ecology: error: \"%s\": \"%s(%s)\": unexpected domain for light model.", e->processfname, p->name, prm);
  }

  if (p->prms->n>1){
    prm = p->prms->se[1]->s;
    ws->pig = prm[0];
  }else{
    ws->pig = 'G';   // default is Guassian.
  }

  int is_geog = 0;
  if (!e->pre_build)
    is_geog = i_is_geographic(e->model);

  eco_write_setup(e,"Is the grid PROJECTION geographic (0 = no; 1= yes): %d\n",is_geog);
  
  if (p->prms->n>2){
    prm = p->prms->se[2]->s;
    ws->moon = prm[0];

    // Over-ride definition of moon if projection is not geographic

    if (is_geog == 0){
      eco_write_setup(e,"Moonlight calculation not undertaken because PROJECTION is not geographic.\n");
      ws->moon = 'X';
    }
    
  }else{
    ws->moon = 'X';   // default is no moolight.
  }
  
  p->workspace = ws;

  /*
   * parameters
   */

  ws->m_l = PhyCellMass(try_parameter_value(e, "PLrad"));
  ws->m_s = PhyCellMass(try_parameter_value(e, "PSrad"));

  ws->xan2chl_l = try_parameter_value(e, "PLxan2chl");
  if (isnan(ws->xan2chl_l))
    ws->xan2chl_l = 0.0;

  ws->xan2chl_s = try_parameter_value(e, "PSxan2chl");
  if (isnan(ws->xan2chl_s))
    ws->xan2chl_s = 0.0;

  ws->rad_l = get_parameter_value(e, "PLrad");
  ws->rad_s = get_parameter_value(e, "PSrad");

  ws->rad_MPB = try_parameter_value(e, "MBrad");
  ws->m_MPB = PhyCellMass(ws->rad_MPB);
  ws->xan2chl_MPB = try_parameter_value(e, "MBxan2chl");
  if (isnan(ws->xan2chl_MPB))
    ws->xan2chl_MPB = 0.0;

  ws->rad_PhyD = try_parameter_value(e, "DFrad");
  ws->m_PhyD = PhyCellMass(ws->rad_PhyD);
  ws->xan2chl_PhyD = try_parameter_value(e, "DFxan2chl");
  if (isnan(ws->xan2chl_PhyD))
    ws->xan2chl_PhyD = 0.0;

  ws->m_Tricho = PhyCellMass(try_parameter_value(e, "Tricho_rad"));
  ws->rad_Tricho = try_parameter_value(e, "Tricho_rad");
  ws->xan2chl_Tricho = try_parameter_value(e, "Trichoxan2chl");
  if (isnan(ws->xan2chl_Tricho))
    ws->xan2chl_Tricho = 0.0;

  ws->bphy = get_parameter_value(e, "bphy");
  ws->NtoCHL = get_parameter_value(e, "NtoCHL");

  /* Get wavelengths for light */

  ws->num_waves = get_parameter_num_values(e, "Light_lambda");
  ws->wave      = get_parameter_value_ptr(e, "Light_lambda");

  /*
   * tracers
   */

  eco_write_setup(e,"Optical code considering: ");

  ws->PhyL_N_i   = e->try_index(tracers, "PhyL_N", e);
  if (ws->PhyL_N_i > -1)
    eco_write_setup(e,"PhyL ");

  ws->PhyS_N_i   = e->try_index(tracers, "PhyS_N", e);
  if (ws->PhyS_N_i > -1)
    eco_write_setup(e,"PhyS ");

  ws->PhyL_I_i   = e->try_index(tracers, "PhyL_I", e);
  ws->PhyS_I_i   = e->try_index(tracers, "PhyS_I", e);
  ws->PhyS_Chl_i = e->try_index(tracers, "PhyS_Chl", e);
  ws->PhyL_Chl_i = e->try_index(tracers, "PhyL_Chl", e);
  ws->DOR_C_i    = e->try_index(tracers, "DOR_C", e);
  ws->EFI_i      = e->find_index(e->tracers, "EFI", e);
  ws->CarbSand_i = e->try_index(e->tracers, "CarbSand", e);
  ws->DetPL_N_i = e->try_index(tracers, "DetPL_N", e);
  ws->DetBL_N_i = e->try_index(tracers, "DetBL_N", e);
  ws->DetR_C_i = e->try_index(tracers, "DetR_C", e);

  ws->MPB_N_i   = e->try_index(tracers, "MPB_N", e);
  if (ws->MPB_N_i > -1)
    eco_write_setup(e,"MPB ");

  ws->MPB_Chl_i = e->try_index(tracers, "MPB_Chl", e);
  ws->MPB_I_i   = e->try_index(tracers, "MPB_I", e);

  ws->PhyD_N_i   = e->try_index(tracers, "PhyD_N", e);
  if (ws->PhyD_N_i > -1)
    eco_write_setup(e,"PhyD ");

  ws->PhyD_Chl_i = e->try_index(tracers, "PhyD_Chl", e);
  ws->PhyD_I_i   = e->try_index(tracers, "PhyD_I", e);

  ws->Tricho_N_i = e->try_index(tracers, "Tricho_N", e);
  if (ws->Tricho_N_i > -1)
    eco_write_setup(e,"Tricho ");

  eco_write_setup(e,"\n");

  ws->Tricho_Chl_i = e->try_index(tracers, "Tricho_Chl", e);
  ws->Tricho_I_i = e->try_index(tracers, "Tricho_I", e);

  ws->PAR_i = e->try_index(tracers, "PAR", e);

  ws->PAR_z_i = -1;

  ws->PAR_z_i = e->try_index(tracers, "PAR_z", e);

  ws->at_440_i   = e->try_index(tracers, "at_440", e);
  ws->bt_550_i = e->try_index(tracers, "bt_550", e);
  ws->Kd_490_i = e->try_index(tracers, "Kd_490", e);
  ws->K_heat_i = e->try_index(tracers, "K_heat", e);

  ws->Ed_440_i   = e->try_index(tracers, "Ed_440", e);
  ws->Ed_490_i   = e->try_index(tracers, "Ed_490", e);
  ws->Ed_550_i   = e->try_index(tracers, "Ed_550", e);
  ws->Ed_670_i   = e->try_index(tracers, "Ed_670", e);
  ws->bb_700_i   = e->try_index(tracers, "bb_700", e);
  

  ws->if_glider = 0;
  if ((ws->Ed_440_i > -1) && (ws->Ed_490_i > -1) && (ws->Ed_550_i > -1) && (ws->Ed_670_i > -1) && ws->bb_700_i > -1){
    ws->if_glider = 1;
    eco_write_setup(e,"Tracer lists contains Ed at 440, 490, 550, and 670, and bb at 700, so do glider calculations. \n");
  }
 
  ws->ap_670_i   = e->try_index(tracers, "ap_670", e);
  ws->Turbidity_i   = e->try_index(tracers, "Turbidity", e);

  ws->Fluorescence_i = e->try_index(tracers, "Fluorescence", e);
  ws->BLP_i = e->try_index(tracers, "BLP", e);

  /* get water mass properties for determining IOPs */

  ws->temp_i = e->find_index(tracers, "temp", e);
  ws->salt_i = e->find_index(tracers, "salt", e);

  ws->SMA_N_i = e->try_index(tracers, "SMA_N", e) ;

  /*
   * common column variables
   */
  ws->cv_lighttop_s_i = find_index_or_add_vector(e->cv_column, "lighttop_s", e, ws->num_waves);

  ws->cv_zenith_i = find_index_or_add(e->cv_column, "zenith", e);

  ws->z_secchi_i = find_index_or_add(e->cv_column, "z_secchi", e);
  ws->E_secchi_i = find_index_or_add(e->cv_column, "E_secchi", e);
  
  ws->bbcell1 = try_parameter_value(e, "bbcell1");
  if (isnan(ws->bbcell1))
    ws->bbcell1 = 4.26e-14;   // Whitmore 2012
  ws->bbcell2 = try_parameter_value(e, "bbcell2");
  if (isnan(ws->bbcell2))
    ws->bbcell2 = 2.028;   // Whitmore 2012

  ws->bbw_ratio = try_parameter_value(e, "bbw_ratio");
  if (isnan(ws->bbw_ratio)){
    ws->bbw_ratio = 0.5;
    eco_write_setup(e,"Code default of clear water backscatter bbw_ratio = %e \n",ws->bbw_ratio);
  }

  if (ws->domain == 'H'||ws->domain == 'C'||ws->domain == 'M'||ws->domain == 'T')
    eco_write_setup(e,"Spectrally-resolved NAP backscattering ratio used. \n");
  else{
    ws->bbnap_ratio = try_parameter_value(e, "bbnap_ratio");
    if (isnan(ws->bbnap_ratio)){
      ws->bbnap_ratio = 0.025;  // Blondeau-Patissier 2009 JGR: 114: C05003.
      eco_write_setup(e,"Code default of bbnap_ratio = %e \n",ws->bbnap_ratio);
    }
  }

  if (ws->domain == 'H'||ws->domain == 'C'||ws->domain == 'G'){
    ws->S_CDOM = try_parameter_value(e, "S_CDOM");
    if (isnan(ws->S_CDOM)){
      ws->S_CDOM = 0.012;  // Blondeau-Patissier 2009 JGR: 114: C05003.
      eco_write_setup(e,"Code default of S_CDOM = %e \n",ws->S_CDOM);
    }
  }
  if (ws->domain == 'T'){
    ws->Scdom_pale = try_parameter_value(e, "Scdom_pale");
    if (isnan(ws->Scdom_pale)){
      ws->Scdom_pale = 0.010;
      eco_write_setup(e,"Code default of Scdom_pale = %e \n",ws->Scdom_pale);
    }
    ws->Scdom_amber = try_parameter_value(e, "Scdom_amber");
    if (isnan(ws->Scdom_amber)){
      ws->Scdom_amber = 0.0125;
      eco_write_setup(e,"Code default of Scdom_amber = %e \n",ws->Scdom_amber);
    }
    ws->Scdom_dark = try_parameter_value(e, "Scdom_dark");
    if (isnan(ws->Scdom_dark)){
      ws->Scdom_dark = 0.015;
      eco_write_setup(e,"Code default of Scdom_dark = %e \n",ws->Scdom_dark);
    }
    ws->a440cdom_pale = try_parameter_value(e, "a440cdom_pale");
    if (isnan(ws->a440cdom_pale)){
      ws->a440cdom_pale = 0.1;
      eco_write_setup(e,"Code default of a440cdom_pale = %e \n",ws->a440cdom_pale);
    }
    ws->a440cdom_amber = try_parameter_value(e, "a440cdom_amber");
    if (isnan(ws->a440cdom_amber)){
      ws->a440cdom_amber = 1.0;
      eco_write_setup(e,"Code default of a440cdom_amber = %e \n",ws->a440cdom_amber);
    }
    ws->a440cdom_dark = try_parameter_value(e, "a440cdom_dark");
    if (isnan(ws->a440cdom_dark)){
      ws->a440cdom_dark = 10.0;
      eco_write_setup(e,"Code default of a440cdom_dark = %e \n",ws->a440cdom_dark);
    }
  }
  if (ws->domain == 'M'){
    ws->Scdom_dark = try_parameter_value(e, "Scdom_dark");
    if (isnan(ws->Scdom_dark)){
      ws->Scdom_dark = 0.015;
      eco_write_setup(e,"Code default of Scdom_dark = %e \n",ws->Scdom_dark);
    }
    ws->a440cdom_dark = try_parameter_value(e, "a440cdom_dark");
    if (isnan(ws->a440cdom_dark)){
      ws->a440cdom_dark = 10.0;
      eco_write_setup(e,"Code default of a440cdom_dark = %e \n",ws->a440cdom_dark);
    }
  }


  /* Extra parameter for DA perturbations */

  ws->SWRscale = try_parameter_value(e, "SWRscale");
  if (isnan(ws->SWRscale)){
    ws->SWRscale = 1.0;
    eco_write_setup(e,"Code default of SWRscale = %e \n",ws->SWRscale);
  }

  ws->FLtoChl = try_parameter_value(e, "FLtoChl");
  if (isnan(ws->FLtoChl)){
    ws->FLtoChl = 1.0;
    eco_write_setup(e,"Code default of FLtoChl = %e \n",ws->FLtoChl);
  }


  /*
   * Output variables
   */
  ws->zenith_i = e->try_index(tracers, "Zenith", e);

}
void ecology_find_rsr_waves(ecology *e);
void light_spectral_wc_postinit(eprocess* p)
{
  double *absorbance;
  double *absorbance_xanth;
  ecology* e = p->ecology;
  workspace* ws = (workspace *)p->workspace;
  int large,w;

  if (e->pre_build) return;

  int is_geog = 0;
  is_geog = i_is_geographic(e->model);

  if (process_present(e,PT_WC,"recom_extras") && is_geog==1 ){
    ws->moon = 'M';
  }
  
  if (ws->moon == 'M'){
    
    eco_write_setup(e,"\nUsing ROLO moonlight model: Kieffer, H. H and T. C. Stone (2005). \n\n");
    
    ws->cv_lunar_zenith_i = find_index_or_add(e->cv_column, "lunar_zenith", e);
    ws->cv_lunar_phase_i = find_index_or_add(e->cv_column, "lunar_phase", e);
    ws->cv_moonlight_i = find_index_or_add(e->cv_column, "moonlight", e);
    ws->cv_moon_fulldisk_i = find_index_or_add(e->cv_column, "moon_fulldisk", e);
    
  }else{
    eco_write_setup(e,"\nNo moonlight forcing \n\n");
  }

  /* 
   * These can only be called when postinit is initiated from the
   * actual ecology_build not pre_build
   * Note: do not call try/find_index on tracers in this function
   */

  if (e->bio_opt==NULL){
    ecology_find_rsr_waves(e);
    e->bio_opt = bio_opt_init(e);
  }

  /* set light indexes so that code isn't looking for them each time step. */

  ws->w_bot_i  = find_index_or_add_vector(e->cv_column, "w_bot", e, e->num_rsr_waves);
  ws->u_surf_i = find_index_or_add_vector(e->cv_column, "u_surf", e, e->num_rsr_waves);

  /* Search for special wavelengths */

  ws->wPAR_bot = -1;
  ws->w440 = -1;
  ws->w470 = -1;
  ws->w490 = -1;
  ws->w550 = -1;
  ws->w590 = -1;
  ws->w670 = -1;
  ws->wPAR_top = -1;
  
  for (w=0; w<ws->num_waves; w++){
    if (ws->wave[w] > 400.0 && ws->wPAR_bot < 0){
      ws->wPAR_bot = w;
      eco_write_setup(e,"Bottom of PAR range is centred on %4.2f nm with edge %4.2f nm \n",ws->wave[w],(ws->wave[w-1]+ws->wave[w])/2.0);
    }
    if (ws->wave[w] == 440.0)
      ws->w440 = w;
    if (ws->wave[w] == 470.0)
      ws->w470 = w;
    if (ws->wave[w] == 490.0)
      ws->w490 = w;
    if (ws->wave[w] == 550.0)
      ws->w550 = w;
    if (ws->wave[w] == 590.0)
      ws->w590 = w;
    if (ws->wave[w] == 670.0)
      ws->w670 = w;
    if (ws->wave[w] > 700.0 && ws->wPAR_top < 0){
      ws->wPAR_top = w-1;
      eco_write_setup(e,"Top of PAR range is centred on %4.2f nm with edge %4.2f nm \n",ws->wave[w-1],(ws->wave[w-1]+ws->wave[w])/2.0);
    }
  }
  
  ws->w488 = -1;
  ws->w678 = -1;
  ws->w700 = -1;
  for (w=0; w<e->num_rsr_waves; w++){
    if (e->rsr_waves[w] == 590.0)
      ws->w590 = w;
    if (e->rsr_waves[w] == 470.0)
      ws->w470 = w;
    if (e->rsr_waves[w] == 488.0)
      ws->w488 = w;
    if (e->rsr_waves[w] == 678.0)
      ws->w678 = w;
    if (e->rsr_waves[w] == 700.0)
      ws->w700 = w;
  }

  if (ws->SMA_N_i > -1){
    
    /* Put in macroalgae initialisations if there is a macroalgae tracer */
    
    ws->MAleafden = get_parameter_value(e, "MAleafden_wc");
    
  }

  ws->KI_s_i = try_index(e->cv_cell, "KI_s", e);
  ws->KI_l_i = try_index(e->cv_cell, "KI_l", e);
  ws->KI_MPB_i = try_index(e->cv_cell, "KI_MPB", e);
  ws->KI_PhyD_i = try_index(e->cv_cell, "KI_PhyD", e);
  ws->KI_Tricho_i = try_index(e->cv_cell, "KI_Tricho", e);
  
  ws->KI_SMA_i = try_index(e->cv_cell, "KI_SMA", e);
  
  ws->yCfac_s_i = try_index(e->cv_cell, "yCfac_s", e);
  ws->yCfac_l_i = try_index(e->cv_cell, "yCfac_l", e);
  ws->yCfac_MPB_i = try_index(e->cv_cell, "yCfac_MPB", e);
  ws->yCfac_PhyD_i = try_index(e->cv_cell, "yCfac_PhyD", e);
  ws->yCfac_Tricho_i = try_index(e->cv_cell, "yCfac_Tricho", e);
  
  /* YC: CALCULATE (PRM = GUAS) OR READ (PRM = HPLC) */
  
  bio_opt_prop *bio = e->bio_opt;
  
  int num_waves = ws->num_waves;
  int num_rsr_waves = e->num_rsr_waves;
  
  ws->yC_s = d_alloc_1d(num_waves);
  ws->yC_l = d_alloc_1d(num_waves);
  ws->yC_Tricho = d_alloc_1d(num_waves);
  ws->yC_MPB = d_alloc_1d(num_waves);
  ws->yC_D = d_alloc_1d(num_waves);

  if (e->num_rsr_waves) {
    ws->yC_s_rsr = d_alloc_1d(num_rsr_waves);
    ws->yC_l_rsr = d_alloc_1d(num_rsr_waves);
    ws->yC_Tricho_rsr = d_alloc_1d(num_rsr_waves);
    ws->yC_MPB_rsr = d_alloc_1d(num_rsr_waves);
    ws->yC_D_rsr = d_alloc_1d(num_rsr_waves);
  } 
  
  if (ws->pig == 'H'){  // HPLC determined absorption coefficients

    eco_write_setup(e,"\nMass-specific absorption coefficients used for microalgae in water column \n");

    for (w=0; w<num_waves; w++){
      ws->yC_s[w] = bio->yC_picoplankton[w];
      ws->yC_l[w] = bio->yC_microplankton[w];
      ws->yC_Tricho[w] = bio->yC_tricho[w];
      ws->yC_MPB[w] = bio->yC_microplankton[w];
      ws->yC_D[w] = bio->yC_microplankton[w];
    }

    for (w=0; w<e->num_rsr_waves; w++){
      ws->yC_s_rsr[w] = bio->yC_picoplankton_rsr[w];
      ws->yC_l_rsr[w] = bio->yC_microplankton_rsr[w];
      ws->yC_Tricho_rsr[w] = bio->yC_tricho_rsr[w];
      ws->yC_MPB_rsr[w] = bio->yC_microplankton_rsr[w];
      ws->yC_D_rsr[w] = bio->yC_microplankton_rsr[w];
    }

    // Overwrite picophytoplankton pigments with rhodomonas_duplex pigments.

    if (ws->domain == 'M'){
      eco_write_setup(e,"\nUsing rhodomonas duplex pigment composition for small phytoplankton \n");
      for (w=0; w<num_waves; w++){
	ws->yC_s[w] = bio->yC_rhodomonas_duplex[w];
      }
      for (w=0; w<e->num_rsr_waves; w++){
	ws->yC_s_rsr[w] = bio->yC_rhodomonas_duplex_rsr[w];
      }
    }
    
  }else{   // Guassian curve determined absorption coefficients 
    
    eco_write_setup(e,"\nGuassian curve determined absorption coefficients in water column \n");
    
    double *yC_chl;
    double *yC_xanth;
    double *yC_chl_rsr;
    double *yC_xanth_rsr;
    
    yC_chl = d_alloc_1d(num_waves);
    yC_xanth = d_alloc_1d(num_waves); 
    absorbwave(yC_chl,1.0, num_waves, ws->wave);
    absorbxanth(yC_xanth,1.0, num_waves, ws->wave);
    
    if (e->num_rsr_waves) {
      yC_chl_rsr = d_alloc_1d(num_rsr_waves);
      yC_xanth_rsr = d_alloc_1d(num_rsr_waves);
      absorbwave(yC_chl_rsr,1.0, num_rsr_waves, e->rsr_waves);
      absorbxanth(yC_xanth_rsr,1.0, num_rsr_waves, e->rsr_waves);
    }
    
    for (w=0; w<num_waves; w++){
      ws->yC_s[w] = yC_chl[w] + ws->xan2chl_s * yC_xanth[w];
      ws->yC_l[w] = yC_chl[w] + ws->xan2chl_l * yC_xanth[w];
      ws->yC_Tricho[w] = yC_chl[w] + ws->xan2chl_Tricho * yC_xanth[w];
      ws->yC_MPB[w] = yC_chl[w] + ws->xan2chl_MPB * yC_xanth[w];
      ws->yC_D[w] = yC_chl[w] + ws->xan2chl_PhyD * yC_xanth[w];
    }
    
    for (w=0; w<num_rsr_waves; w++){
      ws->yC_s_rsr[w] = yC_chl_rsr[w] + ws->xan2chl_s * yC_xanth_rsr[w];
      ws->yC_l_rsr[w] = yC_chl_rsr[w] + ws->xan2chl_l * yC_xanth_rsr[w];
      ws->yC_Tricho_rsr[w] = yC_chl_rsr[w] + ws->xan2chl_Tricho * yC_xanth_rsr[w];
      ws->yC_MPB_rsr[w] = yC_chl_rsr[w] + ws->xan2chl_MPB * yC_xanth_rsr[w];
      ws->yC_D_rsr[w] = yC_chl_rsr[w] + ws->xan2chl_PhyD * yC_xanth_rsr[w];
    }
    
    d_free_1d(yC_chl);
    d_free_1d(yC_xanth);
    if (e->num_rsr_waves) {
      d_free_1d(yC_chl_rsr);
      d_free_1d(yC_xanth_rsr);  
    }  
  }
}

  void light_spectral_wc_destroy(eprocess* p)
{
  workspace* ws = p->workspace;
  d_free_1d(ws->yC_s);
  d_free_1d(ws->yC_l);
  d_free_1d(ws->yC_Tricho);
  d_free_1d(ws->yC_MPB);
  d_free_1d(ws->yC_D);
  free(ws);
}

void light_spectral_wc_precalc(eprocess* p, void* pp)
{
  ecology* e = p->ecology;
  cell* c = (cell*) pp;
  column* col = c->col;
  workspace* ws = p->workspace;
  int num_waves = ws->num_waves;
  double *wave  = ws->wave;
  double* y = c->y;    
  bio_opt_prop *bio = e->bio_opt;
  /* 
   * The column common variables are now vectors so they must be
   * accessed as below. i.e. first get a pointer to the variable and
   * then index into it
   */
  double *lighttop_s = col->cv[ws->cv_lighttop_s_i];
  int w,w2;
  double energy = 0.0;

  if (isnan(lighttop_s[0])) {

    double light = ws->SWRscale * ginterface_getlighttop(col->model, col->b); // *albedo;

    /* check that albedo is zero in the transport file */
    
    for (w=0; w<num_waves; w++) {
      lighttop_s[w] =  bio->landa[w]*light;
    }

    if (ws->moon == 'M'){

      /* Calculate moonlight using USGS ROLO model - see bio_opt.c for spectral coefficients */
      
      /* This model accounts for surface variation in albedo (i.e. position of highlands and seas), as 
	 well as complicated planetary motions and atmospheric effects. */ 
      
      /* variables to be obtained from tide calculations - all angles in radians */
      
      double *cv_moonlight = col->cv[ws->cv_moonlight_i];
      double *cv_lunar_zenith = col->cv[ws->cv_lunar_zenith_i];
      double *cv_lunar_phase = col->cv[ws->cv_lunar_phase_i];
      double *cv_moon_fulldisk = col->cv[ws->cv_moon_fulldisk_i];
      
      // Default values :
      
      double moonlat_obs = 0.0; // B from moon_pos
      double moonlon_obs = 0.0; // L from moon_pos
      double earth_sun_dist = 1496.0e8; // m 
      double moon_earth_dist = 384400.0*1000.0;  // m
      double lunar_angle = 0.0;
      double sun_angle = 0.0;
      double moon_phase = 0.0;
      double lunar_dec = 0.0;
      
      ginterface_moonvars(col->model,col->b+1, &moonlon_obs, &moonlat_obs, &earth_sun_dist,&moon_earth_dist, &lunar_angle, &sun_angle, &moon_phase,&lunar_dec);
      
      // Could have an if statement for moon below horizon
      
      // printf("moon: lat %e, lon %e, E-S dist %e, M-E dist %e, lunar_angle %e, solar angle %e phase %e \n",moonlat_obs,moonlon_obs,earth_sun_dist,moon_earth_dist,lunar_angle,sun_angle,moon_phase);
      
      // Mark to add extra code:
      
      double gmo = fabs(moon_phase-M_PI/2.0);

      // Still not given a good number. 
      
      double moonlon_sun = 7.0/2.0/M_PI; // pi - (solar hour angle - monlon_obs) excluding librations.
      // need to be careful with going past pi.
      
      // **********************************************************************************************************
      // wrap this into a function, in matlab syntax:
      
      // Ak[w] = moon_full_disk_reflectance(moonlat_obs,moonlon_obs,earth_sun_dist,moon_earth_dist,gmo,moonlon_sun,num_waves)
      
      // gmo = gmo * 180.0/M_PI;
      double moonlon_obs_deg = moonlon_obs*180.0/M_PI;
      double moonlat_obs_deg = moonlat_obs*180.0/M_PI;
      double moonlon_sun_deg = 0.0;
      
      double c1 = 0.00034115;
      double c2 = -0.0013425;
      double c3 = 0.00095906;
      double c4 = 0.00066229;
      double p1 = 4.06054*M_PI/180.0;
      double p2 = 12.8802*M_PI/180.0;
      double p3 = -30.5858*M_PI/180.0;
      double p4 = 16.7498*M_PI/180.0;
      
      double Ak[num_waves];
      
      double Ak_nonwave = c1*moonlat_obs_deg + c2*moonlon_obs_deg +c3*moonlon_sun_deg*moonlat_obs_deg + c4*moonlon_sun_deg*moonlon_obs_deg;
      
      for (w=0; w<num_waves; w++){
	Ak[w] = bio->a0[w] + bio->a1[w] * gmo + bio->a2[w] * pow(gmo,2.0) + bio->a3[w] * pow(gmo,3.0);
	Ak[w] += bio->b1[w] * moonlon_sun_deg + bio->b2[w] * pow(moonlon_sun_deg,3.0) + bio->b3[w] * pow(moonlon_sun_deg,5.0);
	Ak[w] += Ak_nonwave;
	Ak[w] += bio->d1[w] * exp(-gmo/p1) + bio->d2[w] * exp(-gmo/p2) + bio->d3[w] * cos((gmo-p3)/p4);
	Ak[w] = exp(Ak[w]);
      }
      
      // *********************************************************************************************************
      
      double cloud = ginterface_get_cloud(col->model,col->b);

      /* Mike's albedo function of cloud cover and hour angle */
 
      /* Get the solar elevation */
      /* hrang *= pi / 180.0; */
      /* dec*=pi/180.0; */
      /* lat*=pi/180.0; */
      // h = sin(moonlat_obs) * sin(lunar_dec) + cos(moonlat_obs) * cos(lunar_dec) * cos(lunar_angle);
      // if (h < 0.0)
      // h = 0.0;
      
      /*-----------------------------------------------------------------*/
      /* Get the albedo of the sea surface */
      // h = asin(h) * 180.0 / pi;
      // iw = (int)(h / 10.0);
      // d1 = (double)iw *10.0;
      // if (iw == 0)
      // alp = ap[tcld][0];
      // else
      // alp = ((ap[tcld][iw + 1] - ap[tcld][iw]) / 10.0) * (h - d1) + ap[tcld][iw];

      // *********************************************************************************************************

      double attenuation_by_clouds = 1.0;  // write interface function to clouds.
      
      double airmass = 1.0; // 1.0/(cos(lunar_zenith)+0.15*pow(93.885-lunar_zenith*180.0/M_PI,-1.253));

      double transmission = airmass;

      /* Use solid angle (6.4177e-5 sr), sunlight (attenuation_by_clouds * 1366.0 * bio->landa[w]), normalised earth / sun / moon 
         distances and above calculated full disk reflectance */
      
      double moonlight = 0.0;
      double moon_tmp = 0.0;
      
      for (w=0; w<num_waves; w++) {
	moon_tmp = Ak[w] * 6.4177e-5 * 1366.0 * bio->landa[w] * pow(earth_sun_dist/1496.0e8 ,2.0) * pow(moon_earth_dist/384400000.0,2.0)/M_PI;
	lighttop_s[w] = lighttop_s[w] + moon_tmp;
	moonlight += moon_tmp;
      }
      /* Correct zenith following solar calcs. */
      
      double lunar_zenith = sin(lunar_dec)*sin(moonlat_obs)+cos(lunar_dec)*cos(moonlat_obs)*cos(lunar_angle);
      
      lunar_zenith = acos(lunar_zenith);
      
      lunar_zenith = (fabs(lunar_zenith)<1.0e-15)? 1.0e-15:lunar_zenith;
      lunar_zenith = ((lunar_zenith-M_PI/2.0)-fabs(lunar_zenith-M_PI/2.0))/2.0+M_PI/2.0;
      
      cv_lunar_zenith[0] = lunar_zenith;
      cv_lunar_phase[0] = moon_phase;
      cv_moon_fulldisk[0] = (lunar_zenith > M_PI/2.0 - 1.0e-15)? 0:moonlight; // zero if below horizon.
      cv_moonlight[0] = moonlight * cos(lunar_zenith) * attenuation_by_clouds * transmission;
      
      // printf("lunar_dec %e, lunar_angle %e, lunar_zenith %e, cos(lunar_zenith) %e \n",lunar_dec,lunar_angle,lunar_zenith,cos(lunar_zenith)); 
    }
  }
  

  /* if night don't do light calculations */

  for (w=ws->wPAR_bot; w <= ws->wPAR_top; w++){
    energy += lighttop_s[w];
  }

  if (energy < 0.0){
    
    for (w=0; w<num_waves; w++) {
      lighttop_s[w] =  0.0;
    }
    if (ws->KI_SMA_i>-1){
      c->cv[ws->KI_SMA_i] = 0.0;
    }
    if (ws->KI_s_i>-1){
      c->cv[ws->KI_s_i] = 0.0;   
      c->cv[ws->yCfac_s_i] = 0.0; 
    }
    if (ws->KI_l_i>-1){
      c->cv[ws->KI_l_i] = 0.0;
      c->cv[ws->yCfac_l_i] = 0.0;
    }
    if (ws->KI_MPB_i>-1){
      c->cv[ws->KI_MPB_i] = 0.0;
      c->cv[ws->yCfac_MPB_i] = 0.0;
    }  
    if (ws->KI_PhyD_i>-1){
      c->cv[ws->KI_PhyD_i] = 0.0;
      c->cv[ws->yCfac_PhyD_i] = 0.0;
    }
    if (ws->KI_Tricho_i>-1){
      c->cv[ws->KI_Tricho_i] = 0.0;
      c->cv[ws->yCfac_Tricho_i] = 0.0;
    }
    if (ws->at_440_i > -1){ // i.e. outputting
      y[ws->at_440_i] = 0.0;
    }
    if (ws->bt_550_i > -1){ // i.e. outputting 
      y[ws->bt_550_i] = 0.0;
    }
    if (ws->Kd_490_i > -1){ // i.e. outputting
      y[ws->Kd_490_i] = 0.0;
    }
    if (ws->PAR_i > -1){ // i.e. outputting
      y[ws->PAR_i] = 0.0;
    }  
    if (ws->PAR_z_i > -1){ // i.e. outputting
      y[ws->PAR_z_i] = 0.0;
    }
    if (ws->K_heat_i > -1){ // i.e. outputting
      y[ws->K_heat_i] = 0.04;  // clear water value. 
    }
    if (ws->if_glider){
      y[ws->Ed_440_i] = 0.0;
      y[ws->Ed_490_i] = 0.0;
      y[ws->Ed_550_i] = 0.0;
      y[ws->Ed_670_i] = 0.0;
      y[ws->bb_700_i] = 0.0;
    }
    return;
   }

  double dz = c->dz_wc;
  double lighttop; 
  double zenith, thetaw;
  double  light_s[num_waves];
  double *zenith_cv = col->cv[ws->cv_zenith_i];
  double *at_s;   /* spectrally-resolved total absorption [m-1] */
  double *bt_s;   /* spectrally-resolved total scattering [m-1] */
  double delenergy = 0.0;
  energy = 0.0; 
  int i;

  /*
   * Allocate local spectral variables
   */

  at_s = d_alloc_1d(num_waves);
  bt_s = d_alloc_1d(num_waves);

  double Kd_s[num_waves];

  double DOC = 0.0;
  double EFI = 0.0;
  double CarbSand = 0.0;
  double DetBL_N = 0.0;
  double DetPL_N = 0.0;
  double DetR_C = 0.0;
  double NAP = 0.0;
  double NAP_noEFI = 0.0;
  double photons = 0.0;

  double Mud_carbonate = 0.0;
  double Mud_mineral = 0.0;
  double Dust = 0.0;
  double FineSed = 0.0;
  double Mud = 0.0;
  
  if (ws->SMA_N_i > -1){
  
    /* Assume that Macroalgae sits at the top of the layer. */
    
    double SMA_N = y[ws->SMA_N_i];
    c->cv[ws->KI_SMA_i] = 0.0;
    
    if (SMA_N > 0.0){  /* macroalgae */
      
      for (w=0; w<num_waves; w++){
	photons += lighttop_s[w] * (1.0 - exp(-SMA_N * c->dz_wc * ws->MAleafden * bio->MA_aAwave[w])) 
	  * wave[w];
	lighttop_s[w] = lighttop_s[w] * exp(-SMA_N * c->dz_wc * ws->MAleafden * bio->MA_aAwave[w]);
      }
      c->cv[ws->KI_SMA_i] = photons * 8.359335857479461e-09 ;
    }
  }

  if (ws->DOR_C_i > -1)
    DOC = y[ws->DOR_C_i];
  if (ws->EFI_i > -1)
    EFI = y[ws->EFI_i];
  if (ws->CarbSand_i > -1)
    CarbSand = y[ws->CarbSand_i];
  if (ws->Mud_i > -1)
    Mud = y[ws->Mud_i];
  if (ws->DetBL_N_i > -1)
    DetBL_N = y[ws->DetBL_N_i];
  if (ws->DetPL_N_i > -1)
    DetPL_N = y[ws->DetPL_N_i];
  if (ws->DetR_C_i > -1)
    DetR_C = y[ws->DetR_C_i];

  if (ws->Mud_carbonate_i > -1)
    Mud_carbonate = y[ws->Mud_carbonate_i];
  if (ws->Mud_mineral_i > -1)
    Mud_mineral = y[ws->Mud_mineral_i];
  if (ws->Dust_i > -1)
    Dust = y[ws->Dust_i];
  if (ws->FineSed_i > -1)
    FineSed = y[ws->FineSed_i];

  /* Need to convert units of detritus to kg m-3. */

  NAP_noEFI = (DetBL_N * atk_W_C + DetPL_N * red_W_C + DetR_C)*1e-6; 
  NAP = EFI + NAP_noEFI;


  /*
   * Initialise zenith on the very first time
   */
  zenith = zenith_cv[0];

  if (isnan(zenith)){

    zenith = ginterface_calc_zenith(e->model, e->t, c->b);
    
    /*  avoid exactly zero zenith which cause trouble with albedo calcs. */
    
    zenith = (fabs(zenith)<1.0e-15)? 1.0e-15:zenith;
    zenith = ((zenith-M_PI/2.0)-fabs(zenith-M_PI/2.0))/2.0+M_PI/2.0;
    
    /* Assign to the zenith column variable */
    
    zenith_cv[0] = zenith;
   
  }
  
  /* Calculate angle of light through water using Snell's law */

  thetaw = asin(sin(zenith)/1.33);

  double gone = e->bio_opt->gone;
  double gtwo = e->bio_opt->gtwo;

  double costhetaw = cos(thetaw);
  
  double sumlight;
  double yCfac;

  double *aA_s_l = NULL;
  double *aA_s_s  = NULL;
  double *aA_s_MPB  = NULL;
  double *aA_s_PhyD = NULL ;
  double *aA_s_Tricho  = NULL;

  if (ws->KI_s_i>-1){
    
    c->cv[ws->KI_s_i] = 0.0;   
    c->cv[ws->yCfac_s_i] = 0.0;
    
    /* if variable Chl need to recalculate absorbance */
    
    double rad = ws->rad_s;
    double vol = 4.0 / 3.0 * M_PI * rad * rad * rad;
    double m = ws->m_s;    
    double Phy_N = y[ws->PhyS_N_i];
    double Phy_Chl = y[ws->PhyS_Chl_i];
    
    double cellnum =  Phy_N / (m * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
    double cellChl = Phy_Chl / (vol * cellnum); /* cellular Chl conc */
    
    // intermediate var
    
    double *absorbance;
    
    absorbance = d_alloc_1d(num_waves);
    aA_s_s = d_alloc_1d(num_waves);
    
    /* weight self-shading coefficient based on photons incident on the top of 
       the layer - multiply by wave[w] to convert energy to photons, but ignore constant 
       (1e9 h c)-1 / Av is constant across all wavelengths */
    
    sumlight = 0.0;
    yCfac = 0.0;
    
    for (w=0; w<num_waves; w++){
      absorbance[w] = ws->yC_s[w] * cellChl;
      yCfac = yCfac + diffaa(absorbance[w],rad)*max(lighttop_s[w]*wave[w],0.00001);
      sumlight = sumlight + max(lighttop_s[w]*wave[w],0.00001);
    }
    
    c->cv[ws->yCfac_s_i] = yCfac/sumlight;
    
    aawave(rad, absorbance, aA_s_s, num_waves,wave);
    
    // free intermediate var
    d_free_1d(absorbance);
  }
  
  if (ws->KI_l_i>-1){
        
    c->cv[ws->KI_l_i] = 0.0;
    c->cv[ws->yCfac_l_i] = 0.0;
    
    /* if variable Chl need to recalculate absorbance */

    double rad = ws->rad_l;
    double vol = 4.0 / 3.0 * M_PI * rad * rad * rad;
    double m   = ws->m_l;
    
    double Phy_N = y[ws->PhyL_N_i];
    double Phy_Chl = y[ws->PhyL_Chl_i];
    
    double cellnum =  Phy_N / (m * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
    double cellChl = Phy_Chl / (vol * cellnum); /* cellular Chl conc */
    
    // intermediate var
    
    double *absorbance;    
    absorbance = d_alloc_1d(num_waves);
    aA_s_l = d_alloc_1d(num_waves);

    /* weight self-shading coefficient based on photons incident on the top of 
       the layer - multiply by wave[w] to convert energy to photons, but ignore constant 
       (1e9 h c)-1 / Av is constant across all wavelengths */
    
    sumlight = 0.0;
    yCfac = 0.0;
    
    for (w=0; w<num_waves; w++){
      absorbance[w] = ws->yC_l[w] * cellChl;
      yCfac = yCfac + diffaa(absorbance[w],rad)*max(lighttop_s[w]*wave[w],0.00001);
      sumlight = sumlight + max(lighttop_s[w]*wave[w],0.00001);
    }
    c->cv[ws->yCfac_l_i] = yCfac/sumlight;
    aawave(rad, absorbance, aA_s_l, num_waves,wave);
        
    // free intermediate var
    d_free_1d(absorbance);
  }

  if (ws->KI_MPB_i>-1){
    
    c->cv[ws->KI_MPB_i] = 0.0;
    c->cv[ws->yCfac_MPB_i] = 0.0;

      /* now do it for microphytobenthos */

    double rad = ws->rad_MPB;
    double vol = 4.0 / 3.0 * M_PI * rad * rad * rad;
    double m   = ws->m_MPB;
    
    double MPB_N = y[ws->MPB_N_i];
    double MPB_Chl = y[ws->MPB_Chl_i];
    
    double cellnum =  MPB_N / (m * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
    double cellChl = MPB_Chl / (vol * cellnum); /* cellular Chl conc */
    
    // intermediate var
    
    double *absorbance;    
    absorbance = d_alloc_1d(num_waves);
    aA_s_MPB = d_alloc_1d(num_waves);

    /* weight self-shading coefficient based on photons incident on the top of 
       the layer - multiply by wave[w] to convert energy to photons, but ignore constant 
       (1e9 h c)-1 / Av which is constant across all wavelengths */
    
    /* if it is dark, assume equal weighting for all wavelengths. */
    
    sumlight = 0.0;
    yCfac = 0.0;
    
    for (w=0; w<num_waves; w++){
      absorbance[w] = ws->yC_MPB[w] * cellChl;
      yCfac = yCfac + diffaa(absorbance[w],rad)*max(lighttop_s[w]*wave[w],0.00001);
      sumlight = sumlight + max(lighttop_s[w]*wave[w],0.00001);
    }
    c->cv[ws->yCfac_MPB_i] = yCfac/sumlight;
    
    aawave(rad, absorbance, aA_s_MPB, num_waves,wave);
    
    // free intermediate var
    d_free_1d(absorbance);
  }
  
  if (ws->KI_PhyD_i>-1){
    
    c->cv[ws->KI_PhyD_i] = 0.0;
    c->cv[ws->yCfac_PhyD_i] = 0.0;
    
    /* now do it for dinoflagellates - use parwaves from above */
    
    double rad = ws->rad_PhyD;
    double vol = 4.0 / 3.0 * M_PI * rad * rad * rad;
    double m   = ws->m_PhyD;
    
    double PhyD_N = y[ws->PhyD_N_i];
    double PhyD_Chl = y[ws->PhyD_Chl_i];
    
    double cellnum =  PhyD_N / (m * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
    double cellChl = PhyD_Chl / (vol * cellnum); /* cellular Chl conc */
    
    // intermediate var
    
    double *absorbance;
    
    absorbance = d_alloc_1d(num_waves);
    aA_s_PhyD = d_alloc_1d(num_waves);

    /* calculate intracellular self-shading coefficients in the Chl a absorption bands */
    
    /* weight self-shading coefficient based on photons incident on the top of 
       the layer - multiply by wave[w] to convert energy to photons, but ignore constant 
       (1e9 h c)-1 / Av is constant across all wavelengths */
    
    sumlight = 0.0;
    yCfac = 0.0;
    
    for (w=0; w<num_waves; w++){
      absorbance[w] = ws->yC_D[w] * cellChl;
      yCfac = yCfac + diffaa(absorbance[w],rad)*max(lighttop_s[w]*wave[w],0.00001);
      sumlight = sumlight + max(lighttop_s[w]*wave[w],0.00001);
    }
    c->cv[ws->yCfac_PhyD_i] = yCfac/sumlight;
    
    aawave(rad, absorbance, aA_s_PhyD, num_waves,wave);
    
    // free intermediate var
    d_free_1d(absorbance);
  }
  
    if (ws->KI_Tricho_i>-1){
      
      c->cv[ws->KI_Tricho_i] = 0.0;
      c->cv[ws->yCfac_Tricho_i] = 0.0;
      
      /* now do it for Trichodesmium - use parwaves from above */
      
      double rad = ws->rad_Tricho;
      double vol = 4.0 / 3.0 * M_PI * rad * rad * rad;
      double m   = ws->m_Tricho;
      
      double Tricho_N = y[ws->Tricho_N_i];
      double Tricho_Chl = y[ws->Tricho_Chl_i];
      
      double cellnum =  Tricho_N / (m * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
      double cellChl = Tricho_Chl / (vol * cellnum); /* cellular Chl conc */
      
      // intermediate var
      
      double *absorbance;
      absorbance = d_alloc_1d(num_waves);
      aA_s_Tricho = d_alloc_1d(num_waves);
      
      /* calculate intracellular self-shading coefficients in the Chl a absorption bands */
      
      /* weight self-shading coefficient based on photons incident on the top of 
	 the layer - multiply by wave[w] to convert energy to photons, but ignore constant 
	 (1e9 h c)-1 / Av is constant across all wavelengths */
      
      sumlight = 0.0;
      yCfac = 0.0;
      
      for (w=0; w<num_waves; w++){
        absorbance[w] = ws->yC_Tricho[w] * cellChl;
	yCfac = yCfac + diffaa(absorbance[w],rad)*max(lighttop_s[w]*wave[w],0.00001);
	sumlight = sumlight + max(lighttop_s[w]*wave[w],0.00001);
      }
      c->cv[ws->yCfac_Tricho_i] = yCfac/sumlight;
      
      aawave(rad, absorbance, aA_s_Tricho, num_waves,wave);
      
      // free intermediate var
      d_free_1d(absorbance);
  }
  
  double acdom443 = e->bio_opt->acdom443star * DOC;

  /* store light at top of the layer in PAR */

  if (ws->PAR_z_i > -1){ /* i.e. outputting */
    double ttmmpp = 0.0;
    for (w=ws->wPAR_bot; w <= ws->wPAR_top; w++){
      ttmmpp += lighttop_s[w] * 8.359335857479461e-09 * wave[w];
    }	
    y[ws->PAR_z_i] = ttmmpp;
  }
  
  for (w=0; w<num_waves; w++){

    lighttop = lighttop_s[w]; // save lighttop_s as lighttop so it can be updated.

    /* 
     * include all calculations for attenuating components, ignoring how
     * it is done in earlier versions
     */

    at_s[w] = bio->kw_s[w] ;  /*  clear water */     

    if (ws->domain == 'G' || ws->domain == 'H'){
      
      /* or alternatively use CDOM vs salinity relationship from Schroeder (2012) 
	 Mar. Poll. Bull 65:210-223 with CDOM spectral slope from 
	 Blondeau-Patissier 2009 JGR: 114: C05003.                 */

      /* added if < 34.145 to avoid extrapolation that may result in negative absorption at above 37, and 
         to avoid deep salt minimums acting like estuarine waters */

      //if (y[ws->salt_i] < 34.145){
      // 	at_s[w] += (1.2336 + (-0.0332 * y[ws->salt_i])) * exp(-ws->S_CDOM*(wave[w]-443.0));
      // }
      //      at_s[w] += 0.01 * exp(-ws->S_CDOM*(wave[w]-443.0)); // low CDOM absorption in marine waters.

      if (ws->cdom_gbr_i > -1){
	//if (y[ws->salt_i] < 36.855){
	  double a443 = (1.2336 + (-0.0332 * (1.0-y[ws->cdom_gbr_i])*36.855));//;+ 0.01; // 0.01 open-ocean a_cdom(443)
	  double S_cdom = 0.0061 * pow(a443,-0.309); // Blondeau-Patissier 2009 JGR: 114: C05003.
	  at_s[w] += a443 * exp(-S_cdom * (wave[w]-443.0));
	  //}else{
	  //at_s[w] += 0.01 * exp(-0.010 * (wave[w]-443.0));
	  //}
      }else{
	if (y[ws->salt_i] < 34.145){
	  double a443 = (1.2336 + (-0.0332 * y[ws->salt_i])) + 0.01; // 0.01 open-ocean a_cdom(443)
	  double S_cdom = 0.0061 * pow(a443,-0.309); // Blondeau-Patissier 2009 JGR: 114: C05003.
	  at_s[w] += a443 * exp(-S_cdom * (wave[w]-443.0));
	}else{
	  at_s[w] += 0.01 * exp(-0.010 * (wave[w]-443.0));
	}
      }
    }
    if (ws->domain == 'C'){
      // From cruises in Oct 2017, Mar 2018 - based on 440 nm
      // at_s[w] += (0.2875 - 0.0059 * min(y[ws->salt_i],35.0)) * exp(-ws->S_CDOM*(wave[w]-440.0));
      at_s[w] += (0.2875 - 0.0059 * min(y[ws->salt_i],35.0)) 
	* exp(-(ws->S_CDOM - max(y[ws->salt_i]-30.0,0.0)*0.007/5.0 ) * (wave[w]-440.0));
    }

    if (ws->domain == 'T'){
      at_s[w] += 0.01 * exp(-0.014 * (wave[w]-440.0));
      at_s[w] += ws->a440cdom_pale * min(y[ws->cdom_pale_i],1.0) * exp(- ws->Scdom_pale * (wave[w]-440.0));
      at_s[w] += ws->a440cdom_amber * min(y[ws->cdom_amber_i],1.0) * exp(- ws->Scdom_amber * (wave[w]-440.0));
      at_s[w] += ws->a440cdom_dark * min(y[ws->cdom_dark_i],1.0) * exp(- ws->Scdom_dark * (wave[w]-440.0));
    }
    
    if (ws->domain == 'M'){
      at_s[w] += 0.01 * exp(-0.014 * (wave[w]-440.0));
      at_s[w] += ws->a440cdom_dark * min(y[ws->cdom_dark_i],1.0) * exp(- ws->Scdom_dark * (wave[w]-440.0));
    }

    switch (ws->domain){     /* Still need to think through NAP */
    case 'G' :
      at_s[w] += bio->ad_A[w] * pow(NAP, bio->ad_B[w]);
      at_s[w] += bio->ad_CarbSand_A[w] * pow(CarbSand, bio->ad_CarbSand_B[w]);
      break;
    case 'H' :
      // at_s[w] += bio->ad_A[w] * (Dust+Mud_mineral+FineSed+NAP_noEFI); // Gladstone Harbour mud
      // at_s[w] += bio->aC_CAL2[w] * Mud_carbonate * 1000.0 ; // Calcite
      at_s[w] += bio->LJCOc_aNAP[w] * Mud_carbonate;
      at_s[w] += bio->LJCO_aNAP[w] * (Dust+Mud_mineral+FineSed+NAP_noEFI);
      break;
    case 'C' :
      at_s[w] += bio->aC_ICE2[w] * NAP * 1.0e3;
      break;
    case 'T' :
      at_s[w] += bio->aC_ICE2[w] * NAP * 1.0e3;
      break;
    case 'M' :
      at_s[w] += bio->aC_ICE2[w] * NAP * 1.0e3;
      break;
    }
    
    bt_s[w] =  bio->bw_s[w];  /* Clear water scattering */

    /* Effect of phytoplankton on absorption - no. cells x absorption cross-section */

    /* Effect of phytoplankton on scattering, assuming constant NtoCHL, since 
       scattering depends on particles (as quantified by PhyL_N) rather than 
       pigment content - may require further thought  */

    if (ws->KI_l_i > -1){
      at_s[w] += (y[ws->PhyL_N_i] / (ws->m_l*red_A_N * 1000.0 * MW_Nitr) ) * aA_s_l[w];
      bt_s[w] += ws->bphy * y[ws->PhyL_N_i] / ws->NtoCHL; 
    }
    if (ws->KI_s_i > -1){
      at_s[w] += (y[ws->PhyS_N_i] / (ws->m_s*red_A_N * 1000.0 * MW_Nitr) ) * aA_s_s[w];
      bt_s[w] += ws->bphy * y[ws->PhyS_N_i] / ws->NtoCHL ;
    }

    /* Note that phytoplankton only add to absorption and scattering if KI_MPB_i, which will only 
       be defined if they are growing */

    if (ws->KI_MPB_i > -1){
      at_s[w] += (y[ws->MPB_N_i] / (ws->m_MPB*red_A_N * 1000.0 * MW_Nitr) ) * aA_s_MPB[w];
      bt_s[w] += ws->bphy * y[ws->MPB_N_i] / ws->NtoCHL ;
    }
    if (ws->KI_PhyD_i > -1){
      at_s[w] += (y[ws->PhyD_N_i] / (ws->m_PhyD*red_A_N * 1000.0 * MW_Nitr) ) * aA_s_PhyD[w];
      bt_s[w] += ws->bphy * y[ws->PhyD_N_i] / ws->NtoCHL ;
    }
    if (ws->KI_Tricho_i > -1){
      at_s[w] += (y[ws->Tricho_N_i] / (ws->m_Tricho*red_A_N * 1000.0 * MW_Nitr) ) * aA_s_Tricho[w];
      bt_s[w] += ws->bphy * y[ws->Tricho_N_i] / ws->NtoCHL ;
    }

    switch (ws->domain){
    case 'G':
     bt_s[w] += bio->bbp_A[w] * pow(NAP, bio->bbp_B[w]); /* NAP scattering */ 
     bt_s[w] += bio->bbp_CarbSand_A[w] * pow(NAP, bio->bbp_CarbSand_B[w]); /* NAP scattering */ 
     break;
    case 'H' :
      //bt_s[w] += bio->bC_CAL2[w] * Mud_carbonate * 1.0e3;
      //bt_s[w] += bio->bbp_A[w] * (Dust+Mud_mineral+FineSed+NAP_noEFI);
      bt_s[w] += bio->LJCOc_bNAP[w] * Mud_carbonate; /* NAP backscattering */
      bt_s[w] += bio->LJCO_bNAP[w] * (Dust + Mud_mineral+FineSed+NAP_noEFI);
      break;
    case 'C' :
      bt_s[w] += bio->bC_ICE2[w] * NAP * 1.0e3;
      break;
    case 'T' :
      bt_s[w] += bio->bC_ICE2[w] * NAP * 1.0e3;
      break;
    case 'M' :
      bt_s[w] += bio->bC_ICE2[w] * NAP * 1.0e3;
     break;
    }
    /*  
     * include the effect of azimuth angle and scattering on
     * pathlength through layer:
     */
    
    double adlen = sqrt(1.0 + (gone*costhetaw-gtwo) * bt_s[w] / at_s[w]) / costhetaw;

    Kd_s[w] = at_s[w] *  adlen;   /* vertical attenuation */

    lighttop_s[w] = lighttop * exp(-Kd_s[w] * dz );

    light_s[w] = (lighttop - lighttop_s[w] ) / (Kd_s[w] * dz );


    /* Also, Escalar = Edown * adlen */

    /* 8.359335857479461e-9 = (1e9 h c)-1 / Av 
     where h = Plank const, c = speed of light, Av = Avogadro const, 1e9 nm m-1  */
    if (ws->KI_s_i > -1){
      c->cv[ws->KI_s_i] += light_s[w] * 8.359335857479461e-09*wave[w] * 
	aA_s_s[w] * adlen;
    }
    if (ws->KI_l_i > -1){
      c->cv[ws->KI_l_i] += light_s[w] * 8.359335857479461e-09*wave[w] * 
	aA_s_l[w] * adlen;
    }
    if (ws->KI_MPB_i > -1){
      c->cv[ws->KI_MPB_i] += light_s[w] * 8.359335857479461e-09*wave[w] * 
	aA_s_MPB[w] * adlen;
    }
    if (ws->KI_PhyD_i > -1){
      c->cv[ws->KI_PhyD_i] += light_s[w] * 8.359335857479461e-09*wave[w] * 
	aA_s_PhyD[w] * adlen;
    }
    if (ws->KI_Tricho_i > -1){
      c->cv[ws->KI_Tricho_i] += light_s[w] * 8.359335857479461e-09*wave[w] * 
	aA_s_Tricho[w] * adlen;
    }
    
    /* Integrated energy loss across layer for PAR: [J s-1 m-2 layer thick-1] */

    if (wave[w] <= 700.0){ // only PAR light heating through water column.
      energy += lighttop; 
      delenergy += (lighttop-lighttop_s[w]);
    }
  } // finish w loop.

  /* Calculate a VERTICAL attenuation coefficient for heat */

  if (ws->K_heat_i > -1){
    if (delenergy > 0.0001 && fabs(energy-delenergy) > 1e-5) { 
      y[ws->K_heat_i] = - log((energy-delenergy)/energy) / dz;
    }
    else{    /* if dark set to default */
      y[ws->K_heat_i] = 0.04;
    }
  } 
  
  /* output IOPs */

  if (ws->at_440_i > -1){ /* i.e. outputting */
    y[ws->at_440_i] = at_s[ws->w440];
  }
  
  if (ws->bt_550_i > -1){ /* i.e. outputting */
    y[ws->bt_550_i] = bt_s[ws->w550];
  }

  if (ws->Kd_490_i > -1){ /* i.e. outputting */
    y[ws->Kd_490_i] = Kd_s[ws->w490];
  }

  if (ws->ap_670_i > -1){ /* total absorption - water at 670 nm */
    y[ws->ap_670_i] = at_s[ws->w670] - bio->kw_s[ws->w670];
  }
  if (ws->if_glider){   // need to divide by waveband width:  W m-2 nm-1
      y[ws->Ed_440_i] = lighttop_s[ws->w440]/10.0;
      y[ws->Ed_490_i] = lighttop_s[ws->w490]/20.0;
      y[ws->Ed_550_i] = lighttop_s[ws->w550]/20.0;
      y[ws->Ed_670_i] = lighttop_s[ws->w670]/20.0;
    }

  /*/ reuse energy variable, but as mol photon m-2 s-1*/
  
  energy = 0.0;
  
  if (ws->PAR_i > -1){ /* i.e. outputting */
    for (w=ws->wPAR_bot; w <= ws->wPAR_top; w++){
      energy += light_s[w] * 8.359335857479461e-09 * wave[w];
    }	
    y[ws->PAR_i] = energy;
  }   
  
  d_free_1d(at_s);
  d_free_1d(bt_s);

  if (aA_s_l != NULL)
        d_free_1d(aA_s_l);
  if (aA_s_s != NULL)
        d_free_1d(aA_s_s);
  if (aA_s_MPB != NULL)
        d_free_1d(aA_s_MPB);
  if (aA_s_PhyD != NULL)
        d_free_1d(aA_s_PhyD);
  if (aA_s_Tricho != NULL)
        d_free_1d(aA_s_Tricho);

}

void light_spectral_wc_postcalc(eprocess* p, void* pp)
{
  ecology* e = p->ecology;
  cell* c = (cell*) pp;
  column* col = c->col;
  workspace* ws = p->workspace;
  double* y = c->y;

  double *lighttop_s = col->cv[ws->cv_lighttop_s_i];
  double nFLH = 0.0;

  /*
   * Section to calculate reflectances - these are done on output fields, with the exception, in transport mode, or salt.
   */

  if (e->num_rsr_waves) {
    
    if (ws->Fluorescence_i > -1  && ws->w470 > -1){

      /* Calculate simulated fluorescence based on spherical, packaged cells absorbing at 470 nm. */

      /* Passive fluorescence (for nFLH) includes photochemical quenching, but active does not (FL) */ 

      double y470 = 0.03;
      double FL = 0.0;
      
      if (ws->KI_s_i>-1){
		
	double rad = ws->rad_s;
	double vol = 4.0 / 3.0 * M_PI * rad * rad * rad;
	double m = ws->m_s;
	
	double Phy_N = y[ws->PhyS_N_i];
	double Phy_Chl = y[ws->PhyS_Chl_i];
	
	double cellnum =  Phy_N / (m * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
	double cellChl = Phy_Chl / (vol * cellnum); /* cellular Chl conc */

	double Phy_I = y[ws->PhyS_I_i];
	double Iquota = max(0.0,min(1.0,(Phy_I / cellnum) / m * red_A_I * 1000.0)); 
	
	double aA_470 = 0.0;
	double temp_470 = 2.0 * ws->yC_s_rsr[ws->w470] * cellChl * rad;
	
	if (temp_470 > 3e-4){
	  aA_470 =  M_PI * rad * rad * (1.0 - 2.0 * (1.0 - (1.0 + temp_470) * exp(-temp_470)) / (temp_470 * temp_470));
	}

	FL = FL + aA_470 / (vol * y470 * cellChl) * Phy_Chl * (1.0 - 0.5 * Iquota) ;

	nFLH += c->cv[ws->KI_s_i] * (1.0 - 0.5 * Iquota) * Iquota * cellnum;
      }

      if (ws->KI_l_i>-1){
		
	double rad = ws->rad_l;
	double vol = 4.0 / 3.0 * M_PI * rad * rad * rad;
	double m = ws->m_l;
	
	double Phy_N = y[ws->PhyL_N_i];
	double Phy_Chl = y[ws->PhyL_Chl_i];
	
	double cellnum =  Phy_N / (m * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
	double cellChl = Phy_Chl / (vol * cellnum); /* cellular Chl conc */

	double Phy_I = y[ws->PhyL_I_i];
	double Iquota = max(0.0,min(1.0,(Phy_I / cellnum) / m * red_A_I * 1000.0));
	
	double aA_470 = 0.0;
	
	double temp_470 = 2.0 * ws->yC_l_rsr[ws->w470] * cellChl * rad;
	
	if (temp_470 > 3e-4){
	  aA_470 =  M_PI * rad * rad * (1.0 - 2.0 * (1.0 - (1.0 + temp_470) * exp(-temp_470)) / (temp_470 * temp_470));
	}

	FL = FL + aA_470 / (vol * y470 * cellChl) * Phy_Chl * (1.0 - 0.5 * Iquota) ;

	nFLH += c->cv[ws->KI_l_i] * (1.0 - 0.5 * Iquota) * Iquota * cellnum;

      }

      if (ws->KI_MPB_i>-1){
		
	double rad = ws->rad_MPB;
	double vol = 4.0 / 3.0 * M_PI * rad * rad * rad;
	double m = ws->m_MPB;
	
	double Phy_N = y[ws->MPB_N_i];
	double Phy_Chl = y[ws->MPB_Chl_i];
	
	double cellnum =  Phy_N / (m * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
	double cellChl = Phy_Chl / (vol * cellnum); /* cellular Chl conc */

	double Phy_I = y[ws->MPB_I_i];
	double Iquota = max(0.0,min(1.0,(Phy_I / cellnum) / m * red_A_I * 1000.0));
	
	double aA_470 = 0.0;
	double temp_470 = 2.0 * ws->yC_l_rsr[ws->w470] * cellChl * rad;
	
	if (temp_470 > 3e-4){
	  aA_470 =  M_PI * rad * rad * (1.0 - 2.0 * (1.0 - (1.0 + temp_470) * exp(-temp_470)) / (temp_470 * temp_470));
	}

	FL = FL + aA_470 / (vol * y470 * cellChl) * Phy_Chl * (1.0 - 0.5 * Iquota) ;

	nFLH += c->cv[ws->KI_MPB_i] * (1.0 - 0.5 * Iquota) * Iquota * cellnum;

      }

      if (ws->KI_Tricho_i>-1){
		
	double rad = ws->rad_Tricho;
	double vol = 4.0 / 3.0 * M_PI * rad * rad * rad;
	double m = ws->m_Tricho;
	
	double Phy_N = y[ws->Tricho_N_i];
	double Phy_Chl = y[ws->Tricho_Chl_i];
	
	double cellnum =  Phy_N / (m * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
	double cellChl = Phy_Chl / (vol * cellnum); /* cellular Chl conc */

	double Phy_I = y[ws->Tricho_I_i];
	double Iquota = max(0.0,min(1.0,(Phy_I / cellnum) / m * red_A_I * 1000.0));
	
	double aA_470 = 0.0;
	double temp_470 = 2.0 * ws->yC_Tricho_rsr[ws->w470] * cellChl * rad;
	
	if (temp_470 > 3e-4){
	  aA_470 =  M_PI * rad * rad * (1.0 - 2.0 * (1.0 - (1.0 + temp_470) * exp(-temp_470)) / (temp_470 * temp_470));
	}

	FL = FL + aA_470 / (vol * y470 * cellChl) * Phy_Chl * (1.0 - 0.5 * Iquota) ;

	nFLH += c->cv[ws->KI_Tricho_i] * (1.0 - 0.5 * Iquota) * Iquota * cellnum;

      }

      if (ws->KI_PhyD_i>-1){
		
	double rad = ws->rad_PhyD;
	double vol = 4.0 / 3.0 * M_PI * rad * rad * rad;
	double m = ws->m_PhyD;
	
	double Phy_N = y[ws->PhyD_N_i];
	double Phy_Chl = y[ws->PhyD_Chl_i];
	
	double cellnum =  Phy_N / (m * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
	double cellChl = Phy_Chl / (vol * cellnum); /* cellular Chl conc */

	double Phy_I = y[ws->PhyD_I_i];
	double Iquota = max(0.0,min(1.0,(Phy_I / cellnum) / m * red_A_I * 1000.0));
	
	double aA_470 = 0.0;
	double temp_470 = 2.0 * ws->yC_Tricho_rsr[ws->w470] * cellChl * rad;
	
	if (temp_470 > 3e-4){
	  aA_470 =  M_PI * rad * rad * (1.0 - 2.0 * (1.0 - (1.0 + temp_470) * exp(-temp_470)) / (temp_470 * temp_470));
	}

	FL = FL + aA_470 / (vol * y470 * cellChl) * Phy_Chl * (1.0 - 0.5 * Iquota) ;

	nFLH += c->cv[ws->KI_PhyD_i] * (1.0 - 0.5 * Iquota) * Iquota * cellnum;

      }

      y[ws->Fluorescence_i] = ws->FLtoChl * FL / 1.333333333334; /* 1.34 for spheres, see Baird et al. (2013) */

    }

    double *at_r;   /* spectrally-resolved total absorption [m-1] */
    double *bt_r;   /* spectrally-resolved total scattering [m-1] */
    double *bb_r;   /* spectrally-resolved backscattering   [m-1] */
    
    double adlen;
    
    double dz = c->dz_wc;
    bio_opt_prop *bio = e->bio_opt;
    
    double gone = e->bio_opt->gone;
    double gtwo = e->bio_opt->gtwo;

    double *z_secchi_cv = col->cv[ws->z_secchi_i];
    double *E_secchi_cv = col->cv[ws->E_secchi_i];
    
    int w2; 
    
    /* Get zenith at postcalc time */
    
    double zenith = ginterface_calc_zenith(e->model, e->t + e->dt, c->b);
    
    if (ws->zenith_i > -1) // this is the 3D zenith output for backward compatability.
      y[ws->zenith_i] = zenith;

    zenith = (fabs(zenith)<1.0e-15)? 1.0e-15:zenith;
    zenith = ((zenith-M_PI/2.0)-fabs(zenith-M_PI/2.0))/2.0+M_PI/2.0;
    
    /* if (zenith > (M_PI/2.0 - 0.01)){
      for (w2=0; w2<e->num_rsr_waves; w2++) {
	col->cv[ws->w_bot_i][w2] = 0.0;
      }
      return;
      } */
    
    double thetaw = asin(sin(zenith)/1.33);
    
    double costhetaw = cos(thetaw);
    
    double *aA_s_l = NULL;
    double *aA_s_s  = NULL;
    double *aA_s_MPB  = NULL;
    double *aA_s_PhyD = NULL ;
    double *aA_s_Tricho  = NULL;
    
    double* w_bot  = col->cv[ws->w_bot_i];
    double* u_surf = col->cv[ws->u_surf_i];
    
    double top,Kd_r,cellnum_l,cellnum_s,cellnum_PhyD,cellnum_Tricho,cellnum_MPB;
    
    double DOC = 0.0;
    double EFI = 0.0;
    double CarbSand = 0.0;
    double DetBL_N = 0.0;
    double DetPL_N = 0.0;
    double DetR_C = 0.0;
    double NAP = 0.0;
    double NAP_noEFI = 0.0;

    double Mud_carbonate = 0.0;
    double Mud_mineral = 0.0;
    double FineSed = 0.0;
    double Dust = 0.0;
    double Mud = 0.0;
    
    if (ws->DOR_C_i > -1)
      DOC = y[ws->DOR_C_i];
    if (ws->EFI_i > -1)
      EFI = y[ws->EFI_i];
    if (ws->CarbSand_i > -1)
      CarbSand = y[ws->CarbSand_i];
    if (ws->DetBL_N_i > -1)
      DetBL_N = y[ws->DetBL_N_i];
    if (ws->DetPL_N_i > -1)
      DetPL_N = y[ws->DetPL_N_i];
    if (ws->DetR_C_i > -1)
      DetR_C = y[ws->DetR_C_i];

    if (ws->Mud_carbonate_i > -1)
      Mud_carbonate = y[ws->Mud_carbonate_i];
    if (ws->Mud_mineral_i > -1)
      Mud_mineral = y[ws->Mud_mineral_i];
    if (ws->FineSed_i > -1)
      FineSed = y[ws->FineSed_i];
    if (ws->Mud_i > -1)
      Mud = y[ws->Mud_i];
    if (ws->Dust_i > -1)
      Dust = y[ws->Dust_i];
    
    NAP_noEFI = (DetBL_N * atk_W_C + DetPL_N * red_W_C + DetR_C)*1e-6;
    NAP = EFI + NAP_noEFI; 

    bb_r = d_alloc_1d(e->num_rsr_waves);
    bt_r = d_alloc_1d(e->num_rsr_waves);
    at_r = d_alloc_1d(e->num_rsr_waves);

    if (isnan(w_bot[0])) {
      for (w2=0; w2<e->num_rsr_waves; w2++){
	w_bot[w2] = 1.0;
	u_surf[w2] = 0.0;
      }
    }

    /* add phytoplankton at reflectance wavelengths */

    if (ws->KI_l_i>-1){
      
      double rad = ws->rad_l;
      double vol = 4.0 / 3.0 * M_PI * rad * rad * rad;
      double m   = ws->m_l;
      
      double Phy_N = y[ws->PhyL_N_i];
      double Phy_Chl = y[ws->PhyL_Chl_i];
      
      cellnum_l =  Phy_N / (m * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
      double cellChl = Phy_Chl / (vol * cellnum_l); /* cellular Chl conc */
      
      double *absorbance;
      
      absorbance = d_alloc_1d(e->num_rsr_waves);

      aA_s_l = d_alloc_1d(e->num_rsr_waves);
      
      for (w2=0; w2<e->num_rsr_waves; w2++){
	absorbance[w2] = ws->yC_l_rsr[w2] * cellChl;
      }
      aawave(rad, absorbance, aA_s_l, e->num_rsr_waves,e->rsr_waves);

      // free intermediate var
      d_free_1d(absorbance);
    }

    if (ws->KI_s_i>-1){
      
      double rad = ws->rad_s;
      double vol = 4.0 / 3.0 * M_PI * rad * rad * rad;
      double m   = ws->m_s;
      
      double Phy_N = y[ws->PhyS_N_i];
      double Phy_Chl = y[ws->PhyS_Chl_i];
      
      cellnum_s =  Phy_N / (m * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
      double cellChl = Phy_Chl / (vol * cellnum_s); /* cellular Chl conc */
      
      double *absorbance;
      
      absorbance = d_alloc_1d(e->num_rsr_waves);

      aA_s_s = d_alloc_1d(e->num_rsr_waves);
      
      for (w2=0; w2<e->num_rsr_waves; w2++){
	absorbance[w2] = ws->yC_s_rsr[w2] * cellChl;
      }
      aawave(rad, absorbance, aA_s_s, e->num_rsr_waves,e->rsr_waves);

      // free intermediate var
      d_free_1d(absorbance);
    }
    if (ws->KI_MPB_i>-1){
      
      double rad = ws->rad_MPB;
      double vol = 4.0 / 3.0 * M_PI * rad * rad * rad;
      double m   = ws->m_MPB;
      
      double Phy_N = y[ws->MPB_N_i];
      double Phy_Chl = y[ws->MPB_Chl_i];
      
      cellnum_MPB =  Phy_N / (m * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
      double cellChl = Phy_Chl / (vol * cellnum_MPB); /* cellular Chl conc */
      
      double *absorbance;
      
      absorbance = d_alloc_1d(e->num_rsr_waves);

      aA_s_MPB = d_alloc_1d(e->num_rsr_waves);
      
      for (w2=0; w2<e->num_rsr_waves; w2++){
	absorbance[w2] = ws->yC_MPB_rsr[w2] * cellChl;
      }
      aawave(rad, absorbance, aA_s_MPB, e->num_rsr_waves,e->rsr_waves);

      // free intermediate var
      d_free_1d(absorbance);
    }

    if (ws->KI_PhyD_i>-1){
      
      double rad = ws->rad_PhyD;
      double vol = 4.0 / 3.0 * M_PI * rad * rad * rad;
      double m   = ws->m_PhyD;
      
      double Phy_N = y[ws->PhyD_N_i];
      double Phy_Chl = y[ws->PhyD_Chl_i];
      
      cellnum_PhyD =  Phy_N / (m * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
      double cellChl = Phy_Chl / (vol * cellnum_PhyD); /* cellular Chl conc */
      
      double *absorbance;
      
      absorbance = d_alloc_1d(e->num_rsr_waves);

      aA_s_PhyD = d_alloc_1d(e->num_rsr_waves);
      
      for (w2=0; w2<e->num_rsr_waves; w2++){
	absorbance[w2] = ws->yC_D_rsr[w2] * cellChl;
      }
      aawave(rad, absorbance, aA_s_PhyD, e->num_rsr_waves,e->rsr_waves);

      // free intermediate var
      d_free_1d(absorbance);
    }

    if (ws->KI_Tricho_i>-1){
      
      double rad = ws->rad_Tricho;
      double vol = 4.0 / 3.0 * M_PI * rad * rad * rad;
      double m   = ws->m_Tricho;
      
      double Phy_N = y[ws->Tricho_N_i];
      double Phy_Chl = y[ws->Tricho_Chl_i];
      
      cellnum_Tricho =  Phy_N / (m * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
      double cellChl = Phy_Chl / (vol * cellnum_Tricho); /* cellular Chl conc */

      double *absorbance;
      
      absorbance = d_alloc_1d(e->num_rsr_waves);

      aA_s_Tricho = d_alloc_1d(e->num_rsr_waves);
      
      for (w2=0; w2<e->num_rsr_waves; w2++){
	absorbance[w2] = ws->yC_Tricho_rsr[w2] * cellChl;
      }
      aawave(rad, absorbance, aA_s_Tricho, e->num_rsr_waves,e->rsr_waves);

      // free intermediate var
      d_free_1d(absorbance);
    }
    // double acdom443 = e->bio_opt->acdom443star * DOC;
    for (w2=0; w2<e->num_rsr_waves; w2++){

      /* need to calculate total absorption and backscatter at reflectance wavelengths */

      at_r[w2] = bio->kw_s2[w2];

      if (ws->domain == 'G'|| ws->domain == 'H'){
	//	if (y[ws->salt_i] < 34.145){
	//  at_r[w2] += (1.2336 + (-0.0332 * y[ws->salt_i])) * exp(-(ws->S_CDOM) * (e->rsr_waves[w2]-443.0)); 
	//}
	//at_r[w2] += 0.01 * exp(-(ws->S_CDOM) * (e->rsr_waves[w2]-443.0));

	if (ws->cdom_gbr_i > -1){
	  // if (y[ws->salt_i] < 36.855){
	    //  double a443 = (1.2336 + (-0.0332 * y[ws->salt_i])) + 0.01; // due to clear water and cdom only.
	    double a443 = (1.2336 + (-0.0332 * (1.0-y[ws->cdom_gbr_i])*36.855));// + 0.01; // 0.01 open-ocean a_cdom(443)
	    double S_cdom = 0.0061 * pow(a443,-0.309);  // Blondeau-Patissier 2009 JGR: 114: C05003.
	    at_r[w2] += a443 * exp(- S_cdom * (e->rsr_waves[w2]-443.0));
	    //}else{
	    // at_r[w2] += 0.01 * exp(-0.010 * (e->rsr_waves[w2]-443.0));
	    //}
	}else{
	  if (y[ws->salt_i] < 34.145){   // 36.855 for 0.01, 34.145 for 0.1
	    double a443 = (1.2336 + (-0.0332 * y[ws->salt_i])) + 0.01; // due to clear water and cdom only.
	    double S_cdom = 0.0061 * pow(a443,-0.309);  // Blondeau-Patissier 2009 JGR: 114: C05003.
	    at_r[w2] += a443 * exp(- S_cdom * (e->rsr_waves[w2]-443.0));
	  }else{
	    at_r[w2] += 0.01 * exp(-0.010 * (e->rsr_waves[w2]-443.0));
	  }
	}
      }
      if (ws->domain == 'C'){
	// From cruises in Oct 2017, Mar 2018
	// at_r[w2] += (0.2875 - 0.0059 * min(y[ws->salt_i],35.0)) * exp(-ws->S_CDOM * (e->rsr_waves[w2]-440.0));
	at_r[w2] += (0.2875 - 0.0059 * min(y[ws->salt_i],35.0)) 
	  * exp(-(ws->S_CDOM - max(y[ws->salt_i]-30.0,0.0)*0.007/5.0 ) * (e->rsr_waves[w2]-440.0));
      }
      if (ws->domain == 'T'){
	at_r[w2] += 0.01 * exp(-0.014 * (e->rsr_waves[w2]-440.0));
	at_r[w2] += ws->a440cdom_pale * min(y[ws->cdom_pale_i],1.0) * exp(- ws->Scdom_pale * (e->rsr_waves[w2]-440.0));
	at_r[w2] += ws->a440cdom_amber * min(y[ws->cdom_amber_i],1.0) * exp(- ws->Scdom_amber * (e->rsr_waves[w2]-440.0));
	at_r[w2] += ws->a440cdom_dark * min(y[ws->cdom_dark_i],1.0) * exp(- ws->Scdom_dark * (e->rsr_waves[w2]-440.0));
      }

      if (ws->domain == 'M'){
	at_r[w2] += 0.01 * exp(-0.014 * (e->rsr_waves[w2]-440.0));
	at_r[w2] += ws->a440cdom_dark * min(y[ws->cdom_dark_i],1.0) * exp(- ws->Scdom_dark * (e->rsr_waves[w2]-440.0));
      }

      bt_r[w2] = bio->bw_s2[w2];  // clear water.
      bb_r[w2] = ws->bbw_ratio * bio->bw_s2[w2];  // clear water.

      switch (ws->domain){     /* STILL NEED TO ADD TOTAL SCATTERING FOR 'G' OPTION */
      case 'G' : 
	at_r[w2] += bio->ad_A2[w2] * pow(NAP, bio->ad_B2[w2]);
	at_r[w2] += bio->ad_CarbSand_A2[w2] * pow(CarbSand, bio->ad_CarbSand_B2[w2]);

	bb_r[w2] += ws->bbnap_ratio * bio->bbp_A2[w2] * pow(NAP, bio->bbp_B2[w2]); /* NAP backscattering */
	bb_r[w2] += ws->bbnap_ratio * bio->bbp_CarbSand_A2[w2] * pow(CarbSand, bio->bbp_CarbSand_B2[w2]); /* NAP backscattering */
      
	break;
      case 'H' :

        // at_r[w2] += bio->aC_CAL1_rsr[w2] * Mud_carbonate * 1000.0 / 4.0;
	// bt_r[w2] += bio->bC_CAL1_rsr[w2] * Mud_carbonate * 1000.0 / 4.0; /* NAP backscattering */
        // bb_r[w2] += bio->B_carb_rsr[w2] * bio->bC_CAL1_rsr[w2] * Mud_carbonate * 1000.0 / 4.0; /* NAP backscattering */

	at_r[w2] += bio->LJCOc_aNAP_rsr[w2] * Mud_carbonate;
	bt_r[w2] += bio->LJCOc_bNAP_rsr[w2] * Mud_carbonate; /* NAP backscattering */
        bb_r[w2] += bio->B_carb_rsr[w2] * bio->LJCOc_bNAP_rsr[w2] * Mud_carbonate; /* NAP backscattering */

	// at_r[w2] += bio->ad_CarbSand_A2[w2] * pow(CarbSand, bio->ad_CarbSand_B2[w2]);
	// bb_r[w2] += ws->bbnap_ratio * bio->bbp_CarbSand_A2[w2] * pow(CarbSand, bio->bbp_CarbSand_B2[w2]); /* NAP backscattering */
	// bt_r[w2] += bio->bbp_CarbSand_A2[w2] * pow(CarbSand, bio->bbp_CarbSand_B2[w2]); /* NAP backscattering */

	//at_r[w2] += bio->aC_KAO1[w2] * (Mud_mineral+FineSed+NAP_noEFI) * 1000.0;
	//* at_r[w2] += bio->ad_A2[w2] * (Dust+Mud_mineral+FineSed+NAP_noEFI) ;
	
        at_r[w2] += bio->LJCO_aNAP_rsr[w2] * (Dust+Mud_mineral+FineSed+NAP_noEFI);
	bt_r[w2] += bio->LJCO_bNAP_rsr[w2] * (Dust + Mud_mineral+FineSed+NAP_noEFI);
	bb_r[w2] += bio->B_terr_rsr[w2] * bio->LJCO_bNAP_rsr[w2] * (Dust + Mud_mineral+FineSed+NAP_noEFI);

	//bb_r[w2] += bio->B_terr_rsr[w2] * bio->bC_KAO1[w2] * (Mud_mineral+FineSed+NAP_noEFI) * 1000.0; /* NAP backscattering */
	//* bt_r[w2] += bio->bbp_A2[w2] * (Dust + Mud_mineral+FineSed+NAP_noEFI);
	//* bb_r[w2] += bio->B_terr_rsr[w2] * bio->bbp_A2[w2] * (Dust + Mud_mineral+FineSed+NAP_noEFI);
	//20nov bt_r[w2] += bio->LJCO_bNAP_rsr[w2] * (Dust + Mud_mineral+FineSed+NAP_noEFI);
	//20nov bb_r[w2] += bio->B_terr_rsr[w2] * bio->LJCO_bNAP_rsr[w2] * (Dust + Mud_mineral+FineSed+NAP_noEFI);
	
	break;
      case 'C' :
	at_r[w2] += bio->aC_ICE2_rsr[w2] * NAP * 1.0e3;
	bt_r[w2] += bio->bC_ICE2_rsr[w2] * NAP * 1.0e3;
	bb_r[w2] += bio->B_terr_rsr[w2] * bio->bC_ICE2_rsr[w2] * NAP * 1.0e3;
	break;
      case 'T' :
	at_r[w2] += bio->aC_ICE2_rsr[w2] * NAP * 1.0e3;
	bt_r[w2] += bio->bC_ICE2_rsr[w2] * NAP * 1.0e3;
	bb_r[w2] += bio->B_terr_rsr[w2] * bio->bC_ICE2_rsr[w2] * NAP * 1.0e3;
	break;
      case 'M' :
	at_r[w2] += bio->aC_ICE2_rsr[w2] * NAP * 1.0e3;
	bt_r[w2] += bio->bC_ICE2_rsr[w2] * NAP * 1.0e3;
	bb_r[w2] += bio->B_terr_rsr[w2] * bio->bC_ICE2_rsr[w2] * NAP * 1.0e3;
       break;
     }

      if (ws->KI_l_i > -1){
	at_r[w2] += cellnum_l * aA_s_l[w2];
	//bt_r[w2] += ws->bbcell1 * pow(ws->rad_l*2e6, ws->bbcell2)*cellnum_l;
	bb_r[w2] += bio->scatfrac[w2] * ws->bbcell1 * pow(ws->rad_l*2e6, ws->bbcell2)*cellnum_l;
	bt_r[w2] += ws->bphy * y[ws->PhyL_N_i] / ws->NtoCHL; 
      }
      if (ws->KI_s_i > -1){
	at_r[w2] += cellnum_s * aA_s_s[w2];
	//bt_r[w2] += ws->bbcell1 * pow(ws->rad_s*2e6, ws->bbcell2)*cellnum_s;
	bb_r[w2] += bio->scatfrac[w2] * ws->bbcell1 * pow(ws->rad_s*2e6, ws->bbcell2)*cellnum_s;
        bt_r[w2] += ws->bphy * y[ws->PhyS_N_i] / ws->NtoCHL ;
      }
      if (ws->KI_MPB_i > -1){
	at_r[w2] +=  cellnum_MPB * aA_s_MPB[w2];
	//bt_r[w2] += ws->bbcell1 * pow(ws->rad_MPB*2e6, ws->bbcell2)*cellnum_MPB;
	bb_r[w2] += bio->scatfrac[w2] * ws->bbcell1 * pow(ws->rad_MPB*2e6, ws->bbcell2)*cellnum_MPB;
	bt_r[w2] += ws->bphy * y[ws->MPB_N_i] / ws->NtoCHL ;
      }
      if (ws->KI_PhyD_i > -1){
	at_r[w2] += cellnum_PhyD * aA_s_PhyD[w2];
	//bt_r[w2] += ws->bbcell1 * pow(ws->rad_PhyD*2e6, ws->bbcell2)*cellnum_PhyD;
	bb_r[w2] += bio->scatfrac[w2] * ws->bbcell1 * pow(ws->rad_PhyD*2e6, ws->bbcell2)*cellnum_PhyD;
	bt_r[w2] += ws->bphy * y[ws->PhyD_N_i] / ws->NtoCHL ;
      }
      if (ws->KI_Tricho_i > -1){
	at_r[w2] += cellnum_Tricho * aA_s_Tricho[w2];
	//bt_r[w2] += ws->bbcell1 * pow(ws->rad_Tricho*2e6, ws->bbcell2)*cellnum_Tricho;
	bb_r[w2] += bio->scatfrac[w2] * ws->bbcell1 * pow(ws->rad_Tricho*2e6, ws->bbcell2)*cellnum_Tricho;
	bt_r[w2] += ws->bphy * y[ws->Tricho_N_i] / ws->NtoCHL ;
      }

      /* output backscatter at 700 nm for glider comparison */

      if (w2 == ws->w700){
	if (ws->bb_700_i > -1)
	  y[ws->bb_700_i] = bb_r[w2];
      }

      /* add Hazen true colour scale Palin (1995) Water and Water Engineering 59, 341-345. ? */ 
      
      /* approx. vertical attentuation coefficient */
      
      adlen = sqrt(1.0 + (gone*costhetaw-gtwo) * bt_r[w2] / at_r[w2]) / costhetaw;
      
      Kd_r = at_r[w2] * adlen;

      /* Now add sources of light to reflectance calculations. The source must be spread across a waveband, which is sensor dependednt */

      double nFLHadd = 0.0;
      double BLadd = 0.0;
      double BLP = 0.0;
      double BLstim = 1.0;
      double iso2upwell = 0.4135669939329333; // 0.5 * 1/(pi*(3/4/pi)^(2/3)) m3/m2

      // I guess we should use cell centred downwelling irradiance, but we haven't calculated it.
      
      if (w2 == ws->w678){ // wavelength of nFLH, wPAR_top is 690.

	/* remember light in is over 20 nm, which is about right */ 
	
	double light678 = lighttop_s[ws->w670] * (12.0/20.0) + lighttop_s[ws->wPAR_top] * (8.0/20.0);
	
	if (light678 > 1.0e-3){  // avoid divide by zero.
	  nFLHadd = nFLH * iso2upwell / (light678 * 8.359335857479461e-09 * 678.0);  // unitless
	}
      }      

      // Add bioluminescence from large phytoplankton at 470 (actually 472 +/- 35 nm for noctiluca Mobely).

      // BLstim is a placeholder to include the effect of stimuli and physiology. At the moment all cells luminescencing at maximum.

      if (w2 == ws->w470){
	
	// BLP = BLstim * cellnum_l * 3.0e10 * (6.63e-34 *3.0e8) / (470.0e-9);  // W m-3 - Mobley 94
	BLP = BLstim * (y[ws->PhyL_Chl_i] / 1000.0) * 1.1 * 6.63e-34 * 3.0e8 *6.022e23 / 470.0e-9 / 86400.0;  // W m-3 - Ondercin 1995
	// printf("ws->BLP_i %d, cellnum_l %e, BLadd %e \n",ws->BLP_i,cellnum_l,BLadd);
	// printf("Chl way %e, cell way %e\n",y[ws->PhyL_Chl_i]*1.1/1000.0,cellnum_l * 3.0e10);
        if (ws->BLP_i > -1){
	  y[ws->BLP_i] = BLP * 1.0e3;  // Ouput as uW m-3
	}
      
	BLadd = iso2upwell * BLP / max(lighttop_s[ws->w470],0.015/10.0); // minimum is 10 percent full moonlight
	// printf("iso2upwell * BLP %e, max(lighttop_s[ws->w470],0.015/10.0) %e, u %e \n",iso2upwell * BLP,lighttop_s[ws->w470],bb_r[w2] / (at_r[w2] + bb_r[w2]));
      }
      
      // w_bot is what is remaining //
      
      top = w_bot[w2];
      w_bot[w2] = top * exp(-2.0 * Kd_r * dz);
      u_surf[w2] += (nFLHadd + BLadd + (bb_r[w2] / (at_r[w2] + bb_r[w2]))) * (top - w_bot[w2]);
    }

    /* do Secchi depth calculation */

    if (ws->w488 > -1){

      /* calculate secchi depth using 488 nm */
      if (isnan(E_secchi_cv[0])) {
	E_secchi_cv[0] = 1.0;    /* Normalised light at top equal 1. */
	z_secchi_cv[0] = 0.0;
      }
      double z_secchi = z_secchi_cv[0];
      double E_secchi = E_secchi_cv[0];
      if (E_secchi > exp(-1.0)){
	adlen = sqrt(1.0 + (gone*costhetaw-gtwo) * bt_r[ws->w488] / at_r[ws->w488]) / costhetaw;
	E_secchi_cv[0] = E_secchi * exp(-at_r[ws->w488] * adlen * dz);
	if (E_secchi_cv[0] < exp(-1.0)){ // this is the step the disk dissappears from view
	  z_secchi_cv[0] = z_secchi + log(exp(-1.0)/E_secchi) / (- at_r[ws->w488] * adlen);
	}else{                           // otherwise still in view.
	  z_secchi_cv[0] = z_secchi + dz;
	}
      } 
    }

    if (ws->Turbidity_i > -1 && ws->w590 > -1){ /* backscatter - clear water backscatter */
      y[ws->Turbidity_i] = 47.02 * (bb_r[ws->w590] - ws->bbw_ratio * bio->bw_s2[ws->w590]) + 0.13;
    }
    
    d_free_1d(at_r);
    d_free_1d(bt_r);
    d_free_1d(bb_r);

    if (aA_s_l != NULL)
      d_free_1d(aA_s_l);
    if (aA_s_s != NULL)
      d_free_1d(aA_s_s);
    if (aA_s_MPB != NULL)
      d_free_1d(aA_s_MPB);
    if (aA_s_PhyD != NULL)
      d_free_1d(aA_s_PhyD);
    if (aA_s_Tricho != NULL)
      d_free_1d(aA_s_Tricho);
  }
  
 }
 
