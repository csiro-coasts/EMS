/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/ecology/ecology_tracer_defaults.c
 *  
 *  Description:
 *  Parameter defaults
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: ecology_tracer_defaults.c 6721 2021-03-29 04:38:30Z bai155 $
 *
 */

#include "hd.h"

/* This is defined in model/hd/ecology.c */
extern const char *ECONAME3D[][2];
extern const char *ECONAME2D[][2];
extern const int NUM_ECO_VARS_3D;
extern const int NUM_ECO_VARS_2D;

/* from ecofunct.c */
extern double PhyCellMass(double r);

/* from sediments.c */
extern void sed_set_tr_att(tracer_info_t *tr, char *key, void *value);


#if defined(HAVE_ECOLOGY_MODULE)
#include "ecology_tracer_defaults.h"
/*-------------------------------------------------------------------*/
/* Sets the ecology tracers standard default attributes              */
/*-------------------------------------------------------------------*/
void eco_defaults_std(tracer_info_t *tracer, char *trname, ecology *e)
{

/* #define SEDIM         0x0001 */
/* #define INTER         0x0002 */
/* #define WATER         0x0004 */
/* #define HYDRO         0x0008 */
/* #define SEDIMENT      0x0010 */
/* #define ECOLOGY       0x0020 */
/* #define WAVE          0x0040 */
/* #define TRACERSTAT    0x0080 */
/* #define PROGNOSTIC    0x0100 */
/* #define DIAGNOSTIC    0x0200 */
/* #define PARAMETER     0x0400 */
/* #define FORCING       0x0800 */

  /* Initial conditions that are independent of parameters */

  double Mic_N = 7.0;
  double Mic_NR = Mic_N/2.0;
  double Mic_PR = Mic_NR/16.0*32.0/14.0;
  double Mic_Chl = Mic_N/7.0; 
  double Mic_I = Mic_NR/14.0*1060.0/16.0;
  
  double SG_N;
  double SGROOT_N;
  double SGH_N;
  double SGHROOT_N;
  double SGP_N;
  double SGPROOT_N;
  double SGD_N;
  double SGDROOT_N;
  double MA_N;
  double CH_N;
  double CS_N;
  double CS_Chl;
  double Mic_N_sed;
  double Mic_NR_sed;
  double Mic_PR_sed;
  double Mic_Chl_sed; 
  double Mic_I_sed;
  double Tricho_N;
  double Tricho_NR;
  double Tricho_PR;
  double Tricho_Chl; 
  double Tricho_I;
  double ZL_N;
  double ZS_N;

  double DetR_C_sed;
  double DetR_N_sed;
  double DetR_P_sed;

  double DOR_C_sed;
  double DOR_N_sed;
  double DOR_P_sed;

  double DetPL_N_sed;
  double DetBL_N_sed;

  double cover;
  double grow;
  double CSm;
  double CSrad;
  double CHpolypden;
  
/* Initial conditions that depend on the parameters chosen */

  if (e){

    /* Benthic plants with 0.2 x 63 % cover */

    eco_write_setup(e,"Reading parameter list to calculated default tracer values \n");
    
    cover = 0.2;

    ////////// Need to avoid if no seagrass processes. ////////////////////

    double SGleafden = try_parameter_value(e, "SGleafden");
    if (isnan(SGleafden)){
      SGleafden = 1.0;
    }
    double SGfrac = try_parameter_value(e, "SGfrac");
    if (isnan(SGfrac)){
      SGfrac = 0.5;
    }
    SG_N = cover / SGleafden;
    SGROOT_N  = SG_N * SGfrac;

    double SGHleafden = try_parameter_value(e, "SGHleafden");
    if (isnan(SGHleafden)){
      SGHleafden = 1.0;
    }
    double SGHfrac = try_parameter_value(e, "SGHfrac");
    if (isnan(SGHfrac)){
      SGHfrac = 0.5;
    }
    SGH_N = cover / SGHleafden;
    SGHROOT_N  = SGH_N * SGHfrac;

    double SGDleafden = try_parameter_value(e, "SGDleafden");
    if (isnan(SGDleafden)){
      SGDleafden = 1.0;
    }
    double SGDfrac = try_parameter_value(e, "SGDfrac");
    if (isnan(SGDfrac)){
      SGDfrac = 0.5;
    }
    SGD_N = cover / SGDleafden;
    SGDROOT_N  = SGD_N * SGDfrac;

    double SGPleafden = try_parameter_value(e, "SGPleafden");
    if (isnan(SGPleafden)){
      SGPleafden = 1.0;
    }
    double SGPfrac = try_parameter_value(e, "SGPfrac");
    if (isnan(SGPfrac)){
      SGPfrac = 0.5;
    }
    SGP_N = cover / SGPleafden;
    SGPROOT_N  = SGP_N * SGPfrac;

    //  Macroalgae

    double MAleafden = try_parameter_value(e, "MAleafden");
    if (isnan(MAleafden)){
      MAleafden = 1.0;
    }
    MA_N = cover / MAleafden;

    /* make symbionts cover 0.1 of the surface area of the corals */

    // cover = 0.2;
    // CSrad = get_parameter_value(e,"CSrad");
    // CSm = PhyCellMass(CSrad);
    // CHpolypden = get_parameter_value(e,"CHpolypden");
    // CH_N = cover / get_parameter_value(e,"CHpolypden");

    // CH_N = cover / CHpolypden;
    // CS_N = 0.1 * (CH_N * CHpolypden) / ((CSm * 16.0 * 1000.0 * 14.01) * CSrad * CSrad * M_PI);

    // Use steady-state values from eReefs
    
    CH_N = 0.15;    
    CS_N = 5.0;

    //    0.5 = ((CS_N / (CSm * red_A_N * 1000.0 * MW_Nitr)) * ws->CSrad * ws->CSrad * M_PI) / (CH_N * ws->CHpolypden*2.0);
    CS_Chl = CS_N/7.0;

    // printf("CH_N %e, CS_N %e, CS_Chl %e \n",CH_N,CS_N,CS_Chl);

    /* Plankton with steady state values at 10 % growth rate */

    grow = 0.1;

    double MBumax = try_parameter_value(e, "MBumax");
    if (isnan(MBumax)){
      MBumax = 0.5;
    }

    double MPB_mQ = try_parameter_value(e, "MPB_mQ");
    if (isnan(MPB_mQ)){
      MPB_mQ = 0.5;
    }

    Mic_N_sed = grow * MBumax / MPB_mQ;
    Mic_NR_sed = Mic_N_sed/2.0;
    Mic_PR_sed = Mic_NR_sed/16.0*32.0/14.0;
    Mic_Chl_sed = Mic_N_sed/7.0; 
    Mic_I_sed = Mic_NR_sed/14.0*1060.0/16.0;

    double Tricho_umax = try_parameter_value(e, "Tricho_umax");
    if (isnan(Tricho_umax)){
      Tricho_umax = 0.5;
    }

    double Tricho_mQ = try_parameter_value(e, "Tricho_mQ");
    if (isnan(Tricho_mQ)){
      Tricho_mQ = 0.5;
    }
    
    Tricho_N = grow * Tricho_umax / Tricho_mQ;
    Tricho_NR = Tricho_N/2.0;
    Tricho_PR = Tricho_NR/16.0*32.0/14.0;
    Tricho_Chl = Tricho_N/7.0; 
    Tricho_I = Tricho_NR/14.0*1060.0/16.0;

    double ZLumax = try_parameter_value(e, "ZLumax");
    if (isnan(ZLumax)){
      ZLumax = 0.5;
    }

    double ZL_mQ = try_parameter_value(e, "ZL_mQ");
    if (isnan(ZL_mQ)){
      ZL_mQ = 0.5;
    }

    double ZSumax = try_parameter_value(e, "ZSumax");
    if (isnan(ZSumax)){
      ZSumax = 0.5;
    }

    double ZS_mQ = try_parameter_value(e, "ZS_mQ");
    if (isnan(ZS_mQ)){
      ZS_mQ = 0.5;
    }

    ZL_N = grow * ZLumax / ZL_mQ;
    ZS_N = grow * ZSumax / ZS_mQ;

    /* Sediment values defaulted to GBR4 - all water column values 0.01 */

    DetR_C_sed = 8000.0;
    DetR_N_sed = 500.0;
    DetR_P_sed = 1100.0;

    DOR_C_sed = 500.0;
    DOR_N_sed = 500.0;
    DOR_P_sed = 500.0;

    DetPL_N_sed = 500.0;
    DetBL_N_sed = 500.0; 

  }else{
    //  e_quit("eco_defaults_std: Failed to read ecology paraneters when setting ecology default tracer initial conditions");
  }
  
  /* 5 - WATER | SEDIMENT ; 2 - BENTHIC */
  
  /* diagn, advect, diff, type, diss, part, in water, in sed, flag, obc, type */
  
  // PROGNOSTIC means model must have them if process asking for them is invoked (i.e. pH).

  /* d - 1: for diagnostic fluxes so divided by time step.
     d - 2: for diagnostic states. */
  
  /* flags for last column (f) ECO_NORESET is 1 */

  /* All tracers */
  eco_def_t eco_def[] = {
    /* name         unit      fillw, fills      d  a  d  t  d  p  iw is f  obc     extended type */ 
    {"Age",       "d",        0.0,    0.0,      0, 1, 1, 5, 1, 0, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"source",    "-",        0.0,    0.0,      2, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"PhyL_N",    "mg N m-3", Mic_N,  0.0,      0, 1, 1, 5, 0, 1, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"PhyS_N",    "mg N m-3", Mic_N,  0.0,      0, 1, 1, 5, 0, 1, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"MPB_N",     "mg N m-3", Mic_N, Mic_N_sed, 0, 1, 1, 5, 0, 1, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"ZooL_N",    "mg N m-3", ZL_N,  0.0,       0, 1, 1, 5, 0, 1, 1, 0, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"ZooS_N",    "mg N m-3", ZS_N,  0.0,       0, 1, 1, 5, 0, 1, 1, 0, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"PhyD_N",    "mg N m-3", Mic_N,  0.0,      0, 1, 1, 5, 0, 1, 1, 0, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"PhyD_C",    "mg C m-3", Mic_N*6.625, 0.0, 0, 1, 1, 5, 0, 1, 1, 0, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"Phy_L_N2",  "mg N m-3", 0.01,  0.01,      0, 1, 1, 5, 0, 1, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"NH4",       "mg N m-3", 0.2,   200.0,     0, 1, 1, 5, 1, 0, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"NO3",       "mg N m-3", 0.1,   500.0,     0, 1, 1, 5, 1, 0, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"DIP",       "mg P m-3", 0.5,   100.0,     0, 1, 1, 5, 1, 0, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"PIP",       "mg P m-3", 0.1,   100.0,     0, 1, 1, 5, 0, 1, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"PIP_Dust",  "mg P m-3", 0.1,   100.0,     0, 1, 1, 5, 0, 1, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"PIPF",      "mg P m-3", 0.0,   100.0,     0, 1, 1, 5, 0, 1, 1, 1, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"PIPI",      "mg P m-3", 0.01,   0.01,     0, 1, 1, 5, 0, 1, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"DIC",       "mg C m-3", 24758.61,24758.61,0, 1, 1, 5, 1, 0, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"DOR_C",     "mg C m-3", 0.01, DOR_C_sed,  0, 1, 1, 5, 1, 0, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"DOR_N",     "mg N m-3", 0.01, DOR_N_sed,  0, 1, 1, 5, 1, 0, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"DOR_P",     "mg P m-3", 0.01, DOR_P_sed,  0, 1, 1, 5, 1, 0, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"DetR_C",    "mg C m-3", 0.01, DetR_C_sed, 0, 1, 1, 5, 0, 1, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"DetR_N",    "mg N m-3", 0.01, DetR_N_sed, 0, 1, 1, 5, 0, 1, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"DetR_P",    "mg P m-3", 0.01, DetR_P_sed, 0, 1, 1, 5, 0, 1, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"DetPL_N",   "mg N m-3", 0.01, DetPL_N_sed,0, 1, 1, 5, 0, 1, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"DetBL_N",   "mg N m-3", 0.0,  DetBL_N_sed,0, 1, 1, 5, 0, 1, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"Oxygen",    "mg O m-3", 6505.0,6505.0,    0, 1, 1, 5, 1, 0, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"COD",       "mg O m-3",     0.0,   0.0,   0, 1, 1, 5, 1, 0, 1, 1, 0, STATIS, ECOLOGY|PROGNOSTIC},
    {"Oxy_sat",   "%",        100.0, 100.0,     2, 0, 0, 5, 0, 0, 1, 1, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"Light",     "W m-2",    0.0,   0.0,       2, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"PAR",     "mol photon m-2 s-1", 0.0,0.0,  2, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"PAR_z","mol photon m-2 s-1", 0.0,0.0,  2, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"K_heat",     "m-1",    0.0,   0.0,        2, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"Kd",        "m-1",      0.0,   0.0,       2, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"Epilight",  "W m-2",    0.0,   0.0,       2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"EpiPAR",  "mol photon m-2 s-1", 0.0, 0.0, 2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"EpiPAR_sg",  "mol photon m-2 d-1", 0.0, 0.0, 2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"Epilightatt","m-1 ",       0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"TN",        "mg N m-3", 0.0,   0.0,       2, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"TP",        "mg P m-3", 0.0,   0.0,       2, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"TC",        "mg C m-3", 0.0,   0.0,       2, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"DIN",       "mg N m-3", 0.0,   0.0,       2, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"Chl_a",     "mg Chl m-3",0.0,  0.0,       2, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"Chl_a_sum", "mg Chl m-3",0.0,  0.0,       2, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"EFI",       "kg m-3",   0.0,   0.0,       2, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"PhyL_N_pr", "mg C m-3 d-1", 0.0,   0.0,   1, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"PhyS_N_pr", "mg C m-3 d-1", 0.0,   0.0,   1, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"ZooL_N_rm", "mg C m-3 d-1", 0.0,   0.0,   1, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"ZooS_N_rm", "mg C m-3 d-1", 0.0,   0.0,   1, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"Den_fl",    "mg N m-2 d-1 per layer", 0.0,   0.0,   1, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"Den_eff",   " ",            0.0,   0.0,   2, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"NH4_pr",    "mg N m-3 d-1", 0.0,   0.0,   1, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"MPB_N_pr",  "mg C m-3 d-1", 0.0,   0.0,   1, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"PhyL_N_gr", "d-1",          0.0,   0.0,   1, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"PhyS_N_gr", "d-1",          0.0,   0.0,   1, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"BOD",       "mg O m-3",     0.0,   0.0,   2, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"EpiBOD",    "mg O m-2",     0.0,   0.0,   2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"MPB_N_gr",  "d-1",          0.0,   0.0,   1, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"ZooL_N_gr", "d-1",          0.0,   0.0,   1, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"ZooS_N_gr", "d-1",          0.0,   0.0,   1, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"PhyD_N_gr", "d-1",          0.0,   0.0,   1, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"PhyD_N_pr", "mg C m-3 d-1", 0.0,   0.0,   1, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"Oxy_pr",    "mg O m-3 d-1 ",0.0,   0.0,   1, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PARAMETER},
    {"Nfix",      "mg N m-3 s-1", 0.0,   0.0,   1, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PARAMETER},
    {"Amm_fl",    "mg N m-3 s-1", 0.0,   0.0,   1, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PARAMETER},
    {"Phy_L_N2_fix","mg N m-3 s-1",0.0,  0.0,   1, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PARAMETER},
    {"SG_N",     "g N m-2",  SG_N,      0.0,    0, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"SGH_N",     "g N m-2",  SGH_N,    0.0,    0, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"SGROOT_N",  "g N m-2",  SGROOT_N, 0.0,    0, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"SGHROOT_N", "g N m-2",  SGHROOT_N,0.0,    0, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"SGP_N",     "g N m-2",  SGP_N,    0.0,    0, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"SGD_N",     "g N m-2",  SGD_N,    0.0,    0, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"SGPROOT_N", "g N m-2", SGPROOT_N, 0.0,    0, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"SGDROOT_N", "g N m-2", SGDROOT_N,0.0,    0, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"MA_N",      "g N m-2",  MA_N,     0.0,    0, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"MAG_N",      "g N m-2",  MA_N,     0.0,    0, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"MAR_N",      "g N m-2",  MA_N,     0.0,    0, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"CS_N",      "mg N m-2",  CS_N,    0.0,    0, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC}, 
    {"CS_NR",      "mg N m-2",  CS_N/2.0,    0.0,    0, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"CS_PR",      "mg P m-2",  CS_N/16.0*32.0/14.0,    0.0,    0, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"CS_I",      "mmol photon m-2", CS_N/14.0*1060.0/16.0,    0.0,    0, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"CS_Chl",    "mg Chl m-2",CS_Chl,  0.0,    0, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"CS_Xh",    "mg Xan m-2",CS_Chl/4.0,  0.0,    0, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"CS_Xp",    "mg Xan m-2",CS_Chl/4.0,  0.0,    0, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"CS_Qred",   "umol RCII m-2",CS_Chl*0.002/893.49/3.0,  0.0,    0, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"CS_Qox",    "umol RCII m-2",CS_Chl*0.002/893.49/3.0,  0.0,    0, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"CS_Qi",    "umol RCII m-2",CS_Chl*0.002/893.49/3.0,  0.0,    0, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"CS_RO",    "umol ROS m-2",0.0*CS_Chl/3.0,  0.0,    0, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"CH_N",      "g N m-2",  CH_N,     0.0,    0, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"EpiTN",     "mg N m-2", 0.0,   0.0,       2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"EpiTP",     "mg P m-2", 0.0,   0.0,       2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"EpiTC",     "mg C m-2", 0.0,   0.0,       2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"SG_N_pr",   "g N m-2 d-1", 0.0,   0.0,   1, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"SGH_N_pr",  "g N m-2 d-1", 0.0,   0.0,   1, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"SGP_N_pr",   "g N m-2 d-1", 0.0,   0.0,   1, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"SGD_N_pr",  "g N m-2 d-1", 0.0,   0.0,   1, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"MA_N_pr",   "g N m-2 d-1", 0.0,   0.0,   1, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"MAG_N_pr",   "g N m-2 d-1", 0.0,   0.0,   1, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"MAR_N_pr",   "g N m-2 d-1", 0.0,   0.0,   1, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"SG_N_gr",   "d-1",0.0,0.0,                1, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"SGH_N_gr",  "d-1",0.0,0.0,                1, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"SGP_N_gr",  "d-1",0.0,0.0,                1, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"SGD_N_gr",  "d-1",0.0,0.0,                1, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"MA_N_gr",   "d-1",0.0,0.0,                1, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"MAG_N_gr",   "d-1",0.0,0.0,                1, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"MAR_N_gr",   "d-1",0.0,0.0,                1, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"EpiOxy_pr", "mg O m-2 d-1", 0.0,   0.0,   1, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"Tricho_N",  "mg N m-3", Tricho_N,  0.0,   0, 1, 1, 5, 0, 1, 1, 0, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"Tricho_N_gr","d-1",        0.0,   0.0,    1, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PARAMETER},
    {"Tricho_N_pr","mg C m-3 d-1",0.0,   0.0,   1, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PARAMETER},
    {"Zenith",    "rad",          0.0,   0.0,   2, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PARAMETER},
    {"PhyL_I",   "mmol photon m-3", Mic_I,0.0,  0, 1, 1, 5, 0, 1, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"PhyL_NR",  "mg N m-3",     Mic_NR,  0.0,  0, 1, 1, 5, 0, 1, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"PhyL_PR",  "mg P m-3",     Mic_PR,  0.0,  0, 1, 1, 5, 0, 1, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"PhyL_Chl", "mg Chl m-3 ",  Mic_Chl, 0.0,  0, 1, 1, 5, 0, 1, 0, 0, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"PhyL_sv",  "m s-1 ",    -0.0001157, 0.0,  2, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"PhyS_I",   "mmol photon m-3",Mic_I, 0.0,  0, 1, 1, 5, 0, 1, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"PhyS_NR",   "mg N m-3",    Mic_NR,  0.0,  0, 1, 1, 5, 0, 1, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"PhyS_PR",   "mg P m-3",    Mic_PR,  0.0,  0, 1, 1, 5, 0, 1, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"PhyS_Chl",  "mg Chl m-3 ", Mic_Chl, 0.0,  0, 1, 1, 5, 0, 1, 0, 0, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"MPB_NR",    "mg N m-3",Mic_NR,Mic_NR_sed, 0, 1, 1, 5, 0, 1, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"MPB_PR",    "mg P m-3",Mic_PR,Mic_PR_sed, 0, 1, 1, 5, 0, 1, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"MPB_I", "mmol photon m-3",Mic_I,Mic_I_sed,0, 1, 1, 5, 0, 1, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"MPB_Chl","mg Chl m-3 ",Mic_Chl,Mic_Chl_sed,0, 1, 1, 5, 0, 1, 0, 0, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"PhyD_NR",   "mg N m-3",    Mic_NR,  0.0,  0, 1, 1, 5, 0, 1, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"PhyD_PR",   "mg P m-3",    Mic_PR,  0.0,  0, 1, 1, 5, 0, 1, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"PhyD_Chl",  "mg Chl m-3",  Mic_Chl, 0.0,  0, 1, 1, 5, 0, 1, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"PhyD_I",   "mmol photon m-3", Mic_I,0.0,  0, 1, 1, 5, 0, 1, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"Tricho_NR", "mg N m-3",    Mic_NR,  0.0,  0, 1, 1, 5, 0, 1, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"Tricho_PR", "mg P m-3",    Mic_PR,  0.0,  0, 1, 1, 5, 0, 1, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"Tricho_I","mmol photon m-3", Mic_I, 0.0,  0, 1, 1, 5, 0, 1, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"Tricho_Chl","mg Chl m-3",  Mic_Chl, 0.0,  0, 1, 1, 5, 0, 1, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"Tricho_sv", "m s-1"      ,   0.003, 0.0,  2, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"ZooL_sv",   "m s-1"      ,   0.003, 0.0,  2, 0, 0, 5, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"PH",           "log(mM)", 8.0,   8.2,     0, 0, 0, 5, 0, 0, 1, 1, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"dco2star",  "mol m-3",      0.0,   0.0,   2, 0, 0, 5, 0, 0, 1, 1, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"CO2_starair","not sure",        0.0, 0.0, 2, 0, 0, 5, 0, 0, 1, 1, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"CO2_star","mol m-3", 0.0, 0.0,            2, 0, 0, 5, 0, 0, 1, 1, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"pco2surf","ppmv", 0.0, 0.0,               2, 0, 0, 5, 0, 0, 1, 1, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"dpCO2","ppmv", 0.0, 0.0,                  2, 0, 0, 5, 0, 0, 1, 1, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"CO32",      "mmol m-3", 262.0, 262.0,     2, 0, 0, 5, 0, 0, 1, 1, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"HCO3",      "mmol m-3",1650.0,1650.0,     2, 0, 0, 5, 0, 0, 1, 1, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"CO2_flux",  "mg C m-2 s-1", 0,     0,     2, 0, 0, 5, 0, 0, 1, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"O2_flux",  "mg O m-2 s-1", 0,     0,      1, 0, 0, 5, 0, 0, 1, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"at_440",    "m-1",        0,     0,       2, 0, 0, 5, 0, 0, 1, 1, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"bt_550",    "m-1",        0,     0,       2, 0, 0, 5, 0, 0, 1, 1, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"Kd_490",    "m-1",        0,     0,       2, 0, 0, 5, 0, 0, 1, 1, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"ap_670",    "m-1",        0,     0,       2, 0, 0, 5, 0, 0, 1, 1, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"Turbidity", "NTU",        0,     0,       2, 0, 0, 5, 0, 0, 1, 1, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"Fluorescence", "mg chla m-3",        0,     0,       2, 0, 0, 5, 0, 0, 1, 1, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"Coral_IN_up","mg N m-2 s-1", 0.0,   0.0, 1, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"Coral_ON_up","mg N m-2 s-1", 0.0,   0.0, 1, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"Gnet",   "mg C m-2 s-1", 0.0,   0.0,      1, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"mucus",   "mg N m-2 s-1", 0.0,   0.0,     1, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"dens",   "kg m-3", 1000.0,  1100.0,       0, 0, 0, 5, 0, 0, 0, 0, 1, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"alk",       "mmol m-3", 2398.2,  2398.2,  0, 1, 1, 5, 1, 0, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"omega_ar",  "nil", 3.0,  3.0,             2, 0, 0, 5, 1, 0, 1, 1, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"omega_ca",  "nil", 3.0,  3.0,             2, 0, 0, 5, 1, 0, 1, 1, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"Age",       "d", 0.0, 0.0,                0, 1, 1, 5, 1, 0, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"passive",   "-", 0.0, 0.0,                0, 1, 1, 5, 1, 0, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"source",  "nil", 0.0, 0.0,                0, 0, 0, 5, 1, 0, 1, 1, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"CS_N_pr",   "mg N m-2 d-1", 0.0,   0.0,   1, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"CH_N_pr",  "g N m-2 d-1", 0.0,   0.0,    1, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"SG_shear_mort",   "g N m-2 d-1", 0.0,   0.0,   1, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"SGH_shear_mort",  "g N m-2 d-1", 0.0,   0.0,    1, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"SGP_shear_mort",  "g N m-2 d-1", 0.0,   0.0,    1, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"SGD_shear_mort",  "g N m-2 d-1", 0.0,   0.0,    1, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"CS_bleach",  "d-1", 0.0,   0.0,    1, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"CS_tempfunc",   "-",0.0,0.0,                2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"temp_clim",   "deg C",0.0 ,26.0,            0, 0, 0, 5, 0, 0, 0, 0, 1, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"OC3M",  "mg Chla m-3", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"OC4Me",  "mg Chla m-3", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"OC3V",  "mg Chla m-3", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"TSSM",  "g TSS m-3", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"KD490M",  "m-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"R_412",  "sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"R_443",  "sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"R_490",  "sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"R_488",  "sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"R_531",  "sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"R_547",  "sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"R_667",  "sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"R_678",  "sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"R_748",  "sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"R_470",  "sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"R_555",  "sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"R_645",  "sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"R_590",  "sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"R_410",  "sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"R_486",  "sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"R_551",  "sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"R_671",  "sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"R_745",  "sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"R_510",  "sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"R_640",  "sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"R_400",  "sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"R_560",  "sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"R_620",  "sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"R_665",  "sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"R_681",  "sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"R_710",  "sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},    
    {"R_709",  "sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"R_753",  "sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"R_754",  "sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"R_482",  "sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"R_655",  "sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"nFLH",  "mW cm-2 um-1 sr-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"Secchi", "m", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"Zenith2D", "rad", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"SWR_bot_abs", "-", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"Oxygen_sedflux", "mg O2 m-2 s-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"DIC_sedflux", "mg C m-2 s-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"NH4_sedflux", "mg N m-2 s-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"NO3_sedflux", "mg N m-2 s-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"DIP_sedflux", "mg P m-2 s-1", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"xco2_in_air", "ppmv",  390.0,    390.0,      0, 0, 0, 5, 0, 0, 1, 1, 1, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"Hue", "degrees", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"Moonlight", "W m-2", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"Lunar_zenith", "radians", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"Lunar_phase", "radians", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"Moon_fulldisk", "W m-2", 0.0,   0.0,    2, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"cdom_pale",       "-", 0.0, 0.0, 0, 1, 1, 5, 1, 0, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"cdom_amber",       "-", 0.0, 0.0, 0, 1, 1, 5, 1, 0, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"cdom_dark",       "-", 0.0, 0.0, 0, 1, 1, 5, 1, 0, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"cdom_gbr",       "-", 0.0, 0.0, 0, 1, 1, 5, 1, 0, 1, 1, 0, FILEIN, ECOLOGY|PROGNOSTIC},
    {"FF_N",      "mg N m-2",  200.0,     0.0,    0, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|PROGNOSTIC},
    {"nFF",   "ind. m-2", 0.0,   0.0,   1, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"FF_N_pr",   "mg C m-2 d-1", 0.0,   0.0,   1, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"FF_N_rm",  "mg C m-2 d-1", 0.0,   0.0,    1, 0, 0, 2, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
    {"NULL",      "NULL",     0.0,   0.0,       0, 0, 0, 0, 0, 0, 0, 0, 0, NOGRAD, ECOLOGY|DIAGNOSTIC},
  };
  
/* diagn, advect, diff, type, diss, part, in water, in sed, flag, obc, type */


  /* Particulate attributes */
  eco_def_partic_t eco_def_partic[] = {
    /* name       psize b_dens i_conc f_conc svel   */
    {"PhyL_N",    1e-5, 1e9,   2e8,   2e8, -3.47e-5 },
    {"PhyL_I",    1e-5, 1e9,   2e6,   2e6, -3.47e-5 },
    {"PhyL_NR",   1e-5, 1e9,   2e6,   2e6, -3.47e-5 },
    {"PhyL_PR",   1e-5, 1e9,   2e6,   2e6, -3.47e-5 },
    {"PhyL_Chl",  1e-5, 1e9,   2e6,   2e6, -3.47e-5 },
    {"PhyS_N",    1e-5, 1e9,   2e6,   2e6,   0.0   },
    {"PhyS_I",    1e-5, 1e9,   2e6,   2e6,   0.0   },
    {"PhyS_NR",   1e-5, 1e9,   2e6,   2e6,   0.0   },
    {"PhyS_PR",   1e-5, 1e9,   2e6,   2e6,   0.0   },
    {"PhyS_Chl",  1e-5, 1e9,   2e6,   2e6,   0.0   },
    {"Tricho_N",  1e-5, 1e9,   2e8,   2e8, -5.6e-6 },
    {"Tricho_NR", 1e-5, 1e9,   2e8,   2e8, -5.6e-6 },
    {"Tricho_PR", 1e-5, 1e9,   2e8,   2e8, -5.6e-6 },
    {"Tricho_I",  1e-5, 1e9,   2e8,   2e8, -5.6e-6 },
    {"Tricho_Chl",1e-5, 1e9,   2e8,   2e8, -5.6e-6 },
    {"MPB_N",     1e-5, 1e9,   2e8,   2e8, -1.16e-4 },
    {"MPB_I",     1e-5, 1e9,   2e6,   2e6, -1.16e-4 },
    {"MPB_NR",    1e-5, 1e9,   2e6,   2e6, -1.16e-4 },
    {"MPB_PR",    1e-5, 1e9,   2e6,   2e6, -1.16e-4 },
    {"MPB_Chl",   1e-5, 1e9,   2e6,   2e6, -1.16e-4 },
    {"ZooL_N",    1e-5, 1e9,   2e8,   2e8,   0.0   },
    {"ZooS_N",    1e-5, 1e9,   2e8,   2e8,   0.0   },
    {"PhyD_N",    1e-5, 1e9,   2e8,   2e8,   0.0   },
    {"PhyD_C",    1e-5, 1e9,   2e8,   2e8,   0.0   },
    {"PhyD_I",    1e-5, 1e9,   2e6,   2e6,   0.0   },
    {"PhyD_NR",   1e-5, 1e9,   2e6,   2e6,   0.0   },
    {"PhyD_PR",   1e-5, 1e9,   2e6,   2e6,   0.0   },
    {"PhyD_Chl",  1e-5, 1e9,   2e6,   2e6,   0.0   },
    {"Phy_L_N2",  1e-5, 1e9,   2e6,   2e6,   0.0   },
    {"PIP",       1.0,  1e9,   2e8,   2e8, -2e-4   },
    {"PIPF",    2.5e-5, 1e9,   2e8,   2e8,   0.0   },
    {"PIPI",    2.5e-5, 1e9,   2e8,   2e8,  -9e-4  },
    {"PIP_Dust",2.5e-5, 1e9,   2e8,   2e8,  -1.1574e-05},
    {"DetR_C",    1e-5, 1e9,   2e6,   2e6, -5.78e-05 },
    {"DetR_N",    1e-5, 1e9,   2e8,   2e8, -5.78e-05 },
    {"DetR_P",    1e-5, 1e9,   2e8,   2e8, -5.78e-05 },
    {"DetPL_N",   1e-5, 1e9,   2e8,   2e8, -5.78e-05 },
    {"DetBL_N",   1e-5, 1e9,   2e8,   2e8, -5.78e-05 },
    {"NULL",      0.0,  0.0,   0.0,   0.0,   0.0   },
  };

  /* Initial the tracer attributes */
  init_tracer_atts_eco(tracer, trname, eco_def, "standard", eco_def_partic);
  
}

/* END eco_defaults_standard()                                       */
/*-------------------------------------------------------------------*/

/*
 * Sets the ecology tracers estuary default attributes
 */
void eco_defaults_est(tracer_info_t *tracer, char *trname, ecology *e)
{
  /*
   * Cut and paste the entire contents of the function
   * eco_defaults_std above and modify the defaults as you see fit.
   *
   * If you think that it would be better to use the above as the
   * baseline and only override a few of the tracers, see Farhan
   */
  quit("'ECO_VARS_ATTS estuary' is not implemented as yet. See the 'standard' counterpart function 'eco_defaults_std' in model/lib/ecology/parameter_defaults.c for an example and populate eco_defaults_est accordingly.\n");  
}

/* END eco_defaults_est()                                            */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Sets the tracer attributes for the given tracer, if found         */
/*-------------------------------------------------------------------*/
void init_tracer_atts_eco(tracer_info_t *tracer, char *trname, 
			  eco_def_t *eco_def, const char *eco_def_name,
			  eco_def_partic_t *eco_def_partic)
{
  int i, n, found;
  
  /* Find the tracer number in the list */
  found = 0;
  for (n = 0; strcmp(eco_def[n].name, "NULL") != 0; ++n) {
    if (strcmp(eco_def[n].name, trname) == 0) {
      found = 1;
      break;
    }
  }

  /* Set the values */
  if (found) {
    /* Allocate ecology private data for this tracer */
    trinfo_priv_eco_t *data = get_private_data_eco(tracer);

    /* Whether this is standard, esturine etc ..*/
    strcpy(data->name, eco_def_name);

    // All of these should really be handled by a generic tracer function
    strcpy(tracer->name, trname);
    strcpy(tracer->long_name, trname); // by default

    // Find the long name and overwrite
    for (i=0; i<NUM_ECO_VARS_3D; i++)
      if (strcmp(ECONAME3D[i][0], trname) == 0) {
	strcpy(tracer->long_name, ECONAME3D[i][1]);
	break;
      }
    if (i == NUM_ECO_VARS_3D) {
      // Try 2D name, if reached the end of the loop above
      for (i=0; i<NUM_ECO_VARS_2D; i++)
	if (strcmp(ECONAME2D[i][0], trname) == 0) {
	  strcpy(tracer->long_name, ECONAME2D[i][1]);
	  break;
	}
    }

    tracer->valid_range_wc[0]  =  0.0;
    tracer->valid_range_sed[0] =  0.0;
    tracer->valid_range_wc[1]  =  1e35;
    tracer->valid_range_sed[1] =  1e35;

    /* Set the standard values */
    strcpy(tracer->units, eco_def[n].units);
    tracer->fill_value_wc  = eco_def[n].fill_wc;
    tracer->fill_value_sed = eco_def[n].fill_sed;
    tracer->diagn   = eco_def[n].diagn;
    tracer->advect  = eco_def[n].advect;
    tracer->diffuse = eco_def[n].diffuse;
    tracer->type    = eco_def[n].type;
    tracer->type    |= eco_def[n].ex_type;
    tracer->inwc    = eco_def[n].inwc;
    tracer->insed   = eco_def[n].insed;
    tracer->dissol  = eco_def[n].dissol;
    tracer->partic  = eco_def[n].partic;

    /* Loop and set all particulates */
    if (tracer->partic) {
      int i;
      found = 0;
      for (i = 0; strcmp(eco_def_partic[i].name, "NULL") != 0; ++i) {
	if (strcmp(eco_def_partic[i].name, trname) == 0) {
	  found = 1;
	  break;
	}
      }
      if (!found)
	hd_quit("Particulate attributes not found for ecology tracer '%s'\n",
		tracer->name);
     
      sed_set_tr_att(tracer, "psize",       &eco_def_partic[i].psize);
      sed_set_tr_att(tracer, "b_dens",      &eco_def_partic[i].b_dens);
      sed_set_tr_att(tracer, "i_conc",      &eco_def_partic[i].i_conc);
      sed_set_tr_att(tracer, "f_conc",      &eco_def_partic[i].f_conc);
      sed_set_tr_att(tracer, "svel",        &eco_def_partic[i].svel);
    }

    /* Ecology specific */
    data->obc  = eco_def[n].obc;
    data->flag = eco_def[n].flag;
  } else
    hd_quit("ecology:init_tracer_atts_eco: Can't find '%s' tracer defaults for %s\n", eco_def_name, trname);
}
#endif

