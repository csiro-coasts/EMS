/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: lib/ecology/parameter_defaults.c
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
 *  $Id: parameter_defaults.c 6353 2019-09-20 14:22:57Z bai155 $
 *
 */

#include <math.h>
#include "utils.h"
#include "parameter_info.h"

/* Local functions */
static void assign_string_values(parameter_info *params, int nprm);
static void init_parameter_values(parameter_info *params);

/** Fills in standard parameter defaults
 * @param nprm Number of parameters
 * @param parameters parameter_info array
 */


void eco_params_bgc2p0(parameter_info **params, int *nprm)
{
  int n;
  parameter_info *parameters = *params;

  init_parameter_values(parameters);

  n = 0;

  /* Temperature / Environment parameters */
  parameters[n].name  = "Tref";
  parameters[n].desc  = "Reference temperature";
  parameters[n].sym   = "T_{ref}";
  parameters[n].units = "Deg C";
  parameters[n].value[0] = 20.0;
  parameters[n].ref = " ";
  parameters[n].index = n++;

  parameters[n].name  = "Q10";
  parameters[n].desc  = "Temperature coefficient for rate parameters";
  parameters[n].sym   = "Q10";
  parameters[n].units = "none";
  parameters[n].value[0] = 2.0;
  parameters[n].ref = " ";
  parameters[n].index = n++;

  parameters[n].name  = "TKEeps";
  parameters[n].desc  = "Nominal rate of TKE dissipation in water column";
  parameters[n].sym   = "\\epsilon";
  parameters[n].units = "m2 s-3";
  parameters[n].value[0] = 1.0e-6;
  parameters[n].index = n++;

  parameters[n].name  = "xco2_in_air";
  parameters[n].desc  = "Atmospheric CO2";
  parameters[n].sym   = "p\\mathrm{CO}_2";
  parameters[n].units = "ppmv";
  parameters[n].value[0] = 396.48;
  parameters[n].ref = "Mean 2013 at Mauna Loa: http://co2now.org/current-co2/co2-now/";
  parameters[n].index = n++;

  parameters[n].name  = "N2";
  parameters[n].desc  = "Concentration of dissolved N2";
  parameters[n].units = "mg N m-3";
  parameters[n].sym   = "[\\mathrm{N}_2]_{gas}";
  parameters[n].value[0] = 2000.0;
  parameters[n].ref = "Robson et al. (2013)";
  parameters[n].index = n++;

  /* Light parameters */

  parameters[n].name = "Light_lambda";
  parameters[n].desc = "Wavelengths of light";
  parameters[n].units= "nm";
  /* Free array allocated above */
  d_free_1d(parameters[n].value);
  parameters[n].num_values = 24;
  parameters[n].value = d_alloc_1d(24);
  parameters[n].value[0] = 290.0;    
  parameters[n].value[1] = 310.0;
  parameters[n].value[2] = 330.0;
  parameters[n].value[3] = 350.0;
  parameters[n].value[4] = 370.0;
  parameters[n].value[5] = 390.0;
  parameters[n].value[6] = 410.0;
  parameters[n].value[7] = 430.0;
  parameters[n].value[8] = 440.0;
  parameters[n].value[9] = 450.0;
  parameters[n].value[10] = 470.0;
  parameters[n].value[11] = 490.0;
  parameters[n].value[12] = 510.0;
  parameters[n].value[13] = 530.0;
  parameters[n].value[14] = 550.0;
  parameters[n].value[15] = 570.0;
  parameters[n].value[16] = 590.0;
  parameters[n].value[17] = 610.0;
  parameters[n].value[18] = 630.0;
  parameters[n].value[19] = 650.0;
  parameters[n].value[20] = 670.0;
  parameters[n].value[21] = 690.0;
  parameters[n].value[22] = 710.0;
  parameters[n].value[23] = 800.0;
  parameters[n].ref = "Approx. 20 nm resolution with 10 nm about 440 nm. PAR (400-700) is integral of bands 6-22.";
  parameters[n].index = n++;

  parameters[n].name  = "acdom443star";
  parameters[n].desc  = "DOC-specific absorption of CDOM 443 nm";
  parameters[n].units = "m2 mg C-1";
  parameters[n].sym   = "k_{CDOM,443}";
  parameters[n].value[0] = 0.00013;
  parameters[n].ref = "Based on Feb 2011 GBR satellite data and modelled DOR_C";
  parameters[n].index = n++;

  parameters[n].name  = "bphy";
  parameters[n].desc  = "Chl-specific scattering coef. for microalgae";
  parameters[n].units = "m-1 (mg Chl a m-3)-1"; 
  parameters[n].sym   = "b_{phy}";
  parameters[n].value[0] = 0.2;
  parameters[n].ref = "Typical microalgae value, Kirk (1994) Light and Photosynthesis in Aquatic ecosystems, Table 4.3";
  parameters[n].index = n++;

  parameters[n].name  = "NtoCHL";
  parameters[n].desc  = "Nominal N:Chl a ratio in phytoplankton by weight";
  parameters[n].units = "g N (g Chl a)-1";
  parameters[n].sym   = "R_{N:Chl}";
  parameters[n].value[0] = 7;
  parameters[n].ref = "Represents a C:Chl ratio of 39.25, Baird et al. (2013) Limnol. Oceanogr. 58: 1215-1226.";
  parameters[n].index = n++;

  /* Phytoplankton */

  parameters[n].name  = "PLumax";
  parameters[n].desc  = "Maximum growth rate of PL at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{PL}^{max}";
  parameters[n].value[0] = 1.4;
  parameters[n].ref = " ";
  parameters[n].index = n++;

  parameters[n].name  = "PLrad";
  parameters[n].desc  = "Radius of the large phytoplankton cells";
  parameters[n].units = "m";
  parameters[n].sym   = "r_{PL}";
  parameters[n].value[0] = 4e-06;
  parameters[n].index = n++;

  parameters[n].name  = "PhyL_mL";
  parameters[n].desc  = "Natural (linear) mortality rate, large phytoplankton";
  parameters[n].units = "d-1";
  parameters[n].sym   = "m_{L,PL}";
  parameters[n].value[0] = 0.1;
  parameters[n].index = n++;

  parameters[n].name  = "PhyL_mL_sed";
  parameters[n].desc  = "Natural (linear) mortality rate in sed., large phyto.";
  parameters[n].units = "d-1";
  parameters[n].sym   = "m_{L,PL,sed}";
  parameters[n].value[0] = 10.0;
  parameters[n].index = n++;

  parameters[n].name  = "PSumax";
  parameters[n].desc  = "Maximum growth rate of PS at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{PL}^{max}";
  parameters[n].value[0] = 1.6;
  parameters[n].index = n++;

  parameters[n].name  = "PSrad";
  parameters[n].desc  = "Radius of the small phytoplankton cells";
  parameters[n].units = "m";
  parameters[n].sym   = "r_{PS}";
  parameters[n].value[0] = 1.e-06;
  parameters[n].index = n++;

  parameters[n].name  = "PhyS_mL";
  parameters[n].desc  = "Natural (linear) mortality rate, small phyto.";
  parameters[n].units = "d-1";
  parameters[n].sym   = "m_{L,PS}";
  parameters[n].value[0] = 0.1;
  parameters[n].index = n++;

  parameters[n].name  = "PhyS_mL_sed";
  parameters[n].desc  = "Natural (linear) mortality rate in sed., small phyto.";
  parameters[n].units = "d-1";
  parameters[n].sym   = "m_{L,PS,sed}";
  parameters[n].value[0] = 1.0;
  parameters[n].index = n++;

  parameters[n].name  = "MBumax";
  parameters[n].desc  = "Maximum growth rate of MB at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{MPB}^{max}";
  parameters[n].value[0] = 0.839;
  parameters[n].index = n++;

  parameters[n].name  = "MBrad";
  parameters[n].desc  = "Radius of the MPB cells";
  parameters[n].units = "m";
  parameters[n].sym   = "r_{MPB}";
  parameters[n].value[0] = 1e-05;
  parameters[n].index = n++;

  parameters[n].name  = "MPB_mQ";
  parameters[n].desc  = "Natural (quadratic) mortality rate, MPB (in sed)";
  parameters[n].units = "d-1 (mg N m-3)-1";
  parameters[n].sym   = "m_{Q,MPB}";
  parameters[n].value[0] = 0.0001;
  parameters[n].ref = "At steady-state, at mu = 0.1 d-1, indep. of temp, MPB_N ~ 0.1 / MPB_mQ = 250 mg N m-3";
  parameters[n].index = n++;

  parameters[n].name  = "PSxan2chl";
  parameters[n].desc  = "Ratio of xanthophyll to chl a of PS";
  parameters[n].units = "mg mg-1";
  parameters[n].sym   = "\\Theta_{xan2chl,PS}";
  parameters[n].value[0] = 0.51;
  parameters[n].ref = "CSIRO parameter library: GBR region";
  parameters[n].index = n++;

  parameters[n].name  = "PLxan2chl";
  parameters[n].desc  = "Ratio of xanthophyll to chl a of PL";
  parameters[n].units = "mg mg-1";
  parameters[n].sym   = "\\Theta_{xan2chl,PL}";
  parameters[n].value[0] = 0.81;
  parameters[n].ref = "CSIRO parameter library: GBR region";
  parameters[n].index = n++;

  parameters[n].name  = "MBxan2chl";
  parameters[n].desc  = "Ratio of xanthophyll to chl a of MPB";
  parameters[n].units = "mg mg-1";
  parameters[n].sym   = "\\Theta_{xan2chl,MPB}";
  parameters[n].value[0] = 0.81;
  parameters[n].ref = "CSIRO parameter library: GBR region WC values";
  parameters[n].index = n++;

  /* Trichodesmium */
  parameters[n].name  = "Tricho_umax";
  parameters[n].desc  = "Maximum growth rate of Trichodesmium at Tref ";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{MPB}^{max}";
  parameters[n].value[0] = 0.24;
  parameters[n].index = n++;

  parameters[n].name  = "Tricho_rad";
  parameters[n].desc  = "Radius of Trichodesmium colonies";
  parameters[n].units = "m";
  parameters[n].sym   = "r_{MPB}";
  parameters[n].value[0] = 0.000005;
  parameters[n].index = n++;

  parameters[n].name  = "Tricho_Sh";
  parameters[n].desc  = "Sherwood number for the Tricho dimensionless";
  parameters[n].units = "none";
  parameters[n].sym   = "Sh_{Tricho}";
  parameters[n].value[0] = 1.0;
  parameters[n].index = n++;

  parameters[n].name  = "Tricho_mL";
  parameters[n].desc  = "Linear mortality for Tricho in sediment";
  parameters[n].units = "d-1";
  parameters[n].sym   = "m_{L,Tricho}";
  parameters[n].value[0] = 0.140000;
  parameters[n].index = n++;

  parameters[n].name  = "Tricho_mQ";
  parameters[n].desc  = "Quadratic mortality for Tricho due to phages in wc";
  parameters[n].units = "d-1 (mg N m-3)-1";
  parameters[n].sym   = "m_{Q,Tricho}";
  parameters[n].value[0] = 0.1000;
  parameters[n].ref = "At steady-state, indep. of temp, Tricho_N ~ Tricho_umax / Tricho_mQ = 0.27 / 0.405 = 0.7 mg N m-3 ~ 0.1 mg Chl m-3";
  parameters[n].index = n++;

  parameters[n].name  = "Tricho_crit";
  parameters[n].desc  = "Critical Tricho above which quadratic mortality applies";
  parameters[n].units = "mg N m-3";
  parameters[n].value[0] = 0.0002000;
  parameters[n].ref = "No source";
  parameters[n].index = n++;

  parameters[n].name  = "p_min";
  parameters[n].desc  = "Minimum density of Trichodesmium";
  parameters[n].units = "kg m-3";
  parameters[n].sym   = "\\rho_{min,Tricho}";
  parameters[n].value[0] = 990.000000;
  parameters[n].ref = "No source";
  parameters[n].index = n++;

  parameters[n].name  = "p_max";
  parameters[n].desc  = "Maximum density of Trichodesmium";
  parameters[n].units = "kg m-3";
  parameters[n].sym   = "\\rho_{max,Tricho}";
  parameters[n].value[0] = 1026.000000;
  parameters[n].ref = "No source";
  parameters[n].index = n++;

  parameters[n].name  = "DINcrit";
  parameters[n].desc  = "DIN conc below which Tricho N fixes ";
  parameters[n].units = "mg N m-3";
  parameters[n].sym   = "DIN_{crit}";
  parameters[n].value[0] = 10.0;
  parameters[n].ref = "Lower end of Robson et al., (2013) 4-20 mg N m-3";
  parameters[n].index = n++;

  parameters[n].name  = "Trichoxan2chl";
  parameters[n].desc  = "Ratio of xanthophyll to chl a of Trichodesmium";
  parameters[n].units = "mg mg-1";
  parameters[n].sym   = "\\Theta_{xan2chl,Tricho}";
  parameters[n].value[0] = 0.5;
  parameters[n].ref = "Subramaniam (1999) LO 44:608-617. Actually redder pigment than xanthophyll.";
  parameters[n].index = n++;

  parameters[n].name  = "C2Chlmin";
  parameters[n].desc  = "Minimum carbon to chlorophyll a ratio";
  parameters[n].units = "wt/wt";
  parameters[n].sym   = "\\theta_{min}";
  parameters[n].value[0] = 20.0;
  parameters[n].ref = "From HPLC in Sathyendranath et al., 2009 MEPS 383,73-84";
  parameters[n].index = n++;


  /* Zooplankton */
  parameters[n].name  = "ZSumax";
  parameters[n].desc  = "Maximum growth rate of ZS at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{max}^{ZS}";
  parameters[n].value[0] = 4.0;
  parameters[n].index = n++;

  parameters[n].name  = "ZSrad";
  parameters[n].desc  = "Radius of the small zooplankton cells";
  parameters[n].units = "m";
  parameters[n].sym   = "r_{ZS}";
  parameters[n].value[0] = 5.0e-06;
  parameters[n].index = n++;

  parameters[n].name  = "ZSswim";
  parameters[n].desc  = "Swimming velocity for small zooplankton";
  parameters[n].units = "m s-1";
  parameters[n].sym   = "U_{ZS}";
  parameters[n].value[0] = 2.0e-4;
  parameters[n].index = n++;

  parameters[n].name  = "ZSmeth";
  parameters[n].desc  = "Grazing technique of small zooplankton";
  parameters[n].units = "none";
  parameters[n].stringvalue = "rect";
  parameters[n].index = n++;

  parameters[n].name  = "ZLumax";
  parameters[n].desc  = "Maximum growth rate of ZL at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{max}^{ZL}";
  parameters[n].value[0] = 1.33;
  parameters[n].index = n++;

  parameters[n].name  = "ZLrad";
  parameters[n].desc  = "Radius of the large zooplankton cells";
  parameters[n].units = "m";
  parameters[n].sym   = "r_{ZL}";
  parameters[n].value[0] = 3.20e-04;
  parameters[n].index = n++;

  parameters[n].name  = "ZLswim";
  parameters[n].desc  = "Swimming velocity for large zooplankton";
  parameters[n].units = "m s-1";
  parameters[n].sym   = "U_{ZL}";
  parameters[n].value[0] = 3.0e-3;
  parameters[n].index = n++;

  parameters[n].name  = "ZLmeth";
  parameters[n].desc  = "Grazing technique of large zooplankton";
  parameters[n].units = "none";
  parameters[n].stringvalue = "rect";
  parameters[n].index = n++;

  parameters[n].name  = "ZL_E";
  parameters[n].desc  = "Growth efficiency, large zooplankton";
  parameters[n].units = "none";
  parameters[n].sym   = "E_{ZL}";
  parameters[n].value[0] = 0.426;
  parameters[n].stderr = 0.0179;
  parameters[n].ref = "Baird and Suthers, 2007 from Hansen et al (1997) LO 42: 687-704";
  parameters[n].index = n++;

  parameters[n].name  = "ZS_E";
  parameters[n].desc  = "Growth efficiency, small zooplankton";
  parameters[n].units = "none";
  parameters[n].sym   = "E_{ZS}";
  parameters[n].value[0] = 0.462;
  parameters[n].stderr = 0.0266;
  parameters[n].ref = "Baird and Suthers, 2007 from Hansen et al (1997) LO 42: 687-704";
  parameters[n].index = n++;

  parameters[n].name  = "ZL_mQ";
  parameters[n].desc  = "Natural (quadratic) mortality rate, large zooplankton";
  parameters[n].sym   = "m_{Q,ZL}";
  parameters[n].units = "d-1 (mg N m-3)-1";
  parameters[n].value[0] = 0.012;
  parameters[n].index = n++;

  parameters[n].name  = "ZS_mQ";
  parameters[n].desc  = "Natural (quadratic) mortality rate, small zooplankton";
  parameters[n].units = "d-1 (mg N m-3)-1";
  parameters[n].sym   = "m_{Q,ZS}";
  parameters[n].value[0] = 0.007;
  parameters[n].index = n++;

  parameters[n].name  = "ZL_FDG";
  parameters[n].desc  = "Fraction of growth inefficiency lost to detritus, large zoo.";
  parameters[n].units = "none";
  parameters[n].sym   = "\\gamma_{ZL}";
  parameters[n].value[0] = 0.5;
  parameters[n].index = n++;

  parameters[n].name  = "ZL_FDM";
  parameters[n].desc  = "Fraction of mortality lost to detritus, large zoo.";
  parameters[n].units = "none";
  parameters[n].sym   = "N/A";
  parameters[n].value[0] = 1.0;
  parameters[n].index = n++;

  parameters[n].name  = "ZS_FDG";
  parameters[n].desc  = "Fraction of growth inefficiency lost to detritus, small zoo.";
  parameters[n].units = "none";
  parameters[n].sym   = "\\gamma_{ZS}";
  parameters[n].value[0] = 0.5;
  parameters[n].index = n++;

  parameters[n].name  = "ZS_FDM";
  parameters[n].desc  = "Fraction of mortality lost to detritus, small zooplankton";
  parameters[n].units = "none";
  parameters[n].sym   = "N/A";
  parameters[n].value[0] = 1.0;
  parameters[n].index = n++;

  /* Remineralisation */
  parameters[n].name  = "F_LD_RD";
  parameters[n].desc  = "Fraction of labile detritus converted to refractory detritus";
  parameters[n].units = "none";
  parameters[n].sym   = "\\zeta_{Red}";
  parameters[n].value[0] = 0.19;
  parameters[n].index = n++;

  parameters[n].name  = "F_LD_DOM";
  parameters[n].desc  = "Fraction of labile detritus converted to DOM";
  parameters[n].units = "none";
  parameters[n].sym   = "\\vartheta_{Red}";
  parameters[n].value[0] = 0.1;
  parameters[n].index = n++;

  parameters[n].name  = "F_RD_DOM";
  parameters[n].desc  = "fraction of refractory detritus that breaks down to DOM";
  parameters[n].units = "none";
  parameters[n].sym   = "\\vartheta_{Ref}";
  parameters[n].value[0] = 0.05;
  parameters[n].index = n++;

  parameters[n].name  = "r_DetPL";
  parameters[n].desc  = "Breakdown rate of labile detritus at 106:16:1";
  parameters[n].units = "d-1";
  parameters[n].sym   = "r_{Red}";
  parameters[n].value[0] = 0.04;
  parameters[n].index = n++;

  parameters[n].name  = "r_DetBL";
  parameters[n].desc  = "Breakdown rate of labile detritus at 550:30:1";
  parameters[n].units = "d-1";
  parameters[n].sym   = "r_{Atk}";
  parameters[n].value[0] = 0.001;
  parameters[n].index = n++;

  parameters[n].name  = "r_RD";
  parameters[n].desc  = "Breakdown rate of refractory detritus";
  parameters[n].units = "d-1";
  parameters[n].sym   = "r_{R}";
  parameters[n].value[0] = 0.001;
  parameters[n].index = n++;

  parameters[n].name  = "r_DOM";
  parameters[n].desc  = "Breakdown rate of dissolved organic matter";
  parameters[n].units = "d-1";
  parameters[n].sym   = "r_{O}";
  parameters[n].value[0] = 0.0001;
  parameters[n].index = n++;

  parameters[n].name  = "Plank_resp";
  parameters[n].desc  = "Respiration as a fraction of umax";
  parameters[n].units = "none";
  parameters[n].sym   = "\\phi";
  parameters[n].value[0] = 0.025;
  parameters[n].index = n++;

  /* Sediment parameters */

  parameters[n].name  = "KO_aer";
  parameters[n].desc  = "Oxygen half-saturation for aerobic respiration";
  parameters[n].units = "mg O m-3";
  parameters[n].sym   = "K_{OA}";
  parameters[n].value[0] = 256.0;
  parameters[n].index = n++;

  parameters[n].name  = "r_nit_wc";
  parameters[n].desc  = "Maximum nitrification rate in water column";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\tau_{nit,wc}";
  parameters[n].value[0] = 0.1;
  parameters[n].index = n++;

  parameters[n].name  = "r_nit_sed";
  parameters[n].desc  = "Maximum nitrification rate in water sediment";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\tau_{nit,sed}";
  parameters[n].value[0] = 20.0;
  parameters[n].index = n++;

  parameters[n].name  = "KO_nit";
  parameters[n].desc  = "Oxygen half-saturation for nitrification";
  parameters[n].units = "mg O m-3";
  parameters[n].sym   = "K_{\\mathrm{O}_2,nit}";
  parameters[n].value[0] = 500.0;
  parameters[n].index = n++;

  parameters[n].name  = "Pads_r";
  parameters[n].desc  = "Rate at which P reaches adsorbed/desorbed equilibrium";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\tau_{Pabs}";
  parameters[n].value[0] = 0.04;
  parameters[n].index = n++;

  parameters[n].name  = "Pads_Kwc";
  parameters[n].desc  = "Freundlich Isothermic Const P adsorption to TSS in wc";
  parameters[n].units = "mg P kg TSS-1";
  parameters[n].sym   = "k_{Pads,wc}";
  parameters[n].value[0] = 30.0;
  parameters[n].index = n++;

  parameters[n].name  = "Pads_Ksed";
  parameters[n].desc  = "Freundlich Isothermic Const P adsorption to TSS in sed";
  parameters[n].units = "mg P kg TSS-1";
  parameters[n].sym   = "k_{Pads,sed}";
  parameters[n].value[0] = 74.0;
  parameters[n].index = n++;

  parameters[n].name  = "Pads_KO";
  parameters[n].desc  = "Oxygen half-saturation for P adsorption";
  parameters[n].units = "mg O m-3";
  parameters[n].sym   = "K_{\\mathrm{O}_2,abs}";
  parameters[n].value[0] = 2000.0;
  parameters[n].index = n++;

  parameters[n].name  = "Pads_exp";
  parameters[n].desc  = "Exponent for Freundlich Isotherm";
  parameters[n].units = "none";
  parameters[n].sym   = "N/A";
  parameters[n].value[0] = 1.0;
  parameters[n].index = n++;

  parameters[n].name  = "r_den";
  parameters[n].desc  = "Maximum denitrification rate";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\tau_{denit}";
  parameters[n].value[0] = 1.0;   // 5.0 in B1p9.
  parameters[n].index = n++;

  parameters[n].name  = "KO_den";
  parameters[n].desc  = "Oxygen half-saturation constant for denitrification";
  parameters[n].units = "mg O m-3";
  parameters[n].sym   = "K_{\\mathrm{O}_2,denit}";
  parameters[n].value[0] = 10000.0;
  parameters[n].index = n++;

  parameters[n].name  = "r_immob_PIP";
  parameters[n].desc  = "Rate of conversion of PIP to immobilised PIP";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\tau_{Pimm}";
  parameters[n].value[0] = 0.0012;
  parameters[n].index = n++;

  parameters[n].name  = "EpiDiffCoeff";
  parameters[n].desc  = "Sediment-water diffusion coefficient";
  parameters[n].units = "m2 s-1";
  parameters[n].sym   = "D";
  parameters[n].value[0] = 3e-7;
  parameters[n].index = n++;

  parameters[n].name  = "EpiDiffDz";
  parameters[n].desc  = "Thickness of diffusive layer";
  parameters[n].units = "m";
  parameters[n].sym   = "h";
  parameters[n].value[0] = 0.0065;
  parameters[n].index = n++;
  
  /* Marcroalgae */
  parameters[n].name  = "MAumax";
  parameters[n].desc  = "Maximum growth rate of MA at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{MA}^{max}";
  parameters[n].value[0] = 1.0;
  parameters[n].index = n++;

  parameters[n].name  = "MA_mL";
  parameters[n].desc  = "Natural (linear) mortality rate, macroalgae";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\zeta_{MA}";
  parameters[n].value[0] = 0.01;
  parameters[n].index = n++;

  parameters[n].name  = "MAleafden";
  parameters[n].desc  = "Nitrogen-specific leaf area of macroalgae";
  parameters[n].units = "m2 g N-1";
  parameters[n].sym   = "\\Omega_{MA}";
  parameters[n].value[0] = 1.0;
  parameters[n].index = n++;

  parameters[n].name  = "Benth_resp";
  parameters[n].desc  = "Respiration as a fraction of umax";
  parameters[n].units = "none";
  parameters[n].sym   = "\\phi";
  parameters[n].value[0] = 0.025;
  parameters[n].index = n++;

  /* Seagrass parameters - Zostera */
  parameters[n].name  = "SGumax";
  parameters[n].desc  = "Maximum growth rate of SG at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{SG}^{max}";
  parameters[n].value[0] = 0.4;
  parameters[n].ref = "x2 nighttime, x2 for roots.";
  parameters[n].index = n++;
  
  parameters[n].name  = "SG_KN";
  parameters[n].desc  = "Half-saturation of SG N uptake in SED";
  parameters[n].units = "mg N m-3";
  parameters[n].sym   = "K_{SG,N}";
  parameters[n].value[0] = 420.0;
  parameters[n].ref = "Lee and Dunton (1999) 1204-1215. Table 3 Zostera";
  parameters[n].index = n++;

  parameters[n].name  = "SG_KP";
  parameters[n].desc  = "Half-saturation of SG P uptake in SED";
  parameters[n].units = "mg P m-3";
  parameters[n].sym   = "K_{SG,P}";
  parameters[n].value[0] = 96.0;
  parameters[n].ref = "Gras et al. (2003) Aquatic Botany 76:299-315. Thalassia testudinum.";
  parameters[n].index = n++;

  parameters[n].name  = "SG_mL";
  parameters[n].desc  = "Natural (linear) mortality rate aboveground seagrass";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\zeta_{SG_A}";
  parameters[n].value[0] = 0.03;
  parameters[n].stderr = 0.001;
  parameters[n].ref = "Fourquean et al.( 2003) Chem. Ecol. 19: 373-390.Thalassia leaves with one component decay";
  parameters[n].index = n++;

  parameters[n].name  = "SGROOT_mL";
  parameters[n].desc  = "Natural (linear) mortality rate belowground seagrass";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\zeta_{SG_B}";
  parameters[n].value[0] = 0.004;
  parameters[n].stderr = 0.0002;
  parameters[n].ref = "Fourquean et al. (2003) Chem. Ecol. 19: 373-390. Thalassia roots with one component decay";
  parameters[n].index = n++;

  parameters[n].name  = "SGfrac";
  parameters[n].desc  = "Fraction (target) of SG biomass below-ground";
  parameters[n].units = "-";
  parameters[n].sym   = "f_{below,SG}";
  parameters[n].value[0] = 0.75;
  parameters[n].ref = "Babcock (2015) Zostera capricornii.";
  parameters[n].index = n++;

  parameters[n].name  = "SGtransrate";
  parameters[n].desc  = "Time scale for seagrass translocation";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\tau_{tran,SG}";
  parameters[n].value[0] = 0.0333;
  parameters[n].ref = "Loosely based on Zostera marine Kaldy et al., 2013 MEPS 487:27-39";
  parameters[n].index = n++;

  parameters[n].name  = "SGleafden";
  parameters[n].desc  = "Nitrogen-specific leaf area of seagrass";
  parameters[n].units = "m2 g N-1";
  parameters[n].sym   = "\\Omega_{SG}";
  parameters[n].value[0] = 1.5;
  parameters[n].ref = "Zostera capricornia: leaf dimensions Kemp et al (1987) Mar Ecol. Prog. Ser. 41:79-86.";
  parameters[n].index = n++;

  parameters[n].name  = "SGseedfrac";
  parameters[n].desc  = "Seagrass seed biomass as fraction of 63 % cover";
  parameters[n].units = "-";
  parameters[n].sym   = "f_{seed,SG}";
  parameters[n].value[0] = 0.01;
  parameters[n].ref = "No source";
  parameters[n].index = n++;

  parameters[n].name  = "SGorient";
  parameters[n].desc  = "Sine of nadir Zostera canopy bending angle";
  parameters[n].units = "-";
  parameters[n].sym   = "\\sin \\beta_{blade,SG}";
  parameters[n].value[0] = 0.5;
  parameters[n].ref = "No source";
  parameters[n].index = n++;

  parameters[n].name  = "SGmlr";
  parameters[n].desc  = "Compensation irradiance for Zostera";
  parameters[n].units = "mol m-2";
  parameters[n].sym   = "E_{comp,SG}";
  parameters[n].value[0] = 4.5;
  parameters[n].ref = "Chartrand (2012) Tech report.";
  parameters[n].index = n++;

  parameters[n].name  = "SGrootdepth";
  parameters[n].desc  = "Maximum depth for Zostera roots";
  parameters[n].units = "m";
  parameters[n].sym   = "z_{root,SG}";
  parameters[n].value[0] = -0.15;
  parameters[n].ref = "Roberts (1993) Aust. J. Mar. Fresh. Res. 44:85-100.";
  parameters[n].index = n++;

  parameters[n].name  = "SG_tau_critical";
  parameters[n].desc  = "Critical shear stress for SG loss";
  parameters[n].units = "N m-2";
  parameters[n].sym   = "\\tau_{SG,shear}";
  parameters[n].value[0] = 1.0;
  parameters[n].ref = "NESP project";
  parameters[n].index = n++;

  parameters[n].name  = "SG_tau_time";
  parameters[n].desc  = "Time-scale for critical shear stress for SG loss";
  parameters[n].units = "s";
  parameters[n].sym   = "\\tau_{SG,time}";
  parameters[n].value[0] = 43200.0;
  parameters[n].ref = "NESP project";
  parameters[n].index = n++;

  /* Seagrass parameters - Halophila */
  parameters[n].name  = "SGHumax";
  parameters[n].desc  = "Maximum growth rate of SGH at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{SGH}^{max}";
  parameters[n].value[0] = 0.4;
  parameters[n].ref = "x2 nighttime, x2 for roots.";
  parameters[n].index = n++;

  parameters[n].name  = "SGH_KN";
  parameters[n].desc  = "Half-saturation of SGH N uptake in SED";
  parameters[n].units = "mg N m-3";
  parameters[n].sym   = "K_{SGH,N}";
  parameters[n].value[0] = 420.0;
  parameters[n].index = n++;

  parameters[n].name  = "SGH_KP";
  parameters[n].desc  = "Half-saturation of SGH P uptake in SED";
  parameters[n].units = "mg P m-3";
  parameters[n].sym   = "K_{SGH,P}";
  parameters[n].value[0] = 96.0;
  parameters[n].index = n++;

  parameters[n].name  = "SGHleafden";
  parameters[n].desc  = "Nitrogen-specific leaf area of SGH"; 
  parameters[n].units = "m2 g N-1";
  parameters[n].sym   = "\\Omega_{SGH}";
  parameters[n].value[0] = 1.9;
  parameters[n].ref = "Halophila ovalis: leaf dimensions from Vermaat et al. (1995)";
  parameters[n].index = n++;

  parameters[n].name  = "SGH_mL";
  parameters[n].desc  = "Natural (linear) mortality rate, aboveground SGH";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\zeta_{SGH_A}";
  parameters[n].value[0] = 0.06;
  parameters[n].stderr = 0.001;
  parameters[n].ref = "Fourquean et al.(2003) Chem. Ecol. 19: 373-390.Thalassia leaves with one component decay";
  parameters[n].index = n++;

  parameters[n].name  = "SGHROOT_mL";
  parameters[n].desc  = "Natural (linear) mortality rate, belowground SGH";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\zeta_{SGH_B}";
  parameters[n].value[0] = 0.004;
  parameters[n].stderr = 0.0002;
  parameters[n].ref = "Fourquean et al. (2003) Chem. Ecol. 19: 373-390. Thalassia roots with one component decay";
  parameters[n].index = n++;

  parameters[n].name  = "SGHfrac";
  parameters[n].desc  = "Fraction (target) of SGH biomass below-ground";
  parameters[n].units = "-";
  parameters[n].sym   = "f_{below,SGH}";
  parameters[n].value[0] = 0.5;
  parameters[n].ref = "Babcock 2015, Halophila ovalis";
  parameters[n].index = n++;

  parameters[n].name  = "SGHtransrate";
  parameters[n].desc  = "Time scale for Halophila translocation";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\tau_{tran,SGH}";
  parameters[n].value[0] = 0.0333;
  parameters[n].ref = "Loosely based on Zostera marine Kaldy et al., 2013 MEPS 487:27-39";
  parameters[n].index = n++;

  parameters[n].name  = "SGHseedfrac";
  parameters[n].desc  = "Halophila seed biomass as fraction of 63 % cover";
  parameters[n].units = "-";
  parameters[n].sym   = "f_{seed,SGH}";
  parameters[n].value[0] = 0.01;
  parameters[n].ref = "No source";
  parameters[n].index = n++;

  parameters[n].name  = "SGHorient";
  parameters[n].desc  = "Sine of nadir Halophila canopy bending angle";
  parameters[n].units = "-";
  parameters[n].sym   = "\\sin \\beta_{blade,SGH}";
  parameters[n].value[0] = 1.0;
  parameters[n].ref = "No source";
  parameters[n].index = n++;

  parameters[n].name  = "SGHmlr";
  parameters[n].desc  = "Compensation irradiance for Halophila";
  parameters[n].units = "mol m-2";
  parameters[n].sym   = "E_{comp,SGH}";
  parameters[n].value[0] = 2.0;
  parameters[n].ref = "Longstaff 2003 UQ PhD thesis";
  parameters[n].index = n++;

  parameters[n].name  = "SGHrootdepth";
  parameters[n].desc  = "Maximum depth for Halophila roots";
  parameters[n].units = "m";
  parameters[n].sym   = "z_{root,SGH}";
  parameters[n].value[0] = -0.08;
  parameters[n].ref = "Roberts (1993) Aust. J. Mar. Fresh. Res. 44:85-100.";
  parameters[n].index = n++;

  parameters[n].name  = "SGH_tau_critical";
  parameters[n].desc  = "Critical shear stress for SGH loss";
  parameters[n].units = "N m-2";
  parameters[n].sym   = "\\tau_{SGH,shear}";
  parameters[n].value[0] = 1.0;
  parameters[n].ref = "NESP project";
  parameters[n].index = n++;

  parameters[n].name  = "SGH_tau_time";
  parameters[n].desc  = "Time-scale for critical shear stress for SGH loss";
  parameters[n].units = "s";
  parameters[n].sym   = "\\tau_{SGH,time}";
  parameters[n].value[0] = 43200.0;
  parameters[n].ref = "NESP project";
  parameters[n].index = n++;


  /* Seagrass parameters - Deep */
  parameters[n].name  = "SGDumax";
  parameters[n].desc  = "Maximum growth rate of SGD at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{SGD}^{max}";
  parameters[n].value[0] = 0.4;
  parameters[n].ref = "x2 nighttime, x2 for roots.";
  parameters[n].index = n++;

  parameters[n].name  = "SGD_KN";
  parameters[n].desc  = "Half-saturation of SGD N uptake in SED";
  parameters[n].units = "mg N m-3";
  parameters[n].sym   = "K_{SGD,N}";
  parameters[n].value[0] = 420.0;
  parameters[n].index = n++;

  parameters[n].name  = "SGD_KP";
  parameters[n].desc  = "Half-saturation of SGD P uptake in SED";
  parameters[n].units = "mg P m-3";
  parameters[n].sym   = "K_{SGD,P}";
  parameters[n].value[0] = 96.0;
  parameters[n].index = n++;

  parameters[n].name  = "SGDleafden";
  parameters[n].desc  = "Nitrogen-specific leaf area of SGD"; 
  parameters[n].units = "m2 g N-1";
  parameters[n].sym   = "\\Omega_{SGD}";
  parameters[n].value[0] = 1.9;
  parameters[n].ref = "Halophila ovalis: leaf dimensions from Vermaat et al. (1995)";
  parameters[n].index = n++;

  parameters[n].name  = "SGD_mL";
  parameters[n].desc  = "Natural (linear) mortality rate, aboveground SGD";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\zeta_{SGD_A}";
  parameters[n].value[0] = 0.06;
  parameters[n].stderr = 0.001;
  parameters[n].ref = "NESP project";
  parameters[n].index = n++;

  parameters[n].name  = "SGDROOT_mL";
  parameters[n].desc  = "Natural (linear) mortality rate, belowground SGD";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\zeta_{SGD_B}";
  parameters[n].value[0] = 0.004;
  parameters[n].stderr = 0.00002;
  parameters[n].ref = "NESP project";
  parameters[n].index = n++;

  parameters[n].name  = "SGDfrac";
  parameters[n].desc  = "Fraction (target) of SGD biomass below-ground";
  parameters[n].units = "-";
  parameters[n].sym   = "f_{below,SGD}";
  parameters[n].value[0] = 0.25;
  parameters[n].ref = "Duarte (1999) Aquatic Biol. 65: 159-174, Halophila ovalis.";
  parameters[n].index = n++;

  parameters[n].name  = "SGDtransrate";
  parameters[n].desc  = "Time scale for deep SG translocation";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\tau_{tran,SGD}";
  parameters[n].value[0] = 0.0333;
  parameters[n].ref = "Loosely based on Zostera marine Kaldy et al., 2013 MEPS 487:27-39";
  parameters[n].index = n++;

  parameters[n].name  = "SGDseedfrac";
  parameters[n].desc  = "Deep SG seed biomass as fraction of 63 % cover";
  parameters[n].units = "-";
  parameters[n].sym   = "f_{seed,SGD}";
  parameters[n].value[0] = 0.01;
  parameters[n].ref = "No source";
  parameters[n].index = n++;

  parameters[n].name  = "SGDorient";
  parameters[n].desc  = "Sine of nadir deep SG canopy bending angle";
  parameters[n].units = "-";
  parameters[n].sym   = "\\sin \\beta_{blade,SGD}";
  parameters[n].value[0] = 1.0;
  parameters[n].ref = "No source";
  parameters[n].index = n++;

  parameters[n].name  = "SGDmlr";
  parameters[n].desc  = "Compensation irradiance for deep SG";
  parameters[n].units = "mol m-2";
  parameters[n].sym   = "E_{comp,SGD}";
  parameters[n].value[0] = 1.5;
  parameters[n].ref = "Chartrand (2017) Tech report.";
  parameters[n].index = n++;

  parameters[n].name  = "SGDrootdepth";
  parameters[n].desc  = "Maximum depth for deep SG roots";
  parameters[n].units = "m";
  parameters[n].sym   = "z_{root,SGD}";
  parameters[n].value[0] = -0.05;
  parameters[n].ref = "NESP project";
  parameters[n].index = n++;

  parameters[n].name  = "SGD_tau_critical";
  parameters[n].desc  = "Critical shear stress for deep SG loss";
  parameters[n].units = "N m-2";
  parameters[n].sym   = "\\tau_{SGD,shear}";
  parameters[n].value[0] = 1.0;
  parameters[n].ref = "NESP project";
  parameters[n].index = n++;

  parameters[n].name  = "SGD_tau_time";
  parameters[n].desc  = "Time-scale for shear stress for deep SG loss";
  parameters[n].units = "s";
  parameters[n].sym   = "\\tau_{SGD,time}";
  parameters[n].value[0] = 43200.0;
  parameters[n].ref = "NESP project";
  parameters[n].index = n++;

  /* Corals */
  parameters[n].name  = "dissCaCO3_sed";
  parameters[n].desc  = "net dissolution rate of sediment without coral";
  parameters[n].units = "mmol C m-2 s-1";
  parameters[n].sym   = "d_{sand}";
  parameters[n].value[0] = 0.001;   // 0.007 in B1p9
  parameters[n].index = n++;
  parameters[n].ref = "Anthony et al. (2013), Biogeosciences 10:4897-4909, Fig 5E: -1 2 3 6  mmol m-2 h-1";

  parameters[n].name  = "CHarea";
  parameters[n].desc  = "Grid scale to reef scale ratio";
  parameters[n].units = "m2 m-2";
  parameters[n].sym   = "A_{CH}";
  parameters[n].value[0] = 0.1;
  parameters[n].index = n++;
  parameters[n].ref = "Heron Island for on 4 km model.";

  parameters[n].name  = "CHpolypden";
  parameters[n].desc  = "Nitrogen-specific host area of coral polyp";
  parameters[n].units = "m2 g N-1";
  parameters[n].sym   = "\\Omega_{CH}";
  parameters[n].value[0] = 2.0;
  parameters[n].index = n++;

  parameters[n].name  = "CHumax";
  parameters[n].desc  = "Max. growth rate of Coral at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{CH}^{max}";
  parameters[n].value[0] = 0.05;
  parameters[n].index = n++;

  parameters[n].name  = "CSumax";
  parameters[n].desc  = "Max. growth rate of zooxanthellae at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{CS}^{max}";
  parameters[n].value[0] = 0.4;
  parameters[n].index = n++;

  parameters[n].name  = "CSrad";
  parameters[n].desc  = "Radius of the zooxanthellae ";
  parameters[n].units = "m";
  parameters[n].sym   = "r_{CS}";
  parameters[n].value[0] = 5e-06;
  parameters[n].index = n++;

  parameters[n].name  = "CHmort";
  parameters[n].desc  = "Quadratic mortality rate of coral polyp ";
  parameters[n].units = "(g N m-2)-1 d-1";
  parameters[n].sym   = "\\zeta_{CH}";
  parameters[n].value[0] = 0.01;
  parameters[n].index = n++;

  parameters[n].name  = "CSmort";
  parameters[n].desc  = "Linear mortality rate of zooxanthellae ";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\zeta_{CS}";
  parameters[n].value[0] = 0.04;
  parameters[n].index = n++;

  parameters[n].name  = "CHremin";
  parameters[n].desc  = "Fraction of coral host death translocated.";
  parameters[n].units = "-";
  parameters[n].sym   = "f_{remin}";
  parameters[n].value[0] = 0.5;
  parameters[n].index = n++;

  parameters[n].name  = "Splank";
  parameters[n].desc  = "Rate coefficent for particle uptake by corals";
  parameters[n].units = "m d-1";
  parameters[n].sym   = "S_{part}";
  parameters[n].value[0] = 3.0;
  parameters[n].ref = "Ribes and Atkinson (2007) Coral Reefs 26: 413-421";
  parameters[n].index = n++;

  parameters[n].name  = "k_day_coral";
  parameters[n].desc  = "Maximum daytime coral calcification";
  parameters[n].units = "mmol C m-2 s-1";
  parameters[n].sym   = "k_{day}";
  parameters[n].value[0] = 0.0132;
  parameters[n].ref = "Anthony et al. (2013), Biogeosciences 10:4897-4909, Fig 5A: 50, 50, 35 55 mmol m-2 h-1 for Acropora aspera n=4";
  parameters[n].index = n++;

  parameters[n].name  = "k_night_coral";
  parameters[n].desc  = "Maximum nightime coral calcification";
  parameters[n].units = "mmol C m-2 s-1";
  parameters[n].sym   = "k_{night}";
  parameters[n].value[0] = 0.0069;
  parameters[n].ref = "Anthony et al. (2013), Biogeosciences 10:4897-4909, Fig 5A: 20, 30, 20, 30  mmol m-2 h-1 for Acropora aspera n=4";
  parameters[n].index = n++;

  parameters[n].name  = "dissCaCO3_shelf";
  parameters[n].desc  = "Carbonate sediment dissolution rate on shelf";
  parameters[n].units = "mmol C m-2 s-1";
  parameters[n].sym   = "d_{shelf}";
  parameters[n].value[0] = 0.0001;
  parameters[n].ref = "Cyronak, T. et al., LO 58:131-143. Heron Island study.";
  parameters[n].index = n++;

  parameters[n].name  = "ageing_decay";
  parameters[n].desc  = "Age tracer growth rate per day";
  parameters[n].units = "d d-1";
  parameters[n].sym   = "n/a";
  parameters[n].value[0] = 1.0;
  parameters[n].ref = "EMS manual";
  parameters[n].index = n++;

  parameters[n].name  = "anti_ageing_decay";
  parameters[n].desc  = "Age tracer decay rate per day outside source";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\Phi";
  parameters[n].value[0] = 0.1;
  parameters[n].ref = "EMS manual";
  parameters[n].index = n++;


  /* Assign acutal number of parameters */
  *nprm = n;

  /* Assign string values */
  assign_string_values(parameters, *nprm);
}

void eco_params_bgc3p1(parameter_info **params, int *nprm)
{
  int n;
  parameter_info *parameters = *params;

  init_parameter_values(parameters);

  n = 0;

  /* Temperature / Environment parameters */
  parameters[n].name  = "Tref";
  parameters[n].desc  = "Reference temperature";
  parameters[n].sym   = "T_{ref}";
  parameters[n].units = "Deg C";
  parameters[n].value[0] = 20.0;
  parameters[n].ref = " ";
  parameters[n].index = n++;

  parameters[n].name  = "Q10";
  parameters[n].desc  = "Temperature coefficient for rate parameters";
  parameters[n].sym   = "Q10";
  parameters[n].units = "none";
  parameters[n].value[0] = 2.0;
  parameters[n].ref = " ";
  parameters[n].index = n++;

  parameters[n].name  = "TKEeps";
  parameters[n].desc  = "Nominal rate of TKE dissipation in water column";
  parameters[n].sym   = "\\epsilon";
  parameters[n].units = "m2 s-3";
  parameters[n].value[0] = 1.0e-6;
  parameters[n].index = n++;

  // parameters[n].name  = "xco2_in_air";
  // parameters[n].desc  = "Atmospheric CO2";
  // parameters[n].sym   = "p\\mathrm{CO}_2";
  // parameters[n].units = "ppmv";
  // parameters[n].value[0] = 396.48;
  // parameters[n].ref = "Mean 2013 at Mauna Loa: http://co2now.org/current-co2/co2-now/";
  // parameters[n].index = n++;

  parameters[n].name  = "N2";
  parameters[n].desc  = "Concentration of dissolved N2";
  parameters[n].units = "mg N m-3";
  parameters[n].sym   = "[\\mathrm{N}_2]_{gas}";
  parameters[n].value[0] = 2000.0;
  parameters[n].ref = "Robson et al. (2013)";
  parameters[n].index = n++;

  /* Light parameters */

  parameters[n].name = "Light_lambda";
  parameters[n].desc = "Wavelengths of light";
  parameters[n].units= "nm";
  /* Free array allocated above */
  d_free_1d(parameters[n].value);
  parameters[n].num_values = 24;
  parameters[n].value = d_alloc_1d(24);
  parameters[n].value[0] = 290.0;    
  parameters[n].value[1] = 310.0;
  parameters[n].value[2] = 330.0;
  parameters[n].value[3] = 350.0;
  parameters[n].value[4] = 370.0;
  parameters[n].value[5] = 390.0;
  parameters[n].value[6] = 410.0;
  parameters[n].value[7] = 430.0;
  parameters[n].value[8] = 440.0;
  parameters[n].value[9] = 450.0;
  parameters[n].value[10] = 470.0;
  parameters[n].value[11] = 490.0;
  parameters[n].value[12] = 510.0;
  parameters[n].value[13] = 530.0;
  parameters[n].value[14] = 550.0;
  parameters[n].value[15] = 570.0;
  parameters[n].value[16] = 590.0;
  parameters[n].value[17] = 610.0;
  parameters[n].value[18] = 630.0;
  parameters[n].value[19] = 650.0;
  parameters[n].value[20] = 670.0;
  parameters[n].value[21] = 690.0;
  parameters[n].value[22] = 710.0;
  parameters[n].value[23] = 800.0;
  parameters[n].ref = "Approx. 20 nm resolution with 10 nm about 440 nm. PAR (400-700) is integral of bands 6-22.";
  parameters[n].index = n++;

  parameters[n].name  = "acdom443star";
  parameters[n].desc  = "DOC-specific absorption of CDOM 443 nm";
  parameters[n].units = "m2 mg C-1";
  parameters[n].sym   = "k_{CDOM,443}";
  parameters[n].value[0] = 0.00013;
  parameters[n].ref = "Based on Feb 2011 GBR satellite data and modelled DOR_C";
  parameters[n].index = n++;

  parameters[n].name  = "bphy";
  parameters[n].desc  = "Chl-specific scattering coef. for microalgae";
  parameters[n].units = "m-1 (mg Chl a m-3)-1"; 
  parameters[n].sym   = "b_{phy}";
  parameters[n].value[0] = 0.2;
  parameters[n].ref = "Typical microalgae value, Kirk (1994) Light and Photosynthesis in Aquatic ecosystems, Table 4.3";
  parameters[n].index = n++;

  parameters[n].name  = "NtoCHL";
  parameters[n].desc  = "Nominal N:Chl a ratio in phytoplankton by weight";
  parameters[n].units = "g N (g Chl a)-1";
  parameters[n].sym   = "R_{N:Chl}";
  parameters[n].value[0] = 7;
  parameters[n].ref = "Represents a C:Chl ratio of 39.25, Baird et al. (2013) Limnol. Oceanogr. 58: 1215-1226.";
  parameters[n].index = n++;

  /* Phytoplankton */

  parameters[n].name  = "PLumax";
  parameters[n].desc  = "Maximum growth rate of PL at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{PL}^{max}";
  parameters[n].value[0] = 1.4;
  parameters[n].ref = " ";
  parameters[n].index = n++;

  parameters[n].name  = "PLrad";
  parameters[n].desc  = "Radius of the large phytoplankton cells";
  parameters[n].units = "m";
  parameters[n].sym   = "r_{PL}";
  parameters[n].value[0] = 4e-06;
  parameters[n].index = n++;

  parameters[n].name  = "PhyL_mL";
  parameters[n].desc  = "Natural (linear) mortality rate, large phytoplankton";
  parameters[n].units = "d-1";
  parameters[n].sym   = "m_{L,PL}";
  parameters[n].value[0] = 0.1;
  parameters[n].index = n++;

  parameters[n].name  = "PhyL_mL_sed";
  parameters[n].desc  = "Natural (linear) mortality rate in sed., large phyto.";
  parameters[n].units = "d-1";
  parameters[n].sym   = "m_{L,PL,sed}";
  parameters[n].value[0] = 10.0;
  parameters[n].index = n++;

  parameters[n].name  = "PSumax";
  parameters[n].desc  = "Maximum growth rate of PS at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{PL}^{max}";
  parameters[n].value[0] = 1.6;
  parameters[n].index = n++;

  parameters[n].name  = "PSrad";
  parameters[n].desc  = "Radius of the small phytoplankton cells";
  parameters[n].units = "m";
  parameters[n].sym   = "r_{PS}";
  parameters[n].value[0] = 1.e-06;
  parameters[n].index = n++;

  parameters[n].name  = "PhyS_mL";
  parameters[n].desc  = "Natural (linear) mortality rate, small phyto.";
  parameters[n].units = "d-1";
  parameters[n].sym   = "m_{L,PS}";
  parameters[n].value[0] = 0.1;
  parameters[n].index = n++;

  parameters[n].name  = "PhyS_mL_sed";
  parameters[n].desc  = "Natural (linear) mortality rate in sed., small phyto.";
  parameters[n].units = "d-1";
  parameters[n].sym   = "m_{L,PS,sed}";
  parameters[n].value[0] = 1.0;
  parameters[n].index = n++;

  parameters[n].name  = "MBumax";
  parameters[n].desc  = "Maximum growth rate of MB at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{MPB}^{max}";
  parameters[n].value[0] = 0.839;
  parameters[n].index = n++;

  parameters[n].name  = "MBrad";
  parameters[n].desc  = "Radius of the MPB cells";
  parameters[n].units = "m";
  parameters[n].sym   = "r_{MPB}";
  parameters[n].value[0] = 1e-05;
  parameters[n].index = n++;

  parameters[n].name  = "MPB_mQ";
  parameters[n].desc  = "Natural (quadratic) mortality rate, MPB (in sed)";
  parameters[n].units = "d-1 (mg N m-3)-1";
  parameters[n].sym   = "m_{Q,MPB}";
  parameters[n].value[0] = 0.0001;
  parameters[n].ref = "At steady-state, at mu = 0.1 d-1, indep. of temp, MPB_N ~ 0.1 / MPB_mQ = 250 mg N m-3";
  parameters[n].index = n++;

  parameters[n].name  = "PSxan2chl";
  parameters[n].desc  = "Ratio of xanthophyll to chl a of PS";
  parameters[n].units = "mg mg-1";
  parameters[n].sym   = "\\Theta_{xan2chl,PS}";
  parameters[n].value[0] = 0.51;
  parameters[n].ref = "CSIRO parameter library: GBR region";
  parameters[n].index = n++;

  parameters[n].name  = "PLxan2chl";
  parameters[n].desc  = "Ratio of xanthophyll to chl a of PL";
  parameters[n].units = "mg mg-1";
  parameters[n].sym   = "\\Theta_{xan2chl,PL}";
  parameters[n].value[0] = 0.81;
  parameters[n].ref = "CSIRO parameter library: GBR region";
  parameters[n].index = n++;

  parameters[n].name  = "MBxan2chl";
  parameters[n].desc  = "Ratio of xanthophyll to chl a of MPB";
  parameters[n].units = "mg mg-1";
  parameters[n].sym   = "\\Theta_{xan2chl,MPB}";
  parameters[n].value[0] = 0.81;
  parameters[n].ref = "CSIRO parameter library: GBR region WC values";
  parameters[n].index = n++;

  /* Trichodesmium */
  parameters[n].name  = "Tricho_umax";
  parameters[n].desc  = "Maximum growth rate of Trichodesmium at Tref ";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{MPB}^{max}";
  parameters[n].value[0] = 0.20;
  parameters[n].index = n++;

  parameters[n].name  = "Tricho_rad";
  parameters[n].desc  = "Radius of Trichodesmium colonies";
  parameters[n].units = "m";
  parameters[n].sym   = "r_{MPB}";
  parameters[n].value[0] = 0.000005;
  parameters[n].index = n++;

  parameters[n].name  = "Tricho_Sh";
  parameters[n].desc  = "Sherwood number for the Tricho dimensionless";
  parameters[n].units = "none";
  parameters[n].sym   = "Sh_{Tricho}";
  parameters[n].value[0] = 1.0;
  parameters[n].index = n++;

  parameters[n].name  = "Tricho_mL";
  parameters[n].desc  = "Linear mortality for Tricho in sediment";
  parameters[n].units = "d-1";
  parameters[n].sym   = "m_{L,Tricho}";
  parameters[n].value[0] = 0.10000;
  parameters[n].index = n++;

  parameters[n].name  = "Tricho_mQ";
  parameters[n].desc  = "Quadratic mortality for Tricho due to phages in wc";
  parameters[n].units = "d-1 (mg N m-3)-1";
  parameters[n].sym   = "m_{Q,Tricho}";
  parameters[n].value[0] = 0.1000;
  parameters[n].ref = "At steady-state, indep. of temp, Tricho_N ~ Tricho_umax / Tricho_mQ = 0.27 / 0.405 = 0.7 mg N m-3 ~ 0.1 mg Chl m-3";
  parameters[n].index = n++;

  parameters[n].name  = "Tricho_crit";
  parameters[n].desc  = "Critical Tricho above which quadratic mortality applies";
  parameters[n].units = "mg N m-3";
  parameters[n].value[0] = 0.0002000;
  parameters[n].ref = "No source";
  parameters[n].index = n++;

  parameters[n].name  = "p_min";
  parameters[n].desc  = "Minimum density of Trichodesmium";
  parameters[n].units = "kg m-3";
  parameters[n].sym   = "\\rho_{min,Tricho}";
  parameters[n].value[0] = 900.000000;
  parameters[n].ref = "No source";
  parameters[n].index = n++;

  parameters[n].name  = "p_max";
  parameters[n].desc  = "Maximum density of Trichodesmium";
  parameters[n].units = "kg m-3";
  parameters[n].sym   = "\\rho_{max,Tricho}";
  parameters[n].value[0] = 1050.000000;
  parameters[n].ref = "No source";
  parameters[n].index = n++;

  parameters[n].name  = "DINcrit";
  parameters[n].desc  = "DIN conc below which Tricho N fixes ";
  parameters[n].units = "mg N m-3";
  parameters[n].sym   = "DIN_{crit}";
  parameters[n].value[0] = 10.0;
  parameters[n].ref = "Lower end of Robson et al., (2013) 4-20 mg N m-3";
  parameters[n].index = n++;

  parameters[n].name  = "Trichoxan2chl";
  parameters[n].desc  = "Ratio of xanthophyll to chl a of Trichodesmium";
  parameters[n].units = "mg mg-1";
  parameters[n].sym   = "\\Theta_{xan2chl,Tricho}";
  parameters[n].value[0] = 0.5;
  parameters[n].ref = "Subramaniam (1999) LO 44:608-617. Actually redder pigment than xanthophyll.";
  parameters[n].index = n++;

  parameters[n].name  = "C2Chlmin";
  parameters[n].desc  = "Minimum carbon to chlorophyll a ratio";
  parameters[n].units = "wt/wt";
  parameters[n].sym   = "\\theta_{min}";
  parameters[n].value[0] = 20.0;
  parameters[n].ref = "From HPLC in Sathyendranath et al., 2009 MEPS 383,73-84";
  parameters[n].index = n++;


  /* Zooplankton */
  parameters[n].name  = "ZSumax";
  parameters[n].desc  = "Maximum growth rate of ZS at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{max}^{ZS}";
  parameters[n].value[0] = 4.0;
  parameters[n].index = n++;

  parameters[n].name  = "ZSrad";
  parameters[n].desc  = "Radius of the small zooplankton cells";
  parameters[n].units = "m";
  parameters[n].sym   = "r_{ZS}";
  parameters[n].value[0] = 5.0e-06;
  parameters[n].index = n++;

  parameters[n].name  = "ZSswim";
  parameters[n].desc  = "Swimming velocity for small zooplankton";
  parameters[n].units = "m s-1";
  parameters[n].sym   = "U_{ZS}";
  parameters[n].value[0] = 2.0e-4;
  parameters[n].index = n++;

  parameters[n].name  = "ZSmeth";
  parameters[n].desc  = "Grazing technique of small zooplankton";
  parameters[n].units = "none";
  parameters[n].stringvalue = "rect";
  parameters[n].index = n++;

  parameters[n].name  = "ZLumax";
  parameters[n].desc  = "Maximum growth rate of ZL at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{max}^{ZL}";
  parameters[n].value[0] = 1.33;
  parameters[n].index = n++;

  parameters[n].name  = "ZLrad";
  parameters[n].desc  = "Radius of the large zooplankton cells";
  parameters[n].units = "m";
  parameters[n].sym   = "r_{ZL}";
  parameters[n].value[0] = 3.20e-04;
  parameters[n].index = n++;

  parameters[n].name  = "ZLswim";
  parameters[n].desc  = "Swimming velocity for large zooplankton";
  parameters[n].units = "m s-1";
  parameters[n].sym   = "U_{ZL}";
  parameters[n].value[0] = 3.0e-3;
  parameters[n].index = n++;

  parameters[n].name  = "ZLmeth";
  parameters[n].desc  = "Grazing technique of large zooplankton";
  parameters[n].units = "none";
  parameters[n].stringvalue = "rect";
  parameters[n].index = n++;

  parameters[n].name  = "ZL_E";
  parameters[n].desc  = "Growth efficiency, large zooplankton";
  parameters[n].units = "none";
  parameters[n].sym   = "E_{ZL}";
  parameters[n].value[0] = 0.426;
  parameters[n].stderr = 0.0179;
  parameters[n].ref = "Baird and Suthers, 2007 from Hansen et al (1997) LO 42: 687-704";
  parameters[n].index = n++;

  parameters[n].name  = "ZS_E";
  parameters[n].desc  = "Growth efficiency, small zooplankton";
  parameters[n].units = "none";
  parameters[n].sym   = "E_{ZS}";
  parameters[n].value[0] = 0.462;
  parameters[n].stderr = 0.0266;
  parameters[n].ref = "Baird and Suthers, 2007 from Hansen et al (1997) LO 42: 687-704";
  parameters[n].index = n++;

  parameters[n].name  = "ZL_mQ";
  parameters[n].desc  = "Natural (quadratic) mortality rate, large zooplankton";
  parameters[n].sym   = "m_{Q,ZL}";
  parameters[n].units = "d-1 (mg N m-3)-1";
  parameters[n].value[0] = 0.012;
  parameters[n].index = n++;

  parameters[n].name  = "ZS_mQ";
  parameters[n].desc  = "Natural (quadratic) mortality rate, small zooplankton";
  parameters[n].units = "d-1 (mg N m-3)-1";
  parameters[n].sym   = "m_{Q,ZS}";
  parameters[n].value[0] = 0.020;
  parameters[n].index = n++;

  parameters[n].name  = "ZL_FDG";
  parameters[n].desc  = "Fraction of growth inefficiency lost to detritus, large zoo.";
  parameters[n].units = "none";
  parameters[n].sym   = "\\gamma_{ZL}";
  parameters[n].value[0] = 0.5;
  parameters[n].index = n++;

  parameters[n].name  = "ZL_FDM";
  parameters[n].desc  = "Fraction of mortality lost to detritus, large zoo.";
  parameters[n].units = "none";
  parameters[n].sym   = "N/A";
  parameters[n].value[0] = 1.0;
  parameters[n].index = n++;

  parameters[n].name  = "ZS_FDG";
  parameters[n].desc  = "Fraction of growth inefficiency lost to detritus, small zoo.";
  parameters[n].units = "none";
  parameters[n].sym   = "\\gamma_{ZS}";
  parameters[n].value[0] = 0.5;
  parameters[n].index = n++;

  parameters[n].name  = "ZS_FDM";
  parameters[n].desc  = "Fraction of mortality lost to detritus, small zooplankton";
  parameters[n].units = "none";
  parameters[n].sym   = "N/A";
  parameters[n].value[0] = 1.0;
  parameters[n].index = n++;

  /* Remineralisation */
  parameters[n].name  = "F_LD_RD";
  parameters[n].desc  = "Fraction of labile detritus converted to refractory detritus";
  parameters[n].units = "none";
  parameters[n].sym   = "\\zeta_{Red}";
  parameters[n].value[0] = 0.19;
  parameters[n].index = n++;

  parameters[n].name  = "F_LD_DOM";
  parameters[n].desc  = "Fraction of labile detritus converted to DOM";
  parameters[n].units = "none";
  parameters[n].sym   = "\\vartheta_{Red}";
  parameters[n].value[0] = 0.1;
  parameters[n].index = n++;

  parameters[n].name  = "F_RD_DOM";
  parameters[n].desc  = "fraction of refractory detritus that breaks down to DOM";
  parameters[n].units = "none";
  parameters[n].sym   = "\\vartheta_{Ref}";
  parameters[n].value[0] = 0.05;
  parameters[n].index = n++;

  parameters[n].name  = "r_DetPL";
  parameters[n].desc  = "Breakdown rate of labile detritus at 106:16:1";
  parameters[n].units = "d-1";
  parameters[n].sym   = "r_{Red}";
  parameters[n].value[0] = 0.04;
  parameters[n].index = n++;

  parameters[n].name  = "r_DetBL";
  parameters[n].desc  = "Breakdown rate of labile detritus at 550:30:1";
  parameters[n].units = "d-1";
  parameters[n].sym   = "r_{Atk}";
  parameters[n].value[0] = 0.001;
  parameters[n].index = n++;

  parameters[n].name  = "r_RD";
  parameters[n].desc  = "Breakdown rate of refractory detritus";
  parameters[n].units = "d-1";
  parameters[n].sym   = "r_{R}";
  parameters[n].value[0] = 0.001;
  parameters[n].index = n++;

  parameters[n].name  = "r_DOM";
  parameters[n].desc  = "Breakdown rate of dissolved organic matter";
  parameters[n].units = "d-1";
  parameters[n].sym   = "r_{O}";
  parameters[n].value[0] = 0.0001;
  parameters[n].index = n++;

  parameters[n].name  = "Plank_resp";
  parameters[n].desc  = "Respiration as a fraction of umax";
  parameters[n].units = "none";
  parameters[n].sym   = "\\phi";
  parameters[n].value[0] = 0.025;
  parameters[n].index = n++;

  /* Sediment parameters */

  parameters[n].name  = "KO_aer";
  parameters[n].desc  = "Oxygen half-saturation for aerobic respiration";
  parameters[n].units = "mg O m-3";
  parameters[n].sym   = "K_{OA}";
  parameters[n].value[0] = 256.0;
  parameters[n].index = n++;

  parameters[n].name  = "r_nit_wc";
  parameters[n].desc  = "Maximum nitrification rate in water column";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\tau_{nit,wc}";
  parameters[n].value[0] = 0.1;
  parameters[n].index = n++;

  parameters[n].name  = "r_nit_sed";
  parameters[n].desc  = "Maximum nitrification rate in water sediment";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\tau_{nit,sed}";
  parameters[n].value[0] = 20.0;
  parameters[n].index = n++;

  parameters[n].name  = "KO_nit";
  parameters[n].desc  = "Oxygen half-saturation for nitrification";
  parameters[n].units = "mg O m-3";
  parameters[n].sym   = "K_{\\mathrm{O}_2,nit}";
  parameters[n].value[0] = 500.0;
  parameters[n].index = n++;

  parameters[n].name  = "Pads_r";
  parameters[n].desc  = "Rate at which P reaches adsorbed/desorbed equilibrium";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\tau_{Pabs}";
  parameters[n].value[0] = 0.04;
  parameters[n].index = n++;

  parameters[n].name  = "Pads_Kwc";
  parameters[n].desc  = "Freundlich Isothermic Const P adsorption to TSS in wc";
  parameters[n].units = "mg P kg TSS-1";
  parameters[n].sym   = "k_{Pads,wc}";
  parameters[n].value[0] = 30.0;
  parameters[n].index = n++;

  parameters[n].name  = "Pads_Ksed";
  parameters[n].desc  = "Freundlich Isothermic Const P adsorption to TSS in sed";
  parameters[n].units = "mg P kg TSS-1";
  parameters[n].sym   = "k_{Pads,sed}";
  parameters[n].value[0] = 74.0;
  parameters[n].index = n++;

  parameters[n].name  = "Pads_KO";
  parameters[n].desc  = "Oxygen half-saturation for P adsorption";
  parameters[n].units = "mg O m-3";
  parameters[n].sym   = "K_{\\mathrm{O}_2,abs}";
  parameters[n].value[0] = 2000.0;
  parameters[n].index = n++;

  parameters[n].name  = "Pads_exp";
  parameters[n].desc  = "Exponent for Freundlich Isotherm";
  parameters[n].units = "none";
  parameters[n].sym   = "N/A";
  parameters[n].value[0] = 1.0;
  parameters[n].index = n++;

  parameters[n].name  = "r_den";
  parameters[n].desc  = "Maximum denitrification rate";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\tau_{denit}";
  parameters[n].value[0] = 0.1;  
  parameters[n].index = n++;

  parameters[n].name  = "KO_den";
  parameters[n].desc  = "Oxygen half-saturation constant for denitrification";
  parameters[n].units = "mg O m-3";
  parameters[n].sym   = "K_{\\mathrm{O}_2,denit}";
  parameters[n].value[0] = 10000.0;
  parameters[n].index = n++;

  parameters[n].name  = "r_immob_PIP";
  parameters[n].desc  = "Rate of conversion of PIP to immobilised PIP";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\tau_{Pimm}";
  parameters[n].value[0] = 0.0012;
  parameters[n].index = n++;

  parameters[n].name  = "EpiDiffCoeff";
  parameters[n].desc  = "Sediment-water diffusion coefficient";
  parameters[n].units = "m2 s-1";
  parameters[n].sym   = "D";
  parameters[n].value[0] = 3e-7;
  parameters[n].index = n++;

  parameters[n].name  = "EpiDiffDz";
  parameters[n].desc  = "Thickness of diffusive layer";
  parameters[n].units = "m";
  parameters[n].sym   = "h";
  parameters[n].value[0] = 0.0065;
  parameters[n].index = n++;
  
  /* Marcroalgae */
  parameters[n].name  = "MAumax";
  parameters[n].desc  = "Maximum growth rate of MA at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{MA}^{max}";
  parameters[n].value[0] = 1.0;
  parameters[n].index = n++;

  parameters[n].name  = "MA_mL";
  parameters[n].desc  = "Natural (linear) mortality rate, macroalgae";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\zeta_{MA}";
  parameters[n].value[0] = 0.01;
  parameters[n].index = n++;

  parameters[n].name  = "MAleafden";
  parameters[n].desc  = "Nitrogen-specific leaf area of macroalgae";
  parameters[n].units = "m2 g N-1";
  parameters[n].sym   = "\\Omega_{MA}";
  parameters[n].value[0] = 1.0;
  parameters[n].index = n++;

  parameters[n].name  = "Benth_resp";
  parameters[n].desc  = "Respiration as a fraction of umax";
  parameters[n].units = "none";
  parameters[n].sym   = "\\phi";
  parameters[n].value[0] = 0.025;
  parameters[n].index = n++;

  /* Seagrass parameters - Zostera */
  parameters[n].name  = "SGumax";
  parameters[n].desc  = "Maximum growth rate of SG at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{SG}^{max}";
  parameters[n].value[0] = 0.4;
  parameters[n].ref = "x2 nighttime, x2 for roots.";
  parameters[n].index = n++;
  
  parameters[n].name  = "SG_KN";
  parameters[n].desc  = "Half-saturation of SG N uptake in SED";
  parameters[n].units = "mg N m-3";
  parameters[n].sym   = "K_{SG,N}";
  parameters[n].value[0] = 420.0;
  parameters[n].ref = "Lee and Dunton (1999) 1204-1215. Table 3 Zostera";
  parameters[n].index = n++;

  parameters[n].name  = "SG_KP";
  parameters[n].desc  = "Half-saturation of SG P uptake in SED";
  parameters[n].units = "mg P m-3";
  parameters[n].sym   = "K_{SG,P}";
  parameters[n].value[0] = 96.0;
  parameters[n].ref = "Gras et al. (2003) Aquatic Botany 76:299-315. Thalassia testudinum.";
  parameters[n].index = n++;

  parameters[n].name  = "SG_mL";
  parameters[n].desc  = "Natural (linear) mortality rate aboveground seagrass";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\zeta_{SG_A}";
  parameters[n].value[0] = 0.03;
  parameters[n].stderr = 0.001;
  parameters[n].ref = "Fourquean et al.( 2003) Chem. Ecol. 19: 373-390.Thalassia leaves with one component decay";
  parameters[n].index = n++;

  parameters[n].name  = "SGROOT_mL";
  parameters[n].desc  = "Natural (linear) mortality rate belowground seagrass";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\zeta_{SG_B}";
  parameters[n].value[0] = 0.004;
  parameters[n].stderr = 0.0002;
  parameters[n].ref = "Fourquean et al. (2003) Chem. Ecol. 19: 373-390. Thalassia roots with one component decay";
  parameters[n].index = n++;

  parameters[n].name  = "SGfrac";
  parameters[n].desc  = "Fraction (target) of SG biomass below-ground";
  parameters[n].units = "-";
  parameters[n].sym   = "f_{below,SG}";
  parameters[n].value[0] = 0.75;
  parameters[n].ref = "Babcock (2015) Zostera capricornii.";
  parameters[n].index = n++;

  parameters[n].name  = "SGtransrate";
  parameters[n].desc  = "Time scale for seagrass translocation";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\tau_{tran,SG}";
  parameters[n].value[0] = 0.0333;
  parameters[n].ref = "Loosely based on Zostera marine Kaldy et al., 2013 MEPS 487:27-39";
  parameters[n].index = n++;

  parameters[n].name  = "SGleafden";
  parameters[n].desc  = "Nitrogen-specific leaf area of seagrass";
  parameters[n].units = "m2 g N-1";
  parameters[n].sym   = "\\Omega_{SG}";
  parameters[n].value[0] = 1.5;
  parameters[n].ref = "Zostera capricornia: leaf dimensions Kemp et al (1987) Mar Ecol. Prog. Ser. 41:79-86.";
  parameters[n].index = n++;

  parameters[n].name  = "SGseedfrac";
  parameters[n].desc  = "Seagrass seed biomass as fraction of 63 % cover";
  parameters[n].units = "-";
  parameters[n].sym   = "f_{seed,SG}";
  parameters[n].value[0] = 0.01;
  parameters[n].ref = "No source";
  parameters[n].index = n++;

  parameters[n].name  = "SGorient";
  parameters[n].desc  = "Sine of nadir Zostera canopy bending angle";
  parameters[n].units = "-";
  parameters[n].sym   = "\\sin \\beta_{blade,SG}";
  parameters[n].value[0] = 0.5;
  parameters[n].ref = "No source";
  parameters[n].index = n++;

  parameters[n].name  = "SGmlr";
  parameters[n].desc  = "Compensation irradiance for Zostera";
  parameters[n].units = "mol m-2";
  parameters[n].sym   = "E_{comp,SG}";
  parameters[n].value[0] = 4.5;
  parameters[n].ref = "Chartrand (2012) Tech report.";
  parameters[n].index = n++;

  parameters[n].name  = "SGrootdepth";
  parameters[n].desc  = "Maximum depth for Zostera roots";
  parameters[n].units = "m";
  parameters[n].sym   = "z_{root,SG}";
  parameters[n].value[0] = -0.15;
  parameters[n].ref = "Roberts (1993) Aust. J. Mar. Fresh. Res. 44:85-100.";
  parameters[n].index = n++;

  parameters[n].name  = "SG_tau_critical";
  parameters[n].desc  = "Critical shear stress for SG loss";
  parameters[n].units = "N m-2";
  parameters[n].sym   = "\\tau_{SG,shear}";
  parameters[n].value[0] = 1.0;
  parameters[n].ref = "NESP project";
  parameters[n].index = n++;

  parameters[n].name  = "SG_tau_efold";
  parameters[n].desc  = "Time-scale for critical shear stress for SG loss";
  parameters[n].units = "s";
  parameters[n].sym   = "\\tau_{SG,time}";
  parameters[n].value[0] = 43200.0;
  parameters[n].ref = "NESP project";
  parameters[n].index = n++;

  /* Seagrass parameters - Halophila */
  parameters[n].name  = "SGHumax";
  parameters[n].desc  = "Maximum growth rate of SGH at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{SGH}^{max}";
  parameters[n].value[0] = 0.4;
  parameters[n].ref = "x2 nighttime, x2 for roots.";
  parameters[n].index = n++;

  parameters[n].name  = "SGH_KN";
  parameters[n].desc  = "Half-saturation of SGH N uptake in SED";
  parameters[n].units = "mg N m-3";
  parameters[n].sym   = "K_{SGH,N}";
  parameters[n].value[0] = 420.0;
  parameters[n].index = n++;

  parameters[n].name  = "SGH_KP";
  parameters[n].desc  = "Half-saturation of SGH P uptake in SED";
  parameters[n].units = "mg P m-3";
  parameters[n].sym   = "K_{SGH,P}";
  parameters[n].value[0] = 96.0;
  parameters[n].index = n++;

  parameters[n].name  = "SGHleafden";
  parameters[n].desc  = "Nitrogen-specific leaf area of SGH"; 
  parameters[n].units = "m2 g N-1";
  parameters[n].sym   = "\\Omega_{SGH}";
  parameters[n].value[0] = 1.9;
  parameters[n].ref = "Halophila ovalis: leaf dimensions from Vermaat et al. (1995)";
  parameters[n].index = n++;

  parameters[n].name  = "SGH_mL";
  parameters[n].desc  = "Natural (linear) mortality rate, aboveground SGH";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\zeta_{SGH_A}";
  parameters[n].value[0] = 0.06;
  parameters[n].stderr = 0.001;
  parameters[n].ref = "Fourquean et al.(2003) Chem. Ecol. 19: 373-390.Thalassia leaves with one component decay";
  parameters[n].index = n++;

  parameters[n].name  = "SGHROOT_mL";
  parameters[n].desc  = "Natural (linear) mortality rate, belowground SGH";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\zeta_{SGH_B}";
  parameters[n].value[0] = 0.004;
  parameters[n].stderr = 0.0002;
  parameters[n].ref = "Fourquean et al. (2003) Chem. Ecol. 19: 373-390. Thalassia roots with one component decay";
  parameters[n].index = n++;

  parameters[n].name  = "SGHfrac";
  parameters[n].desc  = "Fraction (target) of SGH biomass below-ground";
  parameters[n].units = "-";
  parameters[n].sym   = "f_{below,SGH}";
  parameters[n].value[0] = 0.5;
  parameters[n].ref = "Babcock 2015, Halophila ovalis";
  parameters[n].index = n++;

  parameters[n].name  = "SGHtransrate";
  parameters[n].desc  = "Time scale for Halophila translocation";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\tau_{tran,SGH}";
  parameters[n].value[0] = 0.0333;
  parameters[n].ref = "Loosely based on Zostera marine Kaldy et al., 2013 MEPS 487:27-39";
  parameters[n].index = n++;

  parameters[n].name  = "SGHseedfrac";
  parameters[n].desc  = "Halophila seed biomass as fraction of 63 % cover";
  parameters[n].units = "-";
  parameters[n].sym   = "f_{seed,SGH}";
  parameters[n].value[0] = 0.01;
  parameters[n].ref = "No source";
  parameters[n].index = n++;

  parameters[n].name  = "SGHorient";
  parameters[n].desc  = "Sine of nadir Halophila canopy bending angle";
  parameters[n].units = "-";
  parameters[n].sym   = "\\sin \\beta_{blade,SGH}";
  parameters[n].value[0] = 1.0;
  parameters[n].ref = "No source";
  parameters[n].index = n++;

  parameters[n].name  = "SGHmlr";
  parameters[n].desc  = "Compensation irradiance for Halophila";
  parameters[n].units = "mol m-2";
  parameters[n].sym   = "E_{comp,SGH}";
  parameters[n].value[0] = 2.0;
  parameters[n].ref = "Longstaff 2003 UQ PhD thesis";
  parameters[n].index = n++;

  parameters[n].name  = "SGHrootdepth";
  parameters[n].desc  = "Maximum depth for Halophila roots";
  parameters[n].units = "m";
  parameters[n].sym   = "z_{root,SGH}";
  parameters[n].value[0] = -0.08;
  parameters[n].ref = "Roberts (1993) Aust. J. Mar. Fresh. Res. 44:85-100.";
  parameters[n].index = n++;

  parameters[n].name  = "SGH_tau_critical";
  parameters[n].desc  = "Critical shear stress for SGH loss";
  parameters[n].units = "N m-2";
  parameters[n].sym   = "\\tau_{SGH,shear}";
  parameters[n].value[0] = 1.0;
  parameters[n].ref = "NESP project";
  parameters[n].index = n++;

  parameters[n].name  = "SGH_tau_efold";
  parameters[n].desc  = "Time-scale for critical shear stress for SGH loss";
  parameters[n].units = "s";
  parameters[n].sym   = "\\tau_{SGH,time}";
  parameters[n].value[0] = 43200.0;
  parameters[n].ref = "NESP project";
  parameters[n].index = n++;


  /* Seagrass parameters - Deep */
  parameters[n].name  = "SGDumax";
  parameters[n].desc  = "Maximum growth rate of SGD at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{SGD}^{max}";
  parameters[n].value[0] = 0.4;
  parameters[n].ref = "x2 nighttime, x2 for roots.";
  parameters[n].index = n++;

  parameters[n].name  = "SGD_KN";
  parameters[n].desc  = "Half-saturation of SGD N uptake in SED";
  parameters[n].units = "mg N m-3";
  parameters[n].sym   = "K_{SGD,N}";
  parameters[n].value[0] = 420.0;
  parameters[n].index = n++;

  parameters[n].name  = "SGD_KP";
  parameters[n].desc  = "Half-saturation of SGD P uptake in SED";
  parameters[n].units = "mg P m-3";
  parameters[n].sym   = "K_{SGD,P}";
  parameters[n].value[0] = 96.0;
  parameters[n].index = n++;

  parameters[n].name  = "SGDleafden";
  parameters[n].desc  = "Nitrogen-specific leaf area of SGD"; 
  parameters[n].units = "m2 g N-1";
  parameters[n].sym   = "\\Omega_{SGD}";
  parameters[n].value[0] = 1.9;
  parameters[n].ref = "Halophila ovalis: leaf dimensions from Vermaat et al. (1995)";
  parameters[n].index = n++;

  parameters[n].name  = "SGD_mL";
  parameters[n].desc  = "Natural (linear) mortality rate, aboveground SGD";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\zeta_{SGD_A}";
  parameters[n].value[0] = 0.06;
  parameters[n].stderr = 0.001;
  parameters[n].ref = "NESP project";
  parameters[n].index = n++;

  parameters[n].name  = "SGDROOT_mL";
  parameters[n].desc  = "Natural (linear) mortality rate, belowground SGD";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\zeta_{SGD_B}";
  parameters[n].value[0] = 0.004;
  parameters[n].stderr = 0.00002;
  parameters[n].ref = "NESP project";
  parameters[n].index = n++;

  parameters[n].name  = "SGDfrac";
  parameters[n].desc  = "Fraction (target) of SGD biomass below-ground";
  parameters[n].units = "-";
  parameters[n].sym   = "f_{below,SGD}";
  parameters[n].value[0] = 0.25;
  parameters[n].ref = "Duarte (1999) Aquatic Biol. 65: 159-174, Halophila ovalis.";
  parameters[n].index = n++;

  parameters[n].name  = "SGDtransrate";
  parameters[n].desc  = "Time scale for deep SG translocation";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\tau_{tran,SGD}";
  parameters[n].value[0] = 0.0333;
  parameters[n].ref = "Loosely based on Zostera marine Kaldy et al., 2013 MEPS 487:27-39";
  parameters[n].index = n++;

  parameters[n].name  = "SGDseedfrac";
  parameters[n].desc  = "Deep SG seed biomass as fraction of 63 % cover";
  parameters[n].units = "-";
  parameters[n].sym   = "f_{seed,SGD}";
  parameters[n].value[0] = 0.01;
  parameters[n].ref = "No source";
  parameters[n].index = n++;

  parameters[n].name  = "SGDorient";
  parameters[n].desc  = "Sine of nadir deep SG canopy bending angle";
  parameters[n].units = "-";
  parameters[n].sym   = "\\sin \\beta_{blade,SGD}";
  parameters[n].value[0] = 1.0;
  parameters[n].ref = "No source";
  parameters[n].index = n++;

  parameters[n].name  = "SGDmlr";
  parameters[n].desc  = "Compensation irradiance for deep SG";
  parameters[n].units = "mol m-2";
  parameters[n].sym   = "E_{comp,SGD}";
  parameters[n].value[0] = 1.5;
  parameters[n].ref = "Chartrand (2017) Tech report.";
  parameters[n].index = n++;

  parameters[n].name  = "SGDrootdepth";
  parameters[n].desc  = "Maximum depth for deep SG roots";
  parameters[n].units = "m";
  parameters[n].sym   = "z_{root,SGD}";
  parameters[n].value[0] = -0.05;
  parameters[n].ref = "NESP project";
  parameters[n].index = n++;

  parameters[n].name  = "SGD_tau_critical";
  parameters[n].desc  = "Critical shear stress for deep SG loss";
  parameters[n].units = "N m-2";
  parameters[n].sym   = "\\tau_{SGD,shear}";
  parameters[n].value[0] = 1.0;
  parameters[n].ref = "NESP project";
  parameters[n].index = n++;

  parameters[n].name  = "SGD_tau_efold";
  parameters[n].desc  = "Time-scale for shear stress for deep SG loss";
  parameters[n].units = "s";
  parameters[n].sym   = "\\tau_{SGD,time}";
  parameters[n].value[0] = 43200.0;
  parameters[n].ref = "NESP project";
  parameters[n].index = n++;

  /* Corals */
  parameters[n].name  = "dissCaCO3_sed";
  parameters[n].desc  = "net dissolution rate of sediment without coral";
  parameters[n].units = "mmol C m-2 s-1";
  parameters[n].sym   = "d_{sand}";
  parameters[n].value[0] = 0.001;   // 0.007 in B1p9
  parameters[n].index = n++;
  parameters[n].ref = "Anthony et al. (2013), Biogeosciences 10:4897-4909, Fig 5E: -1 2 3 6  mmol m-2 h-1";

  parameters[n].name  = "CHarea";
  parameters[n].desc  = "Grid scale to reef scale ratio";
  parameters[n].units = "m2 m-2";
  parameters[n].sym   = "A_{CH}";
  parameters[n].value[0] = 0.1;
  parameters[n].index = n++;
  parameters[n].ref = "Heron Island for on 4 km model.";

  parameters[n].name  = "CHpolypden";
  parameters[n].desc  = "Nitrogen-specific host area of coral polyp";
  parameters[n].units = "m2 g N-1";
  parameters[n].sym   = "\\Omega_{CH}";
  parameters[n].value[0] = 2.0;
  parameters[n].index = n++;

  parameters[n].name  = "CHumax";
  parameters[n].desc  = "Max. growth rate of Coral at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{CH}^{max}";
  parameters[n].value[0] = 0.05;
  parameters[n].index = n++;

  parameters[n].name  = "CSumax";
  parameters[n].desc  = "Max. growth rate of zooxanthellae at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{CS}^{max}";
  parameters[n].value[0] = 0.4;
  parameters[n].index = n++;

  parameters[n].name  = "CSrad";
  parameters[n].desc  = "Radius of the zooxanthellae ";
  parameters[n].units = "m";
  parameters[n].sym   = "r_{CS}";
  parameters[n].value[0] = 5e-06;
  parameters[n].index = n++;

  parameters[n].name  = "CHmort";
  parameters[n].desc  = "Quadratic mortality rate of coral polyp ";
  parameters[n].units = "(g N m-2)-1 d-1";
  parameters[n].sym   = "\\zeta_{CH}";
  parameters[n].value[0] = 0.01;
  parameters[n].index = n++;

  parameters[n].name  = "CSmort";
  parameters[n].desc  = "Linear mortality rate of zooxanthellae ";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\zeta_{CS}";
  parameters[n].value[0] = 0.04;
  parameters[n].index = n++;

  parameters[n].name  = "CHremin";
  parameters[n].desc  = "Fraction of coral host death translocated.";
  parameters[n].units = "-";
  parameters[n].sym   = "f_{remin}";
  parameters[n].value[0] = 0.5;
  parameters[n].index = n++;

  parameters[n].name  = "Splank";
  parameters[n].desc  = "Rate coefficent for particle uptake by corals";
  parameters[n].units = "m d-1";
  parameters[n].sym   = "S_{part}";
  parameters[n].value[0] = 3.0;
  parameters[n].ref = "Ribes and Atkinson (2007) Coral Reefs 26: 413-421";
  parameters[n].index = n++;

  parameters[n].name  = "k_day_coral";
  parameters[n].desc  = "Maximum daytime coral calcification";
  parameters[n].units = "mmol C m-2 s-1";
  parameters[n].sym   = "k_{day}";
  parameters[n].value[0] = 0.0132;
  parameters[n].ref = "Anthony et al. (2013), Biogeosciences 10:4897-4909, Fig 5A: 50, 50, 35 55 mmol m-2 h-1 for Acropora aspera n=4";
  parameters[n].index = n++;

  parameters[n].name  = "k_night_coral";
  parameters[n].desc  = "Maximum nightime coral calcification";
  parameters[n].units = "mmol C m-2 s-1";
  parameters[n].sym   = "k_{night}";
  parameters[n].value[0] = 0.0069;
  parameters[n].ref = "Anthony et al. (2013), Biogeosciences 10:4897-4909, Fig 5A: 20, 30, 20, 30  mmol m-2 h-1 for Acropora aspera n=4";
  parameters[n].index = n++;

  parameters[n].name  = "dissCaCO3_shelf";
  parameters[n].desc  = "Carbonate sediment dissolution rate on shelf";
  parameters[n].units = "mmol C m-2 s-1";
  parameters[n].sym   = "d_{shelf}";
  parameters[n].value[0] = 0.0001;
  parameters[n].ref = "Cyronak, T. et al., LO 58:131-143. Heron Island study.";
  parameters[n].index = n++;

  parameters[n].name  = "ageing_decay";
  parameters[n].desc  = "Age tracer growth rate per day";
  parameters[n].units = "d d-1";
  parameters[n].sym   = "n/a";
  parameters[n].value[0] = 1.0;
  parameters[n].ref = "EMS manual";
  parameters[n].index = n++;

  parameters[n].name  = "anti_ageing_decay";
  parameters[n].desc  = "Age tracer decay rate per day outside source";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\Phi";
  parameters[n].value[0] = 0.1;
  parameters[n].ref = "EMS manual";
  parameters[n].index = n++;

  /* New parameters for B3p1 */

  parameters[n].name  = "ROSthreshold";
  parameters[n].desc  = "Bleaching ROS threshold per cell";
  parameters[n].units = "mg ROS cell-1";
  parameters[n].sym   = "\\phi_{ROS}";
  parameters[n].value[0] = 1.418e-14;
  parameters[n].ref = "EMS manual";
  parameters[n].index = n++;

  parameters[n].name  = "Xanth_tau";
  parameters[n].desc  = "Xanthophyll switching rate coefficient";
  parameters[n].units = "s-1";
  parameters[n].sym   = "\\tau_{xan}";
  parameters[n].value[0] = 8.333333e-04;
  parameters[n].ref = "Gustafsson et al., 2013";
  parameters[n].index = n++;

  parameters[n].name  = "chla2rcii";
  parameters[n].desc  = "Ratio of RCII to Chlorophyll a";
  parameters[n].units = "mol RCII g Chl-1";
  parameters[n].sym   = "A_{RCII}";
  parameters[n].value[0] = 2.238413e-06;
  parameters[n].ref = "Suggett et al., 2009";
  parameters[n].index = n++;

  parameters[n].name  = "photon2rcii";
  parameters[n].desc  = "Stoichiometric ratio of RCII units to photons";
  parameters[n].units = "mol RCII mol photon-1";
  parameters[n].sym   = "m_{RCII}";
  parameters[n].value[0] = 0.1e-6;
  parameters[n].ref = "";
  parameters[n].index = n++;

  parameters[n].name  = "r_RD_NtoP";
  parameters[n].desc  = "Scaling of DetP to DOP, relative to N";
  parameters[n].units = "-";
  parameters[n].sym   = "r_{RD_NtoP}";
  parameters[n].value[0] = 2.0;
  parameters[n].ref = "EMS manual";
  parameters[n].index = n++;

  parameters[n].name  = "r_DOM_NtoP";
  parameters[n].desc  = "Scaling of DOM to DIP, relative to N";
  parameters[n].units = "-";
  parameters[n].sym   = "r_{DOM_NtoP}";
  parameters[n].value[0] = 1.5;
  parameters[n].ref = "EMS manual";
  parameters[n].index = n++;

  parameters[n].name  = "Tricho_colrad";
  parameters[n].desc  = "Radius of Trichodesmium colonies";
  parameters[n].units = "m";
  parameters[n].sym   = "r_{Tricho colony}";
  parameters[n].value[0] = 0.000005;
  parameters[n].ref = "EMS manual";
  parameters[n].index = n++;

  parameters[n].name  = "CS_photon2ros";
  parameters[n].desc  = "Stoichiometric coefficient of ROS";
  parameters[n].units = "-";
  parameters[n].sym   = "-";
  parameters[n].value[0] = 7.0e7;
  parameters[n].ref = "EMS manual";
  parameters[n].index = n++;

  parameters[n].name  = "ROSmult";
  parameters[n].desc  = "Linear coefficient of bleaching for above threshold fraction";
  parameters[n].units = "unitless";
  parameters[n].sym   = "k_{CS,ROSfrac}";
  parameters[n].value[0] = 1.0;
  parameters[n].ref = "EMS manual";
  parameters[n].index = n++;

  parameters[n].name  = "CSmaxbleachrate";
  parameters[n].desc  = "Maximum coral bleaching rate";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\zeta_{bleach}";
  parameters[n].value[0] = 1.0;
  parameters[n].ref = "EMS manual";
  parameters[n].index = n++;

  /* Assign acutal number of parameters */
  *nprm = n;

  /* Assign string values */
  assign_string_values(parameters, *nprm);
}


void eco_params_gbr4(parameter_info **params, int *nprm)
{
  int n;
  parameter_info *parameters = *params;

  init_parameter_values(parameters);

  n = 0;

  /* Temperature / Environment parameters */
  parameters[n].name  = "Tref";
  parameters[n].desc  = "Reference temperature";
  parameters[n].sym   = "T_{ref}";
  parameters[n].units = "Deg C";
  parameters[n].value[0] = 20.0;
  parameters[n].ref = " ";
  parameters[n].index = n++;

  parameters[n].name  = "Q10";
  parameters[n].desc  = "Temperature coefficient for rate parameters";
  parameters[n].sym   = "Q10";
  parameters[n].units = "none";
  parameters[n].value[0] = 2.0;
  parameters[n].ref = " ";
  parameters[n].index = n++;

  parameters[n].name  = "TKEeps";
  parameters[n].desc  = "Nominal rate of TKE dissipation in water column";
  parameters[n].sym   = "\\epsilon";
  parameters[n].units = "m2 s-3";
  parameters[n].value[0] = 1.0e-6;
  parameters[n].index = n++;

  parameters[n].name  = "xco2_in_air";
  parameters[n].desc  = "Atmospheric CO2";
  parameters[n].sym   = "p\\mathrm{CO}_2";
  parameters[n].units = "ppmv";
  parameters[n].value[0] = 396.48;
  parameters[n].ref = "Mean 2013 at Mauna Loa: http://co2now.org/current-co2/co2-now/";
  parameters[n].index = n++;

  parameters[n].name  = "N2";
  parameters[n].desc  = "Concentration of dissolved N2";
  parameters[n].units = "mg N m-3";
  parameters[n].sym   = "[\\mathrm{N}_2]_{gas}";
  parameters[n].value[0] = 2000.0;
  parameters[n].ref = "Robson et al. (2013)";
  parameters[n].index = n++;

  /* Light parameters */

  parameters[n].name = "Light_lambda";
  parameters[n].desc = "Wavelengths of light";
  parameters[n].units= "nm";
  /* Free array allocated above */
  d_free_1d(parameters[n].value);
  parameters[n].num_values = 24;
  parameters[n].value = d_alloc_1d(24);
  parameters[n].value[0] = 290.0;    
  parameters[n].value[1] = 310.0;
  parameters[n].value[2] = 330.0;
  parameters[n].value[3] = 350.0;
  parameters[n].value[4] = 370.0;
  parameters[n].value[5] = 390.0;
  parameters[n].value[6] = 410.0;
  parameters[n].value[7] = 430.0;
  parameters[n].value[8] = 440.0;
  parameters[n].value[9] = 450.0;
  parameters[n].value[10] = 470.0;
  parameters[n].value[11] = 490.0;
  parameters[n].value[12] = 510.0;
  parameters[n].value[13] = 530.0;
  parameters[n].value[14] = 550.0;
  parameters[n].value[15] = 570.0;
  parameters[n].value[16] = 590.0;
  parameters[n].value[17] = 610.0;
  parameters[n].value[18] = 630.0;
  parameters[n].value[19] = 650.0;
  parameters[n].value[20] = 670.0;
  parameters[n].value[21] = 690.0;
  parameters[n].value[22] = 710.0;
  parameters[n].value[23] = 800.0;
  parameters[n].ref = "Approx. 20 nm resolution with 10 nm about 440 nm. PAR (400-700) is integral of bands 6-22.";
  parameters[n].index = n++;

  parameters[n].name  = "acdom443star";
  parameters[n].desc  = "DOC-specific absorption of CDOM 443 nm";
  parameters[n].units = "m2 mg C-1";
  parameters[n].sym   = "k_{CDOM,443}";
  parameters[n].value[0] = 0.00013;
  parameters[n].ref = "Based on Feb 2011 GBR satellite data and modelled DOR_C";
  parameters[n].index = n++;

  /* Phytoplankton */

  parameters[n].name  = "PLumax";
  parameters[n].desc  = "Maximum growth rate of PL at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{PL}^{max}";
  parameters[n].value[0] = 1.4;
  parameters[n].ref = " ";
  parameters[n].index = n++;

  parameters[n].name  = "PLrad";
  parameters[n].desc  = "Radius of the large phytoplankton cells";
  parameters[n].units = "m";
  parameters[n].sym   = "r_{PL}";
  parameters[n].value[0] = 4e-06;
  parameters[n].index = n++;

  parameters[n].name  = "PhyL_mL";
  parameters[n].desc  = "Natural (linear) mortality rate, large phytoplankton";
  parameters[n].units = "d-1";
  parameters[n].sym   = "m_{L,PL}";
  parameters[n].value[0] = 0.1;
  parameters[n].index = n++;

  parameters[n].name  = "PhyL_mL_sed";
  parameters[n].desc  = "Natural (linear) mortality rate in sediment, large phytoplankton";
  parameters[n].units = "d-1";
  parameters[n].sym   = "m_{L,PL,sed}";
  parameters[n].value[0] = 1.0;
  parameters[n].index = n++;

  parameters[n].name  = "PSumax";
  parameters[n].desc  = "Maximum growth rate of PS at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{PL}^{max}";
  parameters[n].value[0] = 1.6;
  parameters[n].index = n++;

  parameters[n].name  = "PSrad";
  parameters[n].desc  = "Radius of the small phytoplankton cells";
  parameters[n].units = "m";
  parameters[n].sym   = "r_{PS}";
  parameters[n].value[0] = 1.e-06;
  parameters[n].index = n++;

  parameters[n].name  = "PhyS_mL";
  parameters[n].desc  = "Natural (linear) mortality rate, small phytoplankton";
  parameters[n].units = "d-1";
  parameters[n].sym   = "m_{L,PS}";
  parameters[n].value[0] = 0.1;
  parameters[n].index = n++;

  parameters[n].name  = "PhyS_mL_sed";
  parameters[n].desc  = "Natural (linear) mortality rate in sediment, small phytoplankton";
  parameters[n].units = "d-1";
  parameters[n].sym   = "m_{L,PS,sed}";
  parameters[n].value[0] = 1.0;
  parameters[n].index = n++;

  parameters[n].name  = "MBumax";
  parameters[n].desc  = "Maximum growth rate of MB at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{MPB}^{max}";
  parameters[n].value[0] = 0.839;
  parameters[n].index = n++;

  parameters[n].name  = "MBrad";
  parameters[n].desc  = "Radius of the MPB cells";
  parameters[n].units = "m";
  parameters[n].sym   = "r_{MPB}";
  parameters[n].value[0] = 1e-05;
  parameters[n].index = n++;

  parameters[n].name  = "MPB_mQ";
  parameters[n].desc  = "Natural (quadratic) mortality rate, microphytobenthos, applied in sediment";
  parameters[n].units = "d-1 (mg N m-3)-1";
  parameters[n].sym   = "m_{Q,MPB}";
  parameters[n].value[0] = 0.0001;
  parameters[n].ref = "At steady-state, at mu = 0.1 d-1, indep. of temp, MPB_N ~ 0.1 / MPB_mQ = 250 mg N m-3";
  parameters[n].index = n++;

  parameters[n].name  = "PSxan2chl";
  parameters[n].desc  = "Ratio of xanthophyll to chl a of PS";
  parameters[n].units = "mg mg-1";
  parameters[n].sym   = "\\Theta_{xan2chl,PS}";
  parameters[n].value[0] = 0.51;
  parameters[n].ref = "CSIRO parameter library: GBR region";
  parameters[n].index = n++;

  parameters[n].name  = "PLxan2chl";
  parameters[n].desc  = "Ratio of xanthophyll to chl a of PL";
  parameters[n].units = "mg mg-1";
  parameters[n].sym   = "\\Theta_{xan2chl,PL}";
  parameters[n].value[0] = 0.81;
  parameters[n].ref = "CSIRO parameter library: GBR region";
  parameters[n].index = n++;

  parameters[n].name  = "MBxan2chl";
  parameters[n].desc  = "Ratio of xanthophyll to chl a of MPB";
  parameters[n].units = "mg mg-1";
  parameters[n].sym   = "\\Theta_{xan2chl,MPB}";
  parameters[n].value[0] = 0.81;
  parameters[n].ref = "CSIRO parameter library: GBR region WC values";
  parameters[n].index = n++;

  /* Trichodesmium */
  parameters[n].name  = "Tricho_umax";
  parameters[n].desc  = "Maximum growth rate of Trichodesmium at Tref ";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{MPB}^{max}";
  parameters[n].value[0] = 0.24;
  parameters[n].index = n++;

  parameters[n].name  = "Tricho_rad";
  parameters[n].desc  = "Radius of Trichodesmium colonies";
  parameters[n].units = "m";
  parameters[n].sym   = "r_{MPB}";
  parameters[n].value[0] = 0.000005;
  parameters[n].index = n++;

  parameters[n].name  = "Tricho_Sh";
  parameters[n].desc  = "Sherwood number for the Tricho dimensionless";
  parameters[n].units = "none";
  parameters[n].sym   = "Sh_{Tricho}";
  parameters[n].value[0] = 1.0;
  parameters[n].index = n++;

  parameters[n].name  = "Tricho_mL";
  parameters[n].desc  = "Linear mortality for Tricho in sediment";
  parameters[n].units = "d-1";
  parameters[n].sym   = "m_{L,Tricho}";
  parameters[n].value[0] = 0.140000;
  parameters[n].index = n++;

  parameters[n].name  = "Tricho_mQ";
  parameters[n].desc  = "Quadratic mortality for Tricho due to phages in wc";
  parameters[n].units = "d-1 (mg N m-3)-1";
  parameters[n].sym   = "m_{Q,Tricho}";
  parameters[n].value[0] = 0.2000;
  parameters[n].ref = "At steady-state, indep. of temp, Tricho_N ~ Tricho_umax / Tricho_mQ = 1.0 / 0.15 = 6.67 mg N m-3";
  parameters[n].index = n++;

  parameters[n].name  = "Tricho_crit";
  parameters[n].desc  = "Critical Tricho above which quadratic mortality applies";
  parameters[n].units = "mg N m-3";
  parameters[n].value[0] = 0.0002000;
  parameters[n].ref = "No source";
  parameters[n].index = n++;

  parameters[n].name  = "p_min";
  parameters[n].desc  = "Minimum density of Trichodesmium";
  parameters[n].units = "kg m-3";
  parameters[n].sym   = "\\rho_{min,Tricho}";
  parameters[n].value[0] = 990.000000;
  parameters[n].ref = "No source";
  parameters[n].index = n++;

  parameters[n].name  = "p_max";
  parameters[n].desc  = "Maximum density of Trichodesmium";
  parameters[n].units = "kg m-3";
  parameters[n].sym   = "\\rho_{max,Tricho}";
  parameters[n].value[0] = 1060.000000;
  parameters[n].ref = "No source";
  parameters[n].index = n++;

  parameters[n].name  = "DINcrit";
  parameters[n].desc  = "DIN conc below which Tricho N fixes ";
  parameters[n].units = "mg N m-3";
  parameters[n].sym   = "DIN_{crit}";
  parameters[n].value[0] = 10.0;
  parameters[n].ref = "Lower end of Robson et al., (2013) 4-20 mg N m-3";
  parameters[n].index = n++;

  parameters[n].name  = "Trichoxan2chl";
  parameters[n].desc  = "Ratio of xanthophyll to chl a of Trichodesmium";
  parameters[n].units = "mg mg-1";
  parameters[n].sym   = "\\Theta_{xan2chl,Tricho}";
  parameters[n].value[0] = 0.5;
  parameters[n].ref = "Subramaniam (1999) LO 44:608-617. Actually redder pigment than xanthophyll.";
  parameters[n].index = n++;

  parameters[n].name  = "bphy";
  parameters[n].desc  = "Chl-specific scattering coef. for microalgae";
  parameters[n].units = "m-1 (mg Chl a m-3)-1"; 
  parameters[n].sym   = "b_{phy}";
  parameters[n].value[0] = 0.2;
  parameters[n].ref = "Typical microalgae value, Kirk (1994) Light and Photosynthesis in Aquatic ecosystems, Table 4.3";
  parameters[n].index = n++;

  parameters[n].name  = "NtoCHL";
  parameters[n].desc  = "Nominal N:Chl a ratio in phytoplankton by weight";
  parameters[n].units = "g N (g Chl a)-1";
  parameters[n].sym   = "R_{N:Chl}";
  parameters[n].value[0] = 7;
  parameters[n].ref = "Represents a C:Chl ratio of 39.25, Baird et al. (2013) Limnol. Oceanogr. 58: 1215-1226.";
  parameters[n].index = n++;

  parameters[n].name  = "C2Chlmin";
  parameters[n].desc  = "Minimum carbon to chlorophyll a ratio";
  parameters[n].units = "wt/wt";
  parameters[n].sym   = "\\theta_{min}";
  parameters[n].value[0] = 20.0;
  parameters[n].ref = "From HPLC in Sathyendranath et al., 2009 MEPS 383,73-84";
  parameters[n].index = n++;


  /* Zooplankton */
  parameters[n].name  = "ZSumax";
  parameters[n].desc  = "Maximum growth rate of ZS at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{max}^{ZS}";
  parameters[n].value[0] = 4.0;
  parameters[n].index = n++;

  parameters[n].name  = "ZSrad";
  parameters[n].desc  = "Radius of the small zooplankton cells";
  parameters[n].units = "m";
  parameters[n].sym   = "r_{ZS}";
  parameters[n].value[0] = 5.0e-06;
  parameters[n].index = n++;

  parameters[n].name  = "ZSswim";
  parameters[n].desc  = "Swimming velocity for small zooplankton";
  parameters[n].units = "m s-1";
  parameters[n].sym   = "U_{ZS}";
  parameters[n].value[0] = 2.0e-4;
  parameters[n].index = n++;

  parameters[n].name  = "ZSmeth";
  parameters[n].desc  = "Grazing technique of small zooplankton";
  parameters[n].units = "none";
  parameters[n].stringvalue = "rect";
  parameters[n].index = n++;

  parameters[n].name  = "ZLumax";
  parameters[n].desc  = "Maximum growth rate of ZL at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{max}^{ZL}";
  parameters[n].value[0] = 1.33;
  parameters[n].index = n++;

  parameters[n].name  = "ZLrad";
  parameters[n].desc  = "Radius of the large zooplankton cells";
  parameters[n].units = "m";
  parameters[n].sym   = "r_{ZL}";
  parameters[n].value[0] = 3.20e-04;
  parameters[n].index = n++;

  parameters[n].name  = "ZLswim";
  parameters[n].desc  = "Swimming velocity for large zooplankton";
  parameters[n].units = "m s-1";
  parameters[n].sym   = "U_{ZL}";
  parameters[n].value[0] = 3.0e-3;
  parameters[n].index = n++;

  parameters[n].name  = "ZLmeth";
  parameters[n].desc  = "Grazing technique of large zooplankton";
  parameters[n].units = "none";
  parameters[n].stringvalue = "rect";
  parameters[n].index = n++;

  parameters[n].name  = "ZL_E";
  parameters[n].desc  = "Growth efficiency, large zooplankton";
  parameters[n].units = "none";
  parameters[n].sym   = "E_{ZL}";
  parameters[n].value[0] = 0.426;
  parameters[n].stderr = 0.0179;
  parameters[n].ref = "Baird and Suthers, 2007 from Hansen et al (1997) LO 42: 687-704";
  parameters[n].index = n++;

  parameters[n].name  = "ZS_E";
  parameters[n].desc  = "Growth efficiency, small zooplankton";
  parameters[n].units = "none";
  parameters[n].sym   = "E_{ZS}";
  parameters[n].value[0] = 0.462;
  parameters[n].stderr = 0.0266;
  parameters[n].ref = "Baird and Suthers, 2007 from Hansen et al (1997) LO 42: 687-704";
  parameters[n].index = n++;

  parameters[n].name  = "ZL_mQ";
  parameters[n].desc  = "Natural (quadratic) mortality rate, large zooplankton";
  parameters[n].sym   = "m_{Q,ZL}";
  parameters[n].units = "d-1 (mg N m-3)-1";
  parameters[n].value[0] = 0.012;
  parameters[n].index = n++;

  parameters[n].name  = "ZS_mQ";
  parameters[n].desc  = "Natural (quadratic) mortality rate, small zooplankton";
  parameters[n].units = "d-1 (mg N m-3)-1";
  parameters[n].sym   = "m_{Q,ZS}";
  parameters[n].value[0] = 0.007;
  parameters[n].index = n++;

  parameters[n].name  = "ZL_FDG";
  parameters[n].desc  = "Fraction of growth inefficiency lost to detritus, large zooplankton";
  parameters[n].units = "none";
  parameters[n].sym   = "\\gamma_{ZL}";
  parameters[n].value[0] = 0.5;
  parameters[n].index = n++;

  parameters[n].name  = "ZL_FDM";
  parameters[n].desc  = "Fraction of mortality lost to detritus, large zooplankton";
  parameters[n].units = "none";
  parameters[n].sym   = "N/A";
  parameters[n].value[0] = 1.0;
  parameters[n].index = n++;

  parameters[n].name  = "ZS_FDG";
  parameters[n].desc  = "Fraction of growth inefficiency lost to detritus, small zooplankton";
  parameters[n].units = "none";
  parameters[n].sym   = "\\gamma_{ZS}";
  parameters[n].value[0] = 0.5;
  parameters[n].index = n++;

  parameters[n].name  = "ZS_FDM";
  parameters[n].desc  = "Fraction of mortality lost to detritus, small zooplankton";
  parameters[n].units = "none";
  parameters[n].sym   = "N/A";
  parameters[n].value[0] = 1.0;
  parameters[n].index = n++;

  /* Remineralisation */
  parameters[n].name  = "F_LD_RD";
  parameters[n].desc  = "Fraction of labile detritus converted to refractory detritus";
  parameters[n].units = "none";
  parameters[n].sym   = "\\zeta_{Red}";
  parameters[n].value[0] = 0.19;
  parameters[n].index = n++;

  parameters[n].name  = "F_LD_DOM";
  parameters[n].desc  = "Fraction of labile detritus converted to dissolved organic matter";
  parameters[n].units = "none";
  parameters[n].sym   = "\\vartheta_{Red}";
  parameters[n].value[0] = 0.1;
  parameters[n].index = n++;

  parameters[n].name  = "F_RD_DOM";
  parameters[n].desc  = "fraction of refractory detritus that breaks down to DOM";
  parameters[n].units = "none";
  parameters[n].sym   = "\\vartheta_{Ref}";
  parameters[n].value[0] = 0.05;
  parameters[n].index = n++;

  parameters[n].name  = "r_DetPL";
  parameters[n].desc  = "Breakdown rate of labile detritus at 106:16:1";
  parameters[n].units = "d-1";
  parameters[n].sym   = "r_{Red}";
  parameters[n].value[0] = 0.04;
  parameters[n].index = n++;

  parameters[n].name  = "r_DetBL";
  parameters[n].desc  = "Breakdown rate of labile detritus at 550:30:1";
  parameters[n].units = "d-1";
  parameters[n].sym   = "r_{Atk}";
  parameters[n].value[0] = 0.001;
  parameters[n].index = n++;

  parameters[n].name  = "r_RD";
  parameters[n].desc  = "Breakdown rate of refractory detritus";
  parameters[n].units = "d-1";
  parameters[n].sym   = "r_{R}";
  parameters[n].value[0] = 0.001;
  parameters[n].index = n++;

  parameters[n].name  = "r_DOM";
  parameters[n].desc  = "Breakdown rate of dissolved organic matter";
  parameters[n].units = "d-1";
  parameters[n].sym   = "r_{O}";
  parameters[n].value[0] = 0.0001;
  parameters[n].index = n++;

  parameters[n].name  = "Plank_resp";
  parameters[n].desc  = "Respiration as a fraction of umax";
  parameters[n].units = "none";
  parameters[n].sym   = "\\phi";
  parameters[n].value[0] = 0.025;
  parameters[n].index = n++;

  /* Sediment parameters */

  parameters[n].name  = "KO_aer";
  parameters[n].desc  = "Oxygen half-saturation for aerobic respiration";
  parameters[n].units = "mg O m-3";
  parameters[n].sym   = "K_{OA}";
  parameters[n].value[0] = 256.0;
  parameters[n].index = n++;

  parameters[n].name  = "r_nit_wc";
  parameters[n].desc  = "Maximum nitrification rate in water column";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\tau_{nit,wc}";
  parameters[n].value[0] = 0.1;
  parameters[n].index = n++;

  parameters[n].name  = "r_nit_sed";
  parameters[n].desc  = "Maximum nitrification rate in water sediment";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\tau_{nit,sed}";
  parameters[n].value[0] = 20.0;
  parameters[n].index = n++;

  parameters[n].name  = "KO_nit";
  parameters[n].desc  = "Oxygen half-saturation for nitrification";
  parameters[n].units = "mg O m-3";
  parameters[n].sym   = "K_{\\mathrm{O}_2,nit}";
  parameters[n].value[0] = 500.0;
  parameters[n].index = n++;

  parameters[n].name  = "Pads_r";
  parameters[n].desc  = "Rate at which P reaches adsorbed/desorbed equilibrium";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\tau_{Pabs}";
  parameters[n].value[0] = 0.04;
  parameters[n].index = n++;

  parameters[n].name  = "Pads_Kwc";
  parameters[n].desc  = "Freundlich Isothermic Const P adsorption to TSS in water column";
  parameters[n].units = "mg P kg TSS-1";
  parameters[n].sym   = "k_{Pads,wc}";
  parameters[n].value[0] = 300.0;
  parameters[n].index = n++;

  parameters[n].name  = "Pads_Ksed";
  parameters[n].desc  = "Freundlich Isothermic Const P adsorption to TSS in sediment";
  parameters[n].units = "mg P kg TSS-1";
  parameters[n].sym   = "k_{Pads,sed}";
  parameters[n].value[0] = 74.0;
  parameters[n].index = n++;

  parameters[n].name  = "Pads_KO";
  parameters[n].desc  = "Oxygen half-saturation for P adsorption";
  parameters[n].units = "mg O m-3";
  parameters[n].sym   = "K_{\\mathrm{O}_2,abs}";
  parameters[n].value[0] = 2000.0;
  parameters[n].index = n++;

  parameters[n].name  = "Pads_exp";
  parameters[n].desc  = "Exponent for Freundlich Isotherm";
  parameters[n].units = "none";
  parameters[n].sym   = "N/A";
  parameters[n].value[0] = 1.0;
  parameters[n].index = n++;

  parameters[n].name  = "r_den";
  parameters[n].desc  = "Maximum denitrification rate";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\tau_{denit}";
  parameters[n].value[0] = 5.0;
  parameters[n].index = n++;

  parameters[n].name  = "KO_den";
  parameters[n].desc  = "Oxygen half-saturation constant for denitrification";
  parameters[n].units = "mg O m-3";
  parameters[n].sym   = "K_{\\mathrm{O}_2,denit}";
  parameters[n].value[0] = 10000.0;
  parameters[n].index = n++;

  parameters[n].name  = "r_immob_PIP";
  parameters[n].desc  = "Rate of conversion of PIP to immobilised PIP";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\tau_{Pimm}";
  parameters[n].value[0] = 0.0012;
  parameters[n].index = n++;

  parameters[n].name  = "EpiDiffCoeff";
  parameters[n].desc  = "Sediment-water diffusion coefficient";
  parameters[n].units = "m2 s-1";
  parameters[n].sym   = "D";
  parameters[n].value[0] = 3e-7;
  parameters[n].index = n++;

  parameters[n].name  = "EpiDiffDz";
  parameters[n].desc  = "Thickness of diffusive layer";
  parameters[n].units = "m";
  parameters[n].sym   = "h";
  parameters[n].value[0] = 0.0065;
  parameters[n].index = n++;
  
  /* Marcroalgae */
  parameters[n].name  = "MAumax";
  parameters[n].desc  = "Maximum growth rate of MA at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{MA}^{max}";
  parameters[n].value[0] = 1.0;
  parameters[n].index = n++;

  parameters[n].name  = "MA_mL";
  parameters[n].desc  = "Natural (linear) mortality rate, macroalgae";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\zeta_{MA}";
  parameters[n].value[0] = 0.01;
  parameters[n].index = n++;

  parameters[n].name  = "MAleafden";
  parameters[n].desc  = "Nitrogen-specific leaf area of macroalgae";
  parameters[n].units = "m2 g N-1";
  parameters[n].sym   = "\\Omega_{MA}";
  parameters[n].value[0] = 1.0;
  parameters[n].index = n++;

  parameters[n].name  = "Benth_resp";
  parameters[n].desc  = "Respiration as a fraction of umax";
  parameters[n].units = "none";
  parameters[n].sym   = "\\phi";
  parameters[n].value[0] = 0.025;
  parameters[n].index = n++;

  /* Seagrass parameters - Zostera */
  parameters[n].name  = "SGumax";
  parameters[n].desc  = "Maximum growth rate of SG at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{SG}^{max}";
  parameters[n].value[0] = 0.4;
  parameters[n].ref = "x2 nighttime, x2 for roots.";
  parameters[n].index = n++;
  
  parameters[n].name  = "SG_KN";
  parameters[n].desc  = "Half-saturation of SG N uptake in SED";
  parameters[n].units = "mg N m-3";
  parameters[n].sym   = "K_{SG,N}";
  parameters[n].value[0] = 420.0;
  parameters[n].ref = "Lee and Dunton (1999) 1204-1215. Table 3 Zostera";
  parameters[n].index = n++;

  parameters[n].name  = "SG_KP";
  parameters[n].desc  = "Half-saturation of SG P uptake in SED";
  parameters[n].units = "mg P m-3";
  parameters[n].sym   = "K_{SG,P}";
  parameters[n].value[0] = 96.0;
  parameters[n].ref = "Gras et al. (2003) Aquatic Botany 76:299-315. Thalassia testudinum.";
  parameters[n].index = n++;

  parameters[n].name  = "SG_mL";
  parameters[n].desc  = "Natural (linear) mortality rate aboveground seagrass";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\zeta_{SG_A}";
  parameters[n].value[0] = 0.04;
  parameters[n].stderr = 0.001;
  parameters[n].ref = "Fourquean et al.( 2003) Chem. Ecol. 19: 373-390.Thalassia leaves with one component decay";
  parameters[n].index = n++;

  parameters[n].name  = "SGROOT_mL";
  parameters[n].desc  = "Natural (linear) mortality rate belowground seagrass";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\zeta_{SG_B}";
  parameters[n].value[0] = 0.004;
  parameters[n].stderr = 0.0002;
  parameters[n].ref = "Fourquean et al. (2003) Chem. Ecol. 19: 373-390. Thalassia roots with one component decay";
  parameters[n].index = n++;

  parameters[n].name  = "SGfrac";
  parameters[n].desc  = "Fraction (target) of SG biomass below-ground";
  parameters[n].units = "-";
  parameters[n].sym   = "f_{below,SG}";
  parameters[n].value[0] = 0.5;
  parameters[n].ref = "Duarte (1999) Aquatic Biol. 65: 159-174, Zostera capricornii.";
  parameters[n].index = n++;

  parameters[n].name  = "SGtransrate";
  parameters[n].desc  = "Time scale for seagrass translocation";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\tau_{tran,SG}";
  parameters[n].value[0] = 0.0333;
  parameters[n].ref = "Loosely based on Zostera marine Kaldy et al., 2013 MEPS 487:27-39";
  parameters[n].index = n++;

  parameters[n].name  = "SGleafden";
  parameters[n].desc  = "Nitrogen-specific leaf area of seagrass";
  parameters[n].units = "m2 g N-1";
  parameters[n].sym   = "\\Omega_{SG}";
  parameters[n].value[0] = 1.5;
  parameters[n].ref = "Zostera capricornia: leaf dimensions Kemp et al (1987) Mar Ecol. Prog. Ser. 41:79-86.";
  parameters[n].index = n++;

  parameters[n].name  = "SGseedfrac";
  parameters[n].desc  = "Seagrass seed biomass as fraction of 63 % cover";
  parameters[n].units = "-";
  parameters[n].sym   = "f_{seed,SG}";
  parameters[n].value[0] = 0.01;
  parameters[n].ref = "No source";
  parameters[n].index = n++;

  parameters[n].name  = "SGorient";
  parameters[n].desc  = "Sine of nadir Zostera canopy bending angle";
  parameters[n].units = "-";
  parameters[n].sym   = "\\sin \\beta_{blade,SG}";
  parameters[n].value[0] = 0.5;
  parameters[n].ref = "No source";
  parameters[n].index = n++;

  parameters[n].name  = "SGmlr";
  parameters[n].desc  = "Compensation irradiance for Zostera";
  parameters[n].units = "mol m-2";
  parameters[n].sym   = "E_{comp,SG}";
  parameters[n].value[0] = 4.5;
  parameters[n].ref = "Chartrand (2012) Tech report.";
  parameters[n].index = n++;

  parameters[n].name  = "SGrootdepth";
  parameters[n].desc  = "Maximum depth for Zostera roots";
  parameters[n].units = "m";
  parameters[n].sym   = "z_{root,SG}";
  parameters[n].value[0] = -0.15;
  parameters[n].ref = "Roberts (1993) Aust. J. Mar. Fresh. Res. 44:85-100.";
  parameters[n].index = n++;


  /* Seagrass parameters - Halophila */
  parameters[n].name  = "SGHumax";
  parameters[n].desc  = "Maximum growth rate of SGH at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{SGH}^{max}";
  parameters[n].value[0] = 0.4;
  parameters[n].ref = "x2 nighttime, x2 for roots.";
  parameters[n].index = n++;

  parameters[n].name  = "SGH_KN";
  parameters[n].desc  = "Half-saturation of SGH N uptake in SED";
  parameters[n].units = "mg N m-3";
  parameters[n].sym   = "K_{SGH,N}";
  parameters[n].value[0] = 420.0;
  parameters[n].index = n++;

  parameters[n].name  = "SGH_KP";
  parameters[n].desc  = "Half-saturation of SGH P uptake in SED";
  parameters[n].units = "mg P m-3";
  parameters[n].sym   = "K_{SGH,P}";
  parameters[n].value[0] = 96.0;
  parameters[n].index = n++;

  parameters[n].name  = "SGHleafden";
  parameters[n].desc  = "Nitrogen-specific leaf area of SGH"; 
  parameters[n].units = "m2 g N-1";
  parameters[n].sym   = "\\Omega_{SGH}";
  parameters[n].value[0] = 1.9;
  parameters[n].ref = "Halophila ovalis: leaf dimensions from Vermaat et al. (1995)";
  parameters[n].index = n++;

  parameters[n].name  = "SGH_mL";
  parameters[n].desc  = "Natural (linear) mortality rate, aboveground SGH";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\zeta_{SGH_A}";
  parameters[n].value[0] = 0.08;
  parameters[n].stderr = 0.001;
  parameters[n].ref = "Fourquean et al.(2003) Chem. Ecol. 19: 373-390.Thalassia leaves with one component decay";
  parameters[n].index = n++;

  parameters[n].name  = "SGHROOT_mL";
  parameters[n].desc  = "Natural (linear) mortality rate, belowground SGH";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\zeta_{SGH_B}";
  parameters[n].value[0] = 0.004;
  parameters[n].stderr = 0.0002;
  parameters[n].ref = "Fourquean et al. (2003) Chem. Ecol. 19: 373-390. Thalassia roots with one component decay";
  parameters[n].index = n++;

  parameters[n].name  = "SGHfrac";
  parameters[n].desc  = "Fraction (target) of SGH biomass below-ground";
  parameters[n].units = "-";
  parameters[n].sym   = "f_{below,SGH}";
  parameters[n].value[0] = 0.25;
  parameters[n].ref = "Duarte (1999) Aquatic Biol. 65: 159-174, Halophila ovalis.";
  parameters[n].index = n++;

  parameters[n].name  = "SGHtransrate";
  parameters[n].desc  = "Time scale for seagrass translocation";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\tau_{tran,SGH}";
  parameters[n].value[0] = 0.0333;
  parameters[n].ref = "Loosely based on Zostera marine Kaldy et al., 2013 MEPS 487:27-39";
  parameters[n].index = n++;

  parameters[n].name  = "SGHseedfrac";
  parameters[n].desc  = "Halophila seed biomass as fraction of 63 % cover";
  parameters[n].units = "-";
  parameters[n].sym   = "f_{seed,SGH}";
  parameters[n].value[0] = 0.01;
  parameters[n].ref = "No source";
  parameters[n].index = n++;

  parameters[n].name  = "SGHorient";
  parameters[n].desc  = "Sine of nadir Halophila canopy bending angle";
  parameters[n].units = "-";
  parameters[n].sym   = "\\sin \\beta_{blade,SGH}";
  parameters[n].value[0] = 1.0;
  parameters[n].ref = "No source";
  parameters[n].index = n++;

  parameters[n].name  = "SGHmlr";
  parameters[n].desc  = "Compensation irradiance for Halophila";
  parameters[n].units = "mol m-2";
  parameters[n].sym   = "E_{comp,SGH}";
  parameters[n].value[0] = 2.8;
  parameters[n].ref = "Longstaff 2003 UQ PhD thesis";
  parameters[n].index = n++;

  parameters[n].name  = "SGHrootdepth";
  parameters[n].desc  = "Maximum depth for Halophila roots";
  parameters[n].units = "m";
  parameters[n].sym   = "z_{root,SGH}";
  parameters[n].value[0] = -0.08;
  parameters[n].ref = "Roberts (1993) Aust. J. Mar. Fresh. Res. 44:85-100.";
  parameters[n].index = n++;

  /* Seagrass parameters - Deep */

  parameters[n].name  = "SGDumax";
  parameters[n].desc  = "Maximum growth rate of SGD at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{SGD}^{max}";
  parameters[n].value[0] = 0.4;
  parameters[n].ref = "x2 nighttime, x2 for roots.";
  parameters[n].index = n++;

  parameters[n].name  = "SGD_KN";
  parameters[n].desc  = "Half-saturation of SGD N uptake in SED";
  parameters[n].units = "mg N m-3";
  parameters[n].sym   = "K_{SGD,N}";
  parameters[n].value[0] = 420.0;
  parameters[n].index = n++;

  parameters[n].name  = "SGD_KP";
  parameters[n].desc  = "Half-saturation of SGD P uptake in SED";
  parameters[n].units = "mg P m-3";
  parameters[n].sym   = "K_{SGD,P}";
  parameters[n].value[0] = 96.0;
  parameters[n].index = n++;

  parameters[n].name  = "SGDleafden";
  parameters[n].desc  = "Nitrogen-specific leaf area of SGD"; 
  parameters[n].units = "m2 g N-1";
  parameters[n].sym   = "\\Omega_{SGD}";
  parameters[n].value[0] = 1.9;
  parameters[n].ref = "Halophila ovalis: leaf dimensions from Vermaat et al. (1995)";
  parameters[n].index = n++;

  parameters[n].name  = "SGD_mL";
  parameters[n].desc  = "Natural (linear) mortality rate, aboveground SGD";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\zeta_{SGD_A}";
  parameters[n].value[0] = 0.06;
  parameters[n].stderr = 0.001;
  parameters[n].ref = "NESP project";
  parameters[n].index = n++;

  parameters[n].name  = "SGDROOT_mL";
  parameters[n].desc  = "Natural (linear) mortality rate, belowground SGD";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\zeta_{SGD_B}";
  parameters[n].value[0] = 0.004;
  parameters[n].stderr = 0.00002;
  parameters[n].ref = "NESP project";
  parameters[n].index = n++;

  parameters[n].name  = "SGDfrac";
  parameters[n].desc  = "Fraction (target) of SGD biomass below-ground";
  parameters[n].units = "-";
  parameters[n].sym   = "f_{below,SGD}";
  parameters[n].value[0] = 0.25;
  parameters[n].ref = "Duarte (1999) Aquatic Biol. 65: 159-174, Halophila ovalis.";
  parameters[n].index = n++;

  parameters[n].name  = "SGDtransrate";
  parameters[n].desc  = "Time scale for seagrass translocation";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\tau_{tran,SGD}";
  parameters[n].value[0] = 0.0333;
  parameters[n].ref = "Loosely based on Zostera marine Kaldy et al., 2013 MEPS 487:27-39";
  parameters[n].index = n++;

  parameters[n].name  = "SGDseedfrac";
  parameters[n].desc  = "Halophila seed biomass as fraction of 63 % cover";
  parameters[n].units = "-";
  parameters[n].sym   = "f_{seed,SGD}";
  parameters[n].value[0] = 0.01;
  parameters[n].ref = "No source";
  parameters[n].index = n++;

  parameters[n].name  = "SGDorient";
  parameters[n].desc  = "Sine of nadir Halophila canopy bending angle";
  parameters[n].units = "-";
  parameters[n].sym   = "\\sin \\beta_{blade,SGD}";
  parameters[n].value[0] = 1.0;
  parameters[n].ref = "No source";
  parameters[n].index = n++;

  parameters[n].name  = "SGDmlr";
  parameters[n].desc  = "Compensation irradiance for Halophila";
  parameters[n].units = "mol m-2";
  parameters[n].sym   = "E_{comp,SGD}";
  parameters[n].value[0] = 1.5;
  parameters[n].ref = "NESP project";
  parameters[n].index = n++;

  parameters[n].name  = "SGDrootdepth";
  parameters[n].desc  = "Maximum depth for Halophila roots";
  parameters[n].units = "m";
  parameters[n].sym   = "z_{root,SGD}";
  parameters[n].value[0] = -0.05;
  parameters[n].ref = "NESP project";
  parameters[n].index = n++;

  parameters[n].name  = "SGD_tau_critical";
  parameters[n].desc  = "Critical shear stress for SGD loss";
  parameters[n].units = "N m^{-2}";
  parameters[n].sym   = "\\tau_{SGD,shear}";
  parameters[n].value[0] = 1.0;
  parameters[n].ref = "NESP project";
  parameters[n].index = n++;

  parameters[n].name  = "SGD_tau_time";
  parameters[n].desc  = "Time-scale for critical shear stress for SGD loss";
  parameters[n].units = "s";
  parameters[n].sym   = "\\tau_{SGD,time}";
  parameters[n].value[0] = 43200.0;
  parameters[n].ref = "NESP project";
  parameters[n].index = n++;

  /* Corals */
  parameters[n].name  = "dissCaCO3_sed";
  parameters[n].desc  = "net dissolution rate of sediment without coral";
  parameters[n].units = "mmol C m-2 s-1";
  parameters[n].sym   = "d_{sand}";
  parameters[n].value[0] = 0.0007;
  parameters[n].index = n++;
  parameters[n].ref = "Anthony et al. (2013), Biogeosciences 10:4897-4909, Fig 5E: -1 2 3 6  mmol m-2 h-1";

  parameters[n].name  = "CHarea";
  parameters[n].desc  = "Grid scale to reef scale ratio";
  parameters[n].units = "m2 m-2";
  parameters[n].sym   = "A_{CH}";
  parameters[n].value[0] = 0.1;
  parameters[n].index = n++;
  parameters[n].ref = "Heron Island for on 4 km model.";

  parameters[n].name  = "CHpolypden";
  parameters[n].desc  = "Nitrogen-specific host area of coral polyp";
  parameters[n].units = "m2 g N-1";
  parameters[n].sym   = "\\Omega_{CH}";
  parameters[n].value[0] = 2.0;
  parameters[n].index = n++;

  parameters[n].name  = "CHumax";
  parameters[n].desc  = "Max. growth rate of Coral at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{CH}^{max}";
  parameters[n].value[0] = 0.05;
  parameters[n].index = n++;

  parameters[n].name  = "CSumax";
  parameters[n].desc  = "Max. growth rate of zooxanthellae at Tref";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\mu_{CS}^{max}";
  parameters[n].value[0] = 0.4;
  parameters[n].index = n++;

  parameters[n].name  = "CSrad";
  parameters[n].desc  = "Radius of the zooxanthellae ";
  parameters[n].units = "m";
  parameters[n].sym   = "r_{CS}";
  parameters[n].value[0] = 5e-06;
  parameters[n].index = n++;

  parameters[n].name  = "CHmort";
  parameters[n].desc  = "Quadratic mortality rate of coral polyp ";
  parameters[n].units = "(g N m-2)-1 d-1";
  parameters[n].sym   = "\\zeta_{CH}";
  parameters[n].value[0] = 0.01;
  parameters[n].index = n++;

  parameters[n].name  = "CSmort";
  parameters[n].desc  = "Linear mortality rate of zooxanthellae ";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\zeta_{CS}";
  parameters[n].value[0] = 0.04;
  parameters[n].index = n++;

  parameters[n].name  = "CHremin";
  parameters[n].desc  = "Fraction of coral host death translocated.";
  parameters[n].units = "-";
  parameters[n].sym   = "f_{remin}";
  parameters[n].value[0] = 0.5;
  parameters[n].index = n++;

  parameters[n].name  = "Splank";
  parameters[n].desc  = "Rate coefficent for particle uptake by corals";
  parameters[n].units = "m d-1";
  parameters[n].sym   = "S_{part}";
  parameters[n].value[0] = 3.0;
  parameters[n].ref = "Ribes and Atkinson (2007) Coral Reefs 26: 413-421";
  parameters[n].index = n++;

  parameters[n].name  = "k_day_coral";
  parameters[n].desc  = "Maximum daytime coral calcification";
  parameters[n].units = "mmol C m-2 s-1";
  parameters[n].sym   = "k_{day}";
  parameters[n].value[0] = 0.0132;
  parameters[n].ref = "Anthony et al. (2013), Biogeosciences 10:4897-4909, Fig 5A: 50, 50, 35 55 mmol m-2 h-1 for Acropora aspera n=4";
  parameters[n].index = n++;

  parameters[n].name  = "k_night_coral";
  parameters[n].desc  = "Maximum nightime coral calcification";
  parameters[n].units = "mmol C m-2 s-1";
  parameters[n].sym   = "k_{night}";
  parameters[n].value[0] = 0.0069;
  parameters[n].ref = "Anthony et al. (2013), Biogeosciences 10:4897-4909, Fig 5A: 20, 30, 20, 30  mmol m-2 h-1 for Acropora aspera n=4";
  parameters[n].index = n++;

  parameters[n].name  = "dissCaCO3_shelf";
  parameters[n].desc  = "Carbonate sediment dissolution rate on shelf";
  parameters[n].units = "mmol C m-2 s-1";
  parameters[n].sym   = "d_{shelf}";
  parameters[n].value[0] = 0.0001;
  parameters[n].ref = "Cyronak, T. et al., LO 58:131-143. Heron Island study.";
  parameters[n].index = n++;

  parameters[n].name  = "ageing_decay";
  parameters[n].desc  = "Age tracer growth rate per day";
  parameters[n].units = "d d-1";
  parameters[n].sym   = "n/a";
  parameters[n].value[0] = 1.0;
  parameters[n].ref = "EMS manual";
  parameters[n].index = n++;

  parameters[n].name  = "anti_ageing_decay";
  parameters[n].desc  = "Age tracer decay rate per day outside source";
  parameters[n].units = "d-1";
  parameters[n].sym   = "\\Phi";
  parameters[n].value[0] = 0.1;
  parameters[n].ref = "EMS manual";
  parameters[n].index = n++;


  /* Assign acutal number of parameters */
  *nprm = n;

  /* Assign string values */
  assign_string_values(parameters, *nprm);
}



void eco_params_std(parameter_info **params, int *nprm)
{
  int n;
  parameter_info *parameters = *params;

  init_parameter_values(parameters);

  n = 0;
  
  parameters[n].name  = "ZL_E";
  parameters[n].desc  = "Growth efficiency, large zooplankton";
  parameters[n].units = "none";
  parameters[n].value[0] = 0.38;
  parameters[n].index = n++;
  
  parameters[n].name = "ZS_E";
  parameters[n].desc = "Growth efficiency, small zooplankton";
  parameters[n].units= "none";
  parameters[n].value[0] = 0.38;
  parameters[n].index = n++;

  parameters[n].name = "SG_KN";
  parameters[n].desc = "Half-saturation of SG N uptake in SED";
  parameters[n].units= "mg N m-3";
  parameters[n].value[0] = 5.0;
  parameters[n].index = n++;

  parameters[n].name = "SG_KP";
  parameters[n].desc = "Half-saturation of SG P uptake in SED";
  parameters[n].units= "mg N m-3";
  parameters[n].value[0] = 5.0;
  parameters[n].index = n++;

  parameters[n].name = "PhyL_mL";
  parameters[n].desc = "Natural (linear) mortality rate, large phytoplankton (in sediment)";
  parameters[n].units= "d-1";
  parameters[n].value[0] = 0.14;
  parameters[n].index = n++;

  parameters[n].name = "PhyS_mL";
  parameters[n].desc = "Natural (linear) mortality rate, small phytoplankton (in sediment)";
  parameters[n].units= "d-1";
  parameters[n].value[0] = 0.14;
  parameters[n].index = n++;

  parameters[n].name = "MA_mL";
  parameters[n].desc = "Natural (linear) mortality rate, macroalgae";
  parameters[n].units= "d-1";                                 
  parameters[n].value[0] = 0.01;
  parameters[n].index = n++;

  parameters[n].name = "SG_mL";
  parameters[n].desc = "Natural (linear) mortality rate, seagrass";
  parameters[n].units= "d-1";
  parameters[n].value[0] = 0.00274;                                
  parameters[n].index = n++;

  parameters[n].name = "MPB_mQ";
  parameters[n].desc = "Natural (quadratic) mortality rate, microphytobenthos";
  parameters[n].units= "d-1 (mg N m-3)-1"; 
  parameters[n].value[0] = 0.0003;
  parameters[n].index = n++;

  parameters[n].name = "ZL_mQ";
  parameters[n].desc = "Natural (quadratic) mortality rate, large zooplankton";
  parameters[n].units= "d-1 (mg N m-3)-1";
  parameters[n].value[0] = 0.01;                   
  parameters[n].index = n++;

  parameters[n].name = "ZS_mQ";
  parameters[n].desc = "Natural (quadratic) mortality rate, small zooplankton";
  parameters[n].units= "d-1 (mg N m-3)-1";
  parameters[n].value[0] = 0.02;                    
  parameters[n].index = n++;

  parameters[n].name = "ZL_FDG";
  parameters[n].desc = "Fraction of growth inefficiency lost to detritus, large zooplankton";
  parameters[n].units= "none";
  parameters[n].value[0] = 1;
  parameters[n].index = n++;
  
  parameters[n].name = "ZL_FDM";
  parameters[n].desc = "Fraction of mortality lost to detritus, large zooplankton";
  parameters[n].units= "none";
  parameters[n].value[0] = 1;                              
  parameters[n].index = n++;
  
  parameters[n].name = "ZS_FDG";
  parameters[n].desc = "Fraction of growth inefficiency lost to detritus, small zooplankton";
  parameters[n].units= "none";
  parameters[n].value[0] = 1;                           
  parameters[n].index = n++;
  
  parameters[n].name = "ZS_FDM";
  parameters[n].desc = "Fraction of mortality lost to detritus, small zooplankton";
  parameters[n].units= "none";
  parameters[n].value[0] = 1;                               
  parameters[n].index = n++;
  
  parameters[n].name = "F_LD_RD";
  parameters[n].desc = "Fraction of labile detritus converted to refractory detritus";
  parameters[n].units= "none";
  parameters[n].value[0] = 0.19;                                 
  parameters[n].index = n++;
  
  parameters[n].name = "F_LD_DOM";
  parameters[n].desc = "Fraction of labile detritus converted to dissolved organic matter";
  parameters[n].units= "none";
  parameters[n].value[0] = 0.01;                                           
  parameters[n].index = n++;
  
  parameters[n].name = "NtoCHL";
  parameters[n].desc = "Nitrogen:Chlorophyll A ratio in phytoplankton by weight";
  parameters[n].units= "m-1";
  parameters[n].value[0] = 7;                           
  parameters[n].index = n++;
  
  parameters[n].name = "k_w";
  parameters[n].desc = "Background light attenuation coefficient";
  parameters[n].units= "none";
  parameters[n].value[0] = 0.03;                                                                      
  parameters[n].index = n++;
  
  parameters[n].name = "k_DOR_N";
  parameters[n].desc = "DOR_N-specific light attenuation coefficient";
  parameters[n].units= "m-1 (mg N m-3)-1";
  parameters[n].value[0] = 0.0009;                                      
  parameters[n].index = n++;
  
  parameters[n].name = "k_DetL";
  parameters[n].desc = "Detrital N-specific light attenuation coefficient";
  parameters[n].units= "m-1 (mg N m-3)-1";
  parameters[n].value[0] = 0.0038;           
  parameters[n].index = n++;
  
  parameters[n].name = "Light_lambda";
  parameters[n].desc = "Wavelengths of light";
  parameters[n].units= "nm";
  /* Free array allocated above */
  d_free_1d(parameters[n].value);
  parameters[n].num_values = 23;
  parameters[n].value = d_alloc_1d(23);
  parameters[n].value[0] = 300.0;    
  parameters[n].value[1] = 320.0;
  parameters[n].value[2] = 340.0;
  parameters[n].value[3] = 360.0;
  parameters[n].value[4] = 380.0;
  parameters[n].value[5] = 400.0;
  parameters[n].value[6] = 420.0;
  parameters[n].value[7] = 440.0;
  parameters[n].value[8] = 460.0;
  parameters[n].value[9] = 480.0;
  parameters[n].value[10] = 500.0;
  parameters[n].value[11] = 520.0;
  parameters[n].value[12] = 540.0;
  parameters[n].value[13] = 560.0;
  parameters[n].value[14] = 580.0;
  parameters[n].value[15] = 600.0;
  parameters[n].value[16] = 620.0;
  parameters[n].value[17] = 640.0;
  parameters[n].value[18] = 660.0;
  parameters[n].value[19] = 680.0;
  parameters[n].value[20] = 700.0;
  parameters[n].value[21] = 720.0;
  parameters[n].value[22] = 800.0;
  parameters[n].index = n++;
  
  parameters[n].name = "k_TSS";
  parameters[n].desc = "TSS-specific light attenuation coefficient";
  parameters[n].units= "m-1 (kg m-3)-1";
  parameters[n].value[0] = 30.0;                                     
  parameters[n].index = n++;
  
  parameters[n].name = "k_C_fw";
  parameters[n].desc = "CDOM attentuation coefficient of freshwater";
  parameters[n].units= "m-1";
  parameters[n].value[0] = 4.4;  
  parameters[n].index = n++;
  
  parameters[n].name = "k_SWR_PAR";
  parameters[n].desc = "fraction of incident solar radiation that is PAR";
  parameters[n].units= "none";
  parameters[n].value[0] = 0.43;          
  parameters[n].index = n++;
  
  parameters[n].name = "Q10";
  parameters[n].desc = "Temperature coefficient for rate parameters";
  parameters[n].units= "none";
  parameters[n].value[0] = 2.0;  
  parameters[n].index = n++;
  
  parameters[n].name = "PLumax";
  parameters[n].desc = "Maximum growth rate of PL at Tref";
  parameters[n].units= "d-1";
  parameters[n].value[0] = 1.75;                             
  parameters[n].index = n++;

  parameters[n].name = "PLrad";
  parameters[n].desc = "Radius of the large phytoplankton cells";
  parameters[n].units= "m";
  parameters[n].value[0] = 10e-06;
  parameters[n].index = n++;

  parameters[n].name = "PLabsorb";
  parameters[n].desc = "Absorption coefficient of a PL cell";
  parameters[n].units= "m-1";
  parameters[n].value[0] = 50000.;             
  parameters[n].index = n++;

  parameters[n].name = "PLSh";
  parameters[n].desc = "Sherwood number for the PS dimensionless";
  parameters[n].units= "none";
  parameters[n].value[0] = 1;                                
  parameters[n].index = n++;

  parameters[n].name = "PLtable";
  parameters[n].desc = "Netcdf lookup table";
  parameters[n].units= "none";
  parameters[n].stringvalue = "10plkINP";                
  parameters[n].index = n++;

  parameters[n].name = "PLn";
  parameters[n].desc = "Number of limiting nutrients";
  parameters[n].units= "none";
  parameters[n].value[0] = 3;                                                  
  parameters[n].index = n++;
  
  parameters[n].name = "PSumax";
  parameters[n].desc = "Maximum growth rate of PS at Tref";
  parameters[n].units= "d-1";
  parameters[n].value[0] = 1.25;       
  parameters[n].index = n++;
  
  parameters[n].name = "PSrad";
  parameters[n].desc = "Radius of the small phytoplankton cells";
  parameters[n].units= "m";
  parameters[n].value[0] = 2.5e-06;
  parameters[n].index = n++;
  
  parameters[n].name = "PSabsorb";
  parameters[n].desc = "Absorption coefficient of a PS cell";
  parameters[n].units= "m-1";
  parameters[n].value[0] = 50000;
  parameters[n].index = n++;
  
  parameters[n].name = "PSSh";
  parameters[n].desc = "Sherwood number for the PL dimensionless";
  parameters[n].units= "none";
  parameters[n].value[0] = 1;                             
  parameters[n].index = n++;
  
  parameters[n].name = "PStable";
  parameters[n].desc = "Netcdf lookup table";
  parameters[n].units= "none";
  parameters[n].stringvalue = "10plkINP";        
  parameters[n].index = n++;
  
  parameters[n].name = "PSn";
  parameters[n].desc = "Number of limiting nutrients";
  parameters[n].units= "none";
  parameters[n].value[0] = 3;                                            
  parameters[n].index = n++;

  parameters[n].name = "Tricho_aA";
  parameters[n].desc = "Nitrogen specific absorption cross-section of Tricho";
  parameters[n].units= "m2 mg N-1";
  parameters[n].value[0] = 1e-03;          
  parameters[n].index = n++;
  
  parameters[n].name = "Tricho_mL";
  parameters[n].desc = "Linear mortality for Tricho in sediment";
  parameters[n].units= "d-1";
  parameters[n].value[0] = 0.14;         
  parameters[n].index = n++;

  parameters[n].name = "Tricho_mQ";
  parameters[n].desc = "Quadratic mortality for Tricho due to phages in wc"; // NEEDS ATTENTION but a rough first guess";
  parameters[n].units= "d-1 (mg N m-3)-1";
  parameters[n].value[0] = 0.015;    
  parameters[n].index = n++;

  parameters[n].name = "Tricho_crit";
  parameters[n].desc = "Critical Tricho above which quadratic mortality applies";
  parameters[n].units= "mg N m-3";
  parameters[n].value[0] = 0.0;       
  parameters[n].index = n++;
  
  parameters[n].name = "Tricho_rad";
  parameters[n].desc = "Radius of Trichodesmium colonies";
  parameters[n].units= "m";
  parameters[n].value[0] = 5.0e-06;
  parameters[n].index = n++;

  parameters[n].name = "Tricho_umax";
  parameters[n].desc = "Maximum growth rate of Trichodesmium at Tref";
  parameters[n].units= "d-1";
  parameters[n].value[0] = 1.0;                            
  parameters[n].index = n++;

  parameters[n].name = "p_min";
  parameters[n].desc = "Minimum density of Trichodesmium";
  parameters[n].units= "kg m-3";
  parameters[n].value[0] = 990.0;                              
  parameters[n].index = n++;

  parameters[n].name = "p_max";
  parameters[n].desc = "Maximum density of Trichodesmium";
  parameters[n].units= "kg m-3";
  parameters[n].value[0] = 1060.0;                            
  parameters[n].index = n++;

  parameters[n].name = "Tricho_Sh";
  parameters[n].desc = "Sherwood number for Trichodesmium";
  parameters[n].units= "none";
  parameters[n].value[0] = 1;                                    
  parameters[n].index = n++;
  
  parameters[n].name = "N2";
  parameters[n].desc = "Concentration of dissolved N2";
  parameters[n].units= "mg N m-3";
  parameters[n].value[0] = 2e-4;                                    
  parameters[n].index = n++;

  parameters[n].name = "MBumax";
  parameters[n].desc = "Maximum growth rate of MB at Tref";
  parameters[n].units= "d-1";
  parameters[n].value[0] = 0.35;                                                
  parameters[n].index = n++;
  
  parameters[n].name = "MBrad";
  parameters[n].desc = "Radius of the large phytoplankton cells";
  parameters[n].units= "m";
  parameters[n].value[0] = 1e-05;
  parameters[n].index = n++;
  
  parameters[n].name = "MBabsorb";
  parameters[n].desc = "Absorption coefficient of a MB cell";
  parameters[n].units= "m-1";
  parameters[n].value[0] = 50000;                                           
  parameters[n].index = n++;
  
  parameters[n].name = "MBSh";
  parameters[n].desc = "Sherwood number for the PL dimensionless";
  parameters[n].units= "none";
  parameters[n].value[0] = 1;                           
  parameters[n].index = n++;
  
  parameters[n].name = "MBtable";
  parameters[n].desc = "Netcdf lookup table";
  parameters[n].units= "none";
  parameters[n].stringvalue= "10plkINP";       
  parameters[n].index = n++;
  
  parameters[n].name = "MBn";
  parameters[n].desc = "Number of limiting nutrients";
  parameters[n].units= "none";
  parameters[n].value[0] = 3;                                     
  parameters[n].index = n++;
  
  parameters[n].name = "MAumax";
  parameters[n].desc = "Maximum growth rate of MA at Tref";
  parameters[n].units= "d-1";
  parameters[n].value[0] = 0.02;                                            
  parameters[n].index = n++;
  
  parameters[n].name = "MAaA";
  parameters[n].desc = "Nitrogen specific absorption cross-section of MA";
  parameters[n].units= "m2 mg N-1";
  parameters[n].value[0] = 1e-03;      
  parameters[n].index = n++;
  
  parameters[n].name = "MAtable";
  parameters[n].desc = "Netcdf lookup table";
  parameters[n].units= "none";
  parameters[n].stringvalue = "10benINP";       
  parameters[n].index = n++;
  
  parameters[n].name = "MAn";
  parameters[n].desc = "Number of limiting nutrients";
  parameters[n].units= "none";
  parameters[n].value[0] = 3;                                      
  parameters[n].index = n++;
  
  parameters[n].name = "MAm";
  parameters[n].desc = "Stoichometry coefficient of Phosphorus";
  parameters[n].units= "none";
  parameters[n].value[0] = 2.4e-06;
  parameters[n].index = n++;
  
  parameters[n].name = "SGumax";
  parameters[n].desc = "Maximum growth rate of SG at Tref";
  parameters[n].units= "d-1";
  parameters[n].value[0] = 0.1;                                            
  parameters[n].index = n++;
  
  parameters[n].name = "SGaA";
  parameters[n].desc = "Nitrogen specific absorption cross-section of SG";
  parameters[n].units= "m2 mg N-1";
  parameters[n].value[0] = 1e-05;         
  parameters[n].index = n++;
  
  parameters[n].name = "SGm";
  parameters[n].desc = "Stoichometry coefficient of Phosphorus";
  parameters[n].units= "none";
  parameters[n].value[0] = 2.4e-06;
  parameters[n].index = n++;
  
  parameters[n].name = "ZSumax";
  parameters[n].desc = "Maximum growth rate of ZS at Tref";
  parameters[n].units= "d-1";
  parameters[n].value[0] = 3;                                               
  parameters[n].index = n++;
  
  parameters[n].name = "ZSrad";
  parameters[n].desc = "Radius of the small zooplankton cells";
  parameters[n].units= "m-1";
  parameters[n].value[0] = 12.5e-06;
  parameters[n].index = n++;
  
  parameters[n].name = "ZSswim";
  parameters[n].desc = "Swimming velocity for small zooplankton";
  parameters[n].units= "m s-1";
  parameters[n].value[0] = 2.0e-4;
  parameters[n].index = n++;
  
  parameters[n].name = "ZSmeth";
  parameters[n].desc = "Grazing technique of small zooplankton";
  parameters[n].units= "none";
  parameters[n].stringvalue = "rect";                          
  parameters[n].index = n++;

  parameters[n].name = "ZLumax";
  parameters[n].desc = "Maximum growth rate of ZL at Tref";
  parameters[n].units= "d-1";
  parameters[n].value[0] = 1.0;                                                
  parameters[n].index = n++;

  parameters[n].name = "ZLrad";
  parameters[n].desc = "Radius of the large zooplankton cells";
  parameters[n].units= "m-1";
  parameters[n].value[0] = 5.0e-04;                                             
  parameters[n].index = n++;

  parameters[n].name = "ZLswim";
  parameters[n].desc = "Swimming velocity for large zooplankton";
  parameters[n].units= "m s-1";
  parameters[n].value[0] = 1.5e-3;                                                      
  parameters[n].index = n++;

  parameters[n].name = "ZLmeth";
  parameters[n].desc = "Grazing technique of small zooplankton";
  parameters[n].units= "none";
  parameters[n].stringvalue = "rect";                            
  parameters[n].index = n++;

  parameters[n].name = "TKEeps";
  parameters[n].desc = "TKE dissipation in water column";
  parameters[n].units= "m2s-3";
  parameters[n].value[0] = 1.0e-6;                                      
  parameters[n].index = n++;

  parameters[n].name = "cf";
  parameters[n].desc = "drag coefficient of the benthic surface";
  parameters[n].units= "none";
  parameters[n].value[0] = 0.005;
  parameters[n].index = n++;
  
  parameters[n].name = "Ub";
  parameters[n].desc = "velocity at the top of the ben. bound. layer";
  parameters[n].units= "m s-1";
  parameters[n].value[0] = 0.1;    
  parameters[n].index = n++;

  parameters[n].name = "ks";
  parameters[n].desc = "sand-grain roughness of the benthos";
  parameters[n].units= "m";
  parameters[n].value[0] = 0.1;
  parameters[n].index = n++;

  parameters[n].name = "F_RD_DOM";
  parameters[n].desc = "fraction of refractory detritus that breaks down to DOM";
  parameters[n].units= "none";
  parameters[n].value[0] = 0.05;                        
  parameters[n].index = n++;

  parameters[n].name = "r_floc";
  parameters[n].desc = "rate at which TSS floculates above 10 PSU";
  parameters[n].units= "d-1";
  parameters[n].value[0] = 0.01;                              
  parameters[n].index = n++;

  parameters[n].name = "r_DetPL";
  parameters[n].desc = "Breakdown rate of labile detritus at 106:16:1";
  parameters[n].units= "d-1";
  parameters[n].value[0] = 0.1;    
  parameters[n].index = n++;

  parameters[n].name = "r_DetBL";
  parameters[n].desc = "Breakdown rate of labile detritus at 550:30:1";
  parameters[n].units= "d-1";
  parameters[n].value[0] = 0.1;    
  parameters[n].index = n++;

  parameters[n].name = "r_RD";
  parameters[n].desc = "Breakdown rate of refractory detritus";
  parameters[n].units= "d-1";
  parameters[n].value[0] = 0.0036;                                              
  parameters[n].index = n++;

  parameters[n].name = "r_DOM";
  parameters[n].desc = "Breakdown rate of dissolved organic matter";
  parameters[n].units= "d-1";
  parameters[n].value[0] = 0.00176;                            
  parameters[n].index = n++;

  parameters[n].name = "Tref";
  parameters[n].desc = "Reference temperature";
  parameters[n].units= "Deg C";
  parameters[n].value[0] = 15.0;                      
  parameters[n].index = n++;

  parameters[n].name = "Plank_resp";
  parameters[n].desc = "Respiration as a fraction of umax";
  parameters[n].units= "none";
  parameters[n].value[0] = 0.025;                                          
  parameters[n].index = n++;

  parameters[n].name = "Benth_resp";
  parameters[n].desc = "Respiration as a fraction of umax";
  parameters[n].units= "none";
  parameters[n].value[0] = 0.025;                                          
  parameters[n].index = n++;

  parameters[n].name = "DFumax";
  parameters[n].desc = "Maximum growth rate of dinoflagellate at Tref";
  parameters[n].units= "d-1";
  parameters[n].value[0] = 0.4;      
  parameters[n].index = n++;

  parameters[n].name = "DFrad";
  parameters[n].desc = "Radius of dinoflagellate cells";
  parameters[n].units= "m";
  parameters[n].value[0] = 10.0e-6;                                    
  parameters[n].index = n++;

  parameters[n].name = "DFabsorb";
  parameters[n].desc = "Absorption coefficient of a dinoflagellate cell";
  parameters[n].units= "m-1";
  parameters[n].value[0] = 40000.0;      
  parameters[n].index = n++;

  parameters[n].name = "DFSh";
  parameters[n].desc = "Sherwood number for dinoflagellate";
  parameters[n].units= "none";
  parameters[n].value[0] = 1.0;                                            
  parameters[n].index = n++;

  parameters[n].name = "DFn";
  parameters[n].desc = "Number of limiting nutrients for Dinoflagellate";
  parameters[n].units= "none";
  parameters[n].value[0] = 3;          
  parameters[n].index = n++;

  parameters[n].name = "DFtable";
  parameters[n].desc = "Netcdf lookup table for Dinoflagellate";
  parameters[n].units= "none";
  parameters[n].stringvalue = "10plkINP";                                   
  parameters[n].index = n++;

  parameters[n].name = "NOumax";
  parameters[n].desc = "Maximum growth rate of Nodularia at Tref";
  parameters[n].units= "d-1";
  parameters[n].value[0] = 0.6;
  parameters[n].index = n++;

  parameters[n].name = "NOrad";
  parameters[n].desc = "Radius of Nodularia cells";
  parameters[n].units= "m";
  parameters[n].value[0] = 20.0e-6;                         
  parameters[n].index = n++;

  parameters[n].name = "NOabsorb";
  parameters[n].desc = "Absorption coefficient of a Nodularia cell";
  parameters[n].units= "m-1";
  parameters[n].value[0] = 20000.0;                            
  parameters[n].index = n++;

  parameters[n].name = "NO_TG";
  parameters[n].desc = "Temperature below which growth of Nodularia cells slows down";
  parameters[n].units= "Deg";
  parameters[n].value[0] = 19.0;                                   
  parameters[n].index = n++;

  parameters[n].name = "NO_SG";
  parameters[n].desc = "Salinity above which growth of Nodularia cells slows down";
  parameters[n].units= "PSU";
  parameters[n].value[0] = 25.0;                          
  parameters[n].index = n++;

  parameters[n].name = "NOSh";
  parameters[n].desc = "Sherwood number for Nodularia";
  parameters[n].units= "none";
  parameters[n].value[0] = 1.0;                                  
  parameters[n].index = n++;

  parameters[n].name = "NOn";
  parameters[n].desc = "Number of limiting nutrients for Nodularia";
  parameters[n].units= "none";
  parameters[n].value[0] = 3;   
  parameters[n].index = n++;

  parameters[n].name = "NOtable";
  parameters[n].desc = "NetCDF lookup table for Nodularia";
  parameters[n].units= "none";
  parameters[n].stringvalue = "10plkINP";                       
  parameters[n].index = n++;

  parameters[n].name = "NO_mL";
  parameters[n].desc = "Linear mortality rate, nodularia";
  parameters[n].units= "d-1";
  parameters[n].value[0] = 0.1;                                           
  parameters[n].index = n++;

  parameters[n].name = "NO_mQ";
  parameters[n].desc = "Quadratic mortality rate for nodularia";
  parameters[n].units= "d-1";
  parameters[n].value[0] = 0.0002;
  parameters[n].index = n++;

  parameters[n].name = "DFCtoNvar";
  parameters[n].desc = "Maximal to minimal C:N ratio in Dinoflagellate";
  parameters[n].units= "none";
  parameters[n].value[0] = 1.5;    
  parameters[n].index = n++;

  parameters[n].name = "KO_aer";
  parameters[n].desc = "Oxygen half-saturation for aerobic respiration";
  parameters[n].units= "mg O m-3";
  parameters[n].value[0] = 500.0;  
  parameters[n].index = n++;

  parameters[n].name = "r_nit_wc";
  parameters[n].desc = "Maximal nitrification rate in water column";
  parameters[n].units= "d-1";
  parameters[n].value[0] = 0.1;
  parameters[n].index = n++;

  parameters[n].name = "r_nit_sed";
  parameters[n].desc = "Maximal nitrification rate in water sediment";
  parameters[n].units= "d-1";
  parameters[n].value[0] = 5.0;  
  parameters[n].index = n++;

  parameters[n].name = "KO_nit";
  parameters[n].desc = "Oxygen half-saturation for nitrification";
  parameters[n].units= "mg O m-3";
  parameters[n].value[0] = 500.0;
  parameters[n].index = n++;

  parameters[n].name = "Pads_r";
  parameters[n].desc = "Rate at which P reaches adsorbed/desorbed equilibrium";
  parameters[n].units= "d-1";
  parameters[n].value[0] = 0.04;               
  parameters[n].index = n++;

  parameters[n].name = "Pads_Kwc";
  parameters[n].desc = "Freundlich Isothermic Const P adsorption to TSS in water column";
  parameters[n].units= "mg P kg TSS-1";
  parameters[n].value[0] = 300.0;                                  
  parameters[n].index = n++;

  parameters[n].name = "Pads_Ksed";
  parameters[n].desc = "Freundlich Isothermic Const P adsorption to TSS in sediment";
  parameters[n].units= "mg P kg TSS-1";
  parameters[n].value[0] = 74.0;                           
  parameters[n].index = n++;

  parameters[n].name = "Pads_KO";
  parameters[n].desc = "Oxygen half-saturation for P adsorption";
  parameters[n].units= "mg O m-3";
  parameters[n].value[0] = 2000.0;
  parameters[n].index = n++;

  parameters[n].name = "Pads_exp";
  parameters[n].desc = "Exponent for Freundlich Isotherm";
  parameters[n].units= "none";
  parameters[n].value[0] = 1.0;                                      
  parameters[n].index = n++;

  parameters[n].name = "PD_mL";
  parameters[n].desc = "Linear mortality for dinoflagellate in sediment";
  parameters[n].units= "d-1";
  parameters[n].value[0] = 0.14;       
  parameters[n].index = n++;

  parameters[n].name = "r_den";
  parameters[n].desc = "Maximum denitrification rate";
  parameters[n].units= "d-1";
  parameters[n].value[0] = 40.0;                                
  parameters[n].index = n++;

  parameters[n].name = "KO_den";
  parameters[n].desc = "Oxygen content at 50% denitrification rate";
  parameters[n].units= "mg O m-3";
  parameters[n].value[0] = 10000.0;                            
  parameters[n].index = n++;

  parameters[n].name = "r_floc_sed";
  parameters[n].desc = "Rate of the TSS floculation in sediment";
  parameters[n].units= "d-1";
  parameters[n].value[0] = 0.001;
  parameters[n].index = n++;

  parameters[n].name = "r_bury_TSS";
  parameters[n].desc = "Rate of the TSS burying";
  parameters[n].units= "d-1";
  parameters[n].value[0] = 0.001;                    
  parameters[n].index = n++;

  parameters[n].name = "r_immob_PIP";
  parameters[n].desc = "Rate of conversion of PIP to immobilised PIP";
  parameters[n].units= "d-1";
  parameters[n].value[0] = 0.0012;                               
  parameters[n].index = n++;

  parameters[n].name = "IDF";
  parameters[n].desc = "Saturation light intensity for  dinoflagellates";
  parameters[n].units= "mol photon m-2 s-1";
  parameters[n].value[0] = 1.0e-4;     
  parameters[n].index = n++;

  parameters[n].name = "Fmax_Nit_sed";
  parameters[n].desc = "Maximum nitrification efficiency";
  parameters[n].units= "none";
  parameters[n].value[0] = 1.0;                                           
  parameters[n].index = n++;

  parameters[n].name = "EpiDiffCoeff";
  parameters[n].desc = "Diffusion Coefficient";
  parameters[n].units= "m2s-1";
  parameters[n].value[0] = 3e-9;                     
  parameters[n].index = n++;

  parameters[n].name = "EpiDiffDz";
  parameters[n].desc = "Thickness of diffusive layer";
  parameters[n].units= "m";
  parameters[n].value[0] = 0.0065;                            
  parameters[n].index = n++;

  /* Lyngbya */
  parameters[n].name = " Phy_L_N2_pmax";
  parameters[n].desc = "Maximum growth rate of Lyngbya ";
  parameters[n].units= "d-1";
  parameters[n].value[0] = 0.6;                     
  /* observed values 0.7 to 2.5 */                                
  parameters[n].index = n++;

  parameters[n].name = "Plank_L_N2_mort";
  parameters[n].desc = "Mortality rate for lyngbya";
  parameters[n].units= "d-1";
  parameters[n].value[0] = 0.0002;                          
  parameters[n].index = n++;

  parameters[n].name = "Phy_L_N2_Kdin";
  parameters[n].desc = "DIN half saturation constant for lyngbya";
  parameters[n].units= "mg/L";
  /* other literature values are between 0.04 and 0.125 */
  parameters[n].value[0] = 0.03;
  parameters[n].index = n++;

  parameters[n].name = "Phy_L_N2_KP";
  parameters[n].desc = "DIP half saturation constant for lyngbya";
  parameters[n].units= "mg/L";
  /* other literature values are between 0.0025 and 0.02 */
  parameters[n].value[0] = 0.006;
  parameters[n].index = n++;

  parameters[n].name = "Phy_L_N2_alpha";
  parameters[n].desc = "initial slope of the light function (for L)"; // ie the fractional cost of fixing N and is a real value.";
  parameters[n].units= "none";
  parameters[n].value[0] = 5.7;
  parameters[n].index = n++;

  parameters[n].name = "Plank_L_N2_resp";
  parameters[n].desc = "respira L_N2";
  parameters[n].units= "none";
  parameters[n].value[0] = 0.01;
  parameters[n].index = n++;

  parameters[n].name = "NtoCHL_L_N2";
  parameters[n].desc = "Nitrogen:Chlorophyll A ratio in L_N2  by weight";
  parameters[n].units= "m-1";
  parameters[n].value[0] = 7;
  parameters[n].index = n++;

  parameters[n].name  = "SGleafden";
  parameters[n].desc  = "Nitrogen-specific leaf area of seagrass";
  parameters[n].units = "m2 g N-1";
  parameters[n].value[0] = 4.2;
  parameters[n].ref = "Zostera capricornia: leaf dimensions Kemp et al (1987) Mar Ecol. Prog. Ser. 41:79-86.";
  parameters[n].index = n++;

  parameters[n].name  = "SGfrac";
  parameters[n].desc  = "Fraction (target) of SG biomass below-ground";
  parameters[n].units = "-";
  parameters[n].value[0] = 0.4790;
  parameters[n].ref = "Duarte (1999) Aquatic Biol. 65: 159-174, Zostera capricornii.";
  parameters[n].index = n++;

  parameters[n].name  = "SGHleafden";
  parameters[n].desc  = "Nitrogen-specific leaf area of seagrass"; 
  parameters[n].units = "m2 g N-1";
  parameters[n].value[0] = 1.9;
  parameters[n].ref = "Halophila ovalis: leaf dimensions from Vermaat et al. (1995)";
  parameters[n].index = n++;

  parameters[n].name  = "SGHfrac";
  parameters[n].desc  = "Fraction (target) of SGH biomass below-ground";
  parameters[n].units = "-";
  parameters[n].value[0] = 0.278;
  parameters[n].ref = "Duarte (1999) Aquatic Biol. 65: 159-174, Halophila ovalis.";
  parameters[n].index = n++;

  parameters[n].name  = "MAleafden";
  parameters[n].desc  = "Nitrogen-specific leaf area of macroalgae";
  parameters[n].units = "m2 g N-1";
  parameters[n].value[0] = 2.0;
  parameters[n].index = n++;

  parameters[n].name  = "CHpolypden";
  parameters[n].desc  = "Nitrogen-specific host area of coral polyp";
  parameters[n].units = "m2 g N-1";
  parameters[n].value[0] = 2.0;
  parameters[n].index = n++;

  parameters[n].name  = "CSrad";
  parameters[n].desc  = "Radius of the zooxanthellae ";
  parameters[n].units = "m";
  parameters[n].value[0] = 5e-06;
  parameters[n].index = n++;

  /* Assign acutal number of parameters */
  *nprm = n;

  /* Assign string values */
  assign_string_values(parameters, *nprm);
}


/** Fills in estuarine parameter defaults
 * @param nprm Number of parameters
 * @param parameters parameter_info array
 */
void eco_params_est(parameter_info **params, int *nprm)
{
  /*
   * Cut and paste the entire contents of eco_params_std in here and
   * make the changes as needed.
   *
   * Alternatively, you may just call the above eco_params_std
   * function and then simply 'override' some of the parameters
   */
  e_quit("'biofname estuary' is not implemented as yet.  See the 'standard' counterpart function 'eco_params_std' in model/lib/ecology/parameter_defaults.c for an example and populate eco_params_est accordingly.\n");
}

/* Local helper functions */
static void init_parameter_values(parameter_info *params)
{
  int n;

  for (n = 0; n < MAX_PARAMS; n++) {
    parameter_info *prm = &params[n];
    prm->num_values = 1;
    prm->value = d_alloc_1d(prm->num_values);
    prm->value[0] = NaN;
    prm->ref = "Not attributed";
    prm->sym = "";
  }
}

static void assign_string_values(parameter_info *params, int nprm)
{
  int n;
  for (n = 0; n < nprm; n++) {
    parameter_info* prm = &params[n];
    char buf[MAXSTRLEN];
    char *str = buf;
    int i;
    
    /* Build up the array */
    if (!(isnan(prm->value[0]))) {
      for (i=0; i<prm->num_values; i++) {
	char buf2[MAXSTRLEN];
	sprintf(buf2, "%e ", prm->value[i]);
	sprintf(str, "%s", buf2);
	str += strlen(buf2);
      }
      prm->stringvalue = strdup(buf);
    } else {
      /* We can't leave the literal as it will die when freed */
      prm->stringvalue = strdup(prm->stringvalue);
      prm->value[0] = 0.0;
    }
    
    /*
     * Allocate strings
     */
    // name
    str = strdup(prm->name);
    prm->name = str;
    // desc
    str = strdup(prm->desc);
    prm->desc = str;
    // sym
    str = strdup(prm->sym);
    prm->sym = str;
    // units
    str = strdup(prm->units);
    prm->units = str;
    // ref
    str = strdup(prm->ref);
    prm->ref = str;
  }
}

// EOF
