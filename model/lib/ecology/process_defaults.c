/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: lib/ecology/process_defaults.c
 *  
 *  Description:
 *  Process defaults
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: process_defaults.c 6909 2021-09-23 00:09:58Z riz008 $
 *
 */

#include <stdlib.h>
#include "ecology_internal.h"
#include "ecology.h"
#include "eprocess.h"
#include "einterface.h"
#include "utils.h"

/*
 * Here are the list of all the default processes
 *
 * Alternatively, this could be consolidated into a single nested structure
 *
 */

// "porewater_age" water column
const char *WC_PROCESS_PORE[] = {
  "recom_extras"
};
const int NUM_WC_PROCESS_PORE =
  (int)(sizeof(WC_PROCESS_PORE)/sizeof(char *));

// "optics" epibenthos
const char *EPI_PROCESS_PORE[] = {
  "values_common_epi"
};
const int NUM_EPI_PROCESS_PORE =
  (int)(sizeof(EPI_PROCESS_PORE)/sizeof(char *));

// "optics" sediment
const char *SED_PROCESS_PORE[] = {
  "age_sed"
};
const int NUM_SED_PROCESS_PORE =
  (int)(sizeof(SED_PROCESS_PORE)/sizeof(char *));

// **********************************************************

// "optics" water column
const char *WC_PROCESS_OPT[] = {
  "values_common",
  "light_spectral_wc",
  "recom_extras"
};
const int NUM_WC_PROCESS_OPT =
  (int)(sizeof(WC_PROCESS_OPT)/sizeof(char *));

// "optics" epibenthos
const char *EPI_PROCESS_OPT[] = {
  "values_common_epi",
  "light_spectral_uq_epi"
};
const int NUM_EPI_PROCESS_OPT =
  (int)(sizeof(EPI_PROCESS_OPT)/sizeof(char *));

// "optics" sediment
const char *SED_PROCESS_OPT[] = {
  "values_common",
  "light_spectral_sed",
  "recom_extras"
};
const int NUM_SED_PROCESS_OPT =
  (int)(sizeof(SED_PROCESS_OPT)/sizeof(char *));

// **********************************************************

// "gas" water column
const char *WC_PROCESS_GAS[] = { 
  "tfactor",                 
  "carbon_chemistry_wc",
  "gas_exchange_wc(carbon,oxygen)",
  "recom_extras"
};
const int NUM_WC_PROCESS_GAS =
  (int)(sizeof(WC_PROCESS_GAS)/sizeof(char *));

// "gas" epibenthos
const char *EPI_PROCESS_GAS[] = {
  "tfactor_epi"    // removed "diffusion_epi" */
};
const int NUM_EPI_PROCESS_GAS =
  (int)(sizeof(EPI_PROCESS_GAS)/sizeof(char *));

// "gas" sediment
const char *SED_PROCESS_GAS[] = {
  "tfactor"
};
const int NUM_SED_PROCESS_GAS =
  (int)(sizeof(SED_PROCESS_GAS)/sizeof(char *));

// **********************************************************

// "standard" water column
const char *WC_PROCESS_STD[] = {
  "tfactor",                                               
  "viscosity",                                             
  "moldiff",                                               
  "remineralization",                                      
  "phytoplankton_grow_wc(small)",                          
  "phytoplankton_grow_wc(large)",                          
  "zooplankton_large_grow_wc",                             
  "zooplankton_small_grow_wc",                             
  "zooplankton_mortality_wc(small)",                       
  "zooplankton_mortality_wc(large)",                       
  "nitrification_wc",                                      
  "values_common",
  "light_wc"
};
const int NUM_WC_PROCESS_STD =
  (int)(sizeof(WC_PROCESS_STD)/sizeof(char *));

// "standard" epibenthos
const char *EPI_PROCESS_STD[] = {
  "macroalgae_grow_epi",
  "macroalgae_mortality_epi",
  "seagrass_grow_epi",
  "seagrass_mortality_epi",
  "values_common_epi",
  "light_epi",
  "diffusion_epi"
};
const int NUM_EPI_PROCESS_STD =
  (int)(sizeof(EPI_PROCESS_STD)/sizeof(char *));

// "standard" sediment
const char *SED_PROCESS_STD[] = {
  "tfactor",
  "viscosity",
  "moldiff",
  "remineralization",
  "phytoplankton_mortality_sed(small)",
  "phytoplankton_mortality_sed(large)",
  "nitrification_denitrification_sed",
  "p_adsorption_sed",
  "values_common",
  "light_sed"
};
const int NUM_SED_PROCESS_STD =
  (int)(sizeof(SED_PROCESS_STD)/sizeof(char *));


/*
 * Estuarine
 */
// "estuary" water column
const char *WC_PROCESS_EST[] = {
  "tfactor",                                               
  "viscosity",                                             
  "moldiff",                                               
  "remineralization",                                      
  "phytoplankton_grow_wc(small)",                          
  "phytoplankton_grow_wc(large)",                          
  "zooplankton_large_grow_wc",                             
  "zooplankton_small_grow_wc",                             
  "zooplankton_mortality_wc(small)",                       
  "zooplankton_mortality_wc(large)",                       
  "nitrification_wc",                                      
  "values_common",
  "light_wc"
};
const int NUM_WC_PROCESS_EST =
  (int)(sizeof(WC_PROCESS_EST)/sizeof(char *));

// "estuary" epibenthos
const char *EPI_PROCESS_EST[] = {
  "macroalgae_grow_epi",
  "macroalgae_mortality_epi",
  "seagrass_grow_epi",
  "seagrass_mortality_epi",
  "values_common_epi",
  "light_epi",
  "diffusion_epi"
};
const int NUM_EPI_PROCESS_EST =
  (int)(sizeof(EPI_PROCESS_EST)/sizeof(char *));

// "estuary" sediment
const char *SED_PROCESS_EST[] = {
  "tfactor",
  "viscosity",
  "moldiff",
  "remineralization",
  "phytoplankton_mortality_sed(small)",
  "phytoplankton_mortality_sed(large)",
  "nitrification_denitrification_sed",
  "p_adsorption_sed",
  "values_common",
  "light_sed"
};
const int NUM_SED_PROCESS_EST =
  (int)(sizeof(SED_PROCESS_EST)/sizeof(char *));

/*
 * GBR4 for nesting RECOM in GBR4\.
 */

// "GBR4" water column
const char *WC_PROCESS_GBR4[] = {
  "tfactor",
  "viscosity",
  "moldiff",
  "values_common",
  "remineralization",
  "microphytobenthos_spectral_grow_wc",
  "phytoplankton_spectral_grow_wc(small)",
  "phytoplankton_spectral_grow_wc(large)",
  "trichodesmium_mortality_wc",
  "trichodesmium_grow_wc",
  "phytoplankton_spectral_mortality_wc(small)",
  "phytoplankton_spectral_mortality_wc(large)",
  "zooplankton_mortality_wc(small)",
  "zooplankton_mortality_wc(large)",
  "zooplankton_large_carnivore_spectral_grow_wc",
  "zooplankton_small_spectral_grow_wc",
  "nitrification_wc",
  "p_adsorption_wc",
  "carbon_chemistry_wc",
  "gas_exchange_wc(oxygen,carbon)",
  "massbalance_wc",
  "light_spectral_wc",
  "age_wc",
  "recom_extras"
};
const int NUM_WC_PROCESS_GBR4 =
  (int)(sizeof(WC_PROCESS_GBR4)/sizeof(char *));

// "GBR4" epibenthos
const char *EPI_PROCESS_GBR4[] = {
  "tfactor_epi",
  "values_common_epi",
  "macroalgae_spectral_grow_epi",
  "seagrass_spectral_grow_uq_epi(Zostera)",
  "seagrass_spectral_grow_uq_epi(Halophila)",
  "coral_spectral_grow_epi",
  "macroalgae_mortality_epi",
  "seagrass_spectral_mortality_epi(Zostera)",
  "seagrass_spectral_mortality_epi(Halophila)",
  "massbalance_epi",
  "light_spectral_uq_epi",
  "diffusion_epi"
};
const int NUM_EPI_PROCESS_GBR4 =
  (int)(sizeof(EPI_PROCESS_GBR4)/sizeof(char *));

// "standard" sediment
const char *SED_PROCESS_GBR4[] = {
  "tfactor",
  "viscosity",
  "moldiff",
  "values_common",
  "remineralization",
  "light_sed",
  "microphytobenthos_spectral_grow_sed",
  "microphytobenthos_spectral_mortality_sed",
  "phytoplankton_spectral_mortality_sed(small)",
  "phytoplankton_spectral_mortality_sed(large)",
  "zooplankton_mortality_sed(small)",
  "zooplankton_mortality_sed(large)",
  "trichodesmium_mortality_sed",
  "nitrification_denitrification_sed",
  "p_adsorption_sed",
  "massbalance_sed",
  "recom_extras"
};
const int NUM_SED_PROCESS_GBR4 =
  (int)(sizeof(SED_PROCESS_GBR4)/sizeof(char *));

/*
 * GBR4_CORAL for nesting RECOM in GBR4\.
 */

// "GBR4" water column
const char *WC_PROCESS_GBR4_CORAL[] = {
  "tfactor",
  "viscosity",
  "moldiff",
  "values_common",
  "remineralization",
  "microphytobenthos_spectral_grow_wc",
  "phytoplankton_spectral_grow_wc(small)",
  "phytoplankton_spectral_grow_wc(large)",
  "trichodesmium_mortality_wc",
  "trichodesmium_grow_wc",
  "phytoplankton_spectral_mortality_wc(small)",
  "phytoplankton_spectral_mortality_wc(large)",
  "zooplankton_mortality_wc(small)",
  "zooplankton_mortality_wc(large)",
  "zooplankton_large_carnivore_spectral_grow_wc",
  "zooplankton_small_spectral_grow_wc",
  "nitrification_wc",
  "p_adsorption_wc",
  "carbon_chemistry_wc",
  "gas_exchange_wc(oxygen,carbon)",
  "massbalance_wc",
  "recom_extras",
  "light_spectral_wc(GBR,Guassian)",
  "age_wc"
};
const int NUM_WC_PROCESS_GBR4_CORAL =
  (int)(sizeof(WC_PROCESS_GBR4_CORAL)/sizeof(char *));

// "GBR4" epibenthos
const char *EPI_PROCESS_GBR4_CORAL[] = {
  "tfactor_epi",
  "values_common_epi",
  "macroalgae_spectral_grow_epi",
  "seagrass_spectral_grow_proto_epi(Zostera)",
  "seagrass_spectral_grow_proto_epi(Halophila)",
  "seagrass_spectral_grow_proto_epi(Deep)",
  "coral_spectral_grow_bleach_epi",
  "coral_spectral_carb_epi",
  "macroalgae_mortality_epi",
  "seagrass_spectral_mortality_proto_epi(Zostera)",
  "seagrass_spectral_mortality_proto_epi(Halophila)",
  "seagrass_spectral_mortality_proto_epi(Deep)",
  "massbalance_epi",
  "light_spectral_uq_epi",
  "diffusion_epi"
};
const int NUM_EPI_PROCESS_GBR4_CORAL =
  (int)(sizeof(EPI_PROCESS_GBR4_CORAL)/sizeof(char *));

// "standard" sediment
const char *SED_PROCESS_GBR4_CORAL[] = {
  "tfactor",
  "viscosity",
  "moldiff",
  "values_common",
  "remineralization",
  "light_spectral_sed",
  "microphytobenthos_spectral_grow_sed",
  "carbon_chemistry_wc",
  "microphytobenthos_spectral_mortality_sed",
  "phytoplankton_spectral_mortality_sed(small)",
  "phytoplankton_spectral_mortality_sed(large)",
  "zooplankton_mortality_sed(small)",
  "zooplankton_mortality_sed(large)",
  "trichodesmium_mortality_sed",
  "nitrification_denitrification_sed",
  "p_adsorption_sed",
  "massbalance_sed",
  "recom_extras"
};
const int NUM_SED_PROCESS_GBR4_CORAL =
  (int)(sizeof(SED_PROCESS_GBR4_CORAL)/sizeof(char *));


/*
 * GBR4_SEAGRASS for nesting RECOM in GBR4\.
 */

// "GBR4_SEAGRASS" water column
const char *WC_PROCESS_GBR4_SEAGRASS[] = {
  "tfactor",
  "viscosity",
  "moldiff",
  "values_common",
  "remineralization",
  "microphytobenthos_spectral_grow_wc",
  "phytoplankton_spectral_grow_wc(small)",
  "phytoplankton_spectral_grow_wc(large)",
  "trichodesmium_mortality_wc",
  "trichodesmium_grow_wc",
  "phytoplankton_spectral_mortality_wc(small)",
  "phytoplankton_spectral_mortality_wc(large)",
  "zooplankton_mortality_wc(small)",
  "zooplankton_mortality_wc(large)",
  "zooplankton_large_carnivore_spectral_grow_wc",
  "zooplankton_small_spectral_grow_wc",
  "nitrification_wc",
  "p_adsorption_wc",
  "carbon_chemistry_wc",
  "gas_exchange_wc(oxygen,carbon)",
  "massbalance_wc",
  "light_spectral_wc",
  "age_wc",
  "recom_extras"
};
const int NUM_WC_PROCESS_GBR4_SEAGRASS =
  (int)(sizeof(WC_PROCESS_GBR4_SEAGRASS)/sizeof(char *));

// "GBR4_SEAGRASS" epibenthos
const char *EPI_PROCESS_GBR4_SEAGRASS[] = {
  "tfactor_epi",
  "values_common_epi",
  "macroalgae_spectral_grow_epi",
  "seagrass_spectral_grow_proto_epi(Zostera)",
  "seagrass_spectral_grow_proto_epi(Halophila)",
  "seagrass_spectral_grow_proto_epi(Deep)",
  "coral_spectral_grow_epi",
  "macroalgae_mortality_epi",
  "seagrass_spectral_mortality_proto_epi(Zostera)",
  "seagrass_spectral_mortality_proto_epi(Halophila)",
  "seagrass_spectral_mortality_proto_epi(Deep)",
  "massbalance_epi",
  "light_spectral_uq_epi",
  "diffusion_epi"
};
const int NUM_EPI_PROCESS_GBR4_SEAGRASS =
  (int)(sizeof(EPI_PROCESS_GBR4_SEAGRASS)/sizeof(char *));

// "standard" sediment
const char *SED_PROCESS_GBR4_SEAGRASS[] = {
  "tfactor",
  "viscosity",
  "moldiff",
  "values_common",
  "remineralization",
  "light_sed",
  "microphytobenthos_spectral_grow_sed",
  "carbon_chemistry_wc",
  "microphytobenthos_spectral_mortality_sed",
  "phytoplankton_spectral_mortality_sed(small)",
  "phytoplankton_spectral_mortality_sed(large)",
  "zooplankton_mortality_sed(small)",
  "zooplankton_mortality_sed(large)",
  "trichodesmium_mortality_sed",
  "nitrification_denitrification_sed",
  "p_adsorption_sed",
  "massbalance_sed",
  "recom_extras"
};
const int NUM_SED_PROCESS_GBR4_SEAGRASS =
  (int)(sizeof(SED_PROCESS_GBR4_SEAGRASS)/sizeof(char *));


/*
 * BGC2p0 for nesting RECOM\.
 */

// "BGC2p0" water column
const char *WC_PROCESS_BGC2p0[] = {
  "tfactor",
  "viscosity",
  "moldiff",
  "values_common",
  "remineralization",
  "microphytobenthos_spectral_grow_wc",
  "phytoplankton_spectral_grow_wc(small)",
  "phytoplankton_spectral_grow_wc(large)",
  "trichodesmium_mortality_wc",
  "trichodesmium_grow_wc",
  "phytoplankton_spectral_mortality_wc(small)",
  "phytoplankton_spectral_mortality_wc(large)",
  "zooplankton_mortality_wc(small)",
  "zooplankton_mortality_wc(large)",
  "zooplankton_large_carnivore_spectral_grow_wc",
  "zooplankton_small_spectral_grow_wc",
  "nitrification_wc",
  "p_adsorption_wc",
  "carbon_chemistry_wc",
  "gas_exchange_wc(oxygen,carbon)",
  "massbalance_wc",
  "light_spectral_wc",
  "age_wc",
  "recom_extras"
};
const int NUM_WC_PROCESS_BGC2p0 =
  (int)(sizeof(WC_PROCESS_BGC2p0)/sizeof(char *));

// "BGC2p0" epibenthos
const char *EPI_PROCESS_BGC2p0[] = {
  "tfactor_epi",
  "values_common_epi",
  "macroalgae_spectral_grow_epi",
  "seagrass_spectral_grow_proto_epi(Zostera)",
  "seagrass_spectral_grow_proto_epi(Halophila)",
  "seagrass_spectral_grow_proto_epi(Deep)",
  "coral_spectral_grow_epi",
  "macroalgae_mortality_epi",
  "seagrass_spectral_mortality_proto_epi(Zostera)",
  "seagrass_spectral_mortality_proto_epi(Halophila)",
  "seagrass_spectral_mortality_proto_epi(Deep)",
  "massbalance_epi",
  "light_spectral_uq_epi",
  "diffusion_epi"
};
const int NUM_EPI_PROCESS_BGC2p0 =
  (int)(sizeof(EPI_PROCESS_BGC2p0)/sizeof(char *));

// "standard" sediment
const char *SED_PROCESS_BGC2p0[] = {
  "tfactor",
  "viscosity",
  "moldiff",
  "values_common",
  "remineralization",
  "light_spectral_sed",
  "microphytobenthos_spectral_grow_sed",
  "carbon_chemistry_wc",
  "microphytobenthos_spectral_mortality_sed",
  "phytoplankton_spectral_mortality_sed(small)",
  "phytoplankton_spectral_mortality_sed(large)",
  "zooplankton_mortality_sed(small)",
  "zooplankton_mortality_sed(large)",
  "trichodesmium_mortality_sed",
  "nitrification_denitrification_sed",
  "p_adsorption_sed",
  "massbalance_sed",
  "recom_extras"
};
const int NUM_SED_PROCESS_BGC2p0 =
  (int)(sizeof(SED_PROCESS_BGC2p0)/sizeof(char *));


/*
 * BGC3p1 for nesting RECOM\.
 */

// "BGC3p1" water column
const char *WC_PROCESS_BGC3p1[] = {
  "tfactor",
  "viscosity",
  "moldiff",
  "values_common",
  "remineralization",
  "microphytobenthos_spectral_grow_wc",
  "phytoplankton_spectral_grow_wc(small)",
  "phytoplankton_spectral_grow_wc(large)",
  "trichodesmium_mortality_wc",
  "trichodesmium_spectral_grow_wc",
  "phytoplankton_spectral_mortality_wc(small)",
  "phytoplankton_spectral_mortality_wc(large)",
  "zooplankton_mortality_wc(small)",
  "zooplankton_mortality_wc(large)",
  "zooplankton_large_carnivore_spectral_grow_wc",
  "zooplankton_small_spectral_grow_wc",
  "nitrification_wc",
  "p_adsorption_wc",
  "carbon_chemistry_wc",
  "gas_exchange_wc(carbon,oxygen)",
  "massbalance_wc",
  "light_spectral_wc(H,HPLC)",
  "recom_extras"
};
const int NUM_WC_PROCESS_BGC3p1 =
  (int)(sizeof(WC_PROCESS_BGC3p1)/sizeof(char *));

// "BGC3p1" epibenthos
const char *EPI_PROCESS_BGC3p1[] = {
  "tfactor_epi",
  "values_common_epi",
  "macroalgae_spectral_grow_epi",
  "seagrass_spectral_grow_epi(Zostera)",
  "seagrass_spectral_grow_epi(Halophila)",
  "seagrass_spectral_grow_epi(Deep)",
  "coral_spectral_grow_bleach_epi",
  "coral_spectral_carb_epi(H)",
  "macroalgae_mortality_epi",
  "seagrass_spectral_mortality_proto_epi(Zostera)",
  "seagrass_spectral_mortality_proto_epi(Halophila)",
  "seagrass_spectral_mortality_proto_epi(Deep)",
  "massbalance_epi",
  "light_spectral_uq_epi",
  "diffusion_epi"
};
const int NUM_EPI_PROCESS_BGC3p1 =
  (int)(sizeof(EPI_PROCESS_BGC3p1)/sizeof(char *));

// "standard" sediment
const char *SED_PROCESS_BGC3p1[] = {
  "tfactor",
  "viscosity",
  "moldiff",
  "values_common",
  "remineralization",
  "light_spectral_sed(HPLC)",
  "microphytobenthos_spectral_grow_sed",
  "carbon_chemistry_wc",
  "microphytobenthos_spectral_mortality_sed",
  "phytoplankton_spectral_mortality_sed(small)",
  "phytoplankton_spectral_mortality_sed(large)",
  "zooplankton_mortality_sed(small)",
  "zooplankton_mortality_sed(large)",
  "trichodesmium_mortality_sed",
  "nitrification_denitrification_sed",
  "p_adsorption_sed",
  "age_sed",
  "massbalance_sed",
  "recom_extras"
};
const int NUM_SED_PROCESS_BGC3p1 =
  (int)(sizeof(SED_PROCESS_BGC3p1)/sizeof(char *));


/*
 * Dead_only for nesting RECOM\.
 */

// "Dead_only" water column
const char *WC_PROCESS_DEAD[] = {
  "tfactor",
  "viscosity",
  "moldiff",
  "values_common",
  "remineralization",
  "nitrification_denitrification_anammox",
  "p_adsorption",
  "carbon_chemistry_wc",
  "gas_exchange_wc(carbon,oxygen)",
  "massbalance_wc",
  "recom_extras"
};
const int NUM_WC_PROCESS_DEAD =
  (int)(sizeof(WC_PROCESS_DEAD)/sizeof(char *));

// "DEAD" epibenthos
const char *EPI_PROCESS_DEAD[] = {
  "tfactor_epi",
  "values_common_epi",
  "massbalance_epi",
  "diffusion_epi"
};
const int NUM_EPI_PROCESS_DEAD =
  (int)(sizeof(EPI_PROCESS_DEAD)/sizeof(char *));

// "DEAD" sediment
const char *SED_PROCESS_DEAD[] = {
  "tfactor",
  "viscosity",
  "moldiff",
  "values_common",
  "remineralization",
  "carbon_chemistry_wc",
  "nitrification_denitrification_anammox",
  "p_adsorption",
  "age_sed",
  "massbalance_sed",
  "recom_extras"
};
const int NUM_SED_PROCESS_DEAD =
  (int)(sizeof(SED_PROCESS_DEAD)/sizeof(char *));

/*
 * TASSE1p0 for nesting RECOM\.
 */

// "tasse" water column
const char *WC_PROCESS_TASSE1p0[] = {
  "tfactor",
  "viscosity",
  "moldiff",
  "values_common",
  "remineralization",
  "microphytobenthos_spectral_grow_wc",
  "phytoplankton_spectral_grow_wc(small)",
  "phytoplankton_spectral_grow_wc(large)",
  "dinoflagellate_spectral_grow_wc",
  "phytoplankton_spectral_mortality_wc(small)",
  "phytoplankton_spectral_mortality_wc(large)",
  "dinoflagellate_spectral_mortality_wc",
  "zooplankton_mortality_wc(small)",
  "zooplankton_mortality_wc(large)",
  "zooplankton_large_carnivore_spectral_grow_wc",
  "zooplankton_small_spectral_grow_wc",
  "nitrification_denitrification_anammox",
  "p_adsorption",
  "carbon_chemistry_wc",
  "gas_exchange_wc(carbon,oxygen)",
  "massbalance_wc",
  "light_spectral_wc(T,HPLC)",
  "recom_extras"
};
const int NUM_WC_PROCESS_TASSE1p0 =
  (int)(sizeof(WC_PROCESS_TASSE1p0)/sizeof(char *));

// "TASSE1p0" epibenthos
const char *EPI_PROCESS_TASSE1p0[] = {
  "tfactor_epi",
  "values_common_epi",
  "macroalgae_spectral_grow_epi()",
  "seagrass_spectral_grow_epi(Zostera)",
  "macroalgae_mortality_epi()",
  "seagrass_spectral_mortality_proto_epi(Zostera)",
  "massbalance_epi()",
  "light_spectral_uq_epi()",
  "diffusion_epi",
  "filter_feeder_epi"
};
const int NUM_EPI_PROCESS_TASSE1p0 =
  (int)(sizeof(EPI_PROCESS_TASSE1p0)/sizeof(char *));

// "TASSE1p0" sediment
const char *SED_PROCESS_TASSE1p0[] = {
  "tfactor",
  "viscosity",
  "moldiff",
  "values_common",
  "remineralization",
  "light_spectral_sed(HPLC)",
  "microphytobenthos_spectral_grow_sed",
  "carbon_chemistry_wc()",
  "microphytobenthos_spectral_mortality_sed",
  "phytoplankton_spectral_mortality_sed(small)",
  "phytoplankton_spectral_mortality_sed(large)",
  "zooplankton_mortality_sed(small)",
  "zooplankton_mortality_sed(large)",
  "dinoflagellate_spectral_mortality_sed",
  "nitrification_denitrification_anammox",
  "p_adsorption",
  "age_sed",
  "massbalance_sed",
  "recom_extras"
};
const int NUM_SED_PROCESS_TASSE1p0 =
  (int)(sizeof(SED_PROCESS_TASSE1p0)/sizeof(char *));



/*
 * GBR4_SMALL for nesting RECOM in GBR4\.
 */

// "GBR4_SMALL" water column
const char *WC_PROCESS_GBR4_SMALL[] = {
  "tfactor",
  "viscosity",
  "moldiff",
  "remineralization",
  "microphytobenthos_spectral_grow_wc",
  "phytoplankton_spectral_grow_wc(small)",
  "phytoplankton_spectral_grow_wc(large)",
  "trichodesmium_mortality_wc",
  "trichodesmium_grow_wc",
  "phytoplankton_spectral_mortality_wc(small)",
  "phytoplankton_spectral_mortality_wc(large)",
  "zooplankton_mortality_wc(small)",
  "zooplankton_mortality_wc(large)",
  "zooplankton_large_spectral_grow_wc",
  "zooplankton_small_spectral_grow_wc",
  "nitrification_wc",
  "p_adsorption_wc",
  // "carbon_chemistry_wc",
  //  "gas_exchange_wc(oxygen,carbon)",
  "values_common",
  //  "massbalance_wc",
  "light_spectral_wc",
};
const int NUM_WC_PROCESS_GBR4_SMALL =
  (int)(sizeof(WC_PROCESS_GBR4_SMALL)/sizeof(char *));

// "GBR4_SMALL" epibenthos
const char *EPI_PROCESS_GBR4_SMALL[] = {
  // "macroalgae_spectral_grow_epi",
  // "seagrass_spectral_grow_epi(Zostera)",
  // "coral_spectral_grow_epi",
  // "macroalgae_mortality_epi",
  // "seagrass_spectral_mortality_epi(Zostera)",
  "values_common_epi",
  // "massbalance_epi",
  "light_spectral_epi",
  "diffusion_epi"
};
const int NUM_EPI_PROCESS_GBR4_SMALL =
  (int)(sizeof(EPI_PROCESS_GBR4_SMALL)/sizeof(char *));

// "GBR4_SMALL" sediment
const char *SED_PROCESS_GBR4_SMALL[] = {
  "tfactor",
  "viscosity",
  "moldiff",
  "remineralization",
  "light_sed",
  "microphytobenthos_spectral_grow_sed",
  "microphytobenthos_spectral_mortality_sed",
  "phytoplankton_spectral_mortality_sed(small)",
  "phytoplankton_spectral_mortality_sed(large)",
  "trichodesmium_mortality_sed",
  "nitrification_denitrification_sed",
  "p_adsorption_sed",
  "values_common",
  // "massbalance_sed"
};
const int NUM_SED_PROCESS_GBR4_SMALL =
  (int)(sizeof(SED_PROCESS_GBR4_SMALL)/sizeof(char *));


/*
 * "standard" processes
 */
static void eco_processes_std(int type, const char **procs[], int *nprocs)
{
  /* Key off process type */
  switch(type) {
  case PT_WC:
    *procs  = WC_PROCESS_STD;
    *nprocs = NUM_WC_PROCESS_STD;
    return;
  case PT_EPI:
    *procs  = EPI_PROCESS_STD;
    *nprocs = NUM_EPI_PROCESS_STD;
    return;
  case PT_SED:
    *procs  = SED_PROCESS_STD;
    *nprocs = NUM_SED_PROCESS_STD;
    return;
  case PT_COL:
    *procs = NULL;
    *nprocs = 0;
    return;
  default:
    e_quit("eco_processes_std: Unknown Process Type '%d'\n", type);
  }    
}

/*
 * "porewater_age" processes
 */
static void eco_processes_pore(int type, const char **procs[], int *nprocs)
{
  /* Key off process type */
  switch(type) {
  case PT_WC:
    *procs  = WC_PROCESS_PORE;
    *nprocs = NUM_WC_PROCESS_PORE;
    return;
  case PT_EPI:
    *procs  = EPI_PROCESS_PORE;
    *nprocs = NUM_EPI_PROCESS_PORE;
    return;
  case PT_SED:
    *procs  = SED_PROCESS_PORE;
    *nprocs = NUM_SED_PROCESS_PORE;
    return;
  case PT_COL:
    *procs = NULL;
    *nprocs = 0;
    return;
  default:
    e_quit("eco_processes_opt: Unknown Process Type '%d'\n", type);
  }
}

/*
 * "optics_only" processes
 */
static void eco_processes_opt(int type, const char **procs[], int *nprocs)
{
  /* Key off process type */
  switch(type) {
  case PT_WC:
    *procs  = WC_PROCESS_OPT;
    *nprocs = NUM_WC_PROCESS_OPT;
    return;
  case PT_EPI:
    *procs  = EPI_PROCESS_OPT;
    *nprocs = NUM_EPI_PROCESS_OPT;
    return;
  case PT_SED:
    *procs  = SED_PROCESS_OPT;
    *nprocs = NUM_SED_PROCESS_OPT;
    return;
  case PT_COL:
    *procs = NULL;
    *nprocs = 0;
    return;
  default:
    e_quit("eco_processes_opt: Unknown Process Type '%d'\n", type);
  }
}

/*
 * "gas_only" processes
 */
static void eco_processes_gas(int type, const char **procs[], int *nprocs)
{
  /* Key off process type */
  switch(type) {
  case PT_WC:
    *procs  = WC_PROCESS_GAS;
    *nprocs = NUM_WC_PROCESS_GAS;
    return;
  case PT_EPI:
    *procs  = EPI_PROCESS_GAS;
    *nprocs = NUM_EPI_PROCESS_GAS;
    return;
  case PT_SED:
    *procs  = SED_PROCESS_GAS;
    *nprocs = NUM_SED_PROCESS_GAS;
    return;
  case PT_COL:
    *procs = NULL;
    *nprocs = 0;
    return;
  default:
    e_quit("eco_processes_gas: Unknown Process Type '%d'\n", type);
  }
}
/*
 * "estuarine" processes
 */
static void eco_processes_est(int type, const char **procs[], int *nprocs)
{
  /* Key off process type */
  switch(type) {
  case PT_WC:
    *procs  = WC_PROCESS_EST;
    *nprocs = NUM_WC_PROCESS_EST;
    return;
  case PT_EPI:
    *procs  = EPI_PROCESS_EST;
    *nprocs = NUM_EPI_PROCESS_EST;
    return;
  case PT_SED:
    *procs  = SED_PROCESS_EST;
    *nprocs = NUM_SED_PROCESS_EST;
    return;
  case PT_COL:
    *procs = NULL;
    *nprocs = 0;
    return;
  default:
    e_quit("eco_processes_est: Unknown Process Type '%d'\n", type);
  }
}

/*
 * "GBR4" processes
 */
static void eco_processes_gbr4(int type, const char **procs[], int *nprocs)
{
  /* Key off process type */
  switch(type) {
  case PT_WC:
    *procs  = WC_PROCESS_GBR4;
    *nprocs = NUM_WC_PROCESS_GBR4;
    return;
  case PT_EPI:
    *procs  = EPI_PROCESS_GBR4;
    *nprocs = NUM_EPI_PROCESS_GBR4;
    return;
  case PT_SED:
    *procs  = SED_PROCESS_GBR4;
    *nprocs = NUM_SED_PROCESS_GBR4;
    return;
  case PT_COL:
    *procs = NULL;
    *nprocs = 0;
    return;
  default:
    e_quit("eco_processes_gbr4: Unknown Process Type '%d'\n", type);
  }    
}

/*
 * "GBR4_CORAL" processes
 */
static void eco_processes_gbr4_coral(int type, const char **procs[], int *nprocs)
{
  /* Key off process type */
  switch(type) {
  case PT_WC:
    *procs  = WC_PROCESS_GBR4_CORAL;
    *nprocs = NUM_WC_PROCESS_GBR4_CORAL;
    return;
  case PT_EPI:
    *procs  = EPI_PROCESS_GBR4_CORAL;
    *nprocs = NUM_EPI_PROCESS_GBR4_CORAL;
    return;
  case PT_SED:
    *procs  = SED_PROCESS_GBR4_CORAL;
    *nprocs = NUM_SED_PROCESS_GBR4_CORAL;
    return;
  case PT_COL:
    *procs = NULL;
    *nprocs = 0;
    return;
  default:
    e_quit("eco_processes_gbr4_coral: Unknown Process Type '%d'\n", type);
  }    
}

/*
 * "GBR4_SEAGRASS" processes
 */
static void eco_processes_gbr4_seagrass(int type, const char **procs[], int *nprocs)
{
  /* Key off process type */
  switch(type) {
  case PT_WC:
    *procs  = WC_PROCESS_GBR4_SEAGRASS;
    *nprocs = NUM_WC_PROCESS_GBR4_SEAGRASS;
    return;
  case PT_EPI:
    *procs  = EPI_PROCESS_GBR4_SEAGRASS;
    *nprocs = NUM_EPI_PROCESS_GBR4_SEAGRASS;
    return;
  case PT_SED:
    *procs  = SED_PROCESS_GBR4_SEAGRASS;
    *nprocs = NUM_SED_PROCESS_GBR4_SEAGRASS;
    return;
  case PT_COL:
    *procs = NULL;
    *nprocs = 0;
    return;
  default:
    e_quit("eco_processes_gbr4_seagrass: Unknown Process Type '%d'\n", type);
  }    
}

/*
 * "BGC2p0" processes
 */
static void eco_processes_bgc2p0(int type, const char **procs[], int *nprocs)
{
  /* Key off process type */
  switch(type) {
  case PT_WC:
    *procs  = WC_PROCESS_BGC2p0;
    *nprocs = NUM_WC_PROCESS_BGC2p0;
    return;
  case PT_EPI:
    *procs  = EPI_PROCESS_BGC2p0;
    *nprocs = NUM_EPI_PROCESS_BGC2p0;
    return;
  case PT_SED:
    *procs  = SED_PROCESS_BGC2p0;
    *nprocs = NUM_SED_PROCESS_BGC2p0;
    return;
  case PT_COL:
    *procs = NULL;
    *nprocs = 0;
    return;
  default:
    e_quit("eco_processes_bgc2p0: Unknown Process Type '%d'\n", type);
  }    
}

/*
 * "BGC3p1" processes
 */
static void eco_processes_bgc3p1(int type, const char **procs[], int *nprocs)
{

  /* Key off process type */
  switch(type) {
  case PT_WC:
    *procs  = WC_PROCESS_BGC3p1;
    *nprocs = NUM_WC_PROCESS_BGC3p1;
    return;
  case PT_EPI:
    *procs  = EPI_PROCESS_BGC3p1;
    *nprocs = NUM_EPI_PROCESS_BGC3p1;
    return;
  case PT_SED:
    *procs  = SED_PROCESS_BGC3p1;
    *nprocs = NUM_SED_PROCESS_BGC3p1;
    return;
  case PT_COL:
    *procs = NULL;
    *nprocs = 0;
    return;
  default:
    e_quit("eco_processes_bgc3p1: Unknown Process Type '%d'\n", type);
  }
}

/*
 * "Dead_only" processes
 */
static void eco_processes_dead(int type, const char **procs[], int *nprocs)
{

  /* Key off process type */
  switch(type) {
  case PT_WC:
    *procs  = WC_PROCESS_DEAD;
    *nprocs = NUM_WC_PROCESS_DEAD;
    return;
  case PT_EPI:
    *procs  = EPI_PROCESS_DEAD;
    *nprocs = NUM_EPI_PROCESS_DEAD;
    return;
  case PT_SED:
    *procs  = SED_PROCESS_DEAD;
    *nprocs = NUM_SED_PROCESS_DEAD;
    return;
  case PT_COL:
    *procs = NULL;
    *nprocs = 0;
    return;
  default:
    e_quit("eco_processes_dead_only: Unknown Process Type '%d'\n", type);
  }
}

/*
 * "TASSE1p0" processes
 */
static void eco_processes_TASSE1p0(int type, const char **procs[], int *nprocs)
{

  /* Key off process type */
  switch(type) {
  case PT_WC:
    *procs  = WC_PROCESS_TASSE1p0;
    *nprocs = NUM_WC_PROCESS_TASSE1p0;
    return;
  case PT_EPI:
    *procs  = EPI_PROCESS_TASSE1p0;
    *nprocs = NUM_EPI_PROCESS_TASSE1p0;
    return;
  case PT_SED:
    *procs  = SED_PROCESS_TASSE1p0;
    *nprocs = NUM_SED_PROCESS_TASSE1p0;
    return;
  case PT_COL:
    *procs = NULL;
    *nprocs = 0;
    return;
  default:
    e_quit("eco_processes_TASSE1p0: Unknown Process Type '%d'\n", type);
  }
}

/*
 * "GBR4_SMALL" processes
 */
static void eco_processes_gbr4_small(int type, const char **procs[], int *nprocs)
{
  /* Key off process type */
  switch(type) {
  case PT_WC:
    *procs  = WC_PROCESS_GBR4_SMALL;
    *nprocs = NUM_WC_PROCESS_GBR4_SMALL;
    return;
  case PT_EPI:
    *procs  = EPI_PROCESS_GBR4_SMALL;
    *nprocs = NUM_EPI_PROCESS_GBR4_SMALL;
    return;
  case PT_SED:
    *procs  = SED_PROCESS_GBR4_SMALL;
    *nprocs = NUM_SED_PROCESS_GBR4_SMALL;
    return;
  case PT_COL:
    *procs = NULL;
    *nprocs = 0;
    return;
  default:
    e_quit("eco_processes_gbr4_small: Unknown Process Type '%d'\n", type);
  }    
}

/*-------------------------------------------------------------------*/
/* Dynamically prescribed ecology parameter specification            */
/*-------------------------------------------------------------------*/
// static void eco_processes_auto(EPROCESSTYPE type, const char **procs[], int *nprocs)
static void eco_processes_auto(int type, const char **procs[], int *nprocs)

{
  int i;
  int numwater = 0;
  int numepi = 0;
  int numsed = 0;
  char **water, **sediment, **epibenthos;
  ecology E;
  ecology *e = &E;

  /* 
   * Note: This code is not functionaly correct and will seg fault on execution!
   * (farhan)
   */

  if (numwater) {
    type = PT_WC;
    e->npr[type] = numwater;
    e->processes[type] = calloc(numwater, sizeof(void*));
    water = (char **)malloc(numwater * sizeof(char *));
    for (i = 0; i < numwater; i++) {
      water[i] = (char *)malloc(sizeof(char)*MAXSTRLEN);
      strcpy(water[i], "");
      e->processes[type][i] = eprocess_create(e, type, water[i]);
    }
  }
  if (numsed) {
    type = PT_SED;
    e->npr[type] = numsed;
    e->processes[type] = calloc(numsed, sizeof(void*));
    sediment = (char **)malloc(numsed * sizeof(char *));
    for (i = 0; i < numsed; i++) {
      sediment[i] = (char *)malloc(sizeof(char)*MAXSTRLEN);
      strcpy(sediment[i], "");
      e->processes[type][i] = eprocess_create(e, type, sediment[i]);
    }
  }
  if (numepi) {
    type = PT_EPI;
    e->npr[type] = numepi;
    e->processes[type] = calloc(numepi, sizeof(void*));
    epibenthos = (char **)malloc(numepi * sizeof(char *));
    for (i = 0; i < numepi; i++) {
      epibenthos[i] = (char *)malloc(sizeof(char)*MAXSTRLEN);
      strcpy(epibenthos[i], "");
      e->processes[type][i] = eprocess_create(e, type, epibenthos[i]);
    }
  }

  if (numwater)
    free((char **)water);
  if (numsed)
    free((char **)sediment);
  if (numepi)
    free((char **)epibenthos);
}

/*
 * This function is called from create_processes and serves as a
 * gateway for the different default processes types.
 */
int get_eco_processes(char *name, int type, const char **procs[], int *nprocs)
{
  int i,n;

  struct {
    char *name;
    char *description;
    void (*init) (int type, const char **procs[], int *nprocs);
  } param_list[] = {
    {"dead_only","Non-lving components", eco_processes_dead},
    {"porewater_age","Porewater age tracer", eco_processes_pore},
    {"optics_only","Optical processes only",   eco_processes_opt},
    {"gas_only","Gas exchange processes only",   eco_processes_gas},
    {"standard","Standard ecology processes",   eco_processes_std},
    {"estuary", "Estuarine ecology processes", eco_processes_est},
    {"gbr4", "GBR4 ecology processes", eco_processes_gbr4},
    {"gbr4_seagrass", "GBR4 + seagrass processes", eco_processes_gbr4_seagrass},
    {"gbr4_coral", "GBR4 + coral bleaching processes", eco_processes_gbr4_coral},
    {"gbr4_small", "GBR4 SMALL ecology processes", eco_processes_gbr4_small},
    {"BGC2p0", "BGC 2.0 processes", eco_processes_bgc2p0},
    {"BGC3p1", "BGC 3.1 processes", eco_processes_bgc3p1},
    {"TASSE1p0", "TASSE 1p0 processes", eco_processes_TASSE1p0},

    {NULL, NULL, NULL}
  };
  void (*init) (int type, const char **procs[], int *nprocs) = NULL;

  if (strcmp(name, "auto") == 0) {
    eco_processes_auto(type, procs, nprocs);
  } else {
    for (i = 0; (init == NULL) && param_list[i].name; ++i) {
      if (strcasecmp(name, param_list[i].name) == 0) {
	init = param_list[i].init;
      }
    }
    if (init != NULL) {
      init(type, procs, nprocs);
    } else
      return(1);
  }
  return(0);
}

/*-------------------------------------------------------------------*/
/* Writes the ecological parameters to file                          */
/*-------------------------------------------------------------------*/
void write_eco_process(ecology *e)
{
  int type;
  char buf[MAXSTRLEN];
  FILE *fp;

  sprintf(buf, "processes_%s.prm", e->processfname);
  if ((fp = fopen(buf, "w")) == NULL)
    return;

  for (type = 0; type < N_EPROCESS_TYPES; ++type) {
    int n = e->npr[type];
    int i;

    if (n > 0) {
      fprintf(fp, " %s  {\n", eprocess_type_tags[type]);
      for (i = 0; i < n; ++i)
	fprintf(fp, "    %s\n", e->processes[type][i]->fullname);
      fprintf(fp, "  }\n");
    }
  }

  fclose(fp);
}

/* END write_eco_process()                                           */
/*-------------------------------------------------------------------*/


// END OF FILE
