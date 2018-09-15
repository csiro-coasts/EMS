/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: lib/ecology/allprocesses.c
 *  
 *  Description:
 *  List of all available ecological processes.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: allprocesses.c 5945 2018-09-13 21:56:39Z bai155 $
 *
 */

#include "stringtable.h"
#include "eprocess.h"

/* All process header files.
 */
#include "process_library/diffusion_epi.h"
#include "process_library/dinoflagellate_grow_wc.h"
#include "process_library/dinoflagellate_mortality_sed.h"
#include "process_library/anm_epi.h"
#include "process_library/light_epi.h"
#include "process_library/light_sed.h"
#include "process_library/light_wc.h"
#include "process_library/macroalgae_grow_epi.h"
#include "process_library/macroalgae_mortality_epi.h"
#include "process_library/massbalance_epi.h"
#include "process_library/massbalance_sed.h"
#include "process_library/massbalance_wc.h"
#include "process_library/microphytobenthos_grow_sed.h"
#include "process_library/microphytobenthos_grow_wc.h"
#include "process_library/trichodesmium_grow_wc.h"
#include "process_library/trichodesmium_spectral_grow_wc.h"
#include "process_library/microphytobenthos_mortality_sed.h"
#include "process_library/trichodesmium_mortality_sed.h"
#include "process_library/trichodesmium_mortality_wc.h"
#include "process_library/moldiff.h"
#include "process_library/nitrification_denitrification_sed.h"
#include "process_library/nitrification_wc.h"
#include "process_library/nodularia_grow_wc.h"
#include "process_library/nodularia_mortality_wc.h"
#include "process_library/nodularia_mortality_sed.h"
#include "process_library/oxygen_exchange_wc.h"
#include "process_library/oxygen_exchange_sed.h"
#include "process_library/p_adsorption_sed.h"
#include "process_library/p_adsorption_wc.h"
#include "process_library/phytoplankton_grow_wc.h"
#include "process_library/phytoplankton_mortality_sed.h"
#include "process_library/remineralization.h"
#include "process_library/seagrass_grow_epi.h"
#include "process_library/seagrass_mortality_epi.h"
#include "process_library/tfactor.h"
#include "process_library/values_common.h"
#include "process_library/values_common_epi.h"
#include "process_library/viscosity.h"
#include "process_library/anm_wc.h"
#include "process_library/zooplankton_large_grow_wc.h"
#include "process_library/zooplankton_small_grow_wc.h"
#include "process_library/zooplankton_mortality_wc.h"
#include "process_library/zooplankton_mortality_sed.h"
#include "process_library/chlorophyll_normalized_wc.h"
#include "process_library/dinoflagellate_diel_grow_wc.h"
#include "process_library/dinoflagellate_diel_mortality_sed.h"
#include "process_library/microphytobenthos_diel_grow_sed.h"
#include "process_library/microphytobenthos_diel_grow_wc.h"
#include "process_library/microphytobenthos_diel_mortality_sed.h"
#include "process_library/phytoplankton_diel_grow_wc.h"
#include "process_library/phytoplankton_diel_mortality_sed.h"
#include "process_library/zooplankton_large_diel_grow_wc.h"
#include "process_library/zooplankton_small_diel_grow_wc.h"
#include "process_library/light_spectral_wc.h"
#include "process_library/light_spectral2par_epi.h"
#include "process_library/microphytobenthos_spectral_grow_sed.h"
#include "process_library/microphytobenthos_spectral_grow_wc.h"
#include "process_library/microphytobenthos_spectral_mortality_sed.h"
#include "process_library/phytoplankton_spectral_grow_wc.h"
#include "process_library/phytoplankton_spectral_mortality_sed.h"
#include "process_library/phytoplankton_spectral_mortality_wc.h"
#include "process_library/zooplankton_small_spectral_grow_wc.h"
#include "process_library/zooplankton_large_spectral_grow_wc.h"
#include "process_library/zooplankton_large_carnivore_spectral_grow_wc.h"
#include "process_library/dinoflagellate_spectral_grow_wc.h"
#include "process_library/dinoflagellate_spectral_mortality_sed.h"
#include "process_library/seagrass_spectral_grow_epi.h"
#include "process_library/seagrass_spectral_grow_uq_epi.h"
#include "process_library/seagrass_spectral_grow_proto_epi.h"
#include "process_library/macroalgae_spectral_grow_epi.h"
#include "process_library/light_spectral_epi.h"
#include "process_library/light_spectral_sed.h"
#include "process_library/light_spectral_uq_epi.h"
#include "process_library/light_spectral_proto_epi.h"
#include "process_library/seagrass_spectral_mortality_epi.h"
#include "process_library/seagrass_spectral_mortality_proto_epi.h"
#include "process_library/coral_spectral_grow_epi.h"
#include "process_library/coral_spectral_grow_bleach_epi.h"
#include "process_library/coral_spectral_carb_epi.h"
#include "process_library/age_wc.h"
#include "process_library/macroalgae_grow_wc.h"
#include "process_library/macroalgae_mortality_wc.h"
#include "process_library/salmon_waste.h"
#include "process_library/filter_feeder_wc.h"
#include "process_library/filter_feeder_epi.h"
#include "process_library/variable_parameter.h"
#include "process_library/recom_extras.h"
#include "process_library/carbon_leak_sed.h"



/* Do NOT move the lines below */

/*INCLUDEMARKER*/
#include "process_library/light_optical_wc.h"

/*carbon chemistry Mathieu Mongin 18/05/2010*/
#include "process_library/carbon_chemistry_wc.h"
/*#include "process_library/kd_define.h"*/
#include "process_library/co2_exchange_wc.h"
#include "process_library/coral_growth.h"
#include "process_library/lyngbya_wc.h"
#include "process_library/gas_exchange_wc.h"
#include "process_library/gas_exchange_epi.h"

/* surface growing macroalgae spectral*/

#include "process_library/macroalgae_spectral_grow_wc.h"
#include "process_library/macroalgae_spectral_mortality_wc.h"


/** List of all processes in the process library.
 *
 * Following is the entry definition from "eprocess.h":
 *
 * typedef struct {
 *     char* name;
 *     PROCESSTYPE type;
 *     int nparams;
 *     int delayed;
 *     process_initfn init;
 *     process_initfn postinit;
 *     process_initfn destroy;
 *     process_calcfn precalc;
 *     process_calcfn calc;
 *     process_calcfn postcalc;
 * } process_entry;
 */
eprocess_entry eprocesslist[] = {
    {"diffusion_epi", PT_EPI, 0, 0, diffusion_epi_init, NULL, diffusion_epi_destroy, NULL, diffusion_epi_calc, NULL},
    {"dinoflagellate_grow_wc", PT_WC, 0, 0, dinoflagellate_grow_wc_init, dinoflagellate_grow_wc_postinit, dinoflagellate_grow_wc_destroy, dinoflagellate_grow_wc_precalc, dinoflagellate_grow_wc_calc, dinoflagellate_grow_wc_postcalc},
    {"dinoflagellate_mortality_sed", PT_SED, 0, 0, dinoflagellate_mortality_sed_init, dinoflagellate_mortality_sed_postinit, dinoflagellate_mortality_sed_destroy, dinoflagellate_mortality_sed_precalc, dinoflagellate_mortality_sed_calc, dinoflagellate_mortality_sed_postcalc},
    {"anm_epi", PT_EPI, 0, 0, anm_epi_init, NULL, anm_epi_destroy, NULL, anm_epi_calc, NULL},
    {"light_epi", PT_EPI, 0, 0, light_epi_init, NULL, light_epi_destroy, light_epi_precalc, NULL, NULL},
    {"light_spectral2par_epi", PT_EPI, 0, 0, light_spectral2par_epi_init, NULL, light_spectral2par_epi_destroy, light_spectral2par_epi_precalc, NULL, NULL},
    {"light_sed", PT_SED, 0, 0, light_sed_init, NULL, light_sed_destroy, light_sed_precalc, NULL, NULL},
    {"light_spectral_sed", PT_SED, -1, 0, light_spectral_sed_init, light_spectral_sed_postinit, light_spectral_sed_destroy, light_spectral_sed_precalc, NULL, NULL},
    {"light_wc", PT_WC, 0, 0, light_wc_init, light_wc_postinit, light_wc_destroy, light_wc_precalc, NULL, NULL},
    {"light_spectral_wc", PT_WC, -1, 0, light_spectral_wc_init, light_spectral_wc_postinit, light_spectral_wc_destroy, light_spectral_wc_precalc, NULL, light_spectral_wc_postcalc},
    {"macroalgae_grow_epi", PT_EPI, 0, 0, macroalgae_grow_epi_init, macroalgae_grow_epi_postinit, macroalgae_grow_epi_destroy, macroalgae_grow_epi_precalc, macroalgae_grow_epi_calc, macroalgae_grow_epi_postcalc},
    {"macroalgae_mortality_epi", PT_EPI, 0, 0, macroalgae_mortality_epi_init, NULL, macroalgae_mortality_epi_destroy, macroalgae_mortality_epi_precalc, macroalgae_mortality_epi_calc, macroalgae_mortality_epi_postcalc},
    {"massbalance_epi", PT_EPI, 0, 1, massbalance_epi_init, NULL, massbalance_epi_destroy, massbalance_epi_precalc, NULL, massbalance_epi_postcalc},
    {"massbalance_sed", PT_SED, 0, 0, massbalance_sed_init, NULL, massbalance_sed_destroy, massbalance_sed_precalc, NULL, massbalance_sed_postcalc},
    {"massbalance_wc", PT_WC, 0, 0, massbalance_wc_init, NULL, massbalance_wc_destroy, massbalance_wc_precalc, NULL, massbalance_wc_postcalc},
    {"microphytobenthos_grow_sed", PT_SED, 0, 0, microphytobenthos_grow_sed_init, microphytobenthos_grow_sed_postinit, microphytobenthos_grow_sed_destroy, microphytobenthos_grow_sed_precalc, microphytobenthos_grow_sed_calc, microphytobenthos_grow_sed_postcalc},
    {"microphytobenthos_grow_wc", PT_WC, 0, 0, microphytobenthos_grow_wc_init, microphytobenthos_grow_wc_postinit, microphytobenthos_grow_wc_destroy, microphytobenthos_grow_wc_precalc, microphytobenthos_grow_wc_calc, microphytobenthos_grow_wc_postcalc},
    {"microphytobenthos_mortality_sed", PT_SED, 0, 0, microphytobenthos_mortality_sed_init, NULL, microphytobenthos_mortality_sed_destroy, microphytobenthos_mortality_sed_precalc, microphytobenthos_mortality_sed_calc, microphytobenthos_mortality_sed_postcalc},
    {"moldiff", PT_GEN, 0, 0, moldiff_init, NULL, moldiff_destroy, moldiff_precalc, NULL, NULL},
    {"nitrification_denitrification_sed", PT_SED, 0, 0, nitrification_denitrification_sed_init, nitrification_denitrification_sed_postinit, nitrification_denitrification_sed_destroy, nitrification_denitrification_sed_precalc, nitrification_denitrification_sed_calc, nitrification_denitrification_sed_postcalc},
    {"nitrification_wc", PT_WC, 0, 0, nitrification_wc_init, NULL, nitrification_wc_destroy, nitrification_wc_precalc, nitrification_wc_calc, nitrification_wc_postcalc},
    {"nodularia_grow_wc", PT_WC, 0, 0, nodularia_grow_wc_init, nodularia_grow_wc_postinit, nodularia_grow_wc_destroy, nodularia_grow_wc_precalc, nodularia_grow_wc_calc, nodularia_grow_wc_postcalc},
    {"nodularia_mortality_wc", PT_WC, 0, 0, nodularia_mortality_wc_init, NULL, nodularia_mortality_wc_destroy, nodularia_mortality_wc_precalc, nodularia_mortality_wc_calc, nodularia_mortality_wc_postcalc},
    {"nodularia_mortality_sed", PT_SED, 0, 0, nodularia_mortality_sed_init, nodularia_mortality_sed_postinit, nodularia_mortality_sed_destroy, nodularia_mortality_sed_precalc, nodularia_mortality_sed_calc, nodularia_mortality_sed_postcalc},
    {"oxygen_exchange_wc", PT_WC, 0, 0, oxygen_exchange_wc_init, NULL, oxygen_exchange_wc_destroy, NULL, oxygen_exchange_wc_calc, oxygen_exchange_wc_postcalc},
    {"oxygen_exchange_sed", PT_SED, 0, 0, oxygen_exchange_sed_init, NULL, oxygen_exchange_sed_destroy, NULL, oxygen_exchange_sed_calc, oxygen_exchange_sed_postcalc},
    {"phytoplankton_grow_wc", PT_WC, 1, 0, phytoplankton_grow_wc_init, phytoplankton_grow_wc_postinit, phytoplankton_grow_wc_destroy, phytoplankton_grow_wc_precalc, phytoplankton_grow_wc_calc, phytoplankton_grow_wc_postcalc},
    {"phytoplankton_mortality_sed", PT_SED, 1, 0, phytoplankton_mortality_sed_init, phytoplankton_mortality_sed_postinit, phytoplankton_mortality_sed_destroy, phytoplankton_mortality_sed_precalc, phytoplankton_mortality_sed_calc, phytoplankton_mortality_sed_postcalc},
    {"p_adsorption_sed", PT_SED, 0, 0, p_adsorption_sed_init, p_adsorption_sed_postinit, p_adsorption_sed_destroy, p_adsorption_sed_precalc, p_adsorption_sed_calc, p_adsorption_sed_postcalc},
    {"p_adsorption_wc", PT_WC, 0, 0, p_adsorption_wc_init, p_adsorption_wc_postinit, p_adsorption_wc_destroy, p_adsorption_wc_precalc, p_adsorption_wc_calc, p_adsorption_wc_postcalc},
    {"remineralization", PT_GEN, 0, 0, remineralization_init, remineralization_postinit, remineralization_destroy, remineralization_precalc, remineralization_calc, remineralization_postcalc},
    {"seagrass_grow_epi", PT_EPI, 0, 0, seagrass_grow_epi_init, seagrass_grow_epi_postinit, seagrass_grow_epi_destroy, seagrass_grow_epi_precalc, seagrass_grow_epi_calc, seagrass_grow_epi_postcalc},
    {"seagrass_mortality_epi", PT_EPI, 0, 0, seagrass_mortality_epi_init, NULL, seagrass_mortality_epi_destroy, seagrass_mortality_epi_precalc, seagrass_mortality_epi_calc, seagrass_mortality_epi_postcalc},
    {"tfactor", PT_GEN, 0, 0, tfactor_init, NULL, tfactor_destroy, tfactor_precalc, NULL, NULL},
    {"tfactor_epi", PT_EPI, 0, 0, tfactor_init, NULL, tfactor_destroy, tfactor_precalc, NULL, NULL},
    {"recom_extras", PT_GEN, 0, 0, recom_extras_init, recom_extras_postinit, recom_extras_destroy, NULL, NULL, recom_extras_postcalc},
    {"carbon_leak_sed", PT_SED, 0, 0, carbon_leak_sed_init, NULL, carbon_leak_sed_destroy, NULL, carbon_leak_sed_calc, carbon_leak_sed_postcalc},
    {"values_common", PT_GEN, 0, 0, values_common_init, values_common_post_init, values_common_destroy, values_common_precalc, NULL, values_common_postcalc},
    {"values_common_epi", PT_EPI, 0, 0, values_common_epi_init, NULL, values_common_epi_destroy, values_common_epi_precalc, NULL, values_common_epi_postcalc},
    {"viscosity", PT_GEN, 0, 0, viscosity_init, NULL, viscosity_destroy, viscosity_precalc, NULL, NULL},
    {"anm_wc", PT_WC, 0, 0, anm_wc_init, NULL, anm_wc_destroy, anm_wc_precalc, anm_wc_calc, NULL},
    {"zooplankton_large_grow_wc", PT_WC, 0, 0, zooplankton_large_grow_wc_init, zooplankton_large_grow_wc_postinit, zooplankton_large_grow_wc_destroy, zooplankton_large_grow_wc_precalc, zooplankton_large_grow_wc_calc, zooplankton_large_grow_wc_postcalc},
    {"zooplankton_small_grow_wc", PT_WC, 0, 0, zooplankton_small_grow_wc_init, zooplankton_small_grow_wc_postinit, zooplankton_small_grow_wc_destroy, zooplankton_small_grow_wc_precalc, zooplankton_small_grow_wc_calc, zooplankton_small_grow_wc_postcalc},
    {"chlorophyll_normalized_wc",PT_WC,0,1,chlorophyll_normalized_wc_init,chlorophyll_normalized_wc_postinit,chlorophyll_normalized_wc_destroy,chlorophyll_normalized_wc_precalc,chlorophyll_normalized_wc_calc,chlorophyll_normalized_wc_postcalc},
    {"dinoflagellate_diel_grow_wc", PT_WC, 0, 0, dinoflagellate_diel_grow_wc_init, dinoflagellate_diel_grow_wc_postinit, dinoflagellate_diel_grow_wc_destroy, dinoflagellate_diel_grow_wc_precalc, dinoflagellate_diel_grow_wc_calc, dinoflagellate_diel_grow_wc_postcalc},
    {"dinoflagellate_diel_mortality_sed", PT_SED, 0, 0, dinoflagellate_diel_mortality_sed_init, dinoflagellate_diel_mortality_sed_postinit, dinoflagellate_diel_mortality_sed_destroy, dinoflagellate_diel_mortality_sed_precalc, dinoflagellate_diel_mortality_sed_calc, dinoflagellate_diel_mortality_sed_postcalc},
    {"dinoflagellate_spectral_grow_wc", PT_WC, 0, 0, dinoflagellate_spectral_grow_wc_init, dinoflagellate_spectral_grow_wc_postinit, dinoflagellate_spectral_grow_wc_destroy, dinoflagellate_spectral_grow_wc_precalc, dinoflagellate_spectral_grow_wc_calc, dinoflagellate_spectral_grow_wc_postcalc},
    {"dinoflagellate_spectral_mortality_sed", PT_SED, 0, 0, dinoflagellate_spectral_mortality_sed_init, dinoflagellate_spectral_mortality_sed_postinit, dinoflagellate_spectral_mortality_sed_destroy, dinoflagellate_spectral_mortality_sed_precalc, dinoflagellate_spectral_mortality_sed_calc, dinoflagellate_spectral_mortality_sed_postcalc},
    {"microphytobenthos_diel_grow_sed", PT_SED, 0, 0, microphytobenthos_diel_grow_sed_init, microphytobenthos_diel_grow_sed_postinit, microphytobenthos_diel_grow_sed_destroy, microphytobenthos_diel_grow_sed_precalc, microphytobenthos_diel_grow_sed_calc, microphytobenthos_diel_grow_sed_postcalc},
    {"microphytobenthos_spectral_grow_sed", PT_SED, 0, 0, microphytobenthos_spectral_grow_sed_init, microphytobenthos_spectral_grow_sed_postinit, microphytobenthos_spectral_grow_sed_destroy, microphytobenthos_spectral_grow_sed_precalc, microphytobenthos_spectral_grow_sed_calc, microphytobenthos_spectral_grow_sed_postcalc},
    {"microphytobenthos_diel_grow_wc", PT_WC, 0, 0, microphytobenthos_diel_grow_wc_init, microphytobenthos_diel_grow_wc_postinit, microphytobenthos_diel_grow_wc_destroy, microphytobenthos_diel_grow_wc_precalc, microphytobenthos_diel_grow_wc_calc, microphytobenthos_diel_grow_wc_postcalc},
    {"microphytobenthos_spectral_grow_wc", PT_WC, 0, 0, microphytobenthos_spectral_grow_wc_init, microphytobenthos_spectral_grow_wc_postinit, microphytobenthos_spectral_grow_wc_destroy, microphytobenthos_spectral_grow_wc_precalc, microphytobenthos_spectral_grow_wc_calc, microphytobenthos_spectral_grow_wc_postcalc},
    {"microphytobenthos_diel_mortality_sed", PT_SED, 0, 0, microphytobenthos_diel_mortality_sed_init, NULL, microphytobenthos_diel_mortality_sed_destroy, microphytobenthos_diel_mortality_sed_precalc, microphytobenthos_diel_mortality_sed_calc, microphytobenthos_diel_mortality_sed_postcalc},
    {"microphytobenthos_spectral_mortality_sed", PT_SED, 0, 0, microphytobenthos_spectral_mortality_sed_init, microphytobenthos_spectral_mortality_sed_postinit, microphytobenthos_spectral_mortality_sed_destroy, microphytobenthos_spectral_mortality_sed_precalc, microphytobenthos_spectral_mortality_sed_calc, microphytobenthos_spectral_mortality_sed_postcalc},
    {"phytoplankton_diel_grow_wc", PT_WC, 1, 0, phytoplankton_diel_grow_wc_init, phytoplankton_diel_grow_wc_postinit, phytoplankton_diel_grow_wc_destroy, phytoplankton_diel_grow_wc_precalc, phytoplankton_diel_grow_wc_calc, phytoplankton_diel_grow_wc_postcalc},
    {"phytoplankton_spectral_grow_wc", PT_WC, 1, 0, phytoplankton_spectral_grow_wc_init, phytoplankton_spectral_grow_wc_postinit, phytoplankton_spectral_grow_wc_destroy, phytoplankton_spectral_grow_wc_precalc, phytoplankton_spectral_grow_wc_calc, phytoplankton_spectral_grow_wc_postcalc},
    {"phytoplankton_diel_mortality_sed", PT_SED, 1, 0, phytoplankton_diel_mortality_sed_init, phytoplankton_diel_mortality_sed_postinit, phytoplankton_diel_mortality_sed_destroy, phytoplankton_diel_mortality_sed_precalc, phytoplankton_diel_mortality_sed_calc, phytoplankton_diel_mortality_sed_postcalc},
    {"phytoplankton_spectral_mortality_sed", PT_SED, 1, 0, phytoplankton_spectral_mortality_sed_init, phytoplankton_spectral_mortality_sed_postinit, phytoplankton_spectral_mortality_sed_destroy, phytoplankton_spectral_mortality_sed_precalc, phytoplankton_spectral_mortality_sed_calc, phytoplankton_spectral_mortality_sed_postcalc},
    {"phytoplankton_spectral_mortality_wc", PT_WC, 1, 0, phytoplankton_spectral_mortality_wc_init, phytoplankton_spectral_mortality_wc_postinit, phytoplankton_spectral_mortality_wc_destroy, phytoplankton_spectral_mortality_wc_precalc, phytoplankton_spectral_mortality_wc_calc, phytoplankton_spectral_mortality_wc_postcalc},
    {"zooplankton_large_diel_grow_wc", PT_WC, 0, 0, zooplankton_large_diel_grow_wc_init, zooplankton_large_diel_grow_wc_postinit, zooplankton_large_diel_grow_wc_destroy, zooplankton_large_diel_grow_wc_precalc, zooplankton_large_diel_grow_wc_calc, zooplankton_large_diel_grow_wc_postcalc},
    {"zooplankton_small_diel_grow_wc", PT_WC, 0, 0, zooplankton_small_diel_grow_wc_init, zooplankton_small_diel_grow_wc_postinit, zooplankton_small_diel_grow_wc_destroy, zooplankton_small_diel_grow_wc_precalc, zooplankton_small_diel_grow_wc_calc, zooplankton_small_diel_grow_wc_postcalc},
    {"zooplankton_small_spectral_grow_wc", PT_WC, 0, 0, zooplankton_small_spectral_grow_wc_init, zooplankton_small_spectral_grow_wc_postinit, zooplankton_small_spectral_grow_wc_destroy, zooplankton_small_spectral_grow_wc_precalc, zooplankton_small_spectral_grow_wc_calc, zooplankton_small_spectral_grow_wc_postcalc},
    {"zooplankton_large_spectral_grow_wc", PT_WC, 0, 0, zooplankton_large_spectral_grow_wc_init, zooplankton_large_spectral_grow_wc_postinit, zooplankton_large_spectral_grow_wc_destroy, zooplankton_large_spectral_grow_wc_precalc, zooplankton_large_spectral_grow_wc_calc, zooplankton_large_spectral_grow_wc_postcalc},
    {"zooplankton_large_carnivore_spectral_grow_wc", PT_WC, 0, 0, zooplankton_large_carnivore_spectral_grow_wc_init, zooplankton_large_carnivore_spectral_grow_wc_postinit, zooplankton_large_carnivore_spectral_grow_wc_destroy, zooplankton_large_carnivore_spectral_grow_wc_precalc, zooplankton_large_carnivore_spectral_grow_wc_calc, zooplankton_large_carnivore_spectral_grow_wc_postcalc},
    {"seagrass_spectral_grow_epi", PT_EPI, 1, 0, seagrass_spectral_grow_epi_init, seagrass_spectral_grow_epi_postinit, seagrass_spectral_grow_epi_destroy, seagrass_spectral_grow_epi_precalc, seagrass_spectral_grow_epi_calc, seagrass_spectral_grow_epi_postcalc},
    {"seagrass_spectral_grow_uq_epi", PT_EPI, 1, 0, seagrass_spectral_grow_uq_epi_init, seagrass_spectral_grow_uq_epi_postinit, seagrass_spectral_grow_uq_epi_destroy, seagrass_spectral_grow_uq_epi_precalc, seagrass_spectral_grow_uq_epi_calc, seagrass_spectral_grow_uq_epi_postcalc},
    {"seagrass_spectral_grow_proto_epi", PT_EPI, 1, 0, seagrass_spectral_grow_proto_epi_init, seagrass_spectral_grow_proto_epi_postinit, seagrass_spectral_grow_proto_epi_destroy, seagrass_spectral_grow_proto_epi_precalc, seagrass_spectral_grow_proto_epi_calc, seagrass_spectral_grow_proto_epi_postcalc},
    {"macroalgae_spectral_grow_epi", PT_EPI, 0, 0, macroalgae_spectral_grow_epi_init, macroalgae_spectral_grow_epi_postinit, macroalgae_spectral_grow_epi_destroy, macroalgae_spectral_grow_epi_precalc, macroalgae_spectral_grow_epi_calc, macroalgae_spectral_grow_epi_postcalc},
    {"light_spectral_epi", PT_EPI, 0, 0, light_spectral_epi_init, light_spectral_epi_postinit, light_spectral_epi_destroy, light_spectral_epi_precalc, NULL, NULL},
    {"light_spectral_uq_epi", PT_EPI, 0, 0, light_spectral_uq_epi_init, light_spectral_uq_epi_postinit, light_spectral_uq_epi_destroy, light_spectral_uq_epi_precalc, NULL, light_spectral_uq_epi_postcalc},
    {"light_spectral_proto_epi", PT_EPI, 0, 0, light_spectral_proto_epi_init, light_spectral_proto_epi_postinit, light_spectral_proto_epi_destroy, light_spectral_proto_epi_precalc, NULL, NULL},
    {"seagrass_spectral_mortality_epi", PT_EPI, 1, 0, seagrass_spectral_mortality_epi_init, NULL, seagrass_spectral_mortality_epi_destroy, seagrass_spectral_mortality_epi_precalc, seagrass_spectral_mortality_epi_calc, seagrass_spectral_mortality_epi_postcalc},
    {"seagrass_spectral_mortality_proto_epi", PT_EPI, 1, 0, seagrass_spectral_mortality_proto_epi_init, NULL, seagrass_spectral_mortality_proto_epi_destroy, seagrass_spectral_mortality_proto_epi_precalc, seagrass_spectral_mortality_proto_epi_calc, seagrass_spectral_mortality_proto_epi_postcalc},
    {"coral_spectral_grow_epi", PT_EPI, 0, 0, coral_spectral_grow_epi_init, coral_spectral_grow_epi_postinit, coral_spectral_grow_epi_destroy, coral_spectral_grow_epi_precalc, coral_spectral_grow_epi_calc, coral_spectral_grow_epi_postcalc},
    {"coral_spectral_grow_bleach_epi", PT_EPI, 0, 0, coral_spectral_grow_bleach_epi_init, coral_spectral_grow_bleach_epi_postinit, coral_spectral_grow_bleach_epi_destroy, coral_spectral_grow_bleach_epi_precalc, coral_spectral_grow_bleach_epi_calc, coral_spectral_grow_bleach_epi_postcalc},
    {"coral_spectral_carb_epi", PT_EPI, -1, 0, coral_spectral_carb_epi_init, coral_spectral_carb_epi_postinit, coral_spectral_carb_epi_destroy, coral_spectral_carb_epi_precalc, coral_spectral_carb_epi_calc, coral_spectral_carb_epi_postcalc},
    {"filter_feeder_wc", PT_WC, 0, 0, filter_feeder_wc_init, filter_feeder_wc_postinit, filter_feeder_wc_destroy, filter_feeder_wc_precalc, filter_feeder_wc_calc, filter_feeder_wc_postcalc},
    {"filter_feeder_epi", PT_EPI, 0, 0, filter_feeder_epi_init, filter_feeder_epi_postinit, filter_feeder_epi_destroy, filter_feeder_epi_precalc, filter_feeder_epi_calc, filter_feeder_epi_postcalc},
    {"age_wc", PT_WC, 0, 0, age_wc_init, NULL, age_wc_destroy, age_wc_precalc, age_wc_calc, NULL},
/*    {"light_wc_gradient", PT_WC, 0, 0, light_wc_gradient_init, light_wc_gradient_postinit, light_wc_gradient_destroy, light_wc_gradient_precalc, NULL, NULL},
*/

    /*ENTRYMARKER*/
    {"light_optical_wc",PT_WC,0,0,light_optical_wc_init,light_optical_wc_postinit,light_optical_wc_destroy,light_optical_wc_precalc,light_optical_wc_calc,light_optical_wc_postcalc},


/*  ADD new processes here below this line, ending with a comma */

    {"zooplankton_mortality_wc", PT_WC, 1, 0, zooplankton_mortality_wc_init, NULL, zooplankton_mortality_wc_destroy, zooplankton_mortality_wc_precalc, zooplankton_mortality_wc_calc, zooplankton_mortality_wc_postcalc},
	{"zooplankton_mortality_sed", PT_SED, 1, 0, zooplankton_mortality_sed_init, NULL, zooplankton_mortality_sed_destroy, zooplankton_mortality_sed_precalc, zooplankton_mortality_sed_calc, zooplankton_mortality_sed_postcalc},
    {"trichodesmium_grow_wc", PT_WC, 0, 0, trichodesmium_grow_wc_init, trichodesmium_grow_wc_postinit, trichodesmium_grow_wc_destroy, trichodesmium_grow_wc_precalc, trichodesmium_grow_wc_calc, trichodesmium_grow_wc_postcalc},
{"trichodesmium_spectral_grow_wc", PT_WC, 0, 0, trichodesmium_spectral_grow_wc_init, trichodesmium_spectral_grow_wc_postinit, trichodesmium_spectral_grow_wc_destroy, trichodesmium_spectral_grow_wc_precalc, trichodesmium_spectral_grow_wc_calc, trichodesmium_spectral_grow_wc_postcalc},
    {"trichodesmium_mortality_wc", PT_WC, 0, 0, trichodesmium_mortality_wc_init, trichodesmium_mortality_wc_postinit, trichodesmium_mortality_wc_destroy, trichodesmium_mortality_wc_precalc, trichodesmium_mortality_wc_calc, trichodesmium_mortality_wc_postcalc},
    {"trichodesmium_mortality_sed", PT_SED, 0, 0, trichodesmium_mortality_sed_init, trichodesmium_mortality_sed_postinit, trichodesmium_mortality_sed_destroy, trichodesmium_mortality_sed_precalc, trichodesmium_mortality_sed_calc, trichodesmium_mortality_sed_postcalc},
    {"gas_exchange_wc", PT_WC, 2, 0, gas_exchange_wc_init, NULL, gas_exchange_wc_destroy, gas_exchange_wc_precalc, gas_exchange_wc_calc, gas_exchange_wc_postcalc},
   {"gas_exchange_epi", PT_EPI, 0, 0, gas_exchange_epi_init, NULL, gas_exchange_epi_destroy, gas_exchange_epi_precalc, gas_exchange_epi_calc, gas_exchange_epi_postcalc},
    {"carbon_chemistry_wc", PT_GEN, 0, 0, carbon_chemistry_wc_init, carbon_chemistry_wc_postinit, carbon_chemistry_wc_destroy, carbon_chemistry_wc_precalc, carbon_chemistry_wc_calc, carbon_chemistry_wc_postcalc},
    {"co2_exchange_wc", PT_WC, 0, 0, co2_exchange_wc_init, NULL, co2_exchange_wc_destroy, co2_exchange_wc_precalc, co2_exchange_wc_calc, co2_exchange_wc_postcalc},
    {"coral_growth", PT_EPI, 0, 0, coral_growth_init, coral_growth_postinit,coral_growth_destroy,coral_growth_precalc, coral_growth_calc, coral_growth_postcalc},
    {"lyngbya_wc", PT_WC, 0, 0,lyngbya_wc_init, NULL, lyngbya_wc_destroy, lyngbya_wc_precalc,lyngbya_wc_calc,lyngbya_wc_postcalc},
    {"macroalgae_grow_wc", PT_WC, 0, 0, macroalgae_grow_wc_init, macroalgae_grow_wc_postinit, macroalgae_grow_wc_destroy, macroalgae_grow_wc_precalc, macroalgae_grow_wc_calc, macroalgae_grow_wc_postcalc},
    {"macroalgae_mortality_wc", PT_WC, 0, 0, macroalgae_mortality_wc_init, NULL, macroalgae_mortality_wc_destroy, macroalgae_mortality_wc_precalc, macroalgae_mortality_wc_calc, macroalgae_mortality_wc_postcalc},
    {"salmon_waste", PT_WC, 0, 0, salmon_waste_init, NULL, salmon_waste_destroy, salmon_waste_precalc, NULL, salmon_waste_postcalc},

    /*add macroalgae growing on the water column*/

    {"macroalgae_spectral_grow_wc", PT_WC, 0, 0, macroalgae_spectral_grow_wc_init, macroalgae_spectral_grow_wc_postinit, macroalgae_spectral_grow_wc_destroy, macroalgae_spectral_grow_wc_precalc, macroalgae_spectral_grow_wc_calc, macroalgae_spectral_grow_wc_postcalc},
    {"macroalgae_spectral_mortality_wc", PT_WC, 0, 0, macroalgae_spectral_mortality_wc_init, NULL, macroalgae_spectral_mortality_wc_destroy, macroalgae_spectral_mortality_wc_precalc, macroalgae_spectral_mortality_wc_calc, macroalgae_spectral_mortality_wc_postcalc},
    {"variable_parameter", PT_WC, 0, 0,variable_parameter_init,NULL,variable_parameter_destroy,variable_parameter_precalc,variable_parameter_calc,variable_parameter_postcalc}


};
int NEPROCESSES = sizeof(eprocesslist) / sizeof(eprocess_entry);

/** List of all tracers in the process library.
 *
 * Following is the entry definition from "eprocess.h":
 *
 * typedef struct {
 *     char* name;
 *     int diagn;
 * } diagnflag_entry;
 *
 * There is no reason to mantain this list other than to ensure consistency
 * between values of diagnostic flag in the host model and the ecology model.
 */
diagnflag_entry diagnflags[] = {
    /*
     * tracers
     */
    {"temp", 0},
    {"salt", 0},
    {"DetPL_N", 0},
    {"DetBL_N", 0},
    {"DetR_N", 0},
    {"DetR_C", 0},
    {"DetR_P", 0},
    {"DOR_N", 0},
    {"DOR_C", 0},
    {"DOR_P", 0},
    {"NH4", 0},
    {"NH4_pr", 1},
    {"NO3", 0},
    {"Oxygen", 0},
    {"Oxy_pr", 1},
    {"DIC", 0},
    {"DIP", 0},
    {"TN", 2},
    {"TC", 2},
    {"TP", 2},
    {"Kd", 2},
    {"PhyD_N", 0},
    {"PhyD_NR", 0},
    {"PhyD_PR", 0},
    {"PhyD_I", 0},
    {"PhyD_Chl", 0},
    {"PhyD_N_pr", 1},
    {"PhyD_N_gr", 1},
    {"PhyD_C", 0},
    {"Light", 2},
    {"Chl_a", 2},
    {"Oxy_sat", 2},
    {"MPB_N", 0},
    {"MPB_NR", 0},
    {"MPB_PR", 0},
    {"MPB_I", 0},
    {"MPB_Chl", 0},
    {"MPB_N_pr", 1},
    {"MPB_N_gr", 1},
    {"Tricho_N", 0},
    {"Tricho_PR", 0},
    {"Tricho_NR", 0},
    {"Tricho_I", 0},
    {"Tricho_Chl", 0},
    {"Tricho_sv", 0},
    {"Tricho_N_pr", 1},
    {"Tricho_N_gr", 1},
    {"Tricho_svel", 2},
    {"PhyN_N", 0},
    {"Nfix", 1},
    {"PhyN_N_pr", 1},
    {"PhyN_N_gr", 1},
    {"PhyS_N", 0},
    {"PhyS_NR", 0},
    {"PhyS_PR", 0},
    {"PhyS_I", 0},
    {"PhyS_Chl", 0},
    {"PhyS_N_pr", 1},
    {"PhyS_N_gr", 1},
    {"PhyL_N", 0},
    {"PhyL_NR", 0},
    {"PhyL_PR", 0},
    {"PhyL_I", 0},
    {"PhyL_Chl", 0},
    {"PhyL_N_pr", 1},
    {"PhyL_N_gr", 1},
    {"ZooL_N", 0},
    {"ZooL_N_gr", 1},
    {"ZooL_N_rm", 1},
    {"ZooS_N", 0},
    {"ZooS_N_gr", 1},
    {"ZooS_N_rm", 1},
    {"TSSUF", 0},
    {"TSSF", 0},
    {"PIP", 0},
    {"DIN", 2},
    {"Den_fl", 1},
    {"Den_eff", 1},
    {"PIPI", 0},
    /*
     * epibenthic variables
     */
    {"MA_N", 0},
    {"MA_N_pr", 1},
    {"MA_N_gr", 1},
    {"Epilight", 2},
    {"Epilightatt", 2},
    {"EpiTN", 2},
    {"EpiTP", 2},
    {"EpiTC", 2},
    {"SG_N", 0},
    {"SG_N_pr", 1},
    {"SG_N_gr", 1},
    {"EpiOxy_pr", 1},
    {"ustrcw_skin", 0},
    {"SGROOT_N", 0},
    {"Gnet", 0},
    {"Epilight_420", 2},
    {"SGH_N", 0},
    {"CS_N", 0},
    {"CS_Chl", 0},

    /*  ADD new tracer entries here above this line */
    {"Mud", 0},
    {"FineSed", 2},
    {"EFI", 2},
    {"P_Prod", 2},
    {"G_grazing", 2},
    {"alk", 0},
    {"PH", 0},
    {"dco2star", 2},
    {"CO23", 0},
    {"HCO3", 0},
    {"omega_ar", 2},
    {"CO2_flux", 2},
    {"at_420", 2},
    {"bt_420", 2},
    {"Zenith", 2},
    {"Light_420", 2},
    {"Kd_420", 2},
    {"Light_600", 2},
    {"Kd_600", 2}
};

int NDIAGNFLAGS = sizeof(diagnflags) / sizeof(diagnflag_entry);
