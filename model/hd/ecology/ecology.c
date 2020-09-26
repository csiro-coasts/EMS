/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/ecology/ecology.c
 *  
 *  Description:
 *  Interfaces the ecology library to SHOC
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: ecology.c 6376 2019-11-06 04:09:47Z bai155 $
 *
 */

#include <stdlib.h>
#include "hd.h"


#if defined(HAVE_ECOLOGY_MODULE)

#include <assert.h>
#include "einterface.h"
#include "ecology_tracer_defaults.h"

#define ECO_MAXNUMARGS 400

double sinterface_get_svel(void* model,char* name);

/*
 * FR: These should be moved into the ecology library
 */

/* Ecology 3D tracers */
const char *ECONAME3D[][2] = {
  {"Age",         "Tracer age"},
  {"source",      "Ageing zone"},
  {"temp_clim",   "Climatological temperature"},
  {"PhyL_N",      "Large Phytoplankton N"},
  {"PhyL_NR",      "Large Phytoplankton N reserve"},
  {"PhyL_PR",      "Large Phytoplankton P reserve"},
  {"PhyL_I",      "Large Phytoplankton I reserve"},
  {"PhyL_Chl",      "Large Phytoplankton chlorophyll"},
  {"PhyS_N",      "Small Phytoplankton N"},
  {"PhyS_NR",      "Small Phytoplankton N reserve"},
  {"PhyS_PR",      "Small Phytoplankton P reserve"},
  {"PhyS_I",      "Small Phytoplankton I reserve"},
  {"PhyS_Chl",      "Small Phytoplankton chlorophyll"},
  {"PhyL_sv",      "diatoms settling velocity"},
  {"MPB_N",       "Microphytobenthos N"},
  {"MPB_NR",       "Microphytobenthos N reserve"},
  {"MPB_PR",       "Microphytobenthos P reserve"},
  {"MPB_I",       "Microphytobenthos light reserve"},
  {"MPB_Chl",       "Microphytobenthos chlorophyll"},
  {"Zenith",       "Solar zenith"},
  {"ZooL_N",      "Large Zooplankton N"},
  {"ZooS_N",      "Small Zooplankton N"},
  {"PhyD_N",      "Dinoflagellate N"},
  {"PhyD_C",      "Dinoflagellate C"},
  {"Tricho_N",       "Trichodesmium Nitrogen"},
  {"Tricho_NR",       "Trichodesmium N reserve"},
  {"Tricho_PR",       "Trichodesmium P reserve"},
  {"Tricho_I",       "Trichodesmium I reserve"},
  {"Tricho_Chl",       "Trichodesmium chlorophyll"},
  {"Tricho_sv",       "Trichodesmium settling velocity"},
  {"ZooL_sv",       "Large Zooplankton settling velocity"},
  {"Phy_L_N2",    "Lyngbya"},
  {"NH4",         "Ammonia"},
  {"NO3",         "Nitrate"},
  {"DIP",         "Dissolved Inorganic Phosphorus"},
  {"DIC",         "Dissolved Inorganic Carbon"},
  {"DOR_C",       "Dissolved Organic Carbon"},
  {"DOR_N",       "Dissolved Organic Nitrogen"},
  {"DOR_P",       "Dissolved Organic Phosphorus"},
  {"PIP",         "Particulate Inorganic Phosphorus"},
  {"PIP_Dust",    "Particulate Inorganic Phosphorus Dust"},
  {"PIPI",        "Immobilised Particulate Inorganic Phosphorus"},
  {"PIPF",        "Flocculated Particulate Inorganic Phosphorus"},
  {"DetR_C",      "Refractory Detrital Carbon"},
  {"DetR_N",      "Refractory Detrital Nitrogen"},
  {"DetR_P",      "Refractory Detrital Phosphorus"},
  {"DetPL_N",     "Labile Detrital Nitrogen Plank"},
  {"DetBL_N",     "Labile Detrital Nitrogen Benthic"},
  {"Oxygen",      "Dissolved Oxygen"},
  {"Oxy_sat",     "Oxygen saturation percent"},
  {"Light",       "Av. light in layer"},
  {"PAR",         "Av. PAR in layer"},
  {"PAR_z",  "Downwelling PAR at top of layer"},
  {"K_heat",      "Vertical attenuation of heat"},
  {"Kd",          "Attenuation coefficient in layer"},
  {"TN",          "Total N"},
  {"TP",          "Total P"},
  {"TC",          "Total C"},
  {"EFI",         "Ecology Fine Inorganics"},
  {"DIN",         "Dissolved Inorganic Nitrogen"},
  {"Chl_a",       "Total Chlorophyll"},
  {"Chl_a_sum",   "Total Chlorophyll"},
  {"PhyL_N_pr",   "Large Phytoplankton net production"},
  {"PhyS_N_pr",   "Small Phytoplankton net production"},
  {"ZooL_N_rm",   "Large Zooplankton removal rate from Large Phytoplankton"},
  {"ZooS_N_rm",   "Small Zooplankton removal from Small Phytoplankton"},
  {"Den_fl",      "Denitrfication flux"},
  {"Den_eff",     "Denitrification effectiveness"},
  {"NH4_pr",      "Ammonia production"},
  {"MPB_N_pr",    "Microphytobenthos net production"},
  {"Tricho_N_pr", "Trichodesmium net production"},
  {"PhyL_N_gr",   "Large Phytoplankton growth rate"},
  {"PhyS_N_gr",   "Small Phytoplankton growth rate"},
  {"BOD",         "Biochemical Oxygen Demand"},
  {"COD",         "Chemical Oxygen Demand"},
  {"MPB_N_gr",    "Microphytobenthos growth rate"},
  {"ZooL_N_gr",   "Large Zooplankton growth rate"},
  {"ZooS_N_gr",   "Small Zooplankton growth rate"},
  {"PhyD_N_gr",   "Dinoflagellate growth rate"},
  {"Tricho_N_gr", "Trichodesmium growth rate"},
  {"Oxy_pr",      "Oxygen production"},
  {"PhyD_N_pr",   "Dinoflagellate net production"},
  {"dens",        "Density"},
  {"alk",         "Total alkalinity"},
  {"Nfix",        "N2 fixation"},
  {"Amm_fl",      "Anammox flux"},
  {"PH",          "PH"},
  {"CO32",        "Carbonate"},
  {"HCO3",        "Bicarbonate"},
  {"CO2_starair", "Atm CO2 (xco2*ff*atmpres)"}, 
  {"CO2_star",    "Ocean surface aquaeous CO2 concentration"},
  {"dco2star",    "Sea-air Delta CO2"},
  {"pco2surf",    "Oceanic pCO2"},
  {"dpCO2",       "Sea-air Delta pCO2"},
  {"omega_ca",    "Calcite saturation state"},
  {"omega_ar",    "Aragonite saturation state"},
  {"CO2_flux",    "Sea-air CO2 flux"},
  {"O2_flux",     "Sea-air O2 flux"},
  {"at_440",      "Absorption at 440 nm"},
  {"bt_550",      "Scattering at 550 nm"},  
  {"Kd_490",      "Vertical attenuation at 490 nm"},
  {"Turbidity",   "Simulated turbidity vs. bp_590 relationship"},
  {"Fluorescence","Simulated Fluorescence"},
  {"ap_670",      "Absorption at 670 nm minus clear water"},
  {"Phy_L_N2_fix","N2 fix rate Lyngbya"},
  {"xco2_in_air","Atmospheric pCO2"},
};

const int NUM_ECO_VARS_3D = ((int)(sizeof(ECONAME3D)/(2*sizeof(char*))));

/* Ecology 2D tracers */
const char *ECONAME2D[][2] = {
  {"Epilight",    "Light intensity above epibenthos"},
  {"Epilightatt", "Light attenuation in epibenthos"}, 
  {"EpiPAR",      "PAR light below epibenthos"},
  {"EpiPAR_sg",   "PAR light above seagrass"},
  {"EpiBOD",      "Total biochemical oxygen demand in epibenthos"},
  {"SG_N",        "Seagrass N"},
  {"SGH_N",       "Halophila N"},
  {"SGROOT_N",    "Seagrass root N"},
  {"SGHROOT_N",   "Halophila root N"},
  {"SGP_N",       "Posidonia N"},
  {"SGD_N",       "Deep seagrass N"},
  {"SGPROOT_N",   "Posidonia root N"},
  {"SGDROOT_N",   "Deep seagrass root N"},
  {"CS_N",        "Coral symbiont N"},
  {"CS_NR",       "Coral symbiont N reserves"},
  {"CS_PR",       "Coral symbiont P reserves"},
  {"CS_I",        "Coral symbiont energy reserves"},
  {"CS_Chl",      "Coral symbiont Chl"},
  {"CS_Xp",       "Coral symbiont Xan. photo."},
  {"CS_Xh",       "Coral symbiont Xan. heat."},
  {"CS_Qred",     "Coral symbiont red. RC"},
  {"CS_Qox",      "Coral symbiont ox. RC"},
  {"CS_Qi",       "Coral symbiont in. RC"},
  {"CS_RO",       "Coral symbiont reactive oxygen"},
  {"CH_N",        "Coral host N"},
  {"MA_N",        "Macroalgae N"},
  {"EpiTN",       "Total N in epibenthos"},
  {"EpiTP",       "Total P in epibenthos"},
  {"EpiTC",       "Total C in epibenthos"},
  {"SG_N_pr",     "Seagrass net production"},
  {"SGH_N_pr",    "Halophila net production"},
  {"MA_N_pr",     "Macroalgae net production"},
  {"SG_N_gr",     "Seagrass growth rate"},
  {"SGH_N_gr",    "Halophila growth rate"},
  {"MA_N_gr",     "Macroalgae growth rate"},
  {"EpiOxy_pr",   "Net Oxygen Production in epibenthos"},
  {"Coral_IN_up", "Coral inorganic uptake"},
  {"Coral_ON_up", "Coral organic uptake"},
  {"Gnet",        "Coral net calcification"},
  {"mucus",       "Coral mucus production"},
  {"CS_N_pr",     "Coral symbiont net production"},
  {"CH_N_pr",     "Coral host net production"},
  {"SG_shear_mort","Seagrass shear stress mort"},
  {"SGH_shear_mort","Halophila shear stress mort"},
  {"SGP_shear_mort","Posidonia shear stress mort"},
  {"SGD_shear_mort","Deep seagrass shear stress mort"},
  {"CS_bleach",      "Coral sym. expulsion rate"},
  {"CS_tempfunc", "Bleach T function"},
  {"OC3M","Modis-OC3M"},
  {"OC4Me","MERIS-OC4Me"},
  {"OC3V","VIIRS-OC3V"},
  {"TSSM","TSS from 645 nm (Petus et al., 2014)"},
  {"KD490M","Modis-KD490"},
  {"R_412","Remote-sensing reflectance @ 412 nm"},
  {"R_443","Remote-sensing reflectance @ 443 nm"},
  {"R_488","Remote-sensing reflectance @ 488 nm"},
  {"R_531","Remote-sensing reflectance @ 531 nm"},
  {"R_547","Remote-sensing reflectance @ 547 nm"},
  {"R_667","Remote-sensing reflectance @ 667 nm"},
  {"R_678","Remote-sensing reflectance @ 678 nm"},
  {"R_748","Remote-sensing reflectance @ 748 nm"},
  {"R_470","Remote-sensing reflectance @ 470 nm"},
  {"R_555","Remote-sensing reflectance @ 555 nm"},
  {"R_645","Remote-sensing reflectance @ 645 nm"},
  {"R_590","Remote-sensing reflectance @ 590 nm"},
  {"R_410","Remote-sensing reflectance @ 410 nm"},
  {"R_486","Remote-sensing reflectance @ 486 nm"},
  {"R_551","Remote-sensing reflectance @ 551 nm"},
  {"R_671","Remote-sensing reflectance @ 671 nm"},
  {"R_745","Remote-sensing reflectance @ 745 nm"},
  {"R_640","Remote-sensing reflectance @ 640 nm"},
  {"R_510","Remote-sensing reflectance @ 510 nm"},
  {"R_400","Remote-sensing reflectance @ 400 nm"},
  {"R_560","Remote-sensing reflectance @ 560 nm"},
  {"R_620","Remote-sensing reflectance @ 620 nm"},
  {"R_665","Remote-sensing reflectance @ 665 nm"},
  {"R_681","Remote-sensing reflectance @ 681 nm"},
  {"R_710","Remote-sensing reflectance @ 710 nm"},
  {"R_709","Remote-sensing reflectance @ 709 nm"},
  {"R_753","Remote-sensing reflectance @ 753 nm"},
  {"R_754","Remote-sensing reflectance @ 754 nm"},
  {"R_482","Remote-sensing reflectance @ 482 nm"},
  {"R_655","Remote-sensing reflectance @ 655 nm"},
  {"nFLH","normalised Fluorescence Line Height"},
  {"Secchi","Secchi depth"},
  {"Zenith2D","Solar zenith"},
  {"SWR_bot_abs","SWR bottom abs. (PAR)"},
  {"Oxygen_sedflux","Sediment-water oxygen flux"},
  {"DIC_sedflux","Sediment-water DIC flux"},
  {"NH4_sedflux","Sediment-water NH4 flux"},
  {"NO3_sedflux","Sediment-water NO3 flux"},
  {"DIP_sedflux","Sediment-water DIP flux"},
  {"Hue","Hue angle"},
  {"Moonlight","Moonlight (PAR-integrated)"},
  {"Lunar_zenith","Moon zenith"},
  {"Lunar_phase","Moon phase"},
  {"Moon_fulldisk","Moonlight brightness (PAR-integrated)"},
};
const int NUM_ECO_VARS_2D = ((int)(sizeof(ECONAME2D)/(2*sizeof(char*))));

/* Tracerstats */
const eco_def_trstat_t eco_def_trstat[] = {
  /* name */
  {"month_Gnet", "Monthly net calcification rate", "Gnet", "run_mean(Gnet)", "30 days"},
  {"month_EpiPAR_sg", "Monthly PAR light above seagrass", "EpiPAR_sg", "run_mean(EpiPAR_sg)", "30 days"},
  {"week_erdepflux_total", "Weekly net deposition rate", "erdepflux_total", "run_mean(erdepflux_total)", "7 days"},
  {"weekly_Coral_IN_up", "Weekly coral inorganic N uptake", "Coral_IN_up", "run_mean(Coral_IN_up)", "7 days"},
  {"daily_Nfix", "daily N2 fixation", "Nfix", "run_mean(Nfix)", "1 day"},
  {"daily_Den_fl", "daily Denitrification flux", "Den_fl", "run_mean(Den_fl)", "1 day"},
  /* Exposure time must be called <tracer_name>_time, like below */
  {"omega_ar_expose", "Aragonite saturation exposure", "omega_ar", 
                                 "exposure(omega_ar:-3:omega_ar_expose_time)", "30 days"},
  {"NULL", "", ""}
};

static void *private_data_copy_eco(void *src);

static int e_ntr;
static int tr_map[MAXNUMVARS];
static int sed_map[MAXNUMVARS];
static int e_nepi;
static int epi_map[MAXNUMVARS];

extern int i_tracername_exists(void* model, char*name);
extern int i_tracername_exists_2d(void* model, char*name);
extern int i_tracername_exists_sed(void* model, char*name);

static int get_eco_flag_num(char *flags);
static char *eco_flag_str(int flag);

/*-------------------------------------------------------------------*/
/* Ecology specific interface routines
void einterface_ecologyinit(void *model, void *_e)
static void einterface_tracermap(void* model, ecology *e, int ntr)
static void einterface_epimap(void* model, ecology *e, int ntr)
int einterface_get_eco_flag(void* model, char* name)
double einterface_gettracersvel(void* model, char* name)
quitfntype einterface_getquitfn(void)
void einterface_log_error(ecology *e, void *model, int b)
void einterface_check_default_tracers(void)
*/

/*-------------------------------------------------------------------*/
/* Generic interface routines                                        *
/*-------------------------------------------------------------------*/
extern int ginterface_get_num_rsr_tracers(void* model);
extern void ginterface_get_rsr_tracers(void* model, int *rtns);
extern int ginterface_getntracers(void* model);
extern int ginterface_gettracerdiagnflag(void* model, char* name);
extern int ginterface_gettracerparticflag(void* model, char* name);
extern char* ginterface_gettracername(void* model, int i);
extern char* ginterface_get2Dtracername(void* model, int i);
extern int ginterface_tracername_exists(void* model, char*name);
extern int ginterface_tracername_exists_epi(void* model, char*name);
extern void ginterface_get_ij(void* model, int col, int *ij);
extern int ginterface_getnepis(void* model);
extern char* ginterface_getepiname(void* model, int i);
extern int ginterface_getepidiagnflag(void* model, char* name);
extern tracer_info_t* ginterface_getepiinfo(void* model, int* n);
extern double ginterface_getmodeltime(void* model);
extern int ginterface_getnumberwclayers(void* model);
extern int ginterface_getnumbersedlayers(void* model);
extern int ginterface_getnumbercolumns_e(void* model);
extern int ginterface_get_max_numbercolumns(void* model);
extern int ginterface_getwctopk(void *model, int b);
extern int ginterface_getwcbotk(void *model, int b);
extern double ginterface_getcellz(void *model, int b, int k);
extern int ginterface_isboundarycolumn(void *model, int b);
extern int ginterface_getsedtopk(void *model, int b);
extern int ginterface_getsedbotk(void *model, int b);
extern double *ginterface_getwccellthicknesses(void *model, int b);
extern double *ginterface_getsedcellthicknesses(void *model, int b);
extern double ginterface_botstress(void *model, int b);
extern double ginterface_get_windspeed(void *model, int b);
extern double ginterface_getlighttop(void *model, int b);
extern double *ginterface_getporosity(void *model, int b);
extern double ginterface_geterosionrate(void *model, int b);
extern double ginterface_getustrcw(void *model, int b);
extern double** ginterface_getwctracers(void* model, int b);
extern double** ginterface_getsedtracers(void* model, int b);
extern double **ginterface_getepivars(void *model, int b);
extern double ginterface_cellarea(void* hmodel, int b);
extern int ginterface_getverbosity(void *model);
extern char *ginterface_gettimeunits(void *model);
extern int ginterface_transport_mode(void);
extern char * ginterface_get_output_path(void);
extern double ginterface_calc_zenith(void *model, double t, int b);
extern int ginterface_get_win_num(void *model);

/*-------------------------------------------------------------------*/
/* Re-directed interface routines. These should be replaced with     */
/* direct calls to the generic interface routines, bypassing         */
/* wrappers where possible.                                          */
/*-------------------------------------------------------------------*/

int einterface_get_num_rsr_tracers(void* model) {
  return ginterface_get_num_rsr_tracers(model);
}
void einterface_get_rsr_tracers(void* model, int *rtns) {
  ginterface_get_rsr_tracers(model, rtns);
}
int einterface_getntracers(void* model) {
  return ginterface_getntracers(model);
}
int einterface_gettracerdiagnflag(void* model, char* name) {
  return ginterface_gettracerdiagnflag(model, name);
}
int einterface_gettracerparticflag(void* model, char* name) {
  return ginterface_gettracerparticflag(model, name);
}
char* einterface_gettracername(void* model, int i) {
  return ginterface_gettracername(model, i);
}
char* einterface_get2Dtracername(void* model, int i) {
  return ginterface_get2Dtracername(model, i);
}
int einterface_tracername_exists(void* model, char*name) {
  return ginterface_tracername_exists(model, name);
}
int einterface_tracername_exists_epi(void* model, char*name) {
  return ginterface_tracername_exists_epi(model, name);
}
void einterface_get_ij(void* model, int col, int *ij) {
  ginterface_get_ij(model, col, ij);
}
int einterface_getnepis(void* model) {
  return ginterface_getnepis(model);
}
char* einterface_getepiname(void* model, int i) {
  return ginterface_getepiname(model, i);
}
int einterface_getepidiagnflag(void* model, char* name) {
  return ginterface_getepidiagnflag(model, name);
}
tracer_info_t* einterface_getepiinfo(void* model, int* n) {
  return ginterface_getepiinfo(model, n);
}
double einterface_getmodeltime(void* model) {
  return ginterface_getmodeltime(model);
}
int  einterface_getnumberwclayers(void* model) {
  return ginterface_getnumberwclayers(model);
}
int  einterface_getnumbersedlayers(void* model) {
  return ginterface_getnumbersedlayers(model);
}
int  einterface_getnumbercolumns(void* model) {
  return ginterface_getnumbercolumns_e(model);
}
int  einterface_get_max_numbercolumns(void* model) {
  return ginterface_get_max_numbercolumns(model);
}
int einterface_getwctopk(void *model, int b) {
  return ginterface_getwctopk(model, b);
}
int einterface_getwcbotk(void *model, int b) {
  return ginterface_getwcbotk(model, b);
}
double einterface_getcellz(void *model, int b, int k) {
  return ginterface_getcellz(model, b, k);
}
int einterface_isboundarycolumn(void *model, int b) {
  return ginterface_isboundarycolumn(model, b);
}
int einterface_getsedtopk(void *model, int b) {
  return ginterface_getsedtopk(model, b);
}
int einterface_getsedbotk(void *model, int b) {
  return ginterface_getsedbotk(model, b);
}
double *einterface_getwccellthicknesses(void *model, int b) {
  return ginterface_getwccellthicknesses(model, b);
}
double *einterface_getsedcellthicknesses(void *model, int b) {
  return ginterface_getsedcellthicknesses(model, b);
}
double einterface_botstress(void *model, int b) {
  return ginterface_botstress(model, b);
}
double einterface_get_windspeed(void *model, int b) {
  return ginterface_get_windspeed(model, b);
}
double einterface_getlighttop(void *model, int b) {
  return ginterface_getlighttop(model, b);
}
double *einterface_getporosity(void *model, int b) {
  return ginterface_getporosity(model, b);
}
double einterface_geterosionrate(void *model, int b) {
  return ginterface_geterosionrate(model, b);
}
double einterface_getustrcw(void *model, int b) {
  return ginterface_getustrcw(model, b);
}
/* Issues with returning pointers; duplacte einterface_ routines for now.
double**  einterface_getwctracers(void* model, int b) {
  return ginterface_getwctracers(model, b);
}
double** einterface_getsedtracers(void* model, int b) {
  return ginterface_getsedtracers(model, b);
}
double **einterface_getepivars(void *model, int b) {
  return ginterface_getepivars(model, b);
}
*/
double einterface_cellarea(void* hmodel, int b) {
  return ginterface_cellarea(hmodel, b);
}
int einterface_getverbosity(void *model) {
  return ginterface_getverbosity(model);
}
char *einterface_gettimeunits(void *model) {
  return ginterface_gettimeunits(model);
}
int einterface_transport_mode(void) {
  return ginterface_transport_mode();
}
char * einterface_get_output_path(void) {
  ginterface_get_output_path();
}
int einterface_get_win_num(void *model) {
  return ginterface_get_win_num(model);
}
double einterface_calc_zenith(void *model, double t, int b) {
  return ginterface_calc_zenith(model, t, b);
}

/*-------------------------------------------------------------------*/


double**  einterface_getwctracers(void* model, int b)
{
  geometry_t *window = (geometry_t *)model;
  window_t *windat = window->windat;
  int nz = window->nz;
  int ntr = e_ntr;
  double **wctr = (double **)calloc(ntr * nz, sizeof(double*));
  int cs2 = window->wincon->s2[b+1];
  int c2 = window->m2d[cs2];
  int cc = window->c2cc[c2];
  int cs = window->nsur_t[cc];
  int cb = window->bot_t[cc];
  int c, k = window->s2k[cb], n;

  /*
   * Loop from the bottom cell up
   */
  for (c = cb; c != cs; c = window->zp1[c]) {
    for (n = 0; n < ntr; ++n) {
      wctr[k * ntr + n] = &windat->tr_wc[tr_map[n]][c];
    }
    k++;
  }
  // The surface cell
  for (n = 0; n < ntr; ++n) {
    wctr[k * ntr + n] = &windat->tr_wc[tr_map[n]][c];
  }

  return wctr;
}

double** einterface_getsedtracers(void* model, int b)
{
  geometry_t *window = (geometry_t *)model;
  window_t *windat = window->windat;
  int nz = window->sednz;
  int ntr = e_ntr;
  double **sedtr = (double **)calloc(ntr * nz, sizeof(double*));
  int c = window->wincon->s2[b + 1];
  int c2=window->m2d[c];
  int k, n;

  for (k = 0; k < nz; ++k)
    for (n = 0; n < ntr; ++n) {
      sedtr[k * ntr + n] = &windat->tr_sed[sed_map[n]][k][c2];
    }
  return sedtr;
}


double **einterface_getepivars(void *model, int b)
{
  geometry_t *window = (geometry_t *)model;
  window_t *windat = window->windat;
  int nepi = e_nepi;
  double **epivar = (double **)calloc(nepi, sizeof(double*));
  int c = window->wincon->s2[b + 1];
  int c2 = window->m2d[c];
  int n;

  for (n = 0; n < nepi; ++n) {
    epivar[n] = &windat->tr_wcS[epi_map[n]][c2];
  }
  return epivar;
}



/*-------------------------------------------------------------------*/

static void einterface_tracermap(void* model, ecology *e, int ntr)
{
  geometry_t *window = (geometry_t *)model;
  win_priv_t *wincon = window->wincon;
  char trname[MAXSTRLEN];
  int tn,n,m;

  n = m = 0;
  for(tn=0; tn<ntr; tn++) {
    strcpy(trname,ecology_gettracername(e, tn));
    /* Water column */
    if((tr_map[n] = tracer_find_index(trname, wincon->ntr,
              wincon->trinfo_3d)) >= 0) {
      n++;
    }
    else
      hd_quit("ecology interface: Can't find type WATER tracer '%s' in parameter file.\n", trname);
    /* Sediments */
    if((sed_map[m] = tracer_find_index(trname, wincon->nsed,
              wincon->trinfo_sed)) >= 0) {
      m++;
    }
    else
      hd_quit("ecology interface: Can't find type SEDIMENT tracer '%s' in parameter file.\n", trname);
  }
}

static void einterface_epimap(void* model, ecology *e, int ntr)
{
  geometry_t *window = (geometry_t *)model;
  win_priv_t *wincon = window->wincon;
  char trname[MAXSTRLEN];
  int tn,n;

  n = 0;
  for(tn=0; tn<ntr; tn++) {
    strcpy(trname,ecology_getepiname(e, tn));
    if((epi_map[n] = tracer_find_index(trname, wincon->ntrS,
               wincon->trinfo_2d)) >= 0) {
      n++;
    }
    else
      hd_warn("ecology : Can't find type BENTHIC tracer '%s' in parameter file.\n", trname);
  }
}

int einterface_get_eco_flag(void* model, char* name)
{

  
  tracer_info_t *tr = i_get_tracer(model, name);
  trinfo_priv_eco_t *data = tr->private_data[TR_PRV_DATA_ECO];
  int flag = ECO_NONE;
    
  if (data)
    flag = data->flag;

  return flag;
}


double einterface_gettracersvel(void* model, char* name)
{
  return(sinterface_get_svel(model,name));

}


quitfntype einterface_getquitfn(void)
{
  return (quitfntype) hd_quit_and_dump;
}


void einterface_ecologyinit(void *model, void *_e)
{
  geometry_t *window = (geometry_t *)model;
  
  
  ecology* e = (ecology*) _e;

  /*
   * move to ecology_step
   */
  /*
    if (wincon->ecodt < windat->dt) {
    hd_warn
    ("Ecomodel timestep = %f < physical timestep = %f;\nsetting ecomodel timestep = physical timestep.\n",
    wincon->ecodt, windat->dt);
    wincon->ecodt = windat->dt;
    }
    
    wincon->eco_timestep_ratio = wincon->ecodt / windat->dt;
    if (wincon->ecodt != windat->dt * wincon->eco_timestep_ratio) {
    wincon->ecodt = windat->dt * wincon->eco_timestep_ratio;
    hd_warn("Set ecomodel timestep to %f seconds.\n", wincon->ecodt);
    }
  */

  /* Ecology requires water column tracers to be both of type WATER */
  /* and SEDIMENT and requires a 1:1 corresponance between the */
  /* indicies for the water column and sediment components. This */
  /* requires an index map to be established. */
  /* Get the number of water column tracers used by ecology. */
  e_ntr = ecology_getntracers(e);
  /* Get the water column tracers required by ecology and the tracer */
  /* index map. */
  einterface_tracermap(window, e, e_ntr);

  /* Get the number of water column tracers used by ecology */
  e_nepi = ecology_getnepis(e);
  /* Get the water column tracers required by ecology and the tracer */
  /* index map. */
  einterface_epimap(window, e, e_nepi);
}


/*
 * Logs the ecology error and optionally quits if error function is
 * set to none
 */
void einterface_log_error(ecology *e, void *model, int b)
{
  geometry_t *window = (geometry_t*)model;
  window_t   *windat = window->windat;
  win_priv_t *wincon = window->wincon;
  int c  = wincon->s2[b+1];
  int c2 = window->m2d[c];
  int eco_nstep = eco_get_nstep(e);

  hd_warn("Ecology library error @ (%d %d) : %5.2f days\n", 
	  window->s2i[c2], window->s2j[c2], windat->days);
  windat->ecoerr[c2] = 100.0 * 
    (eco_nstep * 0.01 * windat->ecoerr[c2] + 1.0) / (eco_nstep+1);
  
  if (wincon->gint_errfcn == LWARN)
    hd_warn(wincon->gint_error[b]);
  if (wincon->gint_errfcn == LFATAL)
    hd_quit(wincon->gint_error[b]);
}


/*
 * Sanity check (called from ecology_pre_build) to check for
 * duplicates in the tracer default list. The stringtable_add will
 * error out when the same name is added twice
 */
void einterface_check_default_tracers(void)
{
  stringtable *tbl;
  int i;

  tbl = stringtable_create("Tracer-defaults-ECONAME3D");
  for (i = 0; i < NUM_ECO_VARS_3D; i++)
    stringtable_add(tbl, (char*)ECONAME3D[i][0], i);
  stringtable_destroy(tbl);

  tbl = stringtable_create("Tracer-defaults-ECONAME2D");
  for (i = 0; i < NUM_ECO_VARS_2D; i++)
    stringtable_add(tbl, (char*)ECONAME2D[i][0], i);
  stringtable_destroy(tbl);
}


/*-------------------------------------------------------------------*/
/* Ecology step                                                      */
/*-------------------------------------------------------------------*/
void eco_step(geometry_t *window)
{
  win_priv_t *wincon = window->wincon;
  window_t   *windat = window->windat;

  /*-----------------------------------------------------------------*/
  /* Do the ecology if required */
#if defined(HAVE_ECOLOGY_MODULE)
  if (wincon->do_eco) {
    /*
     * Redo eco_timestep_ratio calculation as tratio may have altered
     * windat->dt in hd_step[_trans]
     */
    if (wincon->ecodt < windat->dt) {
      hd_warn
	("Ecomodel timestep = %f < physical timestep = %f;\nsetting ecomodel timestep = physical timestep.\n",
	 wincon->ecodt, windat->dt);
      wincon->ecodt = windat->dt;
    }

    wincon->eco_timestep_ratio = wincon->ecodt / windat->dt;
    /* Trike introduces minor rounding errors due to timeunit conversions */
    if (fabs(wincon->ecodt - (windat->dt * wincon->eco_timestep_ratio)) > 1e-3)
      hd_quit("Inconsistency in setting ecology time step as a multiple of physical timestep: eco_dt=%f, phys_dt=%f\n", wincon->ecodt, windat->dt);

    /* Run ecology if its the right time */
    if (!(windat->nstep % wincon->eco_timestep_ratio)) {
      TIMING_SET;
      ecology_step(wincon->e, wincon->ecodt);
      TIMING_DUMP_WIN(3, "   eco_step", window->wn);
    }
  }
#endif

}

/* END eco_step()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets the ecology tracers standard default attributes              */
/*-------------------------------------------------------------------*/
void eco_set_tracer_defaults(tracer_info_t *tracer, char *trname,
			     char *defname, ecology *e)
{
  int i;
  struct {
    char *name;
    char *description;
    void (*init) (tracer_info_t *tracer, char *trname, ecology *e);
  } def_list[] = {
    {"standard","Standard ecology values",  eco_defaults_std},
    {"estuary", "Estuarine ecology values", eco_defaults_est},
    {NULL, NULL, NULL}
  };

  void (*init) (tracer_info_t *tracer, char *trname, ecology *e)= NULL;

  if (defname != NULL) {
    /* Search for the eco defaults scheme name in the above table */
    for (i = 0; (init == NULL) && def_list[i].name; ++i) {
      if (strcasecmp(defname, def_list[i].name) == 0) {
        init = def_list[i].init;
      }
    }
  }
  if (init != NULL) {
    init(tracer, trname, e);
  }
}


/*-------------------------------------------------------------------*/
/* Initialised ecology private data                                  */
/*-------------------------------------------------------------------*/
trinfo_priv_eco_t *get_private_data_eco(tracer_info_t *tr)
{
  trinfo_priv_eco_t *data = tr->private_data[TR_PRV_DATA_ECO];

  /* Allocate memory, if needed */
  if (data == NULL) {
    data = (trinfo_priv_eco_t *)malloc(sizeof(trinfo_priv_eco_t));
    /*
     * Useful for debugging in case the wrong pointer happens 
     * to be passed in 
     */
    data->type = PTR_BGC;
    
    /* Set data and copy function */
    tr->private_data[TR_PRV_DATA_ECO]      = data;
    tr->private_data_copy[TR_PRV_DATA_ECO] = private_data_copy_eco;
  
    sprintf(data->name, "%c", '\0');
    data->obc  = NOGRAD;
    data->flag = ECO_NONE;
  }
  
  return(data);
}


/*-------------------------------------------------------------------*/
/* Copies ecology private data                                       */
/*-------------------------------------------------------------------*/
static void *private_data_copy_eco(void *src)
{
  void *dest = NULL;

  if (src) {
    dest = (trinfo_priv_eco_t *)malloc(sizeof(trinfo_priv_eco_t));
    memcpy(dest, src, sizeof(trinfo_priv_eco_t));
  }
  return(dest);
}

/* END private_data_copy_eco()                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Initialised ecology private data                                  */
/*-------------------------------------------------------------------*/
void init_eco_private_data(trinfo_priv_eco_t *data)
{
  data->flag = ECO_NONE;
}

/* END init_eco_private_data()                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Searches through the tracer list and returns the tracer type      */
/*-------------------------------------------------------------------*/
int get_eco_var_type(char *trname, char *eco_defs)
{
  tracer_info_t tr;
  memset(&tr, 0, sizeof(tracer_info_t));
  eco_set_tracer_defaults(&tr, trname, eco_defs, NULL);
  if (!tr.type)
    hd_quit("ecology:get_eco_var_type: Tracer '%s' is not listed in '%s'\n",
	    trname, eco_defs);

  return(tr.type);
}


/*-------------------------------------------------------------------*/
/* Checks if a trcer is a valid ecology tracer                       */
/*-------------------------------------------------------------------*/
int is_eco_var(char *trname)
{
  int i;

  for (i = 0; i < NUM_ECO_VARS_3D; i++)
    if (strcmp(trname, ECONAME3D[i][0]) == 0) {
      return(1);
    }
  for (i = 0; i < NUM_ECO_VARS_2D; i++)
    if (strcmp(trname, ECONAME2D[i][0]) == 0) {
      return(1);
    }
  // Return false
  return (0);
}

/* END is_eco_tracer()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads the ecology related attributes                              */
/*-------------------------------------------------------------------*/
void set_eco_atts(tracer_info_t *tr, FILE *fp, char *keyname)
{

}

/* END set_eco_atts()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Perform any custom adjustments to the ecology tracer info         */
/*-------------------------------------------------------------------*/
void eco_tracer_custom(master_t *master)
{
  
  
  
  
  

}

/* END eco_tracer_custom()                                           */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Reads the ecology related attributes                              */
/*-------------------------------------------------------------------*/
void eco_read_tr_atts(tracer_info_t *tr, FILE *fp, char *keyname)
{
  trinfo_priv_eco_t *data;
  char buf[MAXSTRLEN], key[MAXSTRLEN];

  /* Check for a ecology variable */
  // if ( !is_eco_var(tr->name) ) return;

  // Allocate memory, if needed
  data = get_private_data_eco(tr);

  sprintf(key, "%s.eco_flag", keyname);
  if (prm_read_char(fp, key, buf))
    data->flag = get_eco_flag_num(buf);

}


/*-------------------------------------------------------------------*/
/* Writes the attributes for ecology tracers                         */
/*-------------------------------------------------------------------*/
void eco_write_tr_atts(tracer_info_t *tr, FILE *fp, int n)
{
  trinfo_priv_eco_t *data = tr->private_data[TR_PRV_DATA_ECO];

  if (data != NULL && data->type == PTR_BGC) {
    char *flags = eco_flag_str(data->flag);
    fprintf(fp, "TRACER%1.1d.eco_flag        %s\n", n, flags);
    free(flags);
  }
}
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Counts the number of valid ecology classes                        */
/*-------------------------------------------------------------------*/
int count_eco_3d(char *eco_vars)
{
  char *vars[ECO_MAXNUMARGS], buf[MAXSTRLEN];
  int m, n, i, nvars = 0;

  if (eco_vars == NULL || strlen(eco_vars) == 0)
    return(0);

  strcpy(buf, eco_vars);
  m = parseline(buf, vars, ECO_MAXNUMARGS);
  for (n = 0; n < m; n++) {
    for (i = 0; i < NUM_ECO_VARS_3D; i++)
      if (strcmp(vars[n], ECONAME3D[i][0]) == 0) {
	nvars++;
      }
  }
  if (nvars == 0) {
    hd_quit("No valid ecology 3D tracers found.\n");
  }
  return(nvars);
}

int count_eco_2d(char *eco_vars)
{
  char *vars[ECO_MAXNUMARGS], buf[MAXSTRLEN];
  int m, n, i, nvars = 0;

  if (eco_vars == NULL || strlen(eco_vars) == 0)
    return(0);

  strcpy(buf, eco_vars);
  m = parseline(buf, vars, ECO_MAXNUMARGS);
  for (n = 0; n < m; n++) {
    for (i = 0; i < NUM_ECO_VARS_2D; i++)
      if (strcmp(vars[n], ECONAME2D[i][0]) == 0) {
	nvars++;
      }
  }
  if (nvars == 0) {
    //  hd_quit("No valid ecology 2D tracers found.\n");
  }
  return(nvars);
}

/* END count_eco()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to initialise the 2D tracers in the master                */
/*-------------------------------------------------------------------*/
static int eco_set_autotracer(FILE *fp, 
			      int   do_eco, 
			      char *eco_vars,  /* list of vars */ 
			      char *eco_defs,  /* name of class */
			      void *e,         /* pre_build eco struct */
			      tracer_info_t *trinfo, 
			      int ntr, int tn,
			      int         necoclass,
			      const char *ecoclass[][2],
			      int trinfo_type)
{
  char *vars[ECO_MAXNUMARGS], buf[MAXSTRLEN];
  int m, n, i;

  /* Set up the auto-tracers */
  if (do_eco && strlen(eco_vars)) {
    strcpy(buf, eco_vars);
    m = parseline(buf, vars, ECO_MAXNUMARGS);
    for (n = 0; n < m; n++) {
      for (i = 0; i < necoclass; i++) {
	if (strcmp(vars[n], ecoclass[i][0]) == 0) {
	  if ((tracer_find_index(vars[n], ntr, trinfo)) < 0) {
	    /* 
	     * Set the default value 
	     */
	    eco_set_tracer_defaults(&trinfo[tn], vars[n], eco_defs,
				    (ecology*)e);
	    /*
	     * Overwrite any att from prm-file
	     */
	    sed_read_tr_atts(&trinfo[tn], fp, vars[n]); // includes data
	    eco_read_tr_atts(&trinfo[tn], fp, vars[n]);
	    trinfo[tn].m = -1; // does not exist in the prm-file
	    tracer_re_read(&trinfo[tn], fp, trinfo_type);
	    /*
	     * Fill in the tracer number
	     */
	    trinfo[tn].n = tn;
	    tn++;
	  }
	}
      }
    }
  }
  return(tn);
}


/*-------------------------------------------------------------------*/
/* Routine to initialise the 2D tracers in the master                */
/*-------------------------------------------------------------------*/
int ecology_autotracer_3d(FILE *fp, int do_eco, char *eco_vars, char *eco_defs,
			  void *e, tracer_info_t *trinfo, int ntr, int tn)
{
  /* ECO 3D */
  tn = eco_set_autotracer(fp,  do_eco, eco_vars, eco_defs, e, trinfo, ntr, tn,
			  NUM_ECO_VARS_3D, ECONAME3D, WATER);
  return(tn);
}

/* END ecology_autotracer_3d()                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to initialise the 2D tracers in the master                */
/*-------------------------------------------------------------------*/
int ecology_autotracer_2d(FILE *fp, int do_eco, char *eco_vars, char *eco_defs,
			  void *e, tracer_info_t *trinfo, int ntr, int tn)
{
  /* ECO 2D */
  tn = eco_set_autotracer(fp,  do_eco, eco_vars, eco_defs, e, trinfo, ntr, tn,
			  NUM_ECO_VARS_2D, ECONAME2D, INTER);
  
  return(tn);
}

/* END ecolocy_autotracer_2d()                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to initialise the 2D tracers in the master                */
/*-------------------------------------------------------------------*/
int ecology_autotracer_sed(FILE *fp, int do_eco, char *eco_vars, char *eco_defs,
			   void *e, tracer_info_t *trinfo, int ntr, int tn)
{
  /* ECO SED */
  tn = eco_set_autotracer(fp,  do_eco, eco_vars, eco_defs, e, trinfo, ntr, tn,
			  NUM_ECO_VARS_3D, ECONAME3D, SEDIM);
  return(tn);
}

/* END ecology_autotracer_sed()                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Set the open boundary condition as specified in the ecology       */
/* private data for sediment autotracers.                            */
/*-------------------------------------------------------------------*/
int eco_get_obc(tracer_info_t *tr)
{
  trinfo_priv_eco_t *data = tr->private_data[TR_PRV_DATA_ECO];
  int ret = NOGRAD;

  if (data != NULL) {
    if (data->type == PTR_BGC && data->obc > 0)
      ret = data->obc;
  }
  return(ret);
}

/* END eco_get_obc()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Writes the sediment autotracers to the setup file                 */
/*-------------------------------------------------------------------*/
int ecology_autotracer_write(master_t *master, FILE *op, int tn)
{
  parameters_t *params = master->params;
  
  int i, n, m, trn;

  if (!params->do_eco || strlen(params->eco_vars) == 0) return(0);

  fprintf(op, "#\n");
  fprintf(op, "# Ecology tracers - 3D\n");
  fprintf(op, "#\n");
  /* 3D water column and sediment autotracers                        */
  n = tn;
  for (m = 0; m < params->atr; m++) {
    tracer_info_t *tracer = &master->trinfo_3d[m];

    trn = -1;
    for (i = 0; i < NUM_ECO_VARS_3D; i++) {
      if (strcmp(tracer->name, ECONAME3D[i][0]) == 0) {
	trn = i;
	break;
      }
    }
    if (trn < 0) continue;

    /* Call the 3d helper function */
    tracer_write_3d(master, op, tracer, n);

    n++;
  }

  fprintf(op, "#\n");
  fprintf(op, "# Ecology tracers - BENTHIC\n");
  fprintf(op, "#\n");
  /* 2D sediment autotracers                                         */
  for (m = 0; m < master->ntrS; m++) {
    tracer_info_t *tracer = &master->trinfo_2d[m];

    trn = -1;
    for (i = 0; i < NUM_ECO_VARS_2D; i++)
      if (strcmp(tracer->name, ECONAME2D[i][0]) == 0) {
	trn = i;
	break;
      }
    if (trn < 0) continue;

    /* Call the 2d helper function */
    tracer_write_2d(master, op, tracer, n);

    n++;
  }
  return(n);
}

/* END ecology_autotracer_write()                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Eco flag number to string                                         */
/*-------------------------------------------------------------------*/
static char *eco_flag_str(int flag)
{
  char buf[MAXSTRLEN];
  if (flag == 0) {
    sprintf(buf, "NONE");
  } else {
    if (flag & ECO_NORESET)
      sprintf(buf, "NORESET ");
  }

  return(strdup(buf));
}
/* END eco_flag_num                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Eco flag string to number                                         */
/*-------------------------------------------------------------------*/
static int eco_flag_num(char *flag)
{
  if ( (flag == NULL) || (strcmp(flag, "0") == 0) )
    return(ECO_NONE);
  if (strcmp(flag, "NORESET") == 0)
    return(ECO_NORESET);
}
/* END eco_flag_num                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Returns number for all the matching flags                         */
/*-------------------------------------------------------------------*/
static int get_eco_flag_num(char *flags)
{
  int eflag = ECO_NONE;
  char *tok;
  tok = strtok(flags, " ");
  eflag |= eco_flag_num(tok);
  while (tok != NULL) {
    tok = strtok(NULL, " ");
    eflag |= eco_flag_num(tok);
  }

  return (eflag);
}

/* END get_eco_flag_num()                                            */
/*-------------------------------------------------------------------*/

/*-------------*/
/* Tracerstats */
/*-------------*/
/* 
 * Returns the number of ecology tracerstats asscociated with this tracer 
 */
int ecology_count_tracerstats(tracer_info_t *trinfo)
{
  int n = 0;
  const eco_def_trstat_t *trs = &eco_def_trstat[n];
  while (strcmp(trs->name, "NULL") != 0) {
    if (strcmp(trinfo->name, trs->trname) == 0) {
      if (strncmp(trs->trstat, "exposure", 8) == 0)
	return(2); // allow for ex_time
      return(1);
    }
    trs = &eco_def_trstat[++n];
  }
  return(0);
}

/* write out tracerstats */
int ecology_write_tracerstat(master_t *master, FILE *fp, 
			     tracer_info_t *trinfo, int tn)
{
  int n = 0;
  const eco_def_trstat_t *trs = &eco_def_trstat[n];
  char buf[MAXSTRLEN];
  while (strcmp(trs->name, "NULL") != 0) {
    if (strcmp(trinfo->name, trs->trname) == 0) {
      fprintf(fp, "TRACER%1.1d.name            %s\n", tn, trs->name);
      fprintf(fp, "TRACER%1.1d.long_name       %s\n", tn, trs->long_name);
      /* Hack for exposure as it's not supported by sediments just yet */
      if (strncmp(trs->trstat, "exposure", 8)) {
	// fprintf(fp, "TRACER%1.1d.type           %s\n", tn, trtypename(trinfo->type, buf));
	if (trinfo->type & INTER)
	  fprintf(fp, "TRACER%1.1d.type            BENTHIC DIAGNOSTIC\n", tn);
	else
	  fprintf(fp, "TRACER%1.1d.type            WATER SEDIMENT DIAGNOSTIC\n", tn);
      }
      fprintf(fp, "TRACER%1.1d.fill_value      0.0\n", tn);
      fprintf(fp, "TRACER%1.1d.units           %s\n",  tn, trinfo->units);
      fprintf(fp, "TRACER%1.1d.advect          0\n", tn);
      fprintf(fp, "TRACER%1.1d.diffuse         0\n", tn);
      fprintf(fp, "TRACER%1.1d.diagn           0\n", tn);
      fprintf(fp, "TRACER%1.1d.tracerstat      %s\n", tn, trs->trstat);
      fprintf(fp, "TRACER%1.1d.dt              %s\n", tn, trs->trdt);
      fprintf(fp, "\n");
      tn++;
      /* Need to explicitly write out the extra exposure time tracer */
      if (strncmp(trs->trstat, "exposure", 8) == 0) {
	sprintf(buf, "%s_time", trs->name);
	fprintf(fp, "TRACER%1.1d.name            %s\n", tn, buf);
	sprintf(buf, "%s time", trs->long_name);
	fprintf(fp, "TRACER%1.1d.long_name       %s\n", tn, buf);
	fprintf(fp, "TRACER%1.1d.fill_value      0.0\n", tn);
	fprintf(fp, "TRACER%1.1d.units           day\n", tn);
	fprintf(fp, "TRACER%1.1d.advect          0\n", tn);
	fprintf(fp, "TRACER%1.1d.diffuse         0\n", tn);
	fprintf(fp, "TRACER%1.1d.diagn           0\n", tn);
	fprintf(fp, "\n");
	tn++;
      }
    }
    trs = &eco_def_trstat[++n];
  }
  return(tn);
}

#endif /* HAVE_ECOLOGY */

// EOF
