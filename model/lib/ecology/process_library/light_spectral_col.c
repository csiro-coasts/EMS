/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  File: model/lib/ecology/process_library/light_spectral_col.c
 *  
 *  Description:
 *
 *  This process initialises the spectrally-resolved optical model, uses mass-specific optical properties to calculate the 
 *  IOPs, and then uses the entire column and benthos to calculate the AOPs. AOPs are calculated using a ray-tracing technique, 
 *  or using the links to the commercial radiative transfer model HYDROLIGHT. This second option requires a HYDROLIGHT licence.
 *  
 *
 *  Baird, M. E., K. Wild-Allen, J. Parslow, M. Mongin, B. Robson, J. Skerratt, F. Rizwi, M. Soja-Wozniak, E. Jones, 
 *      M. Herzfeld, N. Margvelashvili, J. Andrewartha, C. Langlais, M. Adams, N. Cherukuru, S. Hadley, P. Ralph, 
 *      T. Schroeder, A. Steven, U. Rosebrock, L Laiolo, M. Gustafsson, and D. Harrison (2020). CSIRO Environmental 
 *      Modelling Suite (EMS): Scientific description of the optical and biogeochemical models (vB3p0). 
 *      Geoscientific Model Development.13:4503-4553.
 *  
 *  Column based spectrally-resolved optical model. There is no "calc" routine. The column process is undertaken 
 *      after the postcalc of all wc, epi and sed processes.
 *
 *  Notes: 1. Generates four files, optical_setup.nc, optical_setup_matlab_check.m and optical_setup_python_check.py 
 *              (both for diagnostics) and optical_column_out.nc (for detailed optical outputs in one column).
 *         2. Runs the original EMS light model (Beer-Lambert exponential attenuation) every cell.
 *         3. Option to radiative transfer Hydrolight, by specification in init (#define INCLUDE_HYDROLIGHT 1)
 *         4. In column process, no need to add offset for epi variables in y_epi.
 *         5. EMS has both a lunar and solar zenith and azimuth. Hydrolight has only one - solar if it calculates itself, or 
 *              the choice of solar, lunar or starlight if sent to be EMS.
 *  
 *  Changes from B3p1 light model:
 *
 *         1. Use of algal pigments for absorption cross-section of MPB.
 *         2. Specification of reflectance, absorptance (photosynthetic and non-photosynthetic) and transmission
 *              for benthic plants - reflection reduces light through canopy.
 *         3. CDOM absorption only uses CDOM tracers.
 *         4. Reflectance calculated for sensor bands.
 *         5. Spectrally-resolved bioluminescence and fluorescence emission spectra. 
 *         6. Wavelength dependence of air-sea water interface refraction.
 *         7. SWR fraction into reddest band is the remainder of the light (this reduces moonlight significantly, but doesn't 
 *            change PAR).
 *         8. All AOPs done on the same time step as IOPs.
 *         9. IOPs, AOPs and photosynthetic absorption calculated at postcalc time, then absorption applied during next timestep.
 *         10. Added sediment surface reflectance from free symbionts using pigments for absorption.
 *         11. For AOP calculations, default zenith / azimuth to Sunlight, then Moonlight, then Starlight (zenith = 0).
 *         12. Opaqueness limited to PAR to avoid zeros in pigment absorption.
 *         13. K_heat combines top layers if top layer too thin.
 *         14. SWR_bot_abs weights the absorbance to the spectral of light at the bottom.
 *         15. Interpolated Rrs onto finer grid for 2d sensor response curves.
 *         16. Interpolated E_d_z onto finer grid for 3d sensor response curves. 
 * 
 *  Things to do (excluding hydrolight):
 *
 *         1. Tracer list specification of pigment compositions
 *         3. Tracer list specification of fluorescence.
 *         13. Prevent products being calculated with zero Rrs.
 *         26. Give canopy order in the optical_setup.nc
 *         27. Check out-of-bounds reflectance from benthic substrates - some are zero.
 *         30. Coral skeleton reflectance is hardwired - need another tracer attribute.
 *         31. Suspended macroalgae is hardwired - need another tracer attribute.
 *         39. Colour_of_source_in_water (flu and bio) to add attribute to optical_setup.nc
 *         45. Add potassium decay as light source.
 *         47. Add passive fluorescence for all phytoplankton? Even Symbiodinium?
 *         53. Create function to read SWR at t+dt, perhaps using a tracer.
 *         63. Pigments are still in "csiro_siop_library.nc" - should be in csiro_optical_parameters_library.nc.
 *         64. Restart works for column.nc, but it produces extra times between restart time and last old write out.
 *         65. Make solar_azimuth from zenith a function.
 *         66. Should fulldisk be PAR-integrated?
 *         67. values_common_epi is still hardwired with some plant types.
 *         69. Time ranges and frequency for column output file.
 *         70. Code speed up: a. diffaa should do all wavelengths at one.
 *                            b. X x Omega(X) can be outside the wave loop.
 *                            c. m * red * MW_Nitr * 1000 could be done at initialisation.
 *                            d. 3D sensor - demon outside layer loop.
 *                            e. Secchi pathlength outside layer loop.
 *                            f. moon-phase once for all columns.
 *         71. Generalise fine sensor grid so it doesn't have to be 1 nm bandwidth.
 *         72. We don't have a simulated satellite products for Secchi.
 *         73. Light_spectral_col doesn't work for KEYWORD specification of tracers.
 *         74. Likely to be problems with using col->b to identify output columns in fully-coupled version.
 *         75. Put in ems version into optical_setup.nc file attributes.
 *         76. write_date_created(ncid1) only works for shoc, so commented out.
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
#include "einterface.h"
#include "light_spectral_col.h"
#include "constants.h"
#include "netcdf.h"
#include "ncw.h"

// HYDROLIGHT & BIOLUMINESCENCE are not part of the standard EMS release. Check-in with the
// below two lines commented out (PS. using a zero instead of one does not work).

// #define INCLUDE_HYDROLIGHT 1
// #define INCLUDE_BIOLUMINESCENCE 1

#ifdef INCLUDE_BIOLUMINESCENCE
#include "bioluminescence.h"
#endif

#ifdef INCLUDE_HYDROLIGHT

void hydrolight(long ndep, double wl, const double* depth,const double* a, const double* c, const double* bb, const double* s0);
void output_size(long nphi, long nmu);
void set_run(double wl);
void set_sky_irrad(double ed);
void set_geo(long iskyflag, double suntheta, double doy, double hr, double sunphi, double lat, double lon);
void set_atm(double windspd, double cld, double pres, double wsm, double am, double rh, double wv, double vi, double o3);
void set_water(double t, double refr, double s);
void set_bottom(long ibotm, double rflbot);
void set_depth_out(long ndep, const double *depth);

void copy_output_ns(long nmu, long nphi, long ndep, double* zen, double* azi,
                double* rad_sky, double* rad_ref, double* rad_wl, double *rad_dn, double *rad_up,
                double *irrad_sky, double *irrad_ref, double *irrad_wl, double *irrad_dn, double *irrad_up,
                double *sunz, double *suna);
void copy_output(long nmu, long nphi, long ndep, double* zen, double* azi,
                double* rad_sky, double* rad_ref, double* rad_wl, double *rad_dn, double *rad_up,
                double *irrad_sky, double *irrad_ref, double *irrad_wl, double *irrad_dn, double *irrad_up,
                double *sunz, double *suna, double irrad_et);
#endif

void ginterface_moonvars(void *hmodel, int b, double *mlon, double *mlat,
			 double *dist_earth_sun, double *dist_moon_earth,
			 double *lunar_angle, double *sun_angle,
			 double *moon_phase, double *lunar_dec);
void ginterface_get_lat_lon(void *hmodel, int b, double *lat, double *lon);
void ginterface_get_cloud_and_wind(void *hmodel, int b, double *tcld, double *wind);
void ginterface_get_windspeed_and_cloud_and_mslp_and_rh(void *hmodel, int b, double *windspeed, double *tcld, double *mslp, double *rh);
int  ginterface_getnumbersedlayers(void* model);
int  ginterface_getnumberwclayers(void* model);

double ginterface_getlighttop(void *model, int b);
double ginterface_calc_zenith(void *model, double t, int b);
void ginterface_get_ij(void* model, int col, int *ij);

char *ginterface_gettimeunits(void *model); // model time units.
char *ginterface_getoutputtimeunits(void *model); // model output time units.
char *ginterface_get_output_path(void);
int i_get_num_tracers_3d(void* hmodel, int *atr);
int i_get_num_tracers_2d(void* hmodel, int *atr);
void i_get_names_tracers_2d(void* hmodel, char **trname);
void i_get_names_tracers_3d(void* hmodel, char **trname);
void ginterface_get_tracerunits(void* model, char *name, char *units);
double ginterface_get_eta(void* hmodel, int b);

typedef struct {

  /* Environmental tracers */

  int temp_i;
  int salt_i;
  int eta_i;

  /* Optically-active water column tracers */

  int cdom_gbr_i;
  int cdom_pale_i;
  int cdom_amber_i;
  int cdom_dark_i;

  int PhyL_N_i;
  int PhyL_Chl_i;
  int PhyL_yCfac_i;
  int PhyL_kI_i;
  int PhyL_I_i;

  int PhyS_N_i;
  int PhyS_Chl_i;
  int PhyS_yCfac_i;
  int PhyS_kI_i;

  int PhyD_N_i;
  int PhyD_Chl_i;
  int PhyD_yCfac_i;
  int PhyD_kI_i;

  int Tricho_N_i;
  int Tricho_Chl_i;
  int Tricho_yCfac_i;
  int Tricho_kI_i;

  int MPB_N_i;
  int MPB_Chl_i;
  int MPB_yCfac_i;
  int MPB_kI_i;

  int CS_N_free_i;
  int CS_Chl_free_i;
  int CS_Xh_free_i;
  int CS_Xp_free_i;
  
  int CS_free_yCfac_i;
  int CS_free_kI_i;

  int SMA_N_i;
  int SMA_kI_i;

  /* Optically-active epibenthic variables */

  int CH_N_i;
  int CS_N_i;
  int CS_Chl_i;
  int CS_Xp_i;
  int CS_Xh_i;
  int CS_yCfac_i;
  int CS_kI_i;

  /* Optically-active benthic tracers */

  int MPB_N_sed_i;
  int MPB_Chl_sed_i;
  
  /* Microalgae parameters */

  double m_l;
  double rad_l;
  double vol_l;

  double m_s;
  double rad_s;
  double vol_s;

  double m_Tricho;
  double rad_Tricho;
  double vol_Tricho;

  double m_MPB;
  double rad_MPB;
  double vol_MPB;

  double m_PhyD;
  double rad_PhyD;
  double vol_PhyD;

  double m_CS;
  double rad_CS;
  double vol_CS;

  double MAleafden;

  double bbcell1;
  double bbcell2;

  double bbfact_l;
  double bbfact_s;
  double bbfact_MPB;
  double bbfact_Tricho;
  double bbfact_PhyD;
  double bbfact_CS;

  /* Benthic plant parameters */
  
  double CHpolypden;
  double CHarea;

  /* Water column non-spectrally-resolved non-living optical parameters */ 

  double gi;
  double gii;

  double g0;
  double g1;

  double SWRscale;

  double Scdom_pale;
  double Scdom_amber;
  double Scdom_dark;
  double Scdom_ocean;

  double a440cdom_pale;
  double a440cdom_amber;
  double a440cdom_dark;
  double a440cdom_ocean;

  /* Optical parameters requiring rethink */ 

  double bphy;
  double NtoCHL;

  /* Optical state variable outputs */

  int OC4Me_SR_i;
  int OC3M_SR_i;
  int OC3V_SR_i;
  int Hue_SR_i;
  int nFLH_SR_i;
  int TSSM_SR_i;
  int KD490M_SR_i;

  int EpiPAR_sg_i;
  int EpiPAR_i;

  int PAR_hyd_i;

  int SWR_bot_abs_i;
  int K_heat_i;
  int Secchi_i;
  int Zenith_i;
  int Azimuth_i;
  int Cloud_i;
  
  int PAR_i;
  int PAR_z_i;
  int Kd_490_i;
  int at_440_i;
  int bt_550_i;
  int bb_590_i;
  int bb_700_i;
  int BLP_i;
  int PFL_i;
  int Turbidity_i;

  int Kd_PAR_i;

  int PAR_bio_up_i;

  int Moonlight_i;
  int Lunar_zenith_i;
  int Lunar_azimuth_i;
  int Lunar_phase_i;
  int Moon_fulldisk_i;

  int Solar_zenith_i;
  
  /* Special wavelengths */
  
  int wPAR_top;
  int wPAR_bot;
  int w440;
  int w470;
  int w490;
  int w550;
  int w590;
  int w670;
  int w700;

  /* Indices of output */

  int *i_indices;
  int *j_indices;
  int num_cols_out;

  double PARwidth;
  double bls_sum;
  double fls_sum;

  /* wavelengths at which reflectances are calculated */

  int *RRS_i;
  int *wXXX_i;
  int num_rrs_waves;
  double *rrs_wave;

  /* wavelengths at which downwelling irradiance are calculated */

  int *Ed_i;
  int *wYYY_i;
  int num_ed_waves;
  double *ed_wave;

  /* wavelengths at which sensor response are integrated */

  int num_sensor_waves;
  double *sensor_waves;

  /* Spectrally-resolved optical parameters */

  double* landa;
  double* aw;
  double* bw;
  
  double* nbt_p_555;

  double *yC_s;
  double *yC_l;
  double *yC_Tricho;
  double *yC_MPB;
  double *yC_D;
  double *yC_Symbiodinium;

  // Individual pigments.

  double *yC_diatoxanthin;
  double *yC_diadinoxanthin;

  double *bls;   /* Normalised (to spectral peak) bioluminescence emission spectra */
  double *fls;   /* Normalised (to spectral peak) fluorescence emission spectra */

  double *rhoskel;
  double *SMA_absorption;

  /* Only need satellite bands that are put into derived products */

   /* Sentinel_3B spectral response function - 21 bands */
  
  int Sentinel_3B_Band3_i;
  int Sentinel_3B_Band4_i;
  int Sentinel_3B_Band5_i;
  int Sentinel_3B_Band6_i;
  int Sentinel_3B_Band8_i;
  int Sentinel_3B_Band11_i;

  /* MODIS spectral response function - 16 bands */
  
  int MODIS_Band1_i;
  int MODIS_Band2_i;
  int MODIS_Band4_i;
  int MODIS_Band6_i;
  int MODIS_Band9_i;
  int MODIS_Band10_i;
  int MODIS_Band11_i;

 /* VIIRS spectral response function - 10 bands */
  
  int VIIRS_Band2_i;
  int VIIRS_Band3_i;
  int VIIRS_Band4_i;

  /* Optical grid */

  double* bandedge;
  double* wave;
  double* bandwidth;
  int num_wave;
  int num_wc_layers;
  int num_sed_layers;

  char column_out_name[MAXSTRLEN];
  char* timeunits;
  char* outputtimeunits;

  /* Workspace for indices and parameter values of tracer-list specification 
     of optically-active tracers with mass-specific coefficients */

  int ntr_at;
  int* ind_at;
  double** at_star;  //  [ntr_at x num_wave] 

  int ntr_bt;
  int* ind_bt;
  double** bt_star;  //  [ntr_bt x num_wave]

  int ntr_bb;
  int* ind_bb;
  double** bb_star;  //  [ntr_bb x num_wave]

  int ntr_br;
  int* ind_br;
  double** br_star;  //  [ntr_br x num_wave]

  int ntr_sp3;
  int* ind_sp3;
  double** sp3_star;  //  [ntr_sp3 x num_wave]

  int nepi_A;
  int* ind_A;
  double** A_star;  //  [nepi_at x num_wave] 
  int* ind_kI;

  int nepi_R;
  int* ind_R;
  double** R_star;  //  [nepi_bt x num_wave] 
  
  int nepi_T;
  int* ind_T;
  double** T_star;  //  [nepi_bt x num_wave]

  int nepi_sp2;
  int* ind_sp2;
  double** sp2_star;  //  [nepi_sp2 x sensor_num_wave]

  int numMA; // number of macroalgae types
  int numSG; // number of seagrass types

  double* XXXleafden; // [nepi]
  double* XXXorient;  // [nepi]

  // Spectrally-resolved moonlight reflectance parameters.
  
  double *moon_a0;
  double *moon_a1;
  double *moon_a2;
  double *moon_a3;
  double *moon_b1;
  double *moon_b2;
  double *moon_b3;
  double *moon_d1;
  double *moon_d2;
  double *moon_d3;

  double K_heat_minthickness;

} workspace;

/* Functions listed at bottom of this file */

void moonlight(workspace *ws, column* col, double time_in_days, double* Ak, double *sun_angle, double *moon_phase, double *lunar_dec, double *lunar_angle, double *lunar_zenith, double *lunar_azimuth);
void optdatascalarread(ecology *e, workspace *ws, char *source_file, int ncid, const char *varname, double *out);
void optdataspectralread(ecology *e, workspace *ws, char *source_file, int ncid, const char *varname, const char *wavename, const char *codename, double* out);
void optdataspectralread_fine(ecology *e, workspace *ws, char *source_file, int ncid, const char *varname, const char *wavename, const char *codename, double* out);
void scale_sensor_response(ecology *e, workspace *ws, const char *sensorname, double* in, double* out, double *scale); // NOT TESTED PROPERLY.
void ModelRsestoBandRs(ecology *e, workspace *ws, double* response, double* modelR, double *obsR);
void OC4Me(double band3, double band4, double band5, double band6, double *oc4me_out);
void OC3M(double band2, double band4, double band6, double *oc3m_out);
void OC3V(double band3, double band4, double band6, double *oc3v_out);
void Hue3B(double band3, double band4, double band6, double band8, double band11, double *hue3b_out);
void nFLH_modis(double band9, double band10, double band11, double *nFLH_modis_out);
void TSSM(double band1, double *tssm_out);
void KD490M(double band4, double band6, double *kd490m_out);
void index_of_refraction(ecology *e, workspace *ws, double **y, double* n_stw);
void colour_in_seawater(ecology *e, workspace *ws, const char *varname);
void colour_of_bottom_type(ecology *e, workspace *ws, const char *varname);
void passive_fluorescence_per_cell(ecology *e,workspace *ws, double kI, double Iquota, double *fluorescence);
void usgs_moon_reflectance(ecology *e,workspace *ws);
void chl_specific_aggregate_pigment_absorption(ecology *e,workspace *ws);
void ems_swr_airtrans_and_albedo(double tcld, double lunar_dec, double lunar_angle, double lunar_zenith, double lat, double *airtransmission, double *swr_albedo);
void meb_swr_airtrans_and_albedo(double tcld, double lunar_dec, double lunar_angle, double lunar_zenith, double lat, double *airtransmission, double *swr_albedo);
void light_spectral_col_init(eprocess* p)
{
  ecology* e = p->ecology;
  stringtable* tracers = e->tracers;
  stringtable* epis = e->epis;
  workspace* ws = malloc(sizeof(workspace));
  int i, w;

  char string4netcdf[MAXSTRLEN];

  p->workspace = ws;

  // initialise all workspace elements to zero. 

  memset(ws,0,sizeof(ws));

  /* Read i,j output indices from biological parameter file - default  are LCJO in GBR4 */

  // Lucinda_Jetty 266/23 [in netcdf] col->b = 2548, 146.36 E, 18.52 S
                                     // well offshore 329/165, shallow offshore 176/70 //

  // Read in multiple cols for outputing into separate files - hydrolight only applies to first column listed.
  
  ws->num_cols_out = 0;

  if (try_parameter_value(e, "opt_output_i_index")>-1 && try_parameter_value(e, "opt_output_j_index")>-1){

    ws->num_cols_out = get_parameter_num_values(e,"opt_output_i_index");

    if (get_parameter_num_values(e,"opt_output_j_index") != ws->num_cols_out)
      e_quit("light_spectral_col: specified different number of opt_output_i_index & opt_output_j_index \n");
    
    ws->i_indices = i_alloc_1d(ws->num_cols_out);
    ws->j_indices = i_alloc_1d(ws->num_cols_out);
    
    for (w=0; w<ws->num_cols_out; w++) {
      ws->i_indices[w]= (int)get_parameter_value_ptr(e,"opt_output_i_index")[w];
      ws->j_indices[w] = (int)get_parameter_value_ptr(e,"opt_output_j_index")[w];
      eco_write_setup(e,"If wet, outputs of optical column at i_index = %d j_index = %d \n",ws->i_indices[w],ws->j_indices[w]);
    }
  }

  ws->timeunits = ginterface_gettimeunits(e->model);
  ws->outputtimeunits = ginterface_getoutputtimeunits(e->model);

  /* create matlab plot file */

  FILE *fp1;
  fp1 = fopen("optical_setup_matlab_check.m", "w+");
  fprintf(fp1,"%% plotting file generated by CSIRO EMS for plotting spectrally-resolved optical properties.\n");
  fprintf(fp1,"close all;clear all;\n");
  fprintf(fp1,"file = 'optical_setup.nc';\n");
  fprintf(fp1,"wave = ncread(file,'wave_centre');\n\n");
  fprintf(fp1,"wave_sensor = ncread(file,'wave_sensor');\n\n");
  fclose(fp1);

  /* create phyton plot file */

  char cwd[200];
  getcwd(cwd,sizeof(cwd));
    
  FILE *fp2;
  fp2 = fopen("optical_setup_python_check.py", "w+");
  fprintf(fp2,"## plotting file generated by CSIRO EMS for plotting spectrally-resolved optical properties.\n\n");
  
  fprintf(fp2,"import matplotlib\n");
  fprintf(fp2,"import matplotlib.pyplot as plt\n");
  fprintf(fp2,"import xarray as xar\n\n");
  
  fprintf(fp2,"dpath ='%s'\n", cwd);
  fprintf(fp2,"file = f'{dpath}/optical_setup.nc'\n");
  fprintf(fp2,"XR = xar.open_mfdataset(file, decode_times= True)\n");
  fprintf(fp2,"wave = XR.wave_centre[:].values\n\n");
  fclose(fp2);

  /* Set-up optical grid based on values in biological parameter file. */

  ws->num_wave = get_parameter_num_values(e, "Light_lambda");
  ws->wave = get_parameter_value_ptr(e, "Light_lambda");

  /* Set-up sensor response grid with Rrs_fine - only works for 1 nm at the moment. */

  int ii;
  ws->num_sensor_waves = 501;
  ws->sensor_waves = d_alloc_1d(ws->num_sensor_waves);
  ws->sensor_waves[0] = 300.0;
  
  for (ii = 1; ii<ws->num_sensor_waves; ii++) {
    ws->sensor_waves[ii] = ws->sensor_waves[0]+ii;
  }

  /***********************************************************************************************/
  /* This section looks for the specification of optical properties in the tracer list and loads 
     the data into the workspace. */
  
  int dummy;int dummy111;
    
  char **tnames;       /* 3D tracer names                         */
  char **epinames;       /* 2D tracer names                         */
  
  int ntnames = i_get_num_tracers_3d(e->model,&dummy111);
  
  tnames = malloc(ntnames * sizeof(char *));
  for (i = 0; i < ntnames ; i++)
    tnames[i] = malloc(MAXSTRLEN * sizeof(char));
  
  i_get_names_tracers_3d(e->model,tnames);
  
  for (i = 0; i < ntnames; i++) {  // Loop once to get the size of the arrays
    if (einterface_is_optical(e->model,tnames[i])){
      dummy = e->find_index(tracers,tnames[i],e);
    }
  }
  
  int nepinames = i_get_num_tracers_2d(e->model,&dummy111);
  
  epinames = malloc(nepinames * sizeof(char *));
  for (i = 0; i < nepinames; i++)
    epinames[i] = malloc(MAXSTRLEN * sizeof(char));
  
  i_get_names_tracers_2d(e->model,epinames);
  
  for (i = 0; i < nepinames; i++) {  // Loop once to get the size of the arrays
    if (einterface_is_optical2d(e->model,epinames[i])){
      dummy = e->find_index(epis,epinames[i],e);
      }
    }
  
  ws->ntr_at = 0;
  ws->ntr_bt = 0;
  ws->ntr_bb = 0;
  ws->ntr_br = 0;
  ws->ntr_sp3 = 0;

  ws->nepi_A = 0;
  ws->nepi_R = 0;
  ws->nepi_T = 0;
  ws->nepi_sp2 = 0;

  ws->at_star = NULL;
  ws->bt_star = NULL;
  ws->bb_star = NULL;
  ws->br_star = NULL;
  ws->sp3_star = NULL;

  ws->A_star = NULL;
  ws->R_star = NULL;
  ws->T_star = NULL;
  ws->sp2_star = NULL;

  ws->XXXleafden = NULL;
  ws->XXXorient = NULL;

  char trname[MAXSTRLEN];
  char fname[MAXSTRLEN];
  char absorp_name[MAXSTRLEN];
  char scatter_name[MAXSTRLEN];
  char backscat_name[MAXSTRLEN];
  char benreflt_name[MAXSTRLEN];
  char specresp3d_name[MAXSTRLEN];

  trname[0] = '\0';
        
  for (i = 0; i < e->ntr; i++) {  // Loop once to get the size of the arrays

    strcpy(trname,  e->tracernames[i]);

    if (einterface_is_optical(e->model,trname)){ 
    
      if (einterface_get_absorp_name(e->model, trname, absorp_name))
	ws->ntr_at = ws->ntr_at + 1;
      
      if (einterface_get_scatter_name(e->model, trname, scatter_name))
	ws->ntr_bt = ws->ntr_bt + 1;
    
      if (einterface_get_backscat_name(e->model, trname, backscat_name))
	ws->ntr_bb = ws->ntr_bb + 1;
      
      if (einterface_get_benreflt_name(e->model, trname, benreflt_name))
	ws->ntr_br = ws->ntr_br + 1;

      if (einterface_get_specresp3d_name(e->model, trname, specresp3d_name))
	ws->ntr_sp3 = ws->ntr_sp3 + 1;
    }
  }
  
  eco_write_setup(e,"Light_spectral_col: Found %d ats, %d bts, %d bbs, %d brs, %d sp and %d waves \n",ws->ntr_at,ws->ntr_bt,ws->ntr_bb,ws->ntr_br,ws->ntr_sp3,ws->num_wave);

  // Test that all tracer attributes are specified.

  if ((ws->ntr_br != ws->ntr_bb) || (ws->ntr_br != ws->ntr_at) || (ws->ntr_br != ws->ntr_bt))
    e_quit("eco: stopped because number of tracer attributes different for at, bt, bb, benrefl \n");
  
  // Do the same for 2D tracers.

  char epiname[MAXSTRLEN] = "";

  char abtance_name[MAXSTRLEN];
  char refltce_name[MAXSTRLEN];
  char trnmiss_name[MAXSTRLEN];
  char specresp2d_name[MAXSTRLEN];

  for (i = 0; i < e->nepi; i++) {  // Loop once to get the size of the arrays

    strcpy(epiname,  e->epinames[i]);

    if (einterface_is_optical2d(e->model,epiname)){
    
      if (einterface_get_abtance_name(e->model, epiname, abtance_name))
      	ws->nepi_A = ws->nepi_A + 1;
      
      if (einterface_get_refltce_name(e->model, epiname, refltce_name))
       	ws->nepi_R = ws->nepi_R + 1;
    
      if (einterface_get_trnmiss_name(e->model, epiname, trnmiss_name))
      	ws->nepi_T = ws->nepi_T + 1;

      if (einterface_get_specresp2d_name(e->model, epiname, specresp2d_name))
      	ws->nepi_sp2 = ws->nepi_sp2 + 1;
    }
  }

  eco_write_setup(e,"Light_spectral_col: Found %d As, %d Rs, %d Ts, %d sp and %d waves \n",ws->nepi_A,ws->nepi_R,ws->nepi_T,ws->nepi_sp2,ws->num_wave);

  // Now allocate memory if tracer-specified optical-active tracer found.
  
  if (ws->ntr_at){
    ws->ind_at = i_alloc_1d(ws->ntr_at);
    ws->at_star = d_alloc_2d(ws->num_wave,ws->ntr_at);
  }
  if (ws->ntr_bt){
    ws->ind_bt = i_alloc_1d(ws->ntr_bt);
    ws->bt_star = d_alloc_2d(ws->num_wave,ws->ntr_bt);
  }
  if (ws->ntr_bb){
    ws->ind_bb = i_alloc_1d(ws->ntr_bb);
    ws->bb_star = d_alloc_2d(ws->num_wave,ws->ntr_bb);
  }
  if (ws->ntr_br){
    ws->ind_br = i_alloc_1d(ws->ntr_br);
    ws->br_star = d_alloc_2d(ws->num_wave,ws->ntr_br);
  }
  if (ws->ntr_sp3){
    ws->ind_sp3 = i_alloc_1d(ws->ntr_sp3);
    ws->sp3_star = d_alloc_2d(ws->num_sensor_waves,ws->ntr_sp3);
  }

  // Presently assumes all nepi are A.

  if (ws->nepi_A){
    ws->ind_A = i_alloc_1d(ws->nepi_A);
    ws->ind_kI = i_alloc_1d(ws->nepi_A);
    ws->A_star = d_alloc_2d(ws->num_wave,ws->nepi_A);
  }
  if (ws->nepi_R){
    ws->ind_R = i_alloc_1d(ws->nepi_R);
    ws->R_star = d_alloc_2d(ws->num_wave,ws->nepi_R);
  }
  if (ws->nepi_T){
    ws->ind_T = i_alloc_1d(ws->nepi_T);
    ws->T_star = d_alloc_2d(ws->num_wave,ws->nepi_T);
  }
  if (ws->nepi_sp2){
    ws->ind_sp2 = i_alloc_1d(ws->nepi_sp2);
    ws->sp2_star = d_alloc_2d(ws->num_sensor_waves,ws->nepi_sp2);
  }
  if (ws->nepi_A){
    ws->XXXleafden = d_alloc_1d(ws->nepi_A);
  }
  if (ws->nepi_A){
    ws->XXXorient = d_alloc_1d(ws->nepi_A);
  }

  /* End of tracer-list specified optically-active tracers. 
  /***********************************************************************************************/

    /* Optically-active tracers that are hardwired */

  ws->PhyL_N_i = -1; ws->PhyL_Chl_i = -1; ws->PhyL_yCfac_i = -1; ws->PhyL_kI_i = -1;ws->PhyL_I_i = -1;
  ws->PhyS_N_i = -1; ws->PhyS_Chl_i = -1; ws->PhyS_yCfac_i = -1; ws->PhyS_kI_i = -1;
  ws->PhyD_N_i = -1; ws->PhyD_Chl_i = -1; ws->PhyD_yCfac_i = -1; ws->PhyD_kI_i = -1;
  ws->MPB_N_i = -1;  ws->MPB_Chl_i = -1;  ws->MPB_yCfac_i = -1;  ws->MPB_kI_i = -1;
  ws->Tricho_N_i = -1; ws->Tricho_Chl_i = -1; ws->Tricho_yCfac_i = -1; ws->Tricho_kI_i = -1;
  ws->CS_N_free_i = -1;  ws->CS_Chl_free_i = -1;  ws->CS_free_yCfac_i = -1;  ws->CS_free_kI_i = -1;
  ws->CS_Xh_free_i = -1;  ws->CS_Xp_free_i = -1;
  ws->SMA_N_i = -1; ws->SMA_kI_i = -1;
  
  ws->SMA_N_i = e->try_index(tracers, "SMA_N", e) ;
  ws->SMA_kI_i = e->try_index(tracers, "SMA_kI", e) ;
  
  ws->PhyL_N_i = e->try_index(tracers, "PhyL_N", e);
  if (ws->PhyL_N_i > -1){
    ws->PhyL_Chl_i = e->find_index(tracers, "PhyL_Chl", e);
    if (!e->pre_build){
      ws->PhyL_yCfac_i = e->find_index(tracers, "PhyL_yCfac", e);
      ws->PhyL_kI_i = e->find_index(tracers, "PhyL_kI", e);
    }
    ws->PhyL_I_i = e->find_index(tracers, "PhyL_I", e);
  }

  ws->PhyS_N_i = e->try_index(tracers, "PhyS_N", e);
  if (ws->PhyS_N_i > -1){
    ws->PhyS_Chl_i = e->find_index(tracers, "PhyS_Chl", e);
    if (!e->pre_build){
      ws->PhyS_yCfac_i = e->find_index(tracers, "PhyS_yCfac", e);
      ws->PhyS_kI_i = e->find_index(tracers, "PhyS_kI", e);
    }
  }

  ws->PhyD_N_i = e->try_index(tracers, "PhyD_N", e);
  if (ws->PhyD_N_i > -1){
    ws->PhyD_Chl_i = e->try_index(tracers, "PhyD_Chl", e);
    if (!e->pre_build){
      ws->PhyD_yCfac_i = e->try_index(tracers, "PhyD_yCfac", e);
      ws->PhyD_kI_i = e->try_index(tracers, "PhyD_kI", e);
    }
  }
  
  ws->Tricho_N_i = e->try_index(tracers, "Tricho_N", e);
  if (ws->Tricho_N_i > -1){
    ws->Tricho_Chl_i = e->find_index(tracers, "Tricho_Chl", e);
    if (!e->pre_build){
      ws->Tricho_yCfac_i = e->find_index(tracers, "Tricho_yCfac", e);
      ws->Tricho_kI_i = e->find_index(tracers, "Tricho_kI", e);
    }
  }
  
  ws->MPB_N_i = e->try_index(tracers, "MPB_N", e);
  if (ws->MPB_N_i > -1){
    ws->MPB_Chl_i = e->find_index(tracers, "MPB_Chl", e);
    if (!e->pre_build){
      ws->MPB_yCfac_i = e->find_index(tracers, "MPB_yCfac", e);
      ws->MPB_kI_i = e->find_index(tracers, "MPB_kI", e);
    }
  }
    
  ws->CS_N_free_i = e->try_index(tracers, "CS_N_free", e);
  if (ws->CS_N_free_i > -1){
    ws->CS_Chl_free_i = e->try_index(tracers, "CS_Chl_free", e);
    ws->CS_Xh_free_i = e->try_index(tracers, "CS_Xh_free", e);
    ws->CS_Xp_free_i = e->try_index(tracers, "CS_Xp_free", e);
    if (!e->pre_build){
      ws->CS_free_yCfac_i = e->try_index(tracers, "CS_free_yCfac", e);
      ws->CS_free_kI_i = e->try_index(tracers, "CS_free_kI", e);
    }
  }

  ws->MPB_N_sed_i = -1;  ws->MPB_Chl_sed_i = -1;
  ws->MPB_N_sed_i = e->try_index(tracers, "MPB_N", e);
  if (ws->MPB_N_sed_i > -1){
    ws->MPB_Chl_sed_i = e->find_index(tracers, "MPB_Chl", e);
  }
  
  ws->Secchi_i =-1;
  ws->Secchi_i = e->try_index(epis, "Secchi", e);
  
  ws->EpiPAR_sg_i =-1;
  ws->EpiPAR_sg_i = e->try_index(epis, "EpiPAR_sg", e);

  ws->EpiPAR_i =-1;
  ws->EpiPAR_i = e->try_index(epis, "EpiPAR", e);

  ws->Zenith_i = -1;
  ws->Zenith_i = e->try_index(epis, "Zenith2D", e);

  ws->Azimuth_i = -1;
  ws->Azimuth_i = e->try_index(epis, "Azimuth", e);

  ws->Cloud_i = -1;
  ws->Cloud_i = e->try_index(epis, "Cloud", e);

  /* Output for tracers for comparison */

  ws->Kd_490_i = -1;
  ws->Kd_490_i =  e->try_index(tracers, "Kd_490", e);

  ws->Kd_PAR_i = -1;
  ws->Kd_PAR_i =  e->try_index(tracers, "Kd_PAR", e);

  ws->bt_550_i = -1;
  ws->bt_550_i =  e->try_index(tracers, "bt_550", e);

  ws->at_440_i = -1;
  ws->at_440_i =  e->try_index(tracers, "at_440", e);

  ws->bb_590_i = -1;
  ws->bb_590_i =  e->try_index(tracers, "bb_590", e);

  ws->Turbidity_i = -1;
  ws->Turbidity_i   = e->try_index(tracers, "Turbidity", e);
  if (ws->bb_590_i == -1)
    ws->Turbidity_i = -1; // needs bb_590 to calculate

  ws->bb_700_i = -1;
  ws->bb_700_i =  e->try_index(tracers, "bb_700", e);
  
  ws->PAR_i = -1;
  ws->PAR_i =  e->try_index(tracers, "PAR", e);

  ws->PAR_z_i = -1;
  ws->PAR_z_i =  e->try_index(tracers, "PAR_z", e);

  ws->PAR_bio_up_i = -1;
  ws->PAR_bio_up_i =  e->try_index(tracers, "PAR_bio_up", e);

  ws->K_heat_i = -1;
  ws->K_heat_i =  e->try_index(tracers, "K_heat", e);

  ws->BLP_i = -1;
  ws->BLP_i =  e->try_index(tracers, "BLP", e); 

  ws->PFL_i = -1;
  ws->PFL_i =  e->try_index(tracers, "PFL", e);


#ifdef INCLUDE_HYDROLIGHT
  eco_write_setup(e,"Doing HydroLight calculations \n");

  // Default hydrolight options.
      
  double zengrid[10] = {0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,87.5};
  double azigrid[24] = {0.0,15.0,30.0,45.0,60.0,75.0,90.0,105.0,120.0,135.0,150.0,165.0,180.0,195.0,210.0,225.0,240.0,255.0,270.0,285.0,300.0,315.0,330.0,345.0};

  ws->PAR_hyd_i =  e->try_index(tracers, "PAR_hyd", e);
#endif

  ws->Sentinel_3B_Band3_i = -1;
  ws->Sentinel_3B_Band4_i = -1;
  ws->Sentinel_3B_Band5_i = -1;
  ws->Sentinel_3B_Band6_i = -1;
  ws->Sentinel_3B_Band8_i = -1;
  ws->Sentinel_3B_Band11_i = -1;
  ws->OC4Me_SR_i = -1;
  ws->Hue_SR_i = -1;

  if (!e->pre_build){

    ws->Sentinel_3B_Band3_i = e->try_index(epis,"Sentinel_3B_B3",e);
    ws->Sentinel_3B_Band4_i = e->try_index(epis,"Sentinel_3B_B4",e);
    ws->Sentinel_3B_Band5_i = e->try_index(epis,"Sentinel_3B_B5",e);
    ws->Sentinel_3B_Band6_i = e->try_index(epis,"Sentinel_3B_B6",e);
    ws->Sentinel_3B_Band8_i = e->try_index(epis,"Sentinel_3B_B8",e);
    ws->Sentinel_3B_Band11_i = e->try_index(epis,"Sentinel_3B_B11",e);
    ws->OC4Me_SR_i = e->try_index(epis,"OC4Me_SR",e);
    ws->Hue_SR_i = e->try_index(epis,"Hue_S3B_SR",e);
    
    if ((ws->Sentinel_3B_Band3_i > -1) && (ws->Sentinel_3B_Band4_i > -1) && (ws->Sentinel_3B_Band5_i > -1) && (ws->Sentinel_3B_Band6_i > -1) && (ws->OC4Me_SR_i > -1)){
      eco_write_setup(e,"OC4Me is being calculated using Sentinel_3B Bands 3, 4, 5 and 6 \n");
    }else{
      ws->OC4Me_SR_i = -1;
    }
    if ((ws->Sentinel_3B_Band3_i > -1) && (ws->Sentinel_3B_Band4_i > -1) && (ws->Sentinel_3B_Band6_i > -1) && (ws->Sentinel_3B_Band8_i > -1) && (ws->Sentinel_3B_Band11_i > -1) && (ws->Sentinel_3B_Band11_i > -1) && (ws->Hue_SR_i > -1)){
	eco_write_setup(e,"Hue is being calculated using Sentinel_3B Bands 3, 4, 5, 6, 8 and 11 \n");
    }else{
      ws->Hue_SR_i = -1;
    }
  }
  
  ws->MODIS_Band1_i = -1;
  ws->MODIS_Band2_i = -1;
  ws->MODIS_Band4_i = -1; 
  ws->MODIS_Band6_i = -1;
  ws->MODIS_Band9_i = -1;
  ws->MODIS_Band10_i = -1; 
  ws->MODIS_Band11_i = -1;
  ws->OC3M_SR_i = -1;
  ws->nFLH_SR_i = -1;
  ws->TSSM_SR_i = -1;
  ws->KD490M_SR_i = -1;

  if (!e->pre_build){

    ws->MODIS_Band1_i = e->try_index(epis,"MODIS_B1",e);

    ws->MODIS_Band2_i = e->try_index(epis,"MODIS_B2",e);
    ws->MODIS_Band4_i = e->try_index(epis,"MODIS_B4",e);
    ws->MODIS_Band6_i = e->try_index(epis,"MODIS_B6",e);

    ws->MODIS_Band9_i = e->try_index(epis,"MODIS_B9",e);
    ws->MODIS_Band10_i = e->try_index(epis,"MODIS_B10",e);
    ws->MODIS_Band11_i = e->try_index(epis,"MODIS_B11",e);

    ws->OC3M_SR_i = e->try_index(epis,"OC3M_SR",e);
    ws->nFLH_SR_i = e->try_index(epis,"nFLH_SR",e);
    ws->TSSM_SR_i = e->try_index(epis,"TSSM_SR",e);
    ws->KD490M_SR_i = e->try_index(epis,"KD490M_SR",e);

    // reset derived OC products index to-1 if the correct bands are not present.
    
    if ((ws->MODIS_Band2_i > -1) && (ws->MODIS_Band4_i > -1) && (ws->MODIS_Band6_i > -1) && (ws->OC3M_SR_i > -1)){
      eco_write_setup(e,"OC3M is being calculated using MODIS Bands 2, 4 and 6 \n");
    }else{
      ws->OC3M_SR_i = -1;
    }

    if ((ws->MODIS_Band9_i > -1) && (ws->MODIS_Band10_i > -1) && (ws->MODIS_Band11_i > -1) && (ws->nFLH_SR_i > -1)){
      eco_write_setup(e,"nFLH is being calculated using MODIS Bands 9, 10, 11 \n");
    }else{
      ws->nFLH_SR_i = -1;
    }

    if ((ws->MODIS_Band9_i) && (ws->TSSM_SR_i > -1)){
      eco_write_setup(e,"TSSM is being calculated using MODIS Band 1.\n");
    }else{
      ws->TSSM_SR_i = -1;
    }
    
    if ((ws->MODIS_Band4_i) && (ws->MODIS_Band6_i) && (ws->KD490M_SR_i > -1)){
	eco_write_setup(e,"KD490M is being calculated using MODIS Band 4 and 6.\n");
    }else{
      ws->KD490M_SR_i = -1;
    }
  }

  ws->VIIRS_Band2_i = -1;
  ws->VIIRS_Band3_i = -1;
  ws->VIIRS_Band4_i = -1;
  ws->OC3V_SR_i = -1;

  if (!e->pre_build){
    
    ws->VIIRS_Band2_i = e->try_index(epis,"VIIRS_B2",e);
    ws->VIIRS_Band3_i = e->try_index(epis,"VIIRS_B3",e);
    ws->VIIRS_Band4_i = e->try_index(epis,"VIIRS_B4",e);
    ws->OC3V_SR_i = e->try_index(epis,"OC3V_SR",e);
    
    if ((ws->VIIRS_Band2_i > -1) && (ws->VIIRS_Band3_i > -1) && (ws->VIIRS_Band4_i > -1) && (ws->OC3V_SR_i > -1)){
	eco_write_setup(e,"OC3V is being calculated using VIIRS Bands 2, 3 and 4 \n");
    }else{
      ws->OC3V_SR_i = -1;
    }
  }

  ws->Moonlight_i = -1;
  ws->Lunar_zenith_i = -1;
  ws->Lunar_azimuth_i = -1;
  ws->Lunar_phase_i = -1;
  ws->Moon_fulldisk_i = -1;

  ws->Solar_zenith_i = -1;

  ws->Moonlight_i = e->try_index(epis, "Moonlight", e);
  ws->Lunar_zenith_i = e->try_index(epis, "Lunar_zenith", e);
  ws->Lunar_azimuth_i = e->try_index(epis, "Lunar_azimuth", e);
  ws->Lunar_phase_i = e->try_index(epis, "Lunar_phase", e);
  ws->Moon_fulldisk_i = e->try_index(epis, "Moon_fulldisk", e);

  ws->Solar_zenith_i = e->try_index(epis, "Solar_zenith", e);
  
  ws->SWR_bot_abs_i = -1;
  ws->SWR_bot_abs_i = e->try_index(epis, "SWR_bot_abs", e);
  
  /* Epi-benthic tracers */

  ws->CH_N_i = -1;
  ws->CS_N_i = -1;
  ws->CS_Chl_i = -1;
  ws->CS_Xp_i = -1;
  ws->CS_Xh_i = -1;

  ws->CH_N_i = e->try_index(epis, "CH_N", e);
  ws->CS_N_i = e->try_index(epis, "CS_N", e);
  ws->CS_Chl_i = e->try_index(epis, "CS_Chl", e);
  ws->CS_Xp_i = e->try_index(epis, "CS_Xp", e);
  ws->CS_Xh_i = e->try_index(epis, "CS_Xh", e);

  ws->CS_yCfac_i = -1;
  ws->CS_kI_i = -1;
  if (!e->pre_build){
    ws->CS_yCfac_i = e->try_index(epis, "CS_yCfac", e);
    ws->CS_kI_i = e->try_index(epis, "CS_kI", e);
  }
  /* CDOM tracers */

  ws->cdom_pale_i = -1;
  ws->cdom_amber_i = -1;
  ws->cdom_dark_i = -1;
  ws->cdom_gbr_i = -1;

  /* Look for tracers and then find parameters in netcdf file */

  ws->cdom_pale_i =  e->try_index(tracers, "CDOM_pale", e);
  ws->cdom_amber_i =  e->try_index(tracers, "CDOM_amber", e);
  ws->cdom_dark_i =  e->try_index(tracers, "CDOM_dark", e);
  ws->cdom_gbr_i =  e->try_index(tracers, "cdom_gbr", e);

  ws->salt_i =  e->find_index(tracers, "salt", e);
  ws->temp_i =  e->find_index(tracers, "temp", e);

  /* Extra parameter for DA perturbations */
  
  ws->SWRscale = try_parameter_value(e, "SWRscale");
  if (isnan(ws->SWRscale)){
    ws->SWRscale = 1.0;
    eco_write_setup(e,"Code default of SWRscale = %e \n",ws->SWRscale);
  }

  ws->K_heat_minthickness = try_parameter_value(e, "K_heat_minthickness");
  if (isnan(ws->K_heat_minthickness)){
    ws->K_heat_minthickness = 0.10;
    eco_write_setup(e,"Code default of K_heat_minthickness = %e \n",ws->K_heat_minthickness);
  }
  
  /* Microalgae parameters */

  ws->bbcell1 = try_parameter_value(e, "bbcell1");
  if (isnan(ws->bbcell1)){
    ws->bbcell1 = 4.26e-14;   // Whitmore 2012
    eco_write_setup(e,"Code default of phyto. cell backscatter parameter 1, bbcell1 = %e \n",ws->bbcell1);
  }
  ws->bbcell2 = try_parameter_value(e, "bbcell2");
  if (isnan(ws->bbcell2)){
    ws->bbcell2 = 2.028;   // Whitmore 2012
    eco_write_setup(e,"Code default of phyto. cell backscatter parameter 2, bbcell2 = %e \n",ws->bbcell2);
  }
  
  if (ws->PhyL_N_i > -1){
    ws->m_l = PhyCellMass(try_parameter_value(e, "PLrad"));
    ws->rad_l = get_parameter_value(e, "PLrad");
    ws->vol_l = 4.0 / 3.0 * M_PI * ws->rad_l * ws->rad_l * ws->rad_l;
    ws->bbfact_l = ws->bbcell1 * pow(ws->rad_l*2e6, ws->bbcell2);
  }
  
  if (ws->PhyS_N_i > -1){
    ws->m_s = PhyCellMass(try_parameter_value(e, "PSrad"));
    ws->rad_s = get_parameter_value(e, "PSrad");
    ws->vol_s = 4.0 / 3.0 * M_PI * ws->rad_s * ws->rad_s * ws->rad_s;
    ws->bbfact_s = ws->bbcell1 * pow(ws->rad_s*2e6, ws->bbcell2);
  }

  if (ws->MPB_N_i > -1){
    ws->m_MPB = PhyCellMass(try_parameter_value(e, "MBrad"));
    ws->rad_MPB = get_parameter_value(e, "MBrad");
    ws->vol_MPB = 4.0 / 3.0 * M_PI * ws->rad_MPB * ws->rad_MPB * ws->rad_MPB;
    ws->bbfact_MPB = ws->bbcell1 * pow(ws->rad_MPB*2e6, ws->bbcell2);
  }

  if (ws->Tricho_N_i > -1){
    ws->m_Tricho = PhyCellMass(try_parameter_value(e, "Tricho_rad"));
    ws->rad_Tricho = get_parameter_value(e, "Tricho_rad");
    ws->vol_Tricho = 4.0 / 3.0 * M_PI * ws->rad_Tricho * ws->rad_Tricho * ws->rad_Tricho;
    ws->bbfact_Tricho = ws->bbcell1 * pow(ws->rad_Tricho*2e6, ws->bbcell2);
  }

  if (ws->PhyD_N_i > -1){
    ws->m_PhyD = PhyCellMass(try_parameter_value(e, "DFrad"));
    ws->rad_PhyD = get_parameter_value(e, "DFrad");
    ws->vol_PhyD = 4.0 / 3.0 * M_PI * ws->rad_PhyD * ws->rad_PhyD * ws->rad_PhyD;
    ws->bbfact_PhyD = ws->bbcell1 * pow(ws->rad_PhyD*2e6, ws->bbcell2);
  }

  if (ws->CS_N_i > -1){

    /* Only in coral initialisations if there is a coral tracer */

    ws->CHarea = try_parameter_value(e, "CHarea");
    if (isnan(ws->CHarea))
      ws->CHarea = 1.0;
    
    ws->CHpolypden = get_parameter_value(e, "CHpolypden");
    ws->rad_CS = get_parameter_value(e, "CSrad");
    ws->m_CS = PhyCellMass(ws->rad_CS);
    ws->vol_CS = 4.0 / 3.0 * M_PI * ws->rad_CS * ws->rad_CS * ws->rad_CS;
    ws->bbfact_CS = ws->bbcell1 * pow(ws->rad_CS*2e6, ws->bbcell2);
    
  }

  if (ws->SMA_N_i > -1){
    
    /* Put in macroalgae initialisations if there is a macroalgae tracer */
    
    ws->MAleafden = get_parameter_value(e, "MAleafden_wc");
    
  }

  // Initialise parameter values and kI index for tracer-list specifed benthis plants.

  ws->numMA = 0;
  ws->numSG = 0;
  
  char tempstr[MAXSTRLEN];
  char epiname1[MAXSTRLEN] = "";
  int aaa = -1;
  int is_plant;
  
  for (i = 0; i < e->nepi; i++) {  // Loop once to get the size of the arrays
    
    is_plant = 0;
    
    strcpy(epiname1,  e->epinames[i]);
    
    if (einterface_is_optical2d(e->model,epiname1)){
      
      eco_write_setup(e,"Setting paramaters for optical epi tracer %s \n",epiname1);
      
      char tempstr1[strlen(epiname1)-2];
      
      strncpy(tempstr1,epiname1,strlen(epiname1)-2); // remove the "_N"
      
      tempstr1[strlen(epiname1)-2] = '\0';
      
      if ((tempstr1[0] == 'M') && (tempstr1[1] == 'A')){
	ws->numMA += 1;
	is_plant = 1;
      }
      
      if ((tempstr1[0] == 'S') && (tempstr1[1] == 'G')){
	ws->numSG += 1;
	is_plant = 1;
      }
      
      if (is_plant){
	
	aaa += 1;
	
	sprintf(tempstr, "%sleafden",tempstr1);
	ws->XXXleafden[aaa] = get_parameter_value(e, tempstr);
	
	sprintf(tempstr, "%sorient", tempstr1);
	ws->XXXorient[aaa] = try_parameter_value(e, tempstr);
	
	if (isnan(ws->XXXorient[aaa])){
	  ws->XXXorient[aaa] = 1.0;
	  eco_write_setup(e,"In tracer-specified optical properties: Code default of %s = %e \n",tempstr,ws->XXXorient[aaa]);
	}
	sprintf(tempstr, "%s_kI", tempstr1);
	ws->ind_kI[aaa] = find_index(epis,tempstr,e);
	
      }else{   // Since not MA*_N or SG*_N, must be a spectral response to 2D light field.

	// not a grass or weed.
	
      }
    }
  }
  eco_write_setup(e,"Number of macroalgae types found in tracer list: %d \n",ws->numMA);
  eco_write_setup(e,"Number of seagrass types found in tracer list: %d \n",ws->numSG);

  /* Optical parameters requiring rethink */ 

  ws->bphy = try_parameter_value(e, "bphy");
  if (isnan(ws->bphy)){
    ws->bphy = 0.2;
    eco_write_setup(e,"Code default of bphy = %e \n",ws->bphy);
  }
  
  ws->NtoCHL = try_parameter_value(e, "NtoCHL");
  if (isnan(ws->NtoCHL)){
    ws->NtoCHL = 7.0;
    eco_write_setup(e,"Code default of NtoCHL = %e \n",ws->NtoCHL);
  }

  /* Optical properties coming from the bio.prm file */

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

  ws->Scdom_ocean = try_parameter_value(e, "Scdom_ocean");
  if (isnan(ws->Scdom_ocean)){
    ws->Scdom_ocean = 0.014;
    eco_write_setup(e,"Code default of Scdom_ocean = %e \n",ws->Scdom_ocean);
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
  ws->a440cdom_ocean = try_parameter_value(e, "a440cdom_ocean");
  if (isnan(ws->a440cdom_ocean)){
    ws->a440cdom_ocean = 0.01;
    eco_write_setup(e,"Code default of a440cdom_ocean = %e \n",ws->a440cdom_ocean);
  }

  /* Do a sanity check of wavelengths */

  for (w=0; w<ws->num_wave-1; w++){
    if (ws->wave[w+1] <= ws->wave[w])
     e->quitfn("Optical grid in not monotonically increasing. Check biological parameter file.");
  }

  // Allocate and fill in

  ws->bandedge = d_alloc_1d(ws->num_wave+1);
  ws->bandwidth = d_alloc_1d(ws->num_wave);

  /* Edges of wavebands depends on allocation of solar radiation */

  ws->bandedge[0] = 140.0;

  ws->bandedge[ws->num_wave] = 8000.0;
 
  for (w=0; w < ws->num_wave-1; w++){
    ws->bandedge[w+1] = (ws->wave[w]+ws->wave[w+1])/2.0;
    ws->bandwidth[w] = ws->bandedge[w+1] - ws->bandedge[w];
  }
  ws->bandwidth[ws->num_wave-1] = ws->bandedge[ws->num_wave] - ws->bandedge[ws->num_wave-1];
  
  eco_write_setup(e,"*********************************************************** \n");    
  eco_write_setup(e,"Optical grid calculated in light_spectral_col_init.c \n \n");
  eco_write_setup(e,"Centre of waveband (edges) waveband thickness: \n");
  for (w=0; w < ws->num_wave; w++){
    eco_write_setup(e,"%4.2f (%4.2f - %4.2f) %4.2f \n",ws->wave[w],ws->bandedge[w],ws->bandedge[w+1],ws->bandwidth[w]);
  }

   /* Search for special wavelengths */

  ws->wPAR_bot = -1;
  ws->w440 = -1;
  ws->w470 = -1;
  ws->w490 = -1;
  ws->w550 = -1;
  ws->w590 = -1;
  ws->w670 = -1;
  ws->w700 = -1;    
  ws->wPAR_top = -1;

  for (w=0; w<ws->num_wave; w++){
    if (ws->wave[w] > 400.0 && ws->wPAR_bot < 0){
      ws->wPAR_bot = w;
      eco_write_setup(e,"Bottom of PAR range is centred on %4.2f nm with edge %4.2f nm \n",ws->wave[w],ws->bandedge[w]);
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
    if (ws->wave[w] == 700.0)
      ws->w700 = w;
    if (ws->wave[w] > 700.0 && ws->wPAR_top < 0){
      ws->wPAR_top = w-1;
      eco_write_setup(e,"Top of PAR range is centred on %4.2f nm with edge %4.2f nm \n",ws->wave[w-1],ws->bandedge[w]);
    }
  }

  // if no wavelength > 700 nm then set top wave length to ws->num_wave -1

  if (ws->wPAR_top == -1)
    ws->wPAR_top = ws->num_wave-2;

  ws->PARwidth = (ws->bandedge[ws->wPAR_top+1] - ws->bandedge[ws->wPAR_bot]);

  if ((ws->w440 == -1) && (ws->at_440_i > -1))
    e->quitfn("at_440 is in tracer list, but 440 nm is not in the optical grid.");
  if ((ws->w490 == -1) && (ws->Kd_490_i > -1))
    e->quitfn("Kd_490 is in tracer list, but 490 nm is not in the optical grid.");
  if ((ws->w550 == -1) && (ws->bt_550_i > -1))
    e->quitfn("bt_550 is in tracer list, but 550 nm is not in the optical grid.");
  if ((ws->w590 == -1) && (ws->bb_590_i > -1))
    e->quitfn("bb_590 is in tracer list, but 590 nm is not in the optical grid.");
  if ((ws->w700 == -1) && (ws->bb_700_i > -1))
    e->quitfn("bb_700 is in tracer list, but 700 nm is not in the optical grid.");

  /* load variables from netcdf library */

  /*
   * Read an external spectral library
   */

  int ncid,ncid_sat,ret,varid;
  FILE *fp = NULL;

  char *source_file = "csiro_optical_parameters_library.nc";
  char units[2];
  char description[100];
  char symbol[17];
  char sourcedate[30];
  double bphy;

  //char *satellite_source_file = "csiro_spectral_response_library.nc";
  // ncw_open(satellite_source_file, NC_NOWRITE, &ncid_sat);
  
  fp = fopen(source_file, "r");
  
  if (fp == NULL){
    e->quitfn("In light_spectral_col: cannot read %s \n",source_file);
  }

  if (fp != NULL){
    fclose(fp);
    eco_write_setup(e,"Opening spectral library: %s \n",source_file);
    ncw_open(source_file, NC_NOWRITE, &ncid);
    ncw_get_att_text(source_file,ncid,NC_GLOBAL,"Creation date",sourcedate);

    // Create an optical setup file with spectrally-resolved parameters on the model wavelengths.

    int ncid1;
    int waveid,varid1;
    int dims[2];

    nc_create("optical_setup.nc",NC_CLOBBER, &ncid1);

    write_text_att(ncid1, NC_GLOBAL, "title", "CSIRO Environmental Modelling Suite (EMS) Optical setup file");
    write_text_att(ncid1, NC_GLOBAL, "description", "Optical grid and optical parameter values on the optical grid");
    write_date_created(ncid1);

    /* Source file and date created */

    write_text_att(ncid1, NC_GLOBAL, "default_source_file",source_file);
    nc_put_att_text(ncid1, NC_GLOBAL, "source_file_date_created",strlen(sourcedate),sourcedate);
  
    ncw_def_dim("optical_setup.nc",ncid1,"wave_centre",ws->num_wave, &waveid);
    dims[0] = waveid;
    nc_def_var(ncid1,"wave_centre",NC_DOUBLE,1,dims,&varid1);
    nc_put_att_text(ncid1, varid1, "description", 18,"Centre of waveband");
    nc_put_att_text(ncid1, varid1, "units", 2,"nm");
    nc_put_att_text(ncid1, varid1, "symbol", 8,"\\lambda");
    nc_put_att_text(ncid1, varid1, "puv_uom",52,"http://vocab.nerc.ac.uk/collection/P06/current/UXNM/");
    nc_def_var(ncid1,"bandwidth",NC_DOUBLE,1,dims,&varid1);
    nc_put_att_text(ncid1, varid1, "description", 17,"Width of waveband");
    nc_put_att_text(ncid1, varid1, "units", 2,"nm");
    nc_put_att_text(ncid1, varid1, "symbol", 16,"\\Delta \\lambda");
    nc_put_att_text(ncid1, varid1, "puv_uom",52,"http://vocab.nerc.ac.uk/collection/P06/current/UXNM/");

    int waveedgeid;
    ncw_def_dim("optical_setup.nc",ncid1,"wave_edge",ws->num_wave+1, &waveedgeid);
    dims[0] = waveedgeid;
    nc_def_var(ncid1,"wave_edge",NC_DOUBLE,1,dims,&varid1);
    nc_put_att_text(ncid1, varid1, "description", 17,"Edges of waveband");
    nc_put_att_text(ncid1, varid1, "units", 2,"nm");
    nc_put_att_text(ncid1, varid1, "symbol", 8,"\\lambda");
    nc_put_att_text(ncid1, varid1, "puv_uom",52,"http://vocab.nerc.ac.uk/collection/P06/current/UXNM/");
    
    int wavesensorid;
    int dims22[2];
    
    ncw_def_dim("optical_setup.nc",ncid1,"wave_sensor",ws->num_sensor_waves, &wavesensorid);
    dims22[0] = wavesensorid;
    nc_def_var(ncid1,"wave_sensor",NC_DOUBLE,1,dims22,&varid1);
    nc_put_att_text(ncid1, varid1, "description", 30,"Wavelength of interpolated Rrs");
    nc_put_att_text(ncid1, varid1, "units", 2,"nm");
    nc_put_att_text(ncid1, varid1, "symbol", 8,"\\lambda");
    nc_put_att_text(ncid1, varid1, "puv_uom",52,"http://vocab.nerc.ac.uk/collection/P06/current/UXNM/");
    
    int dim_dummy[2];

    nc_def_var(ncid1,"Benthic_canopy_order",NC_CHAR,0,dim_dummy,&varid1);
    nc_put_att_text(ncid1, varid1, "description", 34,"Vertical order of plants in canopy");
    nc_put_att_text(ncid1, varid1, "orientation", 35,"First is higher in the water column");

    nc_def_var(ncid1,"PARbot",NC_DOUBLE,0,dim_dummy,&varid1);
    nc_put_att_text(ncid1, varid1, "description", 28,"Lower wavelength edge of PAR");
    nc_put_att_text(ncid1, varid1, "units", 2,"nm");
    nc_put_att_text(ncid1, varid1, "puv_uom",52,"http://vocab.nerc.ac.uk/collection/P06/current/UXNM/");

    nc_def_var(ncid1,"PARtop",NC_DOUBLE,0,dim_dummy,&varid1);
    nc_put_att_text(ncid1, varid1, "description", 28,"Upper wavelength edge of PAR");
    nc_put_att_text(ncid1, varid1, "units", 2,"nm");
    nc_put_att_text(ncid1, varid1, "puv_uom",52,"http://vocab.nerc.ac.uk/collection/P06/current/UXNM/");

    nc_def_var(ncid1,"PARwidth",NC_DOUBLE,0,dim_dummy,&varid1);
    nc_put_att_text(ncid1, varid1, "description", 18,"Width of PAR range");
    nc_put_att_text(ncid1, varid1, "units", 2,"nm");

    nc_def_var(ncid1,"bbcell1",NC_DOUBLE,0,dim_dummy,&varid1);
    nc_put_att_text(ncid1, varid1, "description", 23,"Cell backscatter coef 1");
    nc_put_att_text(ncid1, varid1, "units", 2,"??");
    
    nc_def_var(ncid1,"bbcell2",NC_DOUBLE,0,dim_dummy,&varid1);
    nc_put_att_text(ncid1, varid1, "description", 26,"Cell backscatter exponent");
    nc_put_att_text(ncid1, varid1, "units", 1,"-");

    nc_def_var(ncid1,"bphy",NC_DOUBLE,0,dim_dummy,&varid1);
    nc_put_att_text(ncid1, varid1, "description", 45,"Chl-specific scattering coef. for microalgae");
    nc_put_att_text(ncid1, varid1, "units", 20,"m-1 (mg Chl a m-3)-1");

    nc_def_var(ncid1,"NtoCHL",NC_DOUBLE,0,dim_dummy,&varid1);
    nc_put_att_text(ncid1, varid1, "description", 39,"N to Chl ratio for backscattering calc.");
    nc_put_att_text(ncid1, varid1, "units", 5,"g g-1");
    

    if (ws->PhyS_N_i > -1){
      nc_def_var(ncid1,"yC_s",NC_DOUBLE,1,dims,&varid1);
      nc_put_att_text(ncid1, varid1, "description", 50,"Small phytoplankton Chl-a mass specific absorbance");
      nc_put_att_text(ncid1, varid1, "units", 20,"m-1 (mg Chl-a m-3)-1");
      nc_def_var(ncid1,"PhyS_rad",NC_DOUBLE,0,dim_dummy,&varid1);
      nc_put_att_text(ncid1, varid1, "description", 26,"Small phytoplankton radius");
      nc_put_att_text(ncid1, varid1, "units", 1,"m");
    }
    if (ws->PhyL_N_i > -1){
      nc_def_var(ncid1,"yC_l",NC_DOUBLE,1,dims,&varid1);
      nc_put_att_text(ncid1, varid1, "description", 50,"Large phytoplankton Chl-a mass specific absorbance");
      nc_put_att_text(ncid1, varid1, "units", 20,"m-1 (mg Chl-a m-3)-1");
      nc_def_var(ncid1,"PhyL_rad",NC_DOUBLE,0,dim_dummy,&varid1);
      nc_put_att_text(ncid1, varid1, "description", 26,"Large phytoplankton radius");
      nc_put_att_text(ncid1, varid1, "units", 1,"m");
    }
    if (ws->Tricho_N_i > -1){
      nc_def_var(ncid1,"yC_Tricho",NC_DOUBLE,1,dims,&varid1);
      nc_put_att_text(ncid1, varid1, "description", 44,"Trichodesmium Chl-a mass specific absorbance");
      nc_put_att_text(ncid1, varid1, "units", 20,"m-1 (mg Chl-a m-3)-1");
      nc_def_var(ncid1,"Tricho_rad",NC_DOUBLE,0,dim_dummy,&varid1);
      nc_put_att_text(ncid1, varid1, "description", 20,"Trichodesmium radius");
      nc_put_att_text(ncid1, varid1, "units", 1,"m");
    }
    if (ws->PhyD_N_i > -1){
      nc_def_var(ncid1,"yC_D",NC_DOUBLE,1,dims,&varid1);
      nc_put_att_text(ncid1, varid1, "description", 45,"Dinoflagellate Chl-a mass specific absorbance");
      nc_put_att_text(ncid1, varid1, "units", 20,"m-1 (mg Chl-a m-3)-1");
      nc_def_var(ncid1,"PhyD_rad",NC_DOUBLE,0,dim_dummy,&varid1);
      nc_put_att_text(ncid1, varid1, "description", 21,"Dinoflagellate radius");
      nc_put_att_text(ncid1, varid1, "units", 1,"m");
    }
    if (ws->MPB_N_i > -1){
      nc_def_var(ncid1,"yC_MPB",NC_DOUBLE,1,dims,&varid1);
      nc_put_att_text(ncid1, varid1, "description", 48,"Microphytobenthos Chl-a mass specific absorbance");
      nc_put_att_text(ncid1, varid1, "units", 20,"m-1 (mg Chl-a m-3)-1");
      nc_def_var(ncid1,"MPB_rad",NC_DOUBLE,0,dim_dummy,&varid1);
      nc_put_att_text(ncid1, varid1, "description", 25,"Microphytobenthos radius");
      nc_put_att_text(ncid1, varid1, "units", 1,"m");
    }
    if (ws->CS_N_i > -1){
      nc_def_var(ncid1,"yC_Symbiodinium",NC_DOUBLE,1,dims,&varid1);
      nc_put_att_text(ncid1, varid1, "description", 43,"Symbiodinium Chl-a mass specific absorbance");
      nc_put_att_text(ncid1, varid1, "units", 20,"m-1 (mg Chl-a m-3)-1");
      nc_def_var(ncid1,"CS_rad",NC_DOUBLE,0,dim_dummy,&varid1);
      nc_put_att_text(ncid1, varid1, "description", 19,"Symbiodinium radius");
      nc_put_att_text(ncid1, varid1, "units", 1,"m");
    }
    
    nc_close(ncid1);

     /* initialise moonlighting spectral reflectance parameters */

    usgs_moon_reflectance(e,ws);
    
    /* Interpolate clear-water properties so they are available for function colour_in_seawater */

    ws->aw = d_alloc_1d(ws->num_wave);
    optdataspectralread(e,ws,source_file,ncid,"aw","aw_wave","aw",ws->aw);
    
    ws->bw = d_alloc_1d(ws->num_wave);
    optdataspectralread(e,ws,source_file,ncid,"btw","btw_wave","bw",ws->bw);

    
    /*************************************************************/
    /* Now call obtain tracer list-specified optical properties. */

      // might have multiple variable use the same optical parameter in the library. So variable name in source != name in optical_setup.nc

      char trname[MAXSTRLEN] = "";
      char fname[MAXSTRLEN];
      char absorp_name[MAXSTRLEN];
      char scatter_name[MAXSTRLEN];
      char backscat_name[MAXSTRLEN];
      char benreflt_name[MAXSTRLEN];
      char specresp3d_name[MAXSTRLEN];
      char setup_absorp_name[MAXSTRLEN];
      char setup_scatter_name[MAXSTRLEN];
      char setup_backscat_name[MAXSTRLEN];
      char setup_benreflt_name[MAXSTRLEN];
      char setup_specresp3d_name[MAXSTRLEN];
      
      int tmp11;
      int aa = -1;
      int bbb = -1;

      for (i = 0; i < e->ntr; i++) {

	// Assume that any optical property in WC has all of at,bt,bb,bottom reflection.

	strcpy(trname,  e->tracernames[i]);

	if (einterface_is_optical(e->model,trname)) {

	  sprintf(absorp_name,"");
	  sprintf(scatter_name,"");
	  sprintf(backscat_name,"");
	  sprintf(benreflt_name,"");
	  sprintf(specresp3d_name,"");
	  sprintf(setup_absorp_name,"");
	  sprintf(setup_scatter_name,"");
	  sprintf(setup_backscat_name,"");
	  sprintf(setup_benreflt_name,"");
	  sprintf(setup_specresp3d_name,"");
	  
	  tmp11 = einterface_get_optical_file(e->model, trname, fname);
	  tmp11 = einterface_get_absorp_name(e->model, trname, absorp_name);
	  tmp11 = einterface_get_scatter_name(e->model, trname, scatter_name);
	  tmp11 = einterface_get_backscat_name(e->model, trname, backscat_name);
	  tmp11 = einterface_get_benreflt_name(e->model, trname, benreflt_name);
	  tmp11 = einterface_get_specresp3d_name(e->model, trname, specresp3d_name);

	  eco_write_setup(e,"Index %d TRACERNAME: (%s) fname (%s) ab name (%s) scatter name (%s) backscat name (%s) benreflt name (%s) specresp3d name (%s) \n",i,trname,fname,absorp_name,scatter_name,backscat_name,benreflt_name,specresp3d_name);
	  	  
	  sprintf(setup_absorp_name, "ap_%s",trname);
	  sprintf(setup_scatter_name, "bp_%s",trname);
	  sprintf(setup_backscat_name, "B_%s",trname);
	  sprintf(setup_benreflt_name, "dR_%s",trname); // dimensionless reflectance
	  sprintf(setup_specresp3d_name, "sR_%s",trname); // normalised sensor response
	  
	  int ncid9;

	  ncw_open(fname, NC_NOWRITE, &ncid9);
	  if (strcmp(absorp_name,"")){ // strcmp returns 0 for the equal!
	    aa = aa + 1;
	    optdataspectralread(e,ws,fname,ncid9,absorp_name," ",setup_absorp_name,ws->at_star[aa]);
	    ws->ind_at[aa] = i;
	    ws->ind_bt[aa] = i;
	    ws->ind_bb[aa] = i;
	    ws->ind_br[aa] = i;
	  }
	  if (strcmp(scatter_name,""))
	    optdataspectralread(e,ws,fname,ncid9,scatter_name," ",setup_scatter_name,ws->bt_star[aa]);
	  if (strcmp(backscat_name,""))
	    optdataspectralread(e,ws,fname,ncid9,backscat_name," ",setup_backscat_name,ws->bb_star[aa]);
	  if (strcmp(benreflt_name,""))
	    optdataspectralread(e,ws,fname,ncid9,benreflt_name," ",setup_benreflt_name,ws->br_star[aa]);
	
	  colour_in_seawater(e,ws,trname);
	  colour_of_bottom_type(e,ws,trname);

	  if (strcmp(specresp3d_name,"")){
	    bbb = bbb + 1;
	    optdataspectralread_fine(e,ws,fname,ncid9,specresp3d_name," ",setup_specresp3d_name,ws->sp3_star[bbb]);
	    ws->ind_sp3[bbb] = i;
	  }
	  ncw_close(fname,ncid9);
	 }
      }

    // Now do epis 

      char epiname[MAXSTRLEN] = "";
      char fname1[MAXSTRLEN];
      char abtance_name[MAXSTRLEN];
      char refltce_name[MAXSTRLEN];
      char trnmiss_name[MAXSTRLEN]; 
      char specresp2d_name[MAXSTRLEN];
      char setup_abtance_name[MAXSTRLEN];
      char setup_refltce_name[MAXSTRLEN];
      char setup_trnmiss_name[MAXSTRLEN];
      char setup_specresp2d_name[MAXSTRLEN];

      aa = -1;
      bbb = -1;
      char tmmp[MAXSTRLEN] = "";
      int is_plant;

      eco_write_setup(e,"Tracer list has specifed the benthic plants in vertical order: \n");

      string4netcdf[0] = '\0';

      for (i = 0; i < e->nepi; i++) {

	strcpy(epiname,  e->epinames[i]);
	is_plant = 0;

	if (einterface_is_optical2d(e->model,epiname)) {

	  char tempstr1[strlen(epiname)-2];
	  strncpy(tempstr1,epiname,strlen(epiname)-2); // remove the "_N"
	  tempstr1[strlen(epiname)-2] = '\0';

	  if ((tempstr1[0] == 'M') && (tempstr1[1] == 'A')){
	    ws->numMA += 1;
	    is_plant = 1;
	  }
	  
	  if ((tempstr1[0] == 'S') && (tempstr1[1] == 'G')){
	    ws->numSG += 1;
	    is_plant = 1;
	  }
	  
	  if (is_plant){

	  aa = aa + 1;

	  strcpy(tmmp,string4netcdf);

	  sprintf(string4netcdf, "%s %s",tmmp,e->epinames[i]);
	   
	  tmp11 = einterface_get_optical2d_file(e->model, epiname, fname1);

	  tmp11 = einterface_get_abtance_name(e->model, epiname, abtance_name);
	  tmp11 = einterface_get_refltce_name(e->model, epiname, refltce_name);
	  tmp11 = einterface_get_trnmiss_name(e->model, epiname, trnmiss_name);

	  eco_write_setup(e,"\tIndex %d EPINAME: %s fname %s ab name %s refltce name %s trnmiss name %s \n",i,epiname,fname1,abtance_name,refltce_name,trnmiss_name);

	  sprintf(setup_abtance_name, "A_%s",epiname);
	  sprintf(setup_refltce_name, "R_%s",epiname);
	  sprintf(setup_trnmiss_name, "T_%s",epiname);
	  
	  // put interpolated into array "values"
	  
	  int ncid9;
	  
	  ncw_open(fname1, NC_NOWRITE, &ncid9);

	  optdataspectralread(e,ws,fname1,ncid9,abtance_name," ",setup_abtance_name,ws->A_star[aa]);
	  optdataspectralread(e,ws,fname1,ncid9,refltce_name," ",setup_refltce_name,ws->R_star[aa]);
	  optdataspectralread(e,ws,fname1,ncid9,trnmiss_name," ",setup_trnmiss_name,ws->T_star[aa]);

	  colour_of_bottom_type(e,ws,epiname);

	  // Now add nitrogen-specific leaf area of seagrass.
	  
	  int varid22;int varid23;
	  int ncid99; int dim_dummy22[2];
	  char leafden_name[MAXSTRLEN];
	  char orient_name[MAXSTRLEN];

	  ncw_open("optical_setup.nc",NC_WRITE, &ncid99);
	  ncredef(ncid99);

	  sprintf(leafden_name, "%sleafden",tempstr1);
	  sprintf(orient_name, "%sorient",tempstr1);
	  
	  nc_def_var(ncid99,leafden_name,NC_DOUBLE,0,dim_dummy22,&varid22);
	  nc_put_att_text(ncid99, varid22, "description", 27,"Nitrogen-specific leaf area");
	  nc_put_att_text(ncid99, varid22, "units", 6,"g-1 m2");
	  
	  nc_def_var(ncid99,orient_name,NC_DOUBLE,0,dim_dummy22,&varid23);
	  nc_put_att_text(ncid99, varid23, "description", 16,"Leaf orientation");
	  nc_put_att_text(ncid99, varid23, "units", 1,"-");

	  ncw_enddef("optical_setup.nc", ncid99);

	  ncw_put_var_double("optical_setup.nc",ncid99,varid22,&ws->XXXleafden[aa]);
	  ncw_put_var_double("optical_setup.nc",ncid99,varid23,&ws->XXXorient[aa]);

	  nc_close(ncid99);
	  
	  ncw_close(fname1,ncid9);

	  ws->ind_A[aa] = i;
	  ws->ind_R[aa] = i;
	  ws->ind_T[aa] = i;
	  }else{ // no is_plant

	    bbb = bbb + 1;
          
	    tmp11 = einterface_get_specresp2d_name(e->model, epiname, specresp2d_name);

	    tmp11 = einterface_get_optical2d_file(e->model, epiname, fname1);

	    sprintf(setup_specresp2d_name, "s_%s",epiname);

	    eco_write_setup(e,"\tIndex %d bbb %d EPINAME: %s fname %s specresp2d name %s setup_name %s \n",i,bbb,epiname,fname1,specresp2d_name,setup_specresp2d_name);

	    int ncid9;
	    ncw_open(fname1, NC_NOWRITE, &ncid9);
	    optdataspectralread_fine(e,ws,fname1,ncid9,specresp2d_name," ",setup_specresp2d_name,ws->sp2_star[bbb]);
	    ncw_close(fname,ncid9);
	    
	    ws->ind_sp2[bbb] = i;

	}
      }
      eco_write_setup(e,"This is also the order in the canopy, ordered top to bottom. \n");
    }

      /* Add coral skeleton parameters to optical_setup.nc */

      int ncid88; int varid88, varid89; int dim_dummy88[2];
      
      ncw_open("optical_setup.nc",NC_WRITE, &ncid88);
      ncredef(ncid88);

      nc_def_var(ncid88,"CHpolypden",NC_DOUBLE,0,dim_dummy88,&varid88);
      nc_put_att_text(ncid88, varid88, "description", 29,"Nitrogen-specific polyp area");
      nc_put_att_text(ncid88, varid88, "units", 6,"g-1 m2");

      nc_def_var(ncid88,"CHarea",NC_DOUBLE,0,dim_dummy88,&varid89);
      nc_put_att_text(ncid88, varid89, "description", 6,"CHarea");
      nc_put_att_text(ncid88, varid89, "units", 1,"-");
      
      ncw_enddef("optical_setup.nc", ncid88);
      
      ncw_put_var_double("optical_setup.nc",ncid88,varid88,&ws->CHpolypden);
      ncw_put_var_double("optical_setup.nc",ncid88,varid89,&ws->CHarea);

      nc_close(ncid88);

    /* End of tracer-list specified optically-active tracers. 
   /***********************************************************************************************/

    /* Conversion of SWR radiation to spectrally-resolved radiation - if landa in netcdf file */
    
    size_t len2;
    size_t ncstart2[2] = {0,0}; size_t nccount2[2] = {2,1};
    int dimid2;

    varid = ncw_var_id(ncid,"landa");
    if (varid > -1){
      dimid2 = ncw_dim_id(ncid, "landa_wave");
      ncw_inq_dimlen(source_file,ncid,dimid2,&len2);
      double nc_landa[2*len2];
      ncw_get_att_text(source_file,ncid,varid,"Units",units);
      ncw_get_att_text(source_file,ncid,varid,"Symbol",symbol);
      ncw_get_att_text(source_file,ncid,varid,"Description",description);
      
      nccount2[1] = len2;
      ncw_get_vara_double(source_file,ncid, varid, ncstart2,nccount2,nc_landa);
      
      eco_write_setup(e,"Read %s: %s [%s]\n",description,symbol,units);
      eco_write_setup(e,"wave [nm] landa \n");

      /* now interpolate onto model wavelengths */

      double t5[len2];
      double t6[len2];
      
      for (w=0; w<len2; w++){
	eco_write_setup(e,"w %d, %4.2f \t %e \n",w,nc_landa[w],nc_landa[w+len2]);
	t5[w] = nc_landa[w];
	t6[w] = nc_landa[w+len2];
      }
      
      double landa[ws->num_wave];
      double landa_edge[ws->num_wave+1];
     
      interp1d(t5, t6, len2, ws->wave,landa, ws->num_wave);

      /* interpolate landa onto edges */

      interp1d(t5, t6, len2, ws->bandedge, &landa_edge[0], ws->num_wave+1);

      /* 
       * Calculate energy within a wave band
       */
      double landa_out[ws->num_wave];

      double sumlanda_out = 0.0;
      for (w=0; w<ws->num_wave-1; w++){
	landa_out[w] = landa[w] * ws->bandwidth[w];
	sumlanda_out += landa_out[w]; 
      }

      /* Put the rest of the light in the top band */

      landa_out[ws->num_wave-1] = 1.0 - sumlanda_out;

      // Allocate and fill in
      ws->landa = d_alloc_1d(ws->num_wave);
      for (w=0; w<ws->num_wave; w++){
	ws->landa[w] = landa_out[w];
      }
      for (w=0; w<ws->num_wave; w++){
	eco_write_setup(e,"w %d, %4.2f \t %4.2f \t %e \n",w,ws->wave[w],ws->bandwidth[w],ws->landa[w]);
      }

      // Now write landa to the output netcdf file:

      int dimid,dims_l[2];
      size_t len_l;
      size_t ncstart_l[2] = {0,0}; size_t nccount_l[2] = {2,1};
      
      ncw_open("optical_setup.nc",NC_WRITE, &ncid1);
      dimid = ncw_dim_id(ncid1,"wave_centre");
      dims_l[0] = dimid;
      ncredef(ncid1);

      ncw_def_var("optical_setup.nc",ncid1,"landa",NC_DOUBLE,1,dims_l,&varid);
      varid = ncw_var_id(ncid1,"landa");
      nc_put_att_text(ncid1, varid, "description", strlen(description),description);
      nc_put_att_text(ncid1, varid, "symbol", strlen(symbol),symbol);
      nc_put_att_text(ncid1, varid, "units", strlen(units),units);
      ncw_enddef("optical_setup.nc", ncid1);
      ncw_put_var_double("optical_setup.nc",ncid1,varid,ws->landa);

      // const size_t nc_st = 0;
      // const size_t nc_co = strlen(string4netcdf);
      varid = ncw_var_id(ncid1,"Benthic_canopy_order");
      //printf("output to netcdf canopy %s \n",string4netcdf);
      //nc_put_vara_uchar(ncid1,varid,nc_st,nc_co,string4netcdf);
      nc_close(ncid1);
    }

  /* Constant optical properties */

    optdatascalarread(e,ws,source_file,ncid,"g0",&ws->g0);
    optdatascalarread(e,ws,source_file,ncid,"g1",&ws->g1);

    /* Spectrally-independent variables - Kirk (1991) LO 36: 455*/ 

    optdatascalarread(e,ws,source_file,ncid,"gi",&ws->gi);
    optdatascalarread(e,ws,source_file,ncid,"gii",&ws->gii);

    /* Spectrally-dependent variables */
    
    ws->bls = d_alloc_1d(ws->num_wave);
    optdataspectralread(e,ws,source_file,ncid,"bls","bls_wavelength","bls",ws->bls);

    ws->fls = d_alloc_1d(ws->num_wave);
    optdataspectralread(e,ws,source_file,ncid,"fls","fls_wavelength","fls",ws->fls);

    // useful to have sum under the curve 

    double bls_sum = 0.0;
    double fls_sum = 0.0;
    for (w=0; w < ws->num_wave; w++){
      bls_sum += ws->bls[w];
      fls_sum += ws->fls[w];
    }
    ws->bls_sum = bls_sum;
    ws->fls_sum = fls_sum;

    /* Pigments */

    if (ws->CS_Xp_i > -1){
      ws->yC_diatoxanthin = d_alloc_1d(ws->num_wave);
      optdataspectralread(e,ws,source_file,ncid,"Diatoxanthin","pigments_wave","yC_diatoxanthin",ws->yC_diatoxanthin);
      ws->yC_diadinoxanthin = d_alloc_1d(ws->num_wave);
      optdataspectralread(e,ws,source_file,ncid,"Diadinoxanthin","pigments_wave","yC_diadinoxanthin",ws->yC_diadinoxanthin);
    }
    
    /* Benthic plants and sediment surface */
    
    ws->rhoskel = d_alloc_1d(ws->num_wave);
    optdataspectralread(e,ws,source_file,ncid,"rho_coral","rho_coral_wave","rho_coral",ws->rhoskel);

    /* Suspended macroalgae */

    if (ws->SMA_N_i > -1){
      ws->SMA_absorption = d_alloc_1d(ws->num_wave);
      int ncidSMA;
      ncw_open("csiro_benthic_plant_reflectance_library.nc", NC_NOWRITE, &ncidSMA);
      optdataspectralread(e,ws,"csiro_benthic_plant_reflectance_library.nc",ncidSMA,"AZostera_tasmanica","DurakoP_wave","AZostera_tasmanica",ws->SMA_absorption);
      nc_close(ncidSMA);
    }

    /* Mass-specific optical properties of inorganic particles */

    ws->nbt_p_555 = d_alloc_1d(ws->num_wave);
    optdataspectralread(e,ws,source_file,ncid,"nbt_p_555","nbt_p_555_wave","nbt_p_555",ws->nbt_p_555);
    
    ncw_close(source_file,ncid);
    // ncw_close(satellite_source_file,ncid_sat);
  }else{
    eco_write_setup(e,"No spectral library file read \n");
    e->quitfn("Could not find the netcdf spectral library \n");
  }

  /* Find nearest wavelength for reflectance wavelengths specifed in tracer list */

  if (!e->pre_build){
    
    ecology_find_rsr_waves(e);
    ecology_find_ed_waves(e);
    ws->num_rrs_waves = e->num_rsr_waves;
    ws->rrs_wave      = e->rsr_waves;

    int w2;

    char buf[MAXSTRLEN];
    ws->RRS_i = NULL;
    if (ws->num_rrs_waves) 
      ws->RRS_i = i_alloc_1d(ws->num_rrs_waves);
    for (w2=0; w2<ws->num_rrs_waves; w2++) {
      int i_tmp;
      sprintf(buf, "R_%d", (int)ws->rrs_wave[w2]);
      i_tmp = e->try_index(epis, buf, e);
      if (i_tmp > -1){
	ws->RRS_i[w2] = i_tmp;
      }
      else
	e->quitfn("light_spectral_col: Could not find %s\n", buf);
    }

    // now find the index of the wave_centre less than RRS value

    int found;
    ws->wXXX_i = NULL;
    if (ws->num_rrs_waves) 
      ws->wXXX_i = i_alloc_1d(ws->num_rrs_waves);
    for (w2=0; w2<ws->num_rrs_waves; w2++) {
      found = 0;
      ws->wXXX_i[w2] = ws->num_wave-1;
      for (w=0;w<ws->num_wave;w++){
	if ((ws->wave[w] >= ws->rrs_wave[w2] ) && (found == 0)){
	  found = 1;
	  ws->wXXX_i[w2] = w;
	}
      }
      eco_write_setup(e,"\nRemote-sensing reflectance: Choosen wave_centres %f and %f nm for interpolation onto %f nm \n",ws->wave[ws->wXXX_i[w2]-1],ws->wave[ws->wXXX_i[w2]],ws->rrs_wave[w2]);
    }

    ws->num_ed_waves = e->num_ed_waves;
    ws->ed_wave      = e->ed_waves;

    ws->Ed_i = NULL;
    if (ws->num_ed_waves) 
      ws->Ed_i = i_alloc_1d(ws->num_ed_waves);
    for (w2=0; w2<ws->num_ed_waves; w2++) {
      int i_tmp;
      sprintf(buf, "Ed_%d", (int)ws->ed_wave[w2]);
      i_tmp = e->try_index(tracers, buf, e);
      if (i_tmp > -1){
	ws->Ed_i[w2] = i_tmp;
    }  
      else
	e->quitfn("light_spectral_col: Could not find %s\n", buf);
    }

    // now find the index of the wave_centre less than Ed value
    
    ws->wYYY_i = NULL;
    if (ws->num_ed_waves) 
      ws->wYYY_i = i_alloc_1d(ws->num_ed_waves);
    for (w2=0; w2<ws->num_ed_waves; w2++) {
      found = 0;
      ws->wYYY_i[w2] = ws->num_wave-1;
      for (w=0;w<ws->num_wave;w++){
	if ((ws->wave[w] >= ws->ed_wave[w2] ) && (found == 0)){
	  found = 1;
	  ws->wYYY_i[w2] = w;
	}
      }
      eco_write_setup(e,"\nDownwelling irradiance: Choosen wave_centres %f and %f nm for interpolation onto %f nm \n",ws->wave[ws->wYYY_i[w2]-1],ws->wave[ws->wYYY_i[w2]],ws->ed_wave[w2]);
    }
    
    ws->yC_s = d_alloc_1d(ws->num_wave);
    ws->yC_l = d_alloc_1d(ws->num_wave);
    ws->yC_Tricho = d_alloc_1d(ws->num_wave);
    ws->yC_MPB = d_alloc_1d(ws->num_wave);
    ws->yC_D = d_alloc_1d(ws->num_wave);
    ws->yC_Symbiodinium = d_alloc_1d(ws->num_wave);
    
    eco_write_setup(e,"\nMass-specific absorption coefficients used for microalgae in water column \n");
    
    chl_specific_aggregate_pigment_absorption(e,ws);
    
    double var1[ws->num_wave];
    for (w=0; w<ws->num_wave; w++){
      var1[w] = ws->wave[w];
    }
    
    nc_open("optical_setup.nc",NC_WRITE, &ncid);

    double temp_var;

    temp_var = ws->bphy;
    varid = ncw_var_id(ncid,"bphy");
    ncw_put_var_double("optical_setup.nc",ncid,varid,&temp_var);

    temp_var = ws->NtoCHL;
    varid = ncw_var_id(ncid,"NtoCHL");
    ncw_put_var_double("optical_setup.nc",ncid,varid,&temp_var);

    temp_var = ws->bbcell1;
    varid = ncw_var_id(ncid,"bbcell1");
    ncw_put_var_double("optical_setup.nc",ncid,varid,&temp_var);

    temp_var = ws->bbcell2;
    varid = ncw_var_id(ncid,"bbcell2");
    ncw_put_var_double("optical_setup.nc",ncid,varid,&temp_var);
    
    temp_var = ws->bandedge[ws->wPAR_bot];
    varid = ncw_var_id(ncid,"PARbot");
    ncw_put_var_double("optical_setup.nc",ncid,varid,&temp_var);

    temp_var = ws->bandedge[ws->wPAR_top+1];
    varid = ncw_var_id(ncid,"PARtop");
    ncw_put_var_double("optical_setup.nc",ncid,varid,&temp_var);

    temp_var = ws->PARwidth;
    varid = ncw_var_id(ncid,"PARwidth");
    ncw_put_var_double("optical_setup.nc",ncid,varid,&temp_var);

    varid = ncw_var_id(ncid,"wave_centre");
    ncw_put_var_double("optical_setup.nc",ncid,varid,ws->wave);

    varid = ncw_var_id(ncid,"wave_sensor");
    ncw_put_var_double("optical_setup.nc",ncid,varid,ws->sensor_waves);

    varid = ncw_var_id(ncid,"wave_edge");
    ncw_put_var_double("optical_setup.nc",ncid,varid,ws->bandedge);
    
    varid = ncw_var_id(ncid,"bandwidth");
    ncw_put_var_double("optical_setup.nc",ncid,varid,ws->bandwidth);
    
    if (ws->PhyS_N_i > -1){
      varid = ncw_var_id(ncid,"yC_s");
      ncw_put_var_double("optical_setup.nc",ncid,varid,ws->yC_s);
      varid = ncw_var_id(ncid,"PhyS_rad");
      ncw_put_var_double("optical_setup.nc",ncid,varid,&ws->rad_s);
    }
    if (ws->PhyL_N_i > -1){
      varid = ncw_var_id(ncid,"yC_l");
      ncw_put_var_double("optical_setup.nc",ncid,varid,ws->yC_l);
      varid = ncw_var_id(ncid,"PhyL_rad");
      ncw_put_var_double("optical_setup.nc",ncid,varid,&ws->rad_l);
    }
    if (ws->PhyD_N_i > -1){
      varid = ncw_var_id(ncid,"yC_D");
      ncw_put_var_double("optical_setup.nc",ncid,varid,ws->yC_D);
      varid = ncw_var_id(ncid,"PhyD_rad");
      ncw_put_var_double("optical_setup.nc",ncid,varid,&ws->rad_PhyD);
    }
    if (ws->MPB_N_i > -1){
      varid = ncw_var_id(ncid,"yC_MPB");
      ncw_put_var_double("optical_setup.nc",ncid,varid,ws->yC_MPB);
      varid = ncw_var_id(ncid,"MPB_rad");
      ncw_put_var_double("optical_setup.nc",ncid,varid,&ws->rad_MPB);
    }
    if (ws->Tricho_N_i > -1){
      varid = ncw_var_id(ncid,"yC_Tricho");
      ncw_put_var_double("optical_setup.nc",ncid,varid,ws->yC_Tricho);
      varid = ncw_var_id(ncid,"Tricho_rad");
      ncw_put_var_double("optical_setup.nc",ncid,varid,&ws->rad_Tricho);
    }
    if (ws->CS_N_i > -1){
      varid = ncw_var_id(ncid,"yC_Symbiodinium");
      ncw_put_var_double("optical_setup.nc",ncid,varid,ws->yC_Symbiodinium);
      varid = ncw_var_id(ncid,"CS_rad");
      ncw_put_var_double("optical_setup.nc",ncid,varid,&ws->rad_CS);
    }
    if (ws->CS_Xp_i > -1){    // this is probably uneccessary - already done by optspectralread 
      varid = ncw_var_id(ncid,"yC_diatoxanthin");
      ncw_put_var_double("optical_setup.nc",ncid,varid,ws->yC_diatoxanthin);
      varid = ncw_var_id(ncid,"yC_diadinoxanthin");
      ncw_put_var_double("optical_setup.nc",ncid,varid,ws->yC_diadinoxanthin);
    }

    ncw_close("optical_setup.nc",ncid);
        
    // Create netcdf file for radiance output - can reference layer numbers after pre-build.
    
    ws->num_wc_layers = ginterface_getnumberwclayers(e->model);
    ws->num_sed_layers = ginterface_getnumbersedlayers(e->model);

    if (!e->pre_build){   // only do once
      // ws->column_out_name = '\0';
      char* output_path = ginterface_get_output_path();
      if (output_path){
	sprintf(ws->column_out_name,"%soptical_column_out.nc",output_path);
      }else{
	sprintf(ws->column_out_name,"optical_column_out.nc");
      }
      eco_write_setup(e,"Name and path for optical column output file is %s \n",ws->column_out_name);
	
      int waveid,zenid,aziid,recid,depthid,depthcentreid,wavesensorid;
      int dims3[2],dims4[2];

      static double air = 1.0e35;
      static double seafloor = -9999.0;

      // Create column netcdf if it doesn't already exist.
      //double *recordtime;
      int ems_restart = 0;
      ems_restart = nc_create(ws->column_out_name,NC_NOCLOBBER, &ncid); // this will create a compressed file: NC_NOCLOBBER|NC_NETCDF4
      eco_write_setup(e,"ems_restart flag: %d \n",ems_restart);
      //if (ems_restart != 0){
	// remove any record dimensions after e-t.
	// nc_inq_varid(ncid,"time",&varid);
	// nc_get_var_double(ncid,varid,recordtime);
      //}
      
      if (ems_restart == 0){
      write_text_att(ncid, NC_GLOBAL, "title", "CSIRO Environmental Modelling Suite (EMS) optical model output.");
      write_text_att(ncid, NC_GLOBAL, "description", "Spectrally-resolved optical properties in a model column.");
      write_text_att(ncid, NC_GLOBAL, "vertical grid", "Depth is relative to the moving surface, so f(space,time).");
      write_date_created(ncid);
      
      nc_def_var(ncid,"i_index",NC_INT,0,0,&varid);
      nc_def_var(ncid,"j_index",NC_INT,0,0,&varid);
     
      int dim_dummy2[2];

      nc_def_var(ncid,"latitude",NC_DOUBLE,0,dim_dummy2,&varid);
      nc_put_att_text(ncid, varid, "description", 8,"Latitude");
      nc_put_att_text(ncid, varid, "units", 13,"degrees_north");

      nc_def_var(ncid,"longitude",NC_DOUBLE,0,dim_dummy2,&varid);
      nc_put_att_text(ncid, varid, "description", 9,"Longitude");
      nc_put_att_text(ncid, varid, "units", 12,"degrees_east");
      
      recid = ncdimdef(ncid, "record", NC_UNLIMITED);
      
      dims3[0] = recid;
      nc_def_var(ncid,"t",NC_DOUBLE,1,dims3,&varid);
      nc_put_att_text(ncid, varid, "units", strlen(ws->outputtimeunits),ws->outputtimeunits);
      nc_put_att_text(ncid, varid, "long_name",4,"Time");
      nc_put_att_text(ncid, varid, "coordinate_type",4,"time");

      nc_def_var(ncid,"col_index",NC_INT,1,dims3,&varid);
      nc_put_att_text(ncid, varid, "description", 10,"EMS column");
      
      nc_def_var(ncid,"use_radtranx",NC_INT,1,dims3,&varid);
      nc_put_att_text(ncid, varid, "description", 35,"Solar radiation used by hydrolight");
      nc_put_att_text(ncid, varid, "Zero", 35,"EMS light used in hydrolight calcs");
      nc_put_att_text(ncid, varid, "One", 40, "RADTRANX light used in hydrolight calcs");

      nc_def_var(ncid,"ems_solar_zenith",NC_DOUBLE,1,dims3,&varid);
      nc_put_att_text(ncid, varid, "description", 23,"EMS solar zenith angle");
      nc_put_att_text(ncid, varid, "units", 7,"radians");
      nc_def_var(ncid,"ems_solar_azimuth",NC_DOUBLE,1,dims3,&varid);
      nc_put_att_text(ncid, varid, "description", 24,"EMS solar azimuth angle");
      nc_put_att_text(ncid, varid, "units", 7,"radians");

      nc_def_var(ncid,"hyd_solar_zenith",NC_DOUBLE,1,dims3,&varid);
      nc_put_att_text(ncid, varid, "description", 23,"HYD solar zenith angle");
      nc_put_att_text(ncid, varid, "units", 7,"radians");
      
      nc_def_var(ncid,"hyd_solar_azimuth",NC_DOUBLE,1,dims3,&varid);
      nc_put_att_text(ncid, varid, "description", 24,"HYD solar azimuth angle");
      nc_put_att_text(ncid, varid, "units", 7,"radians");

      nc_def_var(ncid,"eta",NC_DOUBLE,1,dims3,&varid);
      nc_put_att_text(ncid, varid, "description", 17,"Surface elevation");
      nc_put_att_text(ncid, varid, "units",1,"m");

      nc_def_var(ncid,"wind_speed",NC_DOUBLE,1,dims3,&varid);
      nc_put_att_text(ncid, varid, "description", 10,"Wind speed");
      nc_put_att_text(ncid, varid, "units",5,"m s-1");

      nc_def_var(ncid,"cloud_cover",NC_DOUBLE,1,dims3,&varid);
      nc_put_att_text(ncid, varid, "description", 20,"Total cloud fraction");
      nc_put_att_text(ncid, varid, "units",1,"-");

      nc_def_var(ncid,"mslp",NC_DOUBLE,1,dims3,&varid);
      nc_put_att_text(ncid, varid, "description", 24,"Mean sea level pressure");
      nc_put_att_text(ncid, varid, "units",2,"Pa");

      nc_def_var(ncid,"rh2m",NC_DOUBLE,1,dims3,&varid);
      nc_put_att_text(ncid, varid, "description", 17,"Relative humidity");
      nc_put_att_text(ncid, varid, "units",1,"%");

      if (ws->Sentinel_3B_Band3_i){
	nc_def_var(ncid,"Sentinel_3B_Band3",NC_DOUBLE,1,dims3,&varid);
	nc_put_att_text(ncid, varid, "description", 47,"Remote-sensing reflectance on Sentinel_3B_Band3");
	nc_put_att_text(ncid, varid, "units", 4,"sr-1");
      }
      if (ws->Sentinel_3B_Band4_i){
	nc_def_var(ncid,"Sentinel_3B_Band4",NC_DOUBLE,1,dims3,&varid);
	nc_put_att_text(ncid, varid, "description", 47,"Remote-sensing reflectance on Sentinel_3B_Band4");
	nc_put_att_text(ncid, varid, "units", 4,"sr-1");
      }
      if (ws->Sentinel_3B_Band5_i){
	nc_def_var(ncid,"Sentinel_3B_Band5",NC_DOUBLE,1,dims3,&varid);
	nc_put_att_text(ncid, varid, "description", 47,"Remote-sensing reflectance on Sentinel_3B_Band5");
	nc_put_att_text(ncid, varid, "units", 4,"sr-1");
      }
      if (ws->Sentinel_3B_Band6_i){
	nc_def_var(ncid,"Sentinel_3B_Band6",NC_DOUBLE,1,dims3,&varid);
	nc_put_att_text(ncid, varid, "description", 47,"Remote-sensing reflectance on Sentinel_3B_Band6");
	nc_put_att_text(ncid, varid, "units", 4,"sr-1");
      }
      if (ws->Sentinel_3B_Band8_i){
	nc_def_var(ncid,"Sentinel_3B_Band8",NC_DOUBLE,1,dims3,&varid);
	nc_put_att_text(ncid, varid, "description", 47,"Remote-sensing reflectance on Sentinel_3B_Band8");
	nc_put_att_text(ncid, varid, "units", 4,"sr-1");
      }
      if (ws->Sentinel_3B_Band11_i){
	nc_def_var(ncid,"Sentinel_3B_Band11",NC_DOUBLE,1,dims3,&varid);
	nc_put_att_text(ncid, varid, "description", 47,"Remote-sensing reflectance on Sentinel_3B_Band11");
	nc_put_att_text(ncid, varid, "units", 4,"sr-1");
      }
      
      if (ws->Moonlight_i > -1){
	nc_def_var(ncid,"ems_lunar_zenith",NC_DOUBLE,1,dims3,&varid);
	nc_put_att_text(ncid, varid, "description", 19,"Lunar zenith angle");
	nc_put_att_text(ncid, varid, "units", 7,"radians");
	nc_def_var(ncid,"ems_lunar_azimuth",NC_DOUBLE,1,dims3,&varid);
	nc_put_att_text(ncid, varid, "description", 20,"Lunar azimuth angle");
	nc_put_att_text(ncid, varid, "units", 7,"radians");
	nc_def_var(ncid,"moon_phase",NC_DOUBLE,1,dims3,&varid);
	nc_put_att_text(ncid, varid, "description", 19,"Lunar phase angle");
	nc_put_att_text(ncid, varid, "units", 7,"radians");
	nc_def_var(ncid,"moon_fulldisk",NC_DOUBLE,1,dims3,&varid);
	nc_put_att_text(ncid, varid, "description", 29,"Clear-sky fulldisk moonlight");
	nc_put_att_text(ncid, varid, "units", 5,"W m-2");
	nc_def_var(ncid,"ems_lunar_albedo",NC_DOUBLE,1,dims3,&varid);
	nc_put_att_text(ncid, varid, "description", 16,"EMS lunar albedo");
	nc_put_att_text(ncid, varid, "units", 1,"-");
	nc_def_var(ncid,"ems_lunar_airtrans",NC_DOUBLE,1,dims3,&varid);
	nc_put_att_text(ncid, varid, "description", 41,"EMS lunar atmospheric transmission loss");
	nc_put_att_text(ncid, varid, "units", 1,"-");
	nc_def_var(ncid,"moonlight",NC_DOUBLE,1,dims3,&varid);
	nc_put_att_text(ncid, varid, "description", 25,"Moonlight at sea surface");
	nc_put_att_text(ncid, varid, "units", 5,"W m-2");
      }
      
      ncw_def_dim(ws->column_out_name,ncid,"k_centre",ws->num_wc_layers, &depthcentreid);
      ncw_def_dim(ws->column_out_name,ncid,"k_grid",ws->num_wc_layers+1, &depthid);
      ncw_def_dim(ws->column_out_name,ncid,"wave_centre",ws->num_wave, &waveid);
      
      int dims6[2] = {recid,depthid};
      nc_def_var(ncid,"depth",NC_DOUBLE,2,dims6,&varid);
      nc_put_att_text(ncid, varid, "description", 33,"metres from sea surface of k_grid");
      
      int dims8[3] = {recid,waveid,depthid};
      nc_def_var(ncid,"E_d_ray",NC_DOUBLE,3,dims8,&varid);
      nc_put_att_text(ncid, varid, "description", 53,"Downwelling irradiance, calculated using ray-tracing");
      nc_put_att_text(ncid, varid, "units", 10,"W m-2 nm-1");
      nc_put_att_double(ncid, varid, "air_value", NC_DOUBLE, 1,&air);
      nc_put_att_double(ncid, varid, "seafloor_value", NC_DOUBLE, 1,&seafloor);
      
      int dims19[2] = {recid,waveid};
      nc_def_var(ncid,"Bottom_reflectance",NC_DOUBLE,2,dims19,&varid);
      nc_put_att_text(ncid, varid, "description", 18,"Bottom_reflectance");
      nc_put_att_text(ncid, varid, "units", 1,"-");
      
      nc_def_var(ncid,"Rrs_odw",NC_DOUBLE,2,dims19,&varid);
      nc_put_att_text(ncid, varid, "description", 64,"Remote-sensing reflectance, Lu/Ed, using optical-depth weighting");
      nc_put_att_text(ncid, varid, "units", 4,"sr-1");
      nc_put_att_text(ncid, varid, "puv_uom",52,"http://vocab.nerc.ac.uk/collection/P06/current/PSTR/");

      ncw_def_dim(ws->column_out_name,ncid,"wave_sensor",ws->num_sensor_waves, &wavesensorid);

      int dims22a[2];
      dims22a[0] = wavesensorid;

      nc_def_var(ncid,"wave_sensor",NC_DOUBLE,1,dims22a,&varid);
      nc_put_att_text(ncid, varid, "description", 30,"Wavelength of interpolated Rrs");
      nc_put_att_text(ncid, varid, "units", 2,"nm");
      nc_put_att_text(ncid, varid, "symbol", 8,"\\lambda");
      nc_put_att_text(ncid, varid, "puv_uom",52,"http://vocab.nerc.ac.uk/collection/P06/current/UXNM/");
      
      int dims19a[2] = {recid,wavesensorid};
      
      nc_def_var(ncid,"Rrs_fine",NC_DOUBLE,2,dims19a,&varid);
      nc_put_att_text(ncid, varid, "description", 77,"Remote-sensing reflectance, Lu/Ed, using optical-depth weighting on fine grid");
      nc_put_att_text(ncid, varid, "units", 4,"sr-1");
      nc_put_att_text(ncid, varid, "puv_uom",52,"http://vocab.nerc.ac.uk/collection/P06/current/PSTR/");

      #ifdef INCLUDE_BIOLUMINESCENCE
      nc_def_var(ncid,"Secchi_colour",NC_DOUBLE,2,dims19,&varid);
      nc_put_att_text(ncid, varid, "description", 21,"Coloured Secchi depth");
      nc_put_att_text(ncid, varid, "units", 1,"m");
      #endif
      
      nc_def_var(ncid,"n_stw",NC_DOUBLE,2,dims19,&varid);
      nc_put_att_text(ncid, varid, "description", 31,"Index of refraction at surface");
      nc_put_att_text(ncid, varid, "units", 1,"-");

      if (ws->Moonlight_i > -1){
	nc_def_var(ncid,"Moonlight_SR",NC_DOUBLE,2,dims19,&varid);
	nc_put_att_text(ncid, varid, "description", 53,"Downwelling moonlight spectrum with no atmosphere");
	nc_put_att_text(ncid, varid, "units", 10,"W m-2 nm-1");
      }

#ifdef INCLUDE_BIOLUMINESCENCE
      nc_def_var(ncid,"E_u_z_bio",NC_DOUBLE,3,dims8,&varid);
      nc_put_att_text(ncid, varid, "description", 35,"Upwelling bioluminescent irradiance");
      nc_put_att_text(ncid, varid, "units", 10,"W m-2 nm-1");
      nc_put_att_double(ncid, varid, "air_value", NC_DOUBLE, 1,&air);
      nc_put_att_double(ncid, varid, "seafloor_value", NC_DOUBLE, 1,&seafloor);
#endif
      
#ifdef INCLUDE_HYDROLIGHT
      nc_def_var(ncid,"Rrs_hyd",NC_DOUBLE,2,dims19,&varid);
      nc_put_att_text(ncid, varid, "description", 64,"Remote-sensing reflectance, Lu/Ed, using hydrolight model solut.");
      nc_put_att_text(ncid, varid, "units", 4,"sr-1");
      nc_put_att_text(ncid, varid, "puv_uom",52,"http://vocab.nerc.ac.uk/collection/P06/current/PSTR/");

      nc_def_var(ncid,"E_sky_hyd",NC_DOUBLE,2,dims19,&varid);
      nc_put_att_text(ncid, varid, "description", 35,"Downwelling irradiance from the sky");
      nc_put_att_text(ncid, varid, "units", 10,"W m-2 nm-1");
      nc_put_att_text(ncid, varid, "puv_uom",52,"http://vocab.nerc.ac.uk/collection/P06/current/UWNM/");

      nc_def_var(ncid,"E_ref_hyd",NC_DOUBLE,2,dims19,&varid);
      nc_put_att_text(ncid, varid, "description", 31,"Reflected irradiance to the sky");
      nc_put_att_text(ncid, varid, "units", 10,"W m-2 nm-1");
      nc_put_att_text(ncid, varid, "puv_uom",52,"http://vocab.nerc.ac.uk/collection/P06/current/UWNM/");

      nc_def_var(ncid,"E_wl_hyd",NC_DOUBLE,2,dims19,&varid);
      nc_put_att_text(ncid, varid, "description",24,"Water-leaving irradiance");
      nc_put_att_text(ncid, varid, "units", 10,"W m-2 nm-1");
      nc_put_att_text(ncid, varid, "puv_uom",52,"http://vocab.nerc.ac.uk/collection/P06/current/UWNM/");
      
      nc_def_var(ncid,"hydrolight_zenith",NC_DOUBLE,1,dims3,&varid);
      nc_def_var(ncid,"hydrolight_azimuth",NC_DOUBLE,1,dims3,&varid);
      
      nc_def_var(ncid,"E_d_hyd",NC_DOUBLE,3,dims8,&varid);
      nc_put_att_text(ncid, varid, "description", 32,"Downwelling radiance, hydrolight");
      nc_put_att_text(ncid, varid, "units", 10,"W m-2 nm-1");
      nc_put_att_double(ncid, varid, "air_value", NC_DOUBLE, 1,&air);
      nc_put_att_double(ncid, varid, "seafloor_value", NC_DOUBLE, 1,&seafloor);
      
      nc_def_var(ncid,"E_u_hyd",NC_DOUBLE,3,dims8,&varid);
      nc_put_att_text(ncid, varid, "description", 30,"Upwelling radiance, hydrolight");
      nc_put_att_text(ncid, varid, "units", 10,"W m-2 nm-1");
      nc_put_att_double(ncid, varid, "air_value", NC_DOUBLE, 1,&air);
      nc_put_att_double(ncid, varid, "seafloor_value", NC_DOUBLE, 1,&seafloor);
#endif
      
      int dims9[3] = {recid,waveid,depthcentreid};

      nc_def_var(ncid,"E_o_ray",NC_DOUBLE,3,dims9,&varid);
      nc_put_att_text(ncid, varid, "description", 47,"Scalar irradiance, calculated using ray-tracing");
      nc_put_att_text(ncid, varid, "units", 10,"W m-2 nm-1");
      nc_put_att_double(ncid, varid, "air_value", NC_DOUBLE, 1,&air);
      nc_put_att_double(ncid, varid, "seafloor_value", NC_DOUBLE, 1,&seafloor);
      nc_put_att_text(ncid, varid, "puv_uom",52,"http://vocab.nerc.ac.uk/collection/P06/current/UWNM/");
      
      nc_def_var(ncid,"at",NC_DOUBLE,3,dims9,&varid);
      nc_put_att_text(ncid, varid, "description", 16,"Total absorption");
      nc_put_att_text(ncid, varid, "units",3,"m-1");
      nc_put_att_double(ncid, varid, "air_value", NC_DOUBLE, 1,&air);
      nc_put_att_double(ncid, varid, "seafloor_value", NC_DOUBLE, 1,&seafloor);
      nc_put_att_text(ncid, varid, "puv_uom",52,"http://vocab.nerc.ac.uk/collection/P06/current/UPRM/");

      nc_def_var(ncid,"bt",NC_DOUBLE,3,dims9,&varid);
      nc_put_att_text(ncid, varid, "description", 16,"Total scattering");
      nc_put_att_text(ncid, varid, "units", 3,"m-1");
      nc_put_att_double(ncid, varid, "air_value", NC_DOUBLE, 1,&air);
      nc_put_att_double(ncid, varid, "seafloor_value", NC_DOUBLE, 1,&seafloor);
      nc_put_att_text(ncid, varid, "puv_uom",52,"http://vocab.nerc.ac.uk/collection/P06/current/UPRM/");
      
      nc_def_var(ncid,"bb",NC_DOUBLE,3,dims9,&varid);
      nc_put_att_text(ncid, varid, "description", 20,"Total backscattering");
      nc_put_att_text(ncid, varid, "units", 3,"m-1");
      nc_put_att_double(ncid, varid, "air_value", NC_DOUBLE, 1,&air);
      nc_put_att_double(ncid, varid, "seafloor_value", NC_DOUBLE, 1,&seafloor);
      nc_put_att_text(ncid, varid, "puv_uom",52,"http://vocab.nerc.ac.uk/collection/P06/current/UPRM/");

      nc_def_var(ncid,"blt",NC_DOUBLE,3,dims9,&varid);
      nc_put_att_text(ncid, varid, "description", 18,"Total luminescence");
      nc_put_att_text(ncid, varid, "units", 10,"W m-3 nm-1");
      nc_put_att_double(ncid, varid, "air_value", NC_DOUBLE, 1,&air);
      nc_put_att_double(ncid, varid, "seafloor_value", NC_DOUBLE, 1,&seafloor);
      
      nc_def_var(ncid,"flt",NC_DOUBLE,3,dims9,&varid);
      nc_put_att_text(ncid, varid, "description", 18,"Total fluorescence");
      nc_put_att_text(ncid, varid, "units", 10,"W m-3 nm-1");
      nc_put_att_double(ncid, varid, "air_value", NC_DOUBLE, 1,&air);
      nc_put_att_double(ncid, varid, "seafloor_value", NC_DOUBLE, 1,&seafloor);
      
      dims3[0] = waveid;
      nc_def_var(ncid,"wave_centre",NC_DOUBLE,1,dims3,&varid);
      nc_put_att_text(ncid, varid, "description", 23,"Mean wavelength of band");
      nc_put_att_text(ncid, varid, "units", 2,"nm");
      
      nc_def_var(ncid,"bandwidth",NC_DOUBLE,1,dims3,&varid);
      nc_put_att_text(ncid, varid, "description", 24,"Width between band edges");
      nc_put_att_text(ncid, varid, "units", 2,"nm");
      
#ifdef INCLUDE_HYDROLIGHT
      
      ncw_def_dim(ws->column_out_name,ncid,"zenith_centre",10, &zenid); // HARDWIRED DIMENSION LENGTH
      dims3[0] = zenid;
      nc_def_var(ncid,"zenith_centre",NC_DOUBLE,1,dims3,&varid);
      nc_put_att_text(ncid, varid, "description", 30,"Radiance grid, zenith centre");
      nc_put_att_text(ncid, varid, "units", 7,"degrees");
      ncw_def_dim(ws->column_out_name,ncid,"azimuth_centre",24, &aziid); // HARDWIRED DIMENSION LENGTH
      dims3[0] = aziid;
      nc_def_var(ncid,"azimuth_centre",NC_DOUBLE,1,dims3,&varid);
      nc_put_att_text(ncid, varid, "description", 31,"Radiance grid, azimuth centre");
      nc_put_att_text(ncid, varid, "units", 7,"degrees");
      
      int dims5[5] = {recid,waveid,depthid,zenid,aziid};
      
      nc_def_var(ncid,"L_u",NC_DOUBLE,5,dims5,&varid);
      nc_put_att_text(ncid, varid, "description", 18,"Upwelling radiance");
      nc_put_att_text(ncid, varid, "units", 15,"W m-2 nm-1 sr-1");
      nc_put_att_double(ncid, varid, "air_value", NC_DOUBLE, 1,&air);
      nc_put_att_double(ncid, varid, "seafloor_value", NC_DOUBLE, 1,&seafloor);
      
      nc_def_var(ncid,"L_d",NC_DOUBLE,5,dims5,&varid);
      nc_put_att_text(ncid, varid, "description", 20,"Downwelling radiance");
      nc_put_att_text(ncid, varid, "units", 15,"W m-2 nm-1 sr-1");
      nc_put_att_double(ncid, varid, "air_value", NC_DOUBLE, 1,&air);
      nc_put_att_double(ncid, varid, "seafloor_value", NC_DOUBLE, 1,&seafloor);
      
      int dims7[4] = {recid,waveid,zenid,aziid};
      nc_def_var(ncid,"L_wl",NC_DOUBLE,4,dims7,&varid);
      nc_put_att_text(ncid, varid, "description", 22,"Water-leaving radiance");
      nc_put_att_text(ncid, varid, "units", 15,"W m-2 nm-1 sr-1");
      nc_put_att_double(ncid, varid, "air_value", NC_DOUBLE, 1,&air);
      nc_put_att_double(ncid, varid, "seafloor_value", NC_DOUBLE, 1,&seafloor);
      
      nc_def_var(ncid,"L_sky",NC_DOUBLE,4,dims7,&varid);
      nc_put_att_text(ncid, varid, "description", 12,"Sky radiance");
      nc_put_att_text(ncid, varid, "units", 15,"W m-2 nm-1 sr-1");
      nc_put_att_double(ncid, varid, "air_value", NC_DOUBLE, 1,&air);
      nc_put_att_double(ncid, varid, "seafloor_value", NC_DOUBLE, 1,&seafloor);
      
      nc_def_var(ncid,"L_ref",NC_DOUBLE,4,dims7,&varid);
      nc_put_att_text(ncid, varid, "description", 18,"Reflected radiance");
      nc_put_att_text(ncid, varid, "units", 15,"W m-2 nm-1 sr-1");
      nc_put_att_double(ncid, varid, "air_value", NC_DOUBLE, 1,&air);
      nc_put_att_double(ncid, varid, "seafloor_value", NC_DOUBLE, 1,&seafloor);
      
#endif

      // Add BGC state variables.

      int dims10[3] = {recid,depthcentreid};
      
      nc_def_var(ncid,"temp",NC_DOUBLE,2,dims10,&varid);
      nc_put_att_text(ncid, varid, "description", 11,"Temperature");
      nc_put_att_text(ncid, varid, "units", 9,"degrees C");

      if (ws->K_heat_i > -1){
	nc_def_var(ncid,"K_heat",NC_DOUBLE,2,dims10,&varid);
	nc_put_att_text(ncid, varid, "description", 30,"Vertical attenuation of heat");
	nc_put_att_text(ncid, varid, "units", 3,"m-1");
	nc_put_att_text(ncid, varid, "puv_uom",52,"http://vocab.nerc.ac.uk/collection/P06/current/UPRM/");
      }
      
      if (ws->PhyL_Chl_i > -1){
	nc_def_var(ncid,"PhyL_Chl",NC_DOUBLE,2,dims10,&varid);
	nc_put_att_text(ncid, varid, "description", 25,"Large phytoplankton Chl a");
	nc_put_att_text(ncid, varid, "units", 6,"mg m-3");
	nc_def_var(ncid,"PhyL_N",NC_DOUBLE,2,dims10,&varid);
	nc_put_att_text(ncid, varid, "description", 21,"Large phytoplankton N");
	nc_put_att_text(ncid, varid, "units", 6,"mg m-3");
      }

      if (ws->PhyS_Chl_i > -1){
	nc_def_var(ncid,"PhyS_Chl",NC_DOUBLE,2,dims10,&varid);
	nc_put_att_text(ncid, varid, "description", 25,"Small phytoplankton Chl a");
	nc_put_att_text(ncid, varid, "units", 6,"mg m-3");
	nc_def_var(ncid,"PhyS_N",NC_DOUBLE,2,dims10,&varid);
	nc_put_att_text(ncid, varid, "description", 21,"Small phytoplankton N");
	nc_put_att_text(ncid, varid, "units", 6,"mg m-3");
      }

      if (ws->Tricho_Chl_i > -1){
	nc_def_var(ncid,"Tricho_Chl",NC_DOUBLE,2,dims10,&varid);
	nc_put_att_text(ncid, varid, "description", 19,"Trichodesmium Chl a");
	nc_put_att_text(ncid, varid, "units", 6,"mg m-3");
	nc_def_var(ncid,"Tricho_N",NC_DOUBLE,2,dims10,&varid);
	nc_put_att_text(ncid, varid, "description", 15,"Trichodesmium N");
	nc_put_att_text(ncid, varid, "units", 6,"mg m-3");
      }

      ncw_close(ws->column_out_name,ncid);
      }
      
      nc_open(ws->column_out_name,NC_WRITE, &ncid);
      
      varid = ncw_var_id(ncid,"wave_centre");
      ncw_put_var_double(ws->column_out_name,ncid,varid,ws->wave);

      varid = ncw_var_id(ncid,"wave_sensor");
      ncw_put_var_double(ws->column_out_name,ncid,varid,ws->sensor_waves);
      
      varid = ncw_var_id(ncid,"bandwidth");
      ncw_put_var_double(ws->column_out_name,ncid,varid,ws->bandwidth);
      
#ifdef INCLUDE_HYDROLIGHT
      varid = ncw_var_id(ncid,"zenith_centre");
      ncw_put_var_double(ws->column_out_name,ncid,varid,zengrid);
      
      varid = ncw_var_id(ncid,"azimuth_centre");
      ncw_put_var_double(ws->column_out_name,ncid,varid,azigrid);
#endif
      
      ncw_close(ws->column_out_name,ncid);

      /* copy complete netcdf header file to the number of columns required. */

      char unixcommand[MAXSTRLEN];
      for (w=0; w<ws->num_cols_out; w++){
	if (output_path){
	  sprintf(unixcommand, "cp %soptical_column_out.nc %soptical_column_out_%d.nc",output_path,output_path,w);
	}else{
	  sprintf(unixcommand, "cp optical_column_out.nc optical_column_out_%d.nc",w);
	}
	system(unixcommand);
      }
      sprintf(unixcommand, "rm %soptical_column_out.nc",output_path);
      system(unixcommand);
    }
  }
}

  void light_spectral_col_destroy(eprocess* p)
{
  workspace* ws = p->workspace;

  // freeing memory created in init.
  
  d_free_1d(ws->sensor_waves);
  
  d_free_1d(ws->bandedge);
  d_free_1d(ws->bandwidth);
  d_free_1d(ws->landa);
  d_free_1d(ws->aw);
  d_free_1d(ws->bw);
  d_free_1d(ws->bls);
  d_free_1d(ws->fls);

  d_free_1d(ws->yC_s);
  d_free_1d(ws->yC_l);
  d_free_1d(ws->yC_Tricho);
  d_free_1d(ws->yC_MPB);
  d_free_1d(ws->yC_D);
  d_free_1d(ws->yC_Symbiodinium);

  d_free_1d(ws->rhoskel);
  if (ws->SMA_N_i > -1){
    d_free_1d(ws->SMA_absorption);
  }
  d_free_1d(ws->moon_a0);
  d_free_1d(ws->moon_a1);
  d_free_1d(ws->moon_a2);
  d_free_1d(ws->moon_a3);
  d_free_1d(ws->moon_b1);
  d_free_1d(ws->moon_b2);
  d_free_1d(ws->moon_b3);
  d_free_1d(ws->moon_d1);
  d_free_1d(ws->moon_d2);
  d_free_1d(ws->moon_d3);
  
  if (ws->at_star)
    d_free_2d(ws->at_star);
  if (ws->bt_star)
    d_free_2d(ws->bt_star);
  if (ws->bb_star)
    d_free_2d(ws->bb_star);
  if (ws->br_star)
    d_free_2d(ws->br_star);
  if (ws->sp3_star)
    d_free_2d(ws->sp3_star);
  if (ws->A_star)
    d_free_2d(ws->A_star);
  if (ws->R_star)
    d_free_2d(ws->R_star);
  if (ws->T_star)
    d_free_2d(ws->T_star);
  if (ws->sp2_star)
    d_free_2d(ws->sp2_star);

  if (ws->XXXleafden)
    d_free_1d(ws->XXXleafden);
  if (ws->XXXorient)
    d_free_1d(ws->XXXorient); 

  free(ws);
}


void light_spectral_col_precalc(eprocess* p, void* pp)
{
  ecology* e = p->ecology;
  column *col = (column *) pp;
  workspace* ws = p->workspace;
  double **y = col->y;
  double *y_sed0 = col->y_sed0;
  double *zc = col->zc;
  double *dz = col->dz;
  double *dz_sed = col->dz_sed;
  double *y_epi = col->y_epi;
  int w,n,i;
  
  int ij[2];
  
  ginterface_get_ij(e->model, col->b, ij);

  /* Column to do hydrolight and / or output spectrally-resolved variables. */

  int col_write = 10000000;
  
  int col_write_num = -1;
  
  for (w=0; w<ws->num_cols_out; w++) {
    if (ij[0] == ws->i_indices[w] && ij[1] == ws->j_indices[w]){
      col_write = col->b;
      col_write_num = w;
    }
  }

 /**************************************/
 /* Internal, spectrally-resolved IOPs */
 /**************************************/

  double** at = d_alloc_2d(col->n_wc,ws->num_wave);  // indexing is reversed.
  double** bt = d_alloc_2d(col->n_wc,ws->num_wave);  // indexing is reversed.
  double** bb = d_alloc_2d(col->n_wc,ws->num_wave);  // indexing is reversed.
  double** blt = d_alloc_2d(col->n_wc,ws->num_wave);  // spectrally-resolved luminescence per waveband
  double** flt = d_alloc_2d(col->n_wc,ws->num_wave);  // total fluorescence

  // Clear sea-water

  for (n = 0; n<col->n_wc; n++) { // 0 is top
    for (w = 0; w<ws->num_wave; w++) {
      at[w][n] = ws->aw[w];
      bt[w][n] = ws->bw[w];
      bb[w][n] = 0.5 * ws->bw[w];
      blt[w][n] = 0.0;  //40K decay causes 4e-12 W m-2, peaks at 380-420 nm (Fig. 8, Massa, 2002).
      flt[w][n] = 0.0;
    }
  }

  // Suspended macroalgae

  if (ws->SMA_N_i > -1){
    for (n = 0; n<col->n_wc; n++) { // ocean
      if (y[ws->SMA_N_i][n] > 1.0e-10){
	for (w = 0; w<ws->num_wave; w++){
	  //at[w][n] += -log(1.0 - min(0.9999,exp(-y[ws->SMA_N_i][n] * dz[n] * ws->MAleafden * fabs(ws->SMA_absorption[w]))))/ dz[n];
	  at[w][n] += y[ws->SMA_N_i][n] * ws->MAleafden * fabs(ws->SMA_absorption[w]);
	}
      }
    }
  }
  
  // CDOM absorption

  if (ws->cdom_gbr_i == - 1){ // ocean component only put in if not GBR.
    for (n = 0; n<col->n_wc; n++) { // ocean
      for (w = 0; w<ws->num_wave; w++){
	  at[w][n] += ws->a440cdom_ocean * exp(- ws->Scdom_ocean * (ws->wave[w]-443.0));
      }
    }
  }
  if (ws->cdom_pale_i > -1){
    for (n = 0; n<col->n_wc; n++) {
      for (w = 0; w<ws->num_wave; w++) {
	at[w][n] += ws->a440cdom_pale * min(y[ws->cdom_pale_i][n],1.0) * exp(- ws->Scdom_pale * (ws->wave[w]-440.0));
      }
    }
  }
  if (ws->cdom_amber_i > -1){
    for (n = 0; n<col->n_wc; n++) {
      for (w = 0; w<ws->num_wave; w++) {
	at[w][n] += ws->a440cdom_amber * min(y[ws->cdom_amber_i][n],1.0) * exp(- ws->Scdom_amber * (ws->wave[w]-440.0));
      }
    }
  }
  if (ws->cdom_dark_i > -1){
    for (n = 0; n<col->n_wc; n++) {
      for (w = 0; w<ws->num_wave; w++) {
	at[w][n] += ws->a440cdom_dark * min(y[ws->cdom_dark_i][n],1.0) * exp(- ws->Scdom_dark * (ws->wave[w]-440.0));
      }
    }
  }
  if (ws->cdom_gbr_i > -1){
    double a443,S_cdom;
    for (n = 0; n<col->n_wc; n++) { // Variable S_cdom for GBR and includes ocean component.
      a443 = (1.2336 + (-0.0332 * (1.0-y[ws->cdom_gbr_i][n])*36.855)); // matching Schroeder salinity relationship.
      S_cdom = 0.0061 * pow(a443,-0.309); // Blondeau-Patissier 2009 JGR: 114: C05003.
      for (w = 0; w<ws->num_wave; w++) {
	at[w][n] += a443 * exp(-S_cdom * (ws->wave[w]-443.0));
      }
    }
  }
  
  // Microalgal absorption and scattering

  double aAtemp[ws->num_wave];
  double rad, vol, m, cellnum, cellChl;

  double* absorbance = d_alloc_1d(ws->num_wave);

  double** aA_s_s = d_alloc_2d(col->n_wc,ws->num_wave);
  double** aA_s_l = d_alloc_2d(col->n_wc,ws->num_wave);   
  double** aA_s_Tricho = d_alloc_2d(col->n_wc,ws->num_wave);    
  double** aA_s_MPB = d_alloc_2d(col->n_wc,ws->num_wave);
  double** aA_s_PhyD = d_alloc_2d(col->n_wc,ws->num_wave);
  double** aA_s_CS_free = d_alloc_2d(col->n_wc,ws->num_wave);
 
  double* aA_s_CS = d_alloc_1d(ws->num_wave);
  double* aA_s_MPB_sed = d_alloc_1d(ws->num_wave);
  double* aA_s_CS_free_sed = d_alloc_1d(ws->num_wave);

  double* biolum = d_alloc_1d(ws->num_wave);
  double* fluor = d_alloc_1d(ws->num_wave);

  if (ws->PhyL_N_i>-1){

    double* Phy_N = y[ws->PhyL_N_i];
    double* Phy_Chl = y[ws->PhyL_Chl_i];
    double* Phy_I = y[ws->PhyL_I_i];

    rad = ws->rad_l;
    vol = ws->vol_l;
    m = ws->m_l;
    
    for (n = 0; n<col->n_wc; n++) { // 0 is top  

      cellnum =  Phy_N[n] / (m * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
      cellChl = Phy_Chl[n] / (vol * cellnum); /* cellular Chl conc */    
 
      for (w=0; w<ws->num_wave; w++){
	absorbance[w] = ws->yC_l[w] * cellChl;
      }
      aawave(rad, absorbance, aAtemp, ws->num_wave,ws->wave);

      double Iquota = 16.0 / 1060.0 * (14.01 / 12.01) * Phy_I[n] / Phy_N[n];
      double kI = y[ws->PhyL_kI_i][n];

      #ifdef INCLUDE_BIOLUMINESCENCE
      bioluminescence_per_cell(ws->num_wave,ws->bls_sum,ws->bls,biolum); // at present this a constant array.
      #endif
      
      passive_fluorescence_per_cell(e,ws, kI, Iquota, fluor);

      for (w=0; w<ws->num_wave; w++){
	aA_s_l[w][n] = aAtemp[w];
	at[w][n] += aA_s_l[w][n] * cellnum;
	bt[w][n] += ws->bphy * Phy_N[n] / ws->NtoCHL;
	//bb[w][n] += ws->nbt_p_555[w] * ws->bbcell1 * pow(rad*2e6, ws->bbcell2)*cellnum;
	bb[w][n] += ws->nbt_p_555[w] * ws->bbfact_l * cellnum;
	blt[w][n] += biolum[w] * cellnum;
	flt[w][n] += fluor[w] * cellnum;
      }
    }
  }

  if (ws->PhyS_N_i> -1){

    double* Phy_N = y[ws->PhyS_N_i];
    double* Phy_Chl = y[ws->PhyS_Chl_i];

    rad = ws->rad_s;
    vol = ws->vol_s;
    m = ws->m_s;
    
    for (n = 0; n<col->n_wc; n++) { // 0 is top
      cellnum =  Phy_N[n] / (m * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
      cellChl = Phy_Chl[n] / (vol * cellnum); /* cellular Chl conc */ 

      for (w=0; w<ws->num_wave; w++){
	absorbance[w] = ws->yC_s[w] * cellChl;
	
      }
      aawave(rad, absorbance, aAtemp, ws->num_wave,ws->wave);

      for (w=0; w<ws->num_wave; w++){
	aA_s_s[w][n] = aAtemp[w];
	at[w][n] += aA_s_s[w][n] * cellnum;
	bt[w][n] += ws->bphy * Phy_N[n] / ws->NtoCHL;
	// bb[w][n] += ws->nbt_p_555[w] * ws->bbcell1 * pow(rad*2e6, ws->bbcell2)*cellnum;
	bb[w][n] += ws->nbt_p_555[w] * ws->bbfact_s * cellnum;
      }
    }
  }

  if (ws->Tricho_N_i> -1){

    double* Phy_N = y[ws->Tricho_N_i];
    double* Phy_Chl = y[ws->Tricho_Chl_i];

    rad = ws->rad_Tricho;
    vol =  ws->vol_Tricho;
    m = ws->m_Tricho;
    
    for (n = 0; n<col->n_wc; n++) { // 0 is top

      cellnum =  Phy_N[n] / (m * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
      cellChl = Phy_Chl[n] / (vol * cellnum); /* cellular Chl conc */ 

      for (w=0; w<ws->num_wave; w++){
	absorbance[w] = ws->yC_Tricho[w] * cellChl;	
      }
      
      aawave(rad, absorbance, aAtemp, ws->num_wave,ws->wave);

      for (w=0; w<ws->num_wave; w++){
	aA_s_Tricho[w][n] = aAtemp[w];
	at[w][n] += aA_s_Tricho[w][n] * cellnum;
	bt[w][n] += ws->bphy * Phy_N[n] / ws->NtoCHL;
	// bb[w][n] += ws->nbt_p_555[w] * ws->bbcell1 * pow(rad*2e6, ws->bbcell2)*cellnum;
	bb[w][n] += ws->nbt_p_555[w] * ws->bbfact_Tricho * cellnum;
      }
    }
  }

  if (ws->MPB_N_i> -1){

    double* Phy_N = y[ws->MPB_N_i];
    double* Phy_Chl = y[ws->MPB_Chl_i];

    rad = ws->rad_MPB;
    vol =  ws->vol_MPB;
    m = ws->m_MPB;
    
    for (n = 0; n<col->n_wc; n++) { // 0 is top

      cellnum =  Phy_N[n] / (m * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
      cellChl = Phy_Chl[n] / (vol * cellnum); /* cellular Chl conc */ 

      for (w=0; w<ws->num_wave; w++){
	absorbance[w] = ws->yC_MPB[w] * cellChl;	
      }
      
      aawave(rad, absorbance, aAtemp, ws->num_wave,ws->wave);

      for (w=0; w<ws->num_wave; w++){
	aA_s_MPB[w][n] = aAtemp[w];
	at[w][n] += aA_s_MPB[w][n] * cellnum;
	bt[w][n] += ws->bphy * Phy_N[n] / ws->NtoCHL;
	// bb[w][n] += ws->nbt_p_555[w] * ws->bbcell1 * pow(rad*2e6, ws->bbcell2)*cellnum;
	bb[w][n] += ws->nbt_p_555[w] * ws->bbfact_MPB * cellnum;
      }
    }
  }

  if (ws->PhyD_N_i>-1){

    double* Phy_N = y[ws->PhyD_N_i];
    double* Phy_Chl = y[ws->PhyD_Chl_i];

    rad = ws->rad_PhyD;
    vol = ws->vol_PhyD;
    m = ws->m_PhyD;
    
    for (n = 0; n<col->n_wc; n++) { // 0 is top  

      cellnum =  Phy_N[n] / (m * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
      cellChl = Phy_Chl[n] / (vol * cellnum); /* cellular Chl conc */    
 
      for (w=0; w<ws->num_wave; w++){
	absorbance[w] = ws->yC_D[w] * cellChl;
      }
      aawave(rad, absorbance, aAtemp, ws->num_wave,ws->wave);

      for (w=0; w<ws->num_wave; w++){
	aA_s_PhyD[w][n] = aAtemp[w];
	at[w][n] += aA_s_PhyD[w][n] * cellnum;
	bt[w][n] += ws->bphy * Phy_N[n] / ws->NtoCHL;
	// bb[w][n] += ws->nbt_p_555[w] * ws->bbcell1 * pow(rad*2e6, ws->bbcell2)*cellnum;
	bb[w][n] += ws->nbt_p_555[w] * ws->bbfact_PhyD * cellnum;
      }
    }
  }

  if (ws->CS_N_free_i>-1){   // this was the last one I avoided.

    /* only do using mass-specific absorption coefficients. */

    double* Phy_N = y[ws->CS_N_free_i];
    double* Phy_Chl = y[ws->CS_Chl_free_i];
   
    for (n = 0; n<col->n_wc; n++) { // 0 is top  

      if (Phy_N[n] > 0.0){

	cellnum =  Phy_N[n] / (ws->m_CS * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
	cellChl = Phy_Chl[n] / (ws->vol_CS * cellnum); /* cellular Chl conc */

	double cellXp = y[ws->CS_Xp_free_i][n] / (ws->vol_CS * cellnum); /* cellular Xp conc */
	double cellXh = y[ws->CS_Xh_free_i][n] / (ws->vol_CS * cellnum); /* cellular Xh conc */
	
	for (w=0; w<ws->num_wave; w++){
	  absorbance[w] = cellXh * ws->yC_diatoxanthin[w] + cellXp * ws->yC_diadinoxanthin[w] + cellChl * ws->yC_Symbiodinium[w];
	}
 
	aawave(ws->rad_CS, absorbance, aAtemp, ws->num_wave,ws->wave);

	for (w=0; w<ws->num_wave; w++){
	  aA_s_CS_free[w][n] = aAtemp[w];
	  at[w][n] += aA_s_CS_free[w][n] * cellnum;
	  bt[w][n] += ws->bphy * Phy_N[n] / ws->NtoCHL;
	  // bb[w][n] += ws->nbt_p_555[w] * ws->bbcell1 * pow(rad*2e6, ws->bbcell2)*cellnum;
	  bb[w][n] += ws->nbt_p_555[w] * ws->bbfact_CS * cellnum;
	}
      }
    }
  }
    
    // Non-algal particulate scattering and absorption.

  // From tracer list

  for (i=0; i<ws->ntr_at; i++){
    for (n = 0; n<col->n_wc; n++) { // 0 is top
      for (w=0; w<ws->num_wave; w++){
	at[w][n] += y[ws->ind_at[i]][n] * ws->at_star[i][w];
      }
    }
  }
  for (i=0; i<ws->ntr_bt; i++){
    for (n = 0; n<col->n_wc; n++) { // 0 is top
      for (w=0; w<ws->num_wave; w++){
	bt[w][n] += y[ws->ind_bt[i]][n] * ws->bt_star[i][w];
      }
    }
  }
  for (i=0; i<ws->ntr_bt; i++){
    for (n = 0; n<col->n_wc; n++) { // 0 is top
      for (w=0; w<ws->num_wave; w++){
	bb[w][n] += y[ws->ind_bb[i]][n] * ws->bb_star[i][w] * ws->bt_star[i][w];
      }
    }
  }
  
  double all_sed = 0.0;
  double u_bot[ws->num_wave];

  double f_BP = 0.0;  // sum of benthic plants.

  double f_polyp = 0.0;
  double f_mpb = 0.0;
  double f_zoo_free = 0.0;
  double f_zoo = 0.0;
  double f_skel = 0.0;
  double frac = 0.0;
  
  /*
   * Here we progress from the tallest species down to the
   * sediments, calculating the fraction of light that is still
   * penetrating down after each reflecting thing
   *
   * Note: variable dependent on things like ws->MAleafden should've
   *       been set to zero above in the case where there is no MA
   *       and so on...
   * frac - what is available at that depth.
   * u_bot - addition to the bottom - bb / 
   */
  
  /* B4p0: Impact of leaves on light defined by A(NP) + A + R + T = 1, potentially from one study.  */

  /* B3p1: Did not consider leaf reflection or non-photosynthetic absorption, A + T = 1, and used a seperately determined R.
   * The R is applied to the area on the bottom that substrate takes up. So this is missing the component that is transmitted through the 
   * leaf, and then reflected by a lower layer. So the total bottom reflection will depend slightly more on the taller benthic plants than 
   * it should, and then is a potential mismatch using an oberved A and R from different experiments. */

  for (w = 0; w<ws->num_wave; w++) {
    u_bot[w] = 0.0;
  }

  double f_XXX;
  for (i = 0; i < ws->nepi_R; i++) {
    double BP_N = y_epi[ws->ind_R[i]];
    f_XXX = (1.0 - f_BP) * (1.0 - exp(-BP_N * ws->XXXorient[i] * ws->XXXleafden[i]));
    f_BP += f_XXX;
    for (w = 0; w<ws->num_wave; w++) {
      u_bot[w] += f_XXX * ws->R_star[i][w];
    }
  }
  
  /*
   * Corals : Polyps, zoos then skeleton
   */
  
  // WARNING: DOUBLE DEFINITION OF ABSORBANCE FOR SEDIMENT LAYER MPB AND CS.
  
  if (ws->CS_N_i > -1){

    double CH_N = y_epi[ws->CH_N_i];
    double CS_N = y_epi[ws->CS_N_i];
    double CS_Chl = y_epi[ws->CS_Chl_i];
    
    if (CS_N > 0.0){
      
      double cellnum = CS_N / (ws->m_CS* 1000.0 * red_A_N * MW_Nitr); /* cell m-2 */
      double cellChl = CS_Chl / (ws->vol_CS * cellnum); /* cellular Chl conc */
      
      if (ws->CS_Xp_i > -1){
	
	double cellXp = y_epi[ws->CS_Xp_i] / (ws->vol_CS * cellnum); /* cellular Xp conc */
	double cellXh = y_epi[ws->CS_Xh_i] / (ws->vol_CS * cellnum); /* cellular Xh conc */
	
	/* Add photosynthetic cart to absorption */
	
	for (w=0; w<ws->num_wave; w++){
	  absorbance[w] = cellXp * ws->yC_diadinoxanthin[w] + cellChl * ws->yC_Symbiodinium[w];
	}
      }else{
	absorbwave(absorbance,cellChl,ws->num_wave, ws->wave); // old way.
      }
      aawave(ws->rad_CS, absorbance, aA_s_CS, ws->num_wave, ws->wave);
      
      double bb_zoo,u_zoo;
      
      double SA = ws->CHarea * (1.0-exp(-CH_N * ws->CHpolypden / ws->CHarea)); /* dimensionless */
      
      f_polyp = (1.0 - f_BP) * SA;  // polyps
      
      f_zoo = min(f_polyp, M_PI*M_PI/(2.0*sqrt(3.0)) * cellnum * ws->rad_CS * ws->rad_CS);
      
      f_skel = f_polyp - f_zoo;
      
      /* zoos need to consider bb/(a+bb) */
      
      for (w = 0; w<ws->num_wave; w++) {
	bb_zoo = ws->nbt_p_555[w] * ws->bbcell1 * pow(ws->rad_CS*2e6, ws->bbcell2);
	u_zoo = bb_zoo/(bb_zoo + aA_s_CS[w]);
	u_bot[w] += f_zoo * u_zoo + f_skel * ws->rhoskel[w];
      }
    }
  }

  /*
   * Microphytobenthos perfectly packed on the bottom using multiple pigments.
   *
   */
  
  if (ws->MPB_N_sed_i > -1){
    
    double u_mpb,bb_mpb;
    
    double cellnum_MPB = y_sed0[ws->MPB_N_sed_i] / (ws->m_MPB * 1000.0 * red_A_N * MW_Nitr); /* cell m-2 */
    double cellChl_MPB = y_sed0[ws->MPB_Chl_sed_i] / (ws->vol_MPB * cellnum_MPB); /* cellular Chl conc */
    
    for (w=0; w<ws->num_wave; w++){
      absorbance[w] = ws->yC_MPB[w] * cellChl_MPB;	
    }
    
    aawave(ws->rad_MPB, absorbance, aA_s_MPB_sed, ws->num_wave,ws->wave);
    
    f_mpb = min(1.0 - f_BP - f_polyp, M_PI*M_PI/(2.0*sqrt(3.0)) * cellnum_MPB * ws->rad_MPB * ws->rad_MPB * (col->dz_sed[ws->num_sed_layers-1]));
    
    for (w = 0; w < ws->num_wave; w++) {
      bb_mpb = ws->nbt_p_555[w] * ws->bbcell1 * pow(ws->rad_MPB*2.0e6, ws->bbcell2);
      u_mpb = bb_mpb/(bb_mpb + aA_s_MPB_sed[w]);
      u_bot[w] += f_mpb * u_mpb;
    }
  }

  /*
   * Free-living symbionts perfectly packed on the bottom using multiple pigments.
   */
  
  if (ws->CS_N_free_i > -1){ // only considering bleaching model means.
    
    double u_cs,bb_cs;
    
    double cellnum_CS = y_sed0[ws->CS_N_free_i] / (ws->m_CS * 1000.0 * red_A_N * MW_Nitr); /* cell m-2 */
    double cellChl_CS = y_sed0[ws->CS_Chl_free_i] / (ws->vol_CS * cellnum_CS); /* cellular Chl conc */
    
    double cellXp_sed = y_sed0[ws->CS_Xp_free_i] / (ws->vol_CS * cellnum); /* cellular Xp conc */
    double cellXh_sed = y_sed0[ws->CS_Xh_free_i] / (ws->vol_CS * cellnum); /* cellular Xh conc */
    
    for (w=0; w<ws->num_wave; w++){
      absorbance[w] = cellXh_sed * ws->yC_diatoxanthin[w] + cellXp_sed * ws->yC_diadinoxanthin[w] + cellChl_CS * ws->yC_Symbiodinium[w];
    }
    
    aawave(ws->rad_CS, absorbance, aA_s_CS_free_sed, ws->num_wave,ws->wave);

    f_zoo_free = min(1.0 - f_BP - f_polyp - f_mpb, M_PI * M_PI/(2.0*sqrt(3.0)) * cellnum_CS * ws->rad_CS * ws->rad_CS * (col->dz_sed[ws->num_sed_layers-1]));
    
    for (w = 0; w < ws->num_wave; w++) {
      bb_cs = ws->nbt_p_555[w] * ws->bbcell1 * pow(ws->rad_CS*2.0e6, ws->bbcell2);
      u_cs = bb_cs/(bb_cs + aA_s_CS_free_sed[w]);
      u_bot[w] += f_zoo_free * u_cs;
    }
  }

  /* Now do sediment base */

  frac = 1.0 - f_BP - f_polyp - f_mpb - f_zoo_free;

  /* 
   * Sediment particulates : we seemed to add up only the non-Gravel components
   * also, this looks like it will add detritus: need benthic reflectance with mg m-2 units.
   */

  all_sed = 0.0;
  for (i = 0; i < ws->ntr_br; i++) {
    all_sed += y_sed0[ws->ind_br[i]];
  }
  
  for (i = 0; i < ws->ntr_br; i++) {
    double frac2 = y_sed0[ws->ind_br[i]]/all_sed;
    for (w = 0; w<ws->num_wave; w++) {
      u_bot[w] += frac * frac2 * ws->br_star[i][w];
      
    }
  }
  
  /***********************************************************/
  /* Calculate spectrally-resolved AOPs using original model */
  /***********************************************************/

  /* get atmospheric conditions */

  double longitude,latitude; // for ginterface call
  ginterface_get_lat_lon(col->model,col->b,&latitude,&longitude);

  double wind_speed,tcld,mslp,rh2m;  
  ginterface_get_windspeed_and_cloud_and_mslp_and_rh(col->model,col->b,&wind_speed,&tcld,&mslp,&rh2m);
  
  /* Lunar radiation characteristics */
  /* =============================== */

  double* Ak = d_alloc_1d(ws->num_wave);
  double* moonlight_SR = d_alloc_1d(ws->num_wave);
  
  double sum_moonlight = 0.0;
  double moon_phase,lunar_dec,sun_angle,lunar_angle,lunar_zenith,lunar_azimuth,moon_fulldisk;

  double lunar_airtransmission = 0.0;
  double lunar_albedo = 1.0;
  double time_in_days = (e->t + e->dt)/86400.0;
  
  if (ws->Moonlight_i > -1){

    // function needs to return spectrally-resolved moonlight, fulldisk moonlight, lunar zenith and lunar phase.
    
    moonlight(ws,col,time_in_days,Ak,&sun_angle,&moon_phase,&lunar_dec,&lunar_angle,&lunar_zenith,&lunar_azimuth);
    
    // Calculate attenuation through the atmosphere and reflection at the air-water interface.
    
    if (lunar_zenith < (M_PI/2.0 - 1.0e-15)){
      ems_swr_airtrans_and_albedo(tcld,lunar_dec,lunar_angle,lunar_zenith,latitude,&lunar_airtransmission,&lunar_albedo);
    }
    
    /* Fulldisk PAR-integrated moonlight and orbital characteristics */
    
    for (w=0; w<ws->num_wave; w++) {
      sum_moonlight += Ak[w];
      moonlight_SR[w] = 0.0;
    }

    if (lunar_zenith < (M_PI/2.0 - 1.0e-15)){
      for (w=0; w<ws->num_wave; w++) {
	moonlight_SR[w] = Ak[w] * cos(lunar_zenith) / ws->bandwidth[w];
      }
    }

    moon_fulldisk = (lunar_zenith >= M_PI/2.0 - 1.0e-15)? 0.0:sum_moonlight;
    
    if (ws->Moon_fulldisk_i > -1)
      y_epi[ws->Moon_fulldisk_i] = moon_fulldisk;
    if (ws->Lunar_zenith_i > -1)
      y_epi[ws->Lunar_zenith_i] = lunar_zenith;
    if (ws->Lunar_azimuth_i > -1)
      y_epi[ws->Lunar_azimuth_i] = lunar_azimuth;
    if (ws->Lunar_phase_i > -1)
      y_epi[ws->Lunar_phase_i] = moon_phase;
    
    /* Moonlight hitting the sea surface */
  
    y_epi[ws->Moonlight_i] = (lunar_zenith >= M_PI/2.0 - 1.0e-15)? 0.0:sum_moonlight * cos(lunar_zenith) * lunar_airtransmission * (1.0-lunar_albedo); // zero if below horizon. 
  }

  /* Obtain surface light field from hydrodynamic / transport model */

  double thetaw;
  double starlight = 7.5e-6; // PAR clearsky starlight (3.0e-6/0.4, W m-2) from openoceanopticsbook.info

  double ems_solar_zenith = ginterface_calc_zenith(e->model, e->t+e->dt, col->b);
  
  double light = ws->SWRscale * ginterface_getlighttop(col->model, col->b); // previous timestep - we could add an interface function that allows for specificing a time.

  // Calculate refraction at sea-surface - and need to change costhetaw to wavelength dependepent.
     
  double* n_stw = d_alloc_1d(ws->num_wave);
  double* costhetaw = d_alloc_1d(ws->num_wave);

  // ****************** this should be a function ************************************************************ //
  
  // Calculate solar azimuth angle:

  int yr,mon,nday;
  double hr,sunphi,doy,lon,ndoy;
  double d1, dec, hrang, solar_azimuth;

  solar_azimuth = 0.0;

  tm_time_to_ymd(e->t-36000.0 + e->dt,ws->timeunits, &yr, &mon, &nday); // Hydrolight uses UTC - subtract 10 hours.
  
  dtime(NULL,ws->timeunits, e->t-36000.0 + e->dt, &yr, &doy, &longitude);
  ndoy = (int)doy;
  hr = 24.0 * (doy - (double)ndoy); 
  
  d1 = ndoy * 2.0 * M_PI / 365.0;
  dec = 0.006918 + 0.070257 * sin(d1) - 0.399912 * cos(d1)
    + 0.000907 * sin(2.0 * d1) - 0.006758 * cos(2.0 * d1)
    + 0.00148 * sin(3.0 * d1) - 0.002697 * cos(3.0 * d1);
  
  hrang = (hr - 2.0) * M_PI / 12.0;  // hr +10-12
  
  if (ems_solar_zenith < M_PI/2.0){
    solar_azimuth = - sin(hrang) * cos(dec) / sin(ems_solar_zenith);
    solar_azimuth = asin(solar_azimuth);
  }

// ****************** this should be a function ************************************************************ //
  
  // zenith, azimuth determine by sun, lunar then starlight.
  
  double zenith = ems_solar_zenith;
  double azimuth = solar_azimuth;

  if (fabs(ems_solar_zenith) > M_PI/2.0 - 0.02){  // sun below the horizon    
    if (fabs(lunar_zenith) < M_PI/2.0 -0.02){ // moon above the horizon
      zenith = lunar_zenith;
      azimuth = lunar_azimuth;
    }else{
      zenith = 0.0; // starlight - sun and moon below the horizon.
      azimuth = 0.0; 
    }
  }

  if (ws->Solar_zenith_i > -1)
      y_epi[ws->Solar_zenith_i] = ems_solar_zenith;

  if (ws->Zenith_i > -1)
    y_epi[ws->Zenith_i] = zenith;
  if (ws->Azimuth_i > -1)
    y_epi[ws->Azimuth_i] = azimuth;  

  double ems_lunar_zenith = lunar_zenith;

  // do transformation of ems_xxxx_zenith, but just for the output.

  ems_solar_zenith = (fabs(ems_solar_zenith)<1.0e-15)? 1.0e-15:ems_solar_zenith;
  ems_solar_zenith = ((ems_solar_zenith-M_PI/2.0)-fabs(ems_solar_zenith-M_PI/2.0))/2.0+M_PI/2.0;
  
  ems_lunar_zenith = (fabs(ems_lunar_zenith)<1.0e-15)? 1.0e-15:ems_lunar_zenith;
  ems_lunar_zenith = ((ems_lunar_zenith-M_PI/2.0)-fabs(ems_lunar_zenith-M_PI/2.0))/2.0+M_PI/2.0;

  index_of_refraction(e,ws,y,n_stw);

  for (w = 0; w<ws->num_wave; w++) {
    thetaw = asin(sin(zenith)/n_stw[w]); // old code was: asin(sin(zenith)/1.33);
    costhetaw[w] = cos(thetaw);
  }

 /*  avoid exactly zero zenith which causes trouble with albedo calcs. */
  
  zenith = (fabs(zenith)<1.0e-15)? 1.0e-15:zenith;
  zenith = ((zenith-M_PI/2.0)-fabs(zenith-M_PI/2.0))/2.0+M_PI/2.0;
  
  /* Do ray-tracing calculation of light field (Baird et al., 2020) */

  /* E - W m-2 band-1 */

  double** E_d = d_alloc_2d(col->n_wc,ws->num_wave);   // volume-weighted downwelling irradiance
  double** E_s = d_alloc_2d(col->n_wc,ws->num_wave);   // volume-weighted downwelling irradiance scalar irradiance
  double** Kd = d_alloc_2d(col->n_wc,ws->num_wave);    // vertical attenuation at cell centre
  double** E_d_z = d_alloc_2d(col->n_wc+1,ws->num_wave); // downwelling irradiance on cell top
  
  // Starlight and sum have same spectral composition, moon has a different one.

  for (w = 0; w<ws->num_wave; w++) {
    E_d_z[w][0] = (light + starlight) * ws->landa[w] + Ak[w] * cos(lunar_zenith) * lunar_airtransmission * (1.0-lunar_albedo); // Ak[w] is moonlight.
    for (n = 0; n < col->n_wc-1; n++){ // 0 is top  
      Kd[w][n] = at[w][n] * sqrt(1.0 + (ws->gi * costhetaw[w] - ws->gii) * bt[w][n] / at[w][n]) / costhetaw[w];
      E_d_z[w][n+1] = E_d_z[w][n] * exp(-Kd[w][n]*dz[n]);
      E_d[w][n] = (E_d_z[w][n] - E_d_z[w][n+1])/(Kd[w][n]*dz[n]);
      E_s[w][n] = E_d[w][n] * Kd[w][n] / at[w][n];         
    }
    n = col->n_wc-1;
    Kd[w][n] = at[w][n]  * sqrt(1.0 + (ws->gi * costhetaw[w] - ws->gii) * bt[w][n] / at[w][n]) / costhetaw[w];
    E_d_z[w][n+1] = E_d_z[w][n] * exp(-Kd[w][n]*dz[n]);
    E_d[w][n] = (E_d_z[w][n] - E_d_z[w][n+1])/(Kd[w][n]*dz[n]);
    E_s[w][n] = E_d[w][n] * Kd[w][n] / at[w][n];            
  }
    
  d_free_1d(Ak);
  
  // Now do the light absorption and opaqueness calculations that require E_d_z for all layers at once
  
  double sumlight,yCfac,sumabsorb;

  // NOTE: UNIT CHANGE FOR ABSORPTION OF SINGLE CELLS (x 8.35E-9) IS IN VALUES_COMMON.C TO AVOID VALUE OF 1E-12 BEING RESET.

  // Suspended macroalgae

  if (ws->SMA_N_i > -1){

    double SMA_photons;
    
    for (n = 0; n<col->n_wc; n++) { // ocean
      if (y[ws->SMA_N_i][n] > 0.0){
	SMA_photons = 0.0;
	for (w = 0; w<ws->num_wave; w++){
	  SMA_photons += E_d[w][n] * (1.0 - exp(-y[ws->SMA_N_i][n] * dz[n] * ws->MAleafden * ws->SMA_absorption[w])) * ws->wave[w];
	}
	y[ws->SMA_kI_i][n] = SMA_photons * 8.359335857479461e-09;
      }else{
	y[ws->SMA_kI_i][n] = 0.0;
      }
    }
  }
  
  if (ws->PhyL_N_i> -1){
      
    double rad = ws->rad_l;
    double vol = ws->vol_l;
    double m = ws->m_l;
 
    for (n = 0; n < col->n_wc; n++){ // 0 is top
      
      double Phy_N = y[ws->PhyL_N_i][n];

      if (Phy_N > 1e-9){ // this is necessary since kI and yC calculated here are used next step in growth

	double Phy_Chl = y[ws->PhyL_Chl_i][n];     
	double cellnum =  Phy_N / (m * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
	double cellChl = Phy_Chl / (vol * cellnum); /* cellular Chl conc */
	
	sumlight = 0.0;
	yCfac = 0.0;
	sumabsorb = 0.0;
	
	for (w=0; w<ws->num_wave; w++){
	  absorbance[w] = ws->yC_l[w] * cellChl;
	  sumabsorb += E_s[w][n] * ws->wave[w] * aA_s_l[w][n];
	}
	for (w=ws->wPAR_bot; w<ws->wPAR_top; w++){
	  yCfac = yCfac + diffaa(absorbance[w],rad)*max(E_d_z[w][n]*ws->wave[w],0.00001);
	  sumlight += max(E_d_z[w][n] * ws->wave[w],0.00001);
	}
      	y[ws->PhyL_kI_i][n] = sumabsorb * 1.0e20;
	y[ws->PhyL_yCfac_i][n] = yCfac/sumlight;
      }else{
	y[ws->PhyL_kI_i][n] = 0.0;
	y[ws->PhyL_yCfac_i][n] = 0.0;
      }
    }
  }
  
  if (ws->PhyS_N_i> -1){
    
    double rad = ws->rad_s;
    double vol = ws->vol_s;
    double m = ws->m_s;
    
    for (n = 0; n < col->n_wc; n++){ // 0 is top
      
      double Phy_N = y[ws->PhyS_N_i][n];

      if (Phy_N > 1e-9){ // this is necessary since kI and yC calculated here are used next step in growth
	
	double Phy_Chl = y[ws->PhyS_Chl_i][n];
	double cellnum =  Phy_N / (m * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
	double cellChl = Phy_Chl / (vol * cellnum); /* cellular Chl conc */
	
	sumlight = 0.0;
	yCfac = 0.0;
	sumabsorb = 0.0;

	for (w=0; w<ws->num_wave; w++){
	  absorbance[w] = ws->yC_s[w] * cellChl;
	  sumabsorb += E_s[w][n] * ws->wave[w] * aA_s_s[w][n];
	}
	for (w=ws->wPAR_bot; w<ws->wPAR_top; w++){
	  yCfac = yCfac + diffaa(absorbance[w],rad)*max(E_d_z[w][n]*ws->wave[w],0.00001);
	  sumlight += max(E_d_z[w][n] * ws->wave[w],0.00001);
	}
	y[ws->PhyS_kI_i][n] = sumabsorb * 1.0e20;
	y[ws->PhyS_yCfac_i][n] = yCfac/sumlight;
      }else{
	y[ws->PhyS_kI_i][n] = 0.0;
	y[ws->PhyS_yCfac_i][n] = 0.0;
      }
    }
  }

  if (ws->Tricho_N_i> -1){

    double rad = ws->rad_Tricho;
    double vol = ws->vol_Tricho;
    double m = ws->m_Tricho;
    
    for (n = 0; n < col->n_wc; n++){ // 0 is top
      
      double Phy_N = y[ws->Tricho_N_i][n];

      if (Phy_N > 1e-9){ // this is necessary since kI and yC calculated here are used next step in growth

	double Phy_Chl = y[ws->Tricho_Chl_i][n];      
	double cellnum =  Phy_N / (m * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
	double cellChl = Phy_Chl / (vol * cellnum); /* cellular Chl conc */
	
	sumlight = 0.0;
	yCfac = 0.0;
	sumabsorb = 0.0;

	for (w=0; w<ws->num_wave; w++){
	  absorbance[w] = ws->yC_Tricho[w] * cellChl;
	  sumabsorb += E_s[w][n] * ws->wave[w] * aA_s_Tricho[w][n];
	}
	for (w=ws->wPAR_bot; w<ws->wPAR_top; w++){
	  yCfac = yCfac + diffaa(absorbance[w],rad)*max(E_d_z[w][n]*ws->wave[w],0.00001);
	  sumlight += max(E_d_z[w][n] * ws->wave[w],0.00001);
	}
	y[ws->Tricho_kI_i][n] = sumabsorb * 1.0e20;
	y[ws->Tricho_yCfac_i][n] = yCfac/sumlight;
      }else{
	y[ws->Tricho_kI_i][n] = 0.0;
	y[ws->Tricho_yCfac_i][n] = 0.0;
      }
    }
  }
  if (ws->MPB_N_i> -1){

    double rad = ws->rad_MPB;
    double vol = ws->vol_MPB;
    double m = ws->m_MPB;
    
    for (n = 0; n < col->n_wc; n++){ // 0 is top
      
      double Phy_N = y[ws->MPB_N_i][n];

      if (Phy_N > 1e-9){ // this is necessary since kI and yC calculated here are used next step in growth
	
	double Phy_Chl = y[ws->MPB_Chl_i][n];
	double cellnum =  Phy_N / (m * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
	double cellChl = Phy_Chl / (vol * cellnum); /* cellular Chl conc */
	
	sumlight = 0.0;
	yCfac = 0.0;
	sumabsorb = 0.0;

	for (w=0; w<ws->num_wave; w++){
	  absorbance[w] = ws->yC_MPB[w] * cellChl;
	  sumabsorb += E_s[w][n] * ws->wave[w] * aA_s_MPB[w][n];
	}
	for (w=ws->wPAR_bot; w<ws->wPAR_top; w++){
	  yCfac = yCfac + diffaa(absorbance[w],rad)*max(E_d_z[w][n]*ws->wave[w],0.00001);
	  sumlight += max(E_d_z[w][n] * ws->wave[w],0.00001);
	}
	
	y[ws->MPB_kI_i][n] = sumabsorb * 1.0e20;
	y[ws->MPB_yCfac_i][n] = yCfac/sumlight;
      }else{
	y[ws->MPB_kI_i][n] = 0.0;
	y[ws->MPB_yCfac_i][n] = 0.0;
      }
    }
  }

  if (ws->PhyD_N_i> -1){
      
    double rad = ws->rad_PhyD;
    double vol = ws->vol_PhyD;
    double m = ws->m_PhyD;
    
    for (n = 0; n < col->n_wc; n++){ // 0 is top
      
      double Phy_N = y[ws->PhyD_N_i][n];

      if (Phy_N > 1e-9){ // this is necessary since kI and yC calculated here are used next step in growth

	double Phy_Chl = y[ws->PhyD_Chl_i][n];     
	double cellnum =  Phy_N / (m * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
	double cellChl = Phy_Chl / (vol * cellnum); /* cellular Chl conc */
	
	sumlight = 0.0;
	yCfac = 0.0;
	sumabsorb = 0.0;

	for (w=0; w<ws->num_wave; w++){
	  absorbance[w] = ws->yC_D[w] * cellChl;
	  sumabsorb += E_s[w][n] * ws->wave[w] * aA_s_PhyD[w][n];
	}
	for (w=ws->wPAR_bot; w<ws->wPAR_top; w++){
	  yCfac = yCfac + diffaa(absorbance[w],rad)*max(E_d_z[w][n]*ws->wave[w],0.00001);
	  sumlight += max(E_d_z[w][n] * ws->wave[w],0.00001);
	}	
      	y[ws->PhyD_kI_i][n] = sumabsorb * 1.0e20;
	y[ws->PhyD_yCfac_i][n] = yCfac/sumlight;
      }else{
	y[ws->PhyD_kI_i][n] = 0.0;
	y[ws->PhyD_yCfac_i][n] = 0.0;
      }
    }
  }

  if (ws->CS_N_free_i> -1){
     
    for (n = 0; n < col->n_wc; n++){ // 0 is top
      
      double Phy_N = y[ws->CS_N_free_i][n];

      if (Phy_N > 1.0e-9){ // this is necessary since kI and yC calculated here are used next step in growth

	double Phy_Chl = y[ws->CS_Chl_free_i][n];     
	double cellnum =  Phy_N / (ws->m_CS * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
	double cellChl = Phy_Chl / (ws->vol_CS * cellnum); /* cellular Chl conc */

	double cellXp = y[ws->CS_Xp_free_i][n] / (ws->vol_CS * cellnum); /* cellular Xp conc */
	double cellXh = y[ws->CS_Xh_free_i][n] / (ws->vol_CS * cellnum); /* cellular Xh conc */
	
	sumlight = 0.0;
	yCfac = 0.0;
	sumabsorb = 0.0;

	for (w=0; w<ws->num_wave; w++){
	  absorbance[w] = cellXp * ws->yC_diadinoxanthin[w] + cellChl * ws->yC_Symbiodinium[w];
	}
	aawave(ws->rad_CS, absorbance, aAtemp, ws->num_wave,ws->wave);
	for (w=0; w<ws->num_wave; w++){
	  sumabsorb += E_s[w][n] * ws->wave[w] * aAtemp[w];
	}
	for (w=ws->wPAR_bot; w<ws->wPAR_top; w++){
	  yCfac = yCfac + diffaa(absorbance[w],rad)*max(E_d_z[w][n]*ws->wave[w],0.00001);
	  sumlight += max(E_d_z[w][n] * ws->wave[w],0.00001);
	}
      	y[ws->CS_free_kI_i][n] = sumabsorb * 1.0e20;
	y[ws->CS_free_yCfac_i][n] = yCfac/sumlight;

      }else{
	y[ws->CS_free_kI_i][n] = 0.0;
	y[ws->CS_free_yCfac_i][n] = 0.0;
      }
    }
  }

  // Can now do the upwelling bioluminescence.

#ifdef INCLUDE_BIOLUMINESCENCE
  double** E_u_z_bio = d_alloc_2d(col->n_wc+1,ws->num_wave); // downwelling irradiance on cell top
  bioluminescence_upwelling_irradiance(ws->num_wave, col->n_wc, ws->gii, dz, at, bt, blt, E_u_z_bio);
  if (ws->PAR_bio_up_i > -1){
    double energysum;
    for (n = 0; n<col->n_wc; n++) { // 0 is top
      energysum = 0.0;
      for (w=ws->wPAR_bot; w <= ws->wPAR_top; w++){
	energysum +=  E_u_z_bio[w][n];
      }	
      y[ws->PAR_bio_up_i][n] = energysum;
    }
  }
#endif
  

  // Now pass light through to benthic plants

  double* lighttop_s = d_alloc_1d(ws->num_wave);

  for (w=0; w<ws->num_wave; w++){
    lighttop_s[w] = E_d_z[w][col->n_wc];
  }

  double photons;

  // Do all benthic tracers with OPTICAL flag.
    
  for (i = 0; i < ws->nepi_A; i++) {
    
    // Calculate light above seagrass.
    
    if (i == ws->numMA) {  // counter includes zero, so nunMA is the first SG
      
      if (ws->EpiPAR_sg_i > -1){
	
	double light_above_sg = 0.0;  // quantum-weighted and in photon m-2 d-1
	
	for (w=ws->wPAR_bot; w <= ws->wPAR_top; w++){
	  light_above_sg += lighttop_s[w] * ws->wave[w];
	}
	y_epi[ws->EpiPAR_sg_i] = light_above_sg * 86400.0 * 8.359335857479461e-09;
      }
    }
    
    // Do absorption benthic plants - same for seagrass and seaweeds - but no loss due to reflection.
    
    double BP_N = y_epi[ws->ind_A[i]];
    
    if (BP_N > 0.0){
      
      photons = 0.0;
      
      for (w=0; w<ws->num_wave; w++){
	photons += lighttop_s[w] * (1.0 - exp(-BP_N * ws->XXXorient[i] * ws->XXXleafden[i] * ws->A_star[i][w])) * ws->wave[w];
	
	// Line below accounts for both absorption and reflection - doesn't use T_star.
	
	lighttop_s[w] = lighttop_s[w] * exp(-BP_N * ws->XXXorient[i] * ws->XXXleafden[i] * (ws->A_star[i][w] + ws->R_star[i][w]));  
      }
      y_epi[ws->ind_kI[i]] = photons * 8.359335857479461e-09;
    }else{
      y_epi[ws->ind_kI[i]] = 0.0;
    }
  }
  
  if (ws->EpiPAR_i > -1) {  // why does this give zero?

    double light_above_epi = 0.0;  // quantum-weighted and in photon m-2 d-1
    
    for (w=ws->wPAR_bot; w <= ws->wPAR_top; w++){
      light_above_epi += lighttop_s[w] * ws->wave[w];
    }
    y_epi[ws->EpiPAR_i] = light_above_epi * 86400.0 * 8.359335857479461e-09;
  }

  if (ws->CS_N_i > -1){
    
    double CH_N = y_epi[ws->CH_N_i];
    double CS_N = y_epi[ws->CS_N_i];

    double CS_Chl = y_epi[ws->CS_Chl_i];
    
    y_epi[ws->CS_yCfac_i] = 0.0;
    y_epi[ws->CS_kI_i] = 0.0;
      
    if (CS_N > 1e-10){  /* Zooxanthellae in corals - represented like cells not leaves */
      /* so need to find an average light in the layer like in wc */
      
      /* Need intracellular chlorophyll concentration - since cells aren't advected, 
	 I could have intra. chl conc as the state variable. But leave as is, incase we 
	 want to have the ability of zoothanthellae being expelled into water column */
      
      double cellnum = CS_N / (ws->m_CS * 1000.0 * red_A_N * MW_Nitr); /* cell m-2 */
      double cellChl = CS_Chl / (ws->vol_CS * cellnum); /* cellular Chl conc */
      
      double cellXp = y_epi[ws->CS_Xp_i] / (ws->vol_CS * cellnum); /* cellular Xp conc */
      
      /* Only include photosynthetic carts in calculations */
      
      sumlight = 0.0;
      yCfac = 0.0;
      sumabsorb = 0.0;
      
      for (w=0; w<ws->num_wave; w++){
	absorbance[w] = cellXp * ws->yC_diadinoxanthin[w] + cellChl * ws->yC_Symbiodinium[w];
      }
      
      aawave(ws->rad_CS, absorbance, aA_s_CS, ws->num_wave,ws->wave);
      
      for (w=0; w<ws->num_wave; w++){
	sumabsorb += lighttop_s[w] * ws->wave[w] * aA_s_CS[w];
      }
      for (w=ws->wPAR_bot; w<ws->wPAR_top; w++){
	yCfac += diffaa(absorbance[w],ws->rad_CS)*max(lighttop_s[w]*ws->wave[w],0.00001);
	sumlight += max(lighttop_s[w]*ws->wave[w],0.00001);
      }
      y_epi[ws->CS_yCfac_i] = yCfac/sumlight;   
      y_epi[ws->CS_kI_i] = sumabsorb * 1.0e20;
    }
  }
  // ****************************************************************
  // Sediment light model (excluding reflectance component)
  // ****************************************************************

  // Lets put symbiodinium_free above MPB.

  if (ws->CS_N_free_i> -1){

    double Phy_N_sed = y_sed0[ws->CS_N_free_i];

    if (Phy_N_sed > 1e-9){ // this is necessary since kI and yC calculated here are used next step in growth

      double Phy_Chl_sed = y_sed0[ws->CS_Chl_free_i];
    
      rad = ws->rad_CS;
      vol =  ws->vol_CS;
      m = ws->m_CS;
      
      double cellnum_sed =  Phy_N_sed / (m * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
      double cellChl_sed = Phy_Chl_sed / (vol * cellnum_sed); /* cellular Chl conc */ 

      double cellXp_sed = y_epi[ws->CS_Xp_free_i] / (vol * cellnum_sed); /* cellular Xp conc */
      
      for (w=0; w<ws->num_wave; w++){
	absorbance[w] = cellXp_sed * ws->yC_diadinoxanthin[w] + ws->yC_Symbiodinium[w] * cellChl_sed;	
      }
      aawave(rad, absorbance, aA_s_CS_free_sed, ws->num_wave,ws->wave);
      
      // Now use this to calculate light field and photon absorption
    
      double dz_sed = col->dz_sed[ws->num_sed_layers-1]; 
    
      double lightbot,meanlight,Kd_sed;
    
      photons = 0.0;
      yCfac = 0.0;
    
      for (w=0; w<ws->num_wave; w++){
      	Kd_sed = (ws->aw[w] + cellnum_sed * aA_s_CS_free_sed[w]) * dz_sed;   /* vertical attenuation */
	lightbot = lighttop_s[w] * exp(-Kd_sed);
	meanlight = (lighttop_s[w] - lightbot) / Kd_sed;
	photons += meanlight *  ws->wave[w] * aA_s_CS_free_sed[w];
      }
      for (w=ws->wPAR_bot; w<ws->wPAR_top; w++){
      	yCfac += diffaa(absorbance[w],rad)*max(lighttop_s[w] * ws->wave[w],0.0001);
	sumlight += max(lighttop_s[w] * ws->wave[w],0.0001);
      }
      y_sed0[ws->CS_free_yCfac_i] = yCfac/sumlight;
      y_sed0[ws->CS_free_kI_i] = photons * 1.0e20;
    }else{
      y_sed0[ws->CS_free_yCfac_i] = 0.0;
      y_sed0[ws->CS_free_kI_i] = 0.0;
    }
  }

  // First calculate absorption cross-section of MPB

  if (ws->MPB_N_i> -1){

    double Phy_N_sed = y_sed0[ws->MPB_N_i];

    if (Phy_N_sed > 1e-9){ // this is necessary since kI and yC calculated here are used next step in growth

      double Phy_Chl_sed = y_sed0[ws->MPB_Chl_i];
    
      rad = ws->rad_MPB;
      vol =  ws->vol_MPB;
      m = ws->m_MPB;
      
      double cellnum_sed =  Phy_N_sed / (m * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
      double cellChl_sed = Phy_Chl_sed / (vol * cellnum_sed); /* cellular Chl conc */ 
      
      for (w=0; w<ws->num_wave; w++){
	absorbance[w] = ws->yC_MPB[w] * cellChl_sed;	
      }
      aawave(rad, absorbance, aA_s_MPB_sed, ws->num_wave,ws->wave);
      
      // Now use this to calculate light field and photon absorption
    
      double dz_sed = col->dz_sed[ws->num_sed_layers-1]; 
    
      double lightbot,meanlight,Kd_sed;
    
      photons = 0.0;
      yCfac = 0.0;
    
      for (w=0; w<ws->num_wave; w++){
	  Kd_sed = (ws->aw[w] + cellnum_sed * aA_s_MPB_sed[w]) * dz_sed;   /* vertical attenuation */
	  lightbot = lighttop_s[w] * exp(-Kd_sed);
	  meanlight = (lighttop_s[w] - lightbot) / Kd_sed;
	  photons += meanlight *  ws->wave[w] * aA_s_MPB_sed[w];
      }
      for (w=ws->wPAR_bot; w<ws->wPAR_top; w++){
	  yCfac += diffaa(absorbance[w],rad)*max(lighttop_s[w] * ws->wave[w],0.0001);
	  sumlight += max(lighttop_s[w] * ws->wave[w],0.0001);
      }
      y_sed0[ws->MPB_yCfac_i] = yCfac/sumlight;
      y_sed0[ws->MPB_kI_i] = photons * 1.0e20;
    }else{
      y_sed0[ws->MPB_yCfac_i] = 0.0;
      y_sed0[ws->MPB_kI_i] = 0.0;
    }
  }
  
  d_free_2d(aA_s_l);  
  d_free_2d(aA_s_s);
  d_free_2d(aA_s_Tricho);
  d_free_2d(aA_s_MPB);
  d_free_2d(aA_s_PhyD);
  d_free_2d(aA_s_CS_free);
  d_free_1d(aA_s_MPB_sed);
  d_free_1d(aA_s_CS_free_sed);
  d_free_1d(aA_s_CS);
  d_free_1d(biolum);
  d_free_1d(fluor);
  d_free_1d(lighttop_s);
  d_free_1d(absorbance);

  double* Rrs = d_alloc_1d(ws->num_wave);
  double* Rrs_fine = d_alloc_1d(ws->num_sensor_waves);
  double** E_d_z_fine = d_alloc_2d(col->n_wc+1,ws->num_sensor_waves);
  
  double w_bot[ws->num_wave];
  double u_surf[ws->num_wave];
  double top,rs;

  /* Do ray-tracing and optical-depth weighting reflectance (Baird et al., 2020) */

  for (w = 0; w<ws->num_wave; w++) {
    w_bot[w] = 1.0;
    u_surf[w] = 0.0;
    for (n = 0; n<col->n_wc; n++) { // 0 is top
      top = w_bot[w];
      w_bot[w] = top * exp(-2.0 * Kd[w][n] * dz[n]);
      u_surf[w] += (bb[w][n] / (at[w][n] + bb[w][n])) * (top - w_bot[w]);
    }
  }

  // Add water column and benthic upwelling irradiance below water and convert to upwelling radiance in the air.

  for (w=0; w<ws->num_wave; w++) {
    u_surf[w] = u_surf[w] + u_bot[w] * w_bot[w];
    rs = u_surf[w] * (ws->g0 + ws->g1 * u_surf[w]);
    Rrs[w] = 0.52 * rs / (1.0 - 1.7 * rs);
  }

   // Interpolated Rrs and E_d_z onto Rrs_fine and E_d_z_fine respectively.

  interp1d(ws->wave, Rrs, ws->num_wave,ws->sensor_waves,Rrs_fine,ws->num_sensor_waves);

  // Need to do a layer by layer interpolation for E_d_z
  
  double tmp_in[ws->num_wave];
  double tmp_out[ws->num_sensor_waves];
  
  for (n = 0; n < col->n_wc; n++){ // 0 is top
    for (w=0; w<ws->num_wave; w++) {
      tmp_in[w] = E_d_z[w][n];
    }
    interp1d(ws->wave, tmp_in, ws->num_wave,ws->sensor_waves,tmp_out,ws->num_sensor_waves);
    
    for (w=0; w<ws->num_sensor_waves; w++) {
      E_d_z_fine[w][n] = tmp_out[w];
    }
  }
  
  double bot_light = 0.0;
  if (ws->SWR_bot_abs_i > -1){
    y_epi[ws->SWR_bot_abs_i] = 0.0;
    for (w=ws->wPAR_bot; w <= ws->wPAR_top; w++){
      bot_light += E_d_z[w][col->n_wc];
      y_epi[ws->SWR_bot_abs_i] += (1.0 - u_bot[w]) * E_d_z[w][col->n_wc];
    }
    y_epi[ws->SWR_bot_abs_i] = y_epi[ws->SWR_bot_abs_i] / bot_light;
  }
  
  // Now calculate tracer-list specified sensor response to downwelling light field.

  double tmmmp = 0.0;
  double denom = 0.0;
  
  for (i=0; i<ws->ntr_sp3; i++){
    for (n = 0; n<col->n_wc; n++) { // 0 is top
      tmmmp = 0.0;
      denom = 0.0;
      for (w = 0; w<ws->num_sensor_waves; w++) {
	tmmmp += ws->sp3_star[i][w] * E_d_z_fine[w][n]; // already multiplied by bandwidth.
	denom += ws->sp3_star[i][w] ;  // assuming bandwidth is 1.0; denominator could be done in init.
      }
      y[ws->ind_sp3[i]][n] = (tmmmp/denom);
    }
  }
  
  // Now calculate tracer-list specified sensor response to remote-sensing reflectance.
  
  for (i=0; i<ws->nepi_sp2; i++){
    tmmmp = 0.0;
    denom = 0.0;
    for (w = 0; w<ws->num_sensor_waves; w++) {
      tmmmp += ws->sp2_star[i][w] * Rrs_fine[w];
      denom += ws->sp2_star[i][w];  // denominator could be done in init.
    }
    y_epi[ws->ind_sp2[i]] = (tmmmp/denom);
  }
  
  /*************************************************************/
  /* Calculate spectrally-resolved AOPs using hydrolight model */
  /*************************************************************/
  
  int do_hydrolight = 0;
  int do_radtranx = 0;

  double hyd_zenith = 0.0;
  double hyd_azimuth = 0.0;
  
#ifdef INCLUDE_HYDROLIGHT

  // Hydrolight needs solar zenith to be above the horizon if using RADTRANX

  // if (col->b == ws->col_write && ems_solar_zenith < (M_PI/2.0 - 0.04) && ems_solar_zenith > 0.04) { // chosen column and daytime (col->b == ws->col_write && zenith <  (M_PI/2.0 - 0.2))
  if (col->b == col_write){
    if (ems_solar_zenith < (M_PI/2.0 - 0.04) && ems_solar_zenith > 0.04 && col_write_num == 0) {
      do_hydrolight = 1;
      do_radtranx = 1;
    }else{
      do_hydrolight = 0;
      do_radtranx  = 0;
    }
  }

  // if passing zenith to hydrolight, then do anytime.

  // if ((col->b == col_write) && (use_radtranx == 0)) {
  //  do_hydrolight = 1;
  //}

  // Need to make a decision about how deep to look.

  double maxcalcdepth = 200.0;
  
  long ndep = col->n_wc+1;  // Hydrolight calculating on the interfaces.
  
  // Default hydrolight has 10 x 24:

  int nphi = 24; // no. of azimuth angle calculations.
  int nmu = 10;  // no. of zenith angle calculations.

  int nn = nmu * nphi; // for index calculations.
  int phi,mu;

  double* positive_zc = d_alloc_1d(col->n_wc);

  double* a_z = d_alloc_1d(col->n_wc);
  double* bb_z = d_alloc_1d(col->n_wc);
  double* c_z = d_alloc_1d(col->n_wc);
  double* s_z = d_alloc_1d(col->n_wc);

  double**** L_d_hyd = d_alloc_4d(nphi,nmu,ndep,ws->num_wave);  // indexing is reversed.
  double**** L_u_hyd = d_alloc_4d(nphi,nmu,ndep,ws->num_wave);  // indexing is reversed.

  double*** L_sky = d_alloc_3d(nphi,nmu,ws->num_wave);  // indexing is reversed.
  double*** L_ref = d_alloc_3d(nphi,nmu,ws->num_wave);  // indexing is reversed.
  double*** L_wl = d_alloc_3d(nphi,nmu,ws->num_wave);  // indexing is reversed.

  double** E_d_hyd = d_alloc_2d(ndep,ws->num_wave);  // indexing is reversed.
  double** E_u_hyd = d_alloc_2d(ndep,ws->num_wave);  // indexing is reversed.

  double* zen = d_alloc_1d(nmu);
  double* azi = d_alloc_1d(nphi);

  double* zg = d_alloc_1d(ndep);
 
  double* rad_sky = d_alloc_1d(nn);
  double* rad_ref = d_alloc_1d(nn);
  double* rad_wl = d_alloc_1d(nn);

  double* rad_dn = d_alloc_1d(nn*ndep);  // indexing is reversed.
  double* rad_up = d_alloc_1d(nn*ndep);  // indexing is reversed.

  double irrad_sky = 0.0;
  double irrad_ref = 0.0;
  double irrad_wl = 0.0;

  double irrad_et = 0.0;

  double* Rrs_hyd = d_alloc_1d(ws->num_wave);
  double* Ed_sky = d_alloc_1d(ws->num_wave);
  double* Ed_ref = d_alloc_1d(ws->num_wave);
  double* Ed_wl = d_alloc_1d(ws->num_wave);

  double* irrad_dn = d_alloc_1d(ndep);
  double* irrad_up = d_alloc_1d(ndep);

  if (do_hydrolight){ // chosen column and daytime (col->b == ws->col_write && zenith <  (M_PI/2.0 - 0.2))

    // Hydrolight needs: (1) depth of the interfaces from the surface - zg[ndep], first being 0.
    //                   (2) depth of cell centres from surface  - positive_zc[ndep-1]

    zg[0] = 0.0;
    positive_zc[0] = col->dz[0]/2.0;
    for (n = 1; n<ndep; n++) { // 0 is top  
      zg[n] = zg[n-1] + col->dz[n-1];
    }
    for (n = 1; n<ndep-1; n++) { // 0 is top  
      positive_zc[n] = (zg[n] + zg[n+1])/2.0;
    }
    positive_zc[0] = 0.0;
    positive_zc[col->n_wc-1] += 100.0;
    
    printf("time units %s, e->t %e, ndoy %f, hr %f, solar_azimuth %e lon %e \n",ws->timeunits,e->t,ndoy,hr,-solar_azimuth*180.0/M_PI,longitude);
    
    // hydrolight does zenith relative to incoming solar. Adjust sunphi so it defines azimuth relative to E in clockwise direction.

    /************ wavelength loop **********************************/
    
    for (w = 0; w<ws->num_wave; w++) {
      
      for (n = 0; n<col->n_wc; n++) { // Choose array and remove clear water.
	a_z[n] = at[w][n] - ws->aw[w];
	bb_z[n] = bb[w][n] - ws->bw[w] * 0.5;
	c_z[n] = at[w][n] + bt[w][n] - ws->bw[w];
	s_z[n] = 0.0; // blt[w][n]/ws->bandwidth[w]; // (blt[w][n] + flt[w][n])/ws->bandwidth[w]; // W nm-1 m-3
      }

      printf("Starting hydrolight calculations for wavelength %e at time %e ems zenith %e azimuth %e \n",ws->wave[w],e->t,zenith*180.0/M_PI,(azimuth>0.0?azimuth:azimuth+2.0*M_PI)*180.0/M_PI);

      set_run(ws->wave[w]);
      
      // Try setting up non-wavelength hydrolight calls outside loop.
      if (do_radtranx){
	set_geo(3,zenith,ndoy, hr, -solar_azimuth*180.0/M_PI, latitude, longitude); // provide nday, hr, lat, lon: get zenith, azimuth.
      }else{
	// wants azimuth relative to wind direction ?????
	set_geo(2,zenith*180.0/M_PI,ndoy, hr, (azimuth>0.0?azimuth:azimuth+2.0*M_PI)*180.0/M_PI, latitude, longitude);     
      }

      // mean windspeed references a file with a int between 0 - 15!
      
      set_atm(wind_speed,tcld,mslp/3386.3752577878,wind_speed,-1.,-1.,-1.,-1.,-1.); // incl. Pa -> mm Hg conversion.
      
      set_water(y[ws->temp_i][0],n_stw[w],y[ws->salt_i][0]);

      set_bottom(1, u_bot[w]);    // (0 - no reflect, 1 - lambertian reflect | bottom reflectance)

      set_depth_out(ndep, zg);
      
      //if (ws->use_radtranx==0){
	//set_sky_irrad(E_d_z[w][0]); // Just above the water - overwrites RADTRANX.
      //}else{
	// no call
      //}

      hydrolight(col->n_wc, ws->wave[w], positive_zc, a_z, c_z, bb_z, s_z);

      //if (1==0)
      // solar_irrad(&irrad_et);

      // Question: does it know what ifladsky is based on previous calls?

      double moon_TA = Ak[w] / ws->bandwidth[w];
      
      if (do_radtranx==0){
	copy_output(nmu, nphi, ndep, zen, azi, rad_sky, rad_ref, rad_wl, rad_dn, rad_up,&irrad_sky,&irrad_ref,&irrad_wl,irrad_dn,irrad_up,&hyd_zenith,&hyd_azimuth,moon_TA);
      }else{
	copy_output_ns(nmu, nphi, ndep, zen, azi, rad_sky, rad_ref, rad_wl, rad_dn, rad_up,&irrad_sky,&irrad_ref,&irrad_wl,irrad_dn,irrad_up,&hyd_zenith,&hyd_azimuth);
      }
      
      Rrs_hyd[w] = (irrad_wl/irrad_dn[0])/2.0/M_PI;  // Hydrolight code uses: Rrs = radw/Ed(0):

      hyd_zenith  = hyd_zenith*2.0/M_PI;
      hyd_azimuth = hyd_azimuth*2.0/M_PI;

      Ed_sky[w] = irrad_sky;
      Ed_ref[w] = irrad_ref;
      Ed_wl[w] = irrad_wl;

      // Put irradiances into arrays/, E: is L in sr-1. If so something may need to be added
      // From hydrolight both E and L are per nm.

      for (n = 0; n < ndep; n++) { // 0 is top
	E_d_hyd[w][n] = irrad_dn[n]; // below water.
	E_u_hyd[w][n] = irrad_up[n];
      }
      
      for (phi = 0; phi<nphi; phi++) {
	for (mu = 0; mu<nmu; mu++) {
	  L_sky[w][mu][phi] = rad_sky[mu + phi*nmu];
	  L_ref[w][mu][phi] = rad_ref[mu + phi*nmu];
	  L_wl[w][mu][phi] =   rad_wl[mu + phi*nmu];
	  
	  for (n = 0; n < ndep; n++) { // L [W/m2/nm/sr]
	    L_d_hyd[w][n][mu][phi] = rad_dn[mu + phi*nmu + n*nn];
	    L_u_hyd[w][n][mu][phi] = rad_up[mu + phi*nmu + n*nn];
	  }
	}
      }
      printf("Finished hydrolight loop for wavelength %e nm \n",ws->wave[w]);
    }
    
    /************ wavelength loop **********************************/

    }else{ // not chosen column or night

    hyd_zenith =  0.0;
    hyd_azimuth = 0.0;

    for (w = 0; w<ws->num_wave; w++) {
     
      for (n = 0; n < ndep; n++) { // 0 is top
	E_d_hyd[w][n] = 0.0;
	E_u_hyd[w][n] = 0.0;
      }
      for (phi = 0; phi<nphi; phi++) {
	for (mu = 0; mu<nmu; mu++) {
	  L_sky[w][mu][phi] = 0.0;
	  L_ref[w][mu][phi] = 0.0;
	  L_wl[w][mu][phi] = 0.0;
	  for (n = 0; n < ndep; n++) { // 0 is top
	    L_d_hyd[w][n][mu][phi] = 0.0;
	    L_u_hyd[w][n][mu][phi] = 0.0;
	  }
	}
      }
    }
  }

#endif
    
  // Use sensor responses in satellite algorithms.
 
  if (ws->OC4Me_SR_i > -1){
      OC4Me(y_epi[ws->Sentinel_3B_Band3_i], y_epi[ws->Sentinel_3B_Band4_i],y_epi[ws->Sentinel_3B_Band5_i],y_epi[ws->Sentinel_3B_Band6_i],&y_epi[ws->OC4Me_SR_i]);
  }
  if (ws->Hue_SR_i > -1){    
    Hue3B(y_epi[ws->Sentinel_3B_Band3_i], y_epi[ws->Sentinel_3B_Band4_i], y_epi[ws->Sentinel_3B_Band6_i], y_epi[ws->Sentinel_3B_Band8_i],y_epi[ws->Sentinel_3B_Band11_i],&y_epi[ws->Hue_SR_i]);
  }

  if (ws->OC3M_SR_i > -1){
    OC3M(y_epi[ws->MODIS_Band2_i],y_epi[ws->MODIS_Band4_i],y_epi[ws->MODIS_Band6_i],&y_epi[ws->OC3M_SR_i]);  
  }   
  
  if (ws->nFLH_SR_i > -1){
    nFLH_modis(y_epi[ws->MODIS_Band9_i],y_epi[ws->MODIS_Band10_i],y_epi[ws->MODIS_Band11_i],&y_epi[ws->nFLH_SR_i]); 
  }
  
  if (ws->TSSM_SR_i > -1)
    TSSM(y_epi[ws->MODIS_Band9_i],&y_epi[ws->TSSM_SR_i]);
  
  if (ws->KD490M_SR_i > -1)
    KD490M(y_epi[ws->MODIS_Band4_i],y_epi[ws->MODIS_Band6_i],&y_epi[ws->KD490M_SR_i]);
      
  if (ws->OC3V_SR_i > -1){
    OC3V(y_epi[ws->VIIRS_Band2_i],y_epi[ws->VIIRS_Band3_i],y_epi[ws->VIIRS_Band4_i],&y_epi[ws->OC3V_SR_i]);  
  }
  
  // Secchi depth calculation.

  double z_secchi = 0.0;
  double E_secchi = 1.0;
  double E_secchi_tmp = 0.0;
  double tmp_adlen;
  
  if (ws->Secchi_i > -1){
    
    for (n = 0; n<col->n_wc; n++) { // 0 is top
      if (E_secchi > exp(-1.0)){
	tmp_adlen = at[ws->w490][n] * sqrt(1.0 + (ws->gi * costhetaw[ws->w490] - ws->gii) * bt[ws->w490][n] / at[ws->w490][n]) / costhetaw[ws->w490];
	E_secchi_tmp = E_secchi;
	E_secchi = E_secchi * exp(- tmp_adlen * dz[n]);
	if (E_secchi < exp(-1.0)){ // this is the step the disk disappears from view
	  z_secchi = z_secchi + log(exp(-1.0)/E_secchi_tmp) / (- tmp_adlen);
	}else{                           // otherwise still in view.
	  z_secchi = z_secchi + dz[n];
	}
      }
    }
    y_epi[ws->Secchi_i] = z_secchi;
  }

  /* Output spectrally-integrated variables */
  
  double energy;

  if (ws->PAR_z_i > -1){ /* i.e. outputting */
  
    for (n = 0; n<col->n_wc; n++) { // 0 is top

      energy = 0.0;

      for (w=ws->wPAR_bot; w <= ws->wPAR_top; w++){
	energy +=  E_d_z[w][n] * 8.359335857479461e-09 * ws->wave[w];
      }	
      y[ws->PAR_z_i][n] = energy;
    }
  }

  if (ws->PAR_i > -1){ /* i.e. outputting */
  
    for (n = 0; n<col->n_wc; n++) { // 0 is top

      energy = 0.0;

      for (w=ws->wPAR_bot; w <= ws->wPAR_top; w++){
	energy +=  E_d[w][n] * 8.359335857479461e-09 * ws->wave[w];
      }	
      y[ws->PAR_i][n] = energy;
    }
  }

#ifdef HYDROLIGHT
  if (ws->PAR_hyd_i > -1){ /* i.e. outputting */
  
    for (n = 0; n<col->n_wc; n++) { // 0 is top

      energy = 0.0;
      for (w=ws->wPAR_bot; w <= ws->wPAR_top; w++){
	energy +=  E_d_hyd[w][n] * 8.359335857479461e-09 * ws->wave[w] * ws->bandwidth[w];
      }	
      y[ws->PAR_hyd_i][n] = energy;
    }
  }
#endif

  if (ws->K_heat_i > -1){ /* i.e. outputting */
    if (E_d[ws->w490][0] > 0.0){

      int thickn = 0;
      double thickdz = dz[0];
      while (thickdz < ws->K_heat_minthickness && thickn < col->n_wc-1){
	thickn += 1;
	thickdz += dz[thickn];
      }
     	
      double energy1,delenergy1;
      for (n = thickn; n<col->n_wc; n++) { // 0 is top
	energy1 = 0.0;
	delenergy1 = 0.0;
	for (w=0; w <= ws->wPAR_top; w++){
	  energy1 += E_d_z[w][n];
	  delenergy1 += (E_d_z[w][n]-E_d_z[w][n+1]);
	}
	if (delenergy1 > 0.0001 && fabs(energy1-delenergy1) > 1.0e-5) { 
	  y[ws->K_heat_i][n] = - log((energy1-delenergy1)/energy1) / dz[n];
	}else{
	  y[ws->K_heat_i][n] = 0.04;
	}
      }
      if (thickn > 0 && thickn < col->n_wc){ // this involves writing over thickn so that two the same.	
	energy1 = 0.0;
	delenergy1 = 0.0;
	for (w=0; w <= ws->wPAR_top; w++){
	  energy1 += E_d_z[w][0];
	  delenergy1 += (E_d_z[w][0]-E_d_z[w][thickn+1]);
	}
	if (delenergy1 > 0.0001 && fabs(energy1-delenergy1) > 1.0e-5) {
	  double avKz = - log((energy1-delenergy1)/energy1) / thickdz;
	  for (n = 0;n < thickn + 1; n++){ 
	    y[ws->K_heat_i][n] = avKz;
	  }
	}else{
	  y[ws->K_heat_i][n] = 0.04;
	}
      }
    }else{
      y[ws->K_heat_i][n] = 0.04;
    }
  }

  if (ws->Kd_490_i > -1){ /* i.e. outputting */
    for (n = 0; n<col->n_wc; n++) { // 0 is top 
      y[ws->Kd_490_i][n] = Kd[ws->w490][n];
    }
  }

  if (ws->Kd_PAR_i > -1){ /* i.e. outputting - normalise by light intensity */
    for (n = 0; n<col->n_wc; n++) { // 0 is top
      double kdtemp = 0.0;
      double E_d_z_sum = 0.0;
      for (w=ws->wPAR_bot; w <= ws->wPAR_top; w++){
	  kdtemp += E_d_z[w][n] * Kd[w][n];
	  E_d_z_sum += E_d_z[w][n];
      } 
      y[ws->Kd_PAR_i][n] = kdtemp/E_d_z_sum;
    }
  }

  if (ws->at_440_i > -1){ // i.e. outputting
    for (n = 0; n<col->n_wc; n++) { // 0 is top
      y[ws->at_440_i][n] = at[ws->w440][n];
    }
  }

  if (ws->bt_550_i > -1){ // i.e. outputting 
    for (n = 0; n<col->n_wc; n++) { // 0 is top 
      y[ws->bt_550_i][n] = bt[ws->w550][n];
    }
  }

  if (ws->bb_590_i > -1){ // i.e. outputting 
    for (n = 0; n<col->n_wc; n++) { // 0 is top 
      y[ws->bb_590_i][n] = bb[ws->w590][n];
    }
  }

  if (ws->Turbidity_i > -1){ // i.e. outputting 
    for (n = 0; n<col->n_wc; n++) { // 0 is top 
      y[ws->Turbidity_i][n] = 47.02 * (bb[ws->w590][n] - 0.5 * ws->bw[ws->w590]) + 0.13;
    }
  }

  if (ws->bb_700_i > -1){ // i.e. outputting 
    for (n = 0; n<col->n_wc; n++) { // 0 is top 
      y[ws->bb_700_i][n] = bb[ws->w700][n];
    }
  }

  if (ws->BLP_i > -1){ // i.e. outputting
    double tmpsum;
    for (n = 0; n<col->n_wc; n++) { // 0 is top
      tmpsum = 0.0;
      for (w=0; w < ws->num_wave; w++){
	tmpsum += blt[w][n];
      }
      y[ws->BLP_i][n] = tmpsum;
    }
  }
  
  if (ws->PFL_i > -1){ // i.e. outputting
    double tmpsum1;
    for (n = 0; n<col->n_wc; n++) { // 0 is top
      tmpsum1 = 0.0;
      for (w=0; w < ws->num_wave; w++){
	tmpsum1 += flt[w][n];
      }
      y[ws->PFL_i][n] = tmpsum1;
    }
  }

  // Interpolate reflectances from resolved wave centres onto benthic tracers named R_XXX.

  int w2;
  for (w2=0; w2<ws->num_rrs_waves; w2++){
    if (ws->wave[ws->wXXX_i[w2]] == ws->rrs_wave[w2]){  // No interpolation necessary
      y_epi[ws->RRS_i[w2]] = Rrs[ws->wXXX_i[w2]];
    }else{
      y_epi[ws->RRS_i[w2]] = (Rrs[ws->wXXX_i[w2]] * (ws->rrs_wave[w2] - ws->wave[ws->wXXX_i[w2]-1]) + Rrs[ws->wXXX_i[w2]-1] * (ws->wave[ws->wXXX_i[w2]] - ws->rrs_wave[w2])) / (ws->wave[ws->wXXX_i[w2]] - ws->wave[ws->wXXX_i[w2]-1]);
    }
  }

  // Interpolate downwelling irradiance at the top of the cell from resolved wave centre onto 3D tracers named Ed_YYY.
  
  if (E_d[ws->w490][0] > 0.0){
    for (w2=0; w2<ws->num_ed_waves; w2++){  
      if (ws->wave[ws->wYYY_i[w2]] == ws->ed_wave[w2]){  // No interpolation necessary
	for (n = 0; n<col->n_wc; n++) {
	  y[ws->Ed_i[w2]][n] = E_d_z[ws->wYYY_i[w2]][n]/ws->bandwidth[ws->wYYY_i[w2]];
	} 
      }else{
	for (n = 0; n<col->n_wc; n++) {
	  y[ws->Ed_i[w2]][n] = (E_d_z[ws->wYYY_i[w2]][n] * (ws->ed_wave[w2] - ws->wave[ws->wYYY_i[w2]-1]) + E_d_z[ws->wYYY_i[w2]-1][n] * (ws->wave[ws->wYYY_i[w2]] - ws->ed_wave[w2])) / (ws->wave[ws->wYYY_i[w2]] - ws->wave[ws->wYYY_i[w2]-1]) / ws->bandwidth[ws->wYYY_i[w2]];
	}
      }
    }
  }else{   // No light at 490 at surface
    for (n = 0; n<col->n_wc; n++) { 
      for (w2=0; w2<ws->num_ed_waves; w2++){
	y[ws->Ed_i[w2]][n] = 0.0;
      }
    }
  }
  
  int do_write = 0;
  if (col_write_num > -1) {
    do_write = 1;
  }
  
#ifdef INCLUDE_HYDROLIGHT

  /* Compare Mobely and Beer solutions */

  if (do_hydrolight){
    if (col->b == col_write && zenith < (M_PI/2.0 - 0.02)){
      for (w = 0; w<ws->num_wave; w++) {
	for (n = 0; n<col->n_wc; n++) { // 0 is top
	  printf("  %f nm layer %d, E_d_ray = %e E_d_mobely %e E_u_mobely %e \n", ws->wave[w], n, E_d_z[w][n]/ws->bandwidth[w],E_d_hyd[w][n],E_u_hyd[w][n]);
	}
	printf("Remote-sensing reflectance at %f nm, Rrs_Baird %e, Rrs_Mobley %e \n", ws->wave[w], Rrs[w],Rrs_hyd[w]);
      }
    }
  }

#endif

  if (do_write==1) {

    // calculations that update variable put into the column output, but don't change anything in the tracer list.

    #ifdef INCLUDE_BIOLUMINESCENCE
    double* Secchi_colour = d_alloc_1d(ws->num_wave);
    Secchi_colour_calc(ws->num_wave,col->n_wc,Kd,dz,Secchi_colour);
    #endif
      
    int ncid,varid,recid;
    size_t reclen;

    char ooutput[MAXSTRLEN];
    char* output_path = ginterface_get_output_path();
    sprintf(ooutput,"%soptical_column_out_%d.nc",output_path,col_write_num);
    nc_open(ooutput,NC_WRITE, &ncid);
 
    nc_inq_dimid(ncid,"record",&recid);
    
    nc_inq_dimlen(ncid,recid,&reclen);
    
    varid = ncw_var_id(ncid,"t");

    /* Here call tm_change_time_units() to convert the model time to the output_tunits are */

    double newt = e->t + e->dt; 
    tm_change_time_units(ws->timeunits, ws->outputtimeunits, &newt, 1);
    nc_put_var1_double(ncid,varid,&reclen,&newt);

    varid = ncw_var_id(ncid,"ems_solar_zenith");
    nc_put_var1_double(ncid,varid,&reclen,&zenith);

    varid = ncw_var_id(ncid,"ems_solar_azimuth");
    nc_put_var1_double(ncid,varid,&reclen,&solar_azimuth);

    varid = ncw_var_id(ncid,"hyd_solar_zenith");
    nc_put_var1_double(ncid,varid,&reclen,&hyd_zenith);

    varid = ncw_var_id(ncid,"hyd_solar_azimuth");
    nc_put_var1_double(ncid,varid,&reclen,&hyd_azimuth);
    
    if (ws->Sentinel_3B_Band3_i){
      varid = ncw_var_id(ncid,"Sentinel_3B_Band3");
      nc_put_var1_double(ncid,varid,&reclen,&y_epi[ws->Sentinel_3B_Band3_i]);
    }
    if (ws->Sentinel_3B_Band4_i){
      varid = ncw_var_id(ncid,"Sentinel_3B_Band4");
      nc_put_var1_double(ncid,varid,&reclen,&y_epi[ws->Sentinel_3B_Band4_i]);
    }
    if (ws->Sentinel_3B_Band5_i){
      varid = ncw_var_id(ncid,"Sentinel_3B_Band5");
      nc_put_var1_double(ncid,varid,&reclen,&y_epi[ws->Sentinel_3B_Band5_i]);
    }
    if (ws->Sentinel_3B_Band6_i){
      varid = ncw_var_id(ncid,"Sentinel_3B_Band6");
      nc_put_var1_double(ncid,varid,&reclen,&y_epi[ws->Sentinel_3B_Band6_i]);
    }
    if (ws->Sentinel_3B_Band8_i){
      varid = ncw_var_id(ncid,"Sentinel_3B_Band8");
      nc_put_var1_double(ncid,varid,&reclen,&y_epi[ws->Sentinel_3B_Band8_i]);
    }
    if (ws->Sentinel_3B_Band11_i){
      varid = ncw_var_id(ncid,"Sentinel_3B_Band11");
      nc_put_var1_double(ncid,varid,&reclen,&y_epi[ws->Sentinel_3B_Band11_i]);
    }

    if (ws->Moonlight_i > -1){

      varid = ncw_var_id(ncid,"ems_lunar_zenith");
      nc_put_var1_double(ncid,varid,&reclen,&lunar_zenith);

      varid = ncw_var_id(ncid,"ems_lunar_azimuth");
      nc_put_var1_double(ncid,varid,&reclen,&lunar_azimuth);

      varid = ncw_var_id(ncid,"moon_phase");
      nc_put_var1_double(ncid,varid,&reclen,&moon_phase);

      varid = ncw_var_id(ncid,"moon_fulldisk");
      nc_put_var1_double(ncid,varid,&reclen,&moon_fulldisk);

      varid = ncw_var_id(ncid,"ems_lunar_albedo");
      nc_put_var1_double(ncid,varid,&reclen,&lunar_albedo);

      varid = ncw_var_id(ncid,"ems_lunar_airtrans");
      nc_put_var1_double(ncid,varid,&reclen,&lunar_airtransmission);

      varid = ncw_var_id(ncid,"moonlight");
      nc_put_var1_double(ncid,varid,&reclen,&y_epi[ws->Moonlight_i]);
      
    }

#ifdef INCLUDE_HYDROLIGHT
    if (col_write_num == 0){
      varid = ncw_var_id(ncid,"hydrolight_zenith");
      nc_put_var1_double(ncid,varid,&reclen,&hyd_zenith);

      varid = ncw_var_id(ncid,"hydrolight_azimuth");
      nc_put_var1_double(ncid,varid,&reclen,&hyd_azimuth);
    }
#endif
  
    double etaa = ginterface_get_eta(col->model,col->b); // why is this one time step old ???
    varid = ncw_var_id(ncid,"eta");
    nc_put_var1_double(ncid,varid,&reclen,&etaa);
    
    varid = ncw_var_id(ncid,"wind_speed");
    nc_put_var1_double(ncid,varid,&reclen,&wind_speed);

    varid = ncw_var_id(ncid,"cloud_cover");
    nc_put_var1_double(ncid,varid,&reclen,&tcld);

    varid = ncw_var_id(ncid,"mslp");
    nc_put_var1_double(ncid,varid,&reclen,&mslp);

    varid = ncw_var_id(ncid,"rh2m");
    nc_put_var1_double(ncid,varid,&reclen,&rh2m);

    varid = ncw_var_id(ncid,"col_index");
    nc_put_var1_int(ncid,varid,&reclen,&col->b);

    varid = ncw_var_id(ncid,"use_radtranx");
    nc_put_var1_int(ncid,varid,&reclen,&do_radtranx);

    int countt;

    // need to put hydrolight calculations on vertical grid from surface.

    double* varr = d_alloc_1d(ws->num_wc_layers+1);

    int out_top = ws->num_wc_layers - col->topk_wc;
    int out_bot = ws->num_wc_layers - col->botk_wc + 1;

    for (n = 0; n < out_top; n++) {
      varr[n] = 1.0e35;
    }
    varr[out_top-1] = 0.0;
    for (n = out_top; n < out_bot; n++) {
      varr[n] = varr[n-1] + col->dz[n-out_top];      
    }
    for (n = out_bot; n < ws->num_wc_layers + 1; n++) {
      varr[n] = -9999.0;
    }

    varid = ncw_var_id(ncid,"depth");
    const size_t start[2] = {reclen,0};
    const size_t count[2] = {1,ws->num_wc_layers+1};

    nc_put_vara_double(ncid,varid,start,count,varr);

    const size_t start2[3] = {reclen,0,0};
    const size_t count2[3] = {1,ws->num_wave,ws->num_wc_layers+1};

    double* E_d_z_tmp = d_alloc_1d((ws->num_wc_layers+1) * ws->num_wave);
    
    countt = 0;

    for (w = 0; w<ws->num_wave; w++) {
      for (n = 0; n < out_top; n++) {
	E_d_z_tmp[countt] = 1.0e35;
	countt = countt+1;
      }
      for (n = out_top; n < out_bot; n++) {
	E_d_z_tmp[countt] = E_d_z[w][n-out_top]/ws->bandwidth[w];
	countt = countt+1;
      }
      for (n =  out_bot; n < ws->num_wc_layers+1; n++) {
	E_d_z_tmp[countt] = -9999.0;
	countt = countt+1;
      }
    }

    varid = ncw_var_id(ncid,"E_d_ray");
    nc_put_vara_double(ncid,varid,start2,count2,E_d_z_tmp);

    #ifdef INCLUDE_BIOLUMINESCENCE
    
    double* E_u_z_bio_tmp = d_alloc_1d((ws->num_wc_layers+1) * ws->num_wave);
    
    countt = 0;

    for (w = 0; w<ws->num_wave; w++) {
      for (n = 0; n < out_top; n++) {
	E_u_z_bio_tmp[countt] = 1.0e35;
	countt = countt+1;
      }
      for (n = out_top; n < out_bot; n++) {
	E_u_z_bio_tmp[countt] = E_u_z_bio[w][n-out_top]/ws->bandwidth[w];
	countt = countt+1;
      }
      for (n =  out_bot; n < ws->num_wc_layers+1; n++) {
	E_u_z_bio_tmp[countt] = -9999.0;
	countt = countt+1;
      }
    }

    varid = ncw_var_id(ncid,"E_u_z_bio");
    nc_put_vara_double(ncid,varid,start2,count2,E_u_z_bio_tmp);

    #endif
    
    // now on wave_centre.

    out_top = ws->num_wc_layers - col->topk_wc-1;
    out_bot = ws->num_wc_layers - col->botk_wc;

    // Now need to turn 2d array into 1 d array.

    double* E_s_tmp = d_alloc_1d(ws->num_wc_layers * ws->num_wave);
    double* at_tmp = d_alloc_1d(ws->num_wc_layers * ws->num_wave);
    double* bt_tmp = d_alloc_1d(ws->num_wc_layers * ws->num_wave);
    double* bb_tmp = d_alloc_1d(ws->num_wc_layers * ws->num_wave);
    double* blt_tmp = d_alloc_1d(ws->num_wc_layers * ws->num_wave);
    double* flt_tmp = d_alloc_1d(ws->num_wc_layers * ws->num_wave);

    countt=0;
    for (w = 0; w<ws->num_wave; w++) {
      for (n = 0; n < out_top; n++) {
	E_s_tmp[countt] = 1.0e35;
	at_tmp[countt] = 1.0e35;
	bt_tmp[countt] = 1.0e35;
	bb_tmp[countt] = 1.0e35;
	blt_tmp[countt] = 1.0e35;
	flt_tmp[countt] = 1.0e35;	
	countt = countt+1;
      }
      for (n = out_top; n < out_bot; n++) {
	E_s_tmp[countt] = E_s[w][n-out_top]/ws->bandwidth[w];
	at_tmp[countt] = at[w][n-out_top];
	bt_tmp[countt] = bt[w][n-out_top];
	bb_tmp[countt] = bb[w][n-out_top];
	blt_tmp[countt] = blt[w][n-out_top]/ws->bandwidth[w];
	flt_tmp[countt] = flt[w][n-out_top]/ws->bandwidth[w];
	countt = countt+1;
      }
      for (n = out_bot; n < ws->num_wc_layers; n++) {
	E_s_tmp[countt] = -9999.0;
	at_tmp[countt] = -9999.0;
	bt_tmp[countt] = -9999.0;
	bb_tmp[countt] = -9999.0;
	blt_tmp[countt] = -9999.0;
	flt_tmp[countt] = -9999.0;
	countt = countt+1;
      }
    }

    const size_t start6[3] = {reclen,0,0};
    const size_t count6[3] = {1,ws->num_wave,ws->num_wc_layers};

    varid = ncw_var_id(ncid,"E_o_ray");
    nc_put_vara_double(ncid,varid,start6,count6,E_s_tmp);

    varid = ncw_var_id(ncid,"at");
    nc_put_vara_double(ncid,varid,start6,count6,at_tmp);

    varid = ncw_var_id(ncid,"bt");
    nc_put_vara_double(ncid,varid,start6,count6,bt_tmp);

    varid = ncw_var_id(ncid,"bb");
    nc_put_vara_double(ncid,varid,start6,count6,bb_tmp);

    varid = ncw_var_id(ncid,"blt");
    nc_put_vara_double(ncid,varid,start6,count6,blt_tmp);

    varid = ncw_var_id(ncid,"flt");
    nc_put_vara_double(ncid,varid,start6,count6,flt_tmp);

    varid = ncw_var_id(ncid,"Bottom_reflectance");
    
    const size_t start5[3] = {reclen,0,0};
    const size_t count5[3] = {1,ws->num_wave};
    const size_t count5a[3] = {1,ws->num_sensor_waves};

    nc_put_vara_double(ncid,varid,start5,count5,u_bot);

    varid = ncw_var_id(ncid,"Rrs_odw");
    nc_put_vara_double(ncid,varid,start5,count5,Rrs);
    
    varid = ncw_var_id(ncid,"Rrs_fine");
    nc_put_vara_double(ncid,varid,start5,count5a,Rrs_fine);

    #ifdef INCLUDE_BIOLUMINESCENCE
    varid = ncw_var_id(ncid,"Secchi_colour");
    nc_put_vara_double(ncid,varid,start5,count5,Secchi_colour);
    #endif
    
    varid = ncw_var_id(ncid,"n_stw");
    nc_put_vara_double(ncid,varid,start5,count5,n_stw);

    if (ws->Moonlight_i > -1){
      varid = ncw_var_id(ncid,"Moonlight_SR");
      nc_put_vara_double(ncid,varid,start5,count5,moonlight_SR);
    }

#ifdef INCLUDE_HYDROLIGHT

    double* E_d_hyd_tmp = d_alloc_1d((ws->num_wc_layers+1) * ws->num_wave);
    double* E_u_hyd_tmp = d_alloc_1d((ws->num_wc_layers+1) * ws->num_wave);

    double* L_u_tmp = d_alloc_1d(nn*ws->num_wave*(ws->num_wc_layers+1));
    double* L_d_tmp = d_alloc_1d(nn*ws->num_wave*(ws->num_wc_layers+1));

    double* L_wl_tmp = d_alloc_1d(nn*ws->num_wave);
    double* L_ref_tmp = d_alloc_1d(nn*ws->num_wave);
    double* L_sky_tmp = d_alloc_1d(nn*ws->num_wave);

    if (col_write_num == 0){

    varid = ncw_var_id(ncid,"Rrs_hyd");
    nc_put_vara_double(ncid,varid,start5,count5,Rrs_hyd);

    varid = ncw_var_id(ncid,"E_sky_hyd");
    nc_put_vara_double(ncid,varid,start5,count5,Ed_sky);

    varid = ncw_var_id(ncid,"E_ref_hyd");
    nc_put_vara_double(ncid,varid,start5,count5,Ed_ref);

    varid = ncw_var_id(ncid,"E_wl_hyd");
    nc_put_vara_double(ncid,varid,start5,count5,Ed_wl);

    // revert back to on the interfaces.

    out_top = ws->num_wc_layers - col->topk_wc;
    out_bot = ws->num_wc_layers - col->botk_wc+1; // change from -1, but not checked.
    
    countt = 0;
    for (w = 0; w<ws->num_wave; w++) {
      for (n = 0; n < out_top; n++) {
	E_d_hyd_tmp[countt] = 1.0e35;
	E_u_hyd_tmp[countt] = 1.0e35;
	countt = countt+1;
      }
      for (n = out_top; n < out_bot; n++) {
	E_d_hyd_tmp[countt] = E_d_hyd[w][n-out_top];
	E_u_hyd_tmp[countt] = E_u_hyd[w][n-out_top];
	countt = countt+1;
      }
      for (n = out_bot; n < ws->num_wc_layers+1; n++){
	E_d_hyd_tmp[countt] = -9999.0;
	E_u_hyd_tmp[countt] = -9999.0;
	countt = countt+1;
      }
    }

    varid = ncw_var_id(ncid,"E_d_hyd");
    
    nc_put_vara_double(ncid,varid,start2,count2,E_d_hyd_tmp);

    varid = ncw_var_id(ncid,"E_u_hyd");
    nc_put_vara_double(ncid,varid,start2,count2,E_u_hyd_tmp);

    countt = 0;
    for (w = 0; w<ws->num_wave; w++) {
      for (n = 0; n < out_top; n++) {
	for (mu = 0; mu<nmu; mu++) {
	  for (phi = 0; phi<nphi; phi++) {
	    L_u_tmp[countt] = 1.0e35;
	    L_d_tmp[countt] = 1.0e35;
	    countt = countt+1;
	  }
	}
      }

      for (n = out_top; n < out_bot; n++) {
	for (mu = 0; mu<nmu; mu++) {
	  for (phi = 0; phi<nphi; phi++) {
	    L_u_tmp[countt] = L_u_hyd[w][n-out_top][mu][phi];
	    L_d_tmp[countt] = L_d_hyd[w][n-out_top][mu][phi];
	    countt = countt+1;
	  }
	}
      }

      for (n = out_bot; n < ws->num_wc_layers+1; n++){
	for (mu = 0; mu<nmu; mu++) {
	  for (phi = 0; phi<nphi; phi++) {
	    L_u_tmp[countt] = -9999.0;
	    L_d_tmp[countt] = -9999.0;
	    countt = countt+1;
	  }
	}
      }
    }

    varid = ncw_var_id(ncid,"L_u");
    const size_t start3[5] = {reclen,0,0,0,0};
    const size_t count3[5] = {1,ws->num_wave,ws->num_wc_layers+1,10,24};
    
    nc_put_vara_double(ncid,varid,start3,count3,L_u_tmp);

    varid = ncw_var_id(ncid,"L_d");
    nc_put_vara_double(ncid,varid,start3,count3,L_d_tmp);
    
    countt = 0;

    for (w = 0; w<ws->num_wave; w++) {
      for (mu = 0; mu<nmu; mu++) {
	for (phi = 0; phi<nphi; phi++) {
	  L_wl_tmp[countt] = L_wl[w][mu][phi];
	  L_sky_tmp[countt] = L_sky[w][mu][phi];
	  L_ref_tmp[countt] = L_ref[w][mu][phi];
	  countt = countt+1;
	}
      }
    }
    
    const size_t start4[4] = {reclen,0,0,0};
    const size_t count4[4] = {1,ws->num_wave,10,24};

    varid = ncw_var_id(ncid,"L_wl");
    nc_put_vara_double(ncid,varid,start4,count4,L_wl_tmp);

    varid = ncw_var_id(ncid,"L_sky");
    nc_put_vara_double(ncid,varid,start4,count4,L_sky_tmp);

    varid = ncw_var_id(ncid,"L_ref");
    nc_put_vara_double(ncid,varid,start4,count4,L_ref_tmp);

    }
    
#endif

    // Add BGC state variables - need to put back onto Z-grid.

    const size_t start15[2] = {reclen,0};
    const size_t count15[2] = {1,ws->num_wc_layers};

    double* varr1 = d_alloc_1d(ws->num_wc_layers);

    for (n = 0; n < out_top; n++) {
      varr1[n] = 1.0e35;
    }
    for (n = out_top; n < out_bot; n++) {
      varr1[n] = y[ws->temp_i][n-out_top];      
    }
    for (n = out_bot; n < ws->num_wc_layers; n++) {
      varr1[n] = -9999.0;
    }

    varid = ncw_var_id(ncid,"temp");
    nc_put_vara_double(ncid,varid,start15,count15,varr1);

    if (ws->K_heat_i > -1){
      varid = ncw_var_id(ncid,"K_heat");
      for (n = out_top; n < out_bot; n++) {
	varr1[n] = y[ws->K_heat_i][n-out_top];      
      }
      nc_put_vara_double(ncid,varid,start15,count15,varr1);
    }

    if (ws->PhyL_Chl_i > -1){
      varid = ncw_var_id(ncid,"PhyL_Chl");
      for (n = out_top; n < out_bot; n++) {
	varr1[n] = y[ws->PhyL_Chl_i][n-out_top];      
      }
      nc_put_vara_double(ncid,varid,start15,count15,varr1);
      
      varid = ncw_var_id(ncid,"PhyL_N");
      for (n = out_top; n < out_bot; n++) {
	varr1[n] = y[ws->PhyL_N_i][n-out_top];      
      }
      nc_put_vara_double(ncid,varid,start15,count15,varr1);
    }

    if (ws->PhyS_Chl_i > -1){
      varid = ncw_var_id(ncid,"PhyS_Chl");
      for (n = out_top; n < out_bot; n++) {
	varr1[n] = y[ws->PhyS_Chl_i][n-out_top];      
      }
      nc_put_vara_double(ncid,varid,start15,count15,varr1);

      varid = ncw_var_id(ncid,"PhyS_N");
      for (n = out_top; n < out_bot; n++) {
	varr1[n] = y[ws->PhyS_N_i][n-out_top];      
      }
      nc_put_vara_double(ncid,varid,start15,count15,varr1);
    }

    if (ws->Tricho_Chl_i > -1){
      varid = ncw_var_id(ncid,"Tricho_Chl");
      for (n = out_top; n < out_bot; n++) {
	varr1[n] = y[ws->Tricho_Chl_i][n-out_top];      
      }
      nc_put_vara_double(ncid,varid,start15,count15,varr1);

      varid = ncw_var_id(ncid,"Tricho_N");
      for (n = out_top; n < out_bot; n++) {
	varr1[n] = y[ws->Tricho_N_i][n-out_top];      
      }
      nc_put_vara_double(ncid,varid,start15,count15,varr1);
    }

    d_free_1d(varr1);

    if (reclen < 2){
      varid = ncw_var_id(ncid,"i_index");
      nc_put_var_int(ncid,varid,&ij[0]);
    
      varid = ncw_var_id(ncid,"j_index");
      nc_put_var_int(ncid,varid,&ij[1]);

      varid = ncw_var_id(ncid,"latitude");
      nc_put_var_double(ncid,varid,&latitude);

      varid = ncw_var_id(ncid,"longitude");
      nc_put_var_double(ncid,varid,&longitude);
    }

    nc_close(ncid);

    // free memory used in netcdf write

    d_free_1d(varr);
    d_free_1d(E_d_z_tmp);
    d_free_1d(E_s_tmp);
    d_free_1d(at_tmp);
    d_free_1d(bt_tmp);
    d_free_1d(bb_tmp);
    d_free_1d(blt_tmp);
    d_free_1d(flt_tmp);

#ifdef INCLUDE_BIOLUMINESCENCE    
    d_free_1d(Secchi_colour);
    d_free_1d(E_u_z_bio_tmp);
#endif

#ifdef INCLUDE_HYDROLIGHT
    d_free_1d(E_d_hyd_tmp);
    d_free_1d(E_u_hyd_tmp);
    d_free_1d(L_d_tmp);
    d_free_1d(L_u_tmp);
    d_free_1d(L_wl_tmp);
    d_free_1d(L_sky_tmp);
    d_free_1d(L_ref_tmp);
#endif
    
  }
  
  // free memory from precalc in order of appearance above.

  d_free_2d(at);
  d_free_2d(bt);
  d_free_2d(bb);
  d_free_2d(blt);
  d_free_2d(flt);
  d_free_2d(Kd);
  d_free_2d(E_d_z);
  d_free_2d(E_d_z_fine);
  d_free_2d(E_d);
  d_free_2d(E_s);
  #ifdef INCLUDE_BIOLUMINESCENCE
     d_free_2d(E_u_z_bio);
  #endif
  d_free_1d(Rrs);
  d_free_1d(Rrs_fine);
  d_free_1d(n_stw);
  d_free_1d(costhetaw);
  d_free_1d(moonlight_SR);
  
  // free hydrolight arrays.

#ifdef INCLUDE_HYDROLIGHT
  if (col_write_num == 0){
    d_free_4d(L_d_hyd);
    d_free_4d(L_u_hyd);
    d_free_3d(L_sky);
    d_free_3d(L_ref);
    d_free_3d(L_wl);    
    d_free_2d(E_d_hyd);
    d_free_2d(E_u_hyd);
    
    d_free_1d(zen);
    d_free_1d(azi);
    d_free_1d(rad_dn);
    d_free_1d(rad_up);
    d_free_1d(rad_sky);
    d_free_1d(rad_ref);
    d_free_1d(rad_wl);
    d_free_1d(irrad_dn);
    d_free_1d(irrad_up);
    d_free_1d(Rrs_hyd);
    d_free_1d(Ed_ref);
    d_free_1d(Ed_wl);
    d_free_1d(Ed_sky);
    
    d_free_1d(a_z);
    d_free_1d(bb_z);
    d_free_1d(c_z);
    d_free_1d(s_z);
    d_free_1d(positive_zc);
    d_free_1d(zg);
  }
#endif
    
  }

void light_spectral_col_postcalc(eprocess* p, void* pp)
{
}

void passive_fluorescence_per_cell(ecology *e,workspace *ws, double kI, double Iquota, double *flpercell)
{
  /* calculates spectrally-resolved passive fluorescence per cell per waveband: W waveband-1 cell-1 */

  // calculate based on light flux in last timestep and physiological state.

  // kI is in W m-2 m2 cell-1 nm, which is done in values_common_precalc.

  // but this now needs to be distributed to another wavelength - ah but we lose energy so okay?

  double total = kI * (1.0 - 0.5 * Iquota) * Iquota;

  int w;

  for (w=0; w < ws->num_wave; w++){
    flpercell[w] = (total / ws->wave[w]) * ws->fls[w]/ws->fls_sum;
  }
}

void usgs_moon_reflectance(ecology *e,workspace *ws)
{

  /* Interpolate moon reflectance data from the USGS ROLO model onto spectral grid

   Kieffer, H. H and T. C. Stone (2005) The spectral irradiance of the moon. The Ast. J. 129: 2887-2901. */

  int num_wave = ws->num_wave;

  double moonwave[32] ={350.0,355.1,405.0,412.3,414.4,441.6,465.8,475.0,486.9,544.0,549.1,553.8,665.1,693.1,703.6,745.3,763.7,774.8,865.3,872.6,882.0,928.4,939.3,942.1,1059.5,1243.2,1538.7,1633.6,1981.5,2126.3,2250.9,2383.6};

  double moondata[32][10] = {{-2.67511,-1.78539, 0.50612,-0.25578,0.03744, 0.00981,-0.00322,0.34185, 0.01441,-0.01602},
			    {-2.71924,-1.74298, 0.44523,-0.23315,0.03492, 0.01142,-0.00383,0.33875, 0.01612,-0.00996},
			    {-2.35754,-1.72134, 0.40337,-0.21105,0.03505, 0.01043,-0.00341,0.35235,-0.03818,-0.00006},
			    {-2.34185,-1.74337, 0.42156,-0.21512,0.03141, 0.01364,-0.00472,0.36591,-0.05902, 0.00080},
			    {-2.43367,-1.72184, 0.43600,-0.22675,0.03474, 0.01188,-0.00422,0.35558,-0.03247,-0.00503},
			    {-2.31964,-1.72114, 0.37286,-0.19304,0.03736, 0.01545,-0.00559,0.37935,-0.09562, 0.00970},
			    {-2.35085,-1.66538, 0.41802,-0.22541,0.04274, 0.01127,-0.00439,0.33450,-0.02546,-0.00484},
			    {-2.28999,-1.63180, 0.36193,-0.20381,0.04007, 0.01216,-0.00437,0.33024,-0.03131, 0.00222},
			    {-2.23351,-1.68573, 0.37632,-0.19877,0.03881, 0.01566,-0.00555,0.36590,-0.08945, 0.00678},
			    {-2.13864,-1.60613, 0.27886,-0.16426,0.03833, 0.01189,-0.00390,0.37190,-0.10629, 0.01428},
			    {-2.10782,-1.66736, 0.41697,-0.22026,0.03451, 0.01452,-0.00517,0.36814,-0.09815, 0.00000},
			    {-2.12504,-1.65970, 0.38409,-0.20655,0.04052, 0.01009,-0.00388,0.37206,-0.10745, 0.00347},
			    {-1.88914,-1.58096, 0.30477,-0.17908,0.04415, 0.00983,-0.00389,0.37141,-0.13514, 0.01248},
			    {-1.89410,-1.58509, 0.28080,-0.16427,0.04429, 0.00914,-0.00351,0.39109,-0.17048, 0.01754},
			    {-1.92103,-1.60151, 0.36924,-0.20567,0.04494, 0.00987,-0.00386,0.37155,-0.13989, 0.00412},
			    {-1.86896,-1.57522, 0.33712,-0.19415,0.03967, 0.01318,-0.00464,0.36888,-0.14828, 0.00958},
			    {-1.85258,-1.47181, 0.14377,-0.11589,0.04435, 0.02000,-0.00738,0.39126,-0.16957, 0.03053},
			    {-1.80271,-1.59357, 0.36351,-0.20326,0.04710, 0.01196,-0.00476,0.36908,-0.16182, 0.00830},
			    {-1.74561,-1.58482, 0.35009,-0.19569,0.04142, 0.01612,-0.00550,0.39200,-0.18837, 0.00978},
			    {-1.76779,-1.60345, 0.37974,-0.20625,0.04645, 0.01170,-0.00424,0.39354,-0.19360, 0.00568},
			    {-1.73011,-1.61156, 0.36115,-0.19576,0.04847, 0.01065,-0.00404,0.40714,-0.21499, 0.01146},
			    {-1.75981,-1.45395, 0.13780,-0.11254,0.05000, 0.01476,-0.00513,0.41900,-0.19963, 0.02940},
			    {-1.76245,-1.49892, 0.07956,-0.07546,0.05461, 0.01355,-0.00464,0.47936,-0.29463, 0.04706},
			    {-1.66473,-1.61875, 0.14630,-0.09216,0.04533, 0.03010,-0.01166,0.57275,-0.38204, 0.04902},
			    {-1.59323,-1.71358, 0.50599,-0.25178,0.04906, 0.03178,-0.01138,0.48160,-0.29486, 0.00116},
			    {-1.53594,-1.55214, 0.31479,-0.18178,0.03965, 0.03009,-0.01123,0.49040,-0.30970, 0.01237},
			    {-1.33802,-1.46208, 0.15784,-0.11712,0.04674, 0.01471,-0.00656,0.53831,-0.38432, 0.03473},
			    {-1.34567,-1.46057, 0.23813,-0.15494,0.03883, 0.02280,-0.00877,0.54393,-0.37182, 0.01845},
			    {-1.26203,-1.25138,-0.06569,-0.04005,0.04157, 0.02036,-0.00772,0.49099,-0.36092, 0.04707},
			    {-1.18946,-2.55069, 2.10026,-0.87285,0.03819,-0.00685,-0.00200,0.29239,-0.34784,-0.13444},
			    {-1.04232,-1.46809, 0.43817,-0.24632,0.04893, 0.00617,-0.00259,0.38154,-0.28937,-0.01110},
			    {-1.08403,-1.31032, 0.20323,-0.15863,0.05955,-0.00940, 0.00083,0.36134,-0.28408, 0.01010}};

  /* Now interpolate onto model's wavelengths */

  double dummy[32];int count;int w;

  ws->moon_a0 = d_alloc_1d(num_wave);
  
  for (count=0; count<32; count++){
    dummy[count] = moondata[count][0];
  }
  interp1d(moonwave, dummy, 32, ws->wave, ws->moon_a0, num_wave);
  
  ws->moon_a1 = d_alloc_1d(num_wave);
  for (count=0; count<32; count++){
    dummy[count] = moondata[count][1];
  }
  interp1d(moonwave, dummy, 32, ws->wave, ws->moon_a1, num_wave);

  ws->moon_a2 = d_alloc_1d(num_wave);
  for (count=0; count<32; count++){
    dummy[count] = moondata[count][2];
  }
  interp1d(moonwave, dummy, 32, ws->wave, ws->moon_a2, num_wave);

  ws->moon_a3 = d_alloc_1d(num_wave);
  for (count=0; count<32; count++){
    dummy[count] = moondata[count][3];
  }
  interp1d(moonwave, dummy, 32, ws->wave, ws->moon_a3, num_wave);

  ws->moon_b1 = d_alloc_1d(num_wave);
  for (count=0; count<32; count++){
    dummy[count] = moondata[count][4];
  }
  interp1d(moonwave, dummy, 32, ws->wave, ws->moon_b1, num_wave);

  ws->moon_b2 = d_alloc_1d(num_wave);
  for (count=0; count<32; count++){
    dummy[count] = moondata[count][5];
  }
  interp1d(moonwave, dummy, 32, ws->wave, ws->moon_b2, num_wave);

  ws->moon_b3 = d_alloc_1d(num_wave);
  for (count=0; count<32; count++){
    dummy[count] = moondata[count][6];
  }
  interp1d(moonwave, dummy, 32, ws->wave, ws->moon_b3, num_wave);

  ws->moon_d1 = d_alloc_1d(num_wave);
  for (count=0; count<32; count++){
    dummy[count] = moondata[count][7];
  }
  interp1d(moonwave, dummy, 32, ws->wave, ws->moon_d1, num_wave);

  ws->moon_d2 = d_alloc_1d(num_wave);
  for (count=0; count<32; count++){
    dummy[count] = moondata[count][8];
  }
  interp1d(moonwave, dummy, 32, ws->wave, ws->moon_d2, num_wave);

  ws->moon_d3 = d_alloc_1d(num_wave);
  for (count=0; count<32; count++){
    dummy[count] = moondata[count][9];
  }
  interp1d(moonwave, dummy, 32, ws->wave, ws->moon_d3, num_wave);

  eco_write_setup(e,"\nMoon disk reflectance coefficients calculated in light_spectral_col.c \n \n");
    
    eco_write_setup(e,"wave \t a0 \t\t a1 \t\t a2 \t\t a3 \t\t b1 \t\t b2 \t\t b3 \t\t d1 \t\t d2 \t\t d3 \n");
    for (w=0; w<num_wave; w++){
      eco_write_setup(e,"%4.2f \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \n",ws->wave[w],ws->moon_a0[w],ws->moon_a1[w],ws->moon_a2[w],ws->moon_a3[w],ws->moon_b1[w],ws->moon_b2[w],ws->moon_b3[w],ws->moon_d1[w],ws->moon_d2[w],ws->moon_d3[w]);
    }
    eco_write_setup(e,"\n");

}

void moonlight(workspace *ws, column* col, double time_in_days, double* Ak, double *sun_angle1, double *moon_phase1, double *lunar_dec1, double *lunar_angle1, double *lunar_zenith1, double *lunar_azimuth1)

// Need ws for num_wave, landa and moon.
// Need col for to pass column number and link to model column mapping.
{
  int w;

  // parameters calculated in moonvars

  // Is this the problem should moonlat and moon lon really be zero?
  
  double moonlat_obs = 0.0;
  double moonlon_obs = 0.0;

  double earth_sun_dist; 
  double moon_earth_dist;

  // Function output parameters

  double lunar_angle;
  double sun_angle;
  double moon_phase;
  double lunar_dec;
  double lunar_zenith;
  double lunar_azimuth;

  ginterface_moonvars(col->model,col->b, &moonlon_obs, &moonlat_obs, &earth_sun_dist, &moon_earth_dist, &lunar_angle, &sun_angle, &moon_phase, &lunar_dec);
      
  // Could have an if statement for moon below horizon

  // Calculate moon phase

  moon_phase = fmod(time_in_days - 3658.18,29.530588853)/29.530588853*M_PI;

  // Mark to add extra code:
  
  double gmo = fabs(moon_phase-M_PI/2.0); // in radians
  
  // Still not given a good number. 
  
  double moonlon_sun = 7.0/2.0/M_PI; // pi - (solar hour angle - monlon_obs) excluding librations.
  // need to be careful with going past pi.

  double moonlon_obs_deg = moonlon_obs*180.0/M_PI;
  double moonlat_obs_deg = moonlat_obs*180.0/M_PI;
  double moonlon_sun_deg = 0.0;
  
  double c1 = 0.00034115;
  double c2 = -0.0013425;
  double c3 = 0.00095906;
  double c4 = 0.00066229;
  double p1 = 4.06054*M_PI/180.0; // in radians
  double p2 = 12.8802*M_PI/180.0;
  double p3 = -30.5858*M_PI/180.0;
  double p4 = 16.7498*M_PI/180.0;
  
  double Ak_nonwave = c1*moonlat_obs_deg + c2*moonlon_obs_deg +c3*moonlon_sun_deg*moonlat_obs_deg + c4*moonlon_sun_deg*moonlon_obs_deg;
  
  for (w=0; w<ws->num_wave; w++){
    Ak[w] = ws->moon_a0[w] + ws->moon_a1[w] * gmo + ws->moon_a2[w] * pow(gmo,2.0) + ws->moon_a3[w] * pow(gmo,3.0);
    Ak[w] += ws->moon_b1[w] * moonlon_sun_deg + ws->moon_b2[w] * pow(moonlon_sun_deg,3.0) + ws->moon_b3[w] * pow(moonlon_sun_deg,5.0);
    Ak[w] += Ak_nonwave;
    Ak[w] += ws->moon_d1[w] * exp(-gmo/p1) + ws->moon_d2[w] * exp(-gmo/p2) + ws->moon_d3[w] * cos((gmo-p3)/p4);
    Ak[w] = exp(Ak[w]);
    Ak[w] *= 6.4177e-5 * 1366.0 * ws->landa[w] * pow(earth_sun_dist/1496.0e8 ,2.0) * pow(moon_earth_dist/384400000.0,2.0)/M_PI;
  }

  lunar_zenith = sin(lunar_dec)*sin(moonlat_obs)+cos(lunar_dec)*cos(moonlat_obs)*cos(lunar_angle);
        
  lunar_zenith = acos(lunar_zenith);
      
  lunar_zenith = (fabs(lunar_zenith)<1.0e-15)? 1.0e-15:lunar_zenith;
  lunar_zenith = ((lunar_zenith-M_PI/2.0)-fabs(lunar_zenith-M_PI/2.0))/2.0+M_PI/2.0;

  // Following a north-clockwise convention

  if (lunar_zenith < M_PI/2.0){
    lunar_azimuth = - sin(lunar_angle) * cos(lunar_dec) / sin(lunar_zenith);
    lunar_azimuth = asin(lunar_azimuth);
  }else{
    lunar_azimuth = 0.0;
  }

  *lunar_angle1 = lunar_angle;
  *sun_angle1   = sun_angle;
  *moon_phase1 = moon_phase;
  *lunar_dec1 = lunar_dec;
  *lunar_zenith1 = lunar_zenith;
  *lunar_azimuth1 = lunar_azimuth;
}

void ems_swr_airtrans_and_albedo(double tcld, double lunar_dec, double lunar_angle, double lunar_zenith, double lat, double *airtransmission, double *swr_albedo)

{
/* This routine calculates the fractional SWR loss through the atmosphere and surface albedo following that calculated in heatflux.c */

  // lunar_angle = hour angle in radians.
  
  double albedo,air,airmass;
  double es = 26.0;                       /* Default vapour pressure                     */
  int  oktas = (int)(tcld * 8);
  double h; //  solar elevation.
  double d2;

  double lat1 = lat*M_PI/180.0;  // convert to radians.
  
  h = sin(lat1) * sin(lunar_dec) + cos(lat1) * cos(lunar_dec) * cos(lunar_angle);
  if (h < 0.0)
    h = 0.0;

 double ap[9][10] =
    { {0.25, 0.25, 0.13, 0.08, 0.06, 0.05, 0.045, 0.04, 0.04, 0.04},
      {0.25, 0.25, 0.13, 0.08, 0.06, 0.05, 0.045, 0.04, 0.04, 0.04},
      {0.25, 0.25, 0.13, 0.08, 0.06, 0.05, 0.045, 0.04, 0.04, 0.04},
      {0.25, 0.25, 0.13, 0.08, 0.06, 0.05, 0.045, 0.04, 0.04, 0.04},
      {0.25, 0.25, 0.13, 0.08, 0.06, 0.05, 0.045, 0.04, 0.04, 0.04},
      {0.22, 0.22, 0.12, 0.08, 0.065, 0.055, 0.05, 0.045, 0.045, 0.045},
      {0.19, 0.19, 0.11, 0.08, 0.065, 0.055, 0.05, 0.05, 0.05, 0.05},
      {0.14, 0.14, 0.09, 0.075, 0.068, 0.063, 0.06, 0.06, 0.06, 0.06},
      {0.09, 0.09, 0.075, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07}
    };
  
  airmass = h * h / ((h + 2.7) * 1.0e-3 * es + 1.085 * h + 0.1); // Zilman 1972,

  // Now reduce again for cloudy skies.

  d2 = sin(lat1) * sin(lunar_dec) + cos(lat1) * cos(lunar_dec);
  if (oktas > 2)
       airmass *= 1.0 - 0.62 * tcld + 0.0019 * asin(d2);

  // Albedo of sea-surface changes with cloud cover and zenith angle.
  
  h = asin(h) * 180.0 / M_PI;
  int iw;
  double d1; /* Counters */
  iw = (int)(h / 10.0);
  d1 = (double)iw *10.0;

  if (iw == 0)
    albedo = ap[oktas][0]; // albedo = 25% !
  else
    albedo = ((ap[oktas][iw + 1] - ap[oktas][iw]) / 10.0) * (h - d1) + ap[oktas][iw];

  *airtransmission = airmass;
  *swr_albedo = albedo;  
}
void meb_swr_airtrans_and_albedo(double tcld, double lunar_dec, double lunar_angle, double lunar_zenith, double lat, double *airtransmission, double *swr_albedo)

{
/* This routine calculates the SWR loss through the atmosphere and surface albedo following that calculated in heatflux.c */
/* NOT TESTED */
  // lunar_angle = hour angle in radians.
  
  double albedo1,albedo2,air;
  double es = 26.0;                       /* Default vapour pressure                     */
  int  oktas = (int)(tcld * 8);
  double airmass1,airmass2;
  double h; //  elevation.
  
  h = sin(lat) * sin(lunar_dec) + cos(lat) * cos(lunar_dec) * cos(lunar_angle);
  if (h < 0.0)
    h = 0.0;

  // printf("ems_swr: lat %e lunar_dec %e lunar_angle %e hour angle %e \n",lat,lunar_dec,lunar_angle,h);
  
  double ap[9][10] =
    { {0.25, 0.25, 0.13, 0.08, 0.06, 0.05, 0.045, 0.04, 0.04, 0.04},
      {0.25, 0.25, 0.13, 0.08, 0.06, 0.05, 0.045, 0.04, 0.04, 0.04},
      {0.25, 0.25, 0.13, 0.08, 0.06, 0.05, 0.045, 0.04, 0.04, 0.04},
      {0.25, 0.25, 0.13, 0.08, 0.06, 0.05, 0.045, 0.04, 0.04, 0.04},
      {0.25, 0.25, 0.13, 0.08, 0.06, 0.05, 0.045, 0.04, 0.04, 0.04},
      {0.22, 0.22, 0.12, 0.08, 0.065, 0.055, 0.05, 0.045, 0.045, 0.045},
      {0.19, 0.19, 0.11, 0.08, 0.065, 0.055, 0.05, 0.05, 0.05, 0.05},
      {0.14, 0.14, 0.09, 0.075, 0.068, 0.063, 0.06, 0.06, 0.06, 0.06},
      {0.09, 0.09, 0.075, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07}
    };
  
  airmass1 = h * h / ((h + 2.7) * 1.0e-3 * es + 1.085 * h + 0.1); // Zilman 1972,

  double Z = lunar_zenith * 180.0/M_PI;
  airmass2 = 1.0/(cos(lunar_zenith) + 0.15*pow((93.885 - Z),-1.253)); // Gregg and Carder (1990).

  double thetaW = asin(sin(lunar_zenith)/1.33);

  double AAA = pow(sin(lunar_zenith - thetaW),2.0)/pow(sin(lunar_zenith + thetaW),2.0);
  double BBB = pow(tan(lunar_zenith - thetaW),2.0)/pow(tan(lunar_zenith + thetaW),2.0);
  
  albedo2 = 0.5*(AAA+BBB);

  // Albedo of sea-surface changes with cloud cover and zenith angle.
  
  h = asin(h) * 180.0 / M_PI;
  int iw;
  double d1; /* Counters */
  iw = (int)(h / 10.0);
  d1 = (double)iw *10.0;

  if (iw == 0)
    albedo1 = ap[oktas][0]; // albedo = 25% !
  else
    albedo1 = ((ap[oktas][iw + 1] - ap[oktas][iw]) / 10.0) * (h - d1) + ap[oktas][iw];

  //  printf("ems: oktas %d, lunar_angle %e, zenith %e Zilman %e, albedo %e, h %e \n",oktas,lunar_angle,lunar_zenith,airmass1,albedo1,h);

  // printf("meb: oktas %d, lunar_angle %e, zenith %e Gregg  %e, albedo %e \n",oktas,lunar_angle,lunar_zenith,airmass2,albedo2);

  *airtransmission = airmass2;
  *swr_albedo = albedo2;
}

void ModelRsestoBandRs(ecology *e,workspace *ws, double* response, double* modelR, double *obsR)
{
  /* This function takes model reflectances and spectral response curves and returns band reflectance */ 
  
  int w;
  double sum = 0.0;
  double denom = 0.0;

  for (w=0; w<ws->num_wave; w++){
    sum += modelR[w]*response[w]*ws->bandwidth[w];
    denom += response[w]*ws->bandwidth[w];
  }
  *obsR = sum/denom;
}

void TSSM(double band1, double *tssm_out)
{
  /* TSS algorithm - use local coastal relationship from Petus et al., 2014 */

  *tssm_out = (12450.0*band1*band1 + 666.0*band1 + 0.48)/1000.0;

}

void KD490M(double band4, double band6, double *kd490m_out)
{
  double a_kd2m[5] = {-0.8813, -2.0584, 2.5878, -3.4885, -1.5061};

  double X = log10(band4/band6);

  *kd490m_out  = pow(10.0,a_kd2m[0] + X * (a_kd2m[1] + X * (a_kd2m[2] + X * (a_kd2m[3] + X * a_kd2m[4])))) + 0.0166;
}

void OC4Me(double band3, double band4, double band5, double band6, double *oc4me_out)
{
  /* Calculate OC4Me from Sentinel bands calculated from hyperspectral reflectance and response curves */

  double ratio;
  double a_oc4[5] = {0.3255, -2.7677, 2.4409, -1.1288, -0.4990};
  
  if (band3 > band4 && band3 > band5){
    ratio = log10(band3/band6);
  }else{
    if (band4 > band5)
      ratio = log10(band4/band6);
    else
      ratio = log10(band5/band6);
  }
  *oc4me_out = pow(10.0,a_oc4[0] + ratio *(a_oc4[1] + ratio *(a_oc4[2] + ratio *(a_oc4[3] + ratio * a_oc4[4]))));
}

void OC3M(double band2, double band4, double band6, double *oc3m_out)
{
  /* Calculate OC3M from MODIS bands calculated from hyperspectral reflectance and response curves */

  double ratio;
  double a_oc3[5] = {0.283, -2.753, 1.457, 0.659, -1.403};
    
  if (band2 > band4) 
    ratio = log10(band2/band6);
  else
    ratio = log10(band4/band6);
  
  *oc3m_out = pow(10.0,a_oc3[0] + ratio *(a_oc3[1] + ratio *(a_oc3[2] + ratio *(a_oc3[3] + ratio * a_oc3[4]))));
}

void OC3V(double band2, double band3, double band4, double *oc3v_out)
{
  /* Calculate OC3V from VIIRS bands calculated from hyperspectral reflectance and response curves */

  double ratio;
  double a_oc3[5] = {0.2228,-2.4683,1.5867,-0.4275,-0.7768};
    
  if (band2 > band3) 
    ratio = log10(band2/band4);
  else
    ratio = log10(band3/band4);
  
  *oc3v_out = pow(10.0,a_oc3[0] + ratio *(a_oc3[1] + ratio *(a_oc3[2] + ratio *(a_oc3[3] + ratio * a_oc3[4]))));
}

void Hue3B(double band3, double band4, double band6, double band8, double band11, double *hue3b_out)
{

  double hue;
  
  double x_tris_MSI = 11.756 * band3 +  6.423 * band4 + 53.696 * band6 + 32.028 * band8 + 0.529 * band11;
  double y_tris_MSI =  1.744 * band3 + 22.289 * band4 + 65.702 * band6 + 16.808 * band8 + 0.192 * band11;
  double z_tris_MSI = 62.696 * band3 + 31.101 * band4 +  1.778 * band6 +  0.015 * band8 + 0.000 * band11;
  
  double cx_MSI = x_tris_MSI/(x_tris_MSI+y_tris_MSI+z_tris_MSI+1.0e-6);
  double cy_MSI = y_tris_MSI/(x_tris_MSI+y_tris_MSI+z_tris_MSI+1.0e-6);
  
  double xw_MSI = cx_MSI - 0.333333333;
  double yw_MSI = cy_MSI - 0.333333333;
  
  double atanterm = atan(yw_MSI/(xw_MSI+1.0e-6));
    
  hue = atanterm*180.0/M_PI;

  if (xw_MSI<0.0){
    if (yw_MSI<0.0){
      hue = 180.0 + atanterm*180.0/M_PI;
    }else{
      hue = 180.0 - atanterm*180.0/M_PI;
    }
  }
  *hue3b_out = hue; 
}


void nFLH_modis(double band9, double band10, double band11,double *nFLH_modis_out)
{
  double nFLH;
  
  nFLH = (148.097 / 3.145926535) * (band10 - (70.0/81.0) * band9 - (11.0/81.0) * band11);
 
  *nFLH_modis_out = nFLH;
}

void index_of_refraction(ecology *e,workspace *ws,double **y, double* n_stw)
{
  // Quan and Fry (1995) Appl. Opt. 34: 3477-3480.

  double S = y[ws->salt_i][0];
  double T = y[ws->temp_i][0];

  double n0 = 1.31405;
  double n1 = 1.779e-4;
  double n2 = -1.05e-6;
  double n3 = 1.6e-8;
  double n4 = -2.02e-6;
  double n5 = 15.868;
  double n6 = 0.01155;
  double n7 = -0.00423;
  double n8 = -4382.0;
  double n9 = 1.1455e6;

  double W;
  int w;

  for (w=0; w<ws->num_wave; w++){
    W = ws->wave[w];
    n_stw[w] = n0 + (n1 + n2*T + n3*T*T)*S + n4*T*T + (n5 + n6*S + n7*T)/W + n8/(W*W) + n9/(W*W*W);
  }
}

void optdatascalarread(ecology *e, workspace *ws, char *source_file, int ncid, const char *varname, double *out)
{
  // read an optical constant from database, place in workspace and write to setup file.

  int varid;
  double value;

  char units[30];
  char description[100];
  char symbol[30];

  char *output_file = "optical_setup.nc";
  
  varid = ncw_var_id(ncid,varname);

  if (varid > -1){
    ncw_get_var_double(source_file,ncid, varid,&value);
    ncw_get_att_text(source_file,ncid,varid,"Units",units);
    ncw_get_att_text(source_file,ncid,varid,"Symbol",symbol);
    ncw_get_att_text(source_file,ncid,varid,"Description",description);

     // now open output netcdf and put in variable.

    int ncid1,varid1,dims[2];
      
    ncw_open(output_file,NC_WRITE, &ncid1);
    ncredef(ncid1);
    nc_def_var(ncid1,varname,NC_DOUBLE,0,dims,&varid1);
    nc_put_att_text(ncid1, varid1, "description", strlen(description),description);
    nc_put_att_text(ncid1, varid1, "sourcename", strlen(varname),varname);
    nc_put_att_text(ncid1, varid1, "symbol", strlen(symbol),symbol);
    nc_put_att_text(ncid1, varid1, "units", strlen(units),units);
    char *hhh = "N";
    nc_put_att_text(ncid1, varid1, "spectrally_resolved",strlen(hhh),hhh);
    
    ncw_enddef(output_file,ncid1);
    ncw_put_var_double(output_file,ncid1,varid1,&value);
    nc_close(ncid1);

    eco_write_setup(e,"optdatascalarread: reading %s for variable %s = %e %s.\n",source_file,varname,value,units);
    
    *out = value;
  }
}
 
void colour_in_seawater(ecology *e, workspace *ws, const char *varname)
{

  /* Add tracer attribute of rgb colour and also units in tracer list */

  char *output_file = "optical_setup.nc";
  int ncid,varid1,varid2,varid3,w;
  char output_name[MAXSTRLEN];
  char tracerunits[30] = "";

  ginterface_get_tracerunits(e->model,varname, tracerunits);

  eco_write_setup(e,"Calculating colour in sea water for variable %s \n",varname);
  
  ncw_open("optical_setup.nc",NC_WRITE, &ncid);

  sprintf(output_name, "ap_%s",varname);
  varid1 = ncw_var_id(ncid,output_name);

  sprintf(output_name, "bp_%s",varname);
  varid2 = ncw_var_id(ncid,output_name);

  sprintf(output_name, "B_%s",varname);
  varid3 = ncw_var_id(ncid,output_name);

  eco_write_setup(e,"IOPs varid1 %d, varid2 %d varid3 %d \n",varid1,varid2,varid3);

  if (varid1 > -1){

    int indw[3];
    
    static double conc;double rgbb[3];double rgb[3];double at[3];double bb[3];
    double atm[ws->num_wave];
    double btm[ws->num_wave];
    double Bm[ws->num_wave];

    if ((ws->w470 > -1) && (ws->w550 > -1) && (ws->w670 > -1)){

      indw[2] = ws->w470; indw[1] = ws->w550; indw[0] = ws->w670;

      // Calculate rgb based on MODIS true colour algorithm
 
      double unscaledR[6] = {0.0, 30.0, 60.0, 120.0, 190.0, 255.0};
      double scaledR[6] = {1.0, 110.0, 160.0, 210.0, 240.0, 255.0};
      size_t ncstart[1] = {0}; size_t nccount[1] = {ws->num_wave};

      conc = 1.0;

      ncw_get_vara_double("optical_setup.nc",ncid, varid1, ncstart,nccount,atm);
      ncw_get_vara_double("optical_setup.nc",ncid, varid2, ncstart,nccount,btm);
      ncw_get_vara_double("optical_setup.nc",ncid, varid3, ncstart,nccount,Bm);        

      for (w = 0; w<3; w++) {
	at[w] = ws->aw[indw[w]] + atm[indw[w]] * conc;
	bb[w] = 0.5 * ws->bw[indw[w]] + (btm[indw[w]] * conc) * Bm[indw[w]];
	rgbb[w] = bb[w]/(at[w] + bb[w])*256.0;
      }
      interp1d(unscaledR,scaledR,6,rgbb,rgb,3);
      for (w = 0; w<3; w++) {
	eco_write_setup(e,"mass-specific IOPs: atm %e, btm %e, Bm %e \n",atm[indw[w]],btm[indw[w]],Bm[indw[w]]);
	eco_write_setup(e,"RGB values: at %e, bb %e, indw %d, w %d, rgbb %e rgb %e \n",at[w],bb[w],indw[w],w,rgbb[w],rgb[w]);
      }
      ncredef(ncid);
      nc_put_att_double(ncid, varid1, "rgb_colour_space", NC_DOUBLE, 3,rgb);
      nc_put_att_double(ncid, varid2, "rgb_colour_space", NC_DOUBLE, 3,rgb);
      nc_put_att_double(ncid, varid3, "rgb_colour_space", NC_DOUBLE, 3,rgb);
      nc_put_att_double(ncid, varid1, "mass_concentration_used_for_rgb_colour_space_in_clear_water", NC_DOUBLE, 1,&conc);
      nc_put_att_double(ncid, varid2, "mass_concentration_used_for_rgb_colour_space_in_clear_water", NC_DOUBLE, 1,&conc);
      nc_put_att_double(ncid, varid3, "mass_concentration_used_for_rgb_colour_space_in_clear_water", NC_DOUBLE, 1,&conc);
      nc_put_att_text(ncid, varid1, "tracer_units", strlen(tracerunits),tracerunits);
      nc_put_att_text(ncid, varid2, "tracer_units", strlen(tracerunits),tracerunits);
      nc_put_att_text(ncid, varid3, "tracer_units", strlen(tracerunits),tracerunits);
      char *hhh = "Y";
      nc_put_att_text(ncid, varid1, "spectrally_resolved", strlen(hhh),hhh);
      ncw_enddef("optical_setup.nc", ncid);
    }
  }
  nc_close(ncid);
}

void colour_of_bottom_type(ecology *e, workspace *ws, const char *varname)
{
  char *output_file = "optical_setup.nc";
  int ncid,varid1,varid2,varid3,w;
  char output_name[MAXSTRLEN];

  char tracerunits[30] = "";

  eco_write_setup(e,"Calculating colour of bottom type for variable %s \n",varname);

  ncw_open("optical_setup.nc",NC_WRITE, &ncid);

  varid1 = -1;

  sprintf(output_name, "dR_%s",varname);
  varid1 = ncw_var_id(ncid,output_name);

  eco_write_setup(e,"In colour_of_bottom_type - 3D in top sediment layer:  varname %s, outputname %s varid %d\n",varname,output_name,varid1);

  if (varid1 < 0){   // try again as it may be a epibenthic variable.

    sprintf(output_name, "R_%s",varname);
    varid1 = ncw_var_id(ncid,output_name);

    sprintf(output_name, "T_%s",varname);
    varid2 = ncw_var_id(ncid,output_name);

    sprintf(output_name, "A_%s",varname);
    varid3 = ncw_var_id(ncid,output_name);
    
    eco_write_setup(e,"In colour_of_bottom_type - Epibenthic:  varname %s, outputname %s varid %d\n",varname,output_name,varid1);
    
  }else{
    ginterface_get_tracerunits(e->model,varname, tracerunits); // must be a tracer.
  }

  if (varid1 > -1){

    int indw[3];
    
    double rgbb[3];double rgb[3];
    double rho[ws->num_wave];

    if ((ws->w470 > -1) && (ws->w550 > -1) && (ws->w670 > -1)){

      indw[2] = ws->w470; indw[1] = ws->w550; indw[0] = ws->w670;

      // Calculate rgb based on MODIS true colour algorithm
 
      double unscaledR[6] = {0.0, 30.0, 60.0, 120.0, 190.0, 255.0};
      double scaledR[6] = {1.0, 110.0, 160.0, 210.0, 240.0, 255.0};

      size_t ncstart[1] = {0}; size_t nccount[1] = {ws->num_wave};

      ncw_get_vara_double("optical_setup.nc",ncid, varid1, ncstart,nccount,rho);      

      for (w = 0; w<3; w++) {
	rgbb[w] = rho[indw[w]]*256.0;
      }
      interp1d(unscaledR,scaledR,6,rgbb,rgb,3);
      for (w = 0; w<3; w++) {
	eco_write_setup(e,"RGB values: indw %d, w %d, rgbb %e rgb %e \n",indw[w],w,rgbb[w],rgb[w]);
      }
      ncredef(ncid);
      nc_put_att_double(ncid, varid1, "rgb_colour_space", NC_DOUBLE, 3,rgb);
      nc_put_att_double(ncid, varid2, "rgb_colour_space", NC_DOUBLE, 3,rgb);
      nc_put_att_double(ncid, varid3, "rgb_colour_space", NC_DOUBLE, 3,rgb);
      // nc_put_att_double(ncid, varid1, "rgb_colour_space_in_1m_of_seawater", NC_DOUBLE, 3,rgb);
      // nc_put_att_double(ncid, varid1, "rgb_colour_space_in_10m_of_seawater", NC_DOUBLE, 3,rgb);
      nc_put_att_text(ncid, varid1, "tracer_units", strlen(tracerunits),tracerunits);
      nc_put_att_text(ncid, varid2, "tracer_units", strlen(tracerunits),tracerunits);
      nc_put_att_text(ncid, varid3, "tracer_units", strlen(tracerunits),tracerunits);
      ncw_enddef("optical_setup.nc", ncid);
    }
  }
  nc_close(ncid);
}

void optdataspectralread(ecology *e, workspace *ws, char *source_file, int ncid, const char *varname, const char *wavename, const char *codename, double* out)
{
  // read a spectrally-resolved variable from database, interpolate onto optical grid, place in workspace and write to setup file.

  int dimid,w,varid;
  int dimids[2];
  size_t len;
  size_t ncstart[2] = {0,0}; size_t nccount[2] = {2,1};

  char units[30];
  char description[100];
  char symbol[30];

  char *output_file = "optical_setup.nc";
  
  varid = ncw_var_id(ncid,varname);

  if (varid > -1){

    eco_write_setup(e,"optdataspectralread: reading %s for variable %s \n",source_file,varname);

    nc_inq_vardimid(ncid, varid, dimids);
    nc_inq_dimlen(ncid, dimids[1],&len);

    double nc_var[2*len];
    ncw_get_att_text(source_file,ncid,varid,"Units",units);
    ncw_get_att_text(source_file,ncid,varid,"Symbol",symbol);
    ncw_get_att_text(source_file,ncid,varid,"Description",description);
    nccount[1] = len;
    ncw_get_vara_double(source_file,ncid, varid, ncstart,nccount,nc_var);
    
    /* now interpolate onto model wavelengths */
    
    double t1[len];
    double t2[len];
    
    for (w=0; w<len; w++){
      t1[w] = nc_var[w];
      t2[w] = nc_var[w+len];
    }
    
    interp1d(t1, t2, len, ws->wave, out, ws->num_wave);

    eco_write_setup(e,"Write %s: %s [%s] to %s \n",description,symbol,units,output_file);
    eco_write_setup(e,"wave [nm] %s \n",varname);

    for (w=0; w<ws->num_wave; w++){
      eco_write_setup(e,"%4.2f \t %e \n",ws->wave[w],out[w]);
    }

    // now open output netcdf and put in variable.

      int ncid1,varid1,dims1,dimid1;
      int dims[2];
      
      ncw_open("optical_setup.nc",NC_WRITE, &ncid1);
      dimid1 = ncw_dim_id(ncid1,"wave_centre");
      dims[0] = dimid1;
      ncredef(ncid1);
      ncw_def_var("optical_setup.nc",ncid1,codename,NC_DOUBLE,1,dims,&varid1);
      nc_put_att_text(ncid1, varid1, "description", strlen(description),description);
      nc_put_att_text(ncid1, varid1, "codename", strlen(codename),codename);
      nc_put_att_text(ncid1, varid1, "sourcefile", strlen(source_file),source_file); // would be good to have creation_date of file
      nc_put_att_text(ncid1, varid1, "sourcename", strlen(varname),varname);
      nc_put_att_text(ncid1, varid1, "symbol", strlen(symbol),symbol);
      nc_put_att_text(ncid1, varid1, "units", strlen(units),units);

      // Now do test of error of interpolation onto model optical grid using area under the curve.

      double area_model = 0.0;
      double area_data = 0.0;

      int model_top = -1;
      int model_bot = -1;

      // find top and bottom model wavelengths in the data range for interpolation error calculation.

      for (w=0; w<ws->num_wave; w++){
        if (ws->wave[w] >= t1[0] && model_bot < 0 && ws->wave[w] >= 400.0){
          model_bot = w;
        }
        if ((ws->wave[w] >= t1[len-1] || ws->wave[w] >= 700.0 ) && model_top < 0 ){
          model_top = w-1;
        }
      }

      if (model_top == -1){
	model_top = ws->num_wave-1;
      }

      // Now find data points closest to the choosen model range.

      int data_top = -1;
      int data_bot = -1;

      if (len>2){

	for (w=0; w<len; w++){
	  if (t1[w] >= ws->wave[model_bot] && data_bot < 0){
	    data_bot = w;
	  }
	  if (t1[w] >= ws->wave[model_top] && data_top < 0){
	    data_top = w-1;
	  }
	}
	if (data_top == -1){
	data_top = len-1;
	}
      }else{
	data_top = 1;
	data_bot = 0;
      }

      for (w=model_bot; w <= model_top; w++){
        area_model += out[w]*ws->bandwidth[w];
      }

      static double areanorm;
      areanorm = area_model;
      nc_put_att_double(ncid1, varid1, "area_under_curve_of_model_in_data_range", NC_DOUBLE, 1,&areanorm);
      
      /* Need to put observations on a band */

      for (w=data_bot; w<=data_top; w++){
	area_data += (t1[w+1]-t1[w-1])/2.0 * t2[w];
      }
      static double areaval;
      areaval = area_data;
      nc_put_att_double(ncid1, varid1, "area_under_curve_of_data_in_data_range", NC_DOUBLE, 1,&areaval);
      
      double normalise_area = (t1[data_top]-t1[data_bot])/(ws->wave[model_top]-ws->wave[model_bot]);

      double errorr = 100.0 * (area_model * normalise_area - area_data)/(area_data + 1.0e-16);
      static double errrval;
      errrval = errorr;
      nc_put_att_double(ncid1, varid1, "interpolation_percent_error_in_data_range", NC_DOUBLE, 1,&errrval);
      
      double range[2];
      range[0] = t1[0];
      range[1] = t1[len-1];
      static double ranval[2];
      ranval[0] = range[0];
      ranval[1] = range[1];
      nc_put_att_double(ncid1, varid1, "source_data_range_in_nm", NC_DOUBLE, 2,ranval);
      
      double range2[2];
      static double ranval2[2];
      range2[0] = ws->wave[model_bot];
      range2[1] = ws->wave[model_top];
      ranval2[0] = range2[0];
      ranval2[1] = range2[1];
      nc_put_att_double(ncid1, varid1, "model_range_used_for_interpolation_error_calculation", NC_DOUBLE, 2,ranval2);
      
      double range3[2];
      static double ranval3[2];
      range3[0] = t1[data_bot];
      range3[1] = t1[data_top];
      ranval3[0] = range3[0];
      ranval3[1] = range3[1];
      nc_put_att_double(ncid1, varid1, "data_range_used_for_interpolation_error_calculation", NC_DOUBLE, 2,ranval3);
      
      eco_write_setup(e,"\t Percent error of interpolation onto optical grid (400 - 700 nm) is %e.\n",errorr);

      ncw_enddef("optical_setup.nc", ncid1);
      ncw_put_var_double("optical_setup.nc",ncid1,varid1,out);
      nc_close(ncid1);

      FILE *fp;
      fp = fopen("optical_setup_matlab_check.m", "a");
      fprintf(fp,"figure \n");
      fprintf(fp,"varname = '%s';\n",codename);
      fprintf(fp,"sourcefile = '%s';\n",source_file);
      fprintf(fp,"sourcename = '%s';\n",varname);
      fprintf(fp,"set(gca,'FontSize',14); \n");
      fprintf(fp,"var = ncread(file,varname);\n");
      fprintf(fp,"error = ncreadatt(file,varname,'interpolation_percent_error_in_data_range');\n");
      fprintf(fp,"var1 = ncread(sourcefile,sourcename);\n");
      fprintf(fp,"units = ncreadatt(file,varname,'units');\n");
      fprintf(fp,"desc = ncreadatt(file,varname,'description');\n");
      fprintf(fp,"try \n");
      fprintf(fp," rgb = ncreadatt(file,varname,'rgb_colour_space');\n");
      fprintf(fp,"catch \n");
      fprintf(fp," rgb = [0 0 0]; \n");
      fprintf(fp,"end \n");
      fprintf(fp,"pp = plot(wave,var,'+k');set(pp,'LineWidth',2,'MarkerSize',16);hold on;\n");
      fprintf(fp,"pp = plot(var1(:,1),var1(:,2),'-k');set(pp,'LineWidth',7);hold on;\n");
      fprintf(fp,"pp = plot(var1(:,1),var1(:,2),'-k');set(pp,'LineWidth',5,'Color',rgb/256);hold on;\n");
      fprintf(fp,"xlabel('wavelength [nm]');\n");
      fprintf(fp,"ylabel([varname ,'[',units,']'],'Interpreter','none');\n");
      fprintf(fp,"title([desc,',: Interp. error: ',num2str(error),' %']);\n");
      fprintf(fp,"set(gca,'xlim',[200 1000]);\n\n");
      fprintf(fp,"text(xlim+0.05*diff(xlim),ylim+0.93*diff(ylim),['Source: ',sourcefile],'Interpreter','none');\n"); 
      fclose(fp);

      FILE *fp2;
      fp2 = fopen("optical_setup_python_check.py", "a");
     
      fprintf(fp2,"plt.figure()\n");
      fprintf(fp2,"matplotlib.rc('xtick', labelsize=14)\n");
      fprintf(fp2,"matplotlib.rc('ytick', labelsize=14)\n");
      fprintf(fp2,"var = XR['%s'][:]\n", codename);      
      fprintf(fp2,"sourcename=var.attrs['sourcefile']\n");
      fprintf(fp2,"sourcefile = f'{dpath}/{sourcename}'\n");
      fprintf(fp2,"XR1 = xar.open_mfdataset(sourcefile, decode_times= True)\n");
      fprintf(fp2,"var1 = XR1['%s']\n",varname);
      fprintf(fp2,"try:\n"); 
      fprintf(fp2,"\t rgb=var.attrs['rgb_colour_space']/256\n"); 
      fprintf(fp2,"except:\n");
      fprintf(fp2,"\t rgb = [0, 0, 0];\n");        
      fprintf(fp2,"plt.plot(var1[0,:], var1[1,:],'-k',linewidth=7.0)\n"); 
      fprintf(fp2,"plt.plot(var1[0,:], var1[1,:],'-',linewidth=5.0, color =rgb)\n"); 
      fprintf(fp2,"plt.plot(wave, var,'Xm',linewidth=2.0)\n"); 
      fprintf(fp2,"units=var.attrs['units']\n");  
      fprintf(fp2,"error=var.attrs['interpolation_percent_error_in_data_range']\n");
      fprintf(fp2,"plt.xlabel('wavelength [nm]', size=15)\n");
      fprintf(fp2,"plt.ylabel(f'%s [{units}]', size=15)\n", codename);
      fprintf(fp2,"plt.title(f'Interpolation error: {error}', size=15)\n");
      fprintf(fp2,"plt.gca().set_xlim([200 ,1000])\n");
      fprintf(fp2,"plt.savefig('%s.png', format='png')\n",codename);
      fclose(fp2);
  }else{
    eco_write_setup(e,"optdataspectralread: reading %s but can't find variable %s \n",source_file,varname);
    e_quit("eco:light_spectral_col, optdataspectralread: variable (%s) not in (%s) \n",varname,source_file);
  }
}

void scale_sensor_response(ecology *e, workspace *ws, const char *sensorname, double* in, double* out, double *scale)
{
  
  ///*** THIS FUNCTION HAS NOT BEEN TESTED ***///////////////

  // Scale sensors bands so that sum of responses is 1.00. Magnitude of scaling given. 
  // If scale < 0.5, then derived products using this sensor are unreliable. Whether the product is still calculated 
  //            is a product-specific decision determined in pre_calc. 

  int w;
  double sum = 0.0; 

  for (w=0; w < ws->num_wave; w++){
    sum += in[w];
  }

  *scale = sum;

  for (w=0; w < ws->num_wave; w++){
    out[w] = in[w]/(sum + 1e-10);
  }

  eco_write_setup(e,"scale_sensor_response: adjusting %s by %e \n",sensorname,sum);

}

void chl_specific_aggregate_pigment_absorption(ecology *e,workspace *ws)
{

  // Call pigments that are not required as an unaggregated value

  char *source_file = "csiro_siop_library.nc";
  int ncid;

  ncw_open(source_file, NC_NOWRITE, &ncid);
  
  double* yC_peridinin;
  yC_peridinin = d_alloc_1d(ws->num_wave);
  optdataspectralread(e,ws,source_file,ncid,"Peridinin","Peridinin_wave","Peridinin",yC_peridinin);

  double* yC_chlorophyllc2;
  yC_chlorophyllc2 = d_alloc_1d(ws->num_wave);
  optdataspectralread(e,ws,source_file,ncid,"Chlorophyll\ c2","Chlorophyllc2_wave","Chlorophyll\ c2",yC_chlorophyllc2);
  
  double* yC_chlorophylla;
  yC_chlorophylla = d_alloc_1d(ws->num_wave);
  optdataspectralread(e,ws,source_file,ncid,"Chlorophyll-a","Chlorophyll-a_wave","Chlorophyll-a",yC_chlorophylla);
  
  double* yC_betacarotene;
  yC_betacarotene= d_alloc_1d(ws->num_wave);
  optdataspectralread(e,ws,source_file,ncid,"B\,B-carotene","B\,B-carotene_wave","B\,B-carotene",yC_betacarotene);

  double* yC_zeaxanthin;
  yC_zeaxanthin= d_alloc_1d(ws->num_wave);
  optdataspectralread(e,ws,source_file,ncid,"Zeaxanthin","Zeaxanthin_wave","Zeaxanthin",yC_zeaxanthin);

  double* yC_echinenone;
  yC_echinenone= d_alloc_1d(ws->num_wave);
  optdataspectralread(e,ws,source_file,ncid,"Echinenone","Echinenone_wave","Echinenone",yC_echinenone);

  double* yC_alloxanthin;
  yC_alloxanthin= d_alloc_1d(ws->num_wave);
  optdataspectralread(e,ws,source_file,ncid,"Alloxanthin","Alloxanthin_wave","Alloxanthin",yC_alloxanthin);

  double* yC_fucoxanthin;
  yC_fucoxanthin= d_alloc_1d(ws->num_wave);
  optdataspectralread(e,ws,source_file,ncid,"Fucoxanthin","Fucoxanthin_wave","Fucoxanthin",yC_fucoxanthin);

  double* yC_myxoxanthophyll;
  yC_myxoxanthophyll= d_alloc_1d(ws->num_wave);
  optdataspectralread(e,ws,source_file,ncid,"Myxoxanthophyll","Myxoxanthophyll_wave","Myxoxanthophyll",yC_myxoxanthophyll);

  double* yC_PE;
  yC_PE = d_alloc_1d(ws->num_wave);
  optdataspectralread(e,ws,source_file,ncid,"Phycoerythrin","Phycoerythrin_wave","Phycoerythrin",yC_PE);

  double* yC_PC;
  yC_PC = d_alloc_1d(ws->num_wave);
  optdataspectralread(e,ws,source_file,ncid,"Phycocyanin","Phycocyanin_wave","Phycocyanin",yC_PC);

  ncw_close(source_file,ncid);
  
  // Calculate aggregate_pigment_absorption, eventually from reading netcdf file of mass_specific absorption coefficients.

  int num_wave = ws->num_wave;int w;

  double yC_symbiodinium[num_wave];
  double yC_picoplankton[num_wave];
  double yC_rhodomonas_duplex[num_wave];
  double yC_microplankton[num_wave];
  double yC_tricho[num_wave];
  
  for (w=0; w<num_wave; w++){
    yC_symbiodinium[w]  = yC_chlorophyllc2[w] * 0.1273 + yC_peridinin[w] * 0.4733 + yC_betacarotene[w] * 0.0446 + yC_chlorophylla[w] * 1.0;
  }

  for (w=0; w<num_wave; w++){   // Monika's change:
    yC_picoplankton[w]  = yC_PE[w] * 0.01 + yC_PC[w] * 0.01 + yC_echinenone[w] * 0.00 +  yC_zeaxanthin[w] * 0.46 + yC_betacarotene[w] * 0.09 + yC_chlorophylla[w] * 1.0;
  }
  
  // Lesley - Derwent R - ignore chlorophyllide-a which is a precursor of chlorophyll
  
  for (w=0; w<num_wave; w++){
    yC_rhodomonas_duplex[w]  = (yC_chlorophyllc2[w] * 529.0 + yC_PE[w] * 269.0 + yC_alloxanthin[w] * 568.0 + yC_betacarotene[w] * 125.0 + yC_chlorophylla[w] * 3067.0)/3067.0;
  }

  for (w=0; w<num_wave; w++){
    yC_microplankton[w]  = yC_fucoxanthin[w] * 0.6 + yC_chlorophylla[w] * 1.0;
  }
  
  for (w=0; w<num_wave; w++){
    yC_tricho[w]  = yC_PE[w] * 2.5 + yC_myxoxanthophyll[w] * 0.02 +  yC_zeaxanthin[w] * 0.1 + yC_betacarotene[w] * 0.09 + yC_chlorophylla[w] * 1.0;
  }

  // Now assign to microalgal types.

  for (w=0; w<num_wave; w++){
    ws->yC_s[w] =      max(0.0,yC_picoplankton[w]);
    ws->yC_l[w] =      max(0.0,yC_microplankton[w]);
    ws->yC_Tricho[w] = max(0.0,yC_tricho[w]);
    ws->yC_MPB[w] =    max(0.0,yC_microplankton[w]);
    ws->yC_D[w] =      max(0.0,yC_microplankton[w]);
    ws->yC_Symbiodinium[w] = max(0.0,yC_symbiodinium[w]);
  }
}


void optdataspectralread_fine(ecology *e, workspace *ws, char *source_file, int ncid, const char *varname, const char *wavename, const char *codename, double* out)
{
  // read a spectrally-resolved variable from database, interpolate onto the fine grid, place in workspace and write to setup file.

  int dimid,w,varid;
  int dimids[2];
  size_t len;
  size_t ncstart[2] = {0,0}; size_t nccount[2] = {2,1};

  char units[30];
  char description[100];
  char symbol[30];

  char *output_file = "optical_setup.nc";
  
  varid = ncw_var_id(ncid,varname);

  if (varid > -1){

    eco_write_setup(e,"optdataspectralread_fine: reading %s for variable %s \n",source_file,varname);

    nc_inq_vardimid(ncid, varid, dimids);
    nc_inq_dimlen(ncid, dimids[1],&len);

    double nc_var[2*len];
    ncw_get_att_text(source_file,ncid,varid,"Units",units);
    ncw_get_att_text(source_file,ncid,varid,"Symbol",symbol);
    ncw_get_att_text(source_file,ncid,varid,"Description",description);
    nccount[1] = len;
    ncw_get_vara_double(source_file,ncid, varid, ncstart,nccount,nc_var);
    
    /* now interpolate onto model wavelengths */
    
    double t1[len];
    double t2[len];
    
    for (w=0; w<len; w++){
      t1[w] = nc_var[w];
      t2[w] = nc_var[w+len];
    }
    
    interp1d(t1, t2, len, ws->sensor_waves, out, ws->num_sensor_waves);

    eco_write_setup(e,"Write %s: %s [%s] to %s \n",description,symbol,units,output_file);
    eco_write_setup(e,"wave [nm] %s \n",varname);

    for (w=0; w<ws->num_sensor_waves; w++){
      eco_write_setup(e,"%4.2f \t %e \n",ws->sensor_waves[w],out[w]);
    }

    // now open output netcdf and put in variable.

      int ncid1,varid1,dims1,dimid1;
      int dims[2];
      
      ncw_open("optical_setup.nc",NC_WRITE, &ncid1);
      dimid1 = ncw_dim_id(ncid1,"wave_sensor");
      dims[0] = dimid1;
      ncredef(ncid1);
      ncw_def_var("optical_setup.nc",ncid1,codename,NC_DOUBLE,1,dims,&varid1);
      nc_put_att_text(ncid1, varid1, "description", strlen(description),description);
      nc_put_att_text(ncid1, varid1, "codename", strlen(codename),codename);
      nc_put_att_text(ncid1, varid1, "sourcefile", strlen(source_file),source_file); // would be good to have creation_date of file
      nc_put_att_text(ncid1, varid1, "sourcename", strlen(varname),varname);
      nc_put_att_text(ncid1, varid1, "symbol", strlen(symbol),symbol);
      nc_put_att_text(ncid1, varid1, "units", strlen(units),units);

      // Now do test of error of interpolation onto model optical grid using area under the curve.

      double area_model = 0.0;
      double area_data = 0.0;

      int model_top = -1;
      int model_bot = -1;

      // find top and bottom model wavelengths in the data range for interpolation error calculation.
      // need to revise for sensor waves

      for (w=0; w<ws->num_sensor_waves; w++){
        if (ws->sensor_waves[w] >= t1[0] && model_bot < 0 && ws->sensor_waves[w] >= 400.0){
          model_bot = w;
        }
        if ((ws->sensor_waves[w] >= t1[len-1] || ws->sensor_waves[w] >= 700.0 ) && model_top < 0 ){
          model_top = w-1;
        }
      }

      if (model_top == -1){
	model_top = ws->num_sensor_waves-1;
      }

      // Now find data points closest to the choosen model range.

      int data_top = -1;
      int data_bot = -1;

      for (w=0; w<len; w++){
        if (t1[w] >= ws->sensor_waves[model_bot] && data_bot < 0){
          data_bot = w;
        }
        if (t1[w] >= ws->sensor_waves[model_top] && data_top < 0){
          data_top = w-1;
        }
      }

      if (data_top == -1){
	data_top = len-1;
      }

      // make sure it is defined.

      if (data_top < 0){
	eco_write_setup(e,"light_spectral_col: THIS SHOULD NEVER HAPPEN  %d \n",data_top);
	data_top = data_bot;
      }

      for (w=model_bot; w <= model_top; w++){
        area_model += out[w]; /// assuming ws->bandwidth[w] is 1.
      }

      static double areanorm;
      areanorm = area_model;
      nc_put_att_double(ncid1, varid1, "area_under_curve_of_model_in_data_range", NC_DOUBLE, 1,&areanorm);
      
      /* Need to put observations on a band */

      for (w=data_bot; w<=data_top; w++){
	area_data += (t1[w+1]-t1[w-1])/2.0 * t2[w];
      }
      static double areaval;
      areaval = area_data;
      nc_put_att_double(ncid1, varid1, "area_under_curve_of_data_in_data_range", NC_DOUBLE, 1,&areaval);
      
      double normalise_area = (t1[data_top]-t1[data_bot])/(ws->sensor_waves[model_top]-ws->sensor_waves[model_bot]);

      double errorr = 100.0 * (area_model * normalise_area - area_data)/(area_data + 1.0e-16);
      static double errrval;
      errrval = errorr;
      nc_put_att_double(ncid1, varid1, "interpolation_percent_error_in_data_range", NC_DOUBLE, 1,&errrval);
      
      double range[2];
      range[0] = t1[0];
      range[1] = t1[len-1];
      static double ranval[2];
      ranval[0] = range[0];
      ranval[1] = range[1];
      nc_put_att_double(ncid1, varid1, "source_data_range_in_nm", NC_DOUBLE, 2,ranval);
      
      double range2[2];
      static double ranval2[2];
      range2[0] = ws->sensor_waves[model_bot];
      range2[1] = ws->sensor_waves[model_top];
      ranval2[0] = range2[0];
      ranval2[1] = range2[1];
      nc_put_att_double(ncid1, varid1, "model_range_used_for_interpolation_error_calculation", NC_DOUBLE, 2,ranval2);
      
      double range3[2];
      static double ranval3[2];
      range3[0] = t1[data_bot];
      range3[1] = t1[data_top];
      ranval3[0] = range3[0];
      ranval3[1] = range3[1];
      nc_put_att_double(ncid1, varid1, "data_range_used_for_interpolation_error_calculation", NC_DOUBLE, 2,ranval3);
      
      // eco_write_setup(e,"\t Percent error of interpolation onto optical grid (400 - 700 nm) is %e.\n",errorr);

      ncw_enddef("optical_setup.nc", ncid1);
      ncw_put_var_double("optical_setup.nc",ncid1,varid1,out);
      nc_close(ncid1);

      FILE *fp;
      fp = fopen("optical_setup_matlab_check.m", "a");
      fprintf(fp,"figure \n");
      fprintf(fp,"varname = '%s';\n",codename);
      fprintf(fp,"sourcefile = '%s';\n",source_file);
      fprintf(fp,"sourcename = '%s';\n",varname);
      fprintf(fp,"set(gca,'FontSize',14); \n");
      fprintf(fp,"var = ncread(file,varname);\n");
      fprintf(fp,"error = ncreadatt(file,varname,'interpolation_percent_error_in_data_range');\n");
      fprintf(fp,"var1 = ncread(sourcefile,sourcename);\n");
      fprintf(fp,"units = ncreadatt(file,varname,'units');\n");
      fprintf(fp,"desc = ncreadatt(file,varname,'description');\n");
      fprintf(fp,"try \n");
      fprintf(fp," rgb = ncreadatt(file,varname,'rgb_colour_space');\n");
      fprintf(fp,"catch \n");
      fprintf(fp," rgb = [0 0 0]; \n");
      fprintf(fp,"end \n");
      fprintf(fp,"pp = plot(wave_sensor,var,'+k');set(pp,'LineWidth',2,'MarkerSize',16);hold on;\n");
      fprintf(fp,"pp = plot(var1(:,1),var1(:,2),'-k');set(pp,'LineWidth',7);hold on;\n");
      fprintf(fp,"pp = plot(var1(:,1),var1(:,2),'-k');set(pp,'LineWidth',5,'Color',rgb/256);hold on;\n");
      fprintf(fp,"xlabel('wavelength [nm]');\n");
      fprintf(fp,"ylabel([varname ,'[',units,']'],'Interpreter','none');\n");
      fprintf(fp,"title([desc,',: Interp. error: ',num2str(error),' %']);\n");
      fprintf(fp,"set(gca,'xlim',[200 1000]);\n\n");
      fprintf(fp,"text(xlim+0.05*diff(xlim),ylim+0.93*diff(ylim),['Source: ',sourcefile],'Interpreter','none');\n"); 
      fclose(fp);

      FILE *fp2;
      fp2 = fopen("optical_setup_python_check.py", "a");
     
      fprintf(fp2,"plt.figure()\n");
      fprintf(fp2,"matplotlib.rc('xtick', labelsize=14)\n");
      fprintf(fp2,"matplotlib.rc('ytick', labelsize=14)\n");
      fprintf(fp2,"var = XR['%s'][:]\n", codename);      
      fprintf(fp2,"sourcename=var.attrs['sourcefile']\n");
      fprintf(fp2,"sourcefile = f'{dpath}/{sourcename}'\n");
      fprintf(fp2,"XR1 = xar.open_mfdataset(sourcefile, decode_times= True)\n");
      fprintf(fp2,"var1 = XR1['%s']\n",varname);
      fprintf(fp2,"try:\n"); 
      fprintf(fp2,"\t rgb=var.attrs['rgb_colour_space']/256 \n"); 
      fprintf(fp2,"except:\n");
      fprintf(fp2,"\t rgb = [0, 0, 0];\n");        
      fprintf(fp2,"plt.plot(var1[0,:], var1[1,:],'-k',linewidth=7.0)\n"); 
      fprintf(fp2,"plt.plot(var1[0,:], var1[1,:],'-',linewidth=5.0, color =rgb)\n"); 
      fprintf(fp2,"plt.plot(wave, var,'Xm',linewidth=2.0)\n"); 
      fprintf(fp2,"units=var.attrs['units']\n");  
      fprintf(fp2,"error=var.attrs['interpolation_percent_error_in_data_range']\n");
      fprintf(fp2,"plt.xlabel('wavelength [nm]', size=15)\n");
      fprintf(fp2,"plt.ylabel(f'%s [{units}]', size=15)\n", codename);
      fprintf(fp2,"plt.title(f'Interpolation error: {error}', size=15)\n");
      fprintf(fp2,"plt.gca().set_xlim([200 ,1000])\n");
      fprintf(fp2,"plt.savefig('%s.png', format='png')\n",codename);
      fclose(fp2);
  }else{
    eco_write_setup(e,"optdataspectralread: reading %s but can't find variable %s \n",source_file,varname);
    e_quit("eco:light_spectral_col, optdataspectralread: variable (%s) not in (%s) \n",varname,source_file);
  }
}

