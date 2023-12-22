/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: lib/ecology/ecofunct.c
 *  
 *  Description:
 *  Functions to be used in ecological processes --
 *  implementation.
 *  
 *  Notes:          Much of the theory behind these functions can be found at:
 *  http://www.oikos.warwick.ac.uk/ecosystems/ThesisArchive/
 *  baird_thesis.html
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: ecofunct.c 7439 2023-10-27 03:06:27Z bai155 $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#if defined(__sparc)
#include <ieeefp.h>
#endif
#include "ems.h"
#include "netcdf.h"
#include "ncw.h"
#include "utils.h"
#include "constants.h"
#include "ecofunct.h"

#define NMAX 4
#define PLANKGROWRATE_MIN 1.0e-15

// See grow_data.c
size_t extract_grow_from_defaults(int n, char* name,
				  double** ptable, double** prate);

/** Opens a netcdf file and extracts a matrix of the axes and the
 * n-dimensional look-up tables
 */
int extract_grow(int n, char* file_name, double** ptable, double** prate)
{
    int ncid, kay1id, norm_id, rate_id;
    size_t rate_len;
    FILE *fp = NULL;

    // See if file_name is actually a file
    fp = fopen(file_name, "r");
    if (fp != NULL) {
      // It is - so proceed by reading this netcdf file
      fclose(fp);
      
      /*
       * open the appropriate netcdf file which has a standard form
       */
      ncw_open(file_name, NC_NOWRITE, &ncid);
      
      /*
       * determine the length of each dimension
       */
      ncw_inq_dimid(file_name, ncid, "kay1", &kay1id);
      ncw_inq_dimlen(file_name, ncid, kay1id, &rate_len);
      
      *ptable = (double*) calloc((int) pow(rate_len, n) + 1, sizeof(double));
      *prate  = (double*) calloc(rate_len * n + 1, sizeof(double));
      
      /*
       * extract the variables that require interpolation
       */
      ncw_inq_varid(file_name, ncid, "normgrow", &norm_id);
      ncw_inq_varid(file_name, ncid, "rates", &rate_id);
      
      /*
       * place values into a local variable
       */
      ncw_get_var_double(file_name, ncid, norm_id, *ptable);
      ncw_get_var_double(file_name, ncid, rate_id, *prate);
      
      ncw_close(file_name, ncid);
    
    } else {
      // Get from defaults
      rate_len = extract_grow_from_defaults(n, file_name, ptable, prate);
    }
    
    return (int)rate_len;
}

/** Takes all the cell-specific parameters and environmental variables,
 * calculates the physical limits, and sends these to CRgrowth, which returns
 * an approximate growth rate.
 */
double plankgrow(double umax, double aA, double psi, double Sh, double diff[], double m, double I, double conc[], int len, int n, double* grow, double* rate, double Plank_resp)
{
    double rates[NMAX];
    double plankgrowrate;
    int i;

    rates[0] = PlankMaxLiteUp(aA, I);   /* photon cell-1 s-1 */
    for (i = 1; i < n; i++)
        /*
         * mol * cell-1 * s-1
         */
        rates[i] = PlankMaxNutUp(psi, diff[i - 1], conc[i - 1], Sh);
    /*
     * growth rate in s-1
     */
    plankgrowrate = CRgrowth(umax, rates, m, grow, rate, len, n, Plank_resp);

    if (plankgrowrate < PLANKGROWRATE_MIN)
	plankgrowrate = 0.0;

    /*
     * do not exit yet, give the host model a chance to report the location
     */
  if (!finite(plankgrowrate))
    emslog(LWARN, "FYI: ecology: plankgrow(): plankgrowrate = inf or NaN\n");

    return plankgrowrate;
}

/** Takes all the cell-specific parameters and environmental variables,
 * calculates the physical limits, and sends these to CRgrowth, which returns
 * an approximate growth rate.
 *
 * Note the calculation is done for growth per m2. Implicitly we are assuming
 * that the particular plant covers the whole bottom. Unlike cells,
 * stoichiometry, m, is a function of biomass.
 */
double benthgrow(double umax, double cf, double vis, double Ub, double ks, double diff[], double m, double I, double conc[], int len, int n, double* grow, double* rate, double Benth_resp)
{
    double rates[NMAX];
    double benthgrowrate;
    int i;

    rates[0] = I;               /* photons m-2 s-1 */
    for (i = 1; i < n; i++)
        /*
         * mol * m-2 * s-1
         */
        rates[i] = BenthicMaxNutUpSmooth(cf, vis, diff[i - 1], conc[i - 1], Ub);
    /*
     * growth rate in s-1
     */
    benthgrowrate = CRgrowth(umax, rates, m, grow, rate, len, n, Benth_resp);

    /*
     * do not exit yet, give the host model a chance to report the location
     */
  if (!finite(benthgrowrate))
    emslog(LWARN, "FYI: ecology: benthgrow(): benthgrowrate = inf or NaN\n");

    return benthgrowrate;
}

/** Multipicative growth function see Bormans and Webster 1999, JPR 31:581-598.
 * Assumes that rates are unitless fractions, obtained from a rectangular
 * hyperbolic function see RectHyp.
 */
double MPgrowth(double umax, double uptake_rates[], int n)
{
    double output;
    int i;

    output = umax;
    for (i = 0; i < n; i++)
        output *= uptake_rates[i];

    return output;
}

/** Law of Minimum growth function -- see PPB.
 * Assumes that rates are unitless fractions, obtained from a rectangular
 * hyperbolic function see RectHyp.
 */
double LMgrowth(double umax, double uptake_rates[], double stoich[], int n)
{
    double output, temp;
    int i;

    output = 1.0;
    for (i = 0; i < n; i++) {
        temp = uptake_rates[i];
        if (output > temp)
            output = temp;
    }
    output *= umax;

    return output;
}

/** Calculates the maximum growth rate at a particular temperature.
 */
double UmaxatT(double umax, double Q10, double T, double Tref)
{

    return umax * pow(Q10, (T - Tref) / 10.0);
}

/** Calculates natural mortality.
 * Note: Natural mortality is very small in the laboratory, and in the field
 *  as well, I would suggest. It may be used to account for loss processes not
 * explicitly modelled (such as viral cell death, toxicity etc.), but this
 * results in a loss of predictive capabilities. To be used only in the
 * calibration stage.
 * @param Pop population suffering death [cell m-3]
 * @param rate rate coefficient [s-1]
 * @return death rate [cell s-1]
 */
double NaturalMortality(double Pop, double rate)
{
    return Pop * rate;
}

/** Calculates dynamic viscosity [N s m-2]
 *
 * Wolf-Gladrow has an ugly expression for viscosity in fresh water:
 * visfresh = 1.002e-3 * exp((1.1709 * (20.0 - T)
 *            - 1.827e-3 * (T - 20.0) * (T - 20.0)) / (20.0 + 89.93))
 */
double Viscosity(double T, double S)
{
    double visfresh = 0.00123161 - 1.228e-5 * T + 4.009e-8 * T * T;
    double vis35 = visfresh / (0.9508 - 0.0007379 * T); /* viscosity at S =
                                                         * 35 */
    S /= 35.0;
    return vis35 * S + visfresh * (1.0 - S);
}

/** Calculates Molecular Diffusivity.
 * Valid between 0 and 25 C.
 * References: Baird, 1999 PhD thesis, Appendix A.
 *             Wolf-Gladrow and Riebesell(1997) Marine Chem 59:17-34.
 *             Li and Gregory, (1974) Geochim at Cos. Acta 38:703-714.
 * @param D25S0 diffusivity at 298 K C, S = 0 [m2 s-1]
 * @param T temperature [deg C]
 * @param S salinity (ppt)
 * $return Moleculkar Diffusivity [m2 s-1]
 */
double MolDiff(double D25S0, double T, double S)
{
    /*
     * calc viscosity at T, S =0
     */
    double vissalty = Viscosity(T, S);
    double visfresh = Viscosity(T, 0);

    /*
     * calc mol diffusivity at T, S
     */
    double DxxSyy = D25S0 * vissalty * (T + 273.15) / (visfresh * 298.15);

    return DxxSyy;
}

/** Calculates maximum nutrient uptake (per unit area). Based on the
 * dimensionless friction coefficient, cf, and Schmidt No. for a smooth
 * surface. Stm for a smooth surface should be ~2e-5, so if Ub = 0.01m,
 * Ub/Stm~500 (or a boundary layer thickness of 0.05 um.
 */
double BenthicMaxNutUpSmooth(double cf, double vis, double D, double N, double Ub)
{
    /*
     * Schmidt number ( = momentum transport / mass transport ). Note that
     * this is kinematic viscosity (m-2 s-1) not dynamic viscosity (N s
     * m-2).
     */
    double Sc = vis / D;

    /*
     * Smooth surface Stanton number, Baird and Atkinson (1997).
     */
    double Stm = sqrt(cf / 2.0) * (0.0575 * pow(Sc, -0.666667) + 0.1184 / Sc);

    return N * Ub * Stm;
}

/* From Baird and Atkinson, 1997, LO 42: 1685-1693.
 *
 * This is based on a heat transfer correlation. The basis is that mass
 * (or heat transfer) is a function of energy dissipated on a surface, which
 * is quantified using at a particular ks and cf. Stm for coral is 3e-4.
 * Be careful to make sure that Rek > 70, or the surface is not hydraulically
 * rough.
 */
double BenthicMaxNutUpRough(double cf, double vis, double D, double N, double Ub, double ks)
{
    /*
     * friction velocity [m s-1]
     */
    double ustar = Ub * sqrt(cf / 2.0);

    /*
     * Reynold's roughness number [-]
     */
    double Rek = ustar * ks / vis;

    /*
     * Schmidt number [-]
     */
    double Sc = vis / D;

    /*
     * roughness Stanton number
     */
    double Stk = 1.0 / (5.19 * pow(Rek, 0.2) * pow(Sc, 0.44) - 8.48);
    double Stm = (cf / 2.0) / (0.9 + (sqrt(cf / 2.0) / Stk));

    return N * Ub * Stm;
}

/** Blanch et. al., 1998 Aquatic Botany 61:181-205.
 * @param Kb specific absorption coefficient [m2 g(biomass)-1]
 * @param I irradiance mol(photons) [m-2 s-1]
 * @return maxlite [mol(photons) m-2 s-1]
 */
double BenthicMaxLiteUp(double I, double Kb)
{
    return I * Kb;
}

/*  Using a Rectangular Hyperbolae to model uptake rates, as first
 *  proposed by Monod in 1942. Most ecological models set maximum to 1,
 *  so that this term becomes a fraction. A more mechanistic approach set
 *  a maximum uptake rate, but this is rarely used in ecological models.
 */
double RectHyp(double maximum, double halfsat, double conc)
{
    return maximum * conc / (halfsat + conc);
}

/** Calculates maximum nutrient uptake by solving the Laplace equation and
 * accounting for relative fluid motion using Sherwood number [see Baird and
 * Emsley 1999 (JPR 21:85-126 Eq. 15)]
 * @param psi diffusion shape factor [m] = 4*pi*r for a sphere (Table II)
 * @param D molecular diffusivity of nutrient [m2 s-1], see Appendix A of
 * Baird, 1999
 * @param Sh Sherwood No = Advective Flux / Diffusive Flux (Table III)
 * @param N nutrient concentration [mol m-3]
 * @return output Uptake [mol s-1 cell-1]
 */
double PlankMaxNutUp(double psi, double D, double N, double Sh)
{
    return psi * D * N * Sh;
}

/** Calculates maximum light uptake from the absorption cross-section and the
 * incident radiation [see Baird and Emsley 1999 (JPR 21:85-126 Eq. 17)]
 * @param aA abosorption cross-section. Can be found as a function of pigment
 * concentration and shape [m-2] [Table IV]
 * @param I irradiance mol(photons) [m-2 s-1]
 * @return maxlite [mol(photons) s-1 cell-1]
 */
double PlankMaxLiteUp(double aA, double I)
{
    return aA * I;
}

/* Decays each element in labile detritus at a constant rate decay, sending one
 * fraction dissolved inorganic nutrients, the rest to the refractory pool
 */
double LDtoRD(double decay, double frac, int n)
{
    return 0.0;
}

/** Calculates the encounter rate of two spheres with individual swimming and
 * sinking velocities in a fluid with a turbulence intensity specified by a
 * mean dissipation rate of TKE. A discription of the theory can be found in
 * Mark Baird's PhD thesis.
 * @param F - method of determining encounter rates:
 *            'r' - Rectilinear coagulation Jackson (1995)
 *            'c' - Curvilinear coagulation Jackson (1995)
 * @param rP radii of prey [m]
 * @param Pswim swimming velocity of phyto [m s-1]
 * @param Psink ? [m s-1]
 * @param rH radius of predator [m]
 * @param Hswim swimming velocity of zooplankton [m s-1]
 * @param Hsink terminal sinking velocity of zooplankton [m s-1]
 * @param epsilon dissipation of turbulent kinetic energy [m2 s-3]
 * @param viscosity kinematic viscosity [m2 s-1]
 * @param density density [kg m-3]
 * @param Twater temperature of the water [C]
 * @return encounter rate
 */
double phi(char* F, double rP, double Pswim, double Psink, double rH, double Hswim, double Hsink, double epsilon, double viscosity, double density, double Twater)
{
    double Uswim, Ueff, p;
    double phiswim = 0.0;
    double phishear = 0.0;
    char* rect = "rect";
    char* curv = "curv";

		if(rP == 0.0)
			return 0.0;
    /*
     * encounter rate based on diffusion
     */
    double phidiff = (2.0 * BOLT * Twater / (3.0 * density * viscosity)) * (1.0 / rP + 1.0 / rH) * (rP + rH);

    /*
     * Now combine the effects of sinking and swimming on encounter rates.
     * Assume differential sedimentation in z-direction, and swimming
     * velocities of predator and prey are uncorrelated to each other, or
     * the vertical (i.e. sedimentation direction).
     */
    double Usink = Hsink - Psink;       /* along same direction */

    /*
     * effective swimming encounter velocity (Gerritsen and Strickler, 1977)
     */
    if (Hswim > Pswim) {
        Uswim = (Pswim * Pswim + 3.0 * Hswim * Hswim) / (3.0 * Hswim);
    } else
        Uswim = (Hswim * Hswim + 3.0 * Pswim * Pswim) / (3.0 * Pswim);

    /*
     * effective encounter velocity (both swimming and sinking)
     */
    if (Usink > Uswim)
        Ueff = (Uswim * Uswim + 3.0 * Usink * Usink) / (3.0 * Usink);
    else
        Ueff = (Usink * Usink + 3.0 * Uswim * Uswim) / (3.0 * Uswim);

    /*
     * rectilinear coagulation coefficient calculations
     */
    if (strcmp(F, rect) == 0) {
        phiswim = M_PI * (pow((rP + rH), 2)) * Ueff;
        phishear = 1.3 * (pow((epsilon / viscosity), 0.5)) * pow((rP + rH), 3.0);
    }
    /*
     * curvilinear coagulation coeffiecient calculations
     */
    else if (strcmp(F, curv) == 0) {
        p = rP / rH;
        phiswim = 0.5 * M_PI * (pow(rP, 2)) * Ueff;
        phishear = 9.8 * pow((p / (1.0 + 2.0 * p)), 2) * pow((rP + rH), 3) * sqrt(epsilon / viscosity);
    } else
        e_quit("error: ecology: phi(): unknown method for calculating encounter rates\n");

    return phishear + phiswim + phidiff;        /* m3 s-1 */
}

/** Calculates the absorption cross-section [m-2] of a sphere with a specified
 * radius [m] and absorbance [m].
 */
double aa(double rad, double absorb)
{
  double temp = 2.0 * absorb * rad;
  double ans = 0.0;
  /* kirk's equation behaves badly for temp < 3e-4 */ 
  if (temp > 3e-4){
    ans = M_PI * rad * rad * (1.0 - 2.0 * (1.0 - (1.0 + temp) * exp(-temp)) / (temp * temp));
  }
  return ans;
}

void aawave(double rad, double *absorb, double *aA, int num_wave, double *wave)

/* replaced 0.0 with analytical approx. for < temp on 27th Oct 2023 */ 
  
{
  int w;
  double temp;
  for (w=0; w<num_wave; w++){
    temp = 2.0 * absorb[w] * rad;
    aA[w] = M_PI * rad * rad * 2.0 / 3.0 * temp; // approx for < 0.001
      /* Duysens (1956) equation behaves badly for temp < 3e-4 */ 
    if (temp > 0.001){
      aA[w] =  M_PI * rad * rad * (1.0 - 2.0 * (1.0 - (1.0 + temp) * exp(-temp)) / (temp * temp));
    }
  }
}

/** Calculates the absorption cross-section [m-2] of a sphere with a specified
 *   radius [m] and Chlorophyll concentration 
 *
 *      yC - array of spec. absorpt coeff. at discrete wavelength
 *              between 350 and 700 nm.                [m-2 mg(pig)-1]
 *      chla - conc of pigments         [mg(pig) m-3]
 *      wave - array of wavelengths at which aA is calc.        [nm]

 *      Data for absorption cross section without package effect.
 *      from Hoepffner and Sathyendranath (1991).
 */

void absorbwave(double *yC, double chla, int num_wave, double *wave)
{
  int w,i; 
  double v[6][3] = {
    {384., 53.8, 0.037},
    {413., 21.3, 0.012},
    {435., 32.1, 0.039},
    {623., 35.0, 0.005},
    {676., 21.6, 0.020},
    {700., 33.5, 0.002}}; 
  
  for (w=0; w<num_wave; w++) {
    yC[w]=0.0;
    for (i=0; i<6; i++) {
      //  if (wave[w] <= 800.0)
	yC[w] += v[i][2]*chla*exp(-(pow(wave[w]-v[i][0],2.0))/(2.0*v[i][1]*v[i][1]));
    }
  } 
}

void absorbxanth(double *yC, double xanth, int num_wave, double *wave)
{
  int w,i; 
  double v[2][3] = {
    {490., 17.1, 0.0313},
    {532., 22.8, 0.0194}}; 
  
  for (w=0; w<num_wave; w++) {
    yC[w]=0.0;
    for (i=0; i<2; i++) {
      //  if (wave[w] <= 800.0)
	yC[w] += v[i][2]*xanth*exp(-(pow(wave[w]-v[i][0],2.0))/(2.0*v[i][1]*v[i][1]));
    }
  } 
}

void absorbxanth_heat(double *yC, double xanth, int num_wave, double *wave)
{
  int w,i; 
  double v[3][3] = {
    {451., 32.0, 0.0632},
    {464., 8.60, 0.0253},
    {493., 12.0, 0.0464}}; 
  
  for (w=0; w<num_wave; w++) {
    yC[w]=0.0;
    for (i=0; i<3; i++) {
      //  if (wave[w] <= 800.0)
	yC[w] += v[i][2]*xanth*exp(-(pow(wave[w]-v[i][0],2.0))/(2.0*v[i][1]*v[i][1]));
    }
  } 
}

/** Calculates the diffusion shape factor for a sphere.
 */
double psi(double rad)
{
    double psi = 4.0 * M_PI * rad;

    return psi;
}

/* Calculates growth rate after comparing the rates of nutrient uptake, light
 * and the maximum growth rate (relative to the cell's stoichiometry).
 */
double CRgrowth(double umax, double uptake_rates[], double m, double* normgrow, double* rate, int k, int n, double resp)
{
    double output = 0.0;        /* where to put it */
    double frac[NMAX + 1];
    int index[NMAX + 1];
    int i;

    uptake_rates[0] = uptake_rates[0] / (umax * m) - resp;

    for (i = 1; i < n; i++)
        uptake_rates[i] /= umax * m;

    /*
     * exit early with zero growth rate if any supply rates are 0
     */
    for (i = 0; i < n; ++i)
    {
        if (uptake_rates[i] <= 0.0 )/*|| !finite(uptake_rates[i])) FR please take this out - could bury something nasty here*/
            return 0.0;
    }

    for (i = 0; i < n; i++) {
        double maxrate = rate[(i + 1) * k - 1];

        /*
         * set uptake rates, index and frac to maximum uptake rates if they are
         * greater than max rate in lookup table
         */
        if (uptake_rates[i] > maxrate) {
            uptake_rates[i] = maxrate;
            index[i] = k - 1;
            frac[i] = 1.0;
        }
        /*
         * otherwise leave uptake_rates alone, and determine index and frac
         * by a binary search
         */
        else {
            int ilow = 0;
            int ihigh = k - 1;
            int ki = k * i;

            while (ihigh - ilow > 1) {
                int imid = (ilow + ihigh) / 2;

                if (uptake_rates[i] >= rate[imid + ki])
                    ilow = imid;
                else
                    ihigh = imid;
            }
            index[i] = ihigh;
            ki += ihigh;
            frac[i] = (uptake_rates[i] - rate[ki - 1]) / (rate[ki] - rate[ki - 1]);
        }
    }

    if (n == 1) {
        int i0 = index[0];
        double f0 = frac[0];

        output = normgrow[i0 - 1] * (1.0 - f0) + normgrow[i0] * f0;
    } else if (n == 2) {
        int i0 = index[0];
        int i1 = index[1];
        double f0 = frac[0];
        double f1 = frac[1];

        int i01 = i0 - 1;
        int i11 = i1 - 1;

        double f01 = 1.0 - f0;
        double f11 = 1.0 - f1;

        output = f01 * f11 * normgrow[i01 * k + i11]
            + f01 * f1 * normgrow[i01 * k + i1]
            + f0 * f11 * normgrow[i0 * k + i11]
            + f0 * f1 * normgrow[i0 * k + i1];
    } else if (n == 3) {
        int i0 = index[0];
        int i1 = index[1];
        int i2 = index[2];
        double f0 = frac[0];
        double f1 = frac[1];
        double f2 = frac[2];
        int k2 = k * k;

        int i01 = i0 - 1;
        int i11 = i1 - 1;
        int i21 = i2 - 1;
        double f01 = 1.0 - f0;
        double f11 = 1.0 - f1;
        double f21 = 1.0 - f2;

        int i01k2 = i01 * k2;
        int i11k = i11 * k;
        int i0k2 = i0 * k2;
        int i1k = i1 * k;

        output = f01 * f11 * f21 * normgrow[i01k2 + i11k + i21]
            + f01 * f11 * f2 * normgrow[i01k2 + i11k + i2]
            + f01 * f1 * f21 * normgrow[i01k2 + i1k + i21]
            + f01 * f1 * f2 * normgrow[i01k2 + i1k + i2]
            + f0 * f11 * f21 * normgrow[i0k2 + i11k + i21]
            + f0 * f11 * f2 * normgrow[i0k2 + i11k + i2]
            + f0 * f1 * f21 * normgrow[i0k2 + i1k + i21]
            + f0 * f1 * f2 * normgrow[i0k2 + i1k + i2];
    } else if (n == 4) {
        int i0 = index[0];
        int i1 = index[1];
        int i2 = index[2];
        int i3 = index[3];
        double f0 = frac[0];
        double f1 = frac[1];
        double f2 = frac[2];
        double f3 = frac[3];
        int k2 = k * k;
        int k3 = k2 * k;

        int i01 = i0 - 1;
        int i11 = i1 - 1;
        int i21 = i2 - 1;
        int i31 = i3 - 1;
        double f01 = 1.0 - f0;
        double f11 = 1.0 - f1;
        double f21 = 1.0 - f2;
        double f31 = 1.0 - f3;

        int i01k3 = i01 * k3;
        int i0k3 = i0 * k3;
        int i11k2 = i11 * k2;
        int i1k2 = i1 * k2;
        int i21k = i21 * k;
        int i2k = i2 * k;

        int i0111 = i01k3 + i11k2;
        int i2131 = i21k + i31;
        int i213_ = i21k + i3;
        int i2_31 = i2k + i31;
        int i2_3_ = i2k + i3;
        int i011_ = i01k3 + i1k2;
        int i0_11 = i0k3 + i11k2;
        int i0_1_ = i0k3 + i1k2;

        double f0111 = f01 * f11;
        double f2131 = f21 * f31;
        double f213_ = f21 * f3;
        double f2_31 = f2 * f31;
        double f2_3_ = f2 * f3;
        double f011_ = f01 * f1;
        double f0_11 = f0 * f11;
        double f0_1_ = f0 * f1;

        output = f0111 * f2131 * normgrow[i0111 + i2131]
            + f0111 * f213_ * normgrow[i0111 + i213_]
            + f0111 * f2_31 * normgrow[i0111 + i2_31]
            + f0111 * f2_3_ * normgrow[i0111 + i2_3_]
            + f011_ * f2131 * normgrow[i011_ + i2131]
            + f011_ * f213_ * normgrow[i011_ + i213_]
            + f011_ * f2_31 * normgrow[i011_ + i2_31]
            + f011_ * f2_3_ * normgrow[i011_ + i2_3_]
            + f0_11 * f2131 * normgrow[i0_11 + i2131]
            + f0_11 * f213_ * normgrow[i0_11 + i213_]
            + f0_11 * f2_31 * normgrow[i0_11 + i2_31]
            + f0_11 * f2_3_ * normgrow[i0_11 + i2_3_]
            + f0_1_ * f2131 * normgrow[i0_1_ + i2131]
            + f0_1_ * f213_ * normgrow[i0_1_ + i213_]
            + f0_1_ * f2_31 * normgrow[i0_1_ + i2_31]
            + f0_1_ * f2_3_ * normgrow[i0_1_ + i2_3_];
    } else
        e_quit("error: ecology: CRgrowth(): too many dimensions. Expected 4 or less nutrients\n");

    return output * umax;
}

double PhyCellMass(double r)
{
    /*
     * Straile (1997) **** gives molP/cell  KWA
     */
    return 9.14e3 * 4.0 / 3.0 * M_PI * r * r * r / 106.0;
}

double PhyCellChl(double r)
{
    /*
     * Finkel (2001)  **** gives mg Chl m-3  - MB
     */

    double vol = 4.0 / 3.0 * M_PI * r * r * r;

    return 2.0*2.09e7 * pow(1.0e18*vol,-0.310);
}

double diffaa(double yC,double r)
{
    /*
     * Differential of the absorption-cross section of a sphere 
     * w.r.t the absorbance. This represents the incremental value of adding 
     * chlorophyll to increasing absorption.  - MB
     */
  double x = yC*r;

  /*  (1.0-exp(-2.0*x)*(2.0*x+1.0))/(x*x*x) - (2.0*exp(-2.0*x)*(2.0*x+1.0)-2.0*exp(-2.0*x))/(2.0*x*x);*/

  /* algebraic reduction */

  return (1.0-exp(-2.0*x)*(2.0*x*x+2.0*x+1.0))/(x*x*x);

}

double ZooCellMass(double r)
{
    /*
     * Hansen (1997) - Table 1 0.126 g C cm-3 = 10.491e3 mmol C m-3: function returns mol P cell-1
     */
    return 10.5e3 * 4.0 / 3.0 * M_PI * r * r * r / 106.0;
}
