/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/math/underwater.c
 *
 *  \brief standard calculations for underwater properties
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: underwater.c 5833 2018-06-27 00:21:35Z riz008 $
 */

#include <math.h>
#include "ems.h"
#include "underwater.h"

#ifndef GRAV_CONST1
#define GRAV_CONST1 5.2788E-3
#endif

#ifndef GRAV_CONST2
#define GRAV_CONST2 2.36E-5
#endif


#ifndef GRAV_CONST3
#define GRAV_CONST3 9.780318
#endif


/* default sound function 
 */
ssfn speed_of_sound = unesco_speed_of_sound;


/**
 * Calculate the gravity at a location at Latitude 'lat'
 * 
 * @return the Latitude dependend gravity 
 */ 
double gravity(double lat)
{
	double sinTh = sin(lat);
		
	return GRAV_CONST3 * (1 + GRAV_CONST1 *(sinTh * sinTh) + ( GRAV_CONST2 * (sinTh * sinTh * sinTh * sinTh)));	
}




/**
 *  Convert the pressure into depth with given gravity, 
 * from National Physical Laboratory, Teddington Middleesex, UK TW11 0LW,
 * Underwater Acoustics, Technical Guides
 * @param pres in MPa
 * @param grav value of gravity
 * @return deoth in metres
 */ 
double pressure2depth(double pres , double grav)
{
	return (9.72659E2 - (0.2512 * pres * pres) + (2.279E-4 * pres * pres  *pres) - (1.82E-7 * pres * pres * pres * pres)  )/
			(grav + 1.092E-4 * pres);
	
}


/**
 * return the pressure at depth in MPa (mega pascal)
 * @param depth in metres
 * @param grav value of gravity
 * @return pressure in MPa
 */
double depth2pressure(double depth , double grav)
{
	return ( 0.0100818 * depth + 2.465E-8 * depth * depth - 
				1.25E-13 * depth * depth * depth + 
				2.8E-19 * depth * depth * depth * depth)  * 
				(( grav - 2E-5*depth) / (9.80612 - 2E-5 * depth));
	
}


/**
 * Speed of sound calculation based on Wong and Zhu 1995, temperature range 0 to 40 degC
 * salinity 0 to 40 PSU, pressure 0 - 1000 bar .
 *
 * @param temp temperature
 * @param salt salinity
 * @param depth_pres depth or pressure
 * @param asDepth whether depth_pres is depth
 * @param grav value of gravity
 * @return speed of sound m/sec 
 * 
 */ 
double unesco_speed_of_sound (double temp, double salt, double depth_pres, int asDepth, double grav)
{
	double press = depth_pres;
	double Aw = 0;
	double Aw1 = 0;
	double Aw2 = 0;
	double Aw3 = 0;
	double Bw = 0;
	double Cw = 0;
	double Cw1 = 0;
	double Cw2 = 0;
	double Cw3 = 0;
	double Dw = 0;
	double texp = temp;
	
	if(asDepth)
		press = depth2pressure(depth_pres,grav) * 10;/*UR convert to bar */
		
	texp = temp;
	
	Cw = UNESCO_C00 + UNESCO_C01 * texp;
	Aw = UNESCO_A00 + UNESCO_A01 * texp;
	
	Cw1 = UNESCO_C10 + UNESCO_C11 * texp;
	Aw1 = UNESCO_A10 + UNESCO_A11 * texp;
	
	Cw2 = UNESCO_C20 + UNESCO_C21 * texp;
	Aw2 = UNESCO_A20 + UNESCO_A21 * texp;
	
	Cw3 = UNESCO_C30 + UNESCO_C31 * texp;
	Aw3 = UNESCO_A30 + UNESCO_A31 * texp;
	
	texp = texp * temp;
	Cw +=  UNESCO_C02 * texp;
	Aw +=  UNESCO_A02 * texp;
	
	Cw1 += UNESCO_C12 * texp;
	Aw1 += UNESCO_A12 * texp;
	
	Cw2 += UNESCO_C22 * texp;
	Aw2 += UNESCO_A22 * texp;
	
	Cw3 += UNESCO_C32 * texp;
	Aw3 += UNESCO_A32 * texp;
	
	texp = texp * temp;
	Cw +=  UNESCO_C03 * texp;
	Aw +=  UNESCO_A03 * texp;
	
	Cw1 += UNESCO_C13 * texp;
	Aw1 += UNESCO_A13 * texp;
	
	Cw2 += UNESCO_C23 * texp;
	Aw2 += UNESCO_A23 * texp;
	
	texp = texp * temp;
	Cw +=  UNESCO_C04 * texp;
	Aw +=  UNESCO_A04 * texp;
	
	Cw1 += UNESCO_C14 * texp;
	Aw1 += UNESCO_A14 * texp;
	
	Cw2 += UNESCO_C24 * texp;

	texp = texp * temp;
	Cw +=  UNESCO_C05 * texp;
	
	texp = press;
	
	Cw1 = Cw1 * texp;
	Aw1 = Aw1 * texp;
	
	texp *= press;
	Cw2 = Cw2 * texp;
	Aw2 = Aw2 * texp;
	
	texp *= press;
	Cw3 = Cw3 * texp;
	Aw3 = Aw3 * texp;
	
	Aw += Aw1 + Aw2 +Aw3;
	Cw += Cw1 + Cw2 +Cw3;
	
	Bw = UNESCO_B00 + UNESCO_B01 * temp + (UNESCO_B10 + UNESCO_B11 * temp) * press;
	Dw = UNESCO_D00 + UNESCO_D10 * press;
	
	
	return Cw + Aw * salt + Bw * pow(salt,1.5) + Dw * salt * salt;
		
}

/**
 * Speed of sound calculation based on Mackenzie 1981, temperature range 2 - 30 degC
 * salinity 25 to 45 PSU, depth 0 - 8000 m.
 *
 * @param temp temperature
 * @param salt salinity
 * @param depth_pres depth or pressure
 * @param asDepth whether depth_pres is depth
 * @param grav value of gravity
 * @return speed of sound m/sec 
 * 
 */ 
double mackenzie_speed_of_sound (double temp, double salt, double depth_pres, int asDepth, double grav)
{
	double depth = depth_pres;
	if(asDepth == 0)
		depth = pressure2depth(depth_pres,grav);
		
	salt = salt - 35.;
		
return 1448.96 + 4.591* temp - 5.304E-2*temp*temp +
                 2.374E-4 * temp * temp * temp +
                 1.340 * salt + 
                 1.630E-2 * depth + 
                 1.675E-7 * depth * depth - 
                 1.025E-2 * temp * salt -
                 7.139E-13 * temp * depth * depth * depth;
}


double apel_speed_of_sound (double temp, double salt, double depth_pres, int asDepth, double grav)
{

  double co = 1493.0;
  double ao =  3.0;
  double bo = -0.006;
  double go = -0.04;
  double d0 = 1.2;
  double eo = -0.01;
  double ho = 0.0164;
  double To = 10.0;
  double T1 = 18.0;
  double So = 35.0;
  double sound;
  double depth = depth_pres;

  if(asDepth == 0)
    depth = pressure2depth(depth_pres,grav);

  sound = co + ao * (temp - To) +
    bo * (temp - To) * (temp - To) +
    go * (temp - T1) * (temp - T1) +
    d0 * (salt - So) +
    eo * (temp - T1) * (salt - So) +
    ho * fabs(depth);
  return(sound);
}
