/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/underwater.h
 *
 *  \brief Include file for speed of sound routines
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: underwater.h 5834 2018-06-27 00:55:37Z riz008 $
 */

#ifndef _UNDERWATER_H_
#define _UNDERWATER_H_




typedef double (*ssfn) (double temp, double salt, double depth_pres, int asDepth, double grav);

double gravity(double lat);

double pressure2depth(double pres , double grav);
double depth2pressure(double depth , double grav);




/*
 
 Constants for UNESCO SPEED of SOUND calculation
  
#ifdef SS_UNESCO  
 */

#define UNESCO_A00 1.389
#define UNESCO_A01 -1.262E-2
#define UNESCO_A02 7.166E-5
#define UNESCO_A03 2.008E-6
#define UNESCO_A04 -3.21E-8
#define UNESCO_A10 9.4742E-5
#define UNESCO_A11 -1.2583E-5
#define UNESCO_A12 -6.4928E-8
#define UNESCO_A13 1.0515E-8
#define UNESCO_A14 -2.0142E-10
#define UNESCO_A20 -3.9064E-7
#define UNESCO_A21 9.1061E-9
#define UNESCO_A22 -1.6009E-10
#define UNESCO_A23 7.994E-12
#define UNESCO_A30 1.100E-10
#define UNESCO_A31 6.651E-12
#define UNESCO_A32 -3.391E-13

#define UNESCO_B00 -1.922E-2
#define UNESCO_B01 -4.42E-5
#define UNESCO_B10 7.3637E-5
#define UNESCO_B11 1.7950E-7

#define UNESCO_C00 1402.388
#define UNESCO_C01 5.0383
#define UNESCO_C02 -5.81090E-2
#define UNESCO_C03 3.3432E-4
#define UNESCO_C04 -1.47797E-6
#define UNESCO_C05 3.14419E-9
#define UNESCO_C10 0.153563
#define UNESCO_C11 6.8999E-4
#define UNESCO_C12 -8.1829E-6
#define UNESCO_C13 1.3632E-7
#define UNESCO_C14 -6.1260E-10
#define UNESCO_C20 3.1260E-5
#define UNESCO_C21 -1.7111E-6
#define UNESCO_C22 2.5986E-8
#define UNESCO_C23 -2.5353E-10
#define UNESCO_C24 1.0415E-12
#define UNESCO_C30 -9.7729E-9
#define UNESCO_C31 3.8513E-10
#define UNESCO_C32 -2.3654E-12

#define UNESCO_D00 1.727E-3
#define UNESCO_D10 -7.9836E-6



double unesco_speed_of_sound (double temp, double salt, double depth_pres, int asDepth, double grav);


/*
 * 
#endif


Constants for the Mackenzie SPEED of SOUND calculation 

#iddef SS_MACKANZIE
 * 
 */

double mackenzie_speed_of_sound (double temp, double salt, double depth_pres, int asDepth, double grav);



/*
 
#endif
 */
 

/*
 * 
#endif


Constants for the Apel SPEED of SOUND calculation 
(Apel, 1987, Principles of Ocean Physics, p 349)

#iddef SS_APEL
 * 
 */

double apel_speed_of_sound (double temp, double salt, double depth_pres, int asDepth, double grav);

/*
 
#endif
 */

extern ssfn  speed_of_sound;
 
 
#endif /*UNDERWATER_H_*/
