/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: lib/ecology/constants.c
 *  
 *  Description:
 *  Constants used by ecological processe
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: constants.c 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#include "constants.h"

double red_W_C = MW_Carb / MW_Nitr * red_A_C / red_A_N; /* [g(C) g(N)-1] */
double red_W_N = 1.0;           /* [g N g N-1] */
double red_W_P = MW_Phos / MW_Nitr * red_A_P / red_A_N; /* [g(P) g(N)-1] */
double red_W_O = MW_Oxyg / MW_Nitr * red_A_O / red_A_N; /* [g(O2) g(N)-1] */

double atk_W_C = MW_Carb / MW_Nitr * atk_A_C / atk_A_N; /* [g(C) g(N)-1] */
double atk_W_N = 1.0;           /* [g N g N-1] */
double atk_W_P = MW_Phos / MW_Nitr * atk_A_P / atk_A_N; /* [g(P) g(N)-1] */
double atk_W_O = MW_Oxyg / MW_Nitr * atk_A_O / atk_A_N; /* [g(O2) g(N)-1] */

/* Carbon to Oxygen ratio by weight 
 */
double C_O_W = MW_Oxyg / MW_Carb * red_A_O / red_A_C;   /* [g(O2) g(C)-1] */

/* Conversion of milligrams (mg) to moles (mol) 
 */
double mgN2molN = 0.001 / MW_Nitr;      /* mol(N) mg(N)-1 */
double mgP2molP = 0.001 / MW_Phos;      /* mol(P) mg(P)-1 */
