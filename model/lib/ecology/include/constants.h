/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: lib/ecology/include/constants.h
 *  
 *  Description:
 *  Constants used by ecological processes
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: constants.h 5908 2018-08-29 04:27:09Z bai155 $
 *
 */

#if !defined(_CONSTANTS_H)

#define BOLT 1.38066e-23        /* Boltzmann's constant [J K-1] */
#define R 8.31451               /* Universal gas const [J K-1 mol-1] */
#define AV 6.02214e23           /* Avogadro's number [mol-1] */
#define KINEMATIC_VISCOSITY 1.0e-6      /* m2 s-1 */
#define SEC_PER_DAY 86400.0

/* Molecular diffusion coefficients in freshwater T=25 C 
 */
#define DNO3_25_0 19.0e-10      /* [m2 s-1], Li and Gregory, 1974 */
#define DNH4_25_0 19.8e-10      /* [m2 s-1], Li and Gregory, 1974 */
#define DN2_25_0 1.9e-4         /* [m2 s-1] (approx. from Figure in Houghton et al., 1962) */
#define DPO4_25_0 7.34e-10      /* [m2 s-1], Li and Gregory, 1974 */
#define DH4SiO4_25_0 10.73e-10  /* [m2 s-1], Wollast and Garrells, 1971 */
#define DCO2_25_0 19.17e-10     /* [m2 s-1], Wolf-Gladrow and Riebesell,
                                 * 1997 */
#define DFe_25_0 7.19e-10       /* [m2 s-1], Li and Gregory, 1974 */

/* Redfield ratios (Redfield, 1934) for planktonic autotrophs 
 */
#define red_A_I 1060.0          /* [mol(quanta) mol(P)-1] */
#define red_A_C 106.0           /* [mol(C) mol(P)-1] */
#define red_A_N 16.0            /* [mol(N) mol(P)-1] */
#define red_A_P 1.0             /* [mol(P) mol(P)-1] */
#define red_A_O 106.0           /* [mol(O2) mol(P)-1] */

/* Atkinson ratios (Atkinson and Smith, 1983) for benthic autotrophs 
 */
#define atk_A_I 5500.0          /* [mol(quanta) mol(P)-1] */
#define atk_A_C 550.0           /* [mol(C) mol(P)-1] */
#define atk_A_N 30.0            /* [mol(N) mol(P)-1] */
#define atk_A_P 1.0             /* [mol(P) mol(P)-1] */
#define atk_A_O 550.0           /* [mol(O2) mol(P)-1] */

/* Define weight ratios relative to N (because N is a commonly used currency)
 * First need Molecular Weight (MW) of elements - (Atkins, 1994)
 */
#define MW_Carb 12.01           /* [g(C) mol(C)-1] */
#define MW_Nitr 14.01           /* [g(N) mol(N)-1] */
#define MW_Phos 30.97           /* [g(P) mol(P)-1] */
#define MW_Oxyg 32.00           /* [g(O2) mol(O2)-1] */

/*
 * 2 moles of DO lost for every mole of N == 4.57 g (see KWA or JP)
 */
#define NIT_N_0 4.57

extern double red_W_C;          /* [g(C) g(N)-1] */
extern double red_W_N;          /* [g N g N-1] */
extern double red_W_P;          /* [g(P) g(N)-1] */
extern double red_W_O;          /* [g(O2) g(N)-1] */

extern double atk_W_C;          /* [g(C) g(N)-1] */
extern double atk_W_N;          /* [g N g N-1] */
extern double atk_W_P;          /* [g(P) g(N)-1] */
extern double atk_W_O;          /* [g(O2) g(N)-1] */

/* Carbon to Oxygen ratio by weight 
 */
extern double C_O_W;            /* [g(O2) g(C)-1] */

/* Conversion of milligrams (mg) to moles (mol) 
 */
extern double mgN2molN;         /* mol(N) mg(N)-1 */
extern double mgP2molP;         /* mol(P) mg(P)-1 */


#define E2W 217620.0; /* Conversion factor Einstein m-2 s-1 to W m-2 */


#define SWR2PAR 0.43; /* Conversion factor Surface Radiation to PAR */

/* Minimal thicknesses for sediment and water column cells; cells thinner
 * that that are assumed to be empty and skipped.
 */
#define THICKNESS_WC_MIN 0.01
#define THICKNESS_SED_MIN 0.0001

/* KWA #define MASSBALANCE_EPS 1.0e-10 */
#define MASSBALANCE_EPS 1.0e-9
#define TRACERVALUE_EPS  1.0e-12


#define _CONSTANTS_H
#endif
