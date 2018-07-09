/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: lib/ecology/include/ecofunct.h
 *  
 *  Description:
 *  Functions to be used in ecological processes -- header
 *  
 *  Note:          Much of the theory behind these functions can be found at:
 *  http://www.oikos.warwick.ac.uk/ecosystems/ThesisArchive/
 *  baird_thesis.html
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: ecofunct.h 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#if !defined(_ECOFUNCT_H)
#include "constants.h"

double plankgrow(double umax, double aA, double psi, double Sh, double diff[], double m, double I, double conc[], int len, int n, double* grow, double* rate, double Plank_resp);
double benthgrow(double umax, double cf, double vis, double Ub, double ks, double diff[], double m, double I, double conc[], int len, int n, double* grow, double* rate, double Benth_resp);
double PlankMaxNutUp(double psi, double D, double N, double Sh);
double PlankMaxLiteUp(double aA, double I);
double Pgrowth(double umax, double* uptake_rates, double* stoich, double m, int n);
double MolDiff(double D25S0, double T, double S);
double UmaxatT(double umax, double Q10, double T, double Tref);
double Viscosity(double T, double S);
double BenthicMaxNutUpSmooth(double cf, double vis, double D, double N, double Ub);
double BenthicMaxLiteUp(double I, double absorptance);
double phi(char* F, double rP, double Pswim, double Psink, double rH, double Hswim, double Hsink, double epsilon, double viscosity, double density, double Twater);
double NaturalMortality(double Pop, double rate);
double BenthicMaxNutUpRough(double cf, double vis, double D, double N, double Ub, double ks);
double Plookup(double umax, double* uptake_rates, double* stoich, double m, int n);
double CRgrowth(double umax, double uptake_rates[], double m, double* normgrow, double* rate, int klen, int n, double resp);
double MPgrowth(double umax, double uptake_rates[], int n);
double LMgrowth(double umax, double uptake_rates[], double stoich[], int n);
double RectHyp(double maximum, double halfsat, double conc);
double LDtoRD(double decay, double frac, int n);
int extract_grow(int n, char* file_name, double** ptable, double** prate);
double aa(double rad, double absorb);
double psi(double rad);
double PhyCellMass(double r);
double ZooCellMass(double r);
double PhyCellChl(double r);
double diffaa(double yC,double r);
void absorbwave(double *yC, double chla, int num_wave, double *wave);
void absorbxanth(double *yC, double xanth, int num_wave, double *wave);
void absorbxanth_heat(double *yC, double xanth, int num_wave, double *wave);
void aawave(double rad, double *absorb, double *aA, int num_wave, double *wave);

#define _ECOFUNCT_H
#endif
