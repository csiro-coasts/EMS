/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/carbon_chemistry_wc.c
 *  
 *  Description:
 *  
 *  Equilibrium state calculation of carbon chemistry system. This process can be applied in 
 *  either the water column or sediment (where it calculates porewater chemistry).
 *
 *  The numberical scheme is adapted from OCMIP routines, but includes a larger search range for pH for the iterative scheme.
 *
 *  For desription see:
 *  Mongin, M., M. E. Baird, B. Tilbrook, R. J. Matear, A. Lenton, M. Herzfeld, K. A. Wild-Allen, J. Skerratt, N. Margvelashvili, 
 *  B. J. Robson, C. M. Duarte, M. S. M. Gustafsson, P. J. Ralph, A. D. L. Steven (2016). The exposure of the Great Barrier Reef 
 *  to ocean acidification. Nature Communications 7, 10732.
 * 
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: carbon_chemistry_wc.c 6542 2020-05-06 06:27:36Z bai155 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ecology_internal.h"
#include "utils.h"
#include "ecofunct.h"
#include "constants.h"
#include "eprocess.h"
#include "stringtable.h"
#include "cell.h"
#include "column.h"
#include "einterface.h"
#include "carbon_chemistry_wc.h"

    void get_co2_param (double t,double s,double *ff,double *k0,double *k1,
						double *k2,double *kb,double*kp1,double *kp2,double *kp3,
						double *ksi,double *kw,double *ks,double *kf,double *bt,
						double *st,double *ft, double *calcium, double *karag,double *kcal)	;

    void ocmip2_co2calc(double  dic_in,double ta_in,double pt_in,double sit_in,
                        double *htotal,double xco2_in,double atmpres,
						double *co2star,double *co2starair,double *dco2star,
					    double *pco2surf,double *dpCO2,double ff,double k0,
					    double k1,double k2,double kb,double kp1,double kp2,
					    double kp3,double ksi,double kw,double ks,double kf,
					    double bt,double st,double ft);
     
    void  drtsafe(double *h_total3,double x1,double x2, double xacc,
						double ff,double k0,double k1,double k2,double kb,double kp1,
						double kp2,double kp3,double ksi,double kw,
						double ks,double kf,double bt,double st,double ft,
						double dic, double pt, double sit, double ta) ;


    void ta_inter1(double x,double *fn,double *df,double ff,double k0,
						double k1,double k2,double kb,double k1p,double k2p,
						double k3p,double ksi,double kw,double ks,double kf,
						double bt,double st,double ft,double dic,double pt, 
						double sit, double ta);
typedef struct {

  int do_mb;                
  
  /*parameters*/
  
  //  double xco2_air;    /*  atmospheric mole fraction CO2 in dry air (ppmv)*/
  
  /* Tracers */
  
  int CO32_i;   	  /* carbonate */
  int HCO3_i;   	  /*bi carbonate*/
  int PH_i;      	  /* PH*/
  int CO2_starair_i;/* co2starair =  xco2 * ff * atmpres~  
		       xco2 is the  atmospheric mole fraction CO2 in dry air (mol/m^3)*/
  /* ff is the CFC solubility for water-vapor saturated air [mol/(mˆ3 * picoatm)]
     atmpres is the atmospheric pressure*/
  int CO2_star_i	; /* co2star    = ocean surface aquaeous CO2 concentration  (mol/m^3)*/
  int TEMP_i;	      /*temperature in degree c*/
  int SALT_i;	      /* salinity PSU*/
  int DIC_i;	      /*    INPUT total inorganic carbon (mmol/m^3)*/
  int ALK_i;	      /*  total alkalinity (eq/m^3) */  
  int dco2star_i;   /* dco2star   = Delta CO2* (mol/m^3)*/
  int pco2surf_i;	  /* pco2surf   = oceanic pCO2 (ppmv)*/
  int dpCO2_i;	  /* dpco2      = Delta pCO2, i.e, pCO2ocn - pCO2atm (ppmv)*/
  int omega_ca_i;   /* calcite saturation state*/ 
  int omega_ar_i;   /* aragonite saturation state*/ 
  int xco2_air_i;	      /*    INPUT pco2 in the atmosphere*/ 
  double xco2_air;	      /*    INPUT pco2 if constant value in atmosphere*/

} workspace;


/* This is only called once during the lifetime of a process/ecology. */
/*******************************************************************************/
void carbon_chemistry_wc_init(eprocess* p)
/*******************************************************************************/

{
  ecology* e = p->ecology;
  stringtable* tracers = e->tracers;
  workspace* ws = malloc(sizeof(workspace));
  p->workspace = ws;
  
  /*parameters*/
  
  ws->xco2_air = try_parameter_value(e, "xco2_in_air");   

  if (isnan(ws->xco2_air)){
    eco_write_setup(e,"No parameter: xco2_in_air, so assuming a tracer instead \n");
    ws->xco2_air_i = e->find_index(tracers,"xco2_in_air", e);
  }
   
  /* TRACERS - necessary ones */
  
  ws->TEMP_i = e->find_index(tracers, "temp", e);
  ws->SALT_i = e->find_index(tracers, "salt", e);
  ws->DIC_i = e->find_index(tracers,"DIC", e);
  ws->ALK_i = e->find_index(tracers,"alk", e);
  ws->PH_i = e->find_index(tracers,"PH", e);
 
  /* dco2star is diagnostic, but necessary since comes from ocmip2_co2calc, and must 
     be found by gas_exchange_wc.c */

  ws->dco2star_i = e->find_index(tracers,"dco2star", e);
  
  /* TRACERS - diagnostic ones */

  ws->CO2_starair_i = e->try_index(tracers,"CO2_starair", e);
  ws->CO2_star_i = e->try_index(tracers,"CO2_star", e);
  ws->CO32_i = e->try_index(tracers,"CO32", e);
  ws->HCO3_i = e->try_index(tracers,"HCO3", e);
  ws->pco2surf_i = e->try_index(tracers,"pco2surf", e);
  ws->dpCO2_i = e->try_index(tracers,"dpCO2", e);
  ws->omega_ar_i = e->try_index(tracers,"omega_ar", e);
  ws->omega_ca_i = e->try_index(tracers,"omega_ca", e);

}
/*******************************************************************************/
void carbon_chemistry_wc_postinit(eprocess* p)
/*******************************************************************************/

{
  ecology* e = p->ecology;
  workspace* ws = p->workspace;
   
}
/*******************************************************************************/
void carbon_chemistry_wc_destroy(eprocess* p)
/*******************************************************************************/
{
	free(p->workspace);

}
/*******************************************************************************/
void carbon_chemistry_wc_precalc(eprocess* p, void* pp)
/*******************************************************************************/
{
  workspace* ws = p->workspace;
  cell* c = (cell*) pp;
  double* y = c->y;
  double* cv = c->cv;
  double dz_wc = c->dz_wc;

  /*LOCAL DECLARATION*/
  
  double co32 ; 				/*diagnostic */
  double hco3 ; 				/*diagnostic*/
  double PH =y[ws->PH_i];  	/*diagnostic but need to keep the previous time step value to be modified*/
  double co2star; 			/*diagnostic*/
  double co2starair; 			/*diagnostic*/
  double temp = y[ws->TEMP_i];
  double salt = max(20.0,y[ws->SALT_i]);
  double dic_in = y[ws->DIC_i] / 12.01 ;  // mg to mmol
  double ta_in=y[ws->ALK_i];

  double xco2_air_in;

  if (isnan(ws->xco2_air))
	    xco2_air_in=y[ws->xco2_air_i];
  else
    xco2_air_in = ws->xco2_air;


  double dco2star; 			/*diagnostic*/
  double pco2surf;			/* diagnostic*/
  double dpCO2 ;				/*diagnostic*/
  double h_total;				/*diagnostic*/
  double omega_ca;			/*diagnostic*/
  double omega_ar;			/*diagnostic*/
  double omegaca;				/*diagnostic*/
  double omegaar;				/*diagnostic*/
  double calcium;				/*internal parameter*/
  double karag;       		/*internal parameter*/
  double kcal;        		/*internal parameter*/
  
  /* define internal parameters */    
  double pt_in=0.00;    /*  inorganic phosphate (mol/m^3) */
  double sit_in=0.00;   /*  inorganic silicate (mol/m^3) */
  double atmpres=1.0;      /*  atmospheric pressure in atmospheres (1 atm = 1013.25 mbar)*/
  
  
  /* define for co2 air sea fluxes - moved to gas exchange process 

  double scco2;
  double io;      
  double Smallk;
  double U;                                      */
  
  /* declare carbon chemistry constant parameters  as double*/
  double ff_co2,k0_co2,k1_co2,k2_co2,kb_co2,kp1_co2,kp2_co2,kp3_co2,ksi_co2,kw_co2,ks_co2,kf_co2,bt_co2,st_co2,ft_co2;
    
  dic_in =dic_in/1000.0 ;/* from mmol m-3 to mol m-3*/
  
  ta_in=ta_in/1000;
  
  
  /* define H+ concentration */

   if (PH < 4.5)
     PH = 7.0;
   if (PH > 10.5)
    PH = 7.0;

  h_total=pow(10,-PH);
  
  /*********************************************************************/	
  /*call routine that calculate the parameter call then with &ff which point to the adress of *ff any change inthe routine will change the adress of ff hence the value of the double ff.
  /*-----------------------------------------------------------------    */
  get_co2_param(temp,salt,&ff_co2,&k0_co2,&k1_co2,&k2_co2,&kb_co2,&kp1_co2,&kp2_co2,&kp3_co2,&ksi_co2,&kw_co2,&ks_co2,&kf_co2,&bt_co2,&st_co2,&ft_co2,&calcium,&karag,&kcal);

  /*********************************************************************/	
  /*Calculate delta co2* and PH from total alkalinity and total CO2 at 
    temperature (t), salinity (s) and "atmpres" atmosphere total pressure.  */
  /*********************************************************************/	
  /*---------------------------------------------------------------------*/
  ocmip2_co2calc(dic_in,ta_in,pt_in,sit_in,&h_total,xco2_air_in,atmpres,&co2star,&co2starair,&dco2star,&pco2surf,&dpCO2,ff_co2,k0_co2,k1_co2,k2_co2,kb_co2,kp1_co2,kp2_co2,kp3_co2,ksi_co2,kw_co2,ks_co2,kf_co2,bt_co2,st_co2,ft_co2);

  /*********************************************************************/	
  /*Calculate carbonate speciation from  */
  /********************************************************************/	    
  
  hco3=dic_in*k1_co2*h_total/(h_total*h_total+k1_co2*h_total+k1_co2*k2_co2);
  co32=dic_in*k1_co2*k2_co2/(h_total*h_total+k1_co2*h_total+k1_co2*k2_co2);
  
  /*********************************************************************/	
  /*Calculate omega calcite and aragonite saturation */
  /********************************************************************/	    
  omegaca = calcium*co32/kcal;    /*co32 is in mol m-3 so does calcium */
  omegaar = calcium*co32/karag;   /* idem*/
  
  /*********************************************************************/	
  /*update the variables if they are in the .prm file */
  /********************************************************************/
  
  if (ws-> CO2_star_i > -1)
    y[ws-> CO2_star_i]    = co2star;    /* co2 [] in water mol/m-3*/
  if (ws-> CO2_starair_i > -1)
    y[ws-> CO2_starair_i] = co2starair; /* co2 [] in air im mol m-3*/
  y[ws-> dco2star_i]    = dco2star;   /* delta CO2 mol m-3*/
  if (ws-> pco2surf_i > -1)  
    y[ws-> pco2surf_i]    = pco2surf;   /* oceanic pCO2 in ppmv*/
  if (ws-> dpCO2_i > -1)
    y[ws-> dpCO2_i]       = dpCO2;      /* delta pco2 pco2 ocean - po2atm in ppm*/
  if (ws-> CO32_i > -1)
    y[ws-> CO32_i]        = co32*1000.0; /*from mol to mmol m-3*/      
  if (ws-> HCO3_i > -1)    
    y[ws-> HCO3_i]        = hco3*1000.0;  /*from mol to mmol m-3*/
  if (ws-> omega_ca_i > -1)
    y[ws-> omega_ca_i]    = omegaca;
  if (ws-> omega_ar_i > -1){
    if (omegaar > 0.0){ 
      y[ws-> omega_ar_i]  = omegaar;
    }else{
      y[ws->omega_ar_i] = 3.8;
    }
  }
 

  if (-log10(h_total)<0){
    y[ws-> PH_i] =y[ws->PH_i];}
  else{
    y[ws-> PH_i]          = -log10(h_total);
  }  
}

/*******************************************************************************/
void carbon_chemistry_wc_calc(eprocess* p, void* pp)
/*******************************************************************************/
{ 

  workspace* ws = p->workspace;
  intargs* ia = (intargs*) pp;
  cell* c = ((cell*) ia->media);
  double* cv = c->cv;
  double* y = ia->y;
  double* y1 = ia->y1;
  double dz_wc = c->dz_wc;
  
}

void carbon_chemistry_wc_postcalc(eprocess* p, void* pp)
{
  cell* c = ((cell*) pp);
  workspace* ws = p->workspace;
  double* y = c->y;
  
}

/*---------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------*/
void get_co2_param(double t,double s,double *ff_co2,double *k0_co2,double *k1_co2,double *k2_co2,double *kb_co2,double *kp1_co2,double *kp2_co2,double *kp3_co2, double *ksi_co2,double *kw_co2,double *ks_co2,double *kf_co2,double *bt_co2,double *st_co2,double *ft_co2,double *calcium,double *karag,double *kcal)	
/*---------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------*/
/* the routine change the value of *ff=7   which copy the value 7 to the adress pointed by ff  /
   {
   /**/
{
  double tk;
  double tk100;
  double tk1002;
  double invtk;
  double dlogtk;
  double is;
  double is2;
  double sqrtis;
  double s2;
  double sqrts;
  double s15;
  double scl;
  double logf_of_s;
  double lnksp0cal;
  double lnksp0arag;
  /************************************************************************/     
  /*some pre calculation first: on salinity and temperature */
  /************************************************************************/  
  tk        = 273.15 + t   ;             /*temperature in Kelvins*/
  tk100     = tk / 100.0;
  tk1002    = tk100 * tk100;
  invtk     = 1.0 / tk;
  dlogtk    = log(tk);
  
  /*salinity*/
  is        = 19.924 * s /(1000.0 - 1.005 * s);
  is2       = is * is;
  sqrtis   = sqrt(is);
  s2        = s * s;
  sqrts     = sqrt(s);
  s15       = pow(sqrts, 3);
  scl       = s / 1.80655;
  logf_of_s = log(1.0 - 0.001005 * s)    ;
  
  /************************************************************************/     
  /* then some calculation of some coefficient
   ************************************************************************/      
  /*================================================
    f = k0(1-pH2O)*correction term for non-ideality
    Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)
    ==================================================*/
  
  *ff_co2 = exp(-162.8301 + 218.2968 / tk100  +  90.9241 *(dlogtk - log(100)) - 1.47696 * tk1002 +               
		s *(0.025695 - 0.025225 * tk100 + 0.0049867 * tk1002));
  
  
  
  /*======================================================
    k0 from Weiss 1974
    =====================================================*/
  
  *k0_co2 = exp(93.4517/tk100 - 60.2409 + 23.3585 * log(tk100) +           
		s *(0.023517 - 0.023656 * tk100 +0.0047036 * tk1002));
  
  
  /*/!=========================================================
    ! k1 = [H][HCO3]/[H2CO3]
    ! k2 = [H][CO3]/[HCO3]
    !
    ! Millero p.664 (1995) using Mehrbach et al. data on seawater scale 
    !==================================================================*/
  
  *k1_co2 =pow( 10.0 , (-(3670.7 * invtk -62.008 + 9.7944 * dlogtk - 0.0118 * s +0.000116 * s2)));
  *k2_co2 = pow(10.0 , (-(1394.7 * invtk + 4.777 - 0.0184 * s +0.000118 * s2)));
  
  /*!===========================================================
    ! q = [H][BO2]/[HBO2]
    !===========================================================
    ! Millero p.669 (1995) using data from Dickson (1990)
    !*/
  
  *kb_co2 = exp((-8966.90 - 2890.53 * sqrts - 77.942 * s +1.728 * s15 -0.0996 * s2) * invtk +   
		(148.0248 + 137.1942 * sqrts +1.62142 * s) +(-24.4344 - 25.085 * sqrts -     
							     0.2474 * s) * dlogtk + 0.053105 * sqrts * tk);
  
  /*!================================================================
    ! k1p = [H][H2PO4]/[H3PO4]
    !================================================================
    ! DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
    !*/
  
  *kp1_co2 = exp(-4576.752 * invtk + 115.525 - 18.453 * dlogtk +               
		 (-106.736 * invtk + 0.69171) *sqrts +                         
		 (-0.65643 * invtk - 0.01844) * s);
  
  /*!
    ! k2p = [H][HPO4]/[H2PO4]
    !
    ! DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
    !*/
  
  *kp2_co2= exp(-8814.715 * invtk + 172.0883 - 27.927 * dlogtk +               
		(-160.340 * invtk + 1.3566) *sqrts +                         
		(0.37335 * invtk- 0.05778) * s);
  
  /*!
    !-----------------------------------------------------------------------
    ! k3p = [H][PO4]/[HPO4]
    !
    ! DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
    !*/
  
  *kp3_co2 = exp(-3070.75 * invtk - 18.141 +(17.27039 * invtk + 2.81197) *  
		 sqrts + (-44.99486 * invtk - 0.09984) * s);
  
  /*!
    !-----------------------------------------------------------------------
    ! ksi = [H][SiO(OH)3]/[Si(OH)4]
    !
    ! Millero p.671 (1995) using data from Yao and Millero (1995)
    !*/
  
  *ksi_co2 = exp(-8904.2 * invtk + 117.385 - 19.334 * dlogtk +               
		 (-458.79 * invtk+ 3.5913) *sqrtis +                        
		 (188.74 * invtk- 1.5998) *is +                            
		 (-12.1652 * invtk + 0.07871) *  
		 is2 + logf_of_s);
  
  /*!
    !-----------------------------------------------------------------------
    ! kw = [H][OH]
    !
    ! Millero p.670 (1995) using composite data
    !*/
  
  *kw_co2 = exp(-13847.26 * invtk + 148.9652 -23.6521 * dlogtk +               
		(118.67 * invtk - 5.977 +1.0495 * dlogtk) *              
		sqrts - 0.01615 * s);
  
  /*!
    !-----------------------------------------------------------------------
    ! ks = [H][SO4]/[HSO4]
    !
    ! Dickson (1990, J. chem. Thermodynamics 22, 113)
    !*/
  
  *ks_co2 = exp(-4276.1 * invtk + 141.328 -      
		23.093 * dlogtk +                
		(-13856.0 * invtk + 324.57 -     
		 47.986 * dlogtk) *              
		sqrtis +                         
		(35474.0 * invtk - 771.54 +      
		 114.723 * dlogtk) * is - 
		2698.0 * invtk *                 
		pow(sqrtis, 3 )+                    
		1776.0 * invtk * is2 +    
		logf_of_s);
  
  /*!
    !-----------------------------------------------------------------------
    ! kf = [H][F]/[HF]
    !
    ! Dickson and Riley (1979) -- change pH scale to total
    !*/
  
  *kf_co2 = exp(1590.2 * invtk - 12.641 +        
		1.525 * sqrtis +                 
		logf_of_s +                      
		log(1.0 + (0.1400 / 96.062) *           
		    scl / *ks_co2));
  
  
  
  /*C kcal and karagonite
    c ksp0cal = gamma*[Ca]gamma*[CO3]/[CaCO3 calcite]
    c  ksp0arag = gamma*[Ca]gamma*[CO3]/[CaCO3 aragonite]
    c  Mucci (1983)*/
  
  /*c  kcal = [Ca][CO3]/[CaCO3 calcite]*/
  
  lnksp0cal=(-171.9065-0.077993*tk+2839.319/tk+71.595*log10(tk)+(-0.77712+0.0028426*tk+178.34/tk)*pow(s,0.5)-0.07711*s+0.0041249*pow(s,1.5));
  
  lnksp0arag=(-171.945-0.077993*tk+2903.293/tk+71.595*log10(tk)+(-0.068393+0.0017276*tk+88.135/tk)*pow(s,0.5)-0.10018*s+0.0059415*pow(s,1.5));
  
  *kcal=pow(10,lnksp0cal);
  *karag=pow(10,lnksp0arag);
  /*printf("muchi****************************\n");
    printf("kcal is %lf\n",*kcal);
    printf("karag is %lf\n",*karag);
    printf("-log(kcal) is %lf\n",-log10(*kcal));
    printf("-log(karag) is %lf\n",-log10(*karag));
    printf("-log(kcal)2 is %lf\n",-lnksp0cal);
    printf("-log(karag) 2is %lf\n",-lnksp0arag);
    printf("****************************************");*/
  
  /*lnksp0cal = -395.8293 + 6537.773/tk + 71.595 *log(tk) - .17959 * tk;
    lnksp0arag = -395.9180 + 6685.079/tk + 71.595 *log(tk) - .17959 * tk;
    
    *kcal = exp(lnksp0cal +(-1.78938 + 410.64/tk + .0065453 * tk) * sqrt(s) -.17755 * s + .0094979 * pow(s,1.5));
    
    *karag = exp(lnksp0arag + (-.157481 + 202.938/tk + .0039780 * tk) * sqrt(s) -.23067 * s + .0136808 * pow(s,1.5));
    
    */
  
  /*printf("ric ha****************************\n");
    
    printf("kcal is %lf\n",*kcal);
    printf("karag is %lf\n",*karag);
    printf("-log(kcal) is %lf\n",-log10(*kcal));
    printf("-log(karag) is %lf\n",-log10(*karag));*/
  /*!
    !-----------------------------------------------------------------------
    ! Calculate concentrations for borate, sulfate, and fluoride
    !
    ! Uppstrom (1974)
    !*/
  
  *bt_co2= 0.000232 / 10.811 * scl;
  /*
    !
    ! Morris & Riley (1966)
    !*/
  
  *st_co2 = 0.14 / 96.062 * scl;
  
  /*!
    ! Riley (1965)
    !*/
  
  *ft_co2 = 0.000067 / 18.9984 * scl;
  
  /* Millero (1982) mol m-3*/
  *calcium = (.010280*s/35.)/1000.   ;   /*mol per l ~ mmol m-3 to mol m-3 */       
  
}

/*---------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------*/
void ocmip2_co2calc(double  dic_in,double ta_in,double pt_in,double sit_in,
                    double *h_total,double xco2_in,double atmpres,
                    double *co2star,double *co2starair,double *dco2star,double *pco2surf,double *dpCO2,
		    double ff,double k0,double k1,double k2,double kb,double kp1,double kp2,double kp3,double ksi,
                    double kw,double ks,double kf,double bt,double st,double ft)
/*---------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------*
  
  Calculate delta co2* from total alkalinity and total CO2 at
  temperature (t), salinity (s) and "atmpres" atmosphere total pressure.
  It is assumed that init_ocmip2_co2calc has already been called with
  the T and S to calculate the various coefficients.
  
  INPUT
  dic_in     = total inorganic carbon (mol/m^3)  where 1 T = 1 metric ton = 1000 kg
  ta_in      = total alkalinity (eq/m^3) 
  pt_in      = inorganic phosphate (mol/m^3) 
  sit_in     = inorganic silicate (mol/m^3) 
  htotallo   = lower limit of htotal range
  htotalhi   = upper limit of htotal range
  htotal     = H+ concentraion
  xco2_in    = atmospheric mole fraction CO2 in dry air (ppmv) 
  atmpres    = atmospheric pressure in atmospheres (1 atm = 1013.25 mbar)
  carbon chemistry constant= ff,k0,k1,k2,kb,kp1,kp2,kp3,ksi,kw,ks,kf,bt,st,ft
  
  OUTPUT
  co2star    = CO2*water (mol/m^3)
  co2starair = CO2*air (mol/m^3)
  dco2star   = Delta CO2* (mol/m^3)
  pco2surf   = oceanic pCO2 (ppmv)
  dpco2      = Delta pCO2, i.e, pCO2ocn - pCO2atm (ppmv)
  
  IMPORTANT: Some words about units - (JCO, 4/4/1999)
  
  - Models carry tracers in mol/m^3 (on a per volume basis)
  
  - Conversely, this routine, which was written by observationalists
  (C. Sabine and R. Key), passes input arguments in umol/kg  
  (i.e., on a per mass basis)
  
  - I have changed things slightly so that input arguments are
  in mol/m^3,
  
  - Thus, all input concentrations (dic, ta, pt, and st) should be 
  given in mol/m^3; output arguments "co2star" and "co2starair"  
  and "dco2star" are likewise be in mol/m^3.
  
  FILES and PROGRAMS NEEDED: drtsafe, ta_iter_1
*/	
{
  /*declare conversion factor */
  double permil = 1.0 / 1024.5;
  double permeg = 1.e-6;
  double xacc = 1.0e-12;
  
  
  /* declare locale variables*/
  double sit;
  double ta;
  double dic;
  double pt;
  double xco2;
  /*-----------------------------------------------------------------------
    !       Change units from the input of mol/m^3 -> mol/kg:
    !       (1 mol/m^3)  x (1 m^3/1024.5 kg)
    !       where the ocean's mean surface density is 1024.5 kg/m^3
    !       Note: mol/kg are actually what the body of this routine uses 
    !       for calculations.  Units are reconverted back to mol/m^3 at the 
    !       end of this routine.
    !---------------------------------------------------------------------
    !
    !       To convert input in mol/m^3 -> mol/kg 
  */
  
  sit = sit_in * permil;  /*inorganic silicate  */
  ta  = ta_in  * permil;  /*total alkalinity  */
  dic = dic_in * permil;  /*total inorganic carbon  */
  pt  = pt_in  * permil;  /*inorganic phosphate  */
  
  
  /* !
     !---------------------------------------------------------------------
     !       Change units from uatm to atm. That is, atm is what the body of 
     !       this routine uses for calculations.
     !       Note: units are reconverted bac to uatm at END of this routine.
     !---------------------------------------------------------------------
     !
     !       To convert input in uatm -> atm
     !*/
  
  xco2= xco2_in * permeg;
  
  
  /*!
    /*-----------------------------------------------------------------------
    !
    ! Calculate [H+] total when DIC and TA are known at T, S and 1 atm.
    ! The solution converges to err of xacc. The solution must be within
    ! the range x1 to x2.
    !
    ! If DIC and TA are known then either a root finding or iterative method
    ! must be used to calculate htotal. In this case we use the
    ! Newton-Raphson "safe" method taken from "Numerical Recipes"
    ! (function "rtsafe.f" with error trapping removed).
    !
    ! As currently set, this procedure iterates about 12 times. The x1
    ! and x2 values set below will accomodate ANY oceanographic values.
    ! If an initial guess of the pH is known, then the number of
    ! iterations can be reduced to about 5 by narrowing the gap between
    ! x1 and x2. It is recommended that the first few time steps be run
    ! with x1 and x2 set as below. After that, set x1 and x2 to the
    ! previous value of the pH +/- ~0.5. The current setting of xacc will
    ! result in co2star accurate to 3 significant figures (xx.y). Making
    ! xacc bigger will result in faster convergence also, but this is not
    ! recommended (xacc of 10**-9 drops precision to 2 significant
    ! figures).
    !*/
  
  double htotalhi=pow(10,(-3.5+log10( *h_total)));
  double htotallo=pow(10,(+3.5+log10( *h_total)));
  /*double htotalhi=pow(10,-10);
    double htotallo=pow(10,-5);*/
  
     if (+3.5+log10( *h_total)<1) {
       htotallo=pow(10,-2);
       }
       if (-3.5+log10( *h_total)>14) {
       htotallo=pow(10,-12);
       }
  
  /* printf("ph hi is %lf\n", -log10(htotalhi));
     printf("ph low  %lf\n",  -log10(htotallo));*/
  
  
  
  double htotal2;
  double h_total3;
  
  /**/
  drtsafe(&h_total3,htotalhi, htotallo, xacc,ff,k0,k1,k2,kb,kp1,kp2,kp3,ksi,kw,ks,kf,bt,st,ft,dic,pt,sit,ta);  
  /* calculate H+ (h_total3) within the range htotalhi
     and htotallo with the error xacc  */
  
  *h_total=h_total3;
  /*printf("h_total after drtsafe is %lf\n",h_total3*1000000);*/
  /*!
    ! Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2, 
    ! ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
    !*/
  
  htotal2    = *h_total * *h_total;
  *co2star    = dic* htotal2 /(htotal2 +k1 * *h_total +k1* k2);
  
  *co2starair = xco2 * ff * atmpres;
  *dco2star   = *co2starair - *co2star;
  
  
  /*!
    !---------------------------------------------------------------
    !c      Add two output arguments for storing pCO2surf
    !c      Should we be using K0 or ff for the solubility here?
    !---------------------------------------------------------------
    !*/
  
  *pco2surf =*co2star / ff;
  *dpCO2   = *pco2surf - xco2 * atmpres;
  
  /*!
    !----------------------------------------------------------------
    !
    ! Convert units of output arguments
    !      Note: dco2star, co2star and co2starair are calculated in
    !            mol/kg within this routine 
    !      Thus Convert now from mol/kg -> mol/m^3
    !*/
  
  
  *co2star= *co2star / permil;
  
  *dco2star  = *dco2star / permil;
  *co2starair = *co2starair/ permil;
  
  
  /*!
    !      Note: pCO2surf and dpCO2 are calculated in atm above. 
    !      Thus convert now to uatm
    !*/
  
  *pco2surf = *pco2surf/ permeg;
  *dpCO2    = *dpCO2 / permeg;
  
}




/*---------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------*/
void   drtsafe(double *h_total3,double x1,double x2, double xacc,
	       double ff,double k0,double k1,double k2,double kb,double kp1,double kp2,
	       double kp3,double ksi,double kw,
	       double ks,double kf,double bt,double st,double ft,double dic, double pt, double sit, double ta)    
/*---------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*needs k1 k2 k1p k2p  k3p st ks dic sit ksi */ 
/*---------------------------------------------------------------------------*/   
/* function drtsafe:
   INPUT
   x1=lower initial guess for H+ concentration
   x2=upper initial guess for H+ concentration
   xacc= error
   
   OUTPUT
   h_total=H+ concentration
*/
{   
  
  /*define local variables */
  int maxit=100;
  int j,num;
  double fh,swap,xl,xh,dxold,dxx,f,tempp;
  double fl;
  double df;
  double h_total4;
  double dum;
  double dum2;
  /*call the iteration function */
  /*input is x1 (lower range of total H+), output fl which is the calculated value for TA and "df" is the value for dTA/dhtotal*/
  ta_inter1(x1,&fl,&df,ff,k0,k1,k2,kb,kp1,kp2,kp3,ksi,kw,ks,kf,bt,st,ft,dic,pt,sit,ta);
  /* printf("x1 is %lf\n",x1);
     printf("df is %lf\n",df);
     printf("fl is %lf\n",fl);*/
  
  /*input is x2 (upper range of total H+)output fh which is the calculated value for TA and "df" 
    void  ta_inter1(double x,double *fn,double *df
    ,double ff,double k0,double k1,double k2,double kb,double k1p,double k2p
    ,double k3p,double ksi,double kw,double ks,double kf,double bt,double st,double ft,double dic,double pt, double sit, double ta)
    is the value for dTA/dhtotal*/
  ta_inter1(x2,&fh,&df,ff,k0,k1,k2,kb,kp1,kp2,kp3,ksi,kw,ks,kf,bt,st,ft,dic,pt,sit,ta);
  
  /* printf("x2 is %lf\n",x2);
     printf("df is %lf\n",df);
     printf("fh is %lf\n",fh);*/
  
  
  if(fl <0.0)   {
    xl=x1;
    xh=x2;    }
  else          {
    xh=x1;
    xl=x2;
    swap=fl;
    fl=fh;
    fh=swap;
  }
  
  *h_total3=0.5*(x1+x2);
  dum2=-log10(*h_total3);
  
  /* printf("h_total at stsart of drtsafe  %lf\n",*h_total3*1000);
     printf("ph at stsart of drtsafe  %lf\n",dum2);*/
  
  dxold=(x2-x1);
  dxx=dxold;
  
  
  
  /*input is drtsafe (middle range pf total H+)output f which is the calculated value for TA and "df" is the value for dTA/dhtotal*/   
  ta_inter1(*h_total3,&f,&df,ff,k0,k1,k2,kb,kp1,kp2,kp3,ksi,kw,ks,kf,bt,st,ft,dic,pt,sit,ta);
  
  
  for(j=1;j<maxit;j++) 
    {   
      
      /*printf(" iteration %d \n ",j);
	printf("dx is %lf \n",dxx);
	printf("*h_total3 is %lf \n",*h_total3);
	printf("xh is %lf \n",xh);
	printf("x1 is %lf \n",x1);
	printf("x2 is %lf \n",x2);
	printf("abs(x2-x1) is %lf \n",fabs(x2-x1));
	printf("dxold is %lf \n",dxold);
	printf("toto is %lf \n",pow(10,-5)-pow(10,-10));*/
      
      
      
      if (((*h_total3-xh)*df-f)*((*h_total3-xl)*df-f) >= 0.0 || fabs(2.0*f) > fabs(dxold*df)) 
	{
	  dxold=dxx;
	  dxx=0.5*(xh-xl);
	  *h_total3=xl+dxx;
	  if (xl == *h_total3) {
	    dum=-log10(*h_total3);
	    /*   printf("Exiting drtsafe at A on iteration %d ph =  %lf \n",j, dum);*/
	    return;
	    
	  }
	}
      else{
	dxold=dxx;
	dxx=f/df;
	tempp=*h_total3;
	*h_total3=*h_total3-dxx;
	if (tempp ==*h_total3) {
	  dum=-log10(*h_total3);
	  /*  printf("Exiting drtsafe at B on iteration %d ph =  %lf \n",j, dum);*/
	  return;
	}
      }
      /*   printf("dxx is %lf \n",dxx);
	   printf("xacc is %lf \n",xacc);
	   printf("dxx is %lf \n", fabs(dxx));*/
      
      if (fabs(dxx) <xacc) {
	/*dum=log10(-*h_total3);*/
	dum=-log10(*h_total3);
	/* printf("Exiting drtsafe at C on iteration %d ph =  %lf \n",j, dum) */
	   return;
      }
      ta_inter1(*h_total3,&f,&df,ff,k0,k1,k2,kb,kp1,kp2,kp3,ksi,kw,ks,kf,bt,st,ft,dic,pt,sit,ta);
      
      if(f < 0.0) {
	xl=*h_total3;
	fl=f;
      }
      else        {
	xh=*h_total3;
	fh=f;
      }
    }
  if (j == maxit-1)
    printf("Total number of iterations, %d exceeded. \n",maxit);  
} /*end iteration loop*/

  /*---------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------*/    
void  ta_inter1(double x,double *fn,double *df
		,double ff,double k0,double k1,double k2,double kb,double k1p,double k2p
		,double k3p,double ksi,double kw,double ks,double kf,double bt,double st,double ft,double dic,double pt, double sit, double ta)
/*---------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------*/
{
  /*    ! This routine expresses TA as a function of DIC, htotal and constants.
	! It also calculates the derivative of this function with respect to 
	! htotal. It is used in the iterative solution for htotal. In the call
	! "x" is the input value for htotal, "fn" is the calculated value for TA
	! and "df" is the value for dTA/dhtotal
	
	INPUT
	carbon chemistry constant: k1,k2,kp1,kp2,kp3,st,ks,sit_in,ksi
	DIC : dic_in  */
  /*  local variables */
  double x2, x3, k12, k12p, k123p, c, a, a2, da, b, b2, db;
  
  x2=x*x;
  
  x3 = x2*x;
  k12 = k1*k2;
  k12p = k1p*k2p;
  k123p = k12p*k3p;
  c = 1.0 + st/ks;
  a = x3 + k1p*x2 + k12p*x + k123p;
  a2 = a*a;
  da = 3.0*x2 + 2.0*k1p*x + k12p;
  b = x2 + k1*x + k12;
  b2 = b*b;
  db = 2.0*x + k1;
  
  *fn = k1*x*dic/b + 2.0*dic *k12/b +bt /(1.0 + x/kb ) +kw /x +pt *k12p*x/a +                        
    2.0*pt *k123p/a +  sit /(1.0 + x/ksi ) - x/c - st /(1.0 + ks /x/c) -                   
    ft / (1.0 + kf /x) - pt *x3/a - ta ;
  
  *df = ((k1  *dic  *b) - k1  *x*dic  *db)/b2 - 2.0*dic  *k12*db/b2 -bt  /kb  /         
    pow((1.0+x/kb  ),2) -kw  /x2 +(pt  *k12p*(a - x*da))/a2 -2.0*pt  *k123p*da/a2 -                 
    sit  /ksi  / pow((1.0+x/ksi  ),2) -1.0/c +  st  *pow((1.0 + ks  /x/c),(-2))*              
    (ks  /c/x2) + ft  *pow((1.0 + kf  /x),(-2))*kf  /x2 -pt  *x2*(3.0*a-x*da)/a2  ;  
  
  
} /*end function inter1 */ 
