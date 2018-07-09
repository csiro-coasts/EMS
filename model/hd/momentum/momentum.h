/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/momentum/momentum.h
 *  
 *  Description:
 *  Include file for momentum equation routines
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: momentum.h 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

/* Use wet points only code */
#define WET_POINTS_ONLY

/* Sets the formulation for horizontal mixing             */
/* HOR_MIX_METHOD = 1 : Complete formulation with metrics */
/* HOR_MIX_METHOD = 0 : Simple formulation                */
#define HOR_MIX_METHOD 1
