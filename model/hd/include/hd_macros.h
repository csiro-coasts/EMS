/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/include/hd_macros.h
 *  
 *  Description:
 *  Include file for hydro model macros
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: hd_macros.h 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

#ifndef _EMSHD_MACROS_H
#define _EMSHD_MACROS_H

#if !defined(MINVAL)
#define MINVAL (1.e-13)
#endif

#define validcelli(i)	((i>=0) && (i<params->nce1))
#define validcellj(j)	((j>=0) && (j<params->nce2))
#define validfacei(i)	((i>=0) && (i<params->nfe1))
#define validfacej(j)	((j>=0) && (j<params->nfe2))
#define validk(k)	((k>=0) && (k<params->nz))

#define validcell(i,j)	 (validcelli(i) && validcellj(j))
#define validu1face(i,j) (validfacei(i) && validcellj(j))
#define validu2face(i,j) (validcelli(i) && validfacej(j))

#endif
