/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  File: model/lib/grid/include/gridlib_version.h
 *  
 *  Description:  Grid library code version
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: gridlib_version.h 6155 2019-03-05 02:57:41Z riz008 $
 *
 */

#if !defined(_GRIDLIB_VERSION_H)

/* Release verions and getters */
#define GRIDLIB_MAJOR_VERSION 1
#define GRIDLIB_MINOR_VERSION 1
#define GRIDLIB_PATCH_VERSION 0

int get_gridlib_major_vers(void);
int get_gridlib_minor_vers(void);
int get_gridlib_patch_vers(void);

#define _GRIDLIB_VERSION_H
#endif
