/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: lib/ecology/include/ecology_version.h
 *  
 *  Description:
 *  Ecology code version.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: ecology_version.h 6482 2020-03-17 02:53:09Z riz008 $
 *
 */

#if !defined(_ECOLOGY_VERSION_H)

/* Release verions and getters */
#define ECOLOGY_MAJOR_VERSION 1
#define ECOLOGY_MINOR_VERSION 2
#define ECOLOGY_PATCH_VERSION 0

int get_ecology_major_vers(void);
int get_ecology_minor_vers(void);
int get_ecology_patch_vers(void);

#define _ECOLOGY_VERSION_H
#endif
