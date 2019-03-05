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
 *  $Id: ecology_version.h 6156 2019-03-05 03:00:49Z riz008 $
 *
 */

#if !defined(_ECOLOGY_VERSION_H)

/* Release verions and getters */
#define ECOLOGY_MAJOR_VERSION 1
#define ECOLOGY_MINOR_VERSION 1
#define ECOLOGY_PATCH_VERSION 1

int get_ecology_major_vers(void);
int get_ecology_minor_vers(void);
int get_ecology_patch_vers(void);

#define _ECOLOGY_VERSION_H
#endif
