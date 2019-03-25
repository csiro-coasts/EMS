/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/include/hd.h
 *  
 *  Description:
 *  Include file for 3-d non-linear hydrodynamic model
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: hd.h 5915 2018-09-05 03:30:40Z riz008 $
 *
 */

#ifndef _HD_H
#define _HD_H

/* include stdio for all files */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <semaphore.h>
#include "ems.h"

/* ANSI value may not be defined in old stdio.h */
#ifndef SEEK_SET
#define SEEK_SET	1
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#endif

/* include file for MAXSTRLEN */
#include "hd_params.h"

/* scheduler_t */
#include "scheduler.h"

/* Add the grid library */
#include "gridlib.h"

#if defined(HAVE_SEDIMENT_MODULE)
#include "sediments.h"
#endif

#if defined(HAVE_TRACERSTATS_MODULE)
#include "tracerstats.h"
#endif

#if defined(HAVE_WAVE_MODULE)
#include "waves.h"
#endif

/* Sparse arrays */
#include "sparse.h"

/* include file for useful macros */
#include "hd_macros.h"

#if defined(HAVE_ECOLOGY_MODULE)
#include "ecology.h"
#endif

#include "dp.h"

/* JIGSAW grid generation */
#ifdef HAVE_JIGSAWLIB
#include "lib_jigsaw.h"
#endif

/* prototypes */
#include "proto.h"

/* global variables */
#include "globals.h"

/* DO_TIMING guards are inside the header */
#include "timing.h"

#endif                          /* _HD_H */
