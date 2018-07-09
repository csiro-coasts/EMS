/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/debug/printflags.c
 *  
 *  Description:
 *  Print out flags to the specified stream.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: printflags.c 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

#include <math.h>
#include <stdio.h>
#include "hd.h"

void printflags(FILE * out, unsigned long f)
{
  if (f & U1SOLID)
    fprintf(out, "U1SOLID ");
  if (f & U2SOLID)
    fprintf(out, "U2SOLID ");
  if (f & U1OUTSIDE)
    fprintf(out, "U1OUTSIDE ");
  if (f & U2OUTSIDE)
    fprintf(out, "U2OUTSIDE ");
  if (f & U1BDRY)
    fprintf(out, "U1BDRY ");
  if (f & U2BDRY)
    fprintf(out, "U2BDRY ");
  if (f & L_EDGE)
    fprintf(out, "L_EDGE ");
  if (f & R_EDGE)
    fprintf(out, "R_EDGE ");
  if (f & B_EDGE)
    fprintf(out, "B_EDGE ");
  if (f & F_EDGE)
    fprintf(out, "F_EDGE ");
  if (f & ETASPEC)
    fprintf(out, "ETASPEC ");
  if (f & DRY)
    fprintf(out, "DRY (No longer used) ");
  if (f & SOLID)
    fprintf(out, "SOLID ");
  if (f & OUTSIDE)
    fprintf(out, "OUTSIDE ");
  if (f & ALLWATER)
    fprintf(out, "ALLWATER ");
  if (f & TRACERSPEC)
    fprintf(out, "TRACERSPEC ");
  if (f & U1TFLUX)
    fprintf(out, "U1TFLUX ");
  if (f & U2TFLUX)
    fprintf(out, "U2TFLUX ");
  if (f & U1AVZERO)
    fprintf(out, "U1AVZERO ");
  if (f & U2AVZERO)
    fprintf(out, "U2AVZERO ");
  fprintf(out, "\n");
}
