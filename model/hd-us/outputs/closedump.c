/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/outputs/closedump.c
 *  
 *  Description:
 *  Close the netCDF dump file
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: closedump.c 5873 2018-07-06 07:23:48Z riz008 $
 *
 */

#include <stdio.h>
#include <netcdf.h>
#include "hd.h"

void dump_close(int cdfid)
{
  nc_close(cdfid);
}
