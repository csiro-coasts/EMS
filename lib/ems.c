/**
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file ems.c
 *
 *  \brief Top-level EMS library initialisation routines
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: ems.c 5832 2018-06-26 23:49:51Z riz008 $
 */

#include <stdlib.h>
#include "ems.h"

int overwrite_out = 0;

void ems_init(FILE* fp, int id)
{
  /*UR added reading of ems log levels                      *
   * Leave this at the front to ensure proper logging       *
   * and shut  down of of logging and/or other cleaning up  *
   
   * There is no error validating in here, it is expected   *
   * that that has been done before                         */    
   
  char buf[MAXSTRLEN];
  char fbuf[MAXSTRLEN];

  if(prm_read_char(fp,"log_file", buf))
  {
    sprintf(fbuf, "%s", buf);
    if (id > 0)
      sprintf(fbuf, "%s_%d", buf, id);

    log_set_output(fbuf);
        emstag(LINFO,"lib:ems:init","Set log file: %s",buf); 
  }else
  {
    sprintf(fbuf, "runlog");
    if (id > 0)
      sprintf(fbuf, "runlog_%d", id);
    log_set_output(fbuf);
      emstag(LINFO,"lib:ems:init","Set log file to default");
  }
  
  if(prm_read_char(fp, "log_levels", buf))
  {
      log_init(buf);
      emstag(LINFO,"lib:ems:init","Set log levels to %s",buf);
  }else
  {
      log_init(NULL);
      emstag(LINFO,"lib:ems:init","Set log levels to default");
  }
  
  if(prm_read_char(fp, "overwrite_output", buf))
  {
    overwrite_out = atoi(buf);
    emstag(LINFO,"lib:ems:init","Set overwrite output to %s",buf);
  }
}


void ems_clear()
{
  log_end(); 
#ifdef DMALLOC
  dmalloc_verify(0L);
#endif 
}


void ems_usage()
{
  fprintf(stderr, "\nEMS LIB Usage:\n");
  fprintf(stderr, "-l [major,info,error,warn,debug,trace,metric,all]\n");
}


int overwrite_output()
{
  return overwrite_out; 
}

/* Version information */
int get_emslib_major_vers(void)
{
  return(EMSLIB_MAJOR_VERSION);
}

int get_emslib_minor_vers(void)
{
  return(EMSLIB_MINOR_VERSION);
}

int get_emslib_patch_vers(void)
{
  return(EMSLIB_PATCH_VERSION);
}
