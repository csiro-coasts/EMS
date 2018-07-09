/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/misc/timing.c
 *  
 *  Description:
 *  Functions to handle timing profile
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: timing.c 5873 2018-07-06 07:23:48Z riz008 $
 *
 */

#include <stdio.h>
#include "hd.h"

#ifdef DO_TIMING
/*
 * Initialises some stuff and prints the main banner
 */
void init_timing(void)
{
  FILE *fp = stderr;
  
  if (timing_logfile != NULL) {
    // This assigns the global var
    tfp = fopen(timing_logfile, "w");
    if (tfp == NULL) {
      // quitting here is probably ok as this shouldn't be called for
      // a production run
      hd_quit("Unable to open the timing log file '%s'\n", timing_logfile);
    }
    fp = tfp;
  }

  // init timing level from the environment variable
  char *tl_env = getenv("EMS_TIMING_LEVEL");
  if (tl_env != NULL) {
    timing_level = atoi(tl_env);
  }
   
  // Print the banner 
  fprintf(fp,"-------------------\n");
  fprintf(fp,"SHOC Timing profile\n");
  fprintf(fp,"-------------------\n");
  fprintf(fp,"\n");
  fprintf(fp,"Each column of numbers corresponds to another functional level.\n");
  fprintf(fp,"To keep the log file size from getting too large, these\n");
  fprintf(fp,"profiles wrap after %d steps of the 'main' loop.\n", TIMING_NUM_STEPS);
  fprintf(fp,"\n");
  fprintf(fp,"=============================================================\n");

  // Set the position
  tfp_pos = ftell(fp);
}

/*
 * Prints a header on top of loops
 */
void print_timing_header(char *str, void *m)
{
  master_t *master = (master_t *)m;
  int this_time    = (int)master->t;
  static int prev_time = 0;
  int dt = this_time - prev_time;

  fprintf(tfp, "\nSimulation time : %d secs (%.4f days) - %d secs from last\n", 
	  this_time, master->days, (dt==this_time ? 0 : dt)); 
  fprintf(tfp, "---------------\n");
  if (strcmp(str, "main") == 0) {
    fprintf(tfp, "%s (%d)\n", str, timing_counter);
  } else {
    fprintf(tfp, "%s\n", str);
  }
  fprintf(tfp, "===============\n");

  prev_time = this_time;
}
#endif

// EOF
