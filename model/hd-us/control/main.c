/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/control/main.c
 *  
 *  Description:
 *  This file contains the main skeleton
 *  of the 3-d non-linear hydrodynamic
 *  model
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: main.c 6597 2020-09-03 05:29:26Z riz008 $
 *
 */

#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <time.h>
#include <execinfo.h>
#include "hd.h"
#include "tracer.h"
#include "ems_version.h"

#ifdef HAVE_PTHREADS
#include <pthread.h>
#endif

/* Proto types */
extern time_t time(time_t *);
void hd_step_trans(hd_data_t *hd_data, double tstop);
void hd_step_2d(hd_data_t *hd_data, double stop_at);
void hd_step(hd_data_t *hd_data, double stop_at);

void usage(void)
{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "compas -p prmfile [-restart] [-nrt] <options>\n");
  fprintf(stderr, "  Run COMPAS using standard parameter file.\n");
  fprintf(stderr, "  prmfile    :          Standard parameter file\n");
  fprintf(stderr, "  [-restart] :          Start using 'restart.nc' file\n");
  fprintf(stderr, "  [-nrt]     :          Near-real-time operation\n\n");
  fprintf(stderr, "compas -g prmfile dumpfile <options>\n");
  fprintf(stderr, "  Generate initial dump using standard parameter file.\n");
  fprintf(stderr, "  prmfile  :            Standard parameter file\n");
  fprintf(stderr, "  dumpfile :            Initialisation dump file\n\n");
  fprintf(stderr, "compas -p prmfile -cr <options>\n");
  fprintf(stderr, "  Run COMPAS using restarts when model crashes.\n");
  fprintf(stderr, "  prmfile  :            Standard parameter file\n\n");
  fprintf(stderr, "compas -a prmfile <options>\n");
  fprintf(stderr, "  Run COMPAS using autostart parameter file.\n");
  fprintf(stderr, "  prmfile  :            Autostart parameter file\n\n");
  fprintf(stderr, "compas -ag prmfile <options>\n");
  fprintf(stderr, "  Generate initial dump using autostart parameter file.\n");
  fprintf(stderr, "  prmfile  :            Autostart parameter file\n\n");
  fprintf(stderr, "compas -t prmfile <options>\n");
  fprintf(stderr, "  Run COMPAS in the transport mode.\n");
  fprintf(stderr, "  prmfile  :            Transport parameter file\n\n");
  fprintf(stderr, "compas -ps\n");
  fprintf(stderr, "  Generate percentile statistics.\n\n");
  fprintf(stderr, "compas -v\n");
  fprintf(stderr, "  Print compas version information.\n\n");
  fprintf(stderr, "options:\n");
  fprintf(stderr,
          "  -warnings on|off      Enable/disable warning messages (default: on).\n");
  fprintf(stderr,
          "  -diag_log <file>|off  Enable/disable diagnostic log (default: diag.txt).\n");
  fprintf(stderr,
          "  -setup_log <file>|off Enable/disable setup log (default: setup.txt).\n");
  fprintf(stderr,
          "  -window_log on|off    Enable/disable window log (default: off).\n");
  fprintf(stderr, "  -l tag,tag,...    Set library log level [info,main,warn,debug,trace,metric] (the more specialized the level the greater the output)\n");
  fprintf(stderr, "  -debug tag,tag,...    Set debug level\n");
  debug_usage();
  exit(1);
}


/* Version information */
static int get_major_vers(void)
{
  return(COMPAS_MAJOR_VERSION);
}

static int get_minor_vers(void)
{
  return(COMPAS_MINOR_VERSION);
}

static int get_patch_vers(void)
{
  return(COMPAS_PATCH_VERSION);
}

void print_vers(void)
{
  fprintf(stderr, "\nModules\n");
  fprintf(stderr, "-------\n");
  fprintf(stderr, "compas\t\t %d.%d.%d\n", get_major_vers(),
	  get_minor_vers(), get_patch_vers());
  fprintf(stderr, "emslib\t\t %d.%d.%d\n", get_emslib_major_vers(),
	  get_emslib_minor_vers(),get_emslib_patch_vers());
#ifdef HAVE_TRACERSTATS_MODULE
  fprintf(stderr, "tracerstats\t %d.%d.%d\n", get_tracerstats_major_vers(),
    get_tracerstats_minor_vers(),get_tracerstats_patch_vers());
#endif
#ifdef HAVE_WAVE_MODULE
  fprintf(stderr, "waves\t\t %d.%d.%d\n", get_waves_major_vers(),
    get_waves_minor_vers(),get_waves_patch_vers());
#endif
#ifdef HAVE_SEDIMENT_MODULE  
  fprintf(stderr, "sediments\t %d.%d.%d\n", get_sediments_major_vers(),
    get_sediments_minor_vers(),get_sediments_patch_vers());
#endif
#ifdef HAVE_ECOLOGY_MODULE
  fprintf(stderr, "ecology\t\t %d.%d.%d\n", get_ecology_major_vers(),
	  get_ecology_minor_vers(),get_ecology_patch_vers());
#endif
  fprintf(stderr,"\n");
}


void process_args(int argc, char *argv[])
{
  char buf[MAXSTRLEN];
  
  if (argc <= 1)
    usage();
  /*UR store the name of this executable for logging */  
  strcpy(executable, argv[0]);
  
  strcpy(diag_logfile, "diag.txt");
  strcpy(setup_logfile, "setup.txt");
  strcpy(window_geom_logfile, "window_geom.txt");
  strcpy(window_map_logfile, "window_map.txt");
#ifdef DO_TIMING
  strcpy(timing_logfile, "timing.txt");
#endif

  /* process the argument list */
  while (--argc > 0) {
    argv++;


    if (strcmp(*argv, "-p") == 0) {
      /* Standard run using parameter file initialisation           */
      if (argc < 2)
        usage();
      strcpy(prmname, *++argv);
      argc--;

    } else if (strcmp(*argv, "-g") == 0) {
      /* Create input file option                                   */
      if (argc < 3)
        usage();
      strcpy(prmname, *++argv);
      strcpy(oname, *++argv);
      autof = 2;
      argc -= 2;

    } else if (strcmp(*argv, "-a") == 0) {
      /* Standard automated initialisation and run                  */
      if (argc < 2)
        usage();
      strcpy(prmname, *++argv);
      autof = 1;
      argc--;

    } else if (strcmp(*argv, "-ag") == 0) {
      /* Standard automated initialisation; no run (exit after the  */
      /* creation of the input and parameter files).                */
      if (argc < 2)
        usage();
      strcpy(prmname, *++argv);
      autof = 3;
      argc--;

    } else if (strcmp(*argv, "-r") == 0) {
      /* ROAM automated initialisation and run                      */
      if (argc < 2)
        usage();
      strcpy(prmname, *++argv);
      autof = 4;
      argc--;

    } else if (strcmp(*argv, "-rg") == 0) {
      /* ROAM automated initialisation; no run (exit after the      */
      /* creation of the input and parameter files).                */
      if (argc < 2)
        usage();
      strcpy(prmname, *++argv);
      autof = 5;
      argc--;

    } else if (strcmp(*argv, "-rt") == 0) {
      /* ROAM PRE_MARVL - same as -rg plus generates transport file */
      /* This is a preprocessing step in MARVL                      */
      if (argc < 2)
        usage();
      strcpy(prmname, *++argv);
      autof = 10;
      argc--;

    } else if (strcmp(*argv, "-restart") == 0) {
      /* Restart the run */
      if (autof == 0 || autof == 8)
        forced_restart = 1;
      else
        usage();
      /*argc--;*/

    } else if (strcmp(*argv, "-cr") == 0) {
      /* Restart on crashes using parameter file initialisation     */
      if (autof == 0 || autof == 1 || autof == 4)
	crash_restart = 1;
      else
        usage();
      /*argc--;*/

    } else if (strcmp(*argv, "-rs") == 0) {
      /* ROAM automated initialisation with restart; variable       */
      /* initialisation uses data from a previous run.              */
      if (argc < 2)
        usage();
      strcpy(prmname, *++argv);
      autof = 6;
      argc--;

    } else if (strcmp(*argv, "-rso") == 0) {
      /* ROAM automated initialisation with restart; variable       */
      /* initialisation uses data from a previous run except for    */
      /* T and S (read from OFAM).                                  */
      if (argc < 2)
        usage();
      strcpy(prmname, *++argv);
      autof = 7;
      argc--;

    } else if (strcmp(*argv, "-t") == 0) {
      /* Transport option.                                          */
      if (argc < 2)
        usage();
      strcpy(prmname, *++argv);
      autof = 8;
      argc--;

    } else if (strcmp(*argv, "-ps") == 0) {
      /* Percentile calculations.                                   */
      if (argc < 2)
        usage();
      strcpy(prmname, *++argv);
      autof = 9;
      argc--;

    } else if (strcmp(*argv, "-nrt") == 0) {
      /* Near real time option                                      */
      if (autof == 0)
        nrt_restart = 1;
      else
        usage();
      /*argc--;*/

    } else if (strcmp(*argv, "-step") == 0) {
      /* Option to run exactly one time step                         */
      if (autof == 0)
        one_step = 1;
      else
        usage();
      /*argc--;*/

    } else if (strcmp(*argv, "-v") == 0) {
      print_vers();
      exit(1);
    } else if (strcmp(*argv, "-warnings") == 0) {
      if (argc < 2)
        usage();
      warnings = (strcmp(*++argv, "on") == 0);
      argc--;

    } else if (strcmp(*argv, "-diag_log") == 0) {
      if (argc < 2)
        usage();
      if (strcmp(*++argv, "off") != 0) {
        diag_log = 1;
        strcpy(diag_logfile, *argv);
      } else
        diag_log = 0;
      argc--;

    } else if (strcmp(*argv, "-setup_log") == 0) {
      if (argc < 2)
        usage();
      if (strcmp(*++argv, "off") != 0) {
        setup_log = 1;
        strcpy(setup_logfile, *argv);
      } else
        setup_log = 0;
      argc--;

    } else if (strcmp(*argv, "-window_log") == 0) {
      if (argc < 2)
        usage();
      windows_log = strcmp(*++argv, "on") == 0;
      argc--;

    } else if (strcmp(*argv, "-debug") == 0) {
      if (argc < 2)
        usage();
      debug_init(*++argv);
      argc--;

    } else if (strcmp(*argv, "-debug_log") == 0) {
      if (argc < 2)
        usage();
      debug_set_output(*++argv);
      argc--;

    } else if (strcmp(*argv, "-l") == 0) {
      if (argc < 2)
        usage();
      log_init(*++argv);
      argc--;

    }else if (strcmp(*argv, "-stat_log") == 0) {
      if (argc < 2)
        usage();
      if (strcmp(*++argv, "off") != 0) {
        stat_log = 1;
        strcpy(stat_logfile, *argv);
      } else
        stat_log = 0;
      argc--;

    }else if (strcmp(*argv, "-dev") == 0) {
      stat_log = 1;
      debug_init("all");
      strcpy(stat_logfile,"statistics.log");
      setup_log = 1;
      warnings = 1;
      diag_log = 1;
      windows_log =1;
      log_init("trace,debug");
      argc--;
    }else {
      usage();
    }
  }
}

/*
 * Routine to print out the stack trace. This is called from the
 * signal handler for a seg. fault
 */
void print_trace (void)
{
  void *array[30]; // 30 ought to be enough
  int size, i;
  char **strings;
  int k = 0;
  
  size    = backtrace(array, 30);
  strings = backtrace_symbols(array, size);

  hd_error("Segmentation violation detect (simulation time = %.4f days)\n", master->days);
  hd_error("Stack trace:\n", size);

  /*
   *  Skip the first 3 as they are not very useful :
   *    0 - this function
   *    1 - the calling function, i.e. kill_signal_handler
   *    2 - handler in libc
   *    3 - most likely culprit
   */
  for (i = 3; i < size; i++)
    hd_error(" [%d] %s\n", k++, strings[i]);
  
  free(strings);
}

void kill_signal_handler(int i)
{
   switch (i) {
   case SIGQUIT:
      hd_warn("A SIGQUIT signal has been received (model time = %g).\n"
              "Commencing graceful termination.", master->t);
      killed = 1;
      break;

   case SIGTERM:
      hd_warn("A SIGTERM signal has been received (model time = %g).\n"
              "Commencing graceful termination.", master->t);
      killed = 1;
      break;

   case SIGSEGV:
     // We can't reliably exit gracefully so just print out the stack
     // and exit with a non-zero status
     print_trace();
     exit(1);
     break;

   default:
     break;
   }
}

int main(int argc, char *argv[])
{
  time_t now;
  FILE *prmfd;
  int mpi_prov;
  
  killed = 0;
  model_running = 0;

#ifdef HAVE_MPI
  /* Must come early in the program */
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &mpi_prov);
  if (mpi_prov != MPI_THREAD_FUNNELED)
    fprintf(stderr, "MPI init thread error\n");
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif

  /* Capture the QUIT, TERMinate & SEGmentationViolation signals */
  signal(SIGQUIT, kill_signal_handler);
  signal(SIGTERM, kill_signal_handler);
  signal(SIGSEGV, kill_signal_handler);

  /* Constuct the version string */
  /* eg. "v1.1 rev(3393)" */
#if !EMS_IS_RELEASE
  /* Only display revision string in development mode */
  sprintf(version, "%s rev(%s)", EMS_VERSION, SVN_REV_STR);
#else
  sprintf(version, "%s", EMS_VERSION);
#endif
  
  /* print heading */
  now = time(NULL);
  fprintf(stderr, "\t\tCOMPAS\n");
  fprintf(stderr, "EMS Version:\t%s\n", version);

  /* get arguments from the command line or interactively */
  memset(prmname, 0, MAXSTRLEN);
  process_args(argc, argv);
  prm_set_case(0);

  fprintf(stderr, "Run start:\t%s\n", ctime(&now));
  
  /* MPI Info */
#ifdef HAVE_MPI
  fprintf(stderr, "\nMPI_RANK %d for %s(%d)\n", mpi_rank, prmname, getpid());
  fprintf(stderr, "\nMPI support is *only* for debugging purposes\n");
#endif

/* Open/Read in the parameter files and setup a scheduler
 */
  if ((prmfd = fopen(prmname, "r")) == NULL)
    hd_quit("Can't open parameter file '%s'.\n", prmname);
  /* params_read(prmfd); */

  if (autof == 9) calc_perc(prmfd);
  
  /*UR-ADDED to initialise any ems lib specific functions
   * leave this in front */
  ems_init(prmfd, 0);

  INIT_TIMING;

/* 
 * Schedule the events and the main data.
 * Maintain order.
 */
  schedule = sched_init(prmfd, now);
  TIMING_SET;
  hd_data = hd_init(prmfd);
  TIMING_DUMP(0, "hd_init");

  /* 
   * Start the main loop
   */

  if (DEBUG("time")) {
    dlog("time", "Start time = %f days.\n",
         schedule->start_time / 86400.0);
    dlog("time", "Stop time  = %f days.\n", schedule->stop_time / 86400.0);
  }

  model_running = 1;
  /*UR-FIX capture actual start time */
  schedule->exec_start_time = time(NULL);

  /*
   * Set marker for the runlog
   */
  emslog_set_log_offset("start of main simulation loop");

  while (!killed && (schedule->t < schedule->stop_time)) {
    double stop_at;

    stop_at = sched_get_next_event(schedule);

    /*
    if (DEBUG("time"))
      dlog("time", "Model time = %f days.\n", schedule->t / 86400.0);
    */
    emstag(LTRACE,"hd:main:main","Model time = %f days.\n", schedule->t / 86400.0);
    PRINT_TIMING_HEADER("main",master);


    if (hd_data->master->runmode == VEL3D) {
      TIMING_SET;
      hd_step(hd_data, stop_at);
      TIMING_DUMP(0, "hd_step");
    } else if (hd_data->master->runmode == VEL2D)
      hd_step_2d(hd_data, stop_at);
    else if (hd_data->master->runmode & TRANS) {
      TIMING_SET;
      hd_step_trans(hd_data, stop_at);
      TIMING_DUMP(0, "hd_step_trans");
    }

    TIMING_SET;
    sched_set_time(schedule, stop_at);
    TIMING_DUMP(0, "event");

    TIMING_COUNTER;
  }

  // Dump this time point when killed
  if (killed) {
    dump_snapshot(sched_get_even_by_name(schedule, "dumps"), schedule->t);
  }
    
  emstag(LINFO,"hd:main:main","Finishid at model time = %f days.\n", schedule->t / 86400.0);

  model_running = 0;


  /*
   * Clean up data assimilation
   */
#ifdef HAVE_DA
  dassim_end();
#endif

/* Free memory associated with the timeseries.
 */
  timeseries_end();

/* Destroy the grid, and deallocate any memory associated with them.
 */
  master_end(master);

/* Remove the time scheduler.
 */
  sched_end(schedule);

/* Close the parameter file.
 */
  fclose(prmfd);

  // wrap up timing
  END_TIMING;

  debug_end();

  /*UR-ADDED finally release the master */
  master_free(master);

  now = time(NULL);
  fprintf(stderr, "Run end:\t%s\n", ctime(&now));
  /*UR-ADDED to clear any ems lib specific functions */
  ems_clear();

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  exit(0);

}

/*UR-ADDED to comply with ems lib addition
 * @depreciated
 */
int verbosity()
{
  return 0;
}

