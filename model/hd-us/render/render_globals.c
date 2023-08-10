/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  File: render_globals.c
 *
 *  Description:
 *  Model global configuration manager

 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: render_globals.c v1.4.0 rev(7068:7069) Wed Mar 16 13:48:05 2022
 her127 $
 *
 */

#include <stdio.h>
#include "hd.h"
#include "tracer.h"
#include "render.h"
#if defined(HAVE_ECOLOGY_MODULE)
#include "eprocess.h"
#endif

eco_global_t eco_glob[] = {
  {"NULL", "NULL", "NULL", NULL, 0, "NULL", NULL, 0, "NULL", NULL, 0}
};

eco_global_param_t eco_param[] = {
  {"NULL", "NULL", "NULL", NULL, 0}
};

sed_global_param_t sed_param[] = {
  {"NULL", "NULL", "NULL", NULL}
};

eco_global_proc_t eco_proc[] = {
  {"NULL", "NULL", "NULL", 0, "NULL", 0, "NULL", 0}
};

hydro_global_t hydro_param[] = {
  {"NULL", "NULL", NULL, NULL}
};

/*
 * This function is called from create_processes and serves as a
 * gateway for the different default processes types.
 */
int get_rendered_processes(char *name, int type, const char **procs[], int *nprocs)
{
  int i,n;

  struct {
    char *name;
    char *description;
    void (*init) (int type, const char **procs[], int *nprocs);
  } param_list[] = {
    {"NULL", "NULL", NULL}
  };
  void (*init) (int type, const char **procs[], int *nprocs) = NULL;

  for (i = 0; (init == NULL) && param_list[i].name; ++i) {
    if (strcasecmp(name, param_list[i].name) == 0) {
      init = param_list[i].init;
    }
  }
  if (init != NULL) {
    init(type, procs, nprocs);
  } else
    return(1);
  return(0);
}


