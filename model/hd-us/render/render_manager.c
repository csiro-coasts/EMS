/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/ecology/render_manager.c
 *  
 *  Description:
 *  Writes a module of the current ecology & sediment configuration.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: render_manager.c 6721 2021-03-29 04:38:30Z her127 $
 *
 */

#include "hd.h"
#include <assert.h>
#include "render.h"

void write_module_header(FILE *fp, char *name, char *desc);
void write_modules(master_t *master, char *name);
void write_rendered_globals(master_t *master, char *name, int existf,
			    int ntr, int nep, int nsp, int nparam,
			    int nwc, int nsed, int nepi, 
			    char *trinfo_name, char *ecopriv_name,
			    char *sedpriv_name);
void render_eco_process(master_t *master, void *e, FILE *fp, char *name, 
			int *wc, int *sed, int *epi);
void render_eco_types(FILE *fp, char *name);
void render_eco_gateway(master_t *master, FILE *fp, char *name, char *desc, int *remp);
int render_eco_parameters(master_t *master, void *ei, FILE *fp, char *name);
void render_sed_params(FILE *fp, char *name);
void render_hydro_parameters(master_t *master, FILE *fp, char *name);
void render_OBC_parameters(master_t *master, FILE *fp, char *name);
void set_global_ecomodule(master_t *master, char *name);
#if defined(HAVE_ECOLOGY_MODULE)
void write_private_eco(FILE *fp, trinfo_priv_eco_t *data, int mode);
void copy_ecoparameter_info(parameter_info *out, parameter_glob_t *in);
int get_rendered_ecoparams(char *name, parameter_info *parameters[], int *nprm);
#endif
#if defined(HAVE_SEDIMENT_MODULE)
void write_private_sed(FILE *fp, trinfo_priv_sed_t *data, int mode);
void copy_sedparameter_info(sed_params_t *out, sed_params_t *in);
int get_rendered_sedparams(char *name, sediment_t *sediment, int mode);
extern sed_global_param_t sed_param[];
#endif
void copy_hydroparameter_info(parameters_t *params, parameters_t *in);
int get_rendered_hydroparams(char *name, parameters_t *in);
void render_list(master_t *master);
void render_dump(master_t *master);
void ecosed_config_write(master_t *master, char *name);
extern hydro_global_t hydro_param[];
extern eco_global_t eco_glob[];
extern eco_global_proc_t eco_proc[];
extern eco_global_param_t eco_param[];


/*------------------------------------------------------------------*/
/* Entry routine for render management                              */
/*------------------------------------------------------------------*/
void render_manage(master_t *master)
{
  if (master->renderopts & R_LIST) render_list(master);

  if (master->renderopts & R_DUMP) render_dump(master);

  if (strlen(master->rendername) || strlen(master->renderrem))
    write_modules(master, master->rendername);
}

/* END render_manage()                                              */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Writes available rendered configurations                         */
void render_list(master_t *master)
{
  FILE *fp;
  char buf[MAXSTRLEN];
  int n, m;

  /*----------------------------------------------------------------*/
  /* List available configurations if required                      */
  if (strlen(master->rendername)) 
    sprintf(buf, "%s_list.txt", master->rendername);
  else
    sprintf(buf, "config_list.txt");
  fp = fopen(buf, "w");
  for (m = 0; strcmp(hydro_param[m].name, "NULL") != 0; ++m) 
    fprintf(fp, "HYDRO%d %s : %s\n", n, hydro_param[n].name, 
	    hydro_param[n].desc);
#if defined(HAVE_ECOLOGY_MODULE)
  for (n = 0; strcmp(eco_glob[n].name, "NULL") != 0; ++n)
    fprintf(fp, "ECO%d %s : %s\n", n, eco_glob[n].name, 
	    eco_glob[n].desc);
#endif
  if (n == 0 && m == 0) 
    fprintf(fp, "There are currently no configurations in the list.\n");
  fclose(fp);
}

/* END render_list()                                                */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Dumps a rendered configuration to parameter file                 */
void render_dump(master_t *master)
{
  parameters_t *params = master->params;
  dump_data_t *dumpdata = master->dumpdata;
  char buf[MAXSTRLEN];

  if (master->rendertype & R_HYDRO) {
    if (strlen(master->rendername)) 
      sprintf(buf, "%s_hydro", master->rendername);
    else
      sprintf(buf, "hydro");
    params_write(params, dumpdata, buf);
  }

  if (master->rendertype & (R_ECO|R_SED)) {
    if (strlen(master->rendername)) 
      sprintf(buf, "%s_ecosed.prm", master->rendername);
    else
      sprintf(buf, "ecosed.prm");
    ecosed_config_write(master, buf);
  }
}

/* END render_dump()                                                */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Writes the sediment and ecology specification to parameter       */
/* file.                                                            */
/*------------------------------------------------------------------*/
void ecosed_config_write(master_t *master, char *name)
{
  parameters_t *params = master->params;
  FILE *fp;
  int n, ntr;

  fp = fopen(name, "w");

  /* Hydro tracers                                                  */
  ntr = tracer_write(params, fp);

  /* Sediment tracers                                               */
#if defined(HAVE_SEDIMENT_MODULE)
  if (params->do_sed) {
    ntr = sediment_autotracer_write(master, fp, ntr);
    fprintf(fp, "\n");
  }
#endif

  /* Ecology tracers                                                */
#if defined(HAVE_ECOLOGY_MODULE)
  if (params->do_eco) {
    ntr = ecology_autotracer_write(master, fp, ntr);
  }
#endif

  /* Sediment parameters                                            */
#if defined(HAVE_SEDIMENT_MODULE)
  if (params->do_sed) {
    win_priv_t *wincon = hd_data->wincon[1];
    fprintf(fp, "############################################################################\n");
    fprintf(fp, "# Sediment specification\n\n");
    trans_write_sed(params, wincon->sediment, fp);
    fprintf(fp, "\n");
  }
#endif

  /* Ecology parameters and processes                               */
#if defined(HAVE_ECOLOGY_MODULE)
  if (params->do_eco) {
    int type;
    win_priv_t *wincon = hd_data->wincon[1];
    ecology *e = wincon->e;
    parameter_info *prm = e->pinfo;
    fprintf(fp, "############################################################################\n");
    fprintf(fp, "# Ecology specification\n\n");
    trans_write_eco(params->eco_vars, params->eco_defs, params->ecodt, wincon->e, fp);
    fprintf(fp, "\n");

    /* Ecology processes                                            */
    for (type = 0; type < N_EPROCESS_TYPES; ++type) {
      int n = e->npr[type];
      int i, mode;
      if (n > 0) {
	if (type == 1)
	  fprintf(fp, "water\n");
	if (type == 2)
	  fprintf(fp, "sediment\n");
	if (type == 3)
	  fprintf(fp, "epibenthos\n");
	fprintf(fp, "{\n");
	for (i = 0; i < n; ++i) {
	  fprintf(fp, "  %s\n", e->processes[type][i]->fullname);
	}
	fprintf(fp, "}\n\n");
      }
    }

    /* Ecology parameters                                           */
    fprintf(fp, "\nNPARAMETERS %d\n\n", e->nprm);
    for (n = 0; n < e->nprm; n++) {
      fprintf(fp, "PARAMETER%d.name %s\n", n, prm[n].name);
      fprintf(fp, "PARAMETER%d.desc %s\n", n, prm[n].desc);
      fprintf(fp, "PARAMETER%d.units %s\n", n, prm[n].units);
      fprintf(fp, "PARAMETER%d.value %s\n", n, prm[n].stringvalue);
      fprintf(fp, "PARAMETER%d.adjust 0\n\n", n);
    }
  }
#endif

  fclose(fp);
}

/* END ecosed_config_write()                                        */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Writes a C module containing a hydrodynamic, sediment or ecology */
/* configuration.                                                   */
/*------------------------------------------------------------------*/
void write_modules(master_t *master, char *name)
{
  parameters_t *params = master->params;
  FILE *fp;
  char buf[MAXSTRLEN];
  char trinfo_name[MAXSTRLEN];
  char ecopriv_name[MAXSTRLEN];
  char sedpriv_name[MAXSTRLEN];
  int n, m, ntr, nep, nsp;
  int nwc, nsed, nepi;
  int nparam;
  int existf = 0;

  /*----------------------------------------------------------------*/
  /* Check this name doesn't exist                                  */
  if (strlen(name) == 0) {
    existf = 1;
    write_rendered_globals(master, name, existf, 0, 0, 0, 0, 
			   0, 0, 0, NULL, NULL, NULL);
    return;
  } else {
    for (n = 0; strcmp(eco_glob[n].name, "NULL") != 0; ++n) {
      if (strcmp(eco_glob[n].name, name) == 0) {
	hd_warn("write_modules: name '%s' already rendered!\n", name);
	existf = 1;
      }
    }
  }

  /*----------------------------------------------------------------*/
  /* Open the configuration header                                  */
  if (strlen(master->renderpath))
    if (endswith(master->renderpath, "/"))
      sprintf(buf, "%s%s.h", master->renderpath, name);
    else
      sprintf(buf, "%s/%s.h", master->renderpath, name);
  else
    sprintf(buf, "%s.h", name);
  fp = fopen(buf, "w");

  /*----------------------------------------------------------------*/
  /* Write the configuration to a header file.                      */
  sprintf(buf, "Model configuration for %s: %s\n", name, master->renderdesc);
  write_module_header(fp, name, buf);
  fprintf(fp,"#include \"ems.h\"\n\n");

#if defined(HAVE_ECOLOGY_MODULE)
  if (master->rendertype & R_ECO) {
    render_eco_process(master, params->pre_eco, fp, name, &nwc, &nsed, &nepi);
    fprintf(fp, "\n");
  }
#endif

#if defined(HAVE_SEDIMENT_MODULE)
  if (master->rendertype & R_SED) {
    render_sed_params(fp, name);
    fprintf(fp, "\n");
  }
#endif

#if defined(HAVE_ECOLOGY_MODULE)
  if (master->rendertype & R_ECO) {
    nparam = render_eco_parameters(master, params->pre_eco, fp, name);
    fprintf(fp, "\n");
  }
#endif

  /*----------------------------------------------------------------*/
  /* Count the ecology tracers                                      */
  ntr = 0;
  if (master->rendertype & (R_ECO|R_SED)) {
    for (n = 0; n < master->ntr; n++) {
      tracer_info_t *tr = &master->trinfo_3d[n];
      if (tr->type & (ECOLOGY|SEDIMENT)) ntr++;
    }
    for (n = 0; n < master->ntrS; n++) {
      tracer_info_t *tr = &master->trinfo_2d[n];
    if (tr->type & (ECOLOGY|SEDIMENT)) ntr++;
    }
    fprintf(fp, "\nint NTR_%s = %d;\n", name, ntr);
    
    /* Write the tracer attributes                                  */
    sprintf(trinfo_name, "tracerlist_%s", name);
    fprintf(fp, "tracer_info_t %s[] = {\n", trinfo_name);
    for (n = 0; n < master->ntr; n++) {
      tracer_info_t *tr = &master->trinfo_3d[n];
      m = (n == ntr) ? 0 : 1;
      if (tr->type & (ECOLOGY|SEDIMENT)) write_auto_atts(fp, tr, m);
    }
    for (n = 0; n < master->ntrS; n++) {
      tracer_info_t *tr = &master->trinfo_2d[n];
      m = (n == ntr) ? 0 : 1;
      if (tr->type & (ECOLOGY|SEDIMENT)) write_auto_atts(fp, tr, m);
    }
    fprintf(fp, "};\n");
  }

  /*----------------------------------------------------------------*/
  /* Count the ecology tracers with private ecology data            */
#if defined(HAVE_ECOLOGY_MODULE)
  nep = 0;
  if (master->rendertype & R_ECO) {
    for (n = 0; n < master->ntr; n++) {
      tracer_info_t *tr = &master->trinfo_3d[n];
      trinfo_priv_eco_t *data = tr->private_data[TR_PRV_DATA_ECO];
      if (data) nep++;
    }
    for (n = 0; n < master->ntrS; n++) {
      tracer_info_t *tr = &master->trinfo_2d[n];
      trinfo_priv_eco_t *data = tr->private_data[TR_PRV_DATA_ECO];
      if (data) nep++;
    }
    fprintf(fp, "\n#if defined(HAVE_ECOLOGY_MODULE)\n");
    fprintf(fp, "\nint NPRIV_ECO_%s = %d;\n", name, nep);
    /* Write the private data                                       */
    sprintf(ecopriv_name, "private_eco_%s", name);
    fprintf(fp, "trinfo_priv_eco_t %s[] = {\n", ecopriv_name);
    for (n = 0; n < master->ntr; n++) {
      tracer_info_t *tr = &master->trinfo_3d[n];
      trinfo_priv_eco_t *data = tr->private_data[TR_PRV_DATA_ECO];
      if (data && master->rendertype & R_ECO) {
	m = (n == nep) ? 0 : 1;
	strcpy(data->trname, tr->name);
	write_private_eco(fp, data, m);
      }
    }
    for (n = 0; n < master->ntrS; n++) {
      tracer_info_t *tr = &master->trinfo_2d[n];
      trinfo_priv_eco_t *data = tr->private_data[TR_PRV_DATA_ECO];
      if (data && master->rendertype & R_ECO) {
	m = (n == nep) ? 0 : 1;
	strcpy(data->trname, tr->name);
	write_private_eco(fp, data, m);
      }
    }
    fprintf(fp, "};\n");
    fprintf(fp, "#endif\n");
  }
#endif

  /*----------------------------------------------------------------*/
  /* Count the ecology tracers with private sediment data           */
#if defined(HAVE_SEDIMENT_MODULE)
  nsp = 0;
  if (master->rendertype & (R_ECO|R_SED)) {
    for (n = 0; n < master->ntr; n++) {
      tracer_info_t *tr = &master->trinfo_3d[n];
      trinfo_priv_sed_t *data = tr->private_data[TR_PRV_DATA_SED];
      if (data) nsp++;
    }
    for (n = 0; n < master->ntrS; n++) {
      tracer_info_t *tr = &master->trinfo_2d[n];
      trinfo_priv_sed_t *data = tr->private_data[TR_PRV_DATA_SED];
      if (data) nsp++;
    }
    fprintf(fp, "\n#if defined(HAVE_SEDIMENT_MODULE)\n");
    fprintf(fp, "\nint NPRIV_SED_%s = %d;\n", name, nsp);
    /* Write the private data                                       */
    sprintf(sedpriv_name, "private_sed_%s", name);
    fprintf(fp, "trinfo_priv_sed_t %s[] = {\n", sedpriv_name);
    for (n = 0; n < master->ntr; n++) {
      tracer_info_t *tr = &master->trinfo_3d[n];
      trinfo_priv_sed_t *data = tr->private_data[TR_PRV_DATA_SED];
      if (data && master->rendertype & (R_ECO|R_SED)) {
	m = (n == nsp) ? 0 : 1;
	strcpy(data->trname, tr->name);
	write_private_sed(fp, data, m);
      }
    }
    for (n = 0; n < master->ntrS; n++) {
      tracer_info_t *tr = &master->trinfo_2d[n];
      trinfo_priv_sed_t *data = tr->private_data[TR_PRV_DATA_SED];
      if (data && master->rendertype & (R_ECO|R_SED)) {
	m = (n == nsp) ? 0 : 1;
	strcpy(data->trname, tr->name);
	write_private_sed(fp, data, m);
      }
    }
    fprintf(fp, "};\n");
    fprintf(fp, "#endif\n");
  }
#endif

  /* Write the hydrodynamic configuration                           */
  if (master->rendertype & R_HYDRO) {
    render_hydro_parameters(master, fp, name);
    render_OBC_parameters(master, fp, name);
  }
  fclose(fp);

  write_rendered_globals(master, name, existf, ntr, nep, nsp, nparam, 
			 nwc, nsed, nepi, trinfo_name, ecopriv_name,
			 sedpriv_name);
}

/* END write_modules()                                              */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Writes the global rendered configurations to a C module.         */
/*------------------------------------------------------------------*/
void write_rendered_globals(master_t *master, 
			    char *name, /* Name of config */
			    int existf, /* cofig already exists */
			    int ntr,    /* # tracers */
			    int nep,    /* tracers with eco priv data */
			    int nsp,    /* tracers with sediment priv data */
			    int nparam, /* # eco parameters */
			    int nwc,    /* # eco WC processes */
			    int nsed,   /* # eco sediment processes */
			    int nepi,   /* # eco epibenthos processes */
			    char *trinfo_name, /* Name of tracer trinfo */
			    char *ecopriv_name, /* Name of eco priv data */
			    char *sedpriv_name  /* Name of sed priv data */
			    )
{
  FILE *fp;
  char buf[MAXSTRLEN];
  char desc[MAXSTRLEN];
  char wcproc_name[MAXSTRLEN];
  char sedproc_name[MAXSTRLEN];
  char epiproc_name[MAXSTRLEN];
  char param_name[MAXSTRLEN];
  char *rem[MAXSTRLEN];
  int hasname = (strlen(name) == 0) ? 0 : 1;
  int n, m, nr, remf, *reme, *remp, *remh;

  /*----------------------------------------------------------------*/
  /* Open the configuration header                                  */
  if (strlen(master->renderpath))
    if (endswith(master->renderpath, "/"))
      sprintf(buf, "%srender_globals.c", master->renderpath);
    else
      sprintf(buf, "%s/render_globals.c", master->renderpath);
  else
    sprintf(buf, "render_globals.c");
  fp = fopen(buf, "w+");
  if (strlen(master->renderdesc))
    strcpy(desc, master->renderdesc);
  else
    strcpy(desc, name);

  /*----------------------------------------------------------------*/
  /* Any configurations to remove?                                  */
  if (strlen(master->renderrem))
    nr = parseline(master->renderrem, rem, MAXNUMARGS);
  else
    nr = 0;
  m = 1;
  for (n = 0; strcmp(eco_glob[n].name, "NULL") != 0; ++n)
    m++;
  reme = i_alloc_1d(m);
  for (n = 0; strcmp(eco_glob[n].name, "NULL") != 0; ++n) {
    reme[n] = 0;
    for (m = 0; m < nr; m++) {
      if (strcmp(eco_glob[n].name, rem[m]) == 0) {
	reme[n] = 1;
	break;
      }
    }
  }
  m = 1;
  for (n = 0; strcmp(eco_proc[n].name, "NULL") != 0; ++n)
    m++;
  remp = i_alloc_1d(m);
  for (n = 0; strcmp(eco_proc[n].name, "NULL") != 0; ++n) {
    remp[n] = 0;
    for (m = 0; m < nr; m++) {
      if (strcmp(eco_proc[n].name, rem[m]) == 0) {
	remp[n] = 1;
	break;
      }
    }
  }
  m = 1;
  for (n = 0; strcmp(hydro_param[n].name, "NULL") != 0; ++n)
    m++;
  remh = i_alloc_1d(m);
  for (n = 0; strcmp(hydro_param[n].name, "NULL") != 0; ++n) {
    remh[n] = 0;
    for (m = 0; m < nr; m++) {
      if (strcmp(hydro_param[n].name, rem[m]) == 0) {
	remh[n] = 1;
	break;
      }
    }
  }

  /*----------------------------------------------------------------*/
  /* Write the headers                                              */
  sprintf(buf, "Model global configuration manager\n");
  write_module_header(fp, "render_globals.c", buf);
  fprintf(fp, "#include <stdio.h>\n");
  fprintf(fp, "#include \"hd.h\"\n");
  fprintf(fp, "#include \"tracer.h\"\n");
  fprintf(fp, "#include \"render.h\"\n");
  fprintf(fp, "#if defined(HAVE_ECOLOGY_MODULE)\n");
  fprintf(fp, "#include \"eprocess.h\"\n");
  fprintf(fp, "#endif\n");
  for (n = 0; strcmp(hydro_param[n].name, "NULL") != 0; ++n) {
    if (strcmp(hydro_param[n].name, name) != 0 && !existf && !remh[n])
      fprintf(fp, "#include \"%s.h\"\n", hydro_param[n].name);
  }
  for (n = 0; strcmp(eco_glob[n].name, "NULL") != 0; ++n) {
    if (strcmp(eco_glob[n].name, name) != 0 && !existf && !reme[n])
      fprintf(fp, "#include \"%s.h\"\n", eco_glob[n].name);
  }
  if (hasname)
    fprintf(fp, "#include \"%s.h\"\n", name);
  fprintf(fp, "\n");

  sprintf(wcproc_name, "WC_PROCESS_%s", name);
  sprintf(sedproc_name, "SED_PROCESS_%s", name);
  sprintf(epiproc_name, "EPI_PROCESS_%s", name);

  /*----------------------------------------------------------------*/
  /* Tracer attributes                                              */
  fprintf(fp, "eco_global_t eco_glob[] = {\n");
  /* Print existing lists                                           */
  for (n = 0; strcmp(eco_glob[n].name, "NULL") != 0; ++n) {
    if (strcmp(eco_glob[n].name, name) != 0 && !existf && !reme[n])
      fprintf(fp, "  {\"%s\", \"%s\", \"%s\", %s, %d, \"%s\", %s, %d, \"%s\", %s, %d},\n",
	      eco_glob[n].name, eco_glob[n].desc,
	      eco_glob[n].trname, eco_glob[n].trname, eco_glob[n].ntr,
	      eco_glob[n].econame, eco_glob[n].econame, eco_glob[n].nep,
	      eco_glob[n].sedname, eco_glob[n].sedname, eco_glob[n].nsp);
  }
  /* Print the list asscociated with 'name' if required             */
  if (hasname && master->rendertype & (R_ECO|R_SED)) {
    fprintf(fp, "  {\"%s\", \"%s\", \"%s\", %s, %d, \"%s\", %s, %d, \"%s\", %s, %d},\n",
	    name, desc, trinfo_name, trinfo_name, ntr,
	    ecopriv_name, ecopriv_name, nep,
	    sedpriv_name, sedpriv_name, nsp);
  }
  fprintf(fp, "  {\"NULL\", \"NULL\", \"NULL\", NULL, 0, \"NULL\", NULL, 0, \"NULL\", NULL, 0}\n");
  fprintf(fp, "};\n\n");

  /*----------------------------------------------------------------*/
  /* Eco parameters                                                 */
  sprintf(param_name, "PARAMETERS_%s", name);
  fprintf(fp, "eco_global_param_t eco_param[] = {\n");
  for (n = 0; strcmp(eco_param[n].name, "NULL") != 0; ++n) {
    remf = 0;
    for (m = 0; m < nr; m++) {
      if (strcmp(eco_param[n].name, rem[m])) {
	remf = 1;
	break;
      }
    }
    if (strcmp(eco_param[n].name, name) != 0 && !existf && !remf)
      fprintf(fp, "  {\"%s\", \"%s\", \"%s\", %s, %d},\n",
	      eco_param[n].name, eco_param[n].desc,
	      eco_param[n].procname, eco_param[n].procname, eco_param[n].nparam);
  }
  if (hasname  && master->rendertype & R_ECO) {
    fprintf(fp, "  {\"%s\", \"%s\", \"%s\", %s, %d},\n", 
	    name, desc, param_name, param_name, nparam);
  }
  fprintf(fp, "  {\"NULL\", \"NULL\", \"NULL\", NULL, 0}\n");
  fprintf(fp, "};\n\n");

  /*----------------------------------------------------------------*/
  /* Sediments parameters                                           */
  sprintf(param_name, "SEDPARAMS_%s", name);
  fprintf(fp, "sed_global_param_t sed_param[] = {\n");
#if defined(HAVE_SEDIMENT_MODULE)
  for (n = 0; strcmp(sed_param[n].name, "NULL") != 0; ++n) {
    remf = 0;
    for (m = 0; m < nr; m++) {
      if (strcmp(sed_param[n].name, rem[m])) {
	remf = 1;
	break;
      }
    }
    if (strcmp(sed_param[n].name, name) != 0 && !existf && !remf)
      fprintf(fp, "  {\"%s\", \"%s\", \"%s\", &%s},\n",
	      sed_param[n].name, sed_param[n].desc, 
	      sed_param[n].procname, sed_param[n].procname);
  }
#endif
  if (hasname && master->rendertype & R_SED) {
    fprintf(fp, "  {\"%s\", \"%s\", \"%s\", &%s},\n", name, desc, 
	    param_name, param_name);
  }
  fprintf(fp, "  {\"NULL\", \"NULL\", \"NULL\", NULL}\n");
  fprintf(fp, "};\n\n");

  /*----------------------------------------------------------------*/
  /* Eco processes                                                  */
  fprintf(fp, "eco_global_proc_t eco_proc[] = {\n");
  for (n = 0; strcmp(eco_proc[n].name, "NULL") != 0; ++n) {
    if (strcmp(eco_proc[n].name, name) != 0 && !existf && !remp[n])
      fprintf(fp, "  {\"%s\", \"%s\", \"%s\", %d, \"%s\", %d, \"%s\", %d},\n",
	      eco_proc[n].name, eco_proc[n].desc,
	      eco_proc[n].wcname, eco_proc[n].nwc_proc,
	      eco_proc[n].sedname, eco_proc[n].nsed_proc,
	      eco_proc[n].epiname, eco_proc[n].nepi_proc);
  }
  if (hasname && master->rendertype & R_ECO) {
    fprintf(fp, "  {\"%s\", \"%s\", \"%s\", %d, \"%s\", %d, \"%s\", %d},\n",
	    name, desc, wcproc_name, nwc,
	    sedproc_name, nsed,
	    epiproc_name, nepi);
  }
  fprintf(fp, "  {\"NULL\", \"NULL\", \"NULL\", 0, \"NULL\", 0, \"NULL\", 0}\n");
  fprintf(fp, "};\n\n");


  /*----------------------------------------------------------------*/
  /* Hydro parameters                                               */
  sprintf(param_name, "HYDRO_%s", name);
  fprintf(fp, "hydro_global_t hydro_param[] = {\n");
  for (n = 0; strcmp(hydro_param[n].name, "NULL") != 0; ++n) {
    if (strcmp(hydro_param[n].name, name) != 0 && !existf && !remh[n])
      fprintf(fp, "  {\"%s\", \"%s\", &HYDRO_%s, HYDRO_%s_obc},\n",
	      hydro_param[n].name, hydro_param[n].desc, hydro_param[n].name, hydro_param[n].name);
  }
  if (hasname && master->rendertype & R_HYDRO) {
    fprintf(fp, "  {\"%s\", \"%s\", &%s, %s_obc},\n", name, desc, param_name, param_name);
  }
  fprintf(fp, "  {\"NULL\", \"NULL\", NULL, NULL}\n");
  fprintf(fp, "};\n\n");

  /* Existing ecology process function for gateway function         */
  for (n = 0; strcmp(eco_proc[n].name, "NULL") != 0; ++n) {
    if (strcmp(eco_proc[n].name, name) != 0 && !existf && !remp[n]) {
      render_eco_types(fp, eco_proc[n].name);
      fprintf(fp, "\n");
    }
  }
  /* New ecology process function for gateway function              */
  if (hasname && master->rendertype & R_ECO) {
    render_eco_types(fp, name);
    fprintf(fp, "\n");
  }
  /* Ecology process gateway function                               */
  render_eco_gateway(master, fp, name, desc, remp);
  fprintf(fp, "\n");

  i_free_1d(reme);
  i_free_1d(remp);
  fclose(fp);
}

/* END write_rendered_globals()                                     */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Sets the hydro parameters from a global rendered configuration.  */
/*------------------------------------------------------------------*/
int get_rendered_hydroparams(char *name, parameters_t *params)
{
  int n, m;
  parameters_t *in;
  open_bdrys_t **open;

  for (n = 0; strcmp(hydro_param[n].name, "NULL") != 0; ++n) {
    if (strcmp(hydro_param[n].name, name) == 0) {
      hd_warn("get_rendered_hydroparams: Using parameter set PARAMETERS_%s\n", name);
      in = hydro_param[n].param;
      open = &hydro_param[n].obc;
      copy_hydroparameter_info(params, in);
      /*
      for (m = 0; m < in->nobc; m++) {
	if (params->open == NULL) params->open[m] = OBC_alloc();
	params->open[m]->ntr = in->ntr;
	copy_OBC_conds(params->open[n], in->open[m], m, NULL);
      }
      */
      return(0);
    }
  }
  return(1);
}

/* END get_rendered_hydroparams()                                   */
/*------------------------------------------------------------------*/


#if defined(HAVE_ECOLOGY_MODULE)

/*------------------------------------------------------------------*/
/* Sets the BGC parameters from a global rendered configuration.    */
/*------------------------------------------------------------------*/
int get_rendered_ecoparams(char *name, parameter_info *parameters[], int *nprm)
{
  int n, m;
  int nparam;
  parameter_glob_t *param;

  for (n = 0; strcmp(eco_param[n].name, "NULL") != 0; ++n) {
    if (strcmp(eco_param[n].name, name) == 0) {
      parameter_info *in = *parameters;
      hd_warn("get_rendered_ecoparams: Using parameter set PARAMETERS_%s\n", name);
      nparam = eco_param[n].nparam;
      param = eco_param[n].param;
      for (m = 0; m < nparam; m++) {
	copy_ecoparameter_info(&in[m], &param[m]);
      }
      *nprm = nparam;
      return(0);
    }
  }
  return(1);
}

/* END get_rendered_ecoparams()                                     */
/*------------------------------------------------------------------*/

#endif

#if defined(HAVE_SEDIMENT_MODULE)

/*------------------------------------------------------------------*/
/* Sets the sediment parameters from a global rendered              */
/* configuration.                                                   */
/*------------------------------------------------------------------*/
int get_rendered_sedparams(char *name, sediment_t *sediment, int mode)
{
  sed_params_t *param = sediment->msparam;
  sed_params_t *gparam;
  int n, m;

  for (n = 0; strcmp(sed_param[n].name, "NULL") != 0; ++n) {
    if (strcmp(sed_param[n].name, name) == 0) {
      hd_warn("get_rendered_sedparams: Using parameter set PARAMETERS_%s\n", name);
      if (mode == 0) return(0);
      gparam = sed_param[n].param;
      copy_sedparameter_info(param, gparam);
      return(0);
    }
  }
  return(1);
}
/* END get_rendered_sedparams()                                     */
/*------------------------------------------------------------------*/

#endif

/*------------------------------------------------------------------*/
/* Writes a module header.                                          */
/*------------------------------------------------------------------*/
void write_module_header(FILE *fp, char *name, char *desc)
{
  long t;

  /* Write the autotracer list                                       */
  time(&t);
  fprintf(fp,"/*\n");
  fprintf(fp," *\n");
  fprintf(fp," *  ENVIRONMENTAL MODELLING SUITE (EMS)\n");
  fprintf(fp," *\n");  
  fprintf(fp," *  File: %s\n", name);
  fprintf(fp," *\n");
  fprintf(fp," *  Description:\n");
  fprintf(fp," *  %s\n", desc);
  fprintf(fp," *  \n");
  fprintf(fp," *  Copyright:\n");
  fprintf(fp," *  Copyright (c) 2018. Commonwealth Scientific and Industrial\n");
  fprintf(fp," *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights\n");
  fprintf(fp," *  reserved. See the license file for disclaimer and full\n");
  fprintf(fp," *  use/redistribution conditions.\n");
  fprintf(fp," *  \n");
  fprintf(fp," *  $Id: %s %s %s her127 $\n", name, version, ctime(&t));
  fprintf(fp," *\n");
  fprintf(fp," */\n\n");

}

/* END write_module_header()                                        */
/*------------------------------------------------------------------*/

#if defined(HAVE_ECOLOGY_MODULE)

/*------------------------------------------------------------------*/
/* Write eco tracer private data                                    */
/*------------------------------------------------------------------*/
void write_private_eco(FILE *fp, trinfo_priv_eco_t *data, int mode)
{
  fprintf(fp, "  {\n");
  fprintf(fp, "    .type = %d,\n", data->type);
  fprintf(fp, "    .name = \"%s\",\n", data->name);
  fprintf(fp, "    .trname = \"%s\",\n", data->trname);
  fprintf(fp, "    .obc = %d,\n", data->obc);
  fprintf(fp, "    .flag = %d,\n", data->flag);
  if (strlen(data->optical_file))
    fprintf(fp, "    .optical_file = %s,\n", data->optical_file);
  if (strlen(data->absorp_name))
    fprintf(fp, "    .absorp_name = %s,\n", data->absorp_name);
  if (strlen(data->scatter_name))
    fprintf(fp, "    .scatter_name = %s,\n", data->scatter_name);
  if (strlen(data->backscat_name))
    fprintf(fp, "    .backscat_name = %s,\n", data->backscat_name);
  if (strlen(data->abtance_name))
    fprintf(fp, "    .abtance_name = %s,\n", data->abtance_name);
  if (strlen(data->refltce_name))
    fprintf(fp, "    .refltce_name = %s,\n", data->refltce_name);
  if (strlen(data->trnmiss_name))
    fprintf(fp, "    .trnmiss_name = %s,\n", data->trnmiss_name);
  if (strlen(data->benreflt_name))
    fprintf(fp, "    .benreflt_name = %s,\n", data->benreflt_name);
  if (mode)
    fprintf(fp, "  },\n");
  else
    fprintf(fp, "  }\n");
}

/* END write_private_eco()                                          */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Eco process selector function.                                   */
/*------------------------------------------------------------------*/
void render_eco_process(master_t *master, void *ei, FILE *fp, char *name,
			int *wc, int *sed, int *epi)
{
  ecology *e = (ecology *)ei;
  int type;

  for (type = 0; type < N_EPROCESS_TYPES; ++type) {
    int n = e->npr[type];
    int i, mode;
    if (n > 0) {
      if (type == 1)
	fprintf(fp, "const char *WC_PROCESS_%s[] = {\n", name);
      if (type == 2)
	fprintf(fp, "const char *SED_PROCESS_%s[] = {\n", name);
      if (type == 3)
	fprintf(fp, "const char *EPI_PROCESS_%s[] = {\n", name);
      for (i = 0; i < n; ++i) {
	fprintf(fp, "    \"%s\"", e->processes[type][i]->fullname);
	if (i == n-1)
	  fprintf(fp, "\n");
	else
	  fprintf(fp, ",\n");
      }
      fprintf(fp, "  };\n");
    }
    if (type == 1) {
      fprintf(fp, "const int NUM_WC_PROCESS_%s =\n", name);
      fprintf(fp, "  (int)(sizeof(WC_PROCESS_%s)/sizeof(char *));\n", name);
      *wc = n;
    }
    if (type == 2) {
      fprintf(fp, "const int NUM_SED_PROCESS_%s =\n", name);
      fprintf(fp, "  (int)(sizeof(SED_PROCESS_%s)/sizeof(char *));\n", name);
      *sed = n;
    }
    if (type == 3) {
      fprintf(fp, "const int NUM_EPI_PROCESS_%s =\n", name);
      fprintf(fp, "  (int)(sizeof(EPI_PROCESS_%s)/sizeof(char *));\n", name);
      *epi = n;
    }
  }
}

/* END render_eco_process()                                         */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Writes the eco parameters to the config header.                  */
/*------------------------------------------------------------------*/
int render_eco_parameters(master_t *master, void *ei, FILE *fp, char *name)

{
  ecology *e = (ecology *)ei;
  parameter_info *prm = e->pinfo;
  int nparams = e->nprm;
  int n, m;

  fprintf(fp, "#if defined(HAVE_ECOLOGY_MODULE)\n");
  fprintf(fp, "\nint NPARAMETERS_%s = %d;\n", name, nparams);
  fprintf(fp, "parameter_glob_t PARAMETERS_%s[] = {\n", name);
  for (n = 0; n < nparams; n++) {
    fprintf(fp, "  {\n");
    fprintf(fp, "    .name = \"%s\",\n", prm[n].name);
    fprintf(fp, "    .desc = \"%s\",\n", prm[n].desc);
    fprintf(fp, "    .units = \"%s\",\n", prm[n].units);
    fprintf(fp, "    .sym = \"%s\",\n", prm[n].sym);
    if (strlen(prm[n].stringvalue))
      fprintf(fp, "    .stringvalue = \"%s\",\n", prm[n].stringvalue);
    fprintf(fp, "    .num_values = %d,\n", prm[n].num_values);
    if (prm[n].num_values) {
      for (m = 0; m < prm[n].num_values; m++)
	fprintf(fp, "    .value[%d] = %e,\n", m, prm[n].value[m]);
    }
    fprintf(fp, "    .stderr = %f,\n", prm[n].stderr);
    fprintf(fp, "    .ref = \"%s\",\n", prm[n].ref);
    fprintf(fp, "    .index = %d,\n", n);
    if (n == nparams-1)
      fprintf(fp, "  }\n");
    else
    fprintf(fp, "  },\n");
  }
  fprintf(fp, "};\n");
  fprintf(fp, "#endif");
  return(nparams);
}

/* END render_eco_parameters()                                      */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Copies ecology parameters to a new data structure                */
/*------------------------------------------------------------------*/
void copy_ecoparameter_info(parameter_info *out, parameter_glob_t *in)
{
  int n;
  out->index = in->index;
  out->name = strdup(in->name);
  out->desc = strdup(in->desc);
  out->sym = strdup(in->sym);
  out->units = strdup(in->units);
  out->stringvalue = strdup(in->stringvalue);
  out->ref = strdup(in->ref);
  out->num_values = in->num_values;
  out->stderr = in->stderr;
  if (out->value) d_free_1d(out->value);
  out->value = d_alloc_1d(out->num_values);
  for (n = 0; n < in->num_values; n++)
    out->value[n] = in->value[n];
}

/* END copy_ecoparameter_info()                                     */
/*------------------------------------------------------------------*/

#endif

/*------------------------------------------------------------------*/
/* Writes the gateway function to select the process list.          */
/*------------------------------------------------------------------*/
void render_eco_gateway(master_t *master, FILE *fp, char *name, char *desc, int *remp)
{
  int n;

  fprintf(fp, "/*\n");
  fprintf(fp, " * This function is called from create_processes and serves as a\n");
  fprintf(fp, " * gateway for the different default processes types.\n");
  fprintf(fp, " */\n");
  fprintf(fp, "int get_rendered_processes(char *name, int type, const char **procs[], int *nprocs)\n");
  fprintf(fp, "{\n");
  fprintf(fp, "  int i,n;\n\n");
  fprintf(fp, "  struct {\n");
  fprintf(fp, "    char *name;\n");
  fprintf(fp, "    char *description;\n");
  fprintf(fp, "    void (*init) (int type, const char **procs[], int *nprocs);\n");
  fprintf(fp, "  } param_list[] = {\n");
  for (n = 0; strcmp(eco_proc[n].name, "NULL") != 0; ++n) {
    if (strcmp(eco_proc[n].name, name) != 0 && !remp[n])
      fprintf(fp, "    {\"%s\",\"%s\", eco_processes_%s},\n", eco_proc[n].name, eco_proc[n].desc, eco_proc[n].name);
  }
  if (strlen(name) && master->rendertype & R_ECO)
    fprintf(fp, "    {\"%s\",\"%s\", eco_processes_%s},\n", name, desc, name);
  fprintf(fp, "    {\"NULL\", \"NULL\", NULL}\n");
  fprintf(fp, "  };\n");
  fprintf(fp, "  void (*init) (int type, const char **procs[], int *nprocs) = NULL;\n\n");
  fprintf(fp, "  for (i = 0; (init == NULL) && param_list[i].name; ++i) {\n");
  fprintf(fp, "    if (strcasecmp(name, param_list[i].name) == 0) {\n");
  fprintf(fp, "      init = param_list[i].init;\n");
  fprintf(fp, "    }\n");
  fprintf(fp, "  }\n");
  fprintf(fp, "  if (init != NULL) {\n");
  fprintf(fp, "    init(type, procs, nprocs);\n");
  fprintf(fp, "  } else\n");
  fprintf(fp, "    return(1);\n");
  fprintf(fp, "  return(0);\n");
  fprintf(fp, "}\n\n");
}

/* END render_eco_gateway()                                         */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Writes an ecology process function for gateway function          */
/*------------------------------------------------------------------*/
void render_eco_types(FILE *fp, char *name)
{

  fprintf(fp, "static void eco_processes_%s(int type, const char **procs[], int *nprocs)\n", name);
  fprintf(fp, "{\n");
  fprintf(fp, "  /* Key off process type */\n");
  fprintf(fp, "  switch(type) {\n");
  fprintf(fp, "  case PT_WC:\n");
  fprintf(fp, "    *procs  = WC_PROCESS_%s;\n", name);
  fprintf(fp, "    *nprocs = NUM_WC_PROCESS_%s;\n", name);
  fprintf(fp, "    return;\n");
  fprintf(fp, "  case PT_EPI:\n");
  fprintf(fp, "    *procs  = EPI_PROCESS_%s;\n", name);
  fprintf(fp, "    *nprocs = NUM_EPI_PROCESS_%s;\n", name);
  fprintf(fp, "    return;\n");
  fprintf(fp, "  case PT_SED:\n");
  fprintf(fp, "    *procs  = SED_PROCESS_%s;\n", name);
  fprintf(fp, "    *nprocs = NUM_SED_PROCESS_%s\n;", name);
  fprintf(fp, "    return;\n");
  fprintf(fp, "  case PT_COL:\n");
  fprintf(fp, "    *procs = NULL;\n");
  fprintf(fp, "    *nprocs = 0;\n");
  fprintf(fp, "    return;\n");
  fprintf(fp, "  default:\n");
  fprintf(fp, "    hd_quit(\"eco_processes_%s: Unknown Process Type \");\n", name);
  fprintf(fp, "  }\n");
  fprintf(fp, "}\n");
}

/* END render_eco_types()                                           */
/*------------------------------------------------------------------*/

#if defined(HAVE_SEDIMENT_MODULE)
#include "sediments.h"

/*------------------------------------------------------------------*/
/* Copies sediment parameters to a new data structure               */
/*------------------------------------------------------------------*/
void copy_sedparameter_info(sed_params_t *out, sed_params_t *in)
{
  int n;
  out->geomorph = in->geomorph;
  out->consolidate = in->consolidate;
  out->minpor_wc = in->minpor_wc;
  out->minpor_sed = in->minpor_sed;
  out->finpor_sed = in->finpor_sed;
  out->mindepth_sedi = in->mindepth_sedi;
  out->consolrate = in->consolrate;
  out->cssmode = in->cssmode;
  out->css = in->css;
  out->css_dep = in->css_dep;
  out->flocmode = in->flocmode;
  out->flocprm1 = in->flocprm1;
  out->flocprm2 = in->flocprm2;
  out->bbl_nonlinear = in->bbl_nonlinear;
  out->calc_ripples = in->calc_ripples;
  out->physriph = in->physriph;
  out->physripl = in->physripl;
  out->bioriph = in->bioriph;
  out->bioripl = in->bioripl;
  out->biodens = in->biodens;
  out->maxbiodepth = in->maxbiodepth;
  out->bi_dissol_kz = in->bi_dissol_kz;
  out->bi_dissol_kz_i = in->bi_dissol_kz_i;
  out->bt_partic_kz = in->bt_partic_kz;
  out->bt_partic_kz_i = in->bt_partic_kz_i;
  out->css_scale = in->css_scale;
  out->biosedprofile = in->biosedprofile;
  out->max_thick_sed = in->max_thick_sed;
  out->quad_bfc = in->quad_bfc;
  out->uf = in->uf;
  out->erflux_scale = in->erflux_scale;
  out->bbl_nonlinear = in->bbl_nonlinear;
  out->hindered_svel_patch = in->hindered_svel_patch;
  out->hindered_svel = in->hindered_svel;
  out->reef_scale_depth = in->reef_scale_depth;
  out->calc_ripples = in->calc_ripples;
  out->verbose_sed = in->verbose_sed;
}

/* END copy_sedparameter_info()                                     */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Write sediment tracer private data                               */
/*------------------------------------------------------------------*/
void write_private_sed(FILE *fp, trinfo_priv_sed_t *data, int mode)
{
  fprintf(fp, "  {\n");
  fprintf(fp, "    .name = \"%s\",\n", data->name);
  fprintf(fp, "    .trname = \"%s\",\n", data->trname);
  fprintf(fp, "    .type = %d,\n", data->type);
  fprintf(fp, "    .cohesive = %d,\n", data->cohesive);
  fprintf(fp, "    .calcvol = %d,\n", data->calcvol);
  fprintf(fp, "    .floc = %d,\n", data->floc);
  fprintf(fp, "    .resuspend = %d,\n", data->resuspend);
  fprintf(fp, "    .deposit = %d,\n", data->deposit);
  fprintf(fp, "    .adsorb = %d,\n", data->adsorb);
  fprintf(fp, "    .obc = 0x%x,\n", data->obc);
  fprintf(fp, "    .psize = %f,\n", data->psize);
  fprintf(fp, "    .b_dens = %f,\n", data->b_dens);
  fprintf(fp, "    .i_conc = %f,\n", data->i_conc);
  fprintf(fp, "    .f_conc = %f,\n", data->f_conc);
  fprintf(fp, "    .css_erosion = %f,\n", data->css_erosion);
  fprintf(fp, "    .css_deposition = %f,\n", data->css_deposition);
  fprintf(fp, "    .svel = %f,\n", data->svel);
  if (strlen(data->svel_name))
    fprintf(fp, "    .svel_name = \"%s\",\n", data->svel_name);
  if (strlen(data->carriername))
  fprintf(fp, "    .adsorbkd = %f,\n", data->adsorbkd);
  fprintf(fp, "    .adsorbrate = %f,\n", data->adsorbrate);
  if (strlen(data->dissolvedname))
  fprintf(fp, "    .dissolvedname = \"%s\",\n", data->dissolvedname);
  if (mode)
    fprintf(fp, "  },\n");
  else
    fprintf(fp, "  }\n");
}

/* END write_private_sed()                                          */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Writes the sediment parameters to the config header.             */
/*------------------------------------------------------------------*/
void render_sed_params(FILE *fp, char *name)
{
  win_priv_t *wincon = hd_data->wincon[1];
  sediment_t *sediment = wincon->sediment;
  sed_params_t *param = sediment->msparam;

  fprintf(fp, "#if defined(HAVE_SEDIMENT_MODULE)\n\n");
  fprintf(fp, "sed_params_t SEDPARAMS_%s = {\n", name);
  fprintf(fp, "    .geomorph = %d,\n", param->geomorph);
  fprintf(fp, "    .consolidate = %d,\n", param->consolidate);
  fprintf(fp, "    .minpor_wc = %f,\n", param->minpor_wc);
  fprintf(fp, "    .minpor_sed = %f,\n", param->minpor_sed);
  fprintf(fp, "    .finpor_sed = %f,\n", param->finpor_sed);
  fprintf(fp, "    .mindepth_sedi = %f,\n", param->mindepth_sedi);
  fprintf(fp, "    .consolrate = %f,\n", param->consolrate);
  fprintf(fp, "    .cssmode = %d,\n", param->cssmode);
  fprintf(fp, "    .css = %f,\n", param->css);
  fprintf(fp, "    .css_dep = %e,\n", param->css_dep);
  fprintf(fp, "    .flocmode = %d,\n", param->flocmode);
  fprintf(fp, "    .flocprm1 = %f,\n", param->flocprm1 );
  fprintf(fp, "    .flocprm2 = %f,\n", param->flocprm2 );
  fprintf(fp, "    .bbl_nonlinear = %d,\n", param->bbl_nonlinear );
  fprintf(fp, "    .calc_ripples = %d,\n", param->calc_ripples);
  if (!param->physriph_spv)
    fprintf(fp, "    .physriph = %f,\n", param->physriph);
  if (!param->physripl_spv)
    fprintf(fp, "    .physripl = %f,\n", param->physripl);
  fprintf(fp, "    .bioriph = %f,\n", param->bioriph);
  fprintf(fp, "    .bioripl = %f,\n", param->bioripl);
  if (!param->biodens_spv)
    fprintf(fp, "    .biodens = %f,\n", param->biodens);
  if (!param->maxbiodepth_spv)
    fprintf(fp, "    .maxbiodepth = %f,\n", param->maxbiodepth);
  if (!param->bi_dissol_kz_spv)
    fprintf(fp, "    .bi_dissol_kz = %e,\n", param->bi_dissol_kz);
  if (!param->bi_dissol_kz_i_spv)
    fprintf(fp, "    .bi_dissol_kz_i = %e,\n", param->bi_dissol_kz_i);
  if (!param->bt_partic_kz_spv)
    fprintf(fp, "    .bt_partic_kz = %e,\n", param->bt_partic_kz);
  if (!param->bt_partic_kz_i_spv)
    fprintf(fp, "    .bt_partic_kz_i = %e,\n", param->bt_partic_kz_i);
  if (!param->css_scale_spv)
    fprintf(fp, "    .css_scale = %f,\n", param->css_scale);
  fprintf(fp, "    .biosedprofile = '%c',\n",  param->biosedprofile);
  fprintf(fp, "    .max_thick_sed = %f,\n",  param->max_thick_sed);
  fprintf(fp, "    .quad_bfc = %f,\n", param->quad_bfc);
  fprintf(fp, "    .uf = %f,\n", param->uf);
  fprintf(fp, "    .erflux_scale = %f,\n", param->erflux_scale);
  fprintf(fp, "    .bbl_nonlinear = %d,\n", param->bbl_nonlinear);
  fprintf(fp, "    .hindered_svel_patch = %d,\n", param->hindered_svel_patch);
  fprintf(fp, "    .hindered_svel = %d,\n", param->hindered_svel);
  fprintf(fp, "    .reef_scale_depth = %f,\n", param->reef_scale_depth);
  fprintf(fp, "    .calc_ripples = %d,\n", param->calc_ripples);
  fprintf(fp, "    .verbose_sed = %d\n", param->verbose_sed);
  fprintf(fp, "};\n");
  fprintf(fp, "#endif\n");
}

/* END render_sed_params()                                          */
/*------------------------------------------------------------------*/

#endif

/*------------------------------------------------------------------*/
/* Sets tracer attributes to those in the global configuration      */
/* specification corresponding to the class edo_defs.               */
/*------------------------------------------------------------------*/
int get_rendered_tracers(FILE *fp, 
			 int do_eco, 
			 char *eco_vars,    /* list of vars         */ 
			 char *eco_defs,    /* name of class        */
			 tracer_info_t *trinfo, 
			 int ntr, int tn,
			 int trinfo_type)
{
  char *vars[ECO_MAXNUMARGS], buf[MAXSTRLEN];
  int m, n, nn, i, j;
  tracer_info_t *tr;
#if defined(HAVE_ECOLOGY_MODULE)
  trinfo_priv_eco_t *ecodata;
#endif
#if defined(HAVE_SEDIMENT_MODULE)
  trinfo_priv_sed_t *seddata;
#endif
  int found;

  if (eco_vars == NULL) return(tn);
  if (strlen(eco_vars) == 0) return(tn);
  /* Loop ovar all tracer names in the list of vars. Note: the      */
  /* tracers list was built in ecology_pre_build() based on the     */
  /* list of processes supplied. This should correspond exactly to  */
  /* the tracers in the global list - exit if it can't find the     */
  /* corresponding tracer.                                          */
  /* First find the global list corresponding to 'eco_defs'         */
  for (n = 0; strcmp(eco_glob[n].name, "NULL") != 0; ++n) {
    if (strcmp(eco_glob[n].name, eco_defs) == 0) {
      hd_warn("get_rendered_tracers: Found global dataset %s.h for tracer attributes.\n",
	      eco_defs);
      tr = eco_glob[n].trlist;
#if defined(HAVE_ECOLOGY_MODULE)
      ecodata = eco_glob[n].eco_priv;
#endif
#if defined(HAVE_SEDIMENT_MODULE)
      seddata = eco_glob[n].sed_priv;
#endif
      /* Set the attributes in the list to the global spec          */
      strcpy(buf, eco_vars);
      m = parseline(buf, vars, ECO_MAXNUMARGS);
      for (nn = 0; nn < m; nn++) {
	found = 0;
	for (i = 0; i < eco_glob[n].ntr; i++) {
	  if (strcmp(vars[nn], tr[i].name) == 0) {
	    found = 1;
	    if (trinfo_type & tr[i].type && (tracer_find_index(vars[nn], ntr, trinfo)) < 0) {
	      tracer_info_t *trg = &tr[i];
	      tracer_info_t *tre = &trinfo[tn];
	      /* Find the ecology private data and copy             */
#if defined(HAVE_ECOLOGY_MODULE)
	      for (j = 0; j < eco_glob[n].nep; j++) {
		if (strcmp(vars[nn], ecodata[j].trname) == 0) {
		  trg->private_data[TR_PRV_DATA_ECO] = &ecodata[j];
		  trg->private_data_copy[TR_PRV_DATA_ECO] = private_data_copy_eco;
		  break;
		} else
		  trg->private_data[TR_PRV_DATA_ECO] = NULL;
	      }
#endif
	      /* Find the sediment private data and copy            */
#if defined(HAVE_SEDIMENT_MODULE)
	      for (j = 0; j < eco_glob[n].nsp; j++) {
		if (strcmp(vars[nn], seddata[j].trname) == 0) {
		  trg->private_data[TR_PRV_DATA_SED] = &seddata[j];
		  trg->private_data_copy[TR_PRV_DATA_SED] = private_data_copy_sed;
		  break;
		} else
		  trg->private_data[TR_PRV_DATA_SED] = NULL;
	      }
#endif
	      /* Copy the tracer attributes                         */
	      tracer_copy(tre, trg);
	      /* Find the sediment private data and copy            */
#if defined(HAVE_SEDIMENT_MODULE)
	      sed_read_tr_atts(&trinfo[tn], fp, vars[nn]); // includes data
#endif
#if defined(HAVE_ECOLOGY_MODULE)
	      eco_read_tr_atts(&trinfo[tn], fp, vars[nn]);
#endif
	      trinfo[tn].m = -1; // does not exist in the prm-file
	      tracer_re_read(&trinfo[tn], fp, trinfo_type);
	      trinfo[tn].n = tn;
	      tn++;
	    }
	  }
	}
	if (!found) hd_quit("eco_global_tracers: Can't find tracer %s in global list %s.h\n",
			    vars[nn], eco_glob[n].name);
      }
    }
  }
  return(tn);
}

/* END get_rendered_tracers()                                       */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Writes the hydro parameters to the config header.                */
/*------------------------------------------------------------------*/
void render_hydro_parameters(master_t *master, FILE *fp, char *name)
{
  parameters_t *params = master->params;
  char buf[MAXSTRLEN];
  int n, m;
  int existf = 0;
  int dof = 0;

  if (dof) {
  /*----------------------------------------------------------------*/
  /* Check this name doesn't exist                                  */
  for (n = 0; strcmp(hydro_param[n].name, "NULL") != 0; ++n) {
    if (strcmp(hydro_param[n].name, name) == 0) {
      hd_warn("render_hydro_parameters: hydro name '%s' already rendered!\n", name);
      existf = 1;
    }
  }

  /*----------------------------------------------------------------*/
  /* Open the configuration header                                  */
  if (strlen(master->renderpath))
    if (endswith(master->renderpath, "/"))
      sprintf(buf, "%s%s.h", master->renderpath, name);
    else
      sprintf(buf, "%s/%s.h", master->renderpath, name);
  else
    sprintf(buf, "%s.h", name);
  fp = fopen(buf, "w");
  }

  /*----------------------------------------------------------------*/
  /* Write the configuration to a header file.                      */
  /*
  sprintf(buf, "Hydrodynamic configuration for %s: %s\n", name, master->renderdesc);
  write_module_header(fp, name, buf);
  fprintf(fp,"#include \"ems.h\"\n\n");
  */
  fprintf(fp, "\nparameters_t HYDRO_%s = {\n", name);
  fprintf(fp, "  .codeheader = \"%s\",\n", params->codeheader);
  fprintf(fp, "  .parameterheader = \"%s\",\n", params->parameterheader);
  fprintf(fp, "  .reference = \"%s\",\n", params->reference);
  fprintf(fp, "  .notes = \"%s\",\n", params->notes);
  fprintf(fp, "  .trl = \"%s\",\n", params->trl);
  fprintf(fp, "  .gridtype = \"%s\",\n", params->gridtype);
  fprintf(fp, "  .gridcode = 0x%x,\n", params->gridcode);
  fprintf(fp, "  .projection = \"%s\",\n", params->projection);
  fprintf(fp, "  .us_type = 0x%x,\n", params->us_type);
  fprintf(fp, "  .opath = \"%s\",\n", params->opath);
  fprintf(fp, "  .trkey = \"%s\",\n", params->trkey);
  fprintf(fp, "  .sequence = \"%s\",\n", params->sequence);
  fprintf(fp, "  .bdrypath = \"%s\",\n", params->bdrypath);
  fprintf(fp, "  .tracerdata = \"%s\",\n", params->tracerdata);
  fprintf(fp, "  .rivldir = \"%s\",\n", params->rivldir);
  fprintf(fp, "  .trans_dt = \"%s\",\n", params->trans_dt);
  fprintf(fp, "  .runnoc = \"%s\",\n", params->runnoc);
  fprintf(fp, "  .runno = %f,\n", params->runno);
  fprintf(fp, "  .runcode = \"%s\",\n", params->runcode);
  fprintf(fp, "  .rev = \"%s\",\n", params->rev);
  fprintf(fp, "  .history = 0x%x,\n", params->history);
  fprintf(fp, "  .start_time = \"%s\",\n", params->start_time);
  fprintf(fp, "  .stop_time = \"%s\",\n", params->stop_time);
  fprintf(fp, "  .lenunit = \"%s\",\n", params->lenunit);
  fprintf(fp, "  .timeunit = \"%s\",\n", params->timeunit);
  fprintf(fp, "  .output_tunit = \"%s\",\n", params->output_tunit);
  fprintf(fp, "  .t = %f,\n", params->t);
  fprintf(fp, "  .ambpress = %f,\n", params->ambpress);
  fprintf(fp, "  .air_dens = %f,\n", params->air_dens);
  fprintf(fp, "  .spec_heat = %f,\n", params->spec_heat);
  fprintf(fp, "  .g = %f,\n", params->g);
  fprintf(fp, "  .grid_dt = %f,\n", params->grid_dt);
  fprintf(fp, "  .iratio = %d,\n", params->iratio);
  fprintf(fp, "  .tratio = %f,\n", params->tratio);
  fprintf(fp, "  .trsplit = %d,\n", params->trsplit);
  fprintf(fp, "  .momsc = 0x%x,\n", params->momsc);
  fprintf(fp, "  .trasc = 0x%x,\n", params->trasc);
  fprintf(fp, "  .ultimate = %d,\n", params->ultimate); \
  fprintf(fp, "  .osl = 0x%x,\n", params->osl);
  fprintf(fp, "  .rkstage = %d,\n", params->rkstage);
  fprintf(fp, "  .smag = \"%s\",\n", params->smag);
  fprintf(fp, "  .smagorinsky = %e,\n", params->smagorinsky);
  fprintf(fp, "  .u1vh = %e,\n", params->u1vh);
  fprintf(fp, "  .u1kh = %e,\n", params->u1kh);
  fprintf(fp, "  .sue1 = %e,\n", params->sue1);
  fprintf(fp, "  .kue1 = %e,\n", params->kue1);
  fprintf(fp, "  .bsue1 = %e,\n", params->bsue1);
  fprintf(fp, "  .bkue1 = %e,\n", params->bkue1);
  fprintf(fp, "  .smag_smooth = %d,\n", params->smag_smooth);
  fprintf(fp, "  .diff_scale = 0x%x,\n", params->diff_scale);
  fprintf(fp, "  .visc_method = 0x%x,\n", params->visc_method);
  fprintf(fp, "  .visc_fact = %e,\n", params->visc_fact);
  fprintf(fp, "  .mixsc = \"%s\",\n", params->mixsc);
  fprintf(fp, "  .s_func = \"%s\",\n", params->s_func);
  fprintf(fp, "  .min_tke = %e,\n", params->min_tke);
  fprintf(fp, "  .min_diss = %e,\n", params->min_diss);
  fprintf(fp, "  .vz0 = %e,\n", params->vz0);
  fprintf(fp, "  .kz0 = %e,\n", params->kz0);
  fprintf(fp, "  .zs = %e,\n", params->zs);
  fprintf(fp, "  .stab = %d,\n", params->stab);
  /*
  fprintf(fp, "  .atr = %d,\n", params->atr);
  fprintf(fp, "  .ntr = %d,\n", params->ntr);
  fprintf(fp, "  .ntrS = %d,\n", params->ntrS);
  fprintf(fp, "  .nsed = %d,\n", params->nsed);
  */
  fprintf(fp, "  .sednz = %d,\n", params->sednz);
  fprintf(fp, "  .autotrpath = \"%s\",\n", params->autotrpath);
  fprintf(fp, "  .cfl = 0x%x,\n", params->cfl);
  fprintf(fp, "  .cfl_dt = \"%s\",\n", params->cfl_dt);
  fprintf(fp, "  .mixlayer = 0x%x,\n", params->mixlayer);
  fprintf(fp, "  .show_layers = %d,\n", params->show_layers);
  fprintf(fp, "  .means = 0x%x,\n", params->means);
  fprintf(fp, "  .means_os = \"%s\",\n", params->means_os);
  fprintf(fp, "  .bathystats = \"%s\",\n", params->bathystats);
  fprintf(fp, "  .particles = \"%s\",\n", params->particles);
  fprintf(fp, "  .addquad = \"%s\",\n", params->addquad);
  fprintf(fp, "  .vorticity = 0x%x,\n", params->vorticity);
  fprintf(fp, "  .numbers = 0x%x,\n", params->numbers);
  fprintf(fp, "  .numbers1 = 0x%x,\n", params->numbers1);
  fprintf(fp, "  .u1_f = 0x%x,\n", params->u1_f);
  fprintf(fp, "  .u1av_f = 0x%x,\n", params->u1av_f);
  fprintf(fp, "  .save_force = 0x%x,\n", params->save_force);
  fprintf(fp, "  .lnm = %f,\n", params->lnm);
  fprintf(fp, "  .trflux = \"%s\",\n", params->trflux);
  fprintf(fp, "  .trfd1 = %d,\n", params->trfd1);
  fprintf(fp, "  .trfd2 = %d,\n", params->trfd2);
  fprintf(fp, "  .trperc = \"%s\",\n", params->trperc);
  fprintf(fp, "  .trpercr = \"%s\",\n", params->trpercr);
  fprintf(fp, "  .trflsh = %d,\n", params->trflsh);
  fprintf(fp, "  .trage = \"%s\",\n", params->trage);
  fprintf(fp, "  .trtend = \"%s\",\n", params->trtend);
  fprintf(fp, "  .ndhw = %d,\n", params->ndhw);
  fprintf(fp, "  .tendf = %d,\n", params->tendf);
  fprintf(fp, "  .thin_merge = %d,\n", params->thin_merge);
  fprintf(fp, "  .sigma = %d,\n", params->sigma);
  fprintf(fp, "  .bathyfill = 0x%x,\n", params->bathyfill);
  fprintf(fp, "  .bvel = %e,\n", params->bvel);
  fprintf(fp, "  .nonlinear = %d,\n", params->nonlinear);
  fprintf(fp, "  .calc_dens = %d,\n", params->calc_dens);
  fprintf(fp, "  .mode2d = %d,\n", params->mode2d);
  fprintf(fp, "  .bmin = %f,\n", params->bmin);
  fprintf(fp, "  .bmax = %f,\n", params->bmax);
  fprintf(fp, "  .mct = \"%s\",\n", params->mct);
  fprintf(fp, "  .smooth = %d,\n", params->smooth);
  fprintf(fp, "  .maxgrad = %f,\n", params->maxgrad);
  fprintf(fp, "  .maxdiff = \"%s\",\n", params->maxdiff);
  fprintf(fp, "  .slipprm = %f,\n", params->slipprm);
  fprintf(fp, "  .nwindows = %d,\n", params->nwindows);
  fprintf(fp, "  .win_reset = %d,\n", params->win_reset);
  fprintf(fp, "  .win_type = 0x%x,\n", params->win_type);
  fprintf(fp, "  .win_block = 0x%x,\n", params->win_block);
  fprintf(fp, "  .win_file = \"%s\",\n", params->win_file);
  for (n = 0; n < params->nwindows; n++) {
    if (params->win_size != NULL) {
      fprintf(fp, "  .win_size[%d] = %f,\n", n, params->win_size[n]);
      fprintf(fp, "  .nwn[n] = %d,\n", params->nwn[n]);
      if (params->nwn != NULL) {
	for (m = 0; m < params->nwn[n]; m++) {
	  fprintf(fp, "  .wnz[%d][%d] = %d,\n", n, m, params->wnx[n][m]);
	  fprintf(fp, "  .wnz[%d][%d] = %d,\n", n, m, params->wny[n][m]);
	}
      }
    }
  }
  fprintf(fp, "  .show_win = %d,\n", params->show_win);
  fprintf(fp, "  .restart_dt = %f,\n", params->restart_dt);
  fprintf(fp, "  .da = %d,\n", params->da);
  fprintf(fp, "  .da_dt = %f,\n", params->da_dt);
  fprintf(fp, "  .da_fcst_dt = %f,\n", params->da_fcst_dt);
  fprintf(fp, "  .prex = %d,\n", params->prex);
  fprintf(fp, "  .ndf = %d,\n", params->ndf);
  fprintf(fp, "  .riverflow = %d,\n", params->riverflow);
  fprintf(fp, "  .meshinfo = %d,\n", params->meshinfo);
  fprintf(fp, "  .swr_type = 0x%x,\n", params->swr_type);
  for (n = 0; n < 6; n++)
    fprintf(fp, "  .swr_ens[%d] = %f,\n", n, params->swr_ens[n]);
  fprintf(fp, "  .u1vhc = \"%s\",\n", params->u1vhc);
  fprintf(fp, "  .u1khc = \"%s\",\n", params->u1khc);
  fprintf(fp, "  .restart_name = \"%s\",\n", params->restart_name);
  fprintf(fp, "  .wind_file = \"%s\",\n", params->wind_file);
  fprintf(fp, "  .dlv0 = %f,\n", params->dlv0);
  fprintf(fp, "  .dlv1 = %f,\n", params->dlv1);
  fprintf(fp, "  .dlc0 = %f,\n", params->dlc0);
  fprintf(fp, "  .dlc1 = %f,\n", params->dlc1);
  fprintf(fp, "  .wind_type = 0x%x,\n", params->wind_type);
  fprintf(fp, "  .stress_fn = 0x%x,\n", params->stress_fn);
  fprintf(fp, "  .wind_scale = %f,\n", params->wind_scale);
  fprintf(fp, "  .wind = \"%s\",\n", params->wind);
  fprintf(fp, "  .bdry_file = \"%s\",\n", params->bdry_file);
  fprintf(fp, "  .ptinname = \"%s\",\n", params->ptinname);
  fprintf(fp, "  .patm = \"%s\",\n", params->patm);
  fprintf(fp, "  .precip = \"%s\",\n", params->precip);
  fprintf(fp, "  .evap = \"%s\",\n", params->evap);
  fprintf(fp, "  .airtemp = \"%s\",\n", params->airtemp);
  fprintf(fp, "  .rh = \"%s\",\n", params->rh);
  fprintf(fp, "  .cloud = \"%s\",\n", params->cloud);
  fprintf(fp, "  .swr = \"%s\",\n", params->swr);
  fprintf(fp, "  .light = \"%s\",\n", params->light);
  fprintf(fp, "  .wetb = \"%s\",\n", params->wetb);
  fprintf(fp, "  .hftemp = \"%s\",\n", params->hftemp);
  fprintf(fp, "  .hf = \"%s\",\n", params->hf);
  fprintf(fp, "  .swr_babs = \"%s\",\n", params->swr_babs);
  fprintf(fp, "  .swr_attn = \"%s\",\n", params->swr_attn);
  fprintf(fp, "  .swr_attn1 = \"%s\",\n", params->swr_attn1);
  fprintf(fp, "  .swr_tran = \"%s\",\n", params->swr_tran);
  fprintf(fp, "  .swr_regions = \"%s\",\n", params->swr_regions);
  fprintf(fp, "  .swreg_dt = %f,\n", params->swreg_dt);
  fprintf(fp, "  .swr_data = \"%s\",\n", params->swr_data);
  fprintf(fp, "  .densname = \"%s\",\n", params->densname);
  fprintf(fp, "  .regions = \"%s\",\n", params->regions);
  fprintf(fp, "  .region_dt = \"%s\",\n", params->region_dt);
  fprintf(fp, "  .region_vars = \"%s\",\n", params->region_vars);
  fprintf(fp, "  .region_mode = \"%s\",\n", params->region_mode);
  fprintf(fp, "  .i_rule = \"%s\",\n", params->i_rule);
  fprintf(fp, "  .cookiecut = \"%s\",\n", params->cookiecut);
  fprintf(fp, "  .imp2df = \"%s\",\n", params->imp2df);
  fprintf(fp, "  .imp3df = \"%s\",\n", params->imp3df);
  fprintf(fp, "  .do_pt = %d,\n", params->do_pt);
  fprintf(fp, "  .do_lag = %d,\n", params->do_lag);
  fprintf(fp, "  .dp_mode = \"%s\",\n", params->dp_mode);
  fprintf(fp, "  .waves = %d,\n", params->waves);
  fprintf(fp, "  .decorr = %d,\n", params->decorr);
  fprintf(fp, "  .decf = 0x%x,\n", params->decf);
  fprintf(fp, "  .monotr = \"%s\",\n", params->monotr);
  fprintf(fp, "  .monomn = %f,\n", params->monomn);
  fprintf(fp, "  .monomx = %f,\n", params->monomx);
  fprintf(fp, "  .orbital = %d,\n", params->orbital);
  fprintf(fp, "  .rampf = 0x%x,\n", params->rampf);
  fprintf(fp, "  .rampstart = %f,\n", params->rampstart);
  fprintf(fp, "  .rampend = %f,\n", params->rampend);
  fprintf(fp, "  .hmin = %f,\n", params->hmin);
  fprintf(fp, "  .uf = %f,\n", params->uf);
  fprintf(fp, "  .quad_bfc = \"%s\",\n", params->quad_bfc);
  fprintf(fp, "  .totals = %d,\n", params->totals);
  fprintf(fp, "  .roammode = 0x%x,\n", params->roammode);
  fprintf(fp, "  .robust = %d,\n", params->robust);
  fprintf(fp, "  .speed = %d,\n", params->speed);
  fprintf(fp, "  .fatal = 0x%x,\n", params->fatal);
  fprintf(fp, "  .avhrr = %d,\n", params->avhrr);
  fprintf(fp, "  .ghrsst = %d,\n", params->ghrsst);
  fprintf(fp, "  .exmapf = %d,\n", params->exmapf);
  fprintf(fp, "  .heatflux = 0x%x,\n", params->heatflux);
  fprintf(fp, "  .saltflux = 0x%x,\n", params->saltflux);
  fprintf(fp, "  .water_type = 0x%x,\n", params->water_type);
  fprintf(fp, "  .tsfile_caching = %d,\n", params->tsfile_caching);
  fprintf(fp, "  .etamax = %f,\n", params->etamax);
  fprintf(fp, "  .velmax = %f,\n", params->velmax);
  fprintf(fp, "  .velmax2d = %f,\n", params->velmax2d);
  fprintf(fp, "  .etadiff = %f,\n", params->etadiff);
  fprintf(fp, "  .wmax = %f,\n", params->wmax);
  fprintf(fp, "  .trsplit = %d,\n", params->trsplit);
  fprintf(fp, "  .eta_ib = %d,\n", params->eta_ib);
  fprintf(fp, "  .etarlx = 0x%x,\n", params->etarlx);
  fprintf(fp, "  .velrlx = 0x%x,\n", params->velrlx);
  fprintf(fp, "  .albedo = %f,\n", params->albedo);
  fprintf(fp, "  .domom = %d,\n", params->domom);
  fprintf(fp, "  .doroms = %d,\n", params->doroms);
  fprintf(fp, "  .doswan = %d,\n", params->doswan);
  fprintf(fp, "  .zref = %f,\n", params->zref);
  fprintf(fp, "  .neutral = %d,\n", params->neutral);
  fprintf(fp, "  .wave_alpha = %f,\n", params->wave_alpha);
  fprintf(fp, "  .wave_hf = %f,\n", params->wave_hf);
  fprintf(fp, "  .wave_b1 = %f,\n", params->wave_b1);
  fprintf(fp, "  .fillf = 0x%x,\n", params->fillf);
  fprintf(fp, "  .filter = 0x%x,\n", params->filter);
  fprintf(fp, "  .trfilter = 0x%x,\n", params->trfilter);
  fprintf(fp, "  .conserve = 0x%x,\n", params->conserve);
  fprintf(fp, "  .lyear = %d,\n", params->lyear);
  fprintf(fp, "  .do_closure = %d,\n", params->do_closure);
  fprintf(fp, "  .rsalt = %d,\n", params->rsalt);
  fprintf(fp, "  .rtemp = %d,\n", params->rtemp);
  fprintf(fp, "  .smooth_VzKz = %d,\n", params->smooth_VzKz);
  fprintf(fp, "  .albedo_l = %f,\n", params->albedo_l);
  fprintf(fp, "  .trout = %d,\n", params->trout);
  fprintf(fp, "  .porusplate = %d,\n", params->porusplate);
  fprintf(fp, "  .sharp_pyc = %d,\n", params->sharp_pyc);
  fprintf(fp, "  .uscf = %d,\n", params->uscf);
  fprintf(fp, "  .tidef = %d,\n", params->tidef);
  fprintf(fp, "  .tidep = %d,\n", params->tidep);
  fprintf(fp, "  .eqt_alpha = %f,\n", params->eqt_alpha);
  fprintf(fp, "  .eqt_beta = %f,\n", params->eqt_beta);
  fprintf(fp, "  .mlat = %e,\n", params->mlat);
  fprintf(fp, "  .mlon = %e,\n", params->mlon);
  fprintf(fp, "  .nprof = \"%s\",\n", params->nprof);
  fprintf(fp, "  .nprof2d = \"%s\",\n", params->nprof2d);
  fprintf(fp, "  .reef_frac = \"%s\",\n", params->reef_frac);
  fprintf(fp, "  .crf = \"%s\",\n", params->crf);
  fprintf(fp, "  .dbi = %d,\n", params->dbi);
  fprintf(fp, "  .dbj = %d,\n", params->dbj);
  fprintf(fp, "  .dbk = %d,\n", params->dbk);
  fprintf(fp, "  .dbgf = 0x%x,\n", params->dbgf);
  fprintf(fp, "  .dbgtime = %f,\n", params->dbgtime);
  fprintf(fp, "  .momfile = \"%s\",\n", params->momfile);
  fprintf(fp, "  .avhrr_path = \"%s\",\n", params->avhrr_path);
  fprintf(fp, "  .ghrsst_type = 0x%x,\n", params->ghrsst_type);
  fprintf(fp, "  .ghrsst_path = \"%s\",\n", params->ghrsst_path);
  fprintf(fp, "  .ghrsst_opt = \"%s\",\n", params->ghrsst_opt);
  fprintf(fp, "  .ghrsst_name = \"%s\",\n", params->ghrsst_name);
  fprintf(fp, "  .ghrsst_dt = \"%s\",\n", params->ghrsst_dt);
  fprintf(fp, "  .trvars = \"%s\",\n", params->trvars);
  fprintf(fp, "  .smooth_v = \"%s\",\n", params->smooth_v);
  fprintf(fp, "  .scale_v = \"%s\",\n", params->scale_v);
  fprintf(fp, "  .orthoweights = \"%s\",\n", params->orthoweights);
  fprintf(fp, "  .nodal_dir = \"%s\",\n", params->nodal_dir);
  fprintf(fp, "  .tide_con_file = \"%s\",\n", params->tide_con_file);
  fprintf(fp, "  .webf = \"%s\",\n", params->webf);
  fprintf(fp, "  .regulate = \"%s\",\n", params->regulate);
  fprintf(fp, "  .regulate_dt = %f,\n", params->regulate_dt);
  fprintf(fp, "  .alert = \"%s\",\n", params->alert);
  fprintf(fp, "  .alertc = \"%s\",\n", params->alertc);
  fprintf(fp, "  .alert_dt = \"%s\",\n", params->alert_dt);
  fprintf(fp, "  .eta_f = %d,\n", params->eta_f);
  fprintf(fp, "  .vel2d_f = %d,\n", params->vel2d_f);
  fprintf(fp, "  .vel3d_f = %d,\n", params->vel3d_f);
  fprintf(fp, "  .wvel_f = %d,\n", params->wvel_f);
  fprintf(fp, "  .tend_f = %d,\n", params->tend_f);
  fprintf(fp, "  .div2d_f = %d,\n", params->div2d_f);
  fprintf(fp, "  .div3d_f = %d,\n", params->div3d_f);
  fprintf(fp, "  .cfl_f = %d,\n", params->cfl_f);
  fprintf(fp, "  .ts_f = %d,\n", params->ts_f);
  fprintf(fp, "  .shear_f = %d,\n", params->shear_f);
  fprintf(fp, "  .hdiff_f = %d,\n", params->hdiff_f);
  fprintf(fp, "  .tide_r = 0x%x,\n", params->tide_r);
  fprintf(fp, "  .parray_inter_botz = %d,\n", params->parray_inter_botz);
  fprintf(fp, "  .compatible = 0x%x,\n", params->compatible);
  fprintf(fp, "  .gint_errfcn = 0x%x,\n", params->gint_errfcn);
  fprintf(fp, "  .nland = %d,\n", params->nland);
  fprintf(fp, "  .nturb = %d,\n", params->nturb);
  fprintf(fp, "  .amax = %f,\n", params->amax);
  fprintf(fp, "  .hmax = %f,\n", params->hmax);
  fprintf(fp, "  .vmax = %f,\n", params->vmax);
  fprintf(fp, "  .btmax = %f,\n", params->btmax);
  fprintf(fp, "  .bcmax = %f,\n", params->bcmax);
  fprintf(fp, "  .cmax = %f,\n", params->cmax);
  fprintf(fp, "  .detamax = %f,\n", params->detamax);
  fprintf(fp, "  .dwmax = %f,\n", params->dwmax);
  fprintf(fp, "  .dtmax = %f,\n", params->dtmax);
  fprintf(fp, "  .dsmax = %f,\n", params->dsmax);
  fprintf(fp, "  .smax = %f,\n", params->smax);
  fprintf(fp, "  .nobc = %d,\n", params->nobc);
  /*fprintf(fp, "  .open = HYDRO_OBC_%s,\n", name);*/

#if defined(HAVE_SEDIMENT_MODULE)
  params->do_sed = 0.0;
  fprintf(fp, "  .do_sed = %d,\n", params->do_sed);
  fprintf(fp, "  .sed_vars = \"%s\",\n", params->sed_vars);
  fprintf(fp, "  .sed_defs = \"%s\",\n", params->sed_defs);
#endif

#if defined(HAVE_ECOLOGY_MODULE)
  fprintf(fp, "  .do_eco = %d,\n", params->do_eco);
  fprintf(fp, "  .eco_vars = \"%s\",\n", params->eco_vars);
  fprintf(fp, "  .eco_defs = \"%s\",\n", params->eco_defs);
#endif

#if defined(HAVE_WAVE_MODULE)
  fprintf(fp, "  .do_wave = 0x%x,\n", params->do_wave);
#endif
  fprintf(fp, "};\n\n");

}

/* END render_hydro_parameters()                                    */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Writes the hydro OBC parameters to the config header.            */
/*------------------------------------------------------------------*/
void render_OBC_parameters(master_t *master, FILE *fp, char *name)
{
  geometry_t *geom = master->geom;
  int nobc = geom->nobc;
  int n, i;

  fprintf(fp, "open_bdrys_t HYDRO_%s_obc[] = {\n", name);
  for (n = 0; n < nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    fprintf(fp, "  {\n");
    fprintf(fp, "    .type = 0x%x,\n", open->type);
    fprintf(fp, "    .id = %d,\n", open->id);
    fprintf(fp, "    .name = \"%s\",\n", open->name);
    fprintf(fp, "    .bcond_nor = 0x%x,\n", open->bcond_nor);
    fprintf(fp, "    .bcond_nor2d = 0x%x,\n", open->bcond_nor2d);
    fprintf(fp, "    .bcond_tan = 0x%x,\n", open->bcond_tan);
    fprintf(fp, "    .bcond_tan2d = 0x%x,\n", open->bcond_tan2d);
    fprintf(fp, "    .bcond_ele = 0x%x,\n", open->bcond_ele);
    fprintf(fp, "    .bcond_w = 0x%x,\n", open->bcond_w);
    fprintf(fp, "    .bcond_Vz = 0x%x,\n", open->bcond_Vz);
    fprintf(fp, "    .bcond_Kz = 0x%x,\n", open->bcond_Kz);
    fprintf(fp, "    .sbcond = %d,\n", open->sbcond);
    fprintf(fp, "    .stagger = 0x%x,\n", open->stagger);
    fprintf(fp, "    .options = 0x%x,\n", open->options);
    fprintf(fp, "    .adjust_flux = %e,\n", open->adjust_flux);
    fprintf(fp, "    .adjust_flux_s = %e,\n", open->adjust_flux_s);
    fprintf(fp, "    .inverse_barometer = %d,\n", open->inverse_barometer);
    fprintf(fp, "    .afr = %e,\n", open->afr);
    fprintf(fp, "    .spf = %e,\n", open->spf);
    fprintf(fp, "    .relax_zone_nor = %d,\n", open->relax_zone_nor);
    fprintf(fp, "    .relax_zone_tan = %d,\n", open->relax_zone_tan);
    fprintf(fp, "    .linear_zone_nor = %d,\n", open->linear_zone_nor);
    fprintf(fp, "    .linear_zone_tan = %d,\n", open->linear_zone_tan);
    fprintf(fp, "    .relax_zone_ele = %d,\n", open->relax_zone_ele);
    fprintf(fp, "    .relax_ele = %d,\n", open->relax_ele);
    fprintf(fp, "    .rnor_b = %e,\n", open->rnor_b);
    fprintf(fp, "    .rnor_i = %e,\n", open->rnor_i);
    fprintf(fp, "    .rtan_b = %e,\n", open->rtan_b);
    fprintf(fp, "    .rtan_i = %e,\n", open->rtan_i);
    fprintf(fp, "    .rele_b = %e,\n", open->rele_b);
    fprintf(fp, "    .rele_i = %e,\n", open->rele_i);
    fprintf(fp, "    .bathycon = %d,\n", open->bathycon);
    fprintf(fp, "    .smooth_z = %d,\n", open->smooth_z);
    fprintf(fp, "    .smooth_n = %d,\n", open->smooth_n);
    fprintf(fp, "    .bout = %d,\n", open->bout);
    fprintf(fp, "    .relax_time = %e,\n", open->relax_time);
    fprintf(fp, "    .relax_timei = %e,\n", open->relax_timei);
    fprintf(fp, "    .sponge_zone = %d,\n", open->sponge_zone);
    fprintf(fp, "    .sponge_zone_h = %e,\n", open->sponge_zone_h);
    fprintf(fp, "    .sponge_f = %e,\n", open->sponge_f);
    fprintf(fp, "    .meanc = %e,\n", open->meanc);
    fprintf(fp, "    .bflux_2d = %e,\n", open->bflux_2d);
    fprintf(fp, "    .bflux_3d = %e,\n", open->bflux_3d);
    fprintf(fp, "    .file_dt = %e,\n", open->file_dt);
    fprintf(fp, "    .file_next = %e,\n", open->file_next);
    fprintf(fp, "    .upmeth = 0x%x,\n", open->upmeth);
    fprintf(fp, "    .rlen = %e,\n", open->rlen);
    fprintf(fp, "    .mindep = %e,\n", open->mindep);
    fprintf(fp, "    .maxdep = %e,\n", open->maxdep);
    fprintf(fp, "    .meandep = %e,\n", open->meandep);
    fprintf(fp, "    .length = %e,\n", open->length);
    fprintf(fp, "    .area = %e,\n", open->area);
    fprintf(fp, "    .ncells = %e,\n", open->ncells);
    fprintf(fp, "    .bhc = %e,\n", open->bhc);
    fprintf(fp, "    .bgz = %d,\n", open->bgz);
    fprintf(fp, "    .v1 = %e,\n", open->v1);
    fprintf(fp, "    .v2 = %e,\n", open->v2);
    fprintf(fp, "    .v3 = %e,\n", open->v3);
    fprintf(fp, "    .bflow = \"%s\",\n", open->bflow);
    fprintf(fp, "    .nzone = \"%s\",\n", open->nzone);
    fprintf(fp, "    .intype = %d,\n", open->intype);
    fprintf(fp, "    .nedges = %d,\n", open->nedges);
    fprintf(fp, "    .minlat = %e,\n", open->minlat);
    fprintf(fp, "    .minlon = %e,\n", open->minlon);
    fprintf(fp, "    .maxlat = %e,\n", open->maxlat);
    fprintf(fp, "    .maxlon = %e,\n", open->maxlon);
    fprintf(fp, "    .elon = %e,\n", open->elon);
    fprintf(fp, "    .elat = %e,\n", open->elat);
    fprintf(fp, "    .slon = %e,\n", open->slon);
    fprintf(fp, "    .slat = %e,\n", open->slat);
    fprintf(fp, "    .mlon = %e,\n", open->mlon);
    fprintf(fp, "    .mlat = %e,\n", open->mlat);
    fprintf(fp, "    .tide_con = \"%s\",\n", open->tide_con);
    fprintf(fp, "    .bstdf = %d,\n", open->bstdf);
    fprintf(fp, "    .nbstd = %d,\n", open->nbstd);
    /*
    fprintf(fp, "    .ntr = %d,\n", open->ntr);
    fprintf(fp, "    .atr = %d,\n", open->atr);
    */
    /*
    for (i = 0; i < open->ntr; i++) {
      fprintf(fp, "    .bcond_tra[%d] = 0x%x,\n", i, open->bcond_tra[i]);
      fprintf(fp, "    .relax_zone_tra[%d] = %d,\n", i, open->relax_zone_tra[i]);
      fprintf(fp, "    .rtra_b[%d] = %e,\n", i, open->rtra_b[i]);
      fprintf(fp, "    .rtra_i[%d] = %e,\n", i, open->rtra_i[i]);
      fprintf(fp, "    .trpc[%d] = %e,\n", i, open->trpc[i]);
      fprintf(fp, "    .clampv[%d] = %e,\n", i, open->clampv[i]);
    }
    */
    if (open->bstd != NULL) {
      for (i = 0; i < NBSTD; i++) {
	fprintf(fp, "    .open->bstd[%d] = \"%s\",\n", i, open->bstd[i]);
      }
    }
    if (n == nobc - 1)
      fprintf(fp, "  }\n");
    else
    fprintf(fp, "  },\n");
  }
  fprintf(fp, "};\n\n");
}

/* END render_OBC_parameters()                                      */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Copies a hydro configuration to a datastructure                  */
/*------------------------------------------------------------------*/
void copy_hydroparameter_info(parameters_t *params, parameters_t *in) {
  int n, m;

  strcpy(params->codeheader, in->codeheader);
  strcpy(params->parameterheader, in->parameterheader);
  strcpy(params->reference, in->reference);
  strcpy(params->reference, in->reference);
  strcpy(params->notes, in->notes);
  strcpy(params->trl, in->trl);
  strcpy(params->opath, in->opath);
  strcpy(params->gridtype, in->gridtype);
  params->gridcode, in->gridcode;
  strcpy(params->projection, in->projection);
  params->us_type = in->us_type;
  strcpy(params->trkey, in->trkey);
  strcpy(params->sequence, in->sequence);
  strcpy(params->bdrypath, in->bdrypath);
  strcpy(params->tracerdata, in->tracerdata);
  strcpy(params->rivldir, in->rivldir);
  strcpy(params->trans_dt, in->trans_dt);
  strcpy(params->runnoc, in->runnoc);
  params->runno = in->runno;
  strcpy(params->runcode, in->runcode);
  strcpy(params->rev, in->rev);
  params->history = in->history;
  /*
  strcpy(params->start_time, in->start_time);
  strcpy(params->stop_time, in->stop_time);
  strcpy(params->timeunit, in->timeunit);
  strcpy(params->output_tunit, in->output_tunit);
  params->t = in->t;
  */
  strcpy(params->lenunit, in->lenunit);
  params->momsc = in->momsc;
  params->trasc = in->trasc;
  params->ultimate = in->ultimate;
  params->osl = in->osl;
  params->rkstage = in->rkstage;
  strcpy(params->smag, in->smag);
  params->smagorinsky = in->smagorinsky;
  params->u1vh = in->u1vh;
  params->u1kh = in->u1kh;
  params->sue1 = in->sue1;
  params->kue1 = in->kue1;
  params->bsue1 = in->bsue1;
  params->bkue1 = in->bkue1;
  params->smag_smooth = in->smag_smooth;
  params->diff_scale = in->diff_scale;
  params->visc_method = in->visc_method;
  params->visc_fact = in->visc_fact;
  strcpy(params->mixsc, in->mixsc);
  strcpy(params->s_func, in->s_func);
  params->min_tke = in->min_tke;
  params->min_diss = in->min_diss;
  params->vz0 = in->vz0;
  params->kz0 = in->kz0;
  params->zs = in->zs;
  params->stab = in->stab;
  /*
  params->atr = in->atr;
  params->ntr = in->ntr;
  params->ntrS = in->ntrS;
  params->nsed = in->nsed;
  */
  params->sednz = in->sednz;
  params->ambpress = in->ambpress;
  params->air_dens = in->air_dens;
  params->spec_heat = in->spec_heat;
  params->g = in->g;
  params->grid_dt = in->grid_dt;
  params->iratio = in->iratio;
  params->tratio = in->tratio;
  params->trsplit = in->trsplit;
  strcpy(params->autotrpath, in->autotrpath);
  params->cfl = in->cfl;
  strcpy(params->cfl_dt, in->cfl_dt);
  params->mixlayer = in->mixlayer;
  params->show_layers = in->show_layers;
  params->means = in->means;
  strcpy(params->means_os, in->means_os);
  strcpy(params->bathystats, in->bathystats);
  strcpy(params->particles, in->particles);
  strcpy(params->addquad, in->addquad);
  params->vorticity = in->vorticity;
  params->numbers = in->numbers;
  params->numbers1 = in->numbers1;
  params->u1_f = in->u1_f;
  params->u1av_f = in->u1av_f;
  params->save_force = in->save_force;
  params->lnm = in->lnm;
  strcpy(params->trflux, in->trflux);
  params->trfd1 = in->trfd1;
  params->trfd2 = in->trfd2;
  strcpy(params->trperc, in->trperc);
  strcpy(params->trpercr, in->trpercr);
  params->trflsh = in->trflsh;
  strcpy(params->trage, in->trage);
  strcpy(params->trtend, in->trtend);
  params->ndhw = in->ndhw;
  params->tendf = in->tendf;
  params->thin_merge = in->thin_merge;
  params->maxgrad = in->maxgrad;
  strcpy(params->maxdiff, in->maxdiff);
  params->bathyfill = in->bathyfill;
  params->bvel = in->bvel;
  params->nonlinear = in->nonlinear;
  params->calc_dens = in->calc_dens;
  params->mode2d = in->mode2d;
  params->bmin = in->bmin;
  params->bmax = in->bmax;
  strcpy(params->mct, in->mct);
  params->maxgrad = in->maxgrad;
  strcpy(params->maxdiff, in->maxdiff);
  params->smooth = in->smooth;
  params->slipprm = in->slipprm;
  params->nwindows = in->nwindows;
  params->win_reset = in->win_reset;
  params->win_type = in->win_type;
  params->win_block = in->win_block;
  strcpy(params->win_file, in->win_file);
  for (n = 0; n < in->nwindows; n++) {
    if (in->win_size) params->win_size[n] = in->win_size[n];
    if (in->nwn) {
      params->nwn[n] = in->nwn[n];
      for (m = 0; m < params->nwn[n]; m++) {
	if (in->wnx) params->wnx[n][m] = in->wnx[n][m];
	if (in->wny) params->wny[n][m] = in->wny[n][m];
      }
    }
  }
  params->show_win = in->show_win;
  params->restart_dt = in->restart_dt;
  params->da = in->da;
  params->da_dt = in->da_dt;
  params->da_fcst_dt = in->da_fcst_dt;
  params->prex = in->prex;
  params->ndf = in->ndf;
  params->riverflow = in->riverflow;
  params->meshinfo = in->meshinfo;
  params->swr_type = in->swr_type;
  for (n = 0; n < 6; n++)
    params->swr_ens[n] = in->swr_ens[n];
  strcpy(params->u1vhc, in->u1vhc);
  strcpy(params->u1khc, in->u1khc);
  strcpy(params->restart_name, in->restart_name);
  strcpy(params->wind_file, in->wind_file);
  strcpy(params->wind, in->wind);
  params->dlv0 = in->dlv0;
  params->dlv1 = in->dlv1;
  params->dlc0 = in->dlc0;
  params->dlc1 = in->dlc1;
  params->wind_type = in->wind_type;
  params->stress_fn = in->stress_fn;
  params->wind_scale = in->wind_scale;
  strcpy(params->bdry_file, in->bdry_file);
  strcpy(params->ptinname, in->ptinname);
  strcpy(params->patm, in->patm);
  strcpy(params->precip, in->precip);
  strcpy(params->evap, in->evap);
  strcpy(params->airtemp, in->airtemp);
  strcpy(params->rh, in->rh);
  strcpy(params->cloud, in->cloud);
  strcpy(params->swr, in->swr);
  strcpy(params->light, in->light);
  strcpy(params->wetb, in->wetb);
  strcpy(params->hftemp, in->hftemp);
  strcpy(params->hf, in->hf);
  strcpy(params->swr_babs, in->swr_babs);
  strcpy(params->swr_attn, in->swr_attn);
  strcpy(params->swr_attn1, in->swr_attn1);
  strcpy(params->swr_tran, in->swr_tran);
  strcpy(params->swr_regions, in->swr_regions);
  params->swreg_dt = in->swreg_dt;
  strcpy(params->swr_data, in->swr_data);
  strcpy(params->densname, in->densname);
  strcpy(params->regions, in->regions);
  strcpy(params->region_dt, in->region_dt);
  strcpy(params->region_vars, in->region_vars);
  strcpy(params->region_mode, in->region_mode);
  strcpy(params->i_rule, in->i_rule);
  strcpy(params->cookiecut, in->cookiecut);
  strcpy(params->imp2df, in->imp2df);
  strcpy(params->imp3df, in->imp3df);
  params->do_pt = in->do_pt;
  params->do_lag = in->do_lag;
  strcpy(params->dp_mode, in->dp_mode);
  params->waves = in->waves;
  params->decorr = in->decorr;
  params->decf = in->decf;
  strcpy(params->monotr, in->monotr);
  params->monomn = in->monomn;
  params->monomx = in->monomx;
  params->orbital = in->orbital;
  params->rampf = in->rampf;
  params->rampstart = in->rampstart;
  params->rampend = in->rampend;
  params->hmin = in->hmin;
  params->uf = in->uf;
  strcpy(params->quad_bfc, in->quad_bfc);
  params->totals = in->totals;
  params->roammode = in->roammode;
  params->robust = in->robust;
  params->speed = in->speed;
  params->fatal = in->fatal;
  params->avhrr = in->avhrr;
  params->ghrsst = in->ghrsst;
  params->exmapf = in->exmapf;
  params->heatflux = in->heatflux;
  params->saltflux = in->saltflux;
  params->water_type = in->water_type;
  params->tsfile_caching = in->tsfile_caching;
  params->etamax = in->etamax;
  params->velmax = in->velmax;
  params->velmax2d = in->velmax2d;
  params->etadiff = in->etadiff;
  params->wmax = in->wmax;
  params->trsplit = in->trsplit;
  params->eta_ib = in->eta_ib;
  params->etarlx = in->etarlx;
  params->velrlx = in->velrlx;
  params->albedo = in->albedo;
  params->domom = in->domom;
  params->doroms = in->doroms;
  params->doswan = in->doswan;
  params->zref = in->zref;
  params->neutral = in->neutral;
  params->wave_alpha = in->wave_alpha;
  params->wave_hf = in->wave_hf;
  params->wave_b1 = in->wave_b1;
  params->fillf = in->fillf;
  params->filter = in->filter;
  params->trfilter = in->trfilter;
  params->conserve = in->conserve;
  params->lyear = in->lyear;
  params->do_closure = in->do_closure;
  params->rsalt = in->rsalt;
  params->rtemp = in->rtemp;
  params->smooth_VzKz = in->smooth_VzKz;
  params->albedo_l = in->albedo_l;
  params->trout = in->trout;
  params->porusplate = in->porusplate;
  params->sharp_pyc = in->sharp_pyc;
  params->uscf = in->uscf;
  params->tidef = in->tidef;
  params->tidep = in->tidep;
  params->eqt_alpha = in->eqt_alpha;
  params->eqt_beta = in->eqt_beta;
  params->mlat = in->mlat;
  params->mlon = in->mlon;
  strcpy(params->nprof, in->nprof);
  strcpy(params->nprof2d, in->nprof2d);
  strcpy(params->reef_frac, in->reef_frac);
  strcpy(params->crf, in->crf);
  params->dbi = in->dbi;
  params->dbj = in->dbj;
  params->dbk = in->dbk;
  params->dbgf = in->dbgf;
  params->dbgtime = in->dbgtime;
  strcpy(params->momfile, in->momfile);
  strcpy(params->avhrr_path, in->avhrr_path);
  params->ghrsst_type = in->ghrsst_type;
  strcpy(params->ghrsst_path, in->ghrsst_path);
  strcpy(params->ghrsst_opt, in->ghrsst_opt);
  strcpy(params->ghrsst_name, in->ghrsst_name);
  strcpy(params->ghrsst_dt, in->ghrsst_dt);
  strcpy(params->trvars, in->trvars);
  strcpy(params->smooth_v, in->smooth_v);
  strcpy(params->scale_v, in->scale_v);
  strcpy(params->orthoweights, in->orthoweights);
  strcpy(params->nodal_dir, in->nodal_dir);
  strcpy(params->tide_con_file, in->tide_con_file);
  params->rampf = in->rampf;
  strcpy(params->webf, in->webf);
  strcpy(params->regulate, in->regulate);
  params->regulate_dt = in->regulate_dt;
  strcpy(params->alert, in->alert);
  strcpy(params->alertc, in->alertc);
  strcpy(params->alert_dt, in->alert_dt);
  params->eta_f = in->eta_f;
  params->vel2d_f = in->vel2d_f;
  params->vel3d_f = in->vel3d_f;
  params->wvel_f = in->wvel_f;
  params->tend_f = in->tend_f;
  params->div2d_f = in->div2d_f;
  params->div3d_f = in->div3d_f;
  params->cfl_f = in->cfl_f;
  params->ts_f = in->ts_f;
  params->shear_f = in->shear_f;
  params->hdiff_f = in->hdiff_f;
  params->tide_r = in->tide_r;
  params->parray_inter_botz = in->parray_inter_botz;
  params->compatible = in->compatible;
  params->gint_errfcn = in->gint_errfcn;
  params->nland = in->nland;
  params->nturb = in->nturb;
  params->amax = in->amax;
  params->hmax = in->hmax;
  params->vmax = in->vmax;
  params->btmax = in->btmax;
  params->bcmax = in->bcmax;
  params->cmax = in->cmax;
  params->detamax = in->detamax;
  params->dwmax = in->dwmax;
  params->dtmax = in->dtmax;
  params->dsmax = in->dsmax;
  params->smax = in->smax;
  strcpy(params->rendername, in->rendername);
  strcpy(params->renderpath, in->renderpath);
  strcpy(params->renderdesc, in->renderdesc);
  strcpy(params->rendertype, in->rendertype);
  strcpy(params->ecosedconfig, in->ecosedconfig);
  params->nobc = in->nobc;

#if defined(HAVE_SEDIMENT_MODULE)
  params->do_sed = in->do_sed;
  strcpy(params->sed_vars, in->sed_vars);
  strcpy(params->sed_defs, in->sed_defs);
#endif

#if defined(HAVE_ECOLOGY_MODULE)
  params->do_eco = in->do_eco;
  strcpy(params->eco_vars, in->eco_vars);
  strcpy(params->eco_defs, in->eco_defs);
#endif

#if defined(HAVE_WAVE_MODULE)
  params->do_wave = in->do_wave;
#endif

}

/* END copy_hydroparameter_info()                                   */
/*------------------------------------------------------------------*/

