/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/tracerstats/tracerstats.c
 *  
 *  Description:
 *  Tracer statistics initialisation and interface routines
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: tracerstats.c 5843 2018-06-29 02:17:55Z riz008 $
 *
 */

#include "tracerstats.h"
#include "string_utils.h"
#include <ctype.h>

/*
 * These flags are set for different stat functions
 */
#define O_1 1   /* Copy args[0] to operation_2d / operation_3d */
#define O_2 2   /* Copy args[1] to tr_2d[on][1] / tr_3d[on][1] */
#define O_3 4   /* Copy args[2] to tr_2d[on][2] / tr_3d[on][2] */
#define O_4 8   /* Copy args[3] to tr_2d[on][3] / tr_3d[on][3] */
#define B_2 16
#define B_3 32
#define OS_1 128
#define OS_2 256
#define OS_B_2 512

#define EX_2D 1  /* Extra 2d tracers: eta                            */
#define EX_3D 6  /* Extra 3d tracers, Kz, Vz, u1kh, u2kh, u1vh, u2vh */

#define MAXSTATFN 7

#define STATFN0 0
#define STATFN1 1
#define STATFN2 2
#define STATFN3 3
#define STATFN4 4
#define STATFN5 5
#define STATFN6 6
#define MAXSECTIONPARAMS 5
#define OUTPUTRESOLUTION 4


static int decode_args(char *p, char *args[], int maxn);
static void corr_3d_pre(trs_t *trs, int* n, void* data, int c);
static void corr_3d_post(trs_t *trs, int* n, void* data, int c);
static void corr_2d_pre(trs_t *trs, int* n, void* data);
static void corr_2d_post(trs_t *trs, int* n, void* data);
static int read_section_coordinates(char* fnamein,  /* input fiilename */
                                      sectioncoord_t* sectiondata,
                                      void* model);
static int* split_ints(char* arg,int* n );
static int count_elements(char* string, int del);
static void init_sectiondata(sectioncoord_t* sectiondata, char* names, char* outscales, char** t3dnames, int m);
static void init_sum(char* names, int* nsumtrs, char** t3dnames, int m);
static int extract_step(void* model, char* st);
void section_scatter(int n, void* origin, void* target);
void section_gather(int n, void* target, void* origin);
void* section_create(void* d);
void fill_stat3dfcn(trs_t *trs, char  *op_3d, char  **tr_3d,
		    lstatfn3d_t *statfn3d, char **trnames, 
		    int m, stat3d_work_arrays* w_arr, int *st_type);

/* Version information */
int get_tracerstats_major_vers(void)
{
  return(TRACERSTATS_MAJOR_VERSION);
}

int get_tracerstats_minor_vers(void)
{
  return(TRACERSTATS_MINOR_VERSION);
}

int get_tracerstats_patch_vers(void)
{
  return(TRACERSTATS_PATCH_VERSION);
}

/*-------------------------------------------------------------------*/
/* Builds and initialises the interface                              */
/*-------------------------------------------------------------------*/
trs_t* trs_build(void* model, FILE *fp)
{
  trs_t* trs = trs_create();

  trs_init(model, trs, fp);
  return trs;
}
/* END trs_build()                                                   */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Allocates memory for the interface structure                      */
/*-------------------------------------------------------------------*/
trs_t* trs_create() {
  trs_t *trs = (trs_t *)malloc(sizeof(trs_t));

  trs->do_tracerstats = 0;
  trs->cols = -1;
  trs->nz = -1;
  trs->sednz = -1;
  trs->ntr = 0;
  trs->ntrS = 0;
  trs->atr = 0;
  trs->atrS = 0;
  trs->nsed = 0;
  trs->e_nstep = NULL;
  trs->o_nstep = NULL;
  trs->new_step = -1;
  trs->news = NULL;
  trs->dt = 0.0;
  trs->use_eta = 0;
  trs->use_etamean = 0;
  trs->use_Kz = 0;
  trs->use_Vz = 0;
  trs->use_fluxe1 = 0;
  trs->use_fluxe2 = 0;
  trs->use_w = 0;
  trs->use_u1 = 0;
  trs->use_u2 = 0;
  trs->use_area_e1_e2 = 0;
  trs->nstep = NULL;
  trs->nprestep = NULL;
  trs->stat_type = NULL;
  trs->stat_type_sed = NULL;
  trs->stat_typeS = NULL;
  trs->statfn3d = NULL;
  trs->statfn2d = NULL;
  trs->fluxfn = NULL;
  trs->domainfn3d = NULL;
  trs->nstatfn3d = 0;
  trs->nstatfn2d = 0;
  trs->nfluxfn = 0;
  trs->ndomainfn3d = 0;
  trs->w_wc.w1 = NULL;
  trs->w_wc.w2 = NULL;
  trs->w_wc.w3 = NULL;
  trs->w_wc.w4 = NULL;
  trs->w_wc.w5 = NULL;
  trs->w_wc.w6 = NULL;
  trs->w_sed.w1 = NULL;
  trs->w_sed.w2 = NULL;
  trs->w_sed.w3 = NULL;
  trs->w_sed.w4 = NULL;
  trs->w1S = NULL;
  trs->w2S = NULL;
  trs->w3S = NULL;
  trs->map_2d = NULL;
  trs->map_3d = NULL;
  trs->map_sed = NULL;
  trs->topk_wc = -1;
  trs->botk_wc = -1;
  trs->topk_sed = -1;
  trs->botk_sed = -1;
  trs->dz_wc = NULL;
  trs->Kz = NULL;
  trs->Vz = NULL;
  trs->u1 = NULL;
  trs->u2 = NULL;
  trs->w = NULL;
  trs->u1flux3d = NULL;
  trs->u2flux3d = NULL;
  trs->area_w = 0.0;
  trs->area_e1 = NULL;
  trs->area_e2 = NULL;
  trs->h1au2 = 0.0;
  trs->h2au1 = 0.0;
  trs->tr_wc = NULL;
  trs->tr_in = NULL;
  trs->tr_sed = NULL;
  trs->trname_3d = NULL;
  trs->trname_2d = NULL;
  trs->trname_sed = NULL;
  trs->tmap_3d = NULL;
  trs->tmap_2d = NULL;
  trs->tmap_sed = NULL;
  trs->depth_wc = -1;
  trs->eta = -1;
  trs->model = NULL;
  return trs;
}
/* END trs_create()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Initialises the tracer statistics structure using information     */
/* from the host model.                                              */
/*-------------------------------------------------------------------*/
void trs_init(void* model, trs_t *trs, FILE *fp) {
  int n, m, i, ntr_tmp, ntrsed_tmp;
  int str;
  int ntr, ntrS,ntrsed;
  int nm, nmS, nm_sed;
  int on, onS, on_sed;
  int stat;
  int tint;
  char **trname_3d;
  char **trname_2d;
  char **trname_sed;
  char buf[MAXSTRLEN], key[MAXSTRLEN], tbuf[MAXSTRLEN];
  char* pbuf[MAXSTRLEN];
  char tok[MAXSTRLEN];
  char *args[256];
  char **operation_3d;
  char **operation_sed;
  char **operation_2d;
  char ***tr_3d  = NULL;
  char ***tr_sed = NULL;
  char ***tr_2d = NULL;
  int *om_2d = NULL;
  int *om_3d = NULL;
  int *om_sed = NULL;
  int nargs;
  int *nstatorder3d ;
  int *nstatorder_sed ;
  int *nstatorder2d ;
  int nflux = 0 ;
  int ndomain = 0 ;
  int fluxdir = -1;
  sectioncoord_t* sectiondata = NULL;
  vdiff_range_t* vdiff = NULL;

  /* allocate memory for the maximum of
   * stat functions
   */
  nstatorder3d   = i_alloc_1d(MAXSTATFN);
  nstatorder_sed = i_alloc_1d(MAXSTATFN);
  nstatorder2d   = i_alloc_1d(MAXSTATFN);
  trs->nprestep  = i_alloc_1d(MAXSTEPS);
  memset(trs->nprestep, 0, MAXSTEPS * sizeof(int));

  /* Allocate memory for the step flags */
  trs->news = i_alloc_1d(MAXSTATFN);
  trs->o_nstep = i_alloc_1d(MAXSTATFN);
  trs->e_nstep = i_alloc_1d(MAXSTATFN);

  for (n = 0; n < MAXSTATFN ; n++)
  {
    nstatorder3d[n]   = 0;
    nstatorder_sed[n] = 0;
    nstatorder2d[n] = 0;
    trs->news[n] = -1;
  }

  /*-----------------------------------------------------------------*/
  /* Set up the grid                                                 */
  /* Get the number of water column layers                           */
  trs->nz = i_get_num_wclayers(model);

  /* Allocate memory for water column layer thickness                */
  trs->dz_wc = d_alloc_1d(trs->nz);


  /* Get the number of sediment layers                               */
  trs->sednz = i_get_num_sedlayers(model);

  /* Get the number of water columns                                 */
  trs->cols = i_get_num_columns(model);

  /* Get the timestep                                                */
  trs->dt = i_get_model_timestep(model);

  /* Get the model start time                                        */
  trs->time = i_get_model_time(model);

  /* Get the 3d tracers names                                        */
  ntr = i_get_num_tracers_3d(model, &trs->atr) + EX_3D;

  if(prm_read_char(fp, "NRTSTATS", buf))
  {
    i =atoi(buf);
  }else
    i = 0;

  /* assign a reference to the parent model */
  trs->model = model;

  trname_3d = malloc((ntr + i) * sizeof(char *));
  operation_3d = malloc((ntr + i) * sizeof(char *));
  for (n = 0; n < ntr; n++) {
    trname_3d[n] = malloc(MAXSTRLEN * sizeof(char));
    operation_3d[n] = malloc(MAXSTRLEN * sizeof(char));
  }
  i_get_names_tracers_3d(model, trname_3d);
  if (ntr || i) {
    tr_3d = c_alloc_3d(MAXSTRLEN, 9, (ntr + i));
    om_3d = i_alloc_1d((ntr + i));
  }

  /* Get the 2d tracers names                                        */
  ntrS = i_get_num_tracers_2d(model, &trs->atrS) + EX_2D;
  trname_2d = malloc(ntrS * sizeof(char *));
  operation_2d = malloc(ntrS * sizeof(char *));
  for (n = 0; n < ntrS; n++) {
    trname_2d[n] = malloc(MAXSTRLEN * sizeof(char));
    operation_2d[n] = malloc(MAXSTRLEN * sizeof(char));
  }
  i_get_names_tracers_2d(model, trname_2d);
  if (ntrS) {
    tr_2d = c_alloc_3d(MAXSTRLEN, 5, ntrS);
    om_2d = i_alloc_1d(ntrS);
  }

  /* Get the sediment tracers names                                        */
  ntrsed = i_get_num_sedtracers(model);
  if(ntrsed)
  {
    trname_sed    = malloc(ntrsed * sizeof(char *));
    operation_sed = malloc(ntrsed * sizeof(char *));
    for (n = 0; n < ntrsed; n++) {
      trname_sed[n]    = malloc(MAXSTRLEN * sizeof(char));
      operation_sed[n] = malloc(MAXSTRLEN * sizeof(char));
    }
    i_get_names_tracers_sed(model, trname_sed);
    tr_sed = c_alloc_3d(MAXSTRLEN, 4, ntrsed);
    om_sed = i_alloc_1d(ntrsed);
  }

  /*-----------------------------------------------------------------*/
  /* Count the number of tracers requiring operations                */
  nm = nmS = nm_sed = 1;
  on = onS = on_sed = 0;

  /* 3D tracers                                                      */
  for (n = 0; n < ntr - EX_3D; n++) {
    stat = 0;
    fluxdir = -1;
    sprintf(key, "TRACER%1.1d.tracerstat", i_get_param_map_3d(model, n));
    if (prm_read_char(fp, key, tbuf)) {
      nargs = decode_args(tbuf, args, 256);
      if (nargs < 2) {
        emstag(LFATAL,"tracerstats:init"," statistics : format <operation>(<tracer>)\n");
        exit(0);
      }
      strcpy(tok, args[0]);
      ntr_tmp = trs->ntr;
      /* Count the number of 3d tracers requiring means              */
      if (strcmp(tok, "mean") == 0) {
        trs->do_tracerstats = 1;
        stat = O_1|O_2;
        /* Read in the averaging period and store temporarily in     */
        /* tr_3d[2].                                                 */
        sprintf(key, "TRACER%1.1d.dt", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][2])) {
          sprintf(key, "STOP_TIME");
          prm_read_char(fp, key, tr_3d[on][2]);
        }
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][3])) {
          sprintf(tr_3d[on][3],"%u", LASTSTEP);
        }

        emstag(LTRACE,"tracerstats:trs_init",
	       "statistics : %s = %s mean of 3d tracer %s\n",
	       trname_3d[n], tr_3d[on][2], args[1]);
        nstatorder3d[STATFN1]++;

	if(ntrsed) {
	  emstag(LTRACE,"tracerstats:trs_init",
		 "statistics : %s = %s mean of sediment tracer %s\n",
		 trname_3d[n], tr_3d[on][2], args[1]);
	  nstatorder_sed[STATFN1]++;
	}
      }else
      if (strcmp(tok, "run_mean") == 0) {
        trs->do_tracerstats = 1;
        stat = O_1|O_2;
        /* Read in the averaging period and store temporarily in     */
        /* tr_3d[2].                                                 */
        sprintf(key, "TRACER%1.1d.dt", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][2])) {
          sprintf(key, "STOP_TIME");
          prm_read_char(fp, key, tr_3d[on][2]);
        }
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][3])) {
          sprintf(tr_3d[on][3],"%u", LASTSTEP);
        }

        emstag(LTRACE,"tracerstats:trs_init",
	       "statistics : %s = %s run_mean of 3d tracer %s\n",
	       trname_3d[n], tr_3d[on][2], args[1]);
        nstatorder3d[STATFN1]++;

	if(ntrsed) {
	  emstag(LTRACE,"tracerstats:trs_init",
		 "statistics : %s = %s run_mean of sediment tracer %s\n",
		 trname_3d[n], tr_3d[on][2], args[1]);
	  nstatorder_sed[STATFN1]++;
	}
      }else

      /* Count the number of 3d tracers requiring integration        */
      if (strcmp(tok, "sum@dt") == 0) {
        trs->do_tracerstats = 1;
        stat = O_1|O_2;
        /* Read in the averaging period and store temporarily in     */
        /* tr_3d[2].                                                 */
        sprintf(key, "TRACER%1.1d.dt", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][2])) {
          sprintf(key, "STOP_TIME");
          prm_read_char(fp, key, tr_3d[on][2]);
        }
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][3])) {
          sprintf(tr_3d[on][3],"%u", LASTSTEP);
        }

        emstag(LTRACE,"tracerstats:trs_init",
	       "statistics : %s = sum every %s of 3d tracer %s\n",
	       trname_3d[n], tr_3d[on][2], args[1]);
        nstatorder3d[STATFN1]++;

	if(ntrsed) {
	  emstag(LTRACE,"tracerstats:trs_init",
		 "statistics : %s = sum every %s of sediment tracer %s\n",
		 trname_3d[n], tr_3d[on][2], args[1]);
	  nstatorder_sed[STATFN1]++;
	}
      }else

      /* Count the number of 3d tracers requiring max */
      if (strcmp(tok, "max") == 0) {
        trs->do_tracerstats = 1;
        stat = O_1|O_2;
	/*
	 * Read in the period over which to max and store temporarily in
	 * tr_3d[2]. If dt is not present, we do it over the whole period.
	 */
        sprintf(key, "TRACER%1.1d.dt", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][2])) {
          sprintf(key, "STOP_TIME");
          prm_read_char(fp, key, tr_3d[on][2]);
        }
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][3])) {
          sprintf(tr_3d[on][3],"%u", LASTSTEP);
        }

        emstag(LTRACE, "tracerstats:trs_init",
	       "statistics : %s = %s max of 3d tracer %s\n",
	       trname_3d[n], args[1]);
        nstatorder3d[STATFN1]++;

	if(ntrsed) {
	  emstag(LTRACE, "tracerstats:trs_init",
		 "statistics : %s = %s max of sediment tracer %s\n",
		 trname_3d[n], args[1]);
	  nstatorder_sed[STATFN1]++;
	}
      }else

      /* Count the number of 3d tracers requiring max */
      if (strcmp(tok, "min") == 0) {
        trs->do_tracerstats = 1;
        stat = O_1|O_2;
	/*
	 * Read in the period over which to min and store temporarily in
	 * tr_3d[2]. If dt is not present, we do it over the whole period.
	 */
        sprintf(key, "TRACER%1.1d.dt", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][2])) {
          sprintf(key, "STOP_TIME");
          prm_read_char(fp, key, tr_3d[on][2]);
        }
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][3])) {
          sprintf(tr_3d[on][3],"%u", LASTSTEP);
        }

        emstag(LTRACE, "tracerstats:trs_init",
	       "statistics : %s = %s min of 3d tracer %s\n",
	       trname_3d[n], args[1]);
        nstatorder3d[STATFN1]++;

	if(ntrsed) {
	  emstag(LTRACE, "tracerstats:trs_init",
		 "statistics : %s = %s min of sediment tracer %s\n",
		 trname_3d[n], args[1]);
	  nstatorder_sed[STATFN1]++;
	}
      }else

      /* Count the number of 3d tracers requiring Degree Heating     */
      /* Weeks (a measure of coral bleaching).                       */
      if (strcmp(tok, "dhw") == 0) {
        trs->do_tracerstats = 1;
        stat = O_1|O_2|O_3;
        sprintf(key, "TRACER%1.1d.dt", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][2]))
	  strcpy(tr_3d[on][2], "7 days");
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][3])) {
          sprintf(tr_3d[on][3],"%u", LASTSTEP);
        }
        emstag(LTRACE,"tracerstats:trs_init","statistics : %s = DHW of 3d tracer %s using %s\n",
           trname_3d[n], args[1], args[2]);
        nstatorder3d[STATFN0]++;

      }else

      /* Count the number of 3d tracers requiring temperature        */
      /* exposure (can be configured to  measure coral bleaching).   */
      if (strcmp(tok, "exposure") == 0) {
        trs->do_tracerstats = 1;
        stat = O_1|O_2|O_3|O_4;
        sprintf(key, "TRACER%1.1d.dt", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][5]))
	  strcpy(tr_3d[on][5], "7 days");

        sprintf(key, "TRACER%1.1d.start", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][6]))
	  sprintf(tr_3d[on][6], "%f seconds", trs->time);
	// Also read the time units
	sprintf(key, "TRACER%1.1d.tunit", i_get_param_map_3d(model, n));
	if (!prm_read_char(fp, key, tr_3d[on][7])) {
	  sprintf(tr_3d[on][7], "since 1990-01-01 00:00:00 +10");
	  emstag(LWARN,"tracerstats:init","Missing parameter tunit for exposure, assuming %s\n", tr_3d[on][7]);
	}

        sprintf(key, "TRACER%1.1d.step", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][3])) {
          sprintf(tr_3d[on][3],"%u", LASTSTEP);

        sprintf(key, "TRACER%1.1d.scale_factor", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][8]))
	  strcpy(tr_3d[on][8], "1 days");
        }
        emstag(LTRACE,"tracerstats:trs_init","statistics : %s = exposure of 3d tracer %s using threshold %s\n",
           trname_3d[n], args[1], args[2]);
        nstatorder3d[STATFN0]++;

      }else

      /* Count the number of 3d tracers requiring ReefTemp exposure  */
      /* (a measure of coral bleaching).                             */
      if (strcmp(tok, "reeftemp") == 0) {
        trs->do_tracerstats = 1;
        stat = O_1|O_2|O_3;

        sprintf(key, "TRACER%1.1d.step", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][3])) {
          sprintf(tr_3d[on][3],"%u", LASTSTEP);

        sprintf(key, "TRACER%1.1d.scale_factor", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][4]))
	  strcpy(tr_3d[on][8], "1 days");
        }
        emstag(LTRACE,"tracerstats:trs_init","statistics : %s = ReefTemp exposure of 3d tracer %s using threshold %s\n",
           trname_3d[n], args[1], args[2]);
        nstatorder3d[STATFN0]++;

      }else

      if (strcmp(tok, "copy") == 0) {
        trs->do_tracerstats = 1;
        stat = O_1|O_2;
        emstag(LTRACE,"tracerstats:trs_init","statistics : %s = copy tracer %s\n",
           trname_3d[n], args[1]);
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][3])) {
          sprintf(tr_3d[on][3],"%u", LASTSTEP);
        }
        nstatorder3d[STATFN1]++;
      } else

      /* Count number of 3d tracers requiring standard deviations    */
      if (strcmp(tok, "stdev") == 0) {
        trs->do_tracerstats = 1;
        stat = O_1|B_2|OS_1;
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][3])) {
          sprintf(tr_3d[on][3],"%u", LASTSTEP);
        }
        emstag(LTRACE,"tracerstats:init",
	       "statistics : %s = standard deviation of 3d tracer %s\n",
	       trname_3d[n], args[1]);
        nstatorder3d[STATFN4]++;

	if(ntrsed) {
	  emstag(LTRACE,"tracerstats:init",
		 "statistics : %s = standard deviation of sediment tracer %s\n",
		 trname_3d[n], args[1]);
	  nstatorder_sed[STATFN4]++;
	}
      } else
      /* Count the number of 3d tracers requiring variances          */
      if (strcmp(tok, "variance") == 0) {
	trs->do_tracerstats = 1;
	stat = O_1|B_2|OS_1;
	sprintf(key, "TRACER%1.1d.step", i_get_param_map_3d(model, n));
	if (!prm_read_char(fp, key, tr_3d[on][3])) {
          sprintf(tr_3d[on][3],"%u", LASTSTEP);
        }
        emstag(LTRACE,"tracerstats:init",
	       "statistics : %s = variance of 3d tracer %s\n",
	       trname_3d[n], args[1]);
        nstatorder3d[STATFN3]++;
	
	if(ntrsed) {
	  emstag(LTRACE,"tracerstats:init",
		 "statistics : %s = variance of sediment tracer %s\n",
		 trname_sed[n], args[1]);
	  nstatorder_sed[STATFN3]++;
	}
      } else
       
      /* Count the number of 3d tracers requiring covariances        */
      if (strcmp(tok, "cov") == 0) {
        trs->do_tracerstats = 1;
        nm = max(nm, 2);
        stat = O_1|B_2|B_3;
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][3])) {
          sprintf(tr_3d[on][3],"%u", LASTSTEP);
        }
        emstag(LTRACE,"tracerstats:init","statistics : %s = covariance of 3d tracer %s and %s\n",
	       trname_3d[n], args[1], args[2]);
        nstatorder3d[STATFN6]++;
      } else 
      /* Count the number of 3d tracers requiring correlation        */
      /* coefficients.                                               */
      if (strcmp(tok, "corr") == 0) {
        trs->do_tracerstats = 1;
        nm = max(nm, 4);
        stat = O_1|B_2|B_3;
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][3])) {
          sprintf(tr_3d[on][3],"%u", LASTSTEP);
        }
        emstag(LTRACE,"tracerstats:init","statistics : %s = correlation coefficient of 3d tracer %s and %s\n",
         trname_3d[n], args[1], args[2]);
        nstatorder3d[STATFN2]++;
        nstatorder3d[STATFN5]++;
      } else
      /* Count the number of 3d tracers requiring RMSE        */
      if (strcmp(tok, "rmse") == 0) {
        trs->do_tracerstats = 1;
        stat = O_1|O_2|O_3;
        emstag(LTRACE,"tracerstats:init","statistics : %s = RMSE of tracer %s and %s\n",
         trname_3d[n], args[1], args[2]);
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][3])) {
          sprintf(tr_3d[on][3],"%u", LASTSTEP);
        }
        nstatorder3d[STATFN0]++;
      } else

      /* Count the number of 3d tracers requiring differences        */
      if (strcmp(tok, "diff") == 0) {
        trs->do_tracerstats = 1;
        stat = O_1|O_2|O_3;
        emstag(LTRACE,"tracerstats:init","statistics : %s = difference of tracer %s - %s\n",
         trname_3d[n], args[1], args[2]);
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][3])) {
          sprintf(tr_3d[on][3],"%u", LASTSTEP);
        }
        nstatorder3d[STATFN0]++;
      } else

      /* Count the number of 3d tracers requiring summing        */
      if (strcmp(tok, "sum") == 0) {
        trs->do_tracerstats = 1;
        /* stat = O_1|O_2|O_3; */
        stat = O_1 | B_2;
        emstag(LTRACE,"tracerstats:init","statistics : %s = sum of tracer %s + %s\n",
         trname_3d[n], args[1], args[2]);
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][3])) {
          sprintf(tr_3d[on][3],"%u", LASTSTEP);
        }
        trs->ntr += count_elements(args[1],(int)',');
        nstatorder3d[STATFN0]++;
      } else

      /* Count the number of 3d tracers requiring u1 fluxes          */
      if (strcmp(tok, "fluxe1") == 0) {
        trs->do_tracerstats = 1;
        trs->use_fluxe1 = 1;
        stat = O_1|O_2;
        emstag(LTRACE,"tracerstats:init","statistics : %s = 3d flux of tracer %s in the e1 direction\n",
         trname_3d[n], args[1]);
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][3])) {
          sprintf(tr_3d[on][3],"%u", LASTSTEP);
        }
        nstatorder3d[STATFN0]++;
      } else
      /* Count the number of 3d tracers requiring u2 fluxes          */
      if (strcmp(tok, "fluxe2") == 0) {
        trs->do_tracerstats = 1;
        trs->use_fluxe2 = 1;
        stat = O_1|O_2;
        emstag(LTRACE,"tracerstats:init","statistics : %s = 3d flux of tracer %s in the e2 direction\n",
             trname_3d[n], args[1]);
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][3])) {
          sprintf(tr_3d[on][3],"%u", LASTSTEP);
        }
        nstatorder3d[STATFN0]++;
      } else
      /* Count the number of 3d tracers requiring w fluxes           */
      if (strcmp(tok, "fluxw") == 0) {
        trs->do_tracerstats = 1;
        trs->use_w = 1;
        stat = O_1|O_2;
        emstag(LDEBUG,"tracerstats:init","statistics : %s = 3d flux of tracer %s in the vertical direction\n",
               trname_3d[n], args[1]);
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][3])) {
          sprintf(tr_3d[on][3],"%u", LASTSTEP);
        }
        nstatorder3d[STATFN0]++;
      } else
      /* Count the number of 3d tracers requiring mean u1 fluxes     */
      if (strcmp(tok, "meanfluxe1") == 0) {
        trs->do_tracerstats = 1;
        trs->use_fluxe1 = 1;
        stat = O_1|O_2;
        sprintf(key, "TRACER%1.1d.dt", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][2])) {
          sprintf(key, "STOP_TIME");
          prm_read_char(fp, key, tr_3d[on][2]);
        }
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][3])) {
          sprintf(tr_3d[on][3],"%u", LASTSTEP);
        }
        emstag(LTRACE,"tracerstats:init","statistics : %s = %s mean e1 flux of 3d tracer %s\n",
               trname_3d[n], tr_3d[on][2], args[1]);
        nflux++;
      } else
      /* Count the number of 3d tracers requiring mean u2 fluxes     */
      if (strcmp(tok, "meanfluxe2") == 0) {
        trs->do_tracerstats = 1;
        trs->use_fluxe2 = 1;
        stat = O_1|O_2;
        sprintf(key, "TRACER%1.1d.dt", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][2])) {
          sprintf(key, "STOP_TIME");
          prm_read_char(fp, key, tr_3d[on][2]);
        }
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][3])) {
          sprintf(tr_3d[on][3],"%u", LASTSTEP);
        }
        emstag(LDEBUG,"tracerstats:init","statistics : %s = %s mean e2 flux of 3d tracer %s\n",
               trname_3d[n], tr_3d[on][2], args[1]);
        nflux++;
      } else
      /* Count the number of 3d tracers requiring mean w fluxes      */
      if (strcmp(tok, "meanfluxw") == 0) {
        trs->do_tracerstats = 1;
        trs->use_w = 1;
        stat = O_1|O_2;
        sprintf(key, "TRACER%1.1d.dt", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][2])) {
          sprintf(key, "STOP_TIME");
          prm_read_char(fp, key, tr_3d[on][2]);
        }
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][3])) {
          sprintf(tr_3d[on][3],"%u", LASTSTEP);
        }
        emstag(LTRACE,"tracerstats:init","statistics : %s = %s mean vertical flux of 3d tracer %s\n",
               trname_3d[n], tr_3d[on][2], args[1]);
        nflux++;
      } else
      /* Count the number of 3d tracers requiring vertical profiles  */
      if (strcmp(tok, "vprof") == 0) {
        trs->do_tracerstats = 1;
        stat = O_1|O_2;
        emstag(LDEBUG,"tracerstats:init","statistics : %s = normalized vertical profile of tracer %s\n",
               trname_3d[n], args[1]);
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][3])) {
          sprintf(tr_3d[on][3],"%u", LASTSTEP);
        }
        nstatorder3d[STATFN0]++;
      } else

      /*UR-ADDED for calculating section fluxes
       */
      if (strcmp(tok, "sectionflux") == 0) {
        if(nargs < 2)
        {
          emstag(LFATAL,"tracerstats:init","Too few arguments for %s :nargs < 2! - exiting!",trname_3d[n]);
          exit(0);
        }
        trs->do_tracerstats = 1;
        stat = B_2;
        sprintf(key, "TRACER%1.1d.data", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][2])) {
          emstag(LFATAL,"tracerstats:init","Missing parameter data for %s, no geometry filename provided! - exiting!",trname_3d[n]);
          exit(1);
        }

        sprintf(key, "TRACER%1.1d.output", i_get_param_map_3d(model, n));
        strcat(tr_3d[on][2]," ");
        if (!prm_read_char(fp, key, buf)) {
          emstag(LFATAL,"tracerstats:init","Missing parameter output for %s, no output filename provided as TRACER#.output - exiting!",trname_3d[n]);
          exit(0);
        }else
             strcat(tr_3d[on][2],buf);

        sprintf(key, "TRACER%1.1d.dt", i_get_param_map_3d(model, n));
        strcat(tr_3d[on][2]," ");
        if (!prm_read_char(fp, key, buf)) {
          sprintf(tbuf,"%f",(i_get_model_timestep(model)*10));
          strcat(tr_3d[on][2],buf);
        }else
        {
          strcat(tr_3d[on][2],"'");
          strcat(tr_3d[on][2],buf);
          strcat(tr_3d[on][2],"'");
        }


        sprintf(key, "TRACER%1.1d.startt", i_get_param_map_3d(model, n));
        strcat(tr_3d[on][2]," ");
        if (!prm_read_char(fp, key, buf)) {
          sprintf(buf,"0.0");
          strcat(tr_3d[on][2],buf);
        }else
        {
          strcat(tr_3d[on][2],"'");
          strcat(tr_3d[on][2],buf);
          strcat(tr_3d[on][2],"'");
        }

        if(nargs < 3)
        {
          emstag(LTRACE,"tracerstats:init","statistics : %s = section flux of 3d tracer %s with geometry file %s",
             trname_3d[n], tr_3d[on][1],tr_3d[on][2]);
          strcat(tr_3d[on][2]," ");
          strcat(tr_3d[on][2],"1");
          trs->use_fluxe1 = 1;
        }else
        {
          strcat(tr_3d[on][2]," ");
          if(strcmp(SSECTION1,args[2]) == 0)
          {
            strcat(tr_3d[on][2],"1");
            trs->use_fluxe1 = 1;
          }else if(strcmp(SSECTION2,args[2]) == 0)
          {
            strcat(tr_3d[on][2],"2");
            trs->use_fluxe2 = 1;
          } else if(strcmp(SSECTIONV,args[2]) == 0)
          {
            strcat(tr_3d[on][2],"0");
            trs->use_w = 1;
          }else
          {
            emstag(LWARN,"tracerstats:init","Section flux direction is not of required value - %s -! Resetting to exchange through E1!",args[2]);
            strcat(tr_3d[on][2],"1");
            trs->use_fluxe1 = 1;
          }
        }

        sprintf(key, "TRACER%1.1d.tscale", i_get_param_map_3d(model, n));
        strcat(tr_3d[on][2]," ");
        if (!prm_read_char(fp, key, buf)) {
          sprintf(buf,"86400");
          strcat(tr_3d[on][2],buf);
          strcat(tr_3d[on][2]," days");
        }else
        {
          strcat(tr_3d[on][2],buf);
          sprintf(key, "TRACER%1.1d.tunit", i_get_param_map_3d(model, n));
          strcat(tr_3d[on][2]," ");
          if (!prm_read_char(fp, key, buf)) {
            emstag(LFATAL,"tracerstats:init","Missing parameter section tunit - providing a time scaling factor requires units for time to be given too! search for %s",key);
            exit(1);
          }else
          {
              strcat(tr_3d[on][2],"'");
              strcat(tr_3d[on][2],buf);
              strcat(tr_3d[on][2],"'");
          }
        }

        sprintf(key, "TRACER%1.1d.outscale", i_get_param_map_3d(model, n));

        if (prm_read_char(fp, key, buf))
        {
          strcat(tr_3d[on][2]," '");
          strcat(tr_3d[on][2],buf);
          strcat(tr_3d[on][2],"'");
        }else
        {
          strcat(tr_3d[on][2],"''");
        }

        sprintf(key, "TRACER%1.1d.outunit", i_get_param_map_3d(model, n));
        if (prm_read_char(fp, key, buf))
        {
          strcat(tr_3d[on][2]," '");
          strcat(tr_3d[on][2],buf);
          strcat(tr_3d[on][2],"'");
        }


        strcpy(operation_3d[on], args[0]);
        strcpy(tr_3d[on][0], trname_3d[n]);
        om_3d[on] = trs->ntr;

        trs->ntr += count_elements(args[1],(int)',');

        sprintf(key, "TRACER%1.1d.step", i_get_param_map_3d(model, n));
        if (!prm_read_char(fp, key, tr_3d[on][3])) {
          sprintf(tr_3d[on][3],"%u", LASTSTEP);
        }
        emstag(LDEBUG,"tracerstats:init","statistics : %s = section flux of 3d tracer: %s with geometry file, output, dt, start-time, direction, tscale, tunit, outscale (outunit) - %s attempting n tracers %d",
               trname_3d[n], args[1], tr_3d[on][2],(trs->ntr - om_3d[on]));

        ndomain++;
      } else {
	/*
	 * If we're here then we've encountered an unknown tracerstat
         *  - probably mispelling in the prm file.
	 */
	emstag(LPANIC, "tracerstats:init", "Unknown 3d tracerstat '%s' found on line '%s'. Exiting ...\n", tok, key);
	exit(1);
      }
      
      if(stat & O_1) {
        strcpy(operation_3d[on], args[0]);
        strcpy(tr_3d[on][0], trname_3d[n]);
        om_3d[on] = ntr_tmp;
        trs->ntr++;
      }

      if(stat & O_2) {
        strcpy(tr_3d[on][1], args[1]);
        trs->ntr++;
      }

      if(stat & O_3) {
        strcpy(tr_3d[on][2], args[2]);
        trs->ntr++;
      }

      if(stat & O_4) {
        strcpy(tr_3d[on][3], args[3]);
        trs->ntr++;
      }

      if(stat & B_2) {
       strcpy(tr_3d[on][1], args[1]);
      }

      if(stat & B_3)
        strcpy(tr_3d[on][2], args[2]);
      
      on++;
    }
  }

  /* Sediment tracers */
  for (n = 0; n < ntrsed; n++) {
    stat = 0;
    fluxdir = -1;
    sprintf(key, "TRACER%1.1d.tracerstat", i_get_param_map_sed(model, n));
    if (prm_read_char(fp, key, tbuf)) {
      nargs = decode_args(tbuf, args, 256);
      if (nargs < 2) {
        emstag(LFATAL,"tracerstats:init"," statistics : format <operation>(<tracer>)\n");
        exit(0);
      }
      strcpy(tok, args[0]);
      ntrsed_tmp = trs->nsed;
      
      /* Count the number of sed tracers requiring means              */
      if (strcmp(tok, "mean") == 0) {
        trs->do_tracerstats = 1;
        stat = OS_1|OS_2;
        /* Read in the averaging period and store temporarily in     */
        /* tr_sed[2].                                                 */
        sprintf(key, "TRACER%1.1d.dt", i_get_param_map_sed(model, n));
        if (!prm_read_char(fp, key, tr_sed[on_sed][2])) {
          sprintf(key, "STOP_TIME");
          prm_read_char(fp, key, tr_sed[on_sed][2]);
        }
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_sed(model, n));
        if (!prm_read_char(fp, key, tr_sed[on_sed][3])) {
          sprintf(tr_sed[on_sed][3],"%u", LASTSTEP);
        }

        emstag(LTRACE,"tracerstats:trs_init",
	       "statistics : %s = %s mean of sed tracer %s\n",
	       trname_sed[n], tr_sed[on_sed][2], args[1]);
        nstatorder_sed[STATFN1]++;

        emstag(LTRACE,"tracerstats:trs_init",
	       "statistics : %s = %s mean of sediment tracer %s\n",
	       trname_sed[n], tr_sed[on_sed][2], args[1]);
        nstatorder_sed[STATFN1]++;
      }else
      if (strcmp(tok, "run_mean") == 0) {
        trs->do_tracerstats = 1;
        stat = OS_1|OS_2;
        /* Read in the averaging period and store temporarily in     */
        /* tr_sed[2].                                                 */
        sprintf(key, "TRACER%1.1d.dt", i_get_param_map_sed(model, n));
        if (!prm_read_char(fp, key, tr_sed[on_sed][2])) {
          sprintf(key, "STOP_TIME");
          prm_read_char(fp, key, tr_sed[on_sed][2]);
        }
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_sed(model, n));
        if (!prm_read_char(fp, key, tr_sed[on_sed][3])) {
          sprintf(tr_sed[on_sed][3],"%u", LASTSTEP);
        }

        emstag(LTRACE,"tracerstats:trs_init",
	       "statistics : %s = %s run_mean of sed tracer %s\n",
	       trname_sed[n], tr_sed[on_sed][2], args[1]);
        nstatorder_sed[STATFN1]++;

        emstag(LTRACE,"tracerstats:trs_init",
	       "statistics : %s = %s run_mean of sediment tracer %s\n",
	       trname_sed[n], tr_sed[on_sed][2], args[1]);
        nstatorder_sed[STATFN1]++;
      }else

      /* Count the number of sed tracers requiring max */
      if (strcmp(tok, "max") == 0) {
        trs->do_tracerstats = 1;
        stat = OS_1|OS_2;
	/*
	 * Read in the period over which to max and store temporarily in
	 * tr_sed[2]. If dt is not present, we do it over the whole period.
	 */
        sprintf(key, "TRACER%1.1d.dt", i_get_param_map_sed(model, n));
        if (!prm_read_char(fp, key, tr_sed[on_sed][2])) {
          sprintf(key, "STOP_TIME");
          prm_read_char(fp, key, tr_sed[on_sed][2]);
        }
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_sed(model, n));
        if (!prm_read_char(fp, key, tr_sed[on_sed][3])) {
          sprintf(tr_sed[on_sed][3],"%u", LASTSTEP);
        }

        emstag(LTRACE, "tracerstats:trs_init",
	       "statistics : %s = %s max of sed tracer %s\n",
	       trname_sed[n], args[1]);
        nstatorder_sed[STATFN1]++;

        emstag(LTRACE, "tracerstats:trs_init",
	       "statistics : %s = %s max of sediment tracer %s\n",
	       trname_sed[n], args[1]);
        nstatorder_sed[STATFN1]++;
      }else

      /* Count the number of sed tracers requiring max */
      if (strcmp(tok, "min") == 0) {
        trs->do_tracerstats = 1;
        stat = OS_1|OS_2;
	/*
	 * Read in the period over which to min and store temporarily in
	 * tr_sed[2]. If dt is not present, we do it over the whole period.
	 */
        sprintf(key, "TRACER%1.1d.dt", i_get_param_map_sed(model, n));
        if (!prm_read_char(fp, key, tr_sed[on_sed][2])) {
          sprintf(key, "STOP_TIME");
          prm_read_char(fp, key, tr_sed[on_sed][2]);
        }
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_sed(model, n));
        if (!prm_read_char(fp, key, tr_sed[on_sed][3])) {
          sprintf(tr_sed[on_sed][3],"%u", LASTSTEP);
        }

        emstag(LTRACE, "tracerstats:trs_init",
	       "statistics : %s = %s min of sed tracer %s\n",
	       trname_sed[n], args[1]);
        nstatorder_sed[STATFN1]++;

        emstag(LTRACE, "tracerstats:trs_init",
	       "statistics : %s = %s min of sediment tracer %s\n",
	       trname_sed[n], args[1]);
        nstatorder_sed[STATFN1]++;
      }else

      /* Count number of sed tracers requiring standard deviations    */
      if (strcmp(tok, "stdev") == 0) {
        trs->do_tracerstats = 1;
        stat = OS_1|OS_B_2;
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_sed(model, n));
        if (!prm_read_char(fp, key, tr_sed[on_sed][3])) {
          sprintf(tr_sed[on_sed][3],"%u", LASTSTEP);
        }
        emstag(LTRACE,"tracerstats:init",
	       "statistics : %s = standard deviation of sed tracer %s\n",
	       trname_sed[n], args[1]);
        nstatorder_sed[STATFN4]++;

        emstag(LTRACE,"tracerstats:init",
	       "statistics : %s = standard deviation of sediment tracer %s\n",
	       trname_sed[n], args[1]);
        nstatorder_sed[STATFN4]++;
      } else
      /* Count the number of sed tracers requiring variances          */
      if (strcmp(tok, "variance") == 0) {
        trs->do_tracerstats = 1;
        stat = OS_1|OS_B_2;
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_sed(model, n));
        if (!prm_read_char(fp, key, tr_sed[on_sed][3])) {
          sprintf(tr_sed[on_sed][3],"%u", LASTSTEP);
        }
        emstag(LTRACE,"tracerstats:init",
	       "statistics : %s = variance of sed tracer %s\n",
	       trname_sed[n], args[1]);
        nstatorder_sed[STATFN3]++;

        emstag(LTRACE,"tracerstats:init",
	       "statistics : %s = variance of sediment tracer %s\n",
	       trname_sed[n], args[1]);
        nstatorder_sed[STATFN3]++;
      } else

      /* Count the number of sed tracers requiring summing        */
      if (strcmp(tok, "sum") == 0) {
        trs->do_tracerstats = 1;
        stat = OS_1 | OS_B_2;
        emstag(LTRACE,"tracerstats:init","statistics : %s = sum of tracer %s + %s\n",
         trname_sed[n], args[1], args[2]);
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_sed(model, n));
        if (!prm_read_char(fp, key, tr_sed[on_sed][3])) {
          sprintf(tr_sed[on_sed][3],"%u", LASTSTEP);
        }
        trs->nsed += count_elements(args[1],(int)',');
        nstatorder_sed[STATFN0]++;
      } else {
	/*
	 * If we're here then we've encountered an unknown tracerstat
	 */
	emstag(LPANIC, "tracerstats:init", "tracerstat: Unknown (or not yet implemented) sediment tracerstat '%s' found. Exiting ...\n", tok);
	exit(1);
      }

      if(stat & OS_1) {
        strcpy(operation_sed[on_sed], args[0]);
        strcpy(tr_sed[on_sed][0], trname_sed[n]);
        om_sed[on_sed] = ntrsed_tmp;
        trs->nsed++;
      }

      if(stat & OS_2) {
        strcpy(tr_sed[on_sed][1], args[1]);
        trs->nsed++;
      }

      if(stat & OS_B_2) {
       strcpy(tr_sed[on_sed][1], args[1]);
      }

      on_sed++;
    }
  }

  /* 2D tracers                                                      */
  emstag(LTRACE,"tracerstats:init"," starting 2d - %d ",(ntrS - EX_2D));
  for (n = 0; n < ntrS - EX_2D; n++) {
    stat = 0;
    sprintf(key, "TRACER%1.1d.tracerstat", i_get_param_map_2d(model, n));
    if (prm_read_char(fp, key, tbuf)) {
      nargs = decode_args(tbuf, args, 256);
      if (nargs < 2) {
       emstag(LFATAL,"tracerstats:init","statistics : format <operation>(<tracer>)\n");
       exit(0);
      }
      strcpy(tok, args[0]);
      ntr_tmp = trs->ntrS;
      /* Count the number of 2d tracers requiring means              */
      if (strcmp(tok, "mean") == 0) {
        trs->do_tracerstats = 1;
        stat = O_1|O_2;
        /* Read in the averaging period and store temporarily in     */
        /* tr_2d[2].                                                 */
        sprintf(key, "TRACER%1.1d.dt", i_get_param_map_2d(model, n));
        if (!prm_read_char(fp, key, tr_2d[onS][2])) {
          sprintf(key, "STOP_TIME");
          prm_read_char(fp, key, tr_2d[onS][2]);
        }
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_2d(model, n));
        if (!prm_read_char(fp, key, tr_2d[onS][3])) {
          sprintf(tr_2d[onS][3],"%u", LASTSTEP);
        }
        emstag(LTRACE,"tracerstats:init","statistics : %s = %s mean of 2d tracer %s\n",
             trname_2d[n], tr_2d[onS][2], args[1]);
        nstatorder2d[STATFN0]++;
      }else

      /* Count the number of 2d tracers requiring running means      */
      if (strcmp(tok, "run_mean") == 0) {
	trs->do_tracerstats = 1;
        stat = O_1|O_2;
        /* Read in the averaging period and store temporarily in     */
        /* tr_2d[2].                                                 */
        sprintf(key, "TRACER%1.1d.dt", i_get_param_map_2d(model, n));
        if (!prm_read_char(fp, key, tr_2d[onS][2])) {
          sprintf(key, "STOP_TIME");
          prm_read_char(fp, key, tr_2d[onS][2]);
        }
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_2d(model, n));
        if (!prm_read_char(fp, key, tr_2d[onS][3])) {
          sprintf(tr_2d[onS][3],"%u", LASTSTEP);
        }
        emstag(LTRACE,"tracerstats:init","statistics : %s = %s run_mean of 2d tracer %s\n",
	       trname_2d[n], tr_2d[onS][2], args[1]);
        nstatorder2d[STATFN0]++;
      }else

      /* Count the number of 2d tracers requiring max */
      if (strcmp(tok, "max") == 0) {
        trs->do_tracerstats = 1;
        stat = O_1|O_2;
        /* Read in the averaging period and store temporarily in     */
        /* tr_2d[2].                                                 */
        sprintf(key, "TRACER%1.1d.dt", i_get_param_map_2d(model, n));
        if (!prm_read_char(fp, key, tr_2d[onS][2])) {
          sprintf(key, "STOP_TIME");
          prm_read_char(fp, key, tr_2d[onS][2]);
        }
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_2d(model, n));
        if (!prm_read_char(fp, key, tr_2d[onS][3])) {
          sprintf(tr_2d[onS][3],"%u", LASTSTEP);
        }
        emstag(LTRACE,
	       "tracerstats:init","statistics : %s = %s max of 2d tracer %s\n",
	       trname_2d[n], tr_2d[onS][2], args[1]);
        nstatorder2d[STATFN0]++;
      }else

      /* Count the number of 2d tracers requiring min */
      if (strcmp(tok, "min") == 0) {
        trs->do_tracerstats = 1;
        stat = O_1|O_2;
        /* Read in the averaging period and store temporarily in     */
        /* tr_2d[2].                                                 */
        sprintf(key, "TRACER%1.1d.dt", i_get_param_map_2d(model, n));
        if (!prm_read_char(fp, key, tr_2d[onS][2])) {
          sprintf(key, "STOP_TIME");
          prm_read_char(fp, key, tr_2d[onS][2]);
        }
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_2d(model, n));
        if (!prm_read_char(fp, key, tr_2d[onS][3])) {
          sprintf(tr_2d[onS][3],"%u", LASTSTEP);
        }
        emstag(LTRACE,
	       "tracerstats:init","statistics : %s = %s min of 2d tracer %s\n",
	       trname_2d[n], tr_2d[onS][2], args[1]);
        nstatorder2d[STATFN0]++;

      }else
      /* Count the number of 2d tracers requiring RMSE               */
      if (strcmp(tok, "rmse") == 0) {
        trs->do_tracerstats = 1;
        stat = O_1|O_2|O_3;
        /* Read in the averaging period and store temporarily in     */
        /* tr_2d[2].                                                 */
        sprintf(key, "TRACER%1.1d.dt", i_get_param_map_2d(model, n));
        if (!prm_read_char(fp, key, tr_2d[onS][2])) {
          sprintf(key, "STOP_TIME");
          prm_read_char(fp, key, tr_2d[onS][2]);
        }
        emstag(LTRACE,"tracerstats:init","statistics : %s = rmse of tracer %s and %s\n",
         trname_2d[n], args[1], args[2]);
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_2d(model, n));
        if (!prm_read_char(fp, key, tr_2d[onS][3])) {
          sprintf(tr_2d[onS][3],"%u", LASTSTEP);
        }
        nstatorder2d[STATFN6]++;
      }else
      /* Count the number of 2d tracers requiring differences        */
      if (strcmp(tok, "diff") == 0) {
        trs->do_tracerstats = 1;
        stat = O_1|O_2|O_3;
        emstag(LTRACE,"tracerstats:init","statistics : %s = difference of tracer %s - %s\n",
         trname_2d[n], args[1], args[2]);
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_2d(model, n));
        if (!prm_read_char(fp, key, tr_2d[onS][3])) {
          sprintf(tr_2d[onS][3],"%u", LASTSTEP);
        }
        nstatorder2d[STATFN6]++;
      }else

      /* Count the number of 2d tracers requiring layer copies   */
      if (strcmp(tok, "copy_layer") == 0) {
        trs->do_tracerstats = 1;
        stat = O_1|OS_1|B_2;
        /* Temporarily the 3D tracer number of the tracer to         */
        /* vertically integrate in tr_2d[2] and increment the        */
        /* number of 3D tracers required.                            */
        emstag(LTRACE,"tracerstats:init","statistics : %s = layer copy of tracer %s\n", trname_2d[n], args[1]);

        /* save the index of the 3d tracer which is filled later */
        sprintf(tr_2d[onS][2], "%d", trs->ntr);
        trs->ntr++;

        emstag(LTRACE,"tracerstats:init","statistics done 1. step : %s = layer copy of tracer %s\n",
         trname_2d[n], args[1]);
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_2d(model, n));
        if (!prm_read_char(fp, key, tr_2d[onS][3])) {
          sprintf(tr_2d[onS][3],"%u", LASTSTEP);
        }
        nstatorder2d[STATFN0]++;
      } else

      /* Count the number of 2d tracers requiring summimg        */
      if (strcmp(tok, "sum") == 0) {
        trs->do_tracerstats = 1;
        stat = O_1| B_2;
        /*|O_2|O_3; */
        emstag(LTRACE,"tracerstats:init","statistics : %s = sum of tracer %s + %s\n",
         trname_2d[n], args[1], args[2]);
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_2d(model, n));
        if (!prm_read_char(fp, key, tr_2d[onS][3])) {
          sprintf(tr_2d[onS][3],"%u", LASTSTEP);
        }
        trs->ntrS += count_elements(args[1],(int)',');
        nstatorder2d[STATFN6]++;
      }else

      /* Count the number of tracers requiring vertical maximums     */
      if (strcmp(tok, "vmax") == 0) {
        trs->do_tracerstats = 1;
        stat = O_1|B_2;
        /* Temporarily the 3D tracer number of the tracer to         */
        /* vertically integrate in tr_2d[2] and increment the        */
        /* number of 3D tracers required.                            */
        emstag(LTRACE,"tracerstats:init","statistics : %s = vertical maximum of tracer %s\n", trname_2d[n], args[1]);

        /* save the index of the 3d tracer which is filled later */
        sprintf(tr_2d[onS][2], "%d", trs->ntr);
        trs->ntr++;
        emstag(LTRACE,"tracerstats:init","statistics done 1. step : %s = vertical maximum of tracer %s\n",
         trname_2d[n], args[1]);
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_2d(model, n));
        if (!prm_read_char(fp, key, tr_2d[onS][3])) {
          sprintf(tr_2d[onS][3],"%u", LASTSTEP);
        }
        nstatorder2d[STATFN0]++;
      }else

      /* Count the number of tracers requiring vertical integrals    */
      if (strcmp(tok, "vint") == 0) {
        trs->do_tracerstats = 1;
        stat = O_1|B_2;
        /* Temporarily the 3D tracer number of the tracer to         */
        /* vertically integrate in tr_2d[2] and increment the        */
        /* number of 3D tracers required.                            */
        emstag(LTRACE,"tracerstats:init","statistics : %s = vertical integral of tracer %s\n", trname_2d[n], args[1]);

        /* save the index of the 3d tracer which is filled later */
        sprintf(tr_2d[onS][2], "%d", trs->ntr);
        trs->ntr++;

        emstag(LTRACE,"tracerstats:init","statistics done 1. step : %s = vertical integral of tracer %s\n",
         trname_2d[n], args[1]);
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_2d(model, n));
        if (!prm_read_char(fp, key, tr_2d[onS][3])) {
          sprintf(tr_2d[onS][3],"%u", LASTSTEP);
        }
        nstatorder2d[STATFN0]++;
      }else

       /* Count the number of tracers requiring vertical mean    */
      if (strcmp(tok, "vmean") == 0) {
        trs->do_tracerstats = 1;
        stat = O_1|B_2;
        /* Temporarily the 3D tracer number of the tracer to         */
        /* vertically integrate in tr_2d[2] and increment the        */
        /* number of 3D tracers required.                            */
        emstag(LTRACE,"tracerstats:init","statistics : %s = vertical mean of tracer %s\n", trname_2d[n], args[1]);

        /* save the index of the 3d tracer which is filled later */
        sprintf(tr_2d[onS][2], "%d", trs->ntr);
        /* increment the 3d tracer list*/
        trs->ntr++;
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_2d(model, n));
        if (!prm_read_char(fp, key, tr_2d[onS][3])) {
          sprintf(tr_2d[onS][3],"%u", LASTSTEP);
        }
        emstag(LTRACE,"tracerstats:init","statistics done 1. step : %s = vertical mean of tracer %s\n",
        trname_2d[n], args[1]);
        nstatorder2d[STATFN0]++;
      }else

      /* for calculating vertical cross calculation
       * to establish stratification
       */
      if (strcmp(tok, "vdiff") == 0) {
        if(nargs < 4)
        {
          emstag(LFATAL,"tracerstats:init","Missing parameters for %s at tracer: %d , required are <tracer name>:<range top>:<range bottom>!",tok,i_get_param_map_3d(model, n));
          exit(0);
        }
        trs->do_tracerstats = 1;
        stat = O_1|B_2;
        emstag(LDEBUG,"tracerstats:init","statistics : %s = vertical diff of tracer %s\n", trname_2d[n], args[1]);

        /* save the index of the 3d tracer which is filled later */


        strcat(tr_2d[onS][2]," ");
        strcat(tr_2d[onS][2], args[2]);
        strcat(tr_2d[onS][2]," ");
        strcat(tr_2d[onS][2],args[3]);
        if(nargs > 4)
        {
          strcat(tr_2d[onS][2]," ");
          strcat(tr_2d[onS][2],args[4]);
        }

        if(nargs > 5)
        {
          int ii;
          strcat(tr_2d[onS][2]," ");
          strcat(tr_2d[onS][2],args[5]);
                  /* find out if we need sediment tracers*/
          ii = extract_step(model,args[5]);
          if(  ii == VDIFFSEDFLUX1 ||
               ii == VDIFFSEDFLUX2 ||
               ii == VDIFFSEDFLUX3)
          {
            sprintf(buf, "%d %s", trs->nsed,tr_2d[onS][2]);
            strcpy(tr_2d[onS][2],buf);
            trs->nsed++;
            if(!trs->sednz) 
    				{
    					emstag(LPANIC,"tracerstats:tracerstats:init","Attempting to calculate vertial diff between sediment and wc,but have no sediment layers! Exiting ...");
    					exit(1);
    				}
            
          }else
          {
            sprintf(buf, "%d %s", trs->ntr,tr_2d[onS][2]);
            strcpy(tr_2d[onS][2],buf);
            trs->ntr++;
          }
        }else
        {
            sprintf(buf, "%d %s", trs->ntr,tr_2d[onS][2]);
            strcpy(tr_2d[onS][2],buf);
        }

        if(nargs > 6)
        {
          strcat(tr_2d[onS][2]," ");
          strcat(tr_2d[onS][2],args[6]);
        }

        sprintf(key, "TRACER%1.1d.dt", i_get_param_map_2d(model, n));
        if (prm_read_char(fp, key, buf))
        {
          strcat(tr_2d[onS][2]," ");
          strcat(tr_2d[onS][2],buf);
        }

        sprintf(key, "TRACER%1.1d.step", i_get_param_map_2d(model, n));
        if (!prm_read_char(fp, key, tr_2d[onS][3])) {
          sprintf(tr_2d[onS][3],"%u", LASTSTEP);
        }
        emstag(LDEBUG,"tracerstats:init","statistics : %s = vertical diff of 3d: tracer %s with layers t/b %s/%s - strict: %s",
               trname_2d[n],  args[1],args[2],args[3],(nargs>4?args[4]:"0"));

        nstatorder2d[STATFN0]++;
      }else

      /* Count number of 2d tracers requiring standard deviations    */
      if (strcmp(tok, "stdev") == 0) {
        trs->do_tracerstats = 1;
        stat = O_1|B_2;
        emstag(LTRACE,"tracerstats:init","statistics : %s = standard deviation of 2d tracer %s\n",
         trname_2d[n], args[1]);
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_2d(model, n));
        if (!prm_read_char(fp, key, tr_2d[onS][3])) {
          sprintf(tr_2d[onS][3],"%u", LASTSTEP);
        }
        nstatorder2d[STATFN3]++;
      }else

      /* Count the number of 2d tracers requiring variances          */
      if (strcmp(tok, "variance") == 0) {
        trs->do_tracerstats = 1;
        stat = O_1|B_2;
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_2d(model, n));
        if (!prm_read_char(fp, key, tr_2d[onS][3])) {
          sprintf(tr_2d[onS][3],"%u", LASTSTEP);
        }
        emstag(LTRACE,"tracerstats:init","statistics : %s = variance of 2d tracer %s\n",
         trname_2d[n], args[1]);
        nstatorder2d[STATFN2]++;
      }else

      /* Count the number of 2d tracers requiring covariances        */
      if (strcmp(tok, "cov") == 0) {
        trs->do_tracerstats = 1;
        nmS = max(nmS, 2);
        stat = O_1|B_2|B_3;
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_2d(model, n));
        if (!prm_read_char(fp, key, tr_2d[onS][3])) {
          sprintf(tr_2d[onS][3],"%u", LASTSTEP);
        }
        emstag(LTRACE,"tracerstats:init","statistics : %s = covariance of 2d tracer %s and %s\n",
         trname_2d[n], args[1], args[2]);
        nstatorder2d[STATFN5]++;
      }else

      /* Count the number of 2d tracers requiring correlation        */
      /* coefficients.                                               */
      if (strcmp(tok, "corr") == 0) {
        trs->do_tracerstats = 1;
        nmS = max(nmS, 4);
        stat = O_1|B_2|B_3;
        sprintf(key, "TRACER%1.1d.step", i_get_param_map_2d(model, n));
        if (!prm_read_char(fp, key, tr_2d[onS][3])) {
          sprintf(tr_2d[onS][3],"%u", LASTSTEP);
        }
        emstag(LTRACE,"tracerstats:init","statistics : %s = correlation coefficient of 2d tracer %s and %s\n",
         trname_2d[n], args[1], args[2]);
        nstatorder2d[STATFN1]++;
        nstatorder2d[STATFN4]++;
      }

      if(stat & O_1) {
        strcpy(operation_2d[onS], args[0]);
        strcpy(tr_2d[onS][0], trname_2d[n]);
        om_2d[onS] = ntr_tmp;
        /* trs->ntrS; */
        trs->ntrS++;
      }
      if(stat & O_2) {
        strcpy(tr_2d[onS][1], args[1]);
        trs->ntrS++;
      }
      if(stat & O_3) {
        strcpy(tr_2d[onS][2], args[2]);
        trs->ntrS++;
      }
      if(stat & B_2)
         strcpy(tr_2d[onS][1], args[1]);

      if(stat & B_3)
        strcpy(tr_2d[onS][2], args[2]);

      if(stat & OS_1)
        strcpy(tr_2d[onS][4], args[2]);
      onS++;
    }
  }

/* Now go through the Non-Tracer stats for the once which do not
 * require a one or two dimensional tracer as output
 * Note: this doesn't make sense for all tracers. */

if(i > 0)
{
  for (n = 0; n < i; n++)
  {
    stat = 0;
    sprintf(key, "RTSTAT%1.1d.name", n);
    if (!prm_read_char(fp, key, tr_3d[on][0]))
    {
      emstag(LWARN,"tracerstats:init"," statistics : format RTSTAT%d.name is missing, require unique name, ignoring... ",n);
      continue;
    }

    sprintf(key, "RTSTAT%1.1d.func", n);
    if (prm_read_char(fp, key, tbuf))
    {
      nargs = decode_args(tbuf, args, 256);
      if (nargs < 2) {
        emstag(LFATAL,"tracerstats:init"," statistics : format RTSTAT#.function <operation>(<tracer>)\n");
        exit(0);
      }
      strcpy(tok, args[0]);
      if (strcmp(tok, "sectionflux") == 0)
      {
        if(nargs < 2)
        {
            emstag(LFATAL,"tracerstats:init","Too few arguments for %s :nargs < 3! - exiting!",trname_3d[n]);
            exit(0);
        }
        trs->do_tracerstats = 1;
        stat = B_2;
        sprintf(key, "RTSTAT%1.1d.data", n);
        if (!prm_read_char(fp, key, tr_3d[on][2])) {
          emstag(LFATAL,"tracerstats:init","Missing parameter data for %s, no geometry filename provided! - exiting!",trname_3d[n]);
          exit(1);
        }

        sprintf(key, "RTSTAT%1.1d.output",n);
        strcat(tr_3d[on][2]," ");
        if (!prm_read_char(fp, key, buf)) {
          emstag(LFATAL,"tracerstats:init","Missing parameter output for %s, no output filename provided as TRACER#.output - exiting!",trname_3d[n]);
          exit(0);
        }else
            strcat(tr_3d[on][2],buf);

        sprintf(key, "RTSTAT%1.1d.dt", n);
        strcat(tr_3d[on][2]," ");
        if (!prm_read_char(fp, key, buf)) {
          sprintf(tbuf,"%f",(i_get_model_timestep(model)*10));
          strcat(tr_3d[on][2],buf);
        }else
        {
          strcat(tr_3d[on][2],"'");
          strcat(tr_3d[on][2],buf);
          strcat(tr_3d[on][2],"'");
        }


        sprintf(key, "RTSTAT%1.1d.startt", n);
        strcat(tr_3d[on][2]," ");
        if (!prm_read_char(fp, key, buf)) {
          sprintf(buf,"0.0");
          strcat(tr_3d[on][2],buf);
        }else
        {
          strcat(tr_3d[on][2],"'");
          strcat(tr_3d[on][2],buf);
          strcat(tr_3d[on][2],"'");
        }

        /* copy the target tracer name */
        strcpy(tr_3d[on][1], args[1]);

        if(nargs < 3)
        {
          emstag(LTRACE,"tracerstats:init","statistics : %s = section flux of 3d tracer %s with geometry file %s",
             trname_3d[n], tr_3d[on][1],tr_3d[on][2]);
          strcat(tr_3d[on][2]," ");
          strcat(tr_3d[on][2],"1");
          trs->use_fluxe1 = 1;
        }else
        {
          strcat(tr_3d[on][2]," ");
          if(strcmp(SSECTION1,args[2]) == 0)
          {
            strcat(tr_3d[on][2],"1");
            trs->use_fluxe1 = 1;
          }else if(strcmp(SSECTION2,args[2]) == 0)
          {
            strcat(tr_3d[on][2],"2");
            trs->use_fluxe2 = 1;
          } else if(strcmp(SSECTIONV,args[2]) == 0)
          {
            strcat(tr_3d[on][2],"0");
            trs->use_w = 1;
          }else
          {
	    // This warning may be redundant if direction is specified
	    // on a per cell basis
            emstag(LWARN,"tracerstats:init","Section flux direction is not of required value - %s -! Resetting to exchange through E1!",args[2]);
            strcat(tr_3d[on][2],"1");
            trs->use_fluxe1 = 1;
          }
        }

        sprintf(key, "RTSTAT%1.1d.tscale",n);
        strcat(tr_3d[on][2]," ");
        if (!prm_read_char(fp, key, buf)) {
          sprintf(buf,"86400");
          strcat(tr_3d[on][2],buf);
          strcat(tr_3d[on][2]," days");
        }else
        {
          strcat(tr_3d[on][2],buf);
          sprintf(key, "RTSTAT%1.1d.tunit", n);
          strcat(tr_3d[on][2]," ");
          if (!prm_read_char(fp, key, buf)) {
            emstag(LFATAL,"tracerstats:init","Missing parameter section tunit - providing a time scaling factor requires units for time to be given too! search for %s",key);
            exit(1);
          }else
          {
              strcat(tr_3d[on][2],"'");
              strcat(tr_3d[on][2],buf);
              strcat(tr_3d[on][2],"'");
          }
        }

	/* param 8 */
	sprintf(key, "RTSTAT%1.1d.mode", n);
        if (prm_read_char(fp, key, buf)) {
	  // Place the value
          strcat(tr_3d[on][2], " ");
          strcat(tr_3d[on][2], buf);
          strcat(tr_3d[on][2], " ");
	} else {	    
	  // Default of integrate
	  strcat(tr_3d[on][2], " integrate ");
	}
	
        /* param 9 */
        sprintf(key, "RTSTAT%1.1d.outscale", n);
        if (prm_read_char(fp, key, buf))
        {
          strcat(tr_3d[on][2]," ");
          strcat(tr_3d[on][2],"'");
          strcat(tr_3d[on][2],buf);
          strcat(tr_3d[on][2],"'");
        }else
        {
          strcat(tr_3d[on][2],"''");
        }

        /* param 10 */
        sprintf(key, "RTSTAT%1.1d.outunit", n);
        if (prm_read_char(fp, key, buf))
        {
          strcat(tr_3d[on][2]," ");
          strcat(tr_3d[on][2],"'");
          strcat(tr_3d[on][2],buf);
          strcat(tr_3d[on][2],"'");
        }

        strcpy(operation_3d[on], args[0]);
        om_3d[on] = trs->ntr;
        trs->ntr += count_elements(args[1],(int)',');

        sprintf(key, "TRACER%1.1d.step", i_get_param_map_2d(model, n));
        sprintf(key, "RTSTAT%1.1d.step", n);
        if (!prm_read_char(fp, key, buf)){
          sprintf(tr_3d[on][3],"%u", LASTSTEP);
        }

        emstag(LDEBUG,"tracerstats:init","RTSTAT statistics : %s = section flux of 3d tracer: %s with geometry file, output, dt, start-time, direction, tscale, tunit, outscale (outunit) - %s attempting n tracers %d",
               tr_3d[on][0], args[1], tr_3d[on][2], (trs->ntr - om_3d[on]));

        ndomain++;
      }
      /* Add more here ...
       */



      if(stat & O_1) {
        strcpy(operation_3d[on], args[0]);
        om_3d[on] = trs->ntr;
        trs->ntr++;
      }


      if(stat & B_2) {
       strcpy(tr_3d[on][1], args[1]);
      }

      if(stat & B_3)
        strcpy(tr_3d[on][2], args[2]);


      on++;
    }else
      emstag(LWARN,"tracerstats:init","No Function defined for RTSTAT entry: %d",n);
  }
}
/* End RTSTAT */
  if(!trs->do_tracerstats)
    return;

  /*-----------------------------------------------------------------*/
  /* Set up the required 3d tracers.                                 */
  /* Allocate memory for 3D tracerstat arrays.                       */
  if(trs->ntr) {
    trs->trname_3d = malloc(trs->ntr * sizeof(char *));
    for (n = 0; n < trs->ntr; n++) {
      trs->trname_3d[n] = malloc(MAXSTRLEN * sizeof(char));
    }
    trs->stat_type = i_alloc_1d(trs->ntr);
    trs->w_wc.w1 = d_alloc_1d(trs->ntr);
    trs->w_wc.w2 = d_alloc_1d(trs->ntr);
    trs->w_wc.w3 = d_alloc_1d(trs->ntr);
    trs->w_wc.w4 = i_alloc_1d(trs->ntr);
    trs->w_wc.w5 = i_alloc_1d(trs->ntr);
    trs->w_wc.w6 = d_alloc_1d(trs->ntr);
    memset(trs->stat_type, 0, trs->ntr * sizeof(int));
    memset(trs->w_wc.w4, 0, trs->ntr * sizeof(int));
    trs->map_3d = i_alloc_2d(nm, trs->ntr);
  }

  /* 
   * setup sediment tracers if required
   */
  if(trs->nsed && trs->sednz) {
    int nt = trs->nsed;
    trs->trname_sed = malloc(nt * sizeof(char *));
    for (n = 0; n < nt; n++) {
      trs->trname_sed[n] = malloc(MAXSTRLEN * sizeof(char));
    }
    trs->stat_type_sed = i_alloc_1d(nt);
    trs->w_sed.w1 = d_alloc_1d(nt);
    trs->w_sed.w2 = d_alloc_1d(nt);
    trs->w_sed.w3 = d_alloc_1d(nt);
    trs->w_sed.w4 = i_alloc_1d(nt);
    memset(trs->stat_type_sed, 0, nt * sizeof(int));
    memset(trs->w_sed.w4, 0, nt * sizeof(int));
    trs->map_sed = i_alloc_2d(nm_sed, nt);
    /* assign memory for sed layer thicknesses */
    trs->dz_sed = d_alloc_1d(trs->sednz);
  }

  trs->nstep = i_alloc_1d(trs->ntr+trs->ntrS+trs->nsed);
  memset(trs->nstep, 0, (trs->ntr+trs->ntrS+trs->nsed) * sizeof(int));

  /*UR moved this up here
   */
  if(is_log_enabled(LDEBUG))
  {
    trs->use_eta = 1;
    trs->use_Kz = 1;
    trs->use_Vz = 1;
    trs->use_fluxe1 = 1;
    trs->use_fluxe2 = 1;
    trs->use_w = 1;
  }

  if(trs->use_Kz)
    trs->Kz = d_alloc_1d(trs->nz);
  if(trs->use_Vz)
    trs->Vz = d_alloc_1d(trs->nz);
  if(trs->use_fluxe1)
    trs->u1flux3d = d_alloc_1d(trs->nz);
  if(trs->use_fluxe2)
    trs->u2flux3d = d_alloc_1d(trs->nz);
  if(trs->use_w)
    trs->w = d_alloc_1d(trs->nz);

  // There are some more at the end of this function as they use flags
  // set by the sectionflux code below
  
  /*
   * Assign how many statfunctions and flux functions we have
   */
  trs->nfluxfn = 0;
  trs->nstatfn3d   = MAXSTATFN;
  trs->nstatfn_sed = MAXSTATFN;
  trs->nstatfn2d   = MAXSTATFN;

  /*
   * Initialise  the function arrays
   */
  trs->statfn3d   = (lstatfn3d_t*) calloc(trs->nstatfn3d, sizeof(lstatfn3d_t));
  trs->statfn_sed = (lstatfn3d_t*) calloc(trs->nstatfn3d, sizeof(lstatfn3d_t));
  trs->statfn2d   = (lstatfn2d_t*) calloc(trs->nstatfn2d, sizeof(lstatfn2d_t));
  if(nflux)
    trs->fluxfn = (fluxfn_t*) malloc(nflux * sizeof(fluxfn_t));

  if(ndomain)
    trs->domainfn3d = (domainfn3d_t*) malloc(ndomain * sizeof(domainfn3d_t));

  for (n = 0; n < trs->nstatfn3d; n++)
  {
    /* trs->statfn[n] = *(lstatfn_t*)malloc(sizeof(lstatfn_t)); */
    if(nstatorder3d[n])
    {
      trs->statfn3d[n].statfn = calloc(nstatorder3d[n],sizeof(statfn3d_t));
      trs->statfn3d[n].nfn =0;
    }
    emstag(LDEBUG,"tracerstats:init","Adding %d at 3d-STATFN%d",nstatorder3d[n],n);

    if(nstatorder_sed[n])
    {
      trs->statfn_sed[n].statfn = calloc(nstatorder_sed[n],sizeof(statfn3d_t));
      trs->statfn_sed[n].nfn =0;
    }
    emstag(LDEBUG,"tracerstats:init","Adding %d at sed-STATFN%d",nstatorder_sed[n],n);

    if(nstatorder2d[n])
    {
      trs->statfn2d[n].statfn = calloc(nstatorder2d[n],sizeof(statfn2d_t));
      trs->statfn2d[n].nfn=0;
    }
    emstag(LDEBUG,"tracerstats:init","Adding %d at 2d-STATFN%d",nstatorder2d[n],n);
  }

  emstag(LDEBUG,"tracerstats:init","Adding %d at FLUXFN",nflux);
  emstag(LDEBUG,"tracerstats:init","Adding %d at DOMAINFN",ndomain);

  /* Initialise 3D WC tracerstat arrays */
  for (n = 0; n < on; n++) {
    lstatfn3d_t* statfn3d = &trs->statfn3d[0];

    m = om_3d[n];
    fill_stat3dfcn(trs, operation_3d[n], tr_3d[n],
		   statfn3d, trs->trname_3d, m, &trs->w_wc, trs->stat_type);
    /*
     * Once we have sediment versions of the ones below working, we
     * can move them up into the more general function above
     */
    if (strcmp(operation_3d[n], "sum@dt") == 0) 
      {
	tm_scale_to_secs(tr_3d[n][2], &trs->w_wc.w1[m]);
	trs->w_wc.w2[m] = 0.0;
	trs->stat_type[m] = SUM_DT;
	strcpy(trs->trname_3d[m],tr_3d[n][0]);
	strcpy(trs->trname_3d[m + 1],tr_3d[n][1]);
	statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].type = SUM_DT;
	statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].statfn = tracer_sum_dt_3d;
	statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].statfn_new = NULL;
	statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].n = m;
	statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].data = NULL;
	statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].step = extract_step(trs->model, tr_3d[n][3]);
	trs->nprestep[statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].step]++;
	
	statfn3d[STATFN1].nfn++;
      } 
    else if (strcmp(operation_3d[n], "dhw") == 0) 
      {
	tm_scale_to_secs(tr_3d[n][2], &trs->w_wc.w1[m]);
	trs->w_wc.w2[m] = 0.0;
	trs->stat_type[m] = DHW;
	strcpy(trs->trname_3d[m],tr_3d[n][0]);
	strcpy(trs->trname_3d[m + 1],tr_3d[n][1]);
	strcpy(trs->trname_3d[m + 2],tr_3d[n][2]);
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].type = DHW;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].statfn = tracer_dhw;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].n = m;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].data = NULL;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].step = extract_step(trs->model, tr_3d[n][3]);
	trs->nprestep[statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].step]++;
	  
	statfn3d[STATFN0].nfn++;
	  
      } 
    else if (strcmp(operation_3d[n], "exposure") == 0) 
      {
	int j;
	double d1, dt, newt = trs->time;
	char str[MAXSTRLEN];
	sprintf(str, "seconds %s", tr_3d[n][7]);
	tm_scale_to_secs(tr_3d[n][8], &trs->w_wc.w6[m]);
	tm_scale_to_secs(tr_3d[n][5], &trs->w_wc.w1[m]);
	tm_scale_to_secs(tr_3d[n][6], &trs->w_wc.w2[m]);
	/* Change model time to match the start time units */
	tm_change_time_units(i_get_model_timeunit(model),
			     str, &newt, 1);
	i = (int)((newt - trs->w_wc.w2[m]) / trs->w_wc.w1[m]);
	/* Note: on old formulation uses a start time, .w2[m], */
	/* so that exposure can be used with restarts (the     */
	/* function tracer_exposureo() in statistics.c).       */
	/* A newer formulation (tracer_exposure() does not     */
	/* require or use this variable.                       */
	trs->w_wc.w2[m] += (double)i * trs->w_wc.w1[m];
	trs->stat_type[m] = EXPOSURE;
	strcpy(trs->trname_3d[m],tr_3d[n][0]);
	strcpy(trs->trname_3d[m + 1],tr_3d[n][1]);
	strcpy(trs->trname_3d[m + 2],tr_3d[n][2]);
	strcpy(trs->trname_3d[m + 3],tr_3d[n][3]);
	trs->w_wc.w5[m] = 0;
	if (trs->trname_3d[m + 2][0] == '+')
	  trs->w_wc.w5[m] |= EX_GT;
	else if (trs->trname_3d[m + 2][0] == '-')
	  trs->w_wc.w5[m] |= EX_LT;
	if (trs->w_wc.w5[m]) {
	  strcpy(buf,trs->trname_3d[m + 2]);
	  for (i = 1; i < strlen(buf); i++)
	    trs->trname_3d[m + 2][i-1] = buf[i];
	  trs->trname_3d[m + 2][i-1] = '\0';
	} else
	  trs->w_wc.w5[m] |= EX_GT;
	if (i_tracername_exists(model, trs->trname_3d[m + 2])) {
	  trs->w_wc.w5[m] |= EX_TR;
	} else {
	  trs->w_wc.w5[m] |= EX_VA;
 	  trs->w_wc.w3[m] = atof(trs->trname_3d[m + 2]);
	}
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].type = EXPOSURE;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].statfn = tracer_exposure;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].statfn_new = NULL;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].n = m;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].data = NULL;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].step = extract_step(trs->model, tr_3d[n][3]);
	trs->nprestep[statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].step]++;
	  
	statfn3d[STATFN0].nfn++;
	  
      } 
    else if (strcmp(operation_3d[n], "reeftemp") == 0) 
      {
	int j;
	double d1, dt, newt = trs->time;
	char str[MAXSTRLEN];
	tm_scale_to_secs(tr_3d[n][4], &trs->w_wc.w6[m]);
	trs->w_wc.w2[m] = tm_time_to_julsecs(i_get_model_timeunit(model));
	trs->stat_type[m] = REEFTEMP;
	strcpy(trs->trname_3d[m],tr_3d[n][0]);
	strcpy(trs->trname_3d[m + 1],tr_3d[n][1]);
	strcpy(trs->trname_3d[m + 2],tr_3d[n][2]);
	trs->w_wc.w5[m] = 0;
	if (i_tracername_exists(model, trs->trname_3d[m + 2])) {
	  trs->w_wc.w5[m] |= EX_TR;
	} else {
	  trs->w_wc.w5[m] |= EX_VA;
 	  trs->w_wc.w3[m] = atof(trs->trname_3d[m + 2]);
	}
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].type = REEFTEMP;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].statfn = tracer_reeftemp;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].statfn_new = NULL;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].n = m;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].data = NULL;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].step = extract_step(trs->model, tr_3d[n][3]);
	trs->nprestep[statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].step]++;
	  
	statfn3d[STATFN0].nfn++;
	  
      } 
    else if (strcmp(operation_3d[n], "cov") == 0) 
      {
	trs->stat_type[m] = COV;
	strcpy(trs->trname_3d[m],tr_3d[n][0]);
	statfn3d[STATFN6].statfn[statfn3d[STATFN6].nfn].type = COV;
	statfn3d[STATFN6].statfn[statfn3d[STATFN6].nfn].statfn = tracer_cov_3d;
	statfn3d[STATFN6].statfn[statfn3d[STATFN6].nfn].statfn_new = NULL;
	statfn3d[STATFN6].statfn[statfn3d[STATFN6].nfn].n = m;
	statfn3d[STATFN6].statfn[statfn3d[STATFN6].nfn].data = NULL;
	statfn3d[STATFN6].statfn[statfn3d[STATFN6].nfn].step = extract_step(trs->model, tr_3d[n][3]);
	trs->nprestep[statfn3d[STATFN6].statfn[statfn3d[STATFN6].nfn].step]++;
	  
	statfn3d[STATFN6].nfn++;
      } 
    else if (strcmp(operation_3d[n], "corr") == 0) 
      {
	trs->stat_type[m] = CORR;
	strcpy(trs->trname_3d[m],tr_3d[n][0]);
	  
	statfn3d[STATFN2].statfn[statfn3d[STATFN2].nfn].type = CORR;
	statfn3d[STATFN2].statfn[statfn3d[STATFN2].nfn].statfn = corr_3d_pre;
	statfn3d[STATFN2].statfn[statfn3d[STATFN2].nfn].statfn_new = NULL;
	statfn3d[STATFN2].statfn[statfn3d[STATFN2].nfn].n = m;
	statfn3d[STATFN2].statfn[statfn3d[STATFN2].nfn].data = NULL;
	statfn3d[STATFN2].statfn[statfn3d[STATFN2].nfn].step = extract_step(trs->model, tr_3d[n][3]);
	trs->nprestep[statfn3d[STATFN2].statfn[statfn3d[STATFN2].nfn].step]++;
	  
	statfn3d[STATFN5].statfn[statfn3d[STATFN5].nfn].type = CORR;
	statfn3d[STATFN5].statfn[statfn3d[STATFN5].nfn].statfn = corr_3d_post;
	statfn3d[STATFN5].statfn[statfn3d[STATFN5].nfn].statfn_new = NULL;
	statfn3d[STATFN5].statfn[statfn3d[STATFN5].nfn].n = m;
	statfn3d[STATFN5].statfn[statfn3d[STATFN5].nfn].data = NULL;
	statfn3d[STATFN5].statfn[statfn3d[STATFN5].nfn].step = extract_step(trs->model, tr_3d[n][3]);
	  
	trs->nprestep[statfn3d[STATFN5].statfn[statfn3d[STATFN5].nfn].step]++;
	  
	statfn3d[STATFN2].nfn++;
	statfn3d[STATFN5].nfn++;
      } 
    else if (strcmp(operation_3d[n], "diff") == 0) 
      {
	trs->stat_type[m] = DIFF;
	strcpy(trs->trname_3d[m],tr_3d[n][0]);
	strcpy(trs->trname_3d[m + 1],tr_3d[n][1]);
	strcpy(trs->trname_3d[m + 2],tr_3d[n][2]);
	  
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].type = DIFF;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].statfn = tracer_diff_3d;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].statfn_new = NULL;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].n = m;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].data = NULL;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].step = extract_step(trs->model, tr_3d[n][3]);
	trs->nprestep[statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].step]++;
	
	statfn3d[STATFN0].nfn++;
      } 
    else if (strcmp(operation_3d[n], "fluxe1") == 0) 
      {
	trs->stat_type[m] = FLUXE1;
	strcpy(trs->trname_3d[m],tr_3d[n][0]);
	strcpy(trs->trname_3d[m + 1],tr_3d[n][1]);
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].type = FLUXE1;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].statfn = tracer_flux3d_e1;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].statfn_new = NULL;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].n = m;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].data = NULL;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].step = extract_step(trs->model, tr_3d[n][3]);
	
	trs->nprestep[statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].step]++;
	statfn3d[STATFN0].nfn++;
      }
    else if (strcmp(operation_3d[n], "fluxe2") == 0) 
      {
	trs->stat_type[m] = FLUXE2;
	strcpy(trs->trname_3d[m],tr_3d[n][0]);
	strcpy(trs->trname_3d[m + 1],tr_3d[n][1]);
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].type = FLUXE2;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].statfn = tracer_flux3d_e2;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].statfn_new = NULL;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].n = m;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].data = NULL;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].step = extract_step(trs->model, tr_3d[n][3]);
	
	trs->nprestep[statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].step]++;
	statfn3d[STATFN0].nfn++;
      } 
    else if (strcmp(operation_3d[n], "fluxw") == 0) 
      {
	trs->stat_type[m] = FLUXW;
	strcpy(trs->trname_3d[m],tr_3d[n][0]);
	strcpy(trs->trname_3d[m + 1],tr_3d[n][1]);
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].type = FLUXW;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].statfn = tracer_flux3d_w;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].statfn_new = NULL;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].n = m;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].data = NULL;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].step = extract_step(trs->model, tr_3d[n][3]);
	
	trs->nprestep[statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].step]++;
	
	statfn3d[STATFN0].nfn++;
      } 
    else if (strcmp(operation_3d[n], "meanfluxe1") == 0) 
      {
	trs->stat_type[m] = MEANFLUXE1;
	strcpy(trs->trname_3d[m],tr_3d[n][0]);
	strcpy(trs->trname_3d[m + 1],tr_3d[n][1]);
	tm_scale_to_secs(tr_3d[n][2], &trs->w_wc.w1[m]);
	trs->w_wc.w2[m] = 0.0;
	
	trs->fluxfn[trs->nfluxfn].type = MEANFLUXE1;
	trs->fluxfn[trs->nfluxfn].fluxfn = tracer_meanflux_3d;
	trs->fluxfn[trs->nfluxfn].data = trs->tr_e1wc;
	trs->fluxfn[trs->nfluxfn].flux = trs->u1flux3d;
	trs->fluxfn[trs->nfluxfn].n = m;
	trs->fluxfn[trs->nfluxfn].step = extract_step(trs->model, tr_3d[n][3]);
	trs->nprestep[trs->fluxfn[trs->nfluxfn].step]++;
	trs->nfluxfn++;
	
	
	emstag(LDEBUG,"tracerstats:init","Initialised : %s = mean e1 flux of 3d tracer %s at step %u",trs->trname_3d[m],trs->trname_3d[m+1],trs->fluxfn[trs->nfluxfn].step);
      } 
    else if (strcmp(operation_3d[n], "meanfluxe2") == 0) 
      {
	trs->stat_type[m] = MEANFLUXE2;
	strcpy(trs->trname_3d[m],tr_3d[n][0]);
	strcpy(trs->trname_3d[m + 1],tr_3d[n][1]);
	tm_scale_to_secs(tr_3d[n][2], &trs->w_wc.w1[m]);
	
	trs->w_wc.w2[m] = 0.0;
	trs->fluxfn[trs->nfluxfn].type = MEANFLUXE2;
	trs->fluxfn[trs->nfluxfn].fluxfn = tracer_meanflux_3d;
	trs->fluxfn[trs->nfluxfn].data = trs->tr_e2wc;
	trs->fluxfn[trs->nfluxfn].flux = trs->u2flux3d;
	trs->fluxfn[trs->nfluxfn].n = m;
	trs->fluxfn[trs->nfluxfn].step = extract_step(trs->model, tr_3d[n][3]);
	trs->nprestep[trs->fluxfn[trs->nfluxfn].step]++;
	
	trs->nfluxfn++;
	emstag(LDEBUG,"tracerstats:init","Initialised : %s = mean e2 flux of 3d tracer %s at step %u",trs->trname_3d[m],trs->trname_3d[m+1],trs->fluxfn[trs->nfluxfn].step);
	
      } 
    else if (strcmp(operation_3d[n], "meanfluxw") == 0) 
      {
	trs->stat_type[m] = MEANFLUXW;
	strcpy(trs->trname_3d[m],tr_3d[n][0]);
	strcpy(trs->trname_3d[m + 1],tr_3d[n][1]);
	tm_scale_to_secs(tr_3d[n][2], &trs->w_wc.w1[m]);
	trs->w_wc.w2[m] = 0.0;
	
	trs->fluxfn[trs->nfluxfn].type = MEANFLUXW;
	trs->fluxfn[trs->nfluxfn].fluxfn = tracer_meanflux_3d;
	trs->fluxfn[trs->nfluxfn].data = NULL;
	trs->fluxfn[trs->nfluxfn].flux = trs->w;
	trs->fluxfn[trs->nfluxfn].n = m;
	trs->fluxfn[trs->nfluxfn].step = extract_step(trs->model, tr_3d[n][3]);
	trs->nprestep[trs->fluxfn[trs->nfluxfn].step]++;
	
	/* mustbe last */
	trs->nfluxfn++;
	
      } 
    else if (strcmp(operation_3d[n], "vprof") == 0) 
      {
	trs->stat_type[m] = VPROF;
	strcpy(trs->trname_3d[m],tr_3d[n][0]);
	strcpy(trs->trname_3d[m + 1],tr_3d[n][1]);
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].type = VPROF;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].statfn = tracer_vprof;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].statfn_new = NULL;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].n = m;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].data = NULL;
	statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].step = extract_step(trs->model, tr_3d[n][3]);
	
	trs->nprestep[statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].step]++;
	
	statfn3d[STATFN0].nfn++;
      } 
    else if (strcmp(operation_3d[n], "sectionflux") == 0) 
      {
	trs->stat_type[m] = SECTFLUX;
	/* assign the name of the target tracer directly,
	 * we don't need this tracer if assigned */
	
	i = parseline(tr_3d[n][2],pbuf,15);
	trs->domainfn3d[trs->ndomainfn3d].domainfn=tracer_section;
	trs->domainfn3d[trs->ndomainfn3d].type = SECTFLUX;
	trs->domainfn3d[trs->ndomainfn3d].n = m;
	trs->domainfn3d[trs->ndomainfn3d].step = extract_step(trs->model, tr_3d[n][3]);
	trs->nprestep[trs->domainfn3d[trs->ndomainfn3d].step]++;
	emstag(LDEBUG,"Tracerstats:init","Section flux, setting up memory for tracer %s and parameters %s ",tr_3d[n][1], tr_3d[n][2]);
	trs->domainfn3d[trs->ndomainfn3d].data = malloc(sizeof(sectioncoord_t));
	sectiondata = (sectioncoord_t*)trs->domainfn3d[trs->ndomainfn3d].data;
	
	if(trs->domainfn3d[trs->ndomainfn3d].data == NULL) {
	  emstag(LERROR,"Tracerstats:init","Section flux, could not assign memory, exiting!");
	  exit(0);
	}
	
	sectiondata->name = (char*) malloc(strlen(tr_3d[n][0])+1);
	strcpy(sectiondata->name,tr_3d[n][0]);
	// Get the mode
	if (strcasecmp(pbuf[7], "average") == 0) {
	  sectiondata->mode = SECTION_AVERAGE;
	} else if (strcasecmp(pbuf[7], "integrate") == 0) {
	  sectiondata->mode = SECTION_INTEGRATE;
	} else {
	  emstag(LERROR,"Tracerstats:init","Section flux init. Unknown section flux mode '%s'",pbuf[7]);
	  exit(1);
	}

	init_sectiondata(sectiondata,tr_3d[n][1],(i > 8 && is_valid(pbuf[8])?pbuf[8]:NULL),trs->trname_3d,m);
	emstag(LDEBUG,"Tracerstats:init","Section flux, extracted %d tracer, first %s ",sectiondata->ntrs, sectiondata->tnames[0]);

	// set area flag
	if (sectiondata->mode == SECTION_AVERAGE) {
	  trs->use_area_e1_e2 = 1;
	  trs->use_u1 = trs->use_u2 = 1;
	}

	sectiondata->fileout = malloc(sizeof(char) * strlen(pbuf[1])+1);
	strcpy(sectiondata->fileout,pbuf[1]);
	
	tm_scale_to_secs(pbuf[2], &(sectiondata->dt));
	tm_scale_to_secs(pbuf[3],&(sectiondata->startt));
	sectiondata->dir = atoi(pbuf[4]);
	
	sectiondata->tscale = atof(pbuf[5]);
	sectiondata->tunit = malloc(sizeof(char) * strlen(pbuf[6])+1);
	strcpy(sectiondata->tunit,pbuf[6]);
	if(i > 9)
	  {
	    sectiondata->outunit = malloc(sizeof(char) * strlen(pbuf[9])+1);
	    strcpy(sectiondata->outunit,pbuf[9]);
	  }else
	    sectiondata->outunit = NULL;
	
	sectiondata->time = trs->time;
	sectiondata->lastdump = trs->time;
	read_section_coordinates(pbuf[0],sectiondata,trs->model);
	emstag(LDEBUG,"Tracerstats:init","Section flux, translated coordiantes...");
	trstat_add_domain_function(sectiondata->name,trs->ndomainfn3d,2,section_create,section_gather,section_scatter);//,NULL);
	trs->ndomainfn3d++;
	emstag(LDEBUG,"Tracerstats:init","Section flux calculating for %s, every - %.0f sec, to file: %s, direction %d, starttime %.0f tscale %.1f tunit %s outscale %s ",tr_3d[n][1],
	       sectiondata->dt,sectiondata->fileout,
	       sectiondata->dir,sectiondata->startt,
	       sectiondata->tscale,sectiondata->tunit,
	       (i > 8?pbuf[8]:"none"));
      }
  }

  /* Initialise SED tracerstat arrays */
  if (trs->sednz) {
    for (n = 0; n < on_sed; n++) {
      m = om_sed[n];
      fill_stat3dfcn(trs, operation_sed[n], tr_sed[n],
		     trs->statfn_sed, trs->trname_sed,
		     m, &trs->w_sed, trs->stat_type_sed);
    }
  }
  
  emstag(LTRACE,"Tracerstats:init","Starting 2d  2.phase - %d",trs->ntrS);

  statfn2d_t* sfn2d;
  /*-----------------------------------------------------------------*/
  /* Set up the required 2d tracers                                  */
  if(trs->ntrS) {
    trs->trname_2d = malloc(trs->ntrS * sizeof(char *));
    for (n = 0; n < trs->ntrS; n++) {
      trs->trname_2d[n] = malloc(MAXSTRLEN * sizeof(char));
    }
    trs->stat_typeS = i_alloc_1d(trs->ntrS);
    trs->w1S = d_alloc_1d(trs->ntrS);
    trs->w2S = d_alloc_1d(trs->ntrS);
    trs->w3S = d_alloc_1d(trs->ntrS);
    memset(trs->stat_typeS, 0, trs->ntrS * sizeof(int));
    trs->map_2d = i_alloc_2d(nmS, trs->ntrS);
  }

  /* Initialise 2D tracerstat arrays                                 */
  for (n = 0; n < onS; n++) {
    m = om_2d[n];
    if (strcmp(operation_2d[n], "vmax") == 0) {
      trs->stat_typeS[m] = VMAX;
      strcpy(trs->trname_2d[m],tr_2d[n][0]);
      sfn2d = &trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn];
      sfn2d->type = VMAX;
      sfn2d->statfn = vertical_max;
      sfn2d->n = m;
      sfn2d->data = malloc(sizeof(int_dummy_t));
      if(sfn2d->data == NULL)
      {
        emstag(LERROR,"Tracerstats:init","Vertical maximum, could not assign memory, exiting!");
        exit(0);
      }
      emstag(LTRACE,"Tracerstats:init","vmax, created - assigning!");
      ((int_dummy_t*)sfn2d->data)->d = atoi(tr_2d[n][2]);

      strcpy(trs->trname_3d[((int_dummy_t*)sfn2d->data)->d],tr_2d[n][1]);
      sfn2d->step = extract_step(model, tr_2d[n][3]);
      trs->nprestep[sfn2d->step]++;
      emstag(LTRACE,"Tracerstats:init","vmax, tracer created on %s",trs->trname_3d[((int_dummy_t*)trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn].data)->d]);

      trs->statfn2d[STATFN0].nfn++;
    }else

    if (strcmp(operation_2d[n], "vmean") == 0) {
      trs->stat_typeS[m] = VMEAN;
      strcpy(trs->trname_2d[m],tr_2d[n][0]);
/*      strcpy(trs->trname_3d[atoi(tr_2d[n][2])],tr_2d[n][1]); */

      sfn2d = &trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn];
      sfn2d->type =VMEAN;
      sfn2d->statfn = vertical_mean;
      sfn2d->n = m;
      sfn2d->data = malloc(sizeof(int_dummy_t));
      if(sfn2d->data == NULL)
      {
        emstag(LERROR,"Tracerstats:init","Vertical integral, could not assign memory, exiting!");
        exit(0);
      }
      emstag(LTRACE,"Tracerstats:init","vmean, created - assigning!");
      ((int_dummy_t*)sfn2d->data)->d = atoi(tr_2d[n][2]);

      strcpy(trs->trname_3d[((int_dummy_t*)sfn2d->data)->d],tr_2d[n][1]);
      sfn2d->step = extract_step(model, tr_2d[n][3]);
      trs->nprestep[sfn2d->step]++;
      emstag(LTRACE,"Tracerstats:init","vmean, tracer created on %s",trs->trname_3d[((int_dummy_t*)trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn].data)->d]);

      trs->statfn2d[STATFN0].nfn++;
    }else

    if (strcmp(operation_2d[n], "vint") == 0) {
      trs->stat_typeS[m] = VINT;
      strcpy(trs->trname_2d[m],tr_2d[n][0]);
      sfn2d = &trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn];
      sfn2d->type = VINT;
      sfn2d->statfn = vertical_integrals;
      sfn2d->n = m;
      sfn2d->data = malloc(sizeof(int_dummy_t));
      if(sfn2d->data == NULL)
      {
        emstag(LERROR,"Tracerstats:init","Vertical integral, could not assign memory, exiting!");
        exit(0);
      }
      emstag(LTRACE,"Tracerstats:init","vint, created - assigning!");
      ((int_dummy_t*)trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn].data)->d = atoi(tr_2d[n][2]);

      strcpy(trs->trname_3d[((int_dummy_t*)sfn2d->data)->d],tr_2d[n][1]);
      sfn2d->step = extract_step(model, tr_2d[n][3]);
      trs->nprestep[sfn2d->step]++;
      emstag(LTRACE,"Tracerstats:init","vint, tracer created on %s",trs->trname_3d[((int_dummy_t*)sfn2d->data)->d]);

      trs->statfn2d[STATFN0].nfn++;
    }
    else

    if (strcmp(operation_2d[n], "copy_layer") == 0) {
      trs->stat_typeS[m] = COPYLAYER;
      strcpy(trs->trname_2d[m],tr_2d[n][0]);
      if (strcmp(tr_2d[n][4], "surface") == 0) 
	trs->w1S[m] = -1.0;
      else if (strcmp(tr_2d[n][4], "bottom") == 0) 
	trs->w1S[m] = -2.0;
      else
	trs->w1S[m] = atof(tr_2d[n][4]);
      sfn2d = &trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn];
      sfn2d->type = COPYLAYER;
      sfn2d->statfn = copy_layer;
      sfn2d->n = m;
      sfn2d->data = malloc(sizeof(int_dummy_t));
      if(sfn2d->data == NULL)
      {
        emstag(LERROR,"Tracerstats:init","Layer copy, could not assign memory, exiting!");
        exit(0);
      }
      emstag(LTRACE,"Tracerstats:init","copy_layer, created - assigning!");
      ((int_dummy_t*)trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn].data)->d = atoi(tr_2d[n][2]);

      strcpy(trs->trname_3d[((int_dummy_t*)sfn2d->data)->d],tr_2d[n][1]);
      sfn2d->step = extract_step(model, tr_2d[n][3]);
      trs->nprestep[sfn2d->step]++;
      emstag(LTRACE,"Tracerstats:init","vint, tracer created on %s",trs->trname_3d[((int_dummy_t*)sfn2d->data)->d]);

      trs->statfn2d[STATFN0].nfn++;
    }
    else

    if (strcmp(operation_2d[n], "mean") == 0) {
      tm_scale_to_secs(tr_2d[n][2], &trs->w1S[m]);
      trs->w2S[m] = 0.0;
      trs->stat_typeS[m] = MEAN;
      strcpy(trs->trname_2d[m],tr_2d[n][0]);
      strcpy(trs->trname_2d[m + 1],tr_2d[n][1]);
      sfn2d = &trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn];
      sfn2d->type = MEAN;
      sfn2d->statfn = tracer_mean_2d;
      sfn2d->n = m;
      sfn2d->data = NULL;
      sfn2d->step = extract_step(model, tr_2d[n][3]);
      trs->nprestep[sfn2d->step]++;

      trs->statfn2d[STATFN0].nfn++;
    }
    else

    if (strcmp(operation_2d[n], "run_mean") == 0) {
      tm_scale_to_secs(tr_2d[n][2], &trs->w1S[m]);
      trs->w2S[m] = 0.0;
      trs->stat_typeS[m] = RUN_MEAN;
      strcpy(trs->trname_2d[m],tr_2d[n][0]);
      strcpy(trs->trname_2d[m + 1],tr_2d[n][1]);
      sfn2d = &trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn];
      sfn2d->type = RUN_MEAN;
      sfn2d->statfn = tracer_run_mean_2d;
      sfn2d->n = m;
      sfn2d->data = NULL;
      sfn2d->step = extract_step(model, tr_2d[n][3]);
      trs->nprestep[sfn2d->step]++;

      trs->statfn2d[STATFN0].nfn++;
    }
    else

    if (strcmp(operation_2d[n], "max") == 0) {
      tm_scale_to_secs(tr_2d[n][2], &trs->w1S[m]);
      trs->w2S[m] = 0.0;
      trs->stat_typeS[m] = MAX;
      strcpy(trs->trname_2d[m],tr_2d[n][0]);
      strcpy(trs->trname_2d[m + 1],tr_2d[n][1]);
      sfn2d = &trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn];
      sfn2d->type = MAX;
      sfn2d->statfn = tracer_max_2d;
      sfn2d->n = m;
      sfn2d->data = NULL;
      sfn2d->step = extract_step(model, tr_2d[n][3]);
      trs->nprestep[sfn2d->step]++;

      trs->statfn2d[STATFN0].nfn++;
    }
    else

    if (strcmp(operation_2d[n], "min") == 0) {
      tm_scale_to_secs(tr_2d[n][2], &trs->w1S[m]);
      trs->w2S[m] = 0.0;
      trs->stat_typeS[m] = MIN;
      strcpy(trs->trname_2d[m],tr_2d[n][0]);
      strcpy(trs->trname_2d[m + 1],tr_2d[n][1]);
      sfn2d = &trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn];
      sfn2d->type = MIN;
      sfn2d->statfn = tracer_min_2d;
      sfn2d->n = m;
      sfn2d->data = NULL;
      sfn2d->step = extract_step(model, tr_2d[n][3]);
      trs->nprestep[sfn2d->step]++;

      trs->statfn2d[STATFN0].nfn++;
    }
    else

    if (strcmp(operation_2d[n], "stdev") == 0) {
      trs->stat_typeS[m] = STDEV;
      strcpy(trs->trname_2d[m],tr_2d[n][0]);

      trs->statfn2d[STATFN3].statfn[trs->statfn2d[STATFN3].nfn].type = STDEV;
      trs->statfn2d[STATFN3].statfn[trs->statfn2d[STATFN3].nfn].statfn = tracer_stdev_2d;
      trs->statfn2d[STATFN3].statfn[trs->statfn2d[STATFN3].nfn].n = m;
      trs->statfn2d[STATFN3].statfn[trs->statfn2d[STATFN3].nfn].data = NULL;
      trs->statfn2d[STATFN3].statfn[trs->statfn2d[STATFN3].nfn].step = extract_step(model, tr_2d[n][3]);
      trs->nprestep[trs->statfn2d[STATFN3].statfn[trs->statfn2d[STATFN3].nfn].step]++;
      trs->statfn2d[STATFN3].nfn++;
    }
    else

    if (strcmp(operation_2d[n], "variance") == 0) {
      trs->stat_typeS[m] = VAR;
      strcpy(trs->trname_2d[m],tr_2d[n][0]);
      trs->statfn2d[STATFN2].statfn[trs->statfn2d[STATFN2].nfn].type = VAR;
      trs->statfn2d[STATFN2].statfn[trs->statfn2d[STATFN2].nfn].statfn = tracer_var_2d;
      trs->statfn2d[STATFN2].statfn[trs->statfn2d[STATFN2].nfn].n = m;
      trs->statfn2d[STATFN2].statfn[trs->statfn2d[STATFN2].nfn].data = NULL;
      trs->statfn2d[STATFN2].statfn[trs->statfn2d[STATFN2].nfn].step = extract_step(model, tr_2d[n][3]);
      trs->nprestep[trs->statfn2d[STATFN2].statfn[trs->statfn2d[STATFN2].nfn].step]++;

      trs->statfn2d[STATFN2].nfn++;
    }
    else

    if (strcmp(operation_2d[n], "cov") == 0) {
      trs->stat_typeS[m] = COV;
      strcpy(trs->trname_2d[m],tr_2d[n][0]);
      trs->statfn2d[STATFN5].statfn[trs->statfn2d[STATFN5].nfn].type = COV;
      trs->statfn2d[STATFN5].statfn[trs->statfn2d[STATFN5].nfn].statfn = tracer_cov_2d;
      trs->statfn2d[STATFN5].statfn[trs->statfn2d[STATFN5].nfn].n = m;
      trs->statfn2d[STATFN5].statfn[trs->statfn2d[STATFN5].nfn].data = NULL;
      trs->statfn2d[STATFN5].statfn[trs->statfn2d[STATFN5].nfn].step = extract_step(model, tr_2d[n][3]);
      trs->nprestep[trs->statfn2d[STATFN5].statfn[trs->statfn2d[STATFN5].nfn].step]++;
      trs->statfn2d[STATFN5].nfn++;
    }
    else

    if (strcmp(operation_2d[n], "corr") == 0) {
      trs->stat_typeS[m] = CORR;
      strcpy(trs->trname_2d[m],tr_2d[n][0]);
      trs->statfn2d[STATFN1].statfn[trs->statfn2d[STATFN1].nfn].type = CORR;
      trs->statfn2d[STATFN1].statfn[trs->statfn2d[STATFN1].nfn].statfn = corr_2d_pre;
      trs->statfn2d[STATFN1].statfn[trs->statfn2d[STATFN1].nfn].n = m;
      trs->statfn2d[STATFN1].statfn[trs->statfn2d[STATFN1].nfn].data = NULL;
      trs->statfn2d[STATFN1].statfn[trs->statfn2d[STATFN1].nfn].step = extract_step(model, tr_2d[n][3]);

      trs->statfn2d[STATFN4].statfn[trs->statfn2d[STATFN4].nfn].type = CORR;
      trs->statfn2d[STATFN4].statfn[trs->statfn2d[STATFN4].nfn].statfn = corr_2d_post;
      trs->statfn2d[STATFN4].statfn[trs->statfn2d[STATFN4].nfn].n = m;
      trs->statfn2d[STATFN4].statfn[trs->statfn2d[STATFN4].nfn].data = NULL;

      trs->statfn2d[STATFN4].statfn[trs->statfn2d[STATFN4].nfn].step = extract_step(model, tr_2d[n][3]);
      trs->nprestep[trs->statfn2d[STATFN4].statfn[trs->statfn2d[STATFN4].nfn].step]++;

      trs->statfn2d[STATFN1].nfn++;
      trs->statfn2d[STATFN4].nfn++;
    }
    else

    if (strcmp(operation_2d[n], "diff") == 0) {
      trs->stat_typeS[m] = DIFF;
      strcpy(trs->trname_2d[m],tr_2d[n][0]);
      strcpy(trs->trname_2d[m + 1],tr_2d[n][1]);
      strcpy(trs->trname_2d[m + 2],tr_2d[n][2]);
      trs->statfn2d[STATFN6].statfn[trs->statfn2d[STATFN6].nfn].type = DIFF;
      trs->statfn2d[STATFN6].statfn[trs->statfn2d[STATFN6].nfn].statfn = tracer_diff_2d;
      trs->statfn2d[STATFN6].statfn[trs->statfn2d[STATFN6].nfn].n = m;
      trs->statfn2d[STATFN6].statfn[trs->statfn2d[STATFN6].nfn].data = NULL;

      trs->statfn2d[STATFN6].statfn[trs->statfn2d[STATFN6].nfn].step = extract_step(model, tr_2d[n][3]);
      trs->nprestep[trs->statfn2d[STATFN6].statfn[trs->statfn2d[STATFN6].nfn].step]++;

      trs->statfn2d[STATFN6].nfn++;
    }
    else

    if (strcmp(operation_2d[n], "rmse") == 0) {
      trs->stat_typeS[m] = RMSE;
      strcpy(trs->trname_2d[m],tr_2d[n][0]);
      strcpy(trs->trname_2d[m + 1],tr_2d[n][1]);
      strcpy(trs->trname_2d[m + 2],tr_2d[n][2]);
      trs->statfn2d[STATFN6].statfn[trs->statfn2d[STATFN6].nfn].type = RMSE;
      trs->statfn2d[STATFN6].statfn[trs->statfn2d[STATFN6].nfn].statfn = tracer_rmse_2d;
      trs->statfn2d[STATFN6].statfn[trs->statfn2d[STATFN6].nfn].n = m;
      trs->statfn2d[STATFN6].statfn[trs->statfn2d[STATFN6].nfn].data = NULL;

      trs->statfn2d[STATFN6].statfn[trs->statfn2d[STATFN6].nfn].step = extract_step(model, tr_2d[n][3]);
      trs->nprestep[trs->statfn2d[STATFN6].statfn[trs->statfn2d[STATFN6].nfn].step]++;

      trs->statfn2d[STATFN6].nfn++;
    }
    else

    if (strcmp(operation_2d[n], "sum") == 0) {
    	int* nsumtrs =malloc(sizeof(int));
      trs->stat_typeS[m] = SUM;
      strcpy(trs->trname_2d[m],tr_2d[n][0]);
      /*
      strcpy(trs->trname_2d[m + 1],tr_2d[n][1]);
      strcpy(trs->trname_2d[m + 2],tr_2d[n][2]);
      */
      trs->statfn2d[STATFN6].statfn[trs->statfn2d[STATFN6].nfn].type = SUM;
      trs->statfn2d[STATFN6].statfn[trs->statfn2d[STATFN6].nfn].statfn = tracer_sum_2d;
      trs->statfn2d[STATFN6].statfn[trs->statfn2d[STATFN6].nfn].n = m;
      trs->statfn2d[STATFN6].statfn[trs->statfn2d[STATFN6].nfn].data = nsumtrs;

      trs->statfn2d[STATFN6].statfn[trs->statfn2d[STATFN6].nfn].step = extract_step(model, tr_2d[n][3]);
      trs->nprestep[trs->statfn2d[STATFN6].statfn[trs->statfn2d[STATFN6].nfn].step]++;
			
      init_sum(tr_2d[n][1],nsumtrs,trs->trname_2d, m+1);
      trs->statfn2d[STATFN6].nfn++;
      //
    }
    else

    if (strcmp(operation_2d[n], "vdiff") == 0) {
      double tt;
      trs->stat_type[m] = TRVDIFF;
      strcpy(trs->trname_2d[m],tr_2d[n][0]);
 /*     strcpy(trs->trname_2d[m + 1],tr_2d[n][1]);  */
      /* get the second and third parameter as well as the index for the 3d tracer back */

      tint = parseline(tr_2d[n][2],pbuf,9);
      trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn].type = TRVDIFF;
      trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn].n = m;
      trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn].data = malloc(sizeof(vdiff_range_t));
      trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn].step = extract_step(model, tr_2d[n][3]);
      trs->nprestep[trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn].step]++;

      vdiff = (vdiff_range_t*) trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn].data;
      if(vdiff == NULL)
      {
        emstag(LERROR,"Tracerstats:init","Vertical diff, could not assign memory, exiting!");
        exit(0);
      }
      emstag(LDEBUG,"Tracerstats:init","vdiff, created - assigning!");
      vdiff->strict = 0;
      vdiff->pairname = NULL;
      if(tint > 3)
        vdiff->strict =atoi(pbuf[3]);

      if(tint > 4)
      {
        vdiff->mode =atoi(pbuf[4]);
      }
      else
        vdiff->mode= VDIFFMODNORMAL;

      vdiff->ntr3d = atoi(pbuf[0]);

      vdiff->botrange = split_ints(pbuf[1], &vdiff->nb);
      vdiff->toprange = split_ints(pbuf[2], &vdiff->nt);

      switch(vdiff->mode)
      {
        case VDIFFMODVOL:
          trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn].statfn = tracer_vertical_strat_vol;
				  strcpy(trs->trname_3d[vdiff->ntr3d],tr_2d[n][1]);
          break;
        case VDIFFMODDIV:
          trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn].statfn = tracer_vertical_strat_div;
  				strcpy(trs->trname_3d[vdiff->ntr3d],tr_2d[n][1]);
          break;
        case VDIFFMODVOLDIV:
          trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn].statfn = tracer_vertical_strat_voldiv;
					strcpy(trs->trname_3d[vdiff->ntr3d],tr_2d[n][1]);
          break;
        case VDIFFSEDFLUX1:
          trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn].statfn = tracer_vertical_sedflux1;
          emstag(LDEBUG,"Tracerstats:init","vdiff, tracer %s created on %s at step %s ",trs->trname_2d[m],trs->trname_3d[vdiff->ntr3d],tr_2d[n][3]);
          strcpy(trs->trname_sed[vdiff->ntr3d],tr_2d[n][1]);

          break;
        case VDIFFSEDFLUX2:
        case VDIFFSEDFLUX3:

          if(tint < 6)
          {
            emstag(LFATAL,"trcerstats:tracerstats:init","vdiff sedflux2, too few arguments, require sedflux1 name for sedflux2 %s!",trs->trname_2d[m]);
            exit(0);
          }else
          {
            vdiff->pairname = malloc(sizeof(char) * strlen(pbuf[5])+1);
            strcpy(vdiff->pairname,pbuf[5]);
          }
          trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn].statfn = tracer_vertical_sedflux2;
          free(vdiff->botrange);
          free(vdiff->toprange);

          strcpy(trs->trname_sed[vdiff->ntr3d],tr_2d[n][1]);
          vdiff->botrange = i_alloc_1d(3);
          if(tint > 6)
          {
            tm_scale_to_secs(pbuf[6],&tt);
            (vdiff->botrange[0]) = (int)tt;
          /*  vdiff->botrange[0] = atoi(pbuf[6]); */
          }
          else
            vdiff->botrange[0] = atoi(pbuf[1]);

          vdiff->botrange[1] = trs->time;
          vdiff->botrange[2] = 0;
          vdiff->nb = -1;
          emstag(LDEBUG,"Tracerstats:init","vdiff, tracer %s created on %s at step %s to pair %s",trs->trname_2d[m],trs->trname_3d[vdiff->ntr3d],tr_2d[n][3],vdiff->pairname);

          break;
        default:
          trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn].statfn = tracer_vertical_strat;
      }
      trs->statfn2d[STATFN0].nfn++;
    }else

    if (strcmp(operation_2d[n], "sedflux") == 0) {
      trs->stat_type[m] = TRVDIFF;
      strcpy(trs->trname_2d[m],tr_2d[n][0]);
 /*     strcpy(trs->trname_2d[m + 1],tr_2d[n][1]);  */
      /* get the second and third parameter as well as the index for the 3d tracer back */

      tint = parseline(tr_2d[n][2],pbuf,9);
      trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn].type = TRVDIFF;
      trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn].n = m;
      trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn].data = malloc(sizeof(vdiff_range_t));
      trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn].step = extract_step(model, tr_2d[n][3]);
      trs->nprestep[trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn].step]++;

      vdiff = (vdiff_range_t*) trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn].data;
      if(vdiff == NULL)
      {
        emstag(LERROR,"Tracerstats:init","Vertical diff, could not assign memory, exiting!");
        exit(0);
      }
      emstag(LDEBUG,"Tracerstats:init","vdiff, created - assigning!");
      vdiff->strict = 0;
      vdiff->pairname = NULL;
      strcpy(trs->trname_sed[vdiff->ntr3d],tr_2d[n][1]);

      if(tint > 3)
        vdiff->strict =atoi(pbuf[3]);

      if(tint > 4)
      {
        vdiff->mode =atoi(pbuf[4]);
      }
      else
        vdiff->mode= VDIFFMODNORMAL;

      vdiff->ntr3d = atoi(pbuf[0]);

      vdiff->botrange = split_ints(pbuf[1], &vdiff->nb);
      vdiff->toprange = split_ints(pbuf[2], &vdiff->nt);

      switch(vdiff->mode)
      {
        case VDIFFSEDFLUX1:
          trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn].statfn = tracer_vertical_sedflux1;
          emstag(LDEBUG,"Tracerstats:init","vdiff, tracer %s created on %s at step %s ",trs->trname_2d[m],trs->trname_3d[vdiff->ntr3d],tr_2d[n][3]);


          break;
        case VDIFFSEDFLUX2:
        case VDIFFSEDFLUX3:
          if(tint < 6)
          {
            emstag(LFATAL,"trcerstats:tracerstats:init","vdiff sedflux2, too few arguments, require sedflux1 name for sedflux2 %s!",trs->trname_2d[m]);
            exit(0);
          }else
          {
            vdiff->pairname = malloc(sizeof(char) * strlen(pbuf[5])+1);
            strcpy(vdiff->pairname,pbuf[5]);
          }
          trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn].statfn = tracer_vertical_sedflux2;
          free(vdiff->botrange);
          free(vdiff->toprange);

          strcpy(trs->trname_sed[vdiff->ntr3d],tr_2d[n][1]);
          vdiff->botrange = i_alloc_1d(3);
          if(tint > 6)
            vdiff->botrange[0] = atoi(pbuf[6]);
          else
            vdiff->botrange[0] = atoi(pbuf[1]);

          vdiff->botrange[1] = trs->time;
          vdiff->botrange[2] = 0;
          vdiff->nb = -1;
          emstag(LDEBUG,"Tracerstats:init","vdiff, tracer %s created on %s at step %s to pair %s",trs->trname_2d[m],trs->trname_3d[vdiff->ntr3d],tr_2d[n][3],vdiff->pairname);

          break;
        default:
          trs->statfn2d[STATFN0].statfn[trs->statfn2d[STATFN0].nfn].statfn = tracer_vertical_strat;
      }
      trs->statfn2d[STATFN0].nfn++;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set up any maps required                                        */
  /* 3d maps                                                         */
  for (n = 0; n < on; n++) {
    m = om_3d[n];
    nm = 0;
    /* Check that tracerstats that require means or standard         */
    /* deviations in their calculations have these operations        */
    /* defined and make maps to the relevant tracers.                */
    /* m this tracer, nm the number of following tracers to consider */
    switch (trs->stat_type[m]) {
    case STDEV:
      strcpy(buf, "standard deviation");
      nm = 1;
      break;
    case VAR:
      strcpy(buf, "variance");
      nm = 1;
      break;
    case COV:
      strcpy(buf, "covariance");
      nm = 2;
      break;
    case CORR:
      strcpy(buf, "correlation coefficient");
      nm = 2;
      break;
    }

    if (nm) {
      for(i = 0; i < nm; i++) {
        trs->map_3d[m][i] = NOT;
        for(str = 0; str < on; str++) {
          if(strcmp(operation_3d[str], "mean") == 0 &&
             strcmp(tr_3d[str][1], tr_3d[n][i + 1]) == 0)
            trs->map_3d[m][i] = om_3d[str];
        }
        if(trs->map_3d[m][i] == NOT) {
          emstag(LFATAL,"tracerstats:init","statistics : %s of tracer %s requires this tracer's mean\n", buf, tr_3d[n][i + 1]);
          exit(0);
        }
        if (strcmp(buf, "correlation coefficient") == 0) {
          trs->map_3d[m][i + 2] = NOT;
          for(str = 0; str < on; str++) {
            if(strcmp(operation_3d[str], "stdev") == 0 &&
               strcmp(tr_3d[str][1], tr_3d[n][i + 1]) == 0)
              trs->map_3d[m][i + 2] = om_3d[str];
          }
          if(trs->map_3d[m][i + 2] == NOT) {
            emstag(LFATAL,"tracerstats:init","statistics : %s of tracer %s requires this tracer's standard deviation\n", buf, tr_3d[n][i + 1]);
            exit(0);
          }
        }
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Sed maps                                                         */
  for (n = 0; n < on_sed; n++) {
    m = om_sed[n];
    nm_sed = 0;
    /* Check that tracerstats that require means or standard         */
    /* deviations in their calculations have these operations        */
    /* defined and make maps to the relevant tracers.                */
    /* m this tracer, nm the number of following tracers to consider */
    switch (trs->stat_type_sed[m]) {
    case STDEV:
      strcpy(buf, "standard deviation");
      nm = 1;
      break;
    case VAR:
      strcpy(buf, "variance");
      nm = 1;
      break;
    }

    if (nm_sed) {
      for(i = 0; i < nm_sed; i++) {
        trs->map_sed[m][i] = NOT;
        for(str = 0; str < on; str++) {
          if(strcmp(operation_sed[str], "mean") == 0 &&
             strcmp(tr_sed[str][1], tr_sed[n][i + 1]) == 0)
            trs->map_sed[m][i] = om_sed[str];
        }
        if(trs->map_sed[m][i] == NOT) {
          emstag(LFATAL,"tracerstats:init","statistics : %s of tracer %s requires this tracer's mean\n", buf, tr_sed[n][i + 1]);
          exit(0);
        }
      }
    }
  }
  
  /* 2d maps                                                         */
  for (n = 0; n < onS; n++) {
    m = om_2d[n];
    nmS = 0;
    /* Check that tracerstats that require means or standard         */
    /* deviations in their calculations have these operations        */
    /* defined and make maps to the relevant tracers.                */
    switch (trs->stat_typeS[m]) {
    case STDEV:
      strcpy(buf, "standard deviation");
      nmS = 1;
      break;
    case VAR:
      strcpy(buf, "variance");
      nmS = 1;
      break;
    case COV:
      strcpy(buf, "covariance");
      nmS = 2;
      break;
    case CORR:
      strcpy(buf, "correlation coefficient");
      nmS = 2;
      break;
    }
    if (nmS) {
      for(i = 0; i < nmS; i++) {
        trs->map_2d[m][i] = NOT;
        for(str = 0; str < onS; str++) {
          if(strcmp(operation_2d[str], "mean") == 0 &&
             strcmp(tr_2d[str][1], tr_2d[n][i + 1]) == 0)
            trs->map_2d[m][i] = om_2d[str];
        }
        if(trs->map_2d[m][i] == NOT) {
          printf("error : statistics : %s of tracer %s requires this tracer's mean\n", buf, tr_2d[n][i + 1]);
          exit(0);
        }
        if (strcmp(buf, "correlation coefficient") == 0) {
          trs->map_2d[m][i + 2] = NOT;
          for(str = 0; str < onS; str++) {
            if(strcmp(operation_2d[str], "stdev") == 0 &&
               strcmp(tr_2d[str][1], tr_2d[n][i + 1]) == 0)
              trs->map_2d[m][i + 2] = om_2d[str];
          }
          if(trs->map_2d[m][i + 2] == NOT) {
            emstag(LFATAL,"tracerstats:init","statistics : %s of tracer %s requires this tracer's standard deviation\n", buf, tr_2d[n][i + 1]);
            exit(0);
          }
        }
      }
    }
    /* Check that tracerstats that require vertical integrals have   */
    /* the corresponding 3D tracer in the tracer list.               */
    /* don't need this anymore
    if (trs->stat_typeS[m] & (VINT)) {
      trs->map_2d[m][0] = NOT;
      for(str = 0; str < trs->ntr; str++) {
        if(strcmp(tr_2d[n][1], trs->trname_3d[str]) == 0)
          trs->map_2d[m][0] = str;
      }
      if(trs->map_2d[m][0] == NOT) {
        emstag(LFATAL,"tracerstats:init","statistics : vertical integral of tracer %s requires tracer %s\n", tr_2d[n][1], tr_2d[n][1]);
        exit(0);
      }
    }
    */
  }
  /*-----------------------------------------------------------------*/
  /* Allocate 3d tracer memory                                          */
  if (trs->ntr) {
    trs->tr_wc = (double ***)p_alloc_2d(trs->nz, trs->ntr);
    trs->tr_e1wc = (double **)d_alloc_2d(trs->nz, trs->ntr);
    trs->tr_e2wc = (double **)d_alloc_2d(trs->nz, trs->ntr);

    trs->tmap_3d = i_get_tmap_3d(model, trs->ntr, trs->trname_3d);

    /* check if all tracers have been found */
    for(i = 0;i< trs->ntr;i++) {
      if(!i_is_valid_tracer(trs->tmap_3d[i]) ) {	
	emstag(LPANIC,"tracerstats:init","Attempting to reference a water tracer which doesn't exist in host model [%s] exiting...",trs->trname_3d[i]);
	exit(1);
      }
    }
  }

  /* verify any sediment and special 2d functions */
  if(trs->nsed && trs->sednz) {
    trs->tr_sed   = (double ***)p_alloc_2d(trs->sednz, trs->nsed);
    trs->tmap_sed = i_get_tmap_sed(model, trs->nsed, trs->trname_sed);
    
    /* check if all tracers have been found */
    for(i = 0;i< trs->nsed;i++) {
      if(!i_is_valid_sed_tracer(trs->tmap_sed[i]) )
	{	
	  emstag(LPANIC,"tracerstats:init","Attempting to reference a sediment tracer which doesn't exist in host model [%s] at %d exiting...",trs->trname_sed[i], i);
	  exit(1);
	}
    }
  }

  /* check that there is a corresponding sedflux stage 1 present
   * we know that the vdiff 2d are at STATFN0 */
	if(nstatorder2d[STATFN0])
	{
	  statfn2d_t* sfn2dl = trs->statfn2d[STATFN0].statfn;
	  vdiff_range_t* vdiffr;
	  for(n = 0; n < trs->statfn2d[STATFN0].nfn; n++)
	  {
	    if(sfn2dl[n].type == TRVDIFF )
	    {
	      vdiffr = (vdiff_range_t*) sfn2dl[n].data;
	      if(vdiffr->mode == VDIFFSEDFLUX2 || vdiffr->mode == VDIFFSEDFLUX3)
	      {
	        vdiffr->nb =-1;
	        for(i = 0; i < trs->statfn2d[STATFN0].nfn; i++)
	        {
	          if(sfn2dl[n].type == TRVDIFF &&
	               ((vdiff_range_t*)sfn2dl[i].data)->mode == VDIFFSEDFLUX1)
	          {
	            if(strcmp(vdiffr->pairname,trs->trname_2d[sfn2dl[i].n]) == 0)
	            {
	              vdiffr->nb = sfn2dl[i].n;
	              if(trs->tmap_sed[vdiffr->ntr3d] != trs->tmap_sed[((vdiff_range_t*)sfn2dl[i].data)->ntr3d])
	              {
	                emstag(LFATAL,"tracerstats:init"," found sedflux pair, but looking atdifferent tracers! %s - %s",trs->trname_2d[sfn2dl[n].n],vdiffr->pairname);
	                exit(0);
	              }
	              break;
	            }
	          }
	        }
	        if(vdiffr->nb == -1)
	        {
	          emstag(LFATAL,"tracerstats:init","Missing stage 1 vediff sedflux for %s, value %s",trs->trname_2d[sfn2dl[n].n],vdiffr->pairname);
	          exit(0);
	        }
	      }
	    }
	  }
	}

  if(trs->ntrS) {
    trs->tr_in = (double**)p_alloc_1d(trs->ntrS);
    trs->tmap_2d = i_get_tmap_2d(model, trs->ntrS, trs->trname_2d);
    /* check if all tracers have been found */
    for(i = 0;i< trs->ntrS;i++) {
      if(!i_is_valid_tracer(trs->tmap_2d[i]) ) {
	emstag(LPANIC,"tracerstats:init","Attempting to reference a 2d tracer which doesn't exist in host model [%s] exiting...",trs->trname_2d[i]);
	exit(1);
      }
    }
  }

  // allocate a few more vectors
  if(trs->use_u1)
    trs->u1 = d_alloc_1d(trs->nz);
  if(trs->use_u2)
    trs->u2 = d_alloc_1d(trs->nz);
  if(trs->use_area_e1_e2) {
    trs->area_e1 = d_alloc_1d(trs->nz);
    trs->area_e2 = d_alloc_1d(trs->nz);
  }

  /* free unused memory */
  if (nstatorder3d)
    i_free_1d(nstatorder3d);
  if (nstatorder2d)
    i_free_1d(nstatorder2d);
  if (nstatorder_sed)
    i_free_1d(nstatorder_sed);

  if(trs->nsed && trs->sednz) {
    free((char **)trname_sed);
    free((char **)operation_sed);
    if (tr_sed)c_free_3d(tr_sed);
    if (om_sed)i_free_1d(om_sed);
  }
  free((char **)trname_3d);
  free((char **)trname_2d);
  free((char **)operation_3d);
  free((char **)operation_2d);
  if (tr_2d)c_free_3d(tr_2d);
  if (tr_3d)c_free_3d(tr_3d);
  if (om_2d)i_free_1d(om_2d);
  if (om_3d)i_free_1d(om_3d);
}

/* END trs_init()                                                    */
/*-------------------------------------------------------------------*/




/*-------------------------------------------------------------------*/
/* Does the interface work. In this case a vertical integral of the  */
/* tracers 'temp' and 'salt' are calculated and written to the 2D    */
/* arrays 'temp_int' and 'salt_int' (if present).                    */
/*-------------------------------------------------------------------*/
void trs_step(void* model, trs_t *trs, int c, int step) {
  int n, i;
  statfn3d_t* sfn3d;
  statfn3d_t* sfn_sed;
  statfn2d_t* sfn2d;
  fluxfn_t ffn;
  domainfn3d_t dfn3d;

  /* 
   * Return if nothing to do for this step
   */
  if(!trs->do_tracerstats || trs->nprestep[step] == 0)
    return;

  /*-----------------------------------------------------------------*/
  /* Get the model attributes                                        */
  /* Get the interface column index                                  */
  i = i_get_interface_counter(model, c);

  /* Get the top index of the water column                           */
  trs->topk_wc = i_get_topk_wc(model, c);

  /* Get the bottom index of the water column                        */
  trs->botk_wc = i_get_botk_wc(model, c);

  /* And for the sediments */
  // Be careful how you use topk_sed as its minus one due to the stat funcs
  trs->topk_sed = i_get_topk_sed(model, c) - 1;
  trs->botk_sed = i_get_botk_sed(model, c);

  /* Get the model step  - figure out if this is the first time
   * this is called for this time step
   * as well as when the model starts  i.e. trs->news[step] = -1  */
  trs->o_nstep[step] = trs->e_nstep[step];
  trs->e_nstep[step] = i_get_nstep(model);
  if (trs->o_nstep[step] != trs->e_nstep[step] || trs->news[step] == -1) {
    for(i = 0; i < trs->ntr+trs->ntrS+trs->nsed;i++)
      trs->nstep[i]++;
    trs->news[step] = 1;
  }else
    trs->news[step] = 0;
  trs->new_step = trs->news[step];

  /* Get the timestep for this step                                  */
  trs->dt = i_get_model_timestep(model);

  /* Get the model time at this step                                  */
  trs->time = i_get_model_time(model);

  /* Get the water column layer thickness                            */
  i_get_dz_wc(model, c, trs->dz_wc);

  /* Get the water depth                                             */
  trs->depth_wc = i_get_depth(model, c);

  /* Get the cell areas for each face */
  trs->area_w = i_get_cellarea_w(model, c);
  if (trs->use_area_e1_e2) {
    i_get_cellarea_e1(model, c, trs->area_e1);
    i_get_cellarea_e2(model, c, trs->area_e2);
  }

  /* Get the surface elevation                                       */
  if(trs->use_eta)
    trs->eta = i_get_eta(model, c);

  /* Get the vertical diffusivity                                    */
  if(trs->use_Kz)
    i_get_Kz_wc(model, c, trs->Kz);

  /* Get the vertical diffusivity                                    */
  if(trs->use_Vz)
    i_get_Vz_wc(model, c, trs->Vz);

  /* Get the e1 flux                                                 */
  if(trs->use_fluxe1)
    i_get_fluxe1_wc(model, c, trs->u1flux3d);

  /* Get the e2 flux                                                 */
  if(trs->use_fluxe2)
    i_get_fluxe2_wc(model, c, trs->u2flux3d);

  /* Get u1 */
  if(trs->use_u1)
    i_get_u1_wc(model, c, trs->u1);
  
  /* Get u2 */
  if(trs->use_u2)
    i_get_u2_wc(model, c, trs->u2);

  /* Get w */
  if(trs->use_w)
    i_get_w_wc(model, c, trs->w);

  /* Read in the water column tracers                                */
  i_get_tracer_wc(model, c, trs->ntr, trs->tmap_3d, trs->tr_wc);
  i_get_tracer_e1(model, c, trs->ntr, trs->tmap_3d, trs->tr_e1wc);
  i_get_tracer_e2(model, c, trs->ntr, trs->tmap_3d, trs->tr_e2wc);

  /* Read in the sediment column tracers                                */
  if(trs->nsed && trs->sednz)
  {
    i_get_tracer_sed(model, c, trs->nsed, trs->tmap_sed, trs->tr_sed);
    i_get_dz_sed(model, c, trs->dz_sed);
  }

  /* Read in the 2D tracers                                          */
  i_get_tracer_2d(model, c, trs->ntrS, trs->tmap_2d, trs->tr_in);

  /*-----------------------------------------------------------------*/
  /* Do the work ...                                                 */
  /*
   * Iterate through the 3d WC functions
   */
  sfn3d = trs->statfn3d[STATFN0].statfn;
  for(n = 0; n < trs->statfn3d[STATFN0].nfn; n++) {
    if(sfn3d[n].step == step) {
      if (sfn3d[n].statfn != NULL) {
	// Old syntax
	sfn3d[n].statfn(trs,&(sfn3d[n].n),sfn3d[n].data,c);
      } else {
	// new
	sfn3d[n].statfn_new(trs, &trs->w_wc, trs->tr_wc, sfn3d[n].n,
			    trs->botk_wc, trs->topk_wc, trs->map_3d,
				    sfn3d[n].data);
      }
    }
  }

  for(n = 0; n < trs->nfluxfn; n++)
    {
      ffn =  trs->fluxfn[n];
      if(ffn.step == step)
	ffn.fluxfn(trs,ffn.flux,ffn.data,&(ffn.n),c);
    }
  
  for(n = 0; n < trs->ndomainfn3d; n++)
    {
      dfn3d =  trs->domainfn3d[n];
      if(dfn3d.step == step)
	dfn3d.domainfn(model,trs,&(dfn3d.n),dfn3d.data,c);
    }
  
  sfn3d = trs->statfn3d[STATFN1].statfn;
  for(n = 0; n < trs->statfn3d[STATFN1].nfn; n++) {
    if(sfn3d[n].step == step) {
      if (sfn3d[n].statfn != NULL) {
	// Old syntax
	sfn3d[n].statfn(trs,&(sfn3d[n].n),sfn3d[n].data,c);
      } else {
	// new
	sfn3d[n].statfn_new(trs, &trs->w_wc, trs->tr_wc, sfn3d[n].n,
			    trs->botk_wc, trs->topk_wc, trs->map_3d, 
			    sfn3d[n].data);
      }
    }
  }
  
  sfn3d = trs->statfn3d[STATFN2].statfn;
  for(n = 0; n < trs->statfn3d[STATFN2].nfn; n++) {
    if(sfn3d[n].step == step) {
      if (sfn3d[n].statfn != NULL) {
	// Old syntax
	sfn3d[n].statfn(trs,&(sfn3d[n].n),sfn3d[n].data,c);
      } else {
	// new
	sfn3d[n].statfn_new(trs, &trs->w_wc, trs->tr_wc, sfn3d[n].n,
			    trs->botk_wc, trs->topk_wc, trs->map_3d, 
			    sfn3d[n].data);
      }
    }
  }
  
  sfn3d = trs->statfn3d[STATFN3].statfn;
  for(n = 0; n < trs->statfn3d[STATFN3].nfn; n++) {
    if(sfn3d[n].step == step) {
      if (sfn3d[n].statfn != NULL) {
	// Old syntax
	sfn3d[n].statfn(trs,&(sfn3d[n].n),sfn3d[n].data,c);
      } else {
	// new
	sfn3d[n].statfn_new(trs, &trs->w_wc, trs->tr_wc, sfn3d[n].n,
			    trs->botk_wc, trs->topk_wc, trs->map_3d, 
			    sfn3d[n].data);
      }
    }
  }
  
  sfn3d = trs->statfn3d[STATFN4].statfn;
  for(n = 0; n < trs->statfn3d[STATFN4].nfn; n++) {
    if(sfn3d[n].step == step) {
      if (sfn3d[n].statfn != NULL) {
	// Old syntax
	sfn3d[n].statfn(trs,&(sfn3d[n].n),sfn3d[n].data,c);
      } else {
	// new
	sfn3d[n].statfn_new(trs, &trs->w_wc, trs->tr_wc, sfn3d[n].n,
			    trs->botk_wc, trs->topk_wc, trs->map_3d, 
			    sfn3d[n].data);
      }
    }
  }
  
  sfn3d = trs->statfn3d[STATFN5].statfn;
  for(n = 0; n < trs->statfn3d[STATFN5].nfn; n++) {
    if(sfn3d[n].step == step) {
      if (sfn3d[n].statfn != NULL) {
	// Old syntax
	sfn3d[n].statfn(trs,&(sfn3d[n].n),sfn3d[n].data,c);
      } else {
	// new
	sfn3d[n].statfn_new(trs, &trs->w_wc, trs->tr_wc, sfn3d[n].n,
			    trs->botk_wc, trs->topk_wc, trs->map_3d, 
			    sfn3d[n].data);
      }
    }
  }
  
  sfn3d = trs->statfn3d[STATFN6].statfn;
  for(n = 0; n < trs->statfn3d[STATFN6].nfn; n++) {
    if(sfn3d[n].step == step) {
      if (sfn3d[n].statfn != NULL) {
	// Old syntax
	sfn3d[n].statfn(trs,&(sfn3d[n].n),sfn3d[n].data,c);
      } else {
	// new
	sfn3d[n].statfn_new(trs, &trs->w_wc, trs->tr_wc, sfn3d[n].n,
			    trs->botk_wc, trs->topk_wc, trs->map_3d, 
			    sfn3d[n].data);
      }
    }
  }
  
  /*
   * Iterate through the 3d sediment functions
   */
  sfn_sed = trs->statfn_sed[STATFN0].statfn;
  for(n = 0; n < trs->statfn_sed[STATFN0].nfn; n++) {
    if(sfn_sed[n].step == step) {
      sfn_sed[n].statfn_new(trs, &trs->w_sed, trs->tr_sed, sfn_sed[n].n,
			    trs->botk_sed, trs->topk_sed, trs->map_sed,
			    sfn_sed[n].data);
    }
  }
  
  sfn_sed = trs->statfn_sed[STATFN1].statfn;
  for(n = 0; n < trs->statfn_sed[STATFN1].nfn; n++) {
    if(sfn_sed[n].step == step) {
      sfn_sed[n].statfn_new(trs, &trs->w_sed, trs->tr_sed, sfn_sed[n].n,
			    trs->botk_sed, trs->topk_sed, trs->map_sed,
			    sfn_sed[n].data);
    }
  }
  
  sfn_sed = trs->statfn_sed[STATFN2].statfn;
  for(n = 0; n < trs->statfn_sed[STATFN2].nfn; n++) {
    if(sfn_sed[n].step == step) {
      sfn_sed[n].statfn_new(trs, &trs->w_sed, trs->tr_sed, sfn_sed[n].n,
			    trs->botk_sed, trs->topk_sed, trs->map_sed,
			    sfn_sed[n].data);
    }
  }
  
  sfn_sed = trs->statfn_sed[STATFN3].statfn;
  for(n = 0; n < trs->statfn_sed[STATFN3].nfn; n++) {
    if(sfn_sed[n].step == step) {
      sfn_sed[n].statfn_new(trs, &trs->w_sed, trs->tr_sed, sfn_sed[n].n,
			    trs->botk_sed, trs->topk_sed, trs->map_sed,
			    sfn_sed[n].data);
    }
  }
  
  sfn_sed = trs->statfn_sed[STATFN4].statfn;
  for(n = 0; n < trs->statfn_sed[STATFN4].nfn; n++) {
    if(sfn_sed[n].step == step) {
      sfn_sed[n].statfn_new(trs, &trs->w_sed, trs->tr_sed, sfn_sed[n].n,
			    trs->botk_sed, trs->topk_sed, trs->map_sed,
			    sfn_sed[n].data);
    }
  }
  
  sfn_sed = trs->statfn_sed[STATFN5].statfn;
  for(n = 0; n < trs->statfn_sed[STATFN5].nfn; n++) {
    if(sfn_sed[n].step == step) {
      sfn_sed[n].statfn_new(trs, &trs->w_sed, trs->tr_sed, sfn_sed[n].n,
			    trs->botk_sed, trs->topk_sed, trs->map_sed,
			    sfn_sed[n].data);
    }
  }
  
  sfn_sed = trs->statfn_sed[STATFN6].statfn;
  for(n = 0; n < trs->statfn_sed[STATFN6].nfn; n++) {
    if(sfn_sed[n].step == step) {
      sfn_sed[n].statfn_new(trs, &trs->w_sed, trs->tr_sed, sfn_sed[n].n,
			    trs->botk_sed, trs->topk_sed, trs->map_sed,
			    sfn_sed[n].data);
    }
  }
  
  /*
   * Iterate through the 2d functions
   */
  
  sfn2d = trs->statfn2d[STATFN0].statfn;
  for(n = 0; n < trs->statfn2d[STATFN0].nfn; n++)
    {
      if(sfn2d[n].step == step)
	sfn2d[n].statfn(trs,&(sfn2d[n].n),sfn2d[n].data);
    }
  
  sfn2d = trs->statfn2d[STATFN1].statfn;
  for(n = 0; n < trs->statfn2d[STATFN1].nfn; n++)
    {
      if(sfn2d[n].step == step)
	sfn2d[n].statfn(trs,&(sfn2d[n].n),sfn2d[n].data);
    }
  
  sfn2d = trs->statfn2d[STATFN2].statfn;
  for(n = 0; n < trs->statfn2d[STATFN2].nfn; n++)
    {
      if(sfn2d[n].step == step)
	sfn2d[n].statfn(trs,&(sfn2d[n].n),sfn2d[n].data);
    }
  
  sfn2d = trs->statfn2d[STATFN3].statfn;
  for(n = 0; n < trs->statfn2d[STATFN3].nfn; n++)
    {
      if(sfn2d[n].step == step)
	sfn2d[n].statfn(trs,&(sfn2d[n].n),sfn2d[n].data);
    }
  
  sfn2d = trs->statfn2d[STATFN4].statfn;
  for(n = 0; n < trs->statfn2d[STATFN4].nfn; n++)
    {
      if(sfn2d[n].step == step)
	sfn2d[n].statfn(trs,&(sfn2d[n].n),sfn2d[n].data);
    }
  
  
  sfn2d = trs->statfn2d[STATFN5].statfn;
  for(n = 0; n < trs->statfn2d[STATFN5].nfn; n++)
    {
      if(sfn2d[n].step == step)
	sfn2d[n].statfn(trs,&(sfn2d[n].n),sfn2d[n].data);
    }
  
  sfn2d = trs->statfn2d[STATFN6].statfn;
  for(n = 0; n < trs->statfn2d[STATFN6].nfn; n++)
    {
  	if(sfn2d[n].step == step)
	  sfn2d[n].statfn(trs,&(sfn2d[n].n),sfn2d[n].data);
    }
}
/* END trs_step()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Read a file providing a list of I,j coordinates                   */
/* Return the number of coordinates read                             */
static int read_section_coordinates(char* fnamein,   /* input fiilename */
				    sectioncoord_t* sectiondata,  /* the section data structure */
				    void* model)   /* the host model*/
{
  FILE* f;
  char** cLines;
  char** outunits;
  char* tline[MAXLINELEN];
  char* token;
  char line[MAXLINELEN];
  int i,j,n,nl,nc,ncc,ni,ni1,ni2=0;
  tracer_info_t* tracer;
  // The number of columns in addition to Time
  int mode_pad = (sectiondata->mode == SECTION_AVERAGE ? SECTION_MODE_PAD : 1);
  // times 4 for conc, squared and 2x for the covariance - see header file
  int num_trs  = (sectiondata->mode == SECTION_AVERAGE ? sectiondata->ntrs*4 : sectiondata->ntrs);
  
  if((f = fopen(fnamein, "r")) == NULL)
  {
    emstag(LFATAL,"tracerstats:read_section_coordinates","Error Opening File - %s",fnamein);
    exit(0);
  }else
    emstag(LMETRIC,"Tracerstats:read_section_coordinates","Opened file %s",fnamein);

  nl = 0;
  outunits = NULL;

  for (n = 0; prm_next_line(line, MAXSTRLEN, f); n++)
    ;
  rewind(f);

  emstag(LTRACE,"tracerstats:read_section_coordinates","read lines %d",n);


  /* allocate memory for the file */
  cLines = c_alloc_2d(MAXSTRLEN,n);

/*  b = fgetc(f); */
  if(n == 0)
  {
    emstag(LFATAL,"tracerstats:read_section_coordinates","Error reading File - %s - NO CONTENT!",fnamein);
    exit(0);
  }

  nl= n;
  nc=ncc=0;

  for (n = 0;  prm_next_line(line, MAXSTRLEN, f); n++)
  {
    if(line != NULL && strlen(line) > 0 &&  line[0] != '#')
    {
      strcpy(cLines[nc],line);
      nc++;
    }
  }
  fclose(f);

  if(nc < 1)
  {
    emstag(LFATAL,"tracerstats:read_section_coordinates","Error reading File - %s - NO Data!",fnamein);
    exit(0);
  }
  sectiondata->data = i_alloc_2d(nc,MAXSECTIONPARAMS);

  for(nl = 0; nl < nc; nl++)
  {
    n = -1;
    ni = parseline(cLines[nl],tline,MAXSECTIONPARAMS);

    // initialise to an invalid value
    sectiondata->data[4][nl] = -1;

    if(ni < 2)
    {
      ncc ++;
     continue;
    }
    ni1 = atoi(tline[0]);
    ni2 = atoi(tline[1]);

    /* assign the i and j     */
    sectiondata->data[0][nl] =  ni1;
    sectiondata->data[1][nl] =  ni2;

    /* use the top and bottom layer index as limits */
    ni1 = i_get_num_wclayers(model)-1;
    ni2 = 0;

    if(ni < 3)
    {
      sectiondata->data[2][nl] = ni1;
      sectiondata->data[3][nl] = ni2;
    } else {
      /*
       * if there are more arguments assign them accordingly, toplayer first
       * match the layer index provided we do allow negative values
       * for the top layer - they will be adjusted in the stat
       * sections function, if it is a vertical exchange - w:
       * subtract it from the top or if that exceeds the bottom, it
       * will be bottom +1 
       * example: top layer 24, bot layer 0 -> i,j,topk,botk
       * =#,#,-100,0 will result in botk = 0,topk = 1, this result in
       * the exchanges the last two water layers 
       *
       * horizontal exchange - e1 & e2:
       * reset to the top layer index
       *
       * 01/10 FR: Added the optional argument for the direction for
       *           each cell. I've allowed for the fact that this
       *           direction can be the 3rd, 4th or 5th argument by
       *           checking for an alpha character
       */
      if (isalpha(*tline[2])) {
	sectiondata->data[4][nl] = (*tline[2] == 'e' ? atoi(&tline[2][1]) : 0);
      } else {
	n = atoi(tline[2]);
	if( n > ni1 )
	  sectiondata->data[2][nl] = ni1;
	else
	  sectiondata->data[2][nl] = n;
      }
      
      if(ni > 3)
      {
	if (isalpha(*tline[3])) {
	  sectiondata->data[4][nl] =  (*tline[3] == 'e' ? atoi(&tline[3][1]):0);
	} else {
	  n = atoi(tline[3]);
	  if(n >= ni1 || n < ni2 )
	    sectiondata->data[3][nl] = ni2;
	  else
	    sectiondata->data[3][nl] = n;
	}
      }
      if (ni > 4) {
	if (isalpha(*tline[4]))
	  sectiondata->data[4][nl] = (*tline[4] == 'e' ? atoi(&tline[4][1]):0);
      }
    }
  }
  emstag(LTRACE,"tracerstats:read_section_coordinates","Read lines %d from file %s ",nc,fnamein);
  sectiondata->nc = (nl - ncc );
  /* initialise the output file */
  f = fopen(sectiondata->fileout,"w");
  if(f == NULL)
  {
    emstag(LFATAL,"tracerstats:tracerstats.c:init_section_out","Failed creating file %s ",sectiondata->fileout);
    exit(0);
  }

  if(sectiondata->outunit== NULL)
    ni =-1;
  else
  {
    ni = 0;
    outunits = malloc(sectiondata->ntrs * sizeof(char*));
    for(i=0 ;i < sectiondata->ntrs;i++)
      outunits[i] = malloc(MAXSTRLEN * sizeof(char));

    token = strtok(sectiondata->outunit,",");
    while(token != NULL)
    {
      strcpy(outunits[ni],token);
      token = strtok(NULL,",");
      ni++;
    }
  }


  /*sectiondata->tscale = 86400.0; */
  /* Write the comment header */
  fprintf(f, "# Time series data from model run - section id: %s \n",sectiondata->name);
  fprintf(f, "# Section profile for tracers, dt %.0f sec, direction %s \n",sectiondata->dt,(sectiondata->dir != SECTIONE1?(sectiondata->dir != SECTIONW?"'e2'":"'w'"):"'e1'"));
  fprintf(f, "# \n");
  fprintf(f, "# Points: i j ktop kbottom\n");
  for(nl = 0; nl < nc; nl++)
    fprintf(f, "# %d %d  %d %d  \n",sectiondata->data[0][nl],sectiondata->data[1][nl],sectiondata->data[2][nl],sectiondata->data[3][nl]);
  fprintf(f, "# \n");

  // Write some sort of a description on what the data represents
  if (sectiondata->mode == SECTION_AVERAGE) {
    fprintf(f, "# SECTION MODE = AVERAGE\n");
    fprintf(f, "# \n");
  } else {
    fprintf(f, "# SECTION MODE = INTEGRATE\n");
    fprintf(f, "# \n");
  }

  fprintf(f,"## COLUMNS %d \n",num_trs+mode_pad);
  fprintf(f,"##\n");
  fprintf(f,"## COLUMN1.name time\n");
  fprintf(f,"## COLUMN1.long_name Time\n");
/* temp */
  fprintf(f,"## COLUMN1.units days since 1990-01-01 00:00:00 +10\n");
/*   fprintf(f,"## COLUMN1.units %s \n",sectiondata->tunit);/ *days since 1990-01-01 00:00:00 +10\n"); */
  fprintf(f,"## COLUMN1.coordinate_type TIME\n");
  fprintf(f,"## COLUMN1.fill_value 0\n");
  fprintf(f,"## COLUMN1.missing_value -999999\n");
  fprintf(f,"##   \n");

  // Write out the velocity header
  if (sectiondata->mode == SECTION_AVERAGE) {
    fprintf(f,"## COLUMN2.name u_vel\n");
    fprintf(f,"## COLUMN2.long_name Mean velocity\n");
    fprintf(f,"## COLUMN2.units m/s\n");
    fprintf(f,"## COLUMN2.fill_value 0\n");
    fprintf(f,"## COLUMN2.missing_value -99999\n");
    fprintf(f,"##\n");
    fprintf(f,"## COLUMN3.name u_abs_vel\n");
    fprintf(f,"## COLUMN3.long_name Mean absolute velocity\n");
    fprintf(f,"## COLUMN3.units m/s\n");
    fprintf(f,"## COLUMN3.fill_value 0\n");
    fprintf(f,"## COLUMN3.missing_value -99999\n");
    fprintf(f,"##\n");
    fprintf(f,"## COLUMN4.name u_sq_vel\n");
    fprintf(f,"## COLUMN4.long_name Mean square velocity\n");
    fprintf(f,"## COLUMN4.units (m/s)^2\n");
    fprintf(f,"## COLUMN4.fill_value 0\n");
    fprintf(f,"## COLUMN4.missing_value -99999\n");
    fprintf(f,"##\n");
    fprintf(f,"## COLUMN5.name area\n");
    fprintf(f,"## COLUMN5.long_name Section area\n");
    fprintf(f,"## COLUMN5.units m^2\n");
    fprintf(f,"## COLUMN5.fill_value 0\n");
    fprintf(f,"## COLUMN5.missing_value -99999\n");
    fprintf(f,"##\n");
  }

  j = 1+mode_pad; // column counter
  for(i = 0 ; i < sectiondata->ntrs  ;i++,j++) {
    tracer = i_get_tracer(model,sectiondata->tnames[i]);
    if(tracer == NULL) {
      fprintf(f,"## COLUMN%d.name %s\n",j, sectiondata->tnames[i]);
      fprintf(f,"## COLUMN%d.long_name %s\n",j, "<not available>");
      emstag(LPANIC,"tracerstats:tracerstats:read_section_coordinates","Failed retrieving tracer %s, does it exist?",sectiondata->tnames[i]);
    } else {
      if (sectiondata->mode == SECTION_AVERAGE) {
	fprintf(f,"## COLUMN%d.name %s_C\n",j, tracer->name);
	fprintf(f,"## COLUMN%d.long_name Mean %s concentration\n",j, tracer->long_name);
      } else {
	fprintf(f,"## COLUMN%d.name %s\n",j, tracer->name);
	fprintf(f,"## COLUMN%d.long_name %s\n",j, tracer->long_name);
      }
    }
    
    fprintf(f,"## COLUMN%d.fill_value 0\n",j);
    fprintf(f,"## COLUMN%d.missing_value -99999\n",j);
    if(ni > i )
      fprintf(f,"## COLUMN%d.units %s\n",j,outunits[i]);
    else if(tracer != NULL)
      fprintf(f,"## COLUMN%d.units %s\n",j,tracer->units);
    
    fprintf(f,"##\n");
    // Handle the extra parts
    if (sectiondata->mode == SECTION_AVERAGE) {
      // squared
      j++;
      if(tracer == NULL) {
	fprintf(f,"## COLUMN%d.name %s_sq\n",j, sectiondata->tnames[i]);
	fprintf(f,"## COLUMN%d.long_name Mean square %s\n",j, "<not available>");
	emstag(LPANIC,"tracerstats:tracerstats:read_section_coordinates","Failed retrieving tracer %s, does it exist?",sectiondata->tnames[i]);
      } else {
	fprintf(f,"## COLUMN%d.name %s_sq\n",j, tracer->name);
	fprintf(f,"## COLUMN%d.long_name Mean square %s\n",j, tracer->long_name);
      }
      fprintf(f,"## COLUMN%d.fill_value 0\n",j);
      fprintf(f,"## COLUMN%d.missing_value -99999\n",j);
      if(ni > i )
	fprintf(f,"## COLUMN%d.units (%s)^2\n",j,outunits[i]);
      else if(tracer != NULL)
	fprintf(f,"## COLUMN%d.units (%s)^2\n",j,tracer->units);
      fprintf(f,"##\n");
      // Conc flux
      j++;
      if(tracer == NULL) {
	fprintf(f,"## COLUMN%d.name %s_F\n",j, sectiondata->tnames[i]);
	fprintf(f,"## COLUMN%d.long_name %s flux\n",j, "<not available>");
	emstag(LPANIC,"tracerstats:tracerstats:read_section_coordinates","Failed retrieving tracer %s, does it exist?",sectiondata->tnames[i]);
      } else {
	fprintf(f,"## COLUMN%d.name %s_F\n",j, tracer->name);
	fprintf(f,"## COLUMN%d.long_name Mean %s flux\n",j, tracer->long_name);
      }
      fprintf(f,"## COLUMN%d.fill_value 0\n",j);
      fprintf(f,"## COLUMN%d.missing_value -99999\n",j);
      if(ni > i )
	fprintf(f,"## COLUMN%d.units %s/m2/s\n",j,outunits[i]);
      else if(tracer != NULL)
	fprintf(f,"## COLUMN%d.units %s/m2/s\n",j,tracer->units);
      fprintf(f,"##\n");
      // U-C covariance
      j++;
      if(tracer == NULL) {
	fprintf(f,"## COLUMN%d.name %s_cv\n",j, sectiondata->tnames[i]);
	fprintf(f,"## COLUMN%d.long_name %s covariance\n",j, "<not available>");
	emstag(LPANIC,"tracerstats:tracerstats:read_section_coordinates","Failed retrieving tracer %s, does it exist?",sectiondata->tnames[i]);
      } else {
	fprintf(f,"## COLUMN%d.name %s_CV\n",j, tracer->name);
	fprintf(f,"## COLUMN%d.long_name Mean %s covariance\n",j, tracer->long_name);
      }
      fprintf(f,"## COLUMN%d.fill_value 0\n",j);
      fprintf(f,"## COLUMN%d.missing_value -99999\n",j);
      fprintf(f,"## COLUMN%d.units none\n",j);
      fprintf(f,"##\n");
    }
  }
  fclose(f);
  
  free(cLines);
  if(ni > -1 )
  {
    for(i = 0 ; i < sectiondata->ntrs  ;i++)
      free(outunits[i]);
    free(outunits);
  }
  emstag(LMETRIC,"tracerstats:read_section_coordinates","Wrote output header lines for file %s ",sectiondata->fileout);
  return (nl - ncc);
}
/* END read_section_coordinates ----------------------------------------*/


static void corr_2d_pre(trs_t *trs, int* n, void* data)
{
  tracer_corr_2d_pre(trs, n, data);
  tracer_cov_2d(trs, n, data);
}


static void corr_2d_post(trs_t *trs, int* n, void* data)
{
  tracer_corr_2d_post(trs, n, data);
}


static void corr_3d_pre(trs_t *trs, int* n, void* data, int c)
{
  tracer_corr_3d_pre(trs, n, data,c);
  tracer_cov_3d(trs, n, data,c);
}


static void corr_3d_post(trs_t *trs, int* n, void* data, int c)
{
  tracer_corr_3d_post(trs, n, data,c);
}


static int decode_args(char *p, char *args[], int maxn)
{
  int len = strlen(p);
  int i;

  for (i = 0; i < len; ++i) {
    if (p[i] == '(' || p[i] == ')' || p[i] == ':')
      p[i] = ' ';
  }

  return parseline(p, args, maxn);
}


/**
 * determine at which step the function should be executed
 *
 * @param st - the string containing either a number between 0 -9 or
 * string the host model will resolve
 * @return step
 */
 static int extract_step(void* model, char* st)
 {
    if(strlen(st) == 1 && isdigit(st[0]))
        return atoi(st);
    else
        return i_get_step(model, st);

 }

/**
 * Extract the int ranges from a mixed string of individual ints comma separated
 * and ranges separated by an underscore
 *  - > moves to lib
 *
 * @param arg - the string containing the content
 * @param vals- reference to the return array of individual ints
 * @return  the number of individual entries found
 * @Depreceated - moves into lib
 */
static int* split_ints(char* arg, int* nv )
{
  char* tok;
  char* tok2;
  char tbuf[MAXSTRLEN];
  char ttbuf[MAXSTRLEN][MAXSTRLEN];
  int lng = strlen(arg);
  int sz=0;
  int nt;
  int idash=0;
  int iind = 0;
  int d,n,nn,nnn;
  int** dashranges;
  int* indvals;
  int* vals;

  /*count how many we have */
  for(n=0;n<lng;n++)
  {
    if(arg[n] == '_')
      idash++;
    else if(arg[n] == ',')
      iind++;
  }
  emstag(LTRACE,"tracerstats:split_ints","Found separators: %d - ranges: %d from %s",iind,idash,arg);

  /* do we have a single value */
  if(iind == 0 && idash == 0)
  {
    vals = i_alloc_1d(1);
    vals[0] = atoi(arg);
    *nv = 1;
    return vals;
  }else if(iind == 0)/* only one range */
  {
    tok = strchr(arg,'_');
    n = atoi(strncpy(tbuf,arg,tok-arg));
    nn = atoi(tok+1);
    if(n==nn) /* silly*/
    {
      vals = i_alloc_1d(1);
      vals[0] = nn;
      *nv = 1;
      return vals;
    }
    if(nn > n)
    {
        vals = i_alloc_1d(nn-n);
        for(nnn=0;n <= nn;n++,nnn++)
        vals[nnn] = n;
    }else
    {
      vals = i_alloc_1d(n-nn);
        for(nnn=0;nn <= n;nn++,nnn++)
        vals[nnn] = nn;
    }

    if(is_log_enabled(LTRACE))
    {
      tbuf[0] = '\0';
      for(nn= 0; nn < nnn ;nn++)
      {
          sprintf(tbuf,"%s %d",tbuf,vals[nn]);
      }
      emstag(LTRACE,"tracerstats:split_ints","Values total of %d as %s ",nnn,tbuf);
    }
    *nv = nnn;
    return vals;
  }else if(idash == 0)/* no range */
  {
      vals = i_alloc_1d(iind+1);
      tok = strtok(arg,",");
      n=0;
      while(tok!= NULL)
      {
        vals[n] = atoi(tok);
        tok = strtok(NULL,",");
        n++;
      }
      *nv = n;
      if(is_log_enabled(LTRACE))
      {
        tbuf[0] = '\0';
        for(nn= 0; nn < n ;nn++)
        {
            sprintf(tbuf,"%s %d",tbuf,vals[nn]);
        }
        emstag(LTRACE,"tracerstats:split_ints","Values total of %d as %s ",n,tbuf);
      }
      return vals;
  }
  /*now we need to do some serious parsing */

  sz = iind + 1 - idash;/* shouldn't be here if sz < 0 */
  dashranges = i_alloc_2d(idash,2);
  indvals = i_alloc_1d(sz);

  n = 0;
  tok2 = strtok(arg,",");
  strcpy(ttbuf[n],tok2);
  n++;

  while(tok2 != NULL)
  {
    tok2 = strtok(NULL,",");
    if(tok2 != NULL)
    {
      strcpy(ttbuf[n],tok2);
      n++;
    }else
      break;
  }

  iind = sz;
  d = 0;
  nnn = 0;
  for(nn=0 ;nn<n;nn++)
  {
    tok = strchr(ttbuf[nn],'_');
    if(tok == NULL)/* single value */
    {
      indvals[nnn]= atoi(ttbuf[nn]);
      nnn++;
    }else
    {
      dashranges[d][1] = atoi(tok+1);
      nt = strlen(tok) - 1;
      ttbuf[nn][nt]='\0';
      dashranges[d][0] = atoi(ttbuf[nn]);
      if(dashranges[d][0] > dashranges[d][1])
      {
        nt = dashranges[d][0];
        dashranges[d][0] = dashranges[d][1];
        dashranges[d][1] = nt;
      }
      sz = sz + abs(dashranges[d][1]-dashranges[d][0])+1;
      d++;
    }
  }
  emstag(LTRACE,"tracerstats:split_ints","Found total vals: %d ",sz);

  vals = i_alloc_1d(sz);
  nnn=0;
  for(n =0;n < d;n++)
  {
    for(; dashranges[n][0] <= dashranges[n][1] ;)
    {
       vals[nnn] =  dashranges[n][0];
       dashranges[n][0]++;
       nnn++;
    }
  }

  for(n =0; n < iind && nnn < sz ;n++)
  {
    vals[nnn] =   indvals[n];
    nnn++;
  }

  tbuf[0] = '\0';
  for(n= 0; n < sz ;n++)
  {
      sprintf(tbuf,"%s %d",tbuf,vals[n]);
  }
  emstag(LTRACE,"tracerstats:split_ints","Values %s ",tbuf);

  free_2d(dashranges);
  free_1d(indvals);

  *nv = sz;
  return vals;
}
/* END split_ints ------------------------------------------------*/


void append_ts(char* filename, double* vals, char** format, int nvals)
{
  FILE* f;
  int i;
  f = fopen(filename,"a");
  if(f == NULL)
  {
    emstag(LERROR,"tracerstats:tracerstats.c:append_ts","Failed writing to file %s ",filename);
    return;
  }
  for(i= 0;i < nvals;i++)
  {
    fprintf(f,format[i],vals[i]);
  }
  fprintf(f,"\n");
  fflush(f);
  fclose(f);
  /*emstag(LTRACE,"tracerstats:tracerstats.c:append_ts","Wrote to file %s nvar %d ",filename,nvals);*/
}


int in_section(void* model, sectioncoord_t* data, int col)
{
  return i_in_window(model,col,data->data[0],data->data[1],data->nc);
}

int* i_get_windowijk(void* model, int col,int* ij)
{
  return i_get_ijk(model,col,ij);
}


static int count_elements(char* string, int del)
{
  int count = 0;
  int i=0;
  if(string == NULL)
    return count;

  while(string[i] != '\0')
  {
    if(i > 0 && string[i] == del)
      count++;
    i++;
  }
  return count+1;
}

static void init_sum(char* names, int* nsumtrs, char** t3dnames, int m)
{
  int i=1;
  char* tok = strtok(names,",");
  strcpy(t3dnames[m], tok);
  
  tok = strtok(NULL,",");
  while(tok != NULL)
    {
      strcpy(t3dnames[m+i],tok);
      tok = strtok(NULL,",");
      i++;
    }
  *nsumtrs=i;
}

// helper function
static double * alloc_results(int num, const char *str)
{
  double *ptr = (double *)malloc(num * sizeof(double));
  if(ptr == NULL) {
    emstag(LFATAL,"tracerstats:tarcerstats:init_sectiondata",
	   "Unable to assign memory for section results %s\n", str);
    exit(0);
  }
  return(ptr);
}

/*
 * This is called once for every window
 */
static void init_sectiondata(sectioncoord_t* sectiondata,
			     char*  names,
			     char*  outscales, 
			     char** t3dnames,
			     int m)
{
  int i=0;
  char* tok = NULL;
  int mode_pad = (sectiondata->mode == SECTION_AVERAGE ? SECTION_MODE_PAD : 1);
  int num_trs  = 0;

  sectiondata->lastdump = 0.0;
  sectiondata->nsteps   = 0;
  sectiondata->ntrs = count_elements(names,(int)',');
  emstag(LDEBUG,"tracerstats:tarcerstats:init_sectiondata",
	 "Considering %d tracers for section %s",
	 sectiondata->ntrs,sectiondata->name);

  num_trs = (sectiondata->mode == SECTION_AVERAGE ? sectiondata->ntrs*4 : sectiondata->ntrs);

  sectiondata->tnames    = malloc(sectiondata->ntrs * sizeof(char*));
  sectiondata->outscales = malloc(sectiondata->ntrs * sizeof(double));
  for (i = 0; i < sectiondata->ntrs; i++)
    sectiondata->tnames[i] = malloc(MAXSTRLEN * sizeof(char));

  sectiondata->trsn = (int *)malloc(sectiondata->ntrs * sizeof(int));
  // Allocate all the results_ arrays
  sectiondata->results = alloc_results(sectiondata->ntrs, "");
  if(sectiondata->mode == SECTION_AVERAGE) {
    sectiondata->results_c       = alloc_results(sectiondata->ntrs, "conc");
    sectiondata->results_sq      = alloc_results(sectiondata->ntrs, "squared");
    sectiondata->results_cusgn_a = alloc_results(sectiondata->ntrs, "cov a");
    sectiondata->results_cusgn_b = alloc_results(sectiondata->ntrs, "cov b");
  } else {
    sectiondata->results_c       = NULL;
    sectiondata->results_sq      = NULL;
    sectiondata->results_cusgn_a = NULL;
    sectiondata->results_cusgn_b = NULL;
  }
  strcpy(sectiondata->tnames[0],strtok(names,","));
  strcpy(t3dnames[m],sectiondata->tnames[0]);
  sectiondata->trsn[0] = m;
  sectiondata->results[0] = 0.0;
  if (sectiondata->mode == SECTION_AVERAGE) {
    sectiondata->results_c[0] = 0.0;
    sectiondata->results_sq[0] = 0.0;
    sectiondata->results_cusgn_a[0] = 0.0;
    sectiondata->results_cusgn_b[0] = 0.0;
  }
  for(i = 1 ;i < sectiondata->ntrs;i++) {
    tok = strtok(NULL,",");
    if(tok == NULL) {
      emstag(LFATAL,"tracerstats:tarcerstats:init_sectiondata","Extracting tracers failed, list for section '%s'does not expand to %d tracers",sectiondata->name,sectiondata->ntrs);
      exit(0);
    }
    strcpy(sectiondata->tnames[i],tok);
    strcpy(t3dnames[m+i],sectiondata->tnames[i]);
    sectiondata->trsn[i] = m+i;
    sectiondata->results[i] = 0.0;
    if (sectiondata->mode == SECTION_AVERAGE) {
      sectiondata->results_c[i] = 0.0;
      sectiondata->results_sq[i] = 0.0;
      sectiondata->results_cusgn_a[i] = 0.0;
      sectiondata->results_cusgn_b[i] = 0.0;
    }
  }
  
  // Initialise velocities
  sectiondata->u = sectiondata->uabs = sectiondata->usq = 0.0;

  if(outscales == NULL) {
    for(i = 0 ;i < sectiondata->ntrs;i++)
      sectiondata->outscales[i] = 1.0;
  } else {
    tok = strtok(outscales,",");
    if(tok != NULL)
      sectiondata->outscales[0] = atof(tok);
    else
      sectiondata->outscales[0] = 1.0;
    for(i = 1 ;i < sectiondata->ntrs;i++) {
      tok = strtok(NULL,",");
      if(tok != NULL)
        sectiondata->outscales[i] = atof(tok);
      else
        sectiondata->outscales[i] = 1.0;
    }
  }

  // plus 1 for time
  sectiondata->format = malloc((num_trs+mode_pad) * sizeof(char*));
  for(i = 0 ;i < num_trs+mode_pad; i++) {
    sectiondata->format[i] = malloc(8 * sizeof(char));
    if (sectiondata->mode == SECTION_AVERAGE) {
      // Average values can be small - maybe we should make this the
      // default anyway? -FR
      if (i) {
	strcpy(sectiondata->format[i],"%.5e  ");
      } else {
	// This suits time a little better
	strcpy(sectiondata->format[i],"%.5f  ");
      }
    } else {
      strcpy(sectiondata->format[i],"%.5f  ");
    }
  }
  
  emstag(LTRACE,"tracerstats:tracerstats:extract_tracer_names","extracted %d tracer names as %s -starting target tracer at %d ",sectiondata->ntrs,names,m);
}

/*
 * Here we create the master section. This function is called only
 * once and is passed in the section from the 1st window only
 */
void* section_create(void* d)
{
  int i, pad;
  sectioncoord_t* sd = (sectioncoord_t*)d;
  sectioncoord_t* section = (sectioncoord_t*) malloc(sizeof(sectioncoord_t));
  int num_trs = 0;

  if(sd != NULL)
  {
    section->dir      = sd->dir;
    section->dt       = sd->dt;
    section->fileout  = (char*) malloc(sizeof(char) * (strlen(sd->fileout)+1));
    strcpy(section->fileout,sd->fileout);

    section->lastdump = sd->lastdump;
    section->time     = sd->time;
    section->ntrs     = sd->ntrs;
    section->startt   = sd->startt;
    section->tscale   = sd->tscale;
    section->mode     = sd->mode;

    num_trs = (section->mode == SECTION_AVERAGE ? section->ntrs*4 : section->ntrs);

    section->tnames = malloc(section->ntrs * sizeof(char *));
    section->results= malloc(section->ntrs * sizeof(double));
    if (section->mode == SECTION_AVERAGE) {
      section->results_c       = malloc(section->ntrs * sizeof(double));
      section->results_sq      = malloc(section->ntrs * sizeof(double));
      section->results_cusgn_a = malloc(section->ntrs * sizeof(double));
      section->results_cusgn_b = malloc(section->ntrs * sizeof(double));
    }
    section->outscales = malloc(section->ntrs * sizeof(double));
    pad = 1; // 1 for time
    if (section->mode == SECTION_AVERAGE)
      pad += SECTION_MODE_PAD-1; // plus 4 for the velocities + area
    section->format = malloc((num_trs + pad) * sizeof(char*));
    
    // format is handled separately
    for (i = 0; i< num_trs+pad; i++) {
      section->format[i] = malloc((strlen(sd->format[i])+1) * sizeof(char));
      strcpy(section->format[i], sd->format[i]);
    }
    // Now for the tracers
    for (i = 0; i < section->ntrs; i++) {
      section->tnames[i] = malloc((strlen(sd->tnames[i])+1) * sizeof(char));
      strcpy(section->tnames[i], sd->tnames[i]);

      section->outscales[i] = sd->outscales[i];
      section->results[i]   = 0.0;
      
      if (section->mode == SECTION_AVERAGE) {
	section->results_c[i]  = 0.0;
	section->results_sq[i] = 0.0;
	section->results_cusgn_a[i] = 0.0;
	section->results_cusgn_b[i] = 0.0;
      }
    }
    // Initialise velocities
    section->u = section->uabs = section->usq = 0.0;
    section->nsteps   = 0;
    section->Tot_area = 0.0;
  }
  emstag(LDEBUG,"trstat:tracerstats:section_create",
	 "cloned sectioncoord succesfully: %s", (sd != NULL?"Yes":"No"));

  return section;
}

/*
 * Gather up all the section bits from multiple windows
 */
void section_gather(int n, void* target, void* origin)
{
  int i;
  trs_t* trs = (trs_t*)origin;
  sectioncoord_t* sectt = (sectioncoord_t*)target;
  sectioncoord_t* sectf = (sectioncoord_t*)trs->domainfn3d[n].data;
  
  // Sum up the data from different sections
  for(i = 0; i < sectf->ntrs;i++) {
    if (sectt->mode == SECTION_AVERAGE) {
      sectt->results[i]    += sectf->results[i];
      sectt->results_c[i]  += sectf->results_c[i];
      sectt->results_sq[i] += sectf->results_sq[i];
      sectt->results_cusgn_a[i] += sectf->results_cusgn_a[i];
      sectt->results_cusgn_b[i] += sectf->results_cusgn_b[i];
    } else {
      sectt->results[i] += sectf->results[i];
    }
    sectf->results[i] = 0.0;
    if (sectf->mode == SECTION_AVERAGE) {
      sectf->results_c[i]       = 0.0;
      sectf->results_sq[i]      = 0.0;
      sectf->results_cusgn_a[i] = 0.0;
      sectf->results_cusgn_b[i] = 0.0;
    }
  }
  // Velocities & area
  if (sectt->mode == SECTION_AVERAGE) {
    sectt->u    += sectf->u;
    sectt->uabs += sectf->uabs;
    sectt->usq  += sectf->usq;
    // reset
    sectf->u = sectf->uabs = sectf->usq = 0.0;
    // Add up the total areas
    sectt->Tot_area += sectf->Tot_area;
    sectf->Tot_area = 0.0;
  }
  
  sectf->lastdump = sectt->lastdump;
  sectt->time = sectf->time;

  // Keep track of timesteps
  sectt->nsteps++;
}

/*
 * If we've gone past the period, then dump to file and reset
 * Divide by the period if mode == 1
 */
void section_scatter(int n, void* origin, void* target)
{
  sectioncoord_t* section = (sectioncoord_t*)origin;
  emstag(LMETRIC,"trstat:tracerstats:section_scatter","view %f at %.2f  start %.2f %s",section->results,(section->time/section->tscale), section->startt,section->name); 
  
  if( (section->time >= section->startt) &&
      (section->time >= (section->lastdump + section->dt)) ) {
    int i,j=0;
    int mode_pad = (section->mode == SECTION_AVERAGE ? SECTION_MODE_PAD : 1);
    int ntrs     = (section->mode == SECTION_AVERAGE ? section->ntrs*4 : section->ntrs);
    int num_cols = ntrs + mode_pad;
    double* tmp  = calloc(num_cols, sizeof(double));
    // Time always comes first
    tmp[j++] = section->time/section->tscale;

    /*
     * NOTE: The order of tmp must be the same as that specified in
     *       the header
     */

    // Add velocities, if (mode)
    if (section->mode == SECTION_AVERAGE) {
      tmp[j++] = section->u    / (section->Tot_area * section->nsteps);
      tmp[j++] = section->uabs / (section->Tot_area * section->nsteps);
      tmp[j++] = section->usq  / (section->Tot_area * section->nsteps);
      tmp[j++] = section->Tot_area;
    }
    // The rest of the tracers
    for(i = 0;i < section->ntrs;i++) {
      if (section->mode == SECTION_AVERAGE) {
	// Average over time
	double C = section->results_c[i]/section->outscales[i] / 
	                            (section->Tot_area * section->nsteps);
	tmp[j++] = C;
	tmp[j++] = section->results_sq[i]/section->outscales[i] / 
	                            (section->Tot_area * section->nsteps);
	tmp[j++] = section->results[i] /section->outscales[i] / 
	              (section->Tot_area * section->nsteps * section->dt);
	// Derive the covariance
	tmp[j++] = (section->results_cusgn_a[i] - 
		    section->results_cusgn_b[i] * C) / 
	  (section->outscales[i] * section->Tot_area * section->nsteps);
      } else {
	// Integration
	tmp[j++] = section->results[i] / section->outscales[i];
      }
    }

    emstag(LTRACE,"trstat:tracerstats:section_scatter","writing %f at %.2f  start %.2f %s",section->results,(section->time/section->tscale), section->startt,section->name);
    // Dump to file
    append_ts(section->fileout, tmp, section->format, num_cols);

    section->lastdump = section->time;
    for(i = 0;i < section->ntrs;i++) {
      section->results[i] = 0.0;
      if (section->mode == SECTION_AVERAGE) {
	section->results_c[i]  = 0.0;
	section->results_sq[i] = 0.0;
	section->results_cusgn_a[i] = 0.0;
	section->results_cusgn_b[i] = 0.0;
      }
    }
    if (section->mode == SECTION_AVERAGE)
      section->u = section->uabs = section->usq = 0.0;

    section->Tot_area = 0.0;

    free(tmp);
    
    // reset number of timesteps
    section->nsteps = 0;
  }
}

/*
 * This is a generic function that allows for either initialising the
 * WC or SED.  Note that not all of the 3d functions have been ported
 * to sediments
 */
void fill_stat3dfcn(trs_t *trs, char  *op_3d, char  **tr_3d,
		    lstatfn3d_t *statfn3d, char **trnames, 
		    int m, stat3d_work_arrays* w_arr, int *st_type)
{
  if (strcmp(op_3d, "mean") == 0) 
    {
      tm_scale_to_secs(tr_3d[2], &w_arr->w1[m]);
      w_arr->w2[m] = 0.0;
      st_type[m] = MEAN;
      strcpy(trnames[m], tr_3d[0]);
      strcpy(trnames[m + 1],tr_3d[1]);
      statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].type = MEAN;
      statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].statfn = NULL;
      statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].statfn_new = tracer_mean_3d;
      statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].n = m;
      statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].data = NULL;
      statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].step = 
	extract_step(trs->model, tr_3d[3]);
      trs->nprestep[statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].step]++;
      
      statfn3d[STATFN1].nfn++;
    }
  else if (strcmp(op_3d, "run_mean") == 0) 
    {
      tm_scale_to_secs(tr_3d[2], &w_arr->w1[m]);
      w_arr->w2[m] = 0.0;
      st_type[m] = RUN_MEAN;
      strcpy(trnames[m], tr_3d[0]);
      strcpy(trnames[m + 1],tr_3d[1]);
      statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].type = RUN_MEAN;
      statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].statfn = NULL;
      // xxx
      statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].statfn_new = tracer_run_mean_3d;
      statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].n = m;
      statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].data = NULL;
      statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].step = 
	extract_step(trs->model, tr_3d[3]);
      trs->nprestep[statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].step]++;
      
      statfn3d[STATFN1].nfn++;
    } 
  else if (strcmp(op_3d, "rmse") == 0) 
    {
      tm_scale_to_secs(tr_3d[2], &w_arr->w1[m]);
      w_arr->w2[m] = 0.0;
      st_type[m] = RMSE;
      strcpy(trnames[m], tr_3d[0]);
      strcpy(trnames[m + 1],tr_3d[1]);
      strcpy(trnames[m + 2],tr_3d[2]);
      statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].type = RMSE;
      statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].statfn = NULL;
      statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].statfn_new = tracer_rmse_3d;
      statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].n = m;
      statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].data = NULL;
      statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].step = 
	extract_step(trs->model, tr_3d[3]);
      trs->nprestep[statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].step]++;
      
      statfn3d[STATFN0].nfn++;
    } 
  else if (strcmp(op_3d, "max") == 0) 
    {
      tm_scale_to_secs(tr_3d[2], &w_arr->w1[m]);
      w_arr->w2[m] = 0.0;
      st_type[m] = MAX;
      strcpy(trnames[m],tr_3d[0]);
      strcpy(trnames[m + 1],tr_3d[1]);

      statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].type = MAX;
      statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].statfn_new = tracer_max_3d;
      statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].statfn = NULL;
      statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].n = m;
      statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].data = NULL;
      statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].step = extract_step(trs->model, tr_3d[3]);
      trs->nprestep[statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].step]++;
      
      statfn3d[STATFN1].nfn++;
    } 
  else if (strcmp(op_3d, "min") == 0) 
    {
      tm_scale_to_secs(tr_3d[2], &w_arr->w1[m]);
      w_arr->w2[m] = 0.0;
      st_type[m] = MIN;
      strcpy(trnames[m],tr_3d[0]);
      strcpy(trnames[m + 1],tr_3d[1]);
      statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].type = MIN;
      statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].statfn_new = tracer_min_3d;
      statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].statfn = NULL;
      statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].n = m;
      statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].data = NULL;
      statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].step = extract_step(trs->model, tr_3d[3]);
      trs->nprestep[statfn3d[STATFN1].statfn[statfn3d[STATFN1].nfn].step]++;
      
      statfn3d[STATFN1].nfn++;
    } 
  else if (strcmp(op_3d, "stdev") == 0) 
    {
      st_type[m] = STDEV;
      strcpy(trnames[m],tr_3d[0]);
      statfn3d[STATFN4].statfn[statfn3d[STATFN4].nfn].type = STDEV;
      statfn3d[STATFN4].statfn[statfn3d[STATFN4].nfn].statfn_new = tracer_stdev_3d;
      statfn3d[STATFN4].statfn[statfn3d[STATFN4].nfn].statfn = NULL;
      statfn3d[STATFN4].statfn[statfn3d[STATFN4].nfn].n = m;
      statfn3d[STATFN4].statfn[statfn3d[STATFN4].nfn].data = NULL;
      statfn3d[STATFN4].statfn[statfn3d[STATFN4].nfn].step = extract_step(trs->model, tr_3d[3]);
      trs->nprestep[statfn3d[STATFN4].statfn[statfn3d[STATFN4].nfn].step]++;
      
      statfn3d[STATFN4].nfn++;
    } 
  else if (strcmp(op_3d, "variance") == 0) 
    {
      st_type[m] = VAR;
      strcpy(trnames[m],tr_3d[0]);
      statfn3d[STATFN3].statfn[statfn3d[STATFN3].nfn].type = VAR;
      statfn3d[STATFN3].statfn[statfn3d[STATFN3].nfn].statfn_new = tracer_var_3d;
      statfn3d[STATFN3].statfn[statfn3d[STATFN3].nfn].statfn = NULL;
      statfn3d[STATFN3].statfn[statfn3d[STATFN3].nfn].n = m;
      statfn3d[STATFN3].statfn[statfn3d[STATFN3].nfn].data = NULL;
      statfn3d[STATFN3].statfn[statfn3d[STATFN3].nfn].step = extract_step(trs->model, tr_3d[3]);
      trs->nprestep[statfn3d[STATFN3].statfn[statfn3d[STATFN3].nfn].step]++;
      
      statfn3d[STATFN3].nfn++;
    } 
  else if (strcmp(op_3d, "sum") == 0) 
    {
      int* nsumtrs =malloc(sizeof(int));
      st_type[m] = SUM;
      strcpy(trnames[m],tr_3d[0]);
      statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].type = SUM;
      statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].statfn_new = tracer_sum_3d;
      statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].statfn = NULL;
      statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].n = m;
      statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].data = nsumtrs;
      statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].step = extract_step(trs->model, tr_3d[3]);
      trs->nprestep[statfn3d[STATFN0].statfn[statfn3d[STATFN0].nfn].step]++;
      init_sum(tr_3d[1],nsumtrs, trnames, m+1);
      statfn3d[STATFN0].nfn++;
    }
}


// EOF
