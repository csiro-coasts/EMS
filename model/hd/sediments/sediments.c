/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/sediments/sediments.c
 *  
 *  Description:
 *  Sediment module interface
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: sediments.c 7484 2024-02-13 23:53:02Z mar644 $
 *
 */

#include <ctype.h>
#include "hd.h"

#if defined(HAVE_OMP)
#include <omp.h>
#endif

#if defined(HAVE_SEDIMENT_MODULE)

char *SEDCLASS[] = {
  "Gravel",
  "Sand",
  "Silt",
  "Clay",
  "Mud",
  "FineSed",
  "Dust",
  "CarbSand",
  "MetalP",
  "Mud-mineral",
  "Mud-carbonate",
  "Sand-mineral",
  "Sand-carbonate",
  "Gravel-mineral",
  "Gravel-carbonate"
};
#define NUM_SED_VARS (int)(sizeof(SEDCLASS)/sizeof(char *))

char *SEDNAME3D[] = {
  "tss",
  "porosity"
};
#define NUM_SED_VARS_3D (int)(sizeof(SEDNAME3D)/sizeof(char *))

char *SEDNAME2D[] = {
  "ustrcw_skin",
  "depth_sedi",
"erdepflux_total_ac"
};
#define NUM_SED_VARS_2D (int)(sizeof(SEDNAME2D)/sizeof(char *))


/*
 * The sediment autotracers are defined using the following definition
 */
typedef struct {
  char *name;
  char *units;
  int type;
  double fill_sed;
  double psize;
  double b_dens;
  double i_conc;

  //2019
  double css_erosion;
  double css_deposition;

  double svel;
  int diagn;
  int advect;
  int diffuse;
  int cohesive;
  int calcvol;
  int floc;
  int resuspend;
  int deposit;
  int obc;
} sed_def_t;

int sinterface_put_ustrcw(void* hmodel, int c, double ustrcw);

// Local functions
/*static void *private_data_copy_sed(void *src);*/
static trinfo_priv_sed_t *get_private_data_sed(tracer_info_t *tr);
static void sed_defaults_std(tracer_info_t *tracer, char *trname);
static void sed_defaults_est(tracer_info_t *tracer, char *trname);
static void sed_defaults_shf(tracer_info_t *tracer, char *trname);
static void sed_defaults_bsc(tracer_info_t *tracer, char *trname);
static void init_tracer_atts_sed(tracer_info_t *tracer, char *trname, 
				 sed_def_t *sed_def, char *sed_def_name);

//2015
int sinterface_getshipfile(FILE* prmfd, char *shipfile);

//2016
double sinterface_erflux_scale(FILE* prmfd);

// 2019 tracer-specific css erosion and css deposition
// (to overwrire default values when present in the tracer specifications)
double sinterface_get_css_erosion(void* model, char *name);
double sinterface_get_css_deposition(void* model, char *name);

//NMY 2024 Amanda
int sinterface_getflocfile(FILE* prmfd, char *flocfile);

/*-------------------------------------------------------------------*/
/* Sediment specific interface routines
int sinterface_getcss_er_val(FILE* prmfd, double *css_er_val)
int sinterface_getcss_er_depth(FILE* prmfd, double *css_er_depth)
void sinterface_putgridz_sed(void* hmodel, int c, double *gridz_sed)
void sinterface_putcellz_sed(void* hmodel, int c, double *cellz_sed)
void sinterface_putgriddz_wc(void* hmodel, int c, double *gridz_wc)
void sinterface_putcellz_wc(void* hmodel, int c, double *cellz_wc)
void sinterface_putdz_wc(void* hmodel, int c, double *dz_wc)
void sinterface_puttopz_wc(void* hmodel, int c, double topz_wc)
void sinterface_putbotz_wc(void* hmodel, int c, double botz_wc)
int sinterface_put_ustrcw(void* hmodel, int c, double ustrcw)
int sinterface_getverbose_sed(FILE* prmfd)
int sinterface_getgeomorph(FILE* prmfd)
int sinterface_getconsolidate(FILE* prmfd)
double sinterface_getfinpor_sed(FILE* prmfd)
double sinterface_getconsolrate(FILE* prmfd)
double sinterface_getmaxthicksed(FILE* prmfd)
int sinterface_getcssmode(FILE* prmfd)
double sinterface_getcss(FILE* prmfd)
double sinterface_getcssdep(FILE* prmfd)
int sinterface_getflocmode(FILE* prmfd)
double sinterface_getflocprm1(FILE* prmfd)
double sinterface_getflocprm2(FILE* prmfd)
int sinterface_getbblnonlinear(FILE* prmfd)
int sinterface_gethindered_svel_patch(FILE* prmfd)
int sinterface_gethindered_svel(FILE* prmfd)
double sinterface_reef_scale_depth(FILE* prmfd)
int sinterface_getcalcripples(FILE* prmfd)
double sinterface_getcssscale(FILE* prmfd)
double sinterface_getphysriph(FILE* prmfd)
double sinterface_getphysripl(FILE* prmfd)
double sinterface_getbioriph(FILE* prmfd)
double sinterface_getbioripl(FILE* prmfd)
double sinterface_getbiodens(FILE* prmfd)
double sinterface_getmaxbiodepth(FILE* prmfd)
double sinterface_getbi_dissol_kz(FILE* prmfd)
double sinterface_getbt_partic_kz(FILE* prmfd)
double sinterface_getbi_dissol_kz_i(FILE* prmfd)
double sinterface_getbt_partic_kz_i(FILE* prmfd)
char sinterface_getbiosedprofile(FILE* prmfd)
double sinterface_getz0(FILE* prmfd)
double sinterface_getquad_bfc(FILE* prmfd)
int sinterface_getshipfile(FILE* prmfd, char *shipfile)
double sinterface_erflux_scale(FILE* prmfd)
FILE* si_getparamfile_tracer(FILE* fp)
FILE* si_getparamfile_sed(FILE *fp)
double sinterface_get_psize(void* model, char *name)
double sinterface_get_b_dens(void* model, char *name)
double sinterface_get_i_conc(void* model, char *name)
double sinterface_get_svel(void* model, char *name)
void sinterface_get_svel_name(void* model, char *name, char *sname)
int sinterface_get_cohesive(void* model, char *name)
int sinterface_get_resuspend(void* model, char *name)
int sinterface_get_deposit(void* model, char *name)
int sinterface_get_floc(void* model, char *name)
int sinterface_get_calcvol(void* model, char *name)
int sinterface_get_adsorb(void* model, char *name)
void sinterface_get_carriername(void* model, char *name,char *carriername )
void sinterface_get_dissolvedname(void* model, char *name, char *dissolvedname)
double sinterface_get_adsorb_kd(void* model, char *name)
double sinterface_get_adsorb_rate(void* model, char *name)
*/

/*-------------------------------------------------------------------*/
/* Generic interface routines                                        *
/*-------------------------------------------------------------------*/
extern void gi_set_errfn_warn(void *hmodel);
extern int ginterface_getdzinit_sed(FILE* prmfd, double *dz_sed, int sednzc);
extern int  ginterface_getnumberwclayers(void* model);
extern int  ginterface_getnumbersedlayers(void* model);
extern int ginterface_getnumbercolumns(void* model);
extern int ginterface_getnumberoftracers(void* hmodel);
extern int ginterface_getnumberofBtracer(void* hmodel);
extern void ginterface_getnameofBtracer(void* hmodel, int n, char *tracername);
extern double ginterface_getvalueofBtracer(void* hmodel, int n, int c);
extern double *ginterface_getpointerBtracer(void* hmodel, int n, int c);
extern int gi_gettracernames(void* hmodel, char **tracername);
extern int gi_getmap2hdsedtracer(void *hmodel, int ns, char *name);
extern int gi_getmap2hdsedtracer(void *hmodel, int ns, char *name);
extern double ginterface_get_fillvalue_wc(void* model, char *name);
extern double ginterface_get_fillvalue_sed(void* model, char *name);
extern int ginterface_get_diagn(void* model, char *name);
extern int ginterface_get_dissol(void* model, char *name);
extern int ginterface_get_partic(void* model, char *name);
extern int ginterface_get_diffuse(void* model, char *name);
extern double ginterface_get_decay(void* model, char *name);
extern void ginterface_get_tracerunits(void* model, char *name, char *units);
extern int ginterface_get_trans_num_omp(void *model);
extern void ginterface_gettheta(void* hmodel, double *theta, int ncol);
extern void gi_set_errfn_warn(void *hmodel);
extern int ginterface_getcelli(void* hmodel, int c);
extern int ginterface_getcellj(void* hmodel, int c);
extern double ginterface_getcellarea(void* hmodel, int c);
extern double ginterface_getwind1(void* hmodel, int c);
extern double ginterface_getwind2(void* hmodel, int c);
extern double ginterface_getwindspeed(void* hmodel, int c);
extern double ginterface_getwave_ub(void* hmodel, int c);
extern double ginterface_getwave_period(void* hmodel, int c);
extern double ginterface_getwave_dir(void* hmodel, int c);
extern int ginterface_gettopk_wc(void* hmodel, int c);
extern int ginterface_getbotk_wc(void* hmodel, int c);
extern double ginterface_getbotz_wc(void* hmodel, int c);
extern double ginterface_gettopz_wc(void* hmodel, int c);
extern void ginterface_getdz_wc(void* hmodel, int c, double *dz_wc);
extern void ginterface_getu1_wc(void* hmodel, int c, double *u1_wc);
extern void ginterface_getu2_wc(void* hmodel, int c, double *u2_wc);
extern void ginterface_getkz_wc(void* hmodel, int c, double *Kz_wc);
extern double ginterface_get_srf_flux(void* hmodel, char *name, int c);
extern double ginterface_getmodeltime(void* hmodel);
extern double ginterface_getmodeltimestep(void* hmodel);
extern void ginterface_getgridz_sed(void* hmodel, int c, double *gridz_sed);
extern int gi_getmap2hdwctracer(void *hmodel, int ns, char *name);
extern int gi_getmap2hdsedtracer(void *hmodel, int ns, char *name);

/*-------------------------------------------------------------------*/
/* Re-directed interface routines. These should be replaced with     */
/* direct calls to the generic interface routines, bypassing         */
/* wrappers where possible.                                          */
/*-------------------------------------------------------------------*/
double sinterface_getmodeltimestep(void* hmodel) {
  return ginterface_getmodeltimestep(hmodel);
}
int  sinterface_getnumberwclayers(void* model) {
  return ginterface_getnumberwclayers(model);
}
int  sinterface_getnumbersedlayers(void* model) {
  return ginterface_getnumbersedlayers(model);
}
int sinterface_getnumbercolumns(void* model) {
  ginterface_getnumbercolumns(model);
}
int sinterface_getnumberoftracers(void* hmodel) {
  return ginterface_getnumberoftracers(hmodel);
}
int sinterface_getdzinit_sed(FILE* prmfd, double *dz_sed, int sednzc) {
  return ginterface_getdzinit_sed(prmfd, dz_sed, sednzc);
}
int sinterface_getnumberofBtracer(void* hmodel) {
  return ginterface_getnumberofBtracer(hmodel);
}
void sinterface_getnameofBtracer(void* hmodel, int n, char *tracername) {
  ginterface_getnameofBtracer(hmodel, n, tracername);
}
double sinterface_getvalueofBtracer(void* hmodel, int n, int c) {
  return ginterface_getvalueofBtracer(hmodel, n, c);
}
double *sinterface_getpointerBtracer(void* hmodel, int n, int c) {
  return ginterface_getpointerBtracer(hmodel, n, c);
}
int si_gettracernames(void* hmodel, char **tracername) {
  return gi_gettracernames(hmodel, tracername);
}
int si_getmap2hdwctracer(void *hmodel, int ns, char *name) {
  return gi_getmap2hdwctracer(hmodel, ns, name);
}
int si_getmap2hdsedtracer(void *hmodel, int ns, char *name) {
  return gi_getmap2hdsedtracer(hmodel, ns, name);
}
double sinterface_get_fillvalue_wc(void* model, char *name) {
  return ginterface_get_fillvalue_wc(model, name);
}
double sinterface_get_fillvalue_sed(void* model, char *name) {
  return ginterface_get_fillvalue_sed(model, name);
}
int sinterface_get_diagn(void* model, char *name) {
  return ginterface_get_diagn(model, name);
}
int sinterface_get_dissol(void* model, char *name) {
  return ginterface_get_dissol(model, name);
}
int sinterface_get_partic(void* model, char *name) {
  return ginterface_get_partic(model, name);
}
int sinterface_get_diffuse(void* model, char *name) {
  return ginterface_get_diffuse(model, name);
}
double sinterface_get_decay(void* model, char *name) {
  return ginterface_get_decay(model, name);
}
void sinterface_get_tracerunits(void* model, char *name, char *units) {
  ginterface_get_tracerunits(model, name, units);
}
#if defined(HAVE_OMP)
int sinterface_get_trans_num_omp(void *model) {
  return ginterface_get_trans_num_omp(model);
}
#endif
void sinterface_gettheta(void* hmodel, double *theta, int ncol) {
  ginterface_gettheta(hmodel, theta, ncol);
}
void si_set_errfn_warn(void *hmodel) {
  gi_set_errfn_warn(hmodel);
}
int sinterface_getcelli(void* hmodel, int c) {
  return ginterface_getcelli(hmodel, c);
}
int sinterface_getcellj(void* hmodel, int c) {
  return ginterface_getcellj(hmodel, c);
}
double sinterface_getcellarea(void* hmodel, int c) {
  return ginterface_getcellarea(hmodel, c);
}
double sinterface_getwind1(void* hmodel, int c) {
  return ginterface_getwind1(hmodel, c);
}
double sinterface_getwind2(void* hmodel, int c) {
  return ginterface_getwind2(hmodel, c);
}
double sinterface_getwindspeed(void* hmodel, int c) {
  return ginterface_getwindspeed(hmodel, c);
}
double sinterface_getwave_ub(void* hmodel, int c) {
  return ginterface_getwave_ub(hmodel, c);
}
double sinterface_getwave_period(void* hmodel, int c) {
  return ginterface_getwave_period(hmodel, c);
}
double sinterface_getwave_dir(void* hmodel, int c) {
  return ginterface_getwave_dir(hmodel, c);
}
int sinterface_gettopk_wc(void* hmodel, int c) {
  return ginterface_gettopk_wc(hmodel, c);
}
int sinterface_getbotk_wc(void* hmodel, int c) {
  return ginterface_getbotk_wc(hmodel, c);
}
double sinterface_getbotz_wc(void* hmodel, int c) {
  return ginterface_getbotz_wc(hmodel, c);
}
double sinterface_gettopz_wc(void* hmodel, int c) {
  return ginterface_gettopz_wc(hmodel, c);
}
void sinterface_getdz_wc(void* hmodel, int c, double *dz_wc) {
  ginterface_getdz_wc(hmodel, c, dz_wc);
}
void sinterface_getu1_wc(void* hmodel, int c, double *u1_wc) {
  ginterface_getu1_wc(hmodel, c, u1_wc);
}
void sinterface_getu2_wc(void* hmodel, int c, double *u2_wc) {
  ginterface_getu2_wc(hmodel, c, u2_wc);
}
void sinterface_getkz_wc(void* hmodel, int c, double *Kz_wc) {
  ginterface_getkz_wc(hmodel, c, Kz_wc);
}
double sinterface_get_srf_flux(void* hmodel, char *name, int c) {
  return ginterface_get_srf_flux(hmodel, name, c);
}
double sinterface_getmodeltime(void* hmodel) {
  return ginterface_getmodeltime(hmodel);
}
void sinterface_getgridz_sed(void* hmodel, int c, double *gridz_sed) {
  ginterface_getgridz_sed(hmodel, c, gridz_sed);
}

// 2021
char* _itoa(int value, char* str, int base);
void strreverse(char* begin, char* end);

static int log_DVM = 0;
/*-------------------------------------------------------------------*/
/* Sediment transport step                               */
/*-------------------------------------------------------------------*/
void sed_step(geometry_t *window)
{
  win_priv_t *wincon = window->wincon;
  window_t   *windat = window->windat;
  double *nu = wincon->w4;
  double *nv = wincon->w5;
  double *nw = wincon->w10;
  int cc, c, nomp = 1;
  int sed_nstep;
  int ncols = wincon->vcs2;
    
#if defined(HAVE_OMP)
  nomp = wincon->trans_num_omp;
#endif

  /*---------------------------*/
  /* Do the sediment transport */
  /*---------------------------*/

  /* Mask sediment cells to exclude */
  process_cell_mask(window, window->csed, window->ncsed);

  /* Calculate cell centered velocities */
  vel_center(window, windat, wincon, nu, nv, nw);

  // Set global
  wincon->sediment->msparam->ncol = ncols;

  /* Do this outside the loop - pulled out from hd2sed */
  sed_nstep = (wincon->sediment->msparam->nstep += 1);
  wincon->sediment->msparam->t = ginterface_getmodeltime(window);
  //     fprintf(stderr,"%f \n",  param->t);

  /*
   * Loop over all wet columns
   *
   * Note: This should be in sync with the number of columns in ecology
   *       see einterface_getnumbercolumns
   */
#if defined(HAVE_OMP)
#pragma omp parallel for private(c)
#endif
  for (cc = 1; cc <= ncols; cc++) {
    sed_column_t *sm = wincon->sediment->sm;
    c = wincon->s2[cc];
    
    /* Exclude process points if required */
    if (wincon->c2[window->m2d[c]]) continue;
    
    /* Re-allocate column structure, if parallel */
#if defined(HAVE_OMP)
    if (nomp > 1)
      sm = wincon->sediment->smn[omp_get_thread_num()];
#endif

    /* Set column index */
    sm->col_number = cc;

    /* Reset error state for this column */
    i_set_error(window, cc, LNONE, NULL);
    
    /* Data flow from 3-d hydromodel to 1-d sediment */
    hd2sed(wincon->sediment, sm, c);
    
    /* Calculate bottom boundary layer */
    if (wincon->gint_errorf[cc] == LNONE)
      bbl(wincon->sediment, sm);
    
    /* Calculate vertical diffusion, settling, resuspension and deposition */
    if (wincon->gint_errorf[cc] == LNONE)
      vert_transport(wincon->sediment, sm);
    
    /* Check limits for instability handling */
    sed_limits(wincon->sediment, sm, c);
    
    if (wincon->gint_errorf[cc] == LNONE)
      /* All good : Data flow from 1-d sediment to 3-d hydromodel */
      sed2hd(wincon->sediment, sm, c);
    else {
      int c2 = window->m2d[c];
      /* Sediment error */
      hd_warn("Sediment library error @ (%d %d) : %5.2f days\n", 
	      window->s2i[c2], window->s2j[c2], windat->days);
      windat->sederr[c2] = 100.0 * ((double)sed_nstep * 0.01 * windat->sederr[c2] + 1.0) /
	(double)(sed_nstep + 1.0);
      if (wincon->gint_errfcn == LWARN)
	hd_warn(wincon->gint_error[cc]);
      if (wincon->gint_errfcn == LFATAL)
	hd_quit(wincon->gint_error[cc]);
    }
  }
}

/**/


		       /*************************** GEOMETRY */



////// 2010

/* critical shear stress varying with depth*/
int sinterface_getcss_er_val(FILE* prmfd, double *css_er_val)
{
    int k;
    int sednz=0;
    double *tmp;
    prm_read_darray(prmfd, "CSS_ER", &css_er_val, &sednz);
    /* change order of numbering in dz_sed[k] */
    tmp = d_alloc_1d(sednz);
    for (k=0;k<sednz;k++)
	tmp[k]=css_er_val[k];
    for (k=0;k<sednz;k++)
	css_er_val[k] = tmp[sednz-1-k];
    d_free_1d(tmp);
    return sednz;
}
int sinterface_getcss_er_depth(FILE* prmfd, double *css_er_depth)
{
    int k;
    int sednz=0;
    double *tmp;
    prm_read_darray(prmfd, "CSS_ER_DEPTH", &css_er_depth, &sednz);
    /* change order of numbering in dz_sed[k] */
    tmp = d_alloc_1d(sednz);
    for (k=0;k<sednz;k++)
	tmp[k]=css_er_depth[k];
    for (k=0;k<sednz;k++)
	css_er_depth[k] = tmp[sednz-1-k];
    d_free_1d(tmp);
    return sednz;
}


/* Write sediment grid (sed2hd) */
void sinterface_putgridz_sed(void* hmodel, int c, double *gridz_sed)
{
    geometry_t* window = (geometry_t*) hmodel;
    int c2=window->m2d[c];
    int bot_k= 0;
    int top_k= window->sednz-1;
    int k;
    for(k=bot_k;k<=top_k+1;k++)
	    window->gridz_sed[k][c2] = gridz_sed[k];
}
void sinterface_putcellz_sed(void* hmodel, int c, double *cellz_sed)
{
    geometry_t* window = (geometry_t*) hmodel;
    int c2=window->m2d[c];
    int bot_k= 0;
    int top_k= window->sednz-1;
    int k;
    for(k=bot_k;k<=top_k;k++)
	    window->cellz_sed[k][c2] = cellz_sed[k];
}

/** geomorph: empty functions */
void sinterface_putgriddz_wc(void* hmodel, int c, double *gridz_wc)
{ }
void sinterface_putcellz_wc(void* hmodel, int c, double *cellz_wc)
{ }
void sinterface_putdz_wc(void* hmodel, int c, double *dz_wc)
{ }
void sinterface_puttopz_wc(void* hmodel, int c, double topz_wc)
{ }
void sinterface_putbotz_wc(void* hmodel, int c, double botz_wc)
{ }

/*************************************** FORCING (hd2sed.c) */


int sinterface_put_ustrcw(void* hmodel, int c, double ustrcw)
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  win_priv_t *wincon=window->wincon;
  int c2;
  if(wincon->waves & ORBITAL) {
  	c2=window->m2d[c];
  	windat->ustrcw[c2] = ustrcw;
        return 1;
  } else
    return 0;
}
/**/


/*************************** SEDIMENT PROCESS PARAMETERS (sed_init.c)*/


/*************************** SEDIMENT PROCESS PARAMETERS (sed_init.c)*/
int sinterface_getverbose_sed(FILE* prmfd)
{
    int variable;
    if( prm_read_int(prmfd,"VERBOSE_SED",&variable) <= 0)
	variable = 0;
    return variable;
}
int sinterface_getgeomorph(FILE* prmfd)
{
    int geomorph;
    if( prm_read_int(prmfd,"GEOMORPH",&geomorph) <= 0)
	geomorph = 0;
    return geomorph;
}
int sinterface_getconsolidate(FILE* prmfd)
{
    int variable;
    if( prm_read_int(prmfd,"CONSOLIDATE",&variable) <= 0)
	variable = 0;
    return variable;
}
double sinterface_getfinpor_sed(FILE* prmfd)
{
    double variable;
    if( prm_read_double(prmfd,"FINPOR_SED",&variable) <= 0)
	variable = 0.4;
    return variable;
}
double sinterface_getconsolrate(FILE* prmfd)
{
    double variable;
    if( prm_read_double(prmfd,"CONSOLRATE",&variable) <= 0)
	variable = 1.e+13;
    return variable;
}
double sinterface_getmaxthicksed(FILE* prmfd)
{
    double variable;
    if( prm_read_double(prmfd,"MAX_THICK_SED",&variable) <= 0)
	variable = 100;
    return variable;
}
int sinterface_getcssmode(FILE* prmfd)
{
    int v=0;
    char buf[MAXSTRLEN];
    if( prm_read_char(prmfd,"CSS_ER_MODE",buf) > 0) {
	if(strcmp(buf,"CSS0")==0)
	    v=0;
	else if (strcmp(buf,"CSS1")==0)
	    v=1;
	else if (strcmp(buf,"CSS2")==0)
	    v=2;
	else if (strcmp(buf,"CSS3")==0)
	    v=3;
	else if (strcmp(buf,"CSS4")==0)
	    v=4;
    }
    else
    {
	v=0;
    }
    return v;
}
double sinterface_getcss(FILE* prmfd)
{
    double v;
    if(prm_read_double(prmfd,"CSS_ER",&v) <= 0.)
	v = 0.2;
    return v;
}
double sinterface_getcssdep(FILE* prmfd)
{
    double v;
    if(prm_read_double(prmfd,"CSS_DEP",&v) <= 0.)
	v = 1.e+13;
    return v;
}
int sinterface_getflocmode(FILE* prmfd)
{
 int v=0;
 char buf[MAXSTRLEN];
 if( prm_read_char(prmfd,"FLOC_MODE",buf) > 0) {
     if(strcmp(buf,"FLOC0")==0)
	 v=0;
     else if (strcmp(buf,"FLOC1")==0)
	 v=1;
     else if (strcmp(buf,"FLOC2")==0)
	 v=2;
     else if (strcmp(buf,"FLOC3")==0)
	 v=3;
     else if (strcmp(buf,"FLOC4")==0)
	 v=4;
     else if (strcmp(buf,"FLOC5")==0) // NMY 2024 Amanda
         v=5;
 }
 else
 {
     v=0;
 }
 return v;
}


// NMY 2024 Amanda
// read floc d50 lookup array
int sinterface_getflocfile(FILE* prmfd, char *flocfile)
{
    int v=0;
    // char buf[MAXSTRLEN];
    if( prm_read_char(prmfd,"FLOC_FILE",flocfile) > 0) {
      v=1;
    }
    else
    {
        v=0;
    }
    return v;
}



double sinterface_getflocprm1(FILE* prmfd)
{
    double v;
    if(prm_read_double(prmfd,"FLOC_PRM1",&v) <= 0.)
	v = -3.4e-4;
    return v;
}
double sinterface_getflocprm2(FILE* prmfd)
{
    double v;
    if(prm_read_double(prmfd,"FLOC_PRM2",&v) <= 0.)
	v = 3.0;
    return v;
}
int sinterface_getbblnonlinear(FILE* prmfd)
{
    int v;
    if(prm_read_int(prmfd,"BBL_NONLINEAR",&v) <= 0.)
	v = 1;
    return v;
}
int sinterface_gethindered_svel_patch(FILE* prmfd)
{
    int v;
    if(prm_read_int(prmfd,"HINDERED_SVEL_PATCH",&v) <= 0.)
	v = 0;
    return v;
}
int sinterface_gethindered_svel(FILE* prmfd)
{
    int v;
    if(prm_read_int(prmfd,"HINDERED_SVEL",&v) <= 0.)
	v = 0;
    return v;
}
double sinterface_reef_scale_depth(FILE* prmfd)
{
    double v;
    if(prm_read_double(prmfd,"REEF_SCALING_DEPTH",&v) <= 0.) {
       if(prm_read_double(prmfd,"REEF_SCALE_DEPTH",&v) <= 0.)
	  v = 0.;
    }
    return v;
}
int sinterface_getcalcripples(FILE* prmfd)
{
    int v;
    if(prm_read_int(prmfd,"CALC_RIPPLES",&v) <= 0.)
	v = 0;
    return v;
}
double sinterface_getcssscale(FILE* prmfd)
{
    double v;
    if(prm_read_double(prmfd,"CSS_SCALE",&v) <= 0.)
	v = 1.;
    return v;
}
double sinterface_getphysriph(FILE* prmfd)
{
    double v;
    if(prm_read_double(prmfd,"PHYSRIPH",&v) <= 0.)
	v = 0.01;
    return v;
}
double sinterface_getphysripl(FILE* prmfd)
{
    double v;
    if(prm_read_double(prmfd,"PHYSRIPL",&v) <= 0.)
	v = 0.5;
    return v;
}
double sinterface_getbioriph(FILE* prmfd)
{
    double v;
    if(prm_read_double(prmfd,"BIORIPH",&v) <= 0.)
	v = 0.;
    return v;
}
double sinterface_getbioripl(FILE* prmfd)
{
    double v;
    if(prm_read_double(prmfd,"BIORIPL",&v) <= 0.)
	v = 0.5;
    return v;
}
double sinterface_getbiodens(FILE* prmfd)
{
    double v;
    if(prm_read_double(prmfd,"BIODENS",&v)<=0.)
	v=0.;
    return v;
}
double sinterface_getmaxbiodepth(FILE* prmfd)
{
    double v;
    if(prm_read_double(prmfd,"MAXBIODEPTH",&v)<=0.)
	v=0.2;
    return v;
}
double sinterface_getbi_dissol_kz(FILE* prmfd)
{
    double v;
    if (prm_read_double(prmfd,"BI_DISSOL_KZ",&v)<=0.)
	v=0.;
    return v;
}
double sinterface_getbt_partic_kz(FILE* prmfd)
{
    double v;
    if(prm_read_double(prmfd,"BT_PARTIC_KZ",&v)<=0.)
	v=0.;
    return v;
}
double sinterface_getbi_dissol_kz_i(FILE* prmfd)
{
    double v;
    if (prm_read_double(prmfd,"BI_DISSOL_KZ_I",&v)<=0.)
	v=0.;
    return v;
}
double sinterface_getbt_partic_kz_i(FILE* prmfd)
{
    double v;
    if(prm_read_double(prmfd,"BT_PARTIC_KZ_I",&v)<=0.)
	v=0.;
    return v;
}
char sinterface_getbiosedprofile(FILE* prmfd)
{
	char v;
  char buf[MAXSTRLEN];
  if( prm_read_char(prmfd, "BIOSEDPROFILE", buf) )
		v = tolower(buf[0]);
  else
        v = 'p'; /* parabolic */
  return v;
}
double sinterface_getz0(FILE* prmfd)
{
    double v;
    if(prm_read_double(prmfd,"Z0_SKIN",&v) <= 0.)
	v = 0.000001;
    return v;
}
double sinterface_getquad_bfc(FILE* prmfd)
{
   double v;
    if(prm_read_double(prmfd,"QBFC",&v) <= 0.)
	v = 0.0002;
    return v;
}
int sinterface_getshipfile(FILE* prmfd, char *shipfile)
{
    int v=0;
    // char buf[MAXSTRLEN];
    if( prm_read_char(prmfd,"SHIP_FILE",shipfile) > 0) {
      v=1;
    }
    else
    {
	v=0;
    }
    return v;
}
double sinterface_erflux_scale(FILE* prmfd)
{
   double v;
    if(prm_read_double(prmfd,"ERFLUX_SCALE",&v) <= 0.)
        v = 1.0;
    return v;
}

/*************************************TRACER ATTRIBUTES (sed_init.c) */



FILE* si_getparamfile_tracer(FILE* fp)
{
	char buf [MAXSTRLEN];
	FILE* fpt;
	if(prm_read_char(fp, "NTRACERS", buf))
		return fp;
	else if (prm_read_char(fp, "TRACERFILE", buf))
	{
		fpt= fopen(buf,"r");
		if(fpt == NULL)
			hd_quit("Cannot open tracer file '%s'",buf);
	  return fpt;
	}
	hd_warn("Cannot find tracer file? ");
	return fp;
	/*return NULL;*/
  
}

FILE* si_getparamfile_sed(FILE *fp)
{
	char buf [MAXSTRLEN];
	FILE* fpt;
	if(prm_read_char(fp, "NSEDLAYERS", buf))
		return fp;
	else if (prm_read_char(fp, "SEDFILE", buf))
	{
		fpt= fopen(buf,"r");
		if(fpt == NULL)
			hd_quit("Cannot open sediment file '%s'",buf);
	  return fpt;
	}
	hd_quit("Cannot find sediment file? ");
	return NULL;
}


double sinterface_get_psize(void* model, char *name)
{
    geometry_t* window = (geometry_t*) model;
    win_priv_t *wincon=window->wincon;
    tracer_info_t *tr = i_get_tracer(model, name);
    trinfo_priv_sed_t *data = tr->private_data[TR_PRV_DATA_SED];
    
    double v=0;
    if (data)
      v = data->psize;
    return v;
}

double sinterface_get_b_dens(void* model, char *name)
{
    geometry_t* window = (geometry_t*) model;
    win_priv_t *wincon=window->wincon;
    tracer_info_t *tr = i_get_tracer(model, name);
    trinfo_priv_sed_t *data = tr->private_data[TR_PRV_DATA_SED];
    
    double v=0;
    if (data)
      v = data->b_dens;
    return v;
}

double sinterface_get_i_conc(void* model, char *name)
{
    geometry_t* window = (geometry_t*) model;
    win_priv_t *wincon=window->wincon;
    tracer_info_t *tr = i_get_tracer(model, name);
    trinfo_priv_sed_t *data = tr->private_data[TR_PRV_DATA_SED];
    
    double v=0;
    if (data)
      v = data->i_conc;
    return v;
}


double sinterface_get_css_erosion(void* model, char *name)
{
    geometry_t* window = (geometry_t*) model;
    win_priv_t *wincon=window->wincon;
    tracer_info_t *tr = i_get_tracer(model, name);
    trinfo_priv_sed_t *data = tr->private_data[TR_PRV_DATA_SED];

    double v=-1;
    if (data)
      v = data->css_erosion;
    return v;
}

double sinterface_get_css_deposition(void* model, char *name)
{
    geometry_t* window = (geometry_t*) model;
    win_priv_t *wincon=window->wincon;
    tracer_info_t *tr = i_get_tracer(model, name);
    trinfo_priv_sed_t *data = tr->private_data[TR_PRV_DATA_SED];

    double v=-1;
    if (data)
      v = data->css_deposition;
    return v;
}


double sinterface_get_svel(void* model, char *name)
{
    geometry_t* window = (geometry_t*) model;
    win_priv_t *wincon=window->wincon;
    tracer_info_t *tr = i_get_tracer(model, name);
    trinfo_priv_sed_t *data = tr->private_data[TR_PRV_DATA_SED];
    
    double v=0;
    if (data)
      v = data->svel;
    return v;
}

void sinterface_get_svel_name(void* model, char *name, char *sname)
{
    geometry_t* window = (geometry_t*) model;
    win_priv_t *wincon=window->wincon;
    tracer_info_t *tr = i_get_tracer(model, name);
    trinfo_priv_sed_t *data = tr->private_data[TR_PRV_DATA_SED];
    
    if (data)
      strcpy(sname, data->svel_name);
}

int sinterface_get_cohesive(void* model, char *name)
{
    geometry_t* window = (geometry_t*) model;
    win_priv_t *wincon=window->wincon;
    tracer_info_t *tr = i_get_tracer(model, name);
    trinfo_priv_sed_t *data = tr->private_data[TR_PRV_DATA_SED];
    
    int v=1;
    if (data)
      v = data->cohesive;
    return v;
}

int sinterface_get_resuspend(void* model, char *name)
{
    geometry_t* window = (geometry_t*) model;
    win_priv_t *wincon=window->wincon;
    tracer_info_t *tr = i_get_tracer(model, name);
    trinfo_priv_sed_t *data = tr->private_data[TR_PRV_DATA_SED];
    
    int v=1;
    if (data)
      v = data->resuspend;
    return v;
}

int sinterface_get_deposit(void* model, char *name)
{
    geometry_t* window = (geometry_t*) model;
    win_priv_t *wincon=window->wincon;
    tracer_info_t *tr = i_get_tracer(model, name);
    trinfo_priv_sed_t *data = tr->private_data[TR_PRV_DATA_SED];
    
    int v=1;
    if (data)
      v = data->deposit;
    return v;
}

int sinterface_get_floc(void* model, char *name)
{
    geometry_t* window = (geometry_t*) model;
    win_priv_t *wincon=window->wincon;
    tracer_info_t *tr = i_get_tracer(model, name);
    trinfo_priv_sed_t *data = tr->private_data[TR_PRV_DATA_SED];
    
    int v=0;
    if (data)
      v = data->floc;
    return v;
}

int sinterface_get_calcvol(void* model, char *name)
{
    geometry_t* window = (geometry_t*) model;
    win_priv_t *wincon=window->wincon;
    tracer_info_t *tr = i_get_tracer(model, name);
    trinfo_priv_sed_t *data = tr->private_data[TR_PRV_DATA_SED];
    int v=0;

    if (data)
      v = data->calcvol;
    return v;
}

int sinterface_get_adsorb(void* model, char *name)
{
    geometry_t* window = (geometry_t*) model;
    win_priv_t *wincon=window->wincon;
    tracer_info_t *tr = i_get_tracer(model, name);
    trinfo_priv_sed_t *data = tr->private_data[TR_PRV_DATA_SED];
    
    int v=0;
    if (data != NULL)
      v = data->adsorb;
    return v;
}

void sinterface_get_carriername(void* model, char *name,char *carriername )
{
    geometry_t* window = (geometry_t*) model;
    win_priv_t *wincon=window->wincon;
    tracer_info_t *tr = i_get_tracer(model, name);
    trinfo_priv_sed_t *data = tr->private_data[TR_PRV_DATA_SED];

    if (data)
      strcpy(carriername, data->carriername);
}

void sinterface_get_dissolvedname(void* model, char *name, char *dissolvedname)
{
    geometry_t* window = (geometry_t*) model;
    win_priv_t *wincon=window->wincon;
    tracer_info_t *tr = i_get_tracer(model, name);
    trinfo_priv_sed_t *data = tr->private_data[TR_PRV_DATA_SED];

    if (data)
      strcpy(dissolvedname, data->dissolvedname);
}

double sinterface_get_adsorb_kd(void* model, char *name)
{
    geometry_t* window = (geometry_t*) model;
    win_priv_t *wincon=window->wincon;
    tracer_info_t *tr = i_get_tracer(model, name);
    trinfo_priv_sed_t *data = tr->private_data[TR_PRV_DATA_SED];
    
    int v=0;
    if (data)
      v = data->adsorbkd;
    return v;
}

double sinterface_get_adsorb_rate(void* model, char *name)
{
    geometry_t* window = (geometry_t*) model;
    win_priv_t *wincon=window->wincon;
    tracer_info_t *tr = i_get_tracer(model, name);
    trinfo_priv_sed_t *data = tr->private_data[TR_PRV_DATA_SED];
    
    int v=0;
    if (data)
      v = data->adsorbrate;
    return v;
}

// 2021
// Read decay_days (e-folding time in days) from prm file
// Default is 0
double sinterface_get_decay_days(FILE* prmfd, void* model, int tr_number)
{
    char keybuf[MAXSTRLEN], buf[MAXSTRLEN];
    double v;

    strcpy(keybuf,"TRACER"); 
    strcat(keybuf, _itoa(tr_number, buf, 10) ); 
    strcat(keybuf, ".decay_days");

    if(prm_read_double(prmfd,keybuf,&v) <= 0.)
        v = 0;
    return v;

}



/************************************** CONCENTRATIONS */


/*hd2sed

Function:
Reads concentrations of tracers in wc

Input:
hmodel is pointer to hd model structure;
c is  number of current column;
n is tracer number as specified in mecosed;
m is tracer number as specified in hd;;

Output:
returns void;
updates ptr_wc[n][k] - a 2d array of pointers to hd tracers.

Note:
k index changes from botk_wc to topk_wc,
where botk_wc and topk_wc are numbers of the bottom and top cells,
respectively, and
botk_wc < topk_wc <= (nwclayers-1);
ptr_wc[n][botk_wc] is address of tracer in the bottom cell;
ptr_wc[n][topk_wc] is address of tracer in the top cell;

If in the hd model cell numbering increases downward
(topk_wc_hd < botk_wc_hd), make sure that in ptr_wc the cell
numbering increases upward, starting from topk_wc_hd in the bottom cell
and to  botk_wc_hd in the top cell.
*/

void si_gettracer_wc(void* hmodel, int c,
			     int n, int m, double ***ptr_wc)
{
    geometry_t* window = (geometry_t*) hmodel;
    window_t *windat = window->windat;

    int c2=window->m2d[c];
    int cc=window->c2cc[c2];
    int botk_wc= window->s2k[window->bot_t[cc]];
    int topk_wc= window->s2k[window->nsur_t[cc]];
    int c3=window->nsur_t[cc];
    int k;
    for(k=topk_wc;k>=botk_wc;k--) {
	ptr_wc[n][k] = &windat->tr_wc[m][c3];
	c3 = window->zm1[c3];
    }
}


/*hd2sed

Function:
Reads tracer concentrations for given sed column

Input:
hmodel is pointer to hd model structure;
c is  number of current column;
n is tracer number as specified in mecosed;
m is tracer number as specified in hd;

Output:
returns void;
updates ptr_sed[n][k] - a 2d array of pointers to hd tracers.

Note:
k index changes from botk_sed to topk_sed,
where botk_sed and topk_sed are numbers of the bottom and top cells,
respectively, and
botk_sed < topk_sed = (numbsedlayers-1);
ptr_sed[n][botk_sed] is address of tracer in the bottom cell;
ptr_sed[n][topk_sed] is address of tracer in the top cell;
*/

void si_gettracer_sed(void* hmodel, int c,
			      int n, int m, double ***ptr_sed)
{
    geometry_t* window = (geometry_t*) hmodel;
    window_t *windat = window->windat;
    int c2=window->m2d[c];
    int bot_k= 0;
    int top_k = window->sednz-1;
    int k;
    for(k=bot_k;k<=top_k;k++)
	ptr_sed[n][k] = &windat->tr_sed[m][k][c2];
}

/***************************************** DIAGNOSTIC OUTPUT */

/* Write diagnostic variables */
/*sed_init
Function:
  get map to 2d diagnostic tracers;
Input:
  hmodel is a pointer to hd model structure;
Output:
  returns void;
  updates n_hripple - hd number of 2d tracer whith name  "hripple";
  updates n_lripple ...
Note:
  Nothing happens if the diagnostic tracer is not specified in the prm file.
*/

void sinterface_getmap2diagtracer_2d(void* hmodel,
int *n_hripple, int *n_lripple, int *n_ustrcw_skin,
int *n_depth_sedi, int *n_dzactive, int *n_erdepflux_total,
int *n_erdepflux_total_ac,int *n_erdepflux_oxygen, int *n_erdepflux_oxygen_ac)
{
    geometry_t* window = (geometry_t*) hmodel;
    window_t *windat = window->windat;
    win_priv_t *wincon=window->wincon;
    int n;
    char name[MAXSTRLEN];

  n=-1;
  *n_hripple= n;
  *n_lripple= n;
  *n_ustrcw_skin= n;
  *n_depth_sedi= n;
  *n_dzactive=n;
  *n_erdepflux_total= n;
  *n_erdepflux_total_ac = n;
  *n_erdepflux_oxygen= n;
  *n_erdepflux_oxygen_ac = n;

    for(n=0; n < windat->ntrS; n++) {
	strcpy(name, wincon->trinfo_2d[n].name);
	if(strcmp(name, "hripple")==0)
	    *n_hripple= n;
	else if(strcmp(name, "lripple")==0)
	    *n_lripple= n;
	else if(strcmp(name, "ustrcw_skin")==0)
	    *n_ustrcw_skin= n;
	else if(strcmp(name, "depth_sedi")==0)
	    *n_depth_sedi= n;
	else if(strcmp(name, "dzactive")==0)
	    *n_dzactive= n;
	else if(strcmp(name, "erdepflux_total")==0)
	    *n_erdepflux_total= n;
	else if(strcmp(name, "erdepflux_total_ac")==0)
	    *n_erdepflux_total_ac = n;
       else if(strcmp(name, "erdepflux_oxygen")==0)
            *n_erdepflux_oxygen = n;
       else if(strcmp(name, "erdepflux_oxygen_ac")==0)
            *n_erdepflux_oxygen_ac = n;
    }
}

/*sed2hd
Function:
  write 2d diagnostic tracers to hd;
Input:
  hmodel is a pointer to hd model structure;
  c is number of current column;
  hripple is calculated ripple hight in column c.
  n_hripple is hd number of a 2d array with name "hripple"
  lripple is ...
  n_lripple is ...
  ...
Output:
  returns void;
  updates 2d hd arrays with numbers n_hripple, n_lripple ...;
Note:
  Nothing happens if the diagnostic tracer is not specified in the prm file.
*/
void sinterface_putdiagtracer_2d(void* hmodel, int c,
double hripple,int n_hripple,  double lripple, int n_lripple,
double ustrcw_skin, int n_ustrcw_skin,
double depth_sedi, int n_depth_sedi, double dzactive,  int n_dzactive,
double erdepflux_total, int n_erdepflux_total, int n_erdepflux_total_ac,
double erdepflux_oxygen,int n_erdepflux_oxygen, int n_erdepflux_oxygen_ac)

{
    geometry_t* window = (geometry_t*) hmodel;
    window_t *windat = window->windat;
    int c2=window->m2d[c];

	if(n_hripple>-1)
	    windat->tr_wcS[n_hripple][c2]= hripple;
	if(n_lripple>-1)
	    windat->tr_wcS[n_lripple][c2]= lripple;
	if(n_ustrcw_skin>-1)
	    windat->tr_wcS[n_ustrcw_skin][c2]= ustrcw_skin;
	if(n_depth_sedi>-1)
	    windat->tr_wcS[n_depth_sedi][c2]= depth_sedi;
	if(n_dzactive>-1)
	    windat->tr_wcS[n_dzactive][c2]= dzactive;
	if(n_erdepflux_total>-1)
	    windat->tr_wcS[n_erdepflux_total][c2] = erdepflux_total;
	if(n_erdepflux_total_ac>-1)
	    windat->tr_wcS[n_erdepflux_total_ac][c2] +=
		windat->dt*erdepflux_total;
        if(n_erdepflux_oxygen>-1)
            windat->tr_wcS[n_erdepflux_oxygen][c2] = erdepflux_oxygen;
        if(n_erdepflux_oxygen_ac>-1)
            windat->tr_wcS[n_erdepflux_oxygen_ac][c2] +=
                windat->dt*erdepflux_oxygen;

}

/************************************************/

void sed_set_grid(geometry_t *geom, geometry_t **window)
{
  int n, cc, c, k;

  return;
  if (geom->sednz) {
    for (n = 1; n <= geom->nwindows; n++) {
      win_priv_t *wincon = window[n]->wincon;
      sediment_t *sediment = wincon->sediment;
      sed_column_t *sm = sediment->sm;

      for (cc = 1; cc <= window[n]->b2_t; c++) {
        c = window[n]->w2_t[cc];
        for (k = 0; k < window[n]->sednz; k++) {
          int kk;               /* change order of numbering to provide
                                   input to ecology */
          kk = window[n]->sednz - 1 - k;
          geom->gridz_sed[k][c] = sm->gridz_sed[kk];
        }
        for (k = 0; k < (window[n]->sednz - 1); k++) {
          int kk;               /* change order of numbering to provide
                                   input to ecology */
          kk = window[n]->sednz - 2 - k;
          geom->cellz_sed[k][c] = sm->cellz_sed[kk];
        }
      }
    }
  }
}


/******************************* Custom routines */

/*
Function:
Custom settling velocities

Input:
hmodel is pointer to hd model
c is hd-index of current column
n is tracer number as specified in mecosed
name is tracer name
topk_wc is number of the top wc cell
botk_wc is number of the bottom wc cell

Output:
svel_wc[k] is 1d array of settling velocities

Note: k-index changes from botk_wc to topk_wc
*/


void sinterface_getsvel_custom(void* hmodel, int c, int n,
		     char *name, double *svel_wc, int topk_wc, int botk_wc)
{
    geometry_t* window = (geometry_t*) hmodel;
    window_t *windat = window->windat;

    int k;
    double max_swim=4/3600.;
    win_priv_t *wincon=window->wincon;
    int c2=window->m2d[c];
    int cc=window->c2cc[c2];
    int c3=window->nsur_t[cc];
    double botz =  window->botz[c2]*wincon->Hn1[c2];
    double topz = windat->eta[c2];
    double depth = topz-botz;
    double dayfrac, zlevel;

    double depth2 = 4;
    if (depth<depth2) depth2=depth;

    if(depth>25.) depth=25.;
    if (strcmp(name,"PhyD_N") == 0 || strcmp(name,"PhyD_C") == 0 )
    {
      dayfrac = fmod(windat->t, 86400)/86400.;
      if (dayfrac >0.25 && dayfrac < 0.75) // day time
	max_swim = fabs(max_swim);
      else //night time
	max_swim = -fabs(max_swim);
      
      zlevel=0.;
      for(k=topk_wc;k>=botk_wc;k--) {
	zlevel = zlevel + 0.5*wincon->dz[c3]*wincon->Hn1[c2];
	
	/*  to activate Dinoflagellate Diel vertical migration code items 1 or 2 below must be activated and
	    item 3 commented out.  The whole code must then be re-compiled into a unique version of shoc.  */
	
	/* 1. Dinoflagellate DVM between 4-25m *  /
	   if(zlevel<depth)
	   svel_wc[k] = -0.1 * fabs(max_swim);
	   else if((zlevel>depth2) && (zlevel<depth))
	   svel_wc[k] = max_swim * (1.-fabs((zlevel-depth2)-0.5*(depth-depth2))/(0.5*(depth-depth2))) ;
	   else
	   svel_wc[k] = 0.1 * fabs(max_swim);
	   
	   if(!log_DVM)
	   {
	   emstag(LINFO,"hd:sediments:sinterface_getsvel_custom","Dinoflagellate vertical migration ON -DVM between 4-25m");
	   log_DVM = 1;
	   }
	*/
	
	/* 2. Dinoflagellate DVM to 0-25m   * /
	   if(zlevel<depth)
	   svel_wc[k] = max_swim * (1.-fabs(zlevel-0.5*depth)/(0.5*depth)) ;
	   else
	   svel_wc[k] = 0.05 * fabs(max_swim);
	   
	   if(!log_DVM)
	   {
	   emstag(LINFO,"hd:sediments:sinterface_getsvel_custom","Dinoflagellate vertical migration ON - DVM to 0-25m");
	   log_DVM = 1;
	   }
	*/
	
	/* 3. No DVM */
	svel_wc[k] = 0.0;
	if(!log_DVM)
	  {
	    emstag(LINFO,"hd:sediments:sinterface_getsvel_custom","Dinoflagellate vertical migration OFF");
	    log_DVM = 1;
	  }
	/* END vertical migration */
	c3 = window->zm1[c3];
      }
      
      svel_wc[botk_wc] =0.;
      svel_wc[topk_wc+1] =0.;
    } 
/* Trichodesmium variable buoyancy */
    else if (strcmp(name,"Tricho_N") == 0 || strcmp(name,"Tricho_NR") == 0 || strcmp(name,"Tricho_PR") == 0 || strcmp(name, "Tricho_I") == 0 || strcmp(name, "Tricho_Chl") == 0) {
      int Tricho_sv_i=-1;
      Tricho_sv_i = gi_getmap2hdwctracer(hmodel, n, "Tricho_sv");

      if (Tricho_sv_i < 0) {
        hd_warn("Diagnostic tracer Tricho_sv not found (required to calculate variable settling velocity for Trichodesmium)\n");
        for(k=topk_wc+1;k>=botk_wc;k--)
           svel_wc[k] = 0.0;
      } else {
        c3=window->nsur_t[cc];
  
        for(k=topk_wc;k>=botk_wc;k--) {
           svel_wc[k] = windat->tr_wc[Tricho_sv_i][c3];
	   c3 = window->zm1[c3];
        }
        svel_wc[botk_wc] =0.;
        svel_wc[topk_wc+1] =0.;
      }
    }
}


/*-------------------------------------------------------------------*/
/* Sets the sediment tracers standard default attributes    
     ./hd/tracers/tracer_info.c       */
/*-------------------------------------------------------------------*/
void sed_set_tracer_defaults(tracer_info_t *tracer, char *trname, char *defname)
{
  int i;
  struct {
    char *name;
    char *description;
    void (*init) (tracer_info_t *tracer, char *trname);
  } def_list[] = {
    {"standard","Standard sediment values",  sed_defaults_std},
    {"estuary", "Estuarine sediment values", sed_defaults_est},
    {"shelf", "Shelf sediment values", sed_defaults_shf},
    {"basic", "Basic sediment values", sed_defaults_bsc},
    {NULL, NULL, NULL}
  };

  void (*init) (tracer_info_t *tracer, char *trname)= NULL;

  if (defname != NULL) {
    /* Search for the sediments scheme name in the above table */
    for (i = 0; (init == NULL) && def_list[i].name; ++i) {
      if (strcasecmp(defname, def_list[i].name) == 0) {
        init = def_list[i].init;
      }
    }
  }
  if (init != NULL) {
    init(tracer, trname);
  }
}

/* END sed_set_tracer_defaults()                                     */
/*-------------------------------------------------------------------*/


#define T1 WATER|SEDIM|SEDIMENT|PROGNOSTIC
#define T2 WATER|SEDIM|SEDIMENT|DIAGNOSTIC
#define T3 INTER|SEDIMENT|DIAGNOSTIC

/*-------------------------------------------------------------------*/
/* Sets the sediment tracers standard default attributes             */
/*                                                                   */
/* Alternatively, may use the expanded format here:                  */
/*   tracer[i].name = "Gravel";                                      */
/*   tracer[i].unit = "kg m-3";                                      */
/*   tracer[i].type = 1;                                             */
/*   tracer[i].fill = 200.0                                          */
/*   etc ...                                                         */
/*-------------------------------------------------------------------*/
void sed_defaults_std(tracer_info_t *tracer, char *trname)
{

 sed_def_t sed_def[] = {
    /* name         unit     type fill_sed   psize    b_dens    i_conc  css_ero css_depo svel dg ad df coh cv fl re de  obc*/
    {"Gravel",      "kg m-3", T1, 200.0,   0.0005,  2.65e+03, 1.2e+03,  -1., -1.,       -0.7,   0, 0, 1, 0, 1, 0, 0, 1, TRCONC|FILEIN},
    {"Sand",        "kg m-3", T1, 400.0,   0.00025, 2.65e+03, 1.2e+03,  -1., -1.,       -2e-3,  0, 0, 1, 0, 1, 0, 0, 1, TRCONC|FILEIN},
    {"Mud",         "kg m-3", T1, 600.0,   0.00025, 2.65e+03, 1.2e+03,  -1., -1.,       -2.e-4, 0, 1, 1, 1, 1, 1, 1, 1, TRCONC|FILEIN},
    {"Gravel-mineral",      "kg m-3", T1, 200.0,   0.0005,  2.65e+03, 1.2e+03,  -1., -1.,       -0.7,   0, 0, 1, 0, 1, 0, 0, 1, TRCONC|FILEIN},
    {"Sand-mineral",        "kg m-3", T1, 400.0,   0.00025, 2.65e+03, 1.2e+03,  -1., -1.,       -2e-3,  0, 0, 1, 0, 1, 0, 0, 1, TRCONC|FILEIN},
    {"Mud-mineral",         "kg m-3", T1, 600.0,   0.00025, 2.65e+03, 1.2e+03,  -1., -1.,       -2.e-4, 0, 1, 1, 1, 1, 1, 1, 1, TRCONC|FILEIN},
    {"Gravel-carbonate",      "kg m-3", T1, 200.0,   0.0005,  2.65e+03, 1.2e+03,  -1., -1.,       -0.7,   0, 0, 1, 0, 1, 0, 0, 1, TRCONC|FILEIN},
    {"Sand-carbonate",        "kg m-3", T1, 400.0,   0.00025, 2.65e+03, 1.2e+03,  -1., -1.,       -2e-3,  0, 0, 1, 0, 1, 0, 0, 1, TRCONC|FILEIN},
    {"Mud-carbonate",         "kg m-3", T1, 600.0,   0.00025, 2.65e+03, 1.2e+03,  -1., -1.,       -2.e-4, 0, 1, 1, 1, 1, 1, 1, 1, TRCONC|FILEIN},
    {"FineSed",     "kg m-3", T1, 0.0,     0.00025, 2.65e+03, 1.2e+03,  -1., -1.,       -2.e-4, 0, 1, 1, 1, 1, 1, 1, 1, TRCONC|FILEIN},
    {"Dust",        "kg m-3", T1, 0.0,     0.00025, 2.65e+03, 1.2e+03,  -1., -1.,       -1.1574e-05, 0, 1, 1, 1, 1, 0, 1, 1, TRCONC|FILEIN},
    {"tss",         "kg m-3", T2, 0.0,     0.0,     0.0,      0.0,      -1., -1.,       0.0,    1, 0, 0, 0, 0, 0, 0, 0, NOTHIN},
    {"porosity",    "-",      T2, 0.5,     0.0,     0.0,      0.0,      -1., -1.,       0.0,    2, 0, 0, 0, 0, 0, 0, 0, NOTHIN},
    {"ustrcw_skin", "m s-1",  T3, 0.0,     0.0,     0.0,      0.0,      -1., -1.,       0.0,    0, 0, 0, 0, 0, 0, 0, 0, NOTHIN},
    {"depth_sedi",  "m    ",  T3, 0.0,     0.0,     0.0,      0.0,      -1., -1.,       0.0,    1, 0, 0, 0, 0, 0, 0, 0, NOTHIN},
    {"erdepflux_total_ac",   "kg m-2",  T3,  0.0,   0.0,      0.0, 0.0, -1., -1.,       0.0,    1, 0, 0, 0, 0, 0, 0, 0, NOTHIN},
    {"NULL",        "NULL",   1, 0.0,      0.0,     0.0,      0.0,      -1., -1.,       0.0,    0, 0, 0, 0, 0, 0, 0, 0, NOTHIN},
  };



  /* Init attributes */
  init_tracer_atts_sed(tracer, trname, sed_def, "standard");
}

/* END sed_defaults_std()                                            */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Sets the sediment tracers estuarine default attributes            */
/*-------------------------------------------------------------------*/
static void sed_defaults_est(tracer_info_t *tracer, char *trname)
{

  sed_def_t sed_def[] = {
    /* name         unit     type fill_sed   psize    b_dens    i_conc css_er css_dep    svel dg ad df coh cv fl re de  obc*/
    {"Gravel",      "kg m-3", T1, 200.0,   0.0005,  2.65e+03, 1.2e+03,  -1.,  -1.,      -0.7,   0, 0, 1, 0, 1, 0, 0, 1, TRCONC|FILEIN},
    {"Sand",        "kg m-3", T1, 400.0,   0.00025, 2.65e+03, 1.2e+03,   -1.,  -1.,     -2e-3,  0, 0, 1, 0, 1, 0, 0, 1, TRCONC|FILEIN},
    {"Mud",         "kg m-3", T1, 600.0,   0.00025, 2.65e+03, 1.2e+03,   -1.,  -1.,     -2.e-4, 0, 1, 1, 1, 1, 1, 1, 1, TRCONC|FILEIN},
    {"CarbSand",    "kg m-3", T1, 0.0,     0.00025, 2.65e+03, 1.2e+03,   -1.,  -1.,     -2.e-4, 0, 1, 1, 1, 1, 1, 1, 1, TRCONC|FILEIN},
    {"FineSed",     "kg m-3", T1, 0.0,     0.00025, 2.65e+03, 1.2e+03,   -1.,  -1.,     -2.e-4, 0, 1, 1, 1, 1, 1, 1, 1, TRCONC|FILEIN},
    {"Dust",        "kg m-3", T1, 0.0,     0.00025, 2.65e+03, 1.2e+03,   -1.,  -1.,     -1.1574e-05, 0, 1, 1, 1, 1, 0, 1, 1, TRCONC|FILEIN},
    {"porosity",    "-",      T2, 0.5,     0.0,     0.0,      0.0,       -1.,  -1.,     0.0,    1, 0, 0, 0, 0, 0, 0, 0, NOTHIN},
    {"tss",         "kg m-3", T2, 0.0,     0.0,     0.0,      0.0,       -1.,  -1.,     0.0,    1, 0, 0, 0, 0, 0, 0, 0, NOTHIN},
    {"ustrcw_skin", "m s-1",  T3, 0.0,     0.0,     0.0,      0.0,       -1.,  -1.,     0.0,    1, 0, 0, 0, 0, 0, 0, 0, NOTHIN},
    {"depth_sedi",   "m   ",  T3, 0.0,     0.0,     0.0,      0.0,       -1.,  -1.,     0.0,    1, 0, 0, 0, 0, 0, 0, 0, NOTHIN},
    {"erdepflux_total_ac",   "kg m-2",  T3, 0.0,    0.0,      0.0, 0.0,  -1.,  -1.,     0.0,    1, 0, 0, 0, 0, 0, 0, 0, NOTHIN},
    {"NULL",        "NULL",   1, 0.0,      0.0,     0.0,      0.0,       -1.,  -1.,     0.0,    0, 0, 0, 0, 0, 0, 0, 0, NOTHIN},
  };



  /* Init attributes */
  init_tracer_atts_sed(tracer, trname, sed_def, "estuarine");
}

/* END sed_defaults_est()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets the sediment tracers shelf default attributes            */
/*-------------------------------------------------------------------*/
static void sed_defaults_shf(tracer_info_t *tracer, char *trname)
{
 sed_def_t sed_def[] = {
    /* name         unit     type fill_sed   psize    b_dens    i_conc  css_er css_dep  svel dg ad df coh cv fl re de  obc*/
    {"Gravel",      "kg m-3", T1, 200.0,   0.0005,  2.65e+03, 1.2e+03,  -1.,  -1.,      -0.7,   0, 0, 1, 0, 1, 0, 0, 1, TRCONC|FILEIN},
    {"Sand",        "kg m-3", T1, 400.0,   0.00025, 2.65e+03, 1.2e+03,  -1.,  -1.,      -2e-3,  0, 0, 1, 0, 1, 0, 0, 1, TRCONC|FILEIN},
    {"Mud",         "kg m-3", T1, 600.0,   0.00025, 2.65e+03, 1.2e+03,  -1.,  -1.,      -2.e-4, 0, 1, 1, 1, 1, 1, 1, 1, TRCONC|FILEIN},
    {"FineSed",     "kg m-3", T1, 0.0,     0.00025, 2.65e+03, 1.2e+03,  -1.,  -1.,      -2.e-4, 0, 1, 1, 1, 1, 1, 1, 1, TRCONC|FILEIN},
    {"Dust",        "kg m-3", T1, 0.0,     0.00025, 2.65e+03, 1.2e+03,  -1.,  -1.,      -1.1574e-05, 0, 1, 1, 1, 1, 0, 1, 1, TRCONC|FILEIN},
    {"tss",         "kg m-3", T2, 0.0,     0.0,     0.0,      0.0,      -1.,  -1.,      0.0,    1, 0, 0, 0, 0, 0, 0, 0, NOGRAD},
    {"porosity",    "-",      T2, 0.5,     0.0,     0.0,      0.0,      -1.,  -1.,      0.0,    1, 0, 0, 0, 0, 0, 0, 0, NOTHIN},
    {"ustrcw_skin", "m s-1",  T3, 0.0,     0.0,     0.0,      0.0,      -1.,  -1.,      0.0,    1, 0, 0, 0, 0, 0, 0, 0, NOTHIN},
    {"depth_sedi",   "m   ",  T3, 0.0,     0.0,     0.0,      0.0,      -1.,  -1.,      0.0,    1, 0, 0, 0, 0, 0, 0, 0, NOTHIN},
    {"erdepflux_total_ac",   "kg m-2",  T3, 0.0,    0.0,      0.0, 0.0, -1.,  -1.,      0.0,    1, 0, 0, 0, 0, 0, 0, 0, NOTHIN},
    {"NULL",        "NULL",   1, 0.0,      0.0,     0.0,      0.0,      -1.,  -1.,      0.0,    0, 0, 0, 0, 0, 0, 0, 0, NOTHIN},
  };

  /* Init attributes */
  init_tracer_atts_sed(tracer, trname, sed_def, "shelf");
}

/* END sed_defaults_est()                                            */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Sets the sediment tracers basic default attributes            */
/*-------------------------------------------------------------------*/
static void sed_defaults_bsc(tracer_info_t *tracer, char *trname)
{

  sed_def_t sed_def[] = {
    /* name         unit     type fill_sed   psize    b_dens    i_conc  css_er css_dep svel    dg ad df coh cv fl re de  obc*/
    {"Gravel",      "kg m-3", T1, 200.0,   0.0005,  2.65e+03, 1.2e+03,  -1.,  -1.,      -0.7,   0, 0, 1, 0, 1, 0, 0, 1, TRCONC|FILEIN},
    {"Sand",        "kg m-3", T1, 400.0,   0.00025, 2.65e+03, 1.2e+03,  -1.,  -1.,      -2e-3,  0, 0, 1, 0, 1, 0, 0, 1, TRCONC|FILEIN},
    {"Mud",         "kg m-3", T1, 600.0,   0.00025, 2.65e+03, 1.2e+03,  -1.,  -1.,      -2.e-4, 0, 1, 1, 1, 1, 1, 1, 1, TRCONC|FILEIN},
    {"FineSed",     "kg m-3", T1, 0.0,     0.00025, 2.65e+03, 1.2e+03,  -1.,  -1.,      -2.e-4, 0, 1, 1, 1, 1, 1, 1, 1, TRCONC|FILEIN},
    {"Dust",        "kg m-3", T1, 0.0,     0.00025, 2.65e+03, 1.2e+03,  -1.,  -1.,      -1.1574e-05, 0, 1, 1, 1, 1, 0, 1, 1, TRCONC|FILEIN},
    {"tss",         "kg m-3", T2, 0.0,     0.0,     0.0,      0.0,      -1.,  -1.,      0.0,    1, 0, 0, 0, 0, 0, 0, 0, NOGRAD},
    {"porosity",    "-",      T2, 0.5,     0.0,     0.0,      0.0,      -1.,  -1.,      0.0,    1, 0, 0, 0, 0, 0, 0, 0, NOTHIN},
    {"ustrcw_skin", "m s-1",  T3, 0.0,     0.0,     0.0,      0.0,      -1.,  -1.,      0.0,    1, 0, 0, 0, 0, 0, 0, 0, NOTHIN},
    {"depth_sedi",   "m   ",  T3, 0.0,     0.0,     0.0,      0.0,      -1.,  -1.,      0.0,    1, 0, 0, 0, 0, 0, 0, 0, NOTHIN},
    {"erdepflux_total_ac",   "kg m-2",  T3, 0.0,    0.0,      0.0, 0.0, -1.,  -1.,      0.0,    1, 0, 0, 0, 0, 0, 0, 0, NOTHIN},
    {"NULL",        "NULL",   1, 0.0,      0.0,     0.0,      0.0,      -1.,  -1.,      0.0,    0, 0, 0, 0, 0, 0, 0, 0, NOTHIN},
  };


  /* Init attributes */
  init_tracer_atts_sed(tracer, trname, sed_def, "basic");
}

/* END sed_defaults_est()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets the tracer attributes for the given tracer, if found         */
/*-------------------------------------------------------------------*/
static void init_tracer_atts_sed(tracer_info_t *tracer, char *trname, 
				 sed_def_t *sed_def, char *sed_def_name)
{
  int n, i, found;
  
  /* Find the tracer number in the list */
  found = 0;
  for (n = 0; strcmp(sed_def[n].name, "NULL") != 0; ++n) {
    if (strcmp(sed_def[n].name, trname) == 0) {
      found = 1;
      break;
    }
  }

  /* Set the values */
  if (found) {
    /* Allocate sediment private data for this tracer */

    trinfo_priv_sed_t *data = get_private_data_sed(tracer);

    /* Whether this is standard, esturine etc ..*/
    strcpy(data->name, sed_def_name);

    /* Same short/long name */
    strcpy(tracer->name, trname);
    strcpy(tracer->long_name, tracer->name);


    /* Fill values */
    tracer->fill_value_wc = 0.0;
    tracer->valid_range_wc[0]  =  -1e35;
    tracer->valid_range_sed[0] =  -1e35;
    tracer->valid_range_wc[1]  =  1e35;
    tracer->valid_range_sed[1] =  1e35;
    /* For particulates */
    if (sed_def[n].psize > 0) {
      tracer->inwc   = 1;
      tracer->insed  = 1;
      tracer->dissol = 0;
      tracer->partic = 1;
      tracer->valid_range_wc[0]  =  0.0;
      tracer->valid_range_sed[0] =  0.0;
    }
    
    /* Set the standard values */
    strcpy(tracer->units, sed_def[n].units);
    tracer->fill_value_sed = sed_def[n].fill_sed;
    tracer->diagn   = sed_def[n].diagn;
    tracer->advect  = sed_def[n].advect;
    tracer->diffuse = sed_def[n].diffuse;
    tracer->m = -1;
    tracer->type = sed_def[n].type;


    /* Set private data */
    data->psize  = sed_def[n].psize;
    data->b_dens = sed_def[n].b_dens;
    data->i_conc = sed_def[n].i_conc;
    data->f_conc = sed_def[n].i_conc;  // by default

   //2019
    data->css_erosion = sed_def[n].css_erosion;
    data->css_deposition = sed_def[n].css_deposition;

    data->svel   = sed_def[n].svel;
    sprintf(data->svel_name, "%c", '\0');
    data->cohesive  = sed_def[n].cohesive;
    data->calcvol   = sed_def[n].calcvol;
    data->floc      = sed_def[n].floc;
    data->resuspend = sed_def[n].resuspend;
    data->deposit   = sed_def[n].deposit;
    data->obc       = sed_def[n].obc;
  }
}

/*-------------------------------------------------------------------*/
/* Initialised sediment private data                                 */
/*-------------------------------------------------------------------*/
static trinfo_priv_sed_t *get_private_data_sed(tracer_info_t *tr)
{
  trinfo_priv_sed_t *data = tr->private_data[TR_PRV_DATA_SED];

  /* Allocate memory, if needed */
  if (data == NULL) {
    data = (trinfo_priv_sed_t *)malloc(sizeof(trinfo_priv_sed_t));
    /*
     * Useful for debugging in case the wrong pointer happens 
     * to be passed in 
     */
    data->type = PTR_SED;

    /* Set data and copy function */
    tr->private_data[TR_PRV_DATA_SED]      = data;
    tr->private_data_copy[TR_PRV_DATA_SED] = private_data_copy_sed;

    sprintf(data->name, "%c", '\0');
    data->psize = 0.0;
    data->b_dens = 0.0;
    data->i_conc = 0.0;
    data->f_conc = 0.0;

  //2019
    data->css_erosion = -1.0;
    data->css_deposition = -1.0;

    data->svel = 0.0;
    sprintf(data->svel_name, "%c", '\0');
    data->cohesive = 1;
    data->calcvol = 0;
    data->floc = 0;
    data->resuspend = 1;
    data->adsorb = 0;
    data->adsorbkd = 0.0;
    data->adsorbrate = 0.0;
    data->deposit = 1;
    data->obc = NOGRAD;
    sprintf(data->carriername, "%c", '\0');
    sprintf(data->dissolvedname, "%c", '\0');
  }

  return(data);
}

/* END get_private_data_sed()                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Checks if a trcer is a valid sediment tracer         
     ./hd/tracers/tracer_info.c          */
/*-------------------------------------------------------------------*/
int is_sed_var(char *trname)
{
  int i;

  for (i = 0; i < NUM_SED_VARS; i++)
    if (strcmp(trname, SEDCLASS[i]) == 0)
      return(1);

  for (i = 0; i < NUM_SED_VARS_3D; i++)
    if (strcmp(trname, SEDNAME3D[i]) == 0)
      return(1);

  for (i = 0; i < NUM_SED_VARS_2D; i++)
    if (strcmp(trname, SEDNAME2D[i]) == 0)
      return(1);

  return (0);
}

/* END is_sed_tracer()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Set a particular private_data attribute                           */
/*-------------------------------------------------------------------*/
void sed_set_tr_att(tracer_info_t *tr, char *key, void *value)
{
  trinfo_priv_sed_t *data = get_private_data_sed(tr);

  /* Switch on key */
  if (strcmp(key, "psize") == 0) 
    {
      data->psize = *(double *)value;
    } 
  else if (strcmp(key, "b_dens") == 0) 
    {
      data->b_dens = *(double *)value;
    }
  else if (strcmp(key, "i_conc") == 0) 
    {
      data->i_conc = *(double *)value;
      data->f_conc = data->i_conc;
    } 
  else if (strcmp(key, "f_conc") == 0) 
    {
      data->f_conc = *(double *)value;
    }
 else if (strcmp(key, "css_erosion") == 0)
    { //2019
      data->css_erosion = *(double *)value;
    }
 else if (strcmp(key, "css_deposition") == 0)
    {
      data->css_deposition = *(double *)value;
    }

  else if (strcmp(key, "svel") == 0) 
    {
      data->svel = *(double *)value;
    } 
  else if (strcmp(key, "svel_name") == 0) 
    {
      sprintf(data->svel_name, "%s", (char*)value);
    } 
  else if (strcmp(key, "cohesive") == 0) 
    {
      data->cohesive = *(int *)value;
    } 
  else if (strcmp(key, "calcvol") == 0) 
    {
      data->calcvol = *(int *)value;
    } 
  else if (strcmp(key, "floc") == 0) 
    {
      data->floc = *(int *)value;
    } 
  else if (strcmp(key, "resuspend") == 0) 
    {
      data->resuspend = *(int *)value;
    } 
  else if (strcmp(key, "deposit") == 0) 
    {
      data->deposit = *(int *)value;
    } 
  else if (strcmp(key, "adsorb") == 0) 
    {
      data->adsorb = *(int *)value;
    } 
  else if (strcmp(key, "adsorbkd") == 0) 
    {
      data->adsorbkd = *(double *)value;
    }
  else if (strcmp(key, "adsorbrate") == 0) 
    {
      data->adsorbrate = *(double *)value;
    }
  else if (strcmp(key, "carriername") == 0) 
    {
      sprintf(data->carriername, "%s", (char *)value);
      // fprintf(stderr,"data->carriername check: %s",data->carriername );

    }
  else if (strcmp(key, "dissolvedname") == 0) 
    {
      sprintf(data->dissolvedname, "%s", (char *)value);
    }
  else
    hd_quit("Unknown sediment tracer attribute '%s'\n", key);

  /* 
   * This is either not an auto-tracer or an auto-tracer with atleast
   * one of its attributes modified
   */
  sprintf(data->name, "user_spec");
}

/*-------------------------------------------------------------------*/
/* Reads the sediment related attributes      
           ./hd/tracers/tracer_info.c              */
/*-------------------------------------------------------------------*/
void sed_read_tr_atts(tracer_info_t *tr, FILE *fp, char *keyname)
{
  trinfo_priv_sed_t *data = get_private_data_sed(tr);
  char buf[MAXSTRLEN], key[MAXSTRLEN];
  double d_val;
  int i_val;

  /* Check for a sediment or ecology variable */
  if (!is_sed_var(tr->name)
#if defined(HAVE_ECOLOGY_MODULE)
      && !is_eco_var(tr->name)
#endif
      ) return;
  
  sprintf(key, "%s.psize", keyname);
  if (prm_read_double(fp, key, &d_val))
    sed_set_tr_att(tr, "psize", &d_val);
  
  sprintf(key, "%s.b_dens", keyname);
  if (prm_read_double(fp, key, &d_val))
    sed_set_tr_att(tr, "b_dens", &d_val);

  sprintf(key, "%s.i_conc", keyname);
  if (prm_read_double(fp, key, &d_val))
    sed_set_tr_att(tr, "i_conc", &d_val);

  sprintf(key, "%s.f_conc", keyname);
  if (prm_read_double(fp, key, &d_val))
    sed_set_tr_att(tr, "f_conc", &d_val);
 
//2019
 sprintf(key, "%s.css_erosion", keyname);
  if (prm_read_double(fp, key, &d_val))
    sed_set_tr_att(tr, "css_erosion", &d_val);

 sprintf(key, "%s.css_deposition", keyname);
  if (prm_read_double(fp, key, &d_val))
    sed_set_tr_att(tr, "css_deposition", &d_val);
 
  sprintf(key, "%s.svel", keyname);
  /* See if its a number */
  if (prm_read_double(fp, key, &d_val))
    sed_set_tr_att(tr, "svel", &d_val);
  else {
    /* Now try char */
    // xxx FR this does not work, prm_read_double will bail out
    if (prm_read_char(fp, key, buf))
      sed_set_tr_att(tr, "svel_name", buf);
  }

  sprintf(key, "%s.cohesive", keyname);
  if (prm_read_int(fp, key, &i_val))
    sed_set_tr_att(tr, "cohesive", &i_val);

  sprintf(key, "%s.calcvol", keyname);
  if (prm_read_int(fp, key, &i_val))
    sed_set_tr_att(tr, "calcvol", &i_val);

  sprintf(key, "%s.floc", keyname);
  if (prm_read_int(fp, key, &i_val))
    sed_set_tr_att(tr, "floc", &i_val);
  
  sprintf(key, "%s.resuspend", keyname);
  if (prm_read_int(fp, key, &i_val))
    sed_set_tr_att(tr, "resuspend", &i_val);

  sprintf(key, "%s.deposit", keyname);
  if (prm_read_int(fp, key, &i_val))
    sed_set_tr_att(tr, "deposit", &i_val);

  sprintf(key, "%s.adsorb", keyname);
  if (prm_read_int(fp, key, &i_val))
    sed_set_tr_att(tr, "adsorb", &i_val);

  sprintf(key, "%s.adsorbkd", keyname);
  if (prm_read_double(fp, key, &d_val))
    sed_set_tr_att(tr, "adsorbkd", &d_val);

  sprintf(key, "%s.adsorbrate", keyname);
  if (prm_read_double(fp, key, &d_val))
    sed_set_tr_att(tr, "adsorbrate", &d_val);

  sprintf(key, "%s.carriername", keyname);
  if (prm_read_char(fp, key, buf)) {
    // fprintf(stderr,"carriername check: %s",buf);
    sed_set_tr_att(tr, "carriername", buf);
  }

  sprintf(key, "%s.dissolvedname", keyname);
  if (prm_read_char(fp, key, buf))
    sed_set_tr_att(tr, "dissolvedname", buf);

  // Cannot read obc as its defined as an int
}

/* END sed_read_tr_atts()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Copies sedimet private data                                       */
/*-------------------------------------------------------------------*/
void *private_data_copy_sed(void *src)
{
  void *dest = NULL;

  if (src) {
    dest = (trinfo_priv_sed_t *)malloc(sizeof(trinfo_priv_sed_t));
    memcpy(dest, src, sizeof(trinfo_priv_sed_t));
  }

  return(dest);
}

/* END private_data_copy_sed()                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Perform any custom adjustments to the sediment tracer info  
  hd/master/hd_init.c    */
/*-------------------------------------------------------------------*/
void sed_tracer_custom(master_t *master)
{
  tracer_info_t *tr2d =  master->trinfo_2d;
  tracer_info_t *tr3d =  master->trinfo_3d;
  tracer_info_t *trsed =  master->trinfo_sed;
  char buf[MAXSTRLEN];
  int n, m;

  /* Check if a spatially varying settling velocity exists           */
  for (n = 0; n < master->atr; n++) {
    strcpy(buf, tr3d[n].name);
    strcat(buf, "_svel");
    for (m = 0 ; m < master->ntrS; m++) {
      if (strcmp(tr2d[m].name, buf) == 0) {
	trinfo_priv_sed_t *data = tr3d[n].private_data[TR_PRV_DATA_SED];
	if (data) {
	  strcpy(data->svel_name, buf);
	}
      }
    }
  }
}

/* END sed_tracer_custom()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Writes the attributes for sediment tracers                 
    /hd/tracers/tracer_info.c            */
/*-------------------------------------------------------------------*/
void sed_write_tr_atts(tracer_info_t *tr, FILE *fp, int n)
{
  trinfo_priv_sed_t *data = tr->private_data[TR_PRV_DATA_SED];

  if (data != NULL && data->type == PTR_SED) {
    fprintf(fp, "TRACER%1.1d.psize           %-4.2e\n", n, data->psize);
    fprintf(fp, "TRACER%1.1d.b_dens          %-4.2e\n", n, data->b_dens);
    fprintf(fp, "TRACER%1.1d.i_conc          %-4.2e\n", n, data->i_conc);
//2019
    fprintf(fp, "TRACER%1.1d.css_erosion          %-4.2e\n", n, data->css_erosion);
    fprintf(fp, "TRACER%1.1d.css_deposition          %-4.2e\n", n, data->css_deposition);

    if (strlen(data->svel_name))
      fprintf(fp, "TRACER%1.1d.svel_name       %s\n", n, data->svel_name);
    fprintf(fp, "TRACER%1.1d.svel            %f\n", n, data->svel);
    fprintf(fp, "TRACER%1.1d.cohesive        %d\n", n, data->cohesive);
    fprintf(fp, "TRACER%1.1d.calcvol         %d\n", n, data->calcvol);
    fprintf(fp, "TRACER%1.1d.floc            %d\n", n, data->floc);
    fprintf(fp, "TRACER%1.1d.resuspend       %d\n", n, data->resuspend);
    fprintf(fp, "TRACER%1.1d.deposit         %d\n", n, data->deposit);
    fprintf(fp, "TRACER%1.1d.adsorb          %d\n", n, data->adsorb);
    fprintf(fp, "TRACER%1.1d.adsorbkd        %-4.2e\n", n, data->adsorbkd);
    fprintf(fp, "TRACER%1.1d.adsorbrate      %-4.2e\n", n, data->adsorbrate);
    if (strlen(data->carriername))
      fprintf(fp, "TRACER%1.1d.carriername     %s\n", n, data->carriername);
    if (strlen(data->dissolvedname))
      fprintf(fp, "TRACER%1.1d.dissolvedname   %s\n", n, data->dissolvedname);
  }
}

/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Coundts the number of valid sediment classes                      */
/*-------------------------------------------------------------------*/
int count_sed_classes(char *sed_vars)
{
  char *vars[MAXSTRLEN * MAXNUMARGS], buf[MAXSTRLEN];
  int m, n, i, nvars = 0;

  if (sed_vars == NULL || strlen(sed_vars) == 0)
    return(0);

  strcpy(buf, sed_vars);
  m = parseline(buf, vars, MAXNUMARGS);
  for (n = 0; n < m; n++) {
    for (i = 0; i < NUM_SED_VARS; i++)
      if (strcmp(vars[n], SEDCLASS[i]) == 0) {
	nvars++;
	break;
      }
  }
  if (nvars == 0) {
    hd_quit("No valid sediment classes found.\n");
  }
  nvars += count_sed_3d();
  return(nvars);
}

int count_sed_3d()
{
  return(NUM_SED_VARS_3D);
}

int count_sed_2d()
{
  return(NUM_SED_VARS_2D);
}

/* END count_sed_classes()                                           */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Helper routine to set autotracer attributes       
    (called locally from routines below)             */
/*-------------------------------------------------------------------*/
static int sed_set_autotracer(FILE *fp, 
			      char *sed_vars,  /* list of vars */ 
			      char *sed_defs,  /* name of class */
			      tracer_info_t *trinfo, 
			      int ntr, int tn,
			      int   nsedclass,
			      char *sedclass[],
			      int trinfo_type)
{
  char *vars[MAXSTRLEN * MAXNUMARGS], buf[MAXSTRLEN];
  int m, n, i;
  int trn = tn;
  
  /* Check if the class is in the global configurations              */
  if (sed_vars == NULL) {
    strcpy(buf, sedclass[0]);
    for (n = 1; n < nsedclass; n++) sprintf(buf, "%s %s", buf, sedclass[n]);
    n = get_rendered_tracers(fp, 1, buf, sed_defs,
			     trinfo, ntr, tn, trinfo_type);
  } else
    n = get_rendered_tracers(fp, 1, sed_vars, sed_defs,
			     trinfo, ntr, tn, trinfo_type);
  if (n != tn) return(n);

  trn = tn;
  if ( sed_vars != NULL ) {
    /* Setup classes as specified */
    strcpy(buf, sed_vars);
    m = parseline(buf, vars, MAXNUMARGS);
    for (n = 0; n < m; n++) {
      for (i = 0; i < nsedclass; i++) {
	if (strcmp(vars[n], sedclass[i]) == 0) {
	  if ((tracer_find_index(sedclass[i], ntr, trinfo)) < 0) {
	    /* 
	     * Set the default value 
	     */
	    sed_set_tracer_defaults(&trinfo[tn], sedclass[i], sed_defs);
	    /*
	     * Overwrite any att from prm-file
	     */
	    sed_read_tr_atts(&trinfo[tn], fp, sedclass[i]);
	    trinfo[tn].m = -1; // does not exist in the prm-file
	    tracer_re_read(&trinfo[tn], fp, trinfo_type);
	    /*
	     * Fill in the tracer number
	     */
	    trinfo[tn].n = tn;
	    tn++;
	  }
	}
      }
    }
  } else {
    /* Add mandatory as per the arguments */
    for (i = 0; i < nsedclass; i++) {
      if ((tracer_find_index(sedclass[i], ntr, trinfo)) < 0) {
	/* 
	 * Set the default value 
	 */
	sed_set_tracer_defaults(&trinfo[tn], sedclass[i], sed_defs);
	/*
	 * Overwrite any att from prm-file
	 */
	sed_read_tr_atts(&trinfo[tn], fp, sedclass[i]);
	tracer_re_read(&trinfo[tn], fp, trinfo_type);
	/*
	 * Fill in the tracer number
	 */
	trinfo[tn].n = tn;
	tn++;
      }
    }
  }
  
  return(tn);
}

/*-------------------------------------------------------------------*/
/* Routine to initialise the 3D tracers in the master     
      ./hd/tracers/load_tracer.c      */
/*-------------------------------------------------------------------*/
int sediment_autotracer_3d(FILE *fp, int do_sed, char *sed_vars, char *sed_defs,
			    tracer_info_t *trinfo, int ntr, int tn)
{
  int i;
  /* SED CLASSES */
  if (do_sed && strlen(sed_vars)) {
    tn = sed_set_autotracer(fp, sed_vars, sed_defs, trinfo,  ntr, tn,
			    NUM_SED_VARS, SEDCLASS, WATER);
    
    /* SED 3D - mandatory*/
    tn = sed_set_autotracer(fp, NULL, sed_defs, trinfo,  ntr, tn,
			    NUM_SED_VARS_3D, SEDNAME3D, WATER);
  }
  return(tn);
}

/* END sediment_autotracer_3d()                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to initialise the 2D tracers in the master  
    ./hd/tracers/load_tracer.c           */
/*-------------------------------------------------------------------*/
int sediment_autotracer_2d(FILE *fp, int do_sed, char *sed_vars, char *sed_defs,
			   tracer_info_t *trinfo, int ntr, int tn)
{
  /* SED 2D - mandatory */
  if (do_sed && strlen(sed_vars)) {
    tn = sed_set_autotracer(fp, NULL, sed_defs, trinfo,  ntr, tn,
			    NUM_SED_VARS_2D, SEDNAME2D, INTER);
  }
  return(tn);
}

/* END sediment_autotracer_2d()                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to initialise the 2D tracers in the master      
     ./hd/tracers/load_tracer.c      */
/*-------------------------------------------------------------------*/
int sediment_autotracer_sed(FILE *fp, int do_sed, char *sed_vars, char *sed_defs,
			    tracer_info_t *trinfo, int ntr, int tn)
{

  /* Temperature and salinity */
  if ((tracer_find_index("salt", ntr, trinfo)) < 0) {
    strcpy(trinfo[tn].name, "salt");
    strcpy(trinfo[tn].long_name, "Salinity");
    strcpy(trinfo[tn].units, "PSU");
    trinfo[tn].valid_range_sed[0] = 0;
    trinfo[tn].valid_range_sed[1] = 40;
    trinfo[tn].fill_value_sed = 35.0;
    trinfo[tn].type = WATER|SEDIM|SEDIMENT|DIAGNOSTIC;
    trinfo[tn].inwc = 1;
    trinfo[tn].dissol = 1;
    trinfo[tn].advect = 1;
    trinfo[tn].diffuse = 1;
    strcpy(trinfo[tn].decay, "0.0");
    trinfo[tn].diagn = 0;
    trinfo[tn].n = tn;
    tn++;
  }
  if ((tracer_find_index("temp", ntr, trinfo)) < 0) {
    strcpy(trinfo[tn].name, "temp");
    strcpy(trinfo[tn].long_name, "Temperature");
    strcpy(trinfo[tn].units, "degrees C");
    trinfo[tn].valid_range_sed[0] = -4;
    trinfo[tn].valid_range_sed[1] = 40;
    trinfo[tn].fill_value_sed = 20.0;
    trinfo[tn].type = WATER|SEDIM|SEDIMENT|DIAGNOSTIC;
    trinfo[tn].inwc = 1;
    trinfo[tn].dissol = 1;
    trinfo[tn].advect = 1;
    trinfo[tn].diffuse = 1;
    strcpy(trinfo[tn].decay, "0.0");
    trinfo[tn].diagn = 0;
    trinfo[tn].n = tn;
    tn++;
  }

  /* SED CLASSES */
  if (do_sed && strlen(sed_vars)) {
    tn = sed_set_autotracer(fp, sed_vars, sed_defs, trinfo,  ntr, tn,
			    NUM_SED_VARS, SEDCLASS, SEDIM);
    
    /* SED 3D - mandatory */
    tn = sed_set_autotracer(fp, NULL, sed_defs, trinfo,  ntr, tn,
			    NUM_SED_VARS_3D, SEDNAME3D, SEDIM);
  }
  return(tn);
}

/* END sediment_autotracer_sed()                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Set the open boundary condition as specified in the sediment      */
/* private data for sediment autotracers.                            */
/*-------------------------------------------------------------------*/
int sed_get_obc(tracer_info_t *tr)
{
  trinfo_priv_sed_t *data = tr->private_data[TR_PRV_DATA_SED];
  int ret = NOGRAD;

  if (data != NULL) {
    if (data->type == PTR_SED && data->obc > 0)
      ret = data->obc;
  }
  return(ret);
}

/* END set_sed_obc()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Writes the sediment autotracers to the setup file      
 ./hd/diagnostics/run_setup.c            */
/*-------------------------------------------------------------------*/
int sediment_autotracer_write(master_t *master, FILE *op, int tn)
{
  parameters_t *params = master->params;
  char key[MAXSTRLEN];
  int i, n, m, trn, sn;

  if (!params->do_sed || strlen(params->sed_vars) == 0) return(tn);

  /* 3D water column and sediment autotracers */
  n = tn;
  for (m = 0; m < params->atr; m++) {
    tracer_info_t *tracer = &master->trinfo_3d[m];

    trn = -1;
    for (i = 0; i < NUM_SED_VARS; i++) {
      if (strcmp(tracer->name, SEDCLASS[i]) == 0) {
	trn = i;
	break;
      }
    }
    if (trn < 0) {
      for (i = 0; i < NUM_SED_VARS_3D; i++) {
	if (strcmp(tracer->name, SEDNAME3D[i]) == 0) {
	trn = i;
	break;
	}
      }
    }
    if (trn < 0) continue;

    /* Call the 3d helper function */
    tracer_write_3d(master, op, tracer, n);

    n++;
  }

  /* 2D sediment autotracers                                         */
  for (m = 0; m < master->ntrS; m++) {
    tracer_info_t *tracer = &master->trinfo_2d[m];

    trn = -1;
    for (i = 0; i < NUM_SED_VARS_2D; i++)
      if (strcmp(tracer->name, SEDNAME2D[i]) == 0) {
	trn = i;
	break;
      }
    if (trn < 0) continue;

    /* Call the 2d helper function */
    tracer_write_2d(master, op, tracer, n);

    n++;
  }
  return(n);
}

/* END sediment_autotracer_write()                                   */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Writes a sediment transport configuration to file   
        ./hd/diagnostics/run_setup.c      */
/*-------------------------------------------------------------------*/
void trans_write_sed(parameters_t *params, sediment_t *sediment, FILE *fp)
{
  sed_params_t *msparam = sediment->msparam;
  int i;
  double d1;
  char tag[MAXSTRLEN];

  fprintf(fp, "############################################################################\n");
  fprintf(fp, "# Sediments\n");

  fprintf(fp, "\n# Sediment layer thicknesses\n");
  fprintf(fp, "NSEDLAYERS      %d\n", params->sednz);
  for (i = 0; i < params->sednz; i++) {
    d1 =  params->gridz_sed[params->sednz - i] - params->gridz_sed[params->sednz - i - 1];
    fprintf(fp, "%5.4f\n", d1);
  }

  if (strlen(params->sed_vars)) {
    fprintf(fp,"\n#Sediment tracers and tracer attribute specification.\n");
    fprintf(fp, "# Avaibale sediment casses are:\n");
    for (i = 0; i < NUM_SED_VARS; i++)
      fprintf(fp, "# %s\n", SEDCLASS[i]);
    fprintf(fp,"SED_VARS        %s", params->sed_vars);
  }
  if (strlen(params->sed_defs))
    fprintf(fp,"\nSED_VARS_ATTS   %s\n", params->sed_defs);

  fprintf(fp, "\nVERBOSE_SED     %d\n", msparam->verbose_sed);
  fprintf(fp, "MAX_THICK_SED   %4.2f\n\n", msparam->max_thick_sed);

  fprintf(fp, "# Activate process\n");
  fprintf(fp, "# Default: 0\n");
  fprintf(fp, "DO_SEDIMENTS    YES\n");
  //  fprintf(fp, "RESUSPEND       YES\n");
  //  fprintf(fp, "DEPOSIT         YES\n");
  fprintf(fp, "CONSOLIDATE     %d\n\n", msparam->consolidate);

  fprintf(fp, "# CRITICAL STRESSES\n");
  fprintf(fp, "# Critical shear stress of cohesive bed erosion N m-2\n");
  fprintf(fp, "# Modes: CSS0, CSS1, CSS2, CSS3\n");
  fprintf(fp, "# Default: css_er = 0.2\n");
  if (msparam->cssmode == 0)
    fprintf(fp, "CSS_ER_MODE     CSS0\n");
  if (msparam->cssmode == 1)
    fprintf(fp, "CSS_ER_MODE     CSS1\n");
  if (msparam->cssmode == 2)
    fprintf(fp, "CSS_ER_MODE     CSS2\n");
  if (msparam->cssmode == 3)
    fprintf(fp, "CSS_ER_MODE     CSS3\n");
  if (msparam->cssmode == 4) {
    fprintf(fp, "CSS_ER_MODE     CSS4\n");
    fprintf(fp, "CSS_ER          %d\n", params->sednz);
    for (i = 0; i < params->sednz; i++) {
         fprintf(fp, "%5.4f\n", msparam->css_er_val[params->sednz-i-1]);
    }
    fprintf(fp, "\n");
    fprintf(fp, "CSS_ER_DEPTH    %d\n", params->sednz);
    for (i = 0; i < params->sednz; i++) {
      fprintf(fp, "%5.4f\n", msparam->css_er_depth[params->sednz-i-1]);
    }
  }
  fprintf(fp, "\n");

  fprintf(fp, "# Critical shear stress of cohesive sediment deposition N m-2\n");
  fprintf(fp, "# Default: CSS_DEP = +inf\n");
  fprintf(fp, "#CSS_DEP	        %4.2e\n\n", msparam->css_dep);

  fprintf(fp, "# FLOCCULATION\n");
  fprintf(fp, "# Modes: FLOC0, FLOC1, FLOC2, FLOC3\n");
  fprintf(fp, "# Default: FLOC0 - no flocculation\n");
  if (msparam->flocmode == 0)
    fprintf(fp, "FLOC_MODE       FLOC0\n");
  if (msparam->flocmode == 1)
    fprintf(fp, "FLOC_MODE       FLOC1\n");
  if (msparam->flocmode == 2)
    fprintf(fp, "FLOC_MODE       FLOC2\n");
  if (msparam->flocmode == 3)
    fprintf(fp, "FLOC_MODE       FLOC3\n");
  if (msparam->flocmode == 4) {
    fprintf(fp, "FLOC_MODE       FLOC4\n");
    fprintf(fp, "FLOC_PRM1       %4.3e\n", msparam->flocprm1);
    fprintf(fp, "FLOC_PRM2       %4.3e\n\n", msparam->flocprm2);
  }
  fprintf(fp, "\n");

  fprintf(fp, "# CONSOLIDATION\n");
  fprintf(fp, "FINPOR_SED      %4.2f\n", msparam->finpor_sed);
  fprintf(fp, "CONSOLRATE      %4.2f\n\n", msparam->consolrate);

  fprintf(fp, "# RIPPLES\n");
  fprintf(fp, "# Bed forms\n");
  fprintf(fp, "# Default: 0\n");
  fprintf(fp, "CALC_RIPPLES    %d\n", msparam->calc_ripples);
  if (msparam->physriph_spv)
    fprintf(fp, "PHYSRIPH spatially varying\n");
  else
    fprintf(fp, "PHYSRIPH        %5.4f\n", msparam->physriph);
  if (msparam->physripl_spv)
    fprintf(fp, "PHYSRIPL spatially varying\n");
  else
    fprintf(fp, "PHYSRIPL        %5.4f\n", msparam->physripl);
  fprintf(fp, "BIORIPH         %5.4f\n", msparam->bioriph);
  fprintf(fp, "BIORIPL         %5.4f\n", msparam->bioripl);
  fprintf(fp, "# Skin roughness\n");
  fprintf(fp, "Z0_SKIN         %5.3e\n\n", msparam->quad_bfc);

  fprintf(fp, "# Maximum depth for biological activity\n");
  if (msparam->maxbiodepth_spv)
    fprintf(fp, "MAXBIODEPTH spatially varying\n");
  else
    fprintf(fp, "MAXBIODEPTH     %4.2f\n\n", msparam->maxbiodepth);

  fprintf(fp, "# Functional form for bioirrigarion and bioturbation activity\n");
  fprintf(fp, "# as a function of depth. Currently can be one of:\n");
  fprintf(fp, "#     constant\n");
  fprintf(fp, "#     linear\n");
  fprintf(fp, "#     parabolic\n");
  fprintf(fp, "#     gaussian\n");
  fprintf(fp, "# Only the first letter of this parameter is significant\n");
  if (msparam->biosedprofile == 'p')
    fprintf(fp, "BIOSEDPROFILE   parabolic\n\n");
  else if (msparam->biosedprofile == 'c')
    fprintf(fp, "BIOSEDPROFILE   constant\n\n");
  else if (msparam->biosedprofile == 'l')
    fprintf(fp, "BIOSEDPROFILE   linear\n\n");
  else if (msparam->biosedprofile == 'g')
    fprintf(fp, "BIOSEDPROFILE   gaussian\n\n");
  else {
    // Make parabolic the default for RECOM
    fprintf(fp, "BIOSEDPROFILE   parabolic\n\n");
  }
  fprintf(fp, "# Biological activity (standard animals per square metre?)\n");
  if (msparam->biodens_spv)
    fprintf(fp, "BIODENS spatially varying\n");
  else
    fprintf(fp, "BIODENS	       %4.2f\n\n", msparam->biodens);

  fprintf(fp, "# Diffusion coefficient for bio-irrigation of sediments.\n");
  fprintf(fp, "# This value is scaled by the amount of biological\n");
  fprintf(fp, "# activity present, and also decreases with depth\n");
  fprintf(fp, "# in the sediment according to some\n");
  fprintf(fp, "# fixed profile. The value here is the value which\n");
  fprintf(fp, "# would apply at zero depth in the sediment.\n");
  fprintf(fp, "# Units are m2 s-1 per animal per m2.\n");
  if (msparam->bi_dissol_kz_spv)
    fprintf(fp, "BI_DISSOL_KZ spatially varying\n");
  else
    fprintf(fp, "BI_DISSOL_KZ    %4.2e\n", msparam->bi_dissol_kz);
  if (msparam->bi_dissol_kz_i_spv)
    fprintf(fp, "BI_DISSOL_KZ_I spatially varying\n");
  else
    fprintf(fp, "BI_DISSOL_KZ_I  %4.2e\n\n", msparam->bi_dissol_kz_i);

  fprintf(fp, "# Diffusion coefficient for bio-turbation of sediments.\n");
  fprintf(fp, "# This value is scaled by the amount of biological\n");
  fprintf(fp, "# activity present, and also decreases with depth\n");
  fprintf(fp, "# in the sediment according to some\n");
  fprintf(fp, "# fixed profile. The value here is the value which\n");
  fprintf(fp, "# would apply at zero depth in the sediment.\n");
  fprintf(fp, "# Units are m2 s-1 per animal per m2.\n");
  if (msparam->bt_partic_kz_spv)
    fprintf(fp, "BT_PARTIC_KZ spatially varying\n");
  else
    fprintf(fp, "BT_PARTIC_KZ    %4.2e\n", msparam->bt_partic_kz);
  if (msparam->bt_partic_kz_i_spv)
    fprintf(fp, "BT_PARTIC_KZ_I spatially varying\n");
  else
    fprintf(fp, "BT_PARTIC_KZ_I  %4.2e\n\n", msparam->bt_partic_kz_i);

  fprintf(fp, "HINDERED_SVEL  %d \n\n", 1);
  fprintf(fp, "REEF_SCALING_DEPTH  %lf \n\n", 70.0);


  fprintf(fp, "############################################################################\n\n");
}

/* END trans_write_sed()                                             */
/*-------------------------------------------------------------------*/
#endif
