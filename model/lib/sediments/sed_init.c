/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/sediments/sed_init.c
 *  
 *  Description:
 *  Initialise sediment module
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: sed_init.c 5848 2018-06-29 05:01:15Z riz008 $
 *
 */

#ifdef __cplusplus
#include "stdafx.h"
extern "C" {
#endif

#include "sediments.h"

#if defined(HAVE_SEDIMENT_MODULE)
/* Initialise the sediment structure.
 */
static void sed_params_init(FILE * prmfd, sediment_t *sediment);
static void sed_tracers_init(FILE * prmfd, sediment_t *sediment);
static void sed_grid_init(FILE * prmfd, sediment_t *sediment, sed_column_t *sm, int write_log);
static void alloc_sed_spatial(sediment_t *sediment);
static void sed_optimization(sediment_t *sediment);

static FILE *get_sed_params(FILE *fp, sediment_t *sediment);
static void sed_params_build(sediment_t *sediment);
static void sed_params_est(sediment_t *sediment);
static void sed_params_std(sediment_t *sediment);
static void sed_params_shf(sediment_t *sediment);
static void sed_params_bsc(sediment_t *sediment);
static void sed_params_auto(sediment_t *sediment);
static void sed_params_write(sediment_t *sediment);
void sed_cleanup(sediment_t *sediment);
double sinterface_getmodeltime(void* hmodel);
double sinterface_getmodeltimestep(void* model);
int sinterface_getnumberwclayers(void* hmodel);
int sinterface_getnumbercolumns(void* hmodel);
int sinterface_getnumbersedlayers(void* hmodel);
int sinterface_getnumberoftracers(void* hmodel);
int si_getmap2hdwctracer(void* hmodel, int n, char *name);
int si_getmap2hdsedtracer(void* hmodel, int n, char *name);
int si_getmap2filetracer(void* hmodel, FILE *fp, int n, char *name);
  //2010
  int sinterface_getcss_er_val(FILE* prmfd, double *css_er_val);
  int sinterface_getcss_er_depth(FILE* prmfd, double *css_er_depth);
//2016
double sinterface_erflux_scale(FILE* prmfd);

/** UR added */
FILE* si_getparamfile_tracer(FILE *fp);
FILE* si_getparamfile_sed(FILE *fp);
/* end UR */
  char* _itoa(int value, char* str, int base);
  void strreverse(char* begin, char* end);

void sinterface_getmap2diagtracer_2d(void* hmodel,
   int *n_hripple, int *n_lripple, int *n_ustrcw_skin,
   int *n_depth_sedi, int *n_dzactive, int *n_erdepflux_total,
   int *n_erdepflux_total_ac);
  /*
void sinterface_getmap2diagtracer_3d(void* hmodel,
   int *n_tss, int *n_svel_floc, int *n_por_sed, int *n_coh_sed);
  */

double sinterface_getquad_bfc(FILE* prmfd);

int sinterface_getverbose_sed(FILE* prmfd);
int sinterface_getgeomorph(FILE* prmfd);
int sinterface_getconsolidate(FILE* prmfd);
double sinterface_getfinpor_sed(FILE* prmfd);
double sinterface_getconsolrate(FILE* prmfd);
int sinterface_getcssmode(FILE* prmfd);
double sinterface_getcss(FILE* prmfd);
double sinterface_getcssdep(FILE* prmfd);
int sinterface_getflocmode(FILE* prmfd);
double sinterface_getflocprm1(FILE* prmfd);
double sinterface_getflocprm2(FILE* prmfd);
int sinterface_getbblnonlinear(FILE* prmfd);
  double sinterface_getcssscale(FILE* prmfd);//Nov12
int sinterface_getcalcripples(FILE* prmfd);
double sinterface_getphysriph(FILE* prmfd);
double sinterface_getbioriph(FILE* prmfd);
double sinterface_getphysripl(FILE* prmfd);
double sinterface_getbioripl(FILE* prmfd);
double sinterface_getbiodens(FILE* prmfd);
double sinterface_getmaxbiodepth(FILE* prmfd);
double sinterface_getbi_dissol_kz(FILE* prmfd);
double sinterface_getbt_partic_kz(FILE* prmfd);

double sinterface_getbi_dissol_kz_i(FILE* prmfd);
double sinterface_getbt_partic_kz_i(FILE* prmfd);

char sinterface_getbiosedprofile(FILE* prmfd);
double sinterface_getz0(FILE* prmfd);
int sinterface_getdzinit_sed(FILE* prmfd, double *dz_sed, int sednzc);

double sinterface_gettheta(void *hmodel, double *theta, int ncol);

void si_set_errfn_warn(void *hmodel);
int si_gettracernames(void* hmodel, char **tracername);

double sinterface_getmaxthicksed(FILE* prmfd);
  
  //2012 
  void sed_tracers_benthic_init(FILE * prmfd, sediment_t *sediment);
  int sinterface_getnumberofBtracer(void* hmodel);
  void sinterface_getnameofBtracer(void* hmodel, int n, char *tracername);

  void sinterface_get_tracerunits(void* model, char *name, char *units);
  double sinterface_get_fillvalue_wc(void* model, char *name);
  double sinterface_get_fillvalue_sed(void* model, char *name);
  double sinterface_get_decay(void* model, char *name);
  double sinterface_get_psize(void* model, char *name);
  double sinterface_get_b_dens(void* model, char *name);
  double sinterface_get_i_conc(void* model, char *name);
  double sinterface_get_svel(void* model, char *name);
  void sinterface_get_svel_name(void* model, char *name, char *sname);
  int sinterface_get_diagn(void* model, char *name);
  int sinterface_get_dissol(void* model, char *name);
  int sinterface_get_partic(void* model, char *name);
  int sinterface_get_adsorb(void* model, char *name);
  int sinterface_get_diffuse(void* model, char *name);
  int sinterface_get_cohesive(void* model, char *name);
  int sinterface_get_floc(void* model, char *name);
  int sinterface_get_resuspend(void* model, char *name);
  int sinterface_get_deposit(void* model, char *name);
  int sinterface_get_calcvol(void* model, char *name);
  int sinterface_get_adsorb_kd(void* model, char *name);
  int sinterface_get_adsorb_rate(void* model, char *name);
  int sinterface_get_carriername(void* model, char *name,char *carriername );
  int sinterface_get_dissolvedname(void* model, char *name, char *dissolvedname);
  int sinterface_gethindered_svel_patch(FILE* prmfd);
  int sinterface_gethindered_svel(FILE* prmfd);
  double sinterface_reef_scale_depth(FILE* prmfd);

int sinterface_getshipfile(FILE* prmfd, char *shipfile);
#if defined(HAVE_OMP)
int sinterface_get_trans_num_omp(void *model);
#endif


/* Version information */
int get_sediments_major_vers(void)
{
  return(SEDIMENTS_MAJOR_VERSION);
}

int get_sediments_minor_vers(void)
{
  return(SEDIMENTS_MINOR_VERSION);
}

int get_sediments_patch_vers(void)
{
  return(SEDIMENTS_PATCH_VERSION);
}
  
sediment_t *sed_init(FILE * prmfd, void *hmodel)
{
  FILE* fd;
  char sedparams[MAXSTRLEN];
  /* Allocate memory for sediment Structure */
  sediment_t *sediment = (sediment_t *)malloc(sizeof(sediment_t));
  if (sediment == NULL) {
      sedtag(LFATAL,"sed:sed_init:sed_init",
        " Could not allocate memory for sediment_t structure.");
      exit(1);
  }
  memset(sediment, 0, sizeof(sediment_t));

  sediment->hmodel = hmodel;

  /* Initialise sediment parameters structure */
  /* MH: 08.2012. included automated parameter specification
  fd = si_getparamfile_sed(prmfd);
  sed_params_init(fd, sediment);
  if(fd != prmfd)
  	fclose(fd);
  */
  if ((fd = get_sed_params(prmfd, sediment))) {
    sed_params_init(fd, sediment);
    if(fd != prmfd)
      fclose(fd);
  }
  /* END MH */
  	
  sed_params_write(sediment);// log file

  /* Initialise sediment tracer structure */
  fd = si_getparamfile_tracer(prmfd);
  sed_tracers_init(fd, sediment);
  if(fd != prmfd)
  	fclose(fd);

  /* Initialise column structure */
  sediment->sm = alloc_sed_column(sediment);
  /* Initialise miscellaneous variables */
  sed_grid_init(prmfd, sediment, sediment->sm, 1);

#if defined(HAVE_OMP)
  /* Initialise number of columns for parallel execution  */
  int i, nomp = sinterface_get_trans_num_omp(hmodel);
  if (nomp > 1) {
    sediment->smn = (sed_column_t **)p_alloc_1d(nomp);
    for (i=0; i<nomp; i++) {
      sediment->smn[i] = alloc_sed_column(sediment);
      /* Initialise miscellaneous variables */
      sed_grid_init(prmfd, sediment, sediment->smn[i], 0);
    }
  }
#endif
  
  /* Initialise internal 3-d variables */
  alloc_sed_spatial(sediment);

  /*UR Apply some optimization for conditional statements */
  sed_optimization(sediment);

  return sediment;
}

/* Close down the window, freeing and dellocate memory, etc.
 */
void sed_cleanup(sediment_t *sediment)
{
  sed_params_t *param = sediment->msparam;
  sed_tracer_t **tracer = &sediment->mstracers;
  sed_column_t *sm = sediment->sm;

  if (param->sednz < 1)
    return;

  d_free_1d(sm->partic_kz);
  d_free_1d(sm->dissol_kz);
  d_free_1d(sm->por_sed);
  d_free_1d(sm->porold_sed);
  d_free_1d(sm->css);
  d_free_1d(sm->cellz_sed);
  d_free_1d(sm->gridz_sed);
  d_free_1d(sm->dz_sed);
  d_free_2d(sm->tr_sed);
  d_free_2d(sm->tr_wc);
  d_free_1d(sm->dz_wc);
  d_free_1d(sm->u1_wc);
  d_free_1d(sm->u2_wc);
  d_free_1d(sm->Kz_wc);
  if (sediment != NULL) {
    free(sm);
    free(tracer);
    free(param);
    free(sediment);
  }
}

/*******************************************************************/

static void sed_params_init(FILE * prmfd, sediment_t *sediment)
{

  int i,j,k;
  char shipfile[MAXSTRLEN],buf[MAXSTRLEN];
  void *hmodel = sediment->hmodel;
  FILE *flog, *fp;

  /* Allocate memory for sed_params_t Structure */
  sed_params_t *param = (sed_params_t *)malloc(sizeof(sed_params_t));
  if (param == NULL) {
      sedtag(LFATAL,"sed:sed_init:sed_tracers_init"," Could not allocate memory for a sediment parameters.");
      exit(1);
  }
  memset(param, 0, sizeof(sed_params_t));
  sediment->msparam = param;

    /* Initialise MecoSedParam Structure */
    /* Data available from hd */
    /*
       param->calcvol_wc  will be defined in MecoSedTracers_init
       param->calcvol_sed  will be defined in MecoSedTracers_init
    */

  param->verbose_sed = sinterface_getverbose_sed(prmfd);

  param->t=sinterface_getmodeltime(hmodel);
  param->dt =  sinterface_getmodeltimestep(hmodel);
  param->nz = sinterface_getnumberwclayers(hmodel)+1;
  /* If there is only one water layer, mecosed will split it in two layers*/
  if(param->nz == 2)
      param->nz = 3;
  param->sednz = sinterface_getnumbersedlayers(hmodel)+1;
  /* Note that in mecosed param->sednz and param->nz are defined as
     number of interfaces (not sediment layers) */
  if(param->sednz < 3) {
      sedtag(LFATAL,"sed:sed_init:sed_tracers_init","At least two sediment layers must be specified");
      exit(1);
  }

  param->nstep=0;
  param->ncol =  sinterface_getnumbercolumns(hmodel);
  param->ntr = sinterface_getnumberoftracers(hmodel);

  sinterface_getmap2diagtracer_2d(hmodel,
          &param->n_hripple,  &param->n_lripple,
          &param->n_ustrcw_skin,
          &param->n_depth_sedi,  &param->n_dzactive,
          &param->n_erdepflux_total,
          &param->n_erdepflux_total_ac);
  /*
  sinterface_getmap2diagtracer_3d(hmodel,
          &param->n_tss,  &param->n_svel_floc,
          &param->n_por_sed,  &param->n_coh_sed);
  */
  si_set_errfn_warn(hmodel);

  param->quad_bfc = sinterface_getquad_bfc(prmfd);
  /* Read data from the parameter file */
  /* param->mindepth_sedi will be defined in sed_grid_init as
     a fraction of the initial sediment-thickness */
  param->geomorph = sinterface_getgeomorph(prmfd);
  param->consolidate =  sinterface_getconsolidate(prmfd);
  param->minpor_wc = 0.1;
  param->minpor_sed = 0.1;
  param->minseddz = 0.000001;
  param->finpor_sed = sinterface_getfinpor_sed( prmfd);
  param->consolrate =  sinterface_getconsolrate( prmfd);
  param->cssmode = sinterface_getcssmode(prmfd);
  param->css = sinterface_getcss(prmfd);
 
 //2010
  param->css_er_val = d_alloc_1d(param->sednz); //2010
  param->css_er_depth = d_alloc_1d(param->sednz); //2010
  if (param->cssmode == 4) {
    i=sinterface_getcss_er_val(prmfd, param->css_er_val);
    if ( i+1 != param->sednz) {
      sedtag(LFATAL,"sed_init", "Number of css values must be equal number of sed layers!");
      exit(1);
    }
    j=sinterface_getcss_er_depth(prmfd, param->css_er_depth);
    if ( j != i) {
      sedtag(LFATAL,"sed_init", "Number of css levels must be equal number of css values!");
      exit(1);
    }
  }

  // ships
  param->ship_N=0;
  param->ship_Kz=0.0;
  k = sinterface_getshipfile(prmfd, shipfile);
  // fprintf(stderr,"shipfile=%s",shipfile);
  if(k) {
   // read from the shipfile 
    //the number of grid-cells on track (ship_C)
    //the number of ships per day (ship_N)
    //time on track for a ship in hours (ship_T)
    // i.j cells comprising ship-track
  fp = fopen(shipfile,"r");
  fscanf(fp,"%s",buf);
  fscanf(fp,"%s",buf);
  fscanf(fp,"%s",buf);
  fscanf(fp,"%s",buf);
  fscanf(fp,"%s",buf);
  fscanf(fp,"%d",&param->ship_C);
  fscanf(fp,"%d",&param->ship_C);
  fscanf(fp,"%d",&param->ship_N);
  fscanf(fp,"%lf",&param->ship_T);
  fscanf(fp,"%lf",&param->ship_Kz);
  // fprintf(stderr,"C=%d N=%d \n", param->ship_C, param->ship_N  );
  param->ship_i = (int *) i_alloc_1d(param->ship_C);
  param->ship_j = (int *) i_alloc_1d(param->ship_C);
  for (i=0;i<param->ship_C;i++) {
     fscanf(fp,"%d",&param->ship_i[i]);
     fscanf(fp,"%d",&param->ship_j[i]);
  }
 fclose(fp);
 }


  

  param->css_dep = sinterface_getcssdep(prmfd);
  param->flocmode = sinterface_getflocmode(prmfd);
  param->flocprm1 = sinterface_getflocprm1( prmfd);
  param->flocprm2 = sinterface_getflocprm2( prmfd);
  param->bbl_nonlinear = sinterface_getbblnonlinear( prmfd);
  param->hindered_svel_patch = sinterface_gethindered_svel_patch( prmfd);
  param->hindered_svel = sinterface_gethindered_svel( prmfd);
  param->reef_scale_depth = sinterface_reef_scale_depth( prmfd);
  param->calc_ripples = sinterface_getcalcripples( prmfd);
  param->physriph = sinterface_getphysriph( prmfd);
  param->css_scale = sinterface_getcssscale( prmfd); //nov12
  param->bioriph = sinterface_getbioriph( prmfd);
  param->physripl = sinterface_getphysripl( prmfd);
  param->bioripl = sinterface_getbioripl( prmfd);
  param->biodens= sinterface_getbiodens( prmfd);
  param->maxbiodepth= sinterface_getmaxbiodepth( prmfd);
  param->bi_dissol_kz= sinterface_getbi_dissol_kz( prmfd);
  param->bt_partic_kz= sinterface_getbt_partic_kz( prmfd);
  param->bi_dissol_kz_i= sinterface_getbi_dissol_kz_i( prmfd);
  param->bt_partic_kz_i= sinterface_getbt_partic_kz_i( prmfd);
  param->biosedprofile = sinterface_getbiosedprofile( prmfd);

  param->max_thick_sed = sinterface_getmaxthicksed(prmfd);//2010

param->erflux_scale = sinterface_erflux_scale(prmfd); //2016

}

/**********************************/

static void sed_tracers_init(FILE * prmfd, sediment_t *sediment)
{
  int i,k,n,nn;
  FILE *flog;
  char carriername[MAXSTRLEN], dissolvedname[MAXSTRLEN];
  char buf[MAXSTRLEN], buf1[MAXSTRLEN];
  sed_params_t *param = sediment->msparam;
  sed_tracer_t **mstracers = &sediment->mstracers;
  void *hmodel = sediment->hmodel;
  int ntr = param->ntr;

  *mstracers = NULL;
  /* Allocate memory for array of sed_tracer_t and fill with zeros */
  *mstracers = (sed_tracer_t *)malloc(sizeof(sed_tracer_t) * ntr);
  memset(*mstracers, 0, sizeof(sed_tracer_t) * ntr);
 	if(!param->ntr)
 	{
 	  sedtag(LFATAL,"sed:sed_init:sed_tracers_init"," No sed tracers, ntr == 0");
  	exit(1);
 	}
/* get tracer names beforehand, as they will be used to
to map mecosed tracers to those in hd module and in prm-file
*/
 {
  char **tracername = (char **)p_alloc_1d(param->ntr);
  int m;
    for(n=0; n < param->ntr; n++)
  tracername[n] = malloc(MAXSTRLEN *sizeof(char));
    m=si_gettracernames(hmodel, tracername);
    if(m!=param->ntr) {
  sedtag(LFATAL,"sed:sed_init:sed_tracers_init"," -> si_gettracernames");
  exit(1);
    }
    for (n = 0; n < ntr; ++n) {
      sed_tracer_t *tr = &sediment->mstracers[n];
      strcpy(tr->name, tracername[n]);
  //	fprintf(stderr,"n=%d, tr->name=%s \n",n, tr->name);
    }
    free(tracername);
  }


  /* Read tracer attributes for all tracers */
  param->calcvol_wc = 0;
  param->calcvol_sed = 0;
  for (n = 0; n < ntr; ++n) {
    sed_tracer_t *tr = &sediment->mstracers[n];
    int nf, col;

    tr->n = n;
    /* maps */
    tr->n_hd_wc = si_getmap2hdwctracer(hmodel, n, tr->name);
    if(tr->n_hd_wc<0) {
      sedtag(LFATAL,"sed:sed_init:sed_tracers_init"," -> si_getmap2hdwctracer, failed to retrieve name: %s",tr->name);
      exit(1);
    }
    tr->n_hd_sed = si_getmap2hdsedtracer(hmodel, n, tr->name);
    if(tr->n_hd_sed<0) {
      sedtag(LFATAL,"sed:sed_init:sed_tracers_init"," -> si_getmap2hdsedtracer, failed to retrieve name: %s",tr->name);
      exit(1);
    }

    tr->inwc = 1;
    tr->insed = 1;
    sinterface_get_tracerunits(hmodel, tr->name, tr->units);
    tr->fill_value_wc = sinterface_get_fillvalue_wc(hmodel, tr->name);
    tr->fill_value_sed = sinterface_get_fillvalue_sed(hmodel, tr->name);
    tr->diagn = sinterface_get_diagn(hmodel, tr->name);
    tr->dissol = sinterface_get_dissol(hmodel, tr->name);
    tr->partic = sinterface_get_partic(hmodel, tr->name);
    tr->adsorb = sinterface_get_adsorb(hmodel, tr->name);

    //NMY2015
    // fprintf(stderr,"sed_init: adsorb=%d, name = %s \n", tr->adsorb,tr->name);

    if(tr->dissol && tr->partic) {
      sedtag(LFATAL,"Sediments:sed_init:sed_tracers_init","Tracer %s must be either partic or dissol \n", tr->name);
      exit(1);
    }
    if(!tr->dissol && !tr->partic) tr->dissol = 1;
    tr->diffuse =  sinterface_get_diffuse(hmodel, tr->name);
    tr->decay = sinterface_get_decay(hmodel, tr->name);
    tr->psize = sinterface_get_psize(hmodel, tr->name);
    tr->b_dens = sinterface_get_b_dens(hmodel, tr->name);
    tr->i_conc = sinterface_get_i_conc(hmodel, tr->name);
    tr->svel = d_alloc_1d(param->ncol);
    for (col = 0; col < param->ncol; col++)
      tr->svel[col] = sinterface_get_svel(hmodel, tr->name);
    /* END MH */

    /* rescale units from mg/m3 to kg/m3 */
    if (tr->partic  && !tr->adsorb && strcmp(tr->units,"kg m-3") != 0 && strcmp(tr->units,"kgm-3") != 0)
    {
      tr->b_dens=tr->b_dens/1.e+6;
      tr->i_conc=tr->i_conc/1.e+6;
      sedtag(LWARN,"Sediments:sed_init:sed_tracers_init", "Particulate tracer %s units are not kg m-3, and are assumed to be mg m-3 \n", tr->name);
      tr->u_scale=1e+6;
    }else
      tr->u_scale = 1;


    /* Read data from prm file */
    if(tr->partic) { /* if the tracer is particulate or adsorbed */
      /* MH: 08.2012. Read tracer attributes via the trinfo private data  */
      tr->cohesive = sinterface_get_cohesive(hmodel, tr->name);
      tr->floc = sinterface_get_floc(hmodel, tr->name);
      tr->resuspend = sinterface_get_resuspend(hmodel, tr->name);
      tr->deposit = sinterface_get_deposit(hmodel, tr->name);
      tr->calcvol_sed = sinterface_get_calcvol(hmodel, tr->name);
      /* END MH */

      if (param->geomorph)
        tr->calcvol_wc = tr->calcvol_sed;
      else
        tr->calcvol_wc = 0;

     /* tr->adsorb = sinterface_getadsorb(prmfd, nf);	*/
      if( tr->calcvol_wc + tr->calcvol_sed ) {
        tr->adsorb=0;
        tr->carriernum = -1;
        tr->dissolvednum = -1;
        tr->adsorbkd = -1.;
        tr->adsorbrate = -1.;
      }
      else
      {
        if(tr->adsorb) { /* if the tracer is adsorbed */
          tr->calcvol_wc  = 0;
          tr->calcvol_sed = 0;
	  tr->carriernum = 0;
	  tr->dissolvednum = 0;
	  /* MH: 08.2012. Read tracer attributes via the trinfo private data */
	  sinterface_get_carriername(hmodel, tr->name, carriername);
	  sinterface_get_dissolvedname(hmodel, tr->name, dissolvedname);

    //NMY2015
    fprintf(stderr,"sed_init: carrier name = %s \n", carriername);

	  /* END MH */
	  for (nn=0; nn<ntr;nn++) {
	    sed_tracer_t *trr = &sediment->mstracers[nn];
	    if ( strcmp(carriername, trr->name) == 0)
	      tr->carriernum = nn;
	    if ( strcmp(dissolvedname, trr->name) == 0)
	      tr->dissolvednum = nn;          
	  }
	  if (tr->carriernum == 0 || tr->dissolvednum == 0) {
	    fprintf(stderr,"ERROR:sed_init:sed_tracers_init: either carriername or dissolvedname missed in prm file \n"); exit(1);
	  }
       /*
	* MH: 08.2012. Read tracer attributes via the trinfo private data 
	*/
	  tr->adsorbkd = sinterface_get_adsorb_kd(hmodel, tr->name);
	  tr->adsorbrate = sinterface_get_adsorb_rate(hmodel, tr->name);
	  /* END MH */
        }else { /* if dissolved num is not specified then
		   there is no sorption/desorption exchange */
          tr->dissolvednum = -1;
          tr->adsorbkd = -1.;
          tr->adsorbrate = -1.;
        }
      }
    }else /* if the tracer is dissolved */
    {
        tr->cohesive=0;
        tr->calcvol_wc=0;
        tr->calcvol_sed=0;
        tr->adsorb=0;
        tr->carriernum = -1;
        tr->dissolvednum = -1;
        tr->adsorbkd = -1.;
        tr->adsorbrate = -1.;
    }

    param->calcvol_wc += tr->calcvol_wc;
    param->calcvol_sed += tr->calcvol_sed;
  }


  if (param->calcvol_sed < 1){
      sedtag(LFATAL,"sed:sed_init:sed_tracers_init","At least one volumetric particulate tracer must be specified \n");
      exit(1);
  }


 // 2012 subregions
  //spatially varying parameters 
  //nov11
 {
    int m, col;
    int np;  //max number of prm
    // read names and total number of benthic (ie 2d) variables
    sed_tracers_benthic_init(prmfd, sediment);
    // set max number of spatially varying prms
    np = param->ntr + 30;
    //allocate mem
    param->prmpointS = (double ***)p_alloc_2d(param->ncol,np);
    param->prmnameS = (char **)p_alloc_1d(np);
    param->prmindexS = i_alloc_1d(np);
    param->trindexS = i_alloc_1d(np);
    for(n=0; n < np; n++)
      param->prmnameS[n] = malloc(MAXSTRLEN *sizeof(char));
    for (i=0;i<np;i++) {
      for (col = 0; col < param->ncol; col++)
	param->prmpointS[i][col] = NULL;
      sprintf(param->prmnameS[i], "%c", '\0');
    }
    // fprintf(stderr,"%s \n",  param->prmnameS[3]);
    //fill up the list of the candidate prm-names and the list of pointers to corresonding values
    // start with settling velocities
    for (n=0; n<param->ntr; n++) {
      sed_tracer_t *tr = &sediment->mstracers[n];
      /* MH: 08/2012. Reference spatially variable parameters by name
      strcpy(buf,"TRACER"); strcat(buf,_itoa(tr->n_file,buf1,10));strcat(buf,".svel");
      */
      sprintf(buf, "%c", '\0');
      sinterface_get_svel_name(hmodel, tr->name, buf);

      // fprintf(stderr, "n=%d tr->name=%s \n", n, buf);

      //NM commented out MH above
      // strcpy(buf,"TRACER"); strcat(buf,_itoa(tr->n_file,buf1,10)); strcat(buf,".svel");

    // NM if we decide to get rid of numbers in prm tracer names then
      // we can specify prmnameS for svel as follows 
      // (no dots in names because dive does not like that) 
      // strcpy(buf,tr->name); strcat(buf,"_svel");

      strcpy(param->prmnameS[n],buf);
      for (col = 0; col < param->ncol; col++)
	param->prmpointS[n][col] = &tr->svel[col];
    }

    // proceeed with other parameters
    i= param->ntr;
    // FR 07-2015: Note there are no assignments to prmpointS[] done
    // here as the spattially varying case is handled by the next loop
    // using the *_spv variables
    strcpy(param->prmnameS[i+0],"PHYSRIPH");
    strcpy(param->prmnameS[i+1],"PHYSRIPL");
    strcpy(param->prmnameS[i+2],"MAXBIODEPTH");
    strcpy(param->prmnameS[i+3],"BIODENS");
    strcpy(param->prmnameS[i+4],"BT_PARTIC_KZ");
    strcpy(param->prmnameS[i+5],"BT_PARTIC_KZ_I");
    strcpy(param->prmnameS[i+6],"BI_DISSOL_KZ");
    strcpy(param->prmnameS[i+7],"BI_DISSOL_KZ_I");
    strcpy(param->prmnameS[i+8],"CSS_SCALE");

    // map spatially varying parameter to the corresponding tracer
    param->nprmS=0;
    for(n=0; n < param->ntrB; n++) {
      for (k=0;k<np;k++) {
	if (strcmp(param->trnameB[n],param->prmnameS[k]) == 0){
	  i=param->nprmS;
	  param->prmindexS[i] = k; // prm index (eg param->prmpointS[k]) 
	  param->trindexS[i] = n; // tracer index  
	  // tracer value at location c is acecssed in hd2sed.c through 
	  // sinterface_getvalueofBtracer(void* hmodel, int n, int c)
	  // and then that value is assigned to  *param->prmpointS[k]
	  param->nprmS += 1; 
	  // Allocate memory pointed to the *_spv variables
	  if (strcmp(param->prmnameS[k],"PHYSRIPH") == 0) {
	    param->physriph_spv = d_alloc_1d(param->ncol);
	    for (col = 0; col < param->ncol; col++)
	      param->prmpointS[k][col] = &param->physriph_spv[col];
	  } else
	    if (strcmp(param->prmnameS[k],"PHYSRIPL") == 0) {
	      param->physripl_spv = d_alloc_1d(param->ncol);
	      for (col = 0; col < param->ncol; col++)
	      param->prmpointS[k][col] = &param->physripl_spv[col];
	    }
	    if (strcmp(param->prmnameS[k],"MAXBIODEPTH") == 0) {
	      param->maxbiodepth_spv = d_alloc_1d(param->ncol);
	      for (col = 0; col < param->ncol; col++)
	      param->prmpointS[k][col] = &param->maxbiodepth_spv[col];
	    }
	    if (strcmp(param->prmnameS[k],"BIODENS") == 0) {
	      param->biodens_spv = d_alloc_1d(param->ncol);
	      for (col = 0; col < param->ncol; col++)
	      param->prmpointS[k][col] = &param->biodens_spv[col];
	    }
	    if (strcmp(param->prmnameS[k],"BT_PARTIC_KZ") == 0) {
	      param->bt_partic_kz_spv = d_alloc_1d(param->ncol);
	      for (col = 0; col < param->ncol; col++)
	      param->prmpointS[k][col] = &param->bt_partic_kz_spv[col];
	    }
	    if (strcmp(param->prmnameS[k],"BT_PARTIC_KZ_I") == 0) {
	      param->bt_partic_kz_i_spv = d_alloc_1d(param->ncol);
	      for (col = 0; col < param->ncol; col++)
	      param->prmpointS[k][col] = &param->bt_partic_kz_i_spv[col];
	    }
	    if (strcmp(param->prmnameS[k],"BI_DISSOL_KZ") == 0) {
	      param->bi_dissol_kz_spv = d_alloc_1d(param->ncol);
	      for (col = 0; col < param->ncol; col++)
	      param->prmpointS[k][col] = &param->bi_dissol_kz_spv[col];
	    }
	    if (strcmp(param->prmnameS[k],"BI_DISSOL_KZ_I") == 0) {
	      param->bi_dissol_kz_i_spv = d_alloc_1d(param->ncol);
	      for (col = 0; col < param->ncol; col++)
	      param->prmpointS[k][col] = &param->bi_dissol_kz_i_spv[col];
	    }
	    if (strcmp(param->prmnameS[k],"CSS_SCALE") == 0) {
	      param->css_scale_spv = d_alloc_1d(param->ncol);
	      for (col = 0; col < param->ncol; col++)
	      param->prmpointS[k][col] = &param->css_scale_spv[col];
	    }
	  
	}
      }
   }
    // print to check
    if(param->nprmS)
      for (m=0;m<param->nprmS;m++){
	k=param->prmindexS[m]; 
	n=param->trindexS[m];
	// FR commented 03/2014 - if its important it can go in the
	//    runlog or setup.txt
	// fprintf(stderr,"spatial prm %s; tarcer number %d \n",param->prmnameS[k],n);
      }
    
 } //end nov11

 //2012 hardsub
 param->hardsub_numb = -1;
 for(n=0; n < param->ntrB; n++) {
   if ( (strcmp(param->trnameB[n],"reef") == 0) ||  
	(strcmp(param->trnameB[n],"Reef") == 0) ){
     param->hardsub_numb = n;
   }
 }


  if (param->verbose_sed) {
    flog = fopen("sedlog.txt","a");

    fprintf(flog, "\n Tracer attributes \n \n");
    for (n = 0; n < ntr; ++n) {
      sed_tracer_t *tr = &sediment->mstracers[n];
      fprintf(flog, "n=%d, n_hd_wc=%d, n_hd_sed=%d, n_file=%d \n",
        n, tr->n_hd_wc,tr->n_hd_sed, tr->n_file);
      fprintf(flog, "n=%d, tr.name=%s, tr.units=%s \n",
        n, tr->name,tr->units);
      fprintf(flog, "n=%d, tr->fill_value_wc=%f, tr->fill_value_sed=%f \n",
        n, tr->fill_value_wc , tr->fill_value_sed );
      fprintf(flog, "n=%d, tr->diagn=%d, tr->inwc=%d, tr->insed=%d,  \n",
        n, tr->diagn, tr->inwc, tr->insed );
      fprintf(flog, "n=%d, tr->dissol=%d, tr->partic=%d \n",
        n, tr->dissol, tr->partic);
      fprintf(flog, "n=%d, tr->diffuse=%d, tr->decay=%f \n",
        n, tr->diffuse, tr->decay);
      fprintf(flog, "n=%d, tr->psize=%f, tr->b_dens=%f,  tr->i_conc=%f \n",
        n, tr->psize , tr->b_dens,  tr->i_conc );
      if (strlen(param->prmnameS[n]))
	fprintf(flog,"n=%d, tr->svel (spatially varying) = %s\n", n, param->prmnameS[n]);
      else
	fprintf(flog,"n=%d, tr->svel=%f\n", n, tr->svel[0]);
      fprintf(flog,"n=%d, tr->cohesive=%d \n",n,tr->cohesive);
      fprintf(flog,"n=%d, tr->floc=%d \n",n,tr->floc);
      fprintf(flog,"n=%d, tr->resuspend=%d \n",n,tr->resuspend);
      fprintf(flog,"n=%d, tr->deposit=%d \n",n,tr->deposit);
      fprintf(flog, "n=%d, tr->calcvol_wc=%d, tr->calcvol_sed=%d \n",
        n, tr->calcvol_wc , tr->calcvol_sed );
      fprintf(flog,
        "n=%d, tr->adsorb=%d, tr->carriernum=%d, tr->dissolvednum=%d \n",
        n, tr->adsorb , tr->carriernum, tr->dissolvednum);
      fprintf(flog, "n=%d, tr->adsorbkd=%f, tr->adsorbrate=%f \n \n",
        n, tr->adsorbkd , tr->adsorbrate);

    }
    fprintf(flog,"param->calcvol_wc=%d \n",param->calcvol_wc);
    fprintf(flog,"param->calcvol_sed=%d \n \n",param->calcvol_sed);

    fclose(flog);
  }


}

/******************************************************/
/******************************************************/
/*2012*/
 //read benthic variables (number and names)
void sed_tracers_benthic_init(FILE * prmfd, sediment_t *sediment)
{
  int i,k,n,ntrB;
  FILE *flog;
  char buf[MAXSTRLEN], buf1[MAXSTRLEN];
  sed_params_t *param = sediment->msparam;
  sed_tracer_t **mstracers = &sediment->mstracers;
  void *hmodel = sediment->hmodel;
  //number of B tracers
  ntrB = sinterface_getnumberofBtracer(hmodel);
  //B tracer names
  param->trnameB = (char **)p_alloc_1d(ntrB);
  for(n=0; n < ntrB; n++) {
    param->trnameB[n] = malloc(MAXSTRLEN *sizeof(char));
    sinterface_getnameofBtracer(hmodel, n, param->trnameB[n]);
    //    fprintf(stderr,"sed_tracers_benthic_init %d %d %s \n",  ntrB, n, param->trnameB[n]);
    
  }
  param->ntrB = ntrB;  
}

  /**/

/****************************************************************/

sed_column_t *alloc_sed_column(sediment_t *sediment)
{
  sed_params_t *msparam = sediment->msparam;
  sed_column_t *sm = (sed_column_t *)malloc(sizeof(sed_column_t));
  if (sm == NULL) {
      sedtag(LFATAL,"sed:sed_init:alloc_sed_columns"," Could not allocate memory for sediment column.");
      exit(1);
  }
  memset(sm, 0, sizeof(sed_column_t));

  /* allocate memory for sed_column_t arrays, and
     for i,j-dependent variables */

  sm->col_number = 0;

  sm->dz_sed = d_alloc_1d(msparam->sednz-1);
  sm->dz_wc = d_alloc_1d(msparam->nz-1);
  /*UR added 5/2006 */
  sm->dzface_sed = d_alloc_1d(msparam->sednz+1);
  sm->dzface_wc = d_alloc_1d(msparam->nz+1);

  sm->dzold_sed = d_alloc_1d(msparam->sednz-1);
  sm->dzold_wc = d_alloc_1d(msparam->nz-1);
  sm->gridz_sed = d_alloc_1d(msparam->sednz);
  sm->gridz_wc = d_alloc_1d(msparam->nz);
  sm->cellz_sed = d_alloc_1d(msparam->sednz-1);
  sm->cellz_wc = d_alloc_1d(msparam->nz-1);
  sm->por_sed = d_alloc_1d(msparam->sednz-1);
  sm->por_wc = d_alloc_1d(msparam->nz-1);
  sm->porold_sed = d_alloc_1d(msparam->sednz-1);
  sm->porold_wc = d_alloc_1d(msparam->nz-1);
  sm->tr_sed = d_alloc_2d(msparam->sednz-1, msparam->ntr);
  sm->tr_wc = d_alloc_2d(msparam->nz-1, msparam->ntr);
  sm->svel_wc = d_alloc_2d(msparam->nz, msparam->ntr);

  sm->ptr_sed = (double ***) p_alloc_2d(msparam->sednz-1, msparam->ntr);
  sm->ptr_wc = (double ***) p_alloc_2d(msparam->nz-1, msparam->ntr);

  sm->tss_wc = d_alloc_1d(msparam->nz-1);
  sm->tss_sed = d_alloc_1d(msparam->sednz-1);
  sm->sigma_sed = d_alloc_1d(msparam->sednz);

  sm->tmp_sed = d_alloc_1d(msparam->sednz);
  sm->tmp1000a = d_alloc_1d(1001);
  sm->tmp1000b = d_alloc_1d(1001);
  sm->tmp1000c = d_alloc_1d(1001);

  sm->dissol_kz = d_alloc_1d(msparam->sednz);
  sm->partic_kz = d_alloc_1d(msparam->sednz);
  sm->svel_consolid = d_alloc_1d(msparam->sednz);
  sm->svel_floc = d_alloc_1d(msparam->nz);
  sm->gridvel_sed = d_alloc_1d(msparam->sednz);
  sm->gridvel_wc = d_alloc_1d(msparam->nz);
  sm->watvel_wc = d_alloc_1d(msparam->nz);
  sm->watvel_sed = d_alloc_1d(msparam->sednz);

  sm->u1_wc = d_alloc_1d(msparam->nz-1);
  sm->u2_wc = d_alloc_1d(msparam->nz-1);
  sm->Kz_wc = d_alloc_1d(msparam->nz);
  sm->diss_wc = d_alloc_1d(msparam->nz);

  sm->css = d_alloc_1d(msparam->sednz);
  sm->coh_sed = d_alloc_1d(msparam->sednz-1);
  sm->coh_wc = d_alloc_1d(msparam->nz-1);
  sm->depflux = d_alloc_1d(msparam->ntr);
  sm->erdepflux_ac = d_alloc_1d(msparam->ntr);
  sm->erdepflux = d_alloc_1d(msparam->ntr);
  sm->bedloadflux = d_alloc_1d(msparam->ntr);
  sm->erflux = d_alloc_1d(msparam->ntr);
  sm->ref_c_eq = d_alloc_1d(msparam->ntr);
  sm->ref_c = d_alloc_1d(msparam->ntr);
  sm->cbnm1 = d_alloc_1d(msparam->ntr);
  sm->cbnm2 = d_alloc_1d(msparam->ntr);
  sm->cbnm3 = d_alloc_1d(msparam->ntr);
  sm->cbnm4 = d_alloc_1d(msparam->ntr);
  sm->cbfilt = d_alloc_1d(msparam->ntr);
  sm->tr_srf_flux = d_alloc_1d(msparam->ntr);

  return(sm);

}

/**********************************************************/

static void sed_grid_init(FILE * prmfd, sediment_t *sediment, sed_column_t *sm, int write_log)
{
  sed_params_t *param = sediment->msparam;
  int k;
  FILE *flog;
  double dhinter;
  int botk_sed = 0;
  int topk_sed = param->sednz-2;

  /* initialise sediment grid:
     Note that sm->gridz_sed[k] will be updated every (i,j) cycle
     in hyd2sed routine. Here we need it just to specify mindepth_sedi
     and sigma_sed */
  k=sinterface_getdzinit_sed(prmfd, sm->dz_sed, param->sednz-1);
  if ( k+1 != param->sednz) {
      sedtag(LFATAL,"sed:sed_init:sed_grid_init", "Wrong number of sediment layers!");
      exit(1);
  }
  sm->gridz_sed[topk_sed+1] = 0.;
  for(k=topk_sed; k>=botk_sed; k--)
      sm->gridz_sed[k] = sm->gridz_sed[k+1] - sm->dz_sed[k];

  // initialise dzactive
  sm->dzactive = sm->dz_sed[topk_sed]; //2011

  /* Specify minimal thickness of the sediment bed and z0_skin*/
  sm->depth_sedi = sm->gridz_sed[topk_sed+1] - sm->gridz_sed[botk_sed];
  // param->mindepth_sedi = 0.001 * sm->depth_sedi;
  param->mindepth_sedi = 2. * sm->dzactive;//Sep2011
  if( param->mindepth_sedi >=  sm->depth_sedi) {
    fprintf(stderr,"sed_init.c:ERROR:  initial thickness of sediments must exceed the minimum thickness of sediments (which is the double thickness of the top sediment layer) \n");
    exit(1);
  }

  sm->z0_skin = sinterface_getz0(prmfd);

  /* Specify sigma levels for sediment bed below active layer.
     Note that the sigma levels are specified only once and then
     do not change during calculations */
  dhinter = (-sm->gridz_sed[botk_sed] + sm->gridz_sed[topk_sed]);
  for(k=botk_sed;k<=topk_sed;k++)
      sm->sigma_sed[k] = -
    (-sm->gridz_sed[k] + sm->gridz_sed[topk_sed])/dhinter;

  /* print section */
  if (param->verbose_sed && write_log) {
      flog = fopen("sedlog.txt","a");
      for(k=0; k <= topk_sed+1; k++)
	fprintf(flog,"sediment_init k = %d, sm->gridz_sed = %f \n",
		k, sm->gridz_sed[k]);
      fprintf(flog,"sm->z0_skin = %f \n", sm->z0_skin);
      for(k=0;k<=topk_sed;k++)
	fprintf(flog," k = %d, sm->sigma_sed = %f \n", k, sm->sigma_sed[k]);
      fclose(flog);
  }

 
}

static void alloc_sed_spatial(sediment_t *sediment)
{
  sed_params_t *msparam = sediment->msparam;
  void *hmodel = sediment->hmodel;
  int c;

  sed_spatial_t *spatial = (sed_spatial_t *)malloc(sizeof(sed_spatial_t));
  if (spatial == NULL) {
      sedtag(LFATAL,"sed:sed_init:alloc_sed_spatial"," Could not allocate memory for sediment spatial variables.");
      exit(1);
  }
  memset(spatial, 0, sizeof(sed_spatial_t));
  sediment->spatial = spatial;

  spatial->theta = d_alloc_1d(msparam->ncol);
  spatial->hripples = d_alloc_1d(msparam->ncol);
  spatial->lripples = d_alloc_1d(msparam->ncol);
  spatial->cbnm1 = d_alloc_2d(msparam->ncol, msparam->ntr);
  spatial->cbnm2 = d_alloc_2d(msparam->ncol, msparam->ntr);
  spatial->cbnm3 = d_alloc_2d(msparam->ncol, msparam->ntr);
  spatial->cbnm4 = d_alloc_2d(msparam->ncol, msparam->ntr);
  spatial->erdeprate = d_alloc_2d(msparam->ncol, msparam->ntr);

  // FR 06-2015: Deferred to hd2sed
  /*  
      for (c=0;c<msparam->ncol;c++)
      {
      spatial->hripples[c] = max(msparam->physriph, msparam->bioriph);
      spatial->lripples[c] = max(msparam->physripl, msparam->bioripl);
      }
  */
  sinterface_gettheta(hmodel, spatial->theta, msparam->ncol);
}


static void sed_optimization(sediment_t *sediment)
{
  int n;
  sed_params_t *param = sediment->msparam;


  if(param->flocmode == 1 && strcmp("salt",(&sediment->mstracers[0])->name) != 0) {
    sedtag(LFATAL,"sed:sed_init:sed_optimization","First tracer name should be salt when floculation mode = 1 - %s",(&sediment->mstracers[0])->name);
    exit(1);
  }

  param->vtransp_por_sed = calloc(param->ntr,sizeof(int));
  param->n_vtransp_por_sed = 0;

  param->vtransp_por_wc = calloc(param->ntr,sizeof(int));
  param->n_vtransp_por_wc = 0;

  param->vtransp_decay = calloc(param->ntr,sizeof(int));
  param->n_vtransp_decay = 0;

  param->vtransp_adsorb = calloc(param->ntr,sizeof(int));
  param->n_vtransp_adsorb = 0;

  param->vtransp_floc = calloc(param->ntr,sizeof(int));
  param->n_vtransp_floc = 0;

  param->vtransp_dissolv = calloc(param->ntr,sizeof(int));
  param->n_vtransp_dissolv = 0;

  param->vtransp_partic = calloc(param->ntr,sizeof(int));
  param->n_vtransp_partic = 0;

  param->vtransp_insed_partic = calloc(param->ntr,sizeof(int));
  param->n_vtransp_insed_partic = 0;

  param->vtransp_nd_insed_partic = calloc(param->ntr,sizeof(int));
  param->n_vtransp_nd_insed_partic = 0;

  /* establish at this time which tracers should have porosity in sediment calculated */
  for (n = 0; n < param->ntr; n++) {
    sed_tracer_t *tracer = &sediment->mstracers[n];
    if (!tracer->diagn && tracer->insed && tracer->partic && !(tracer->adsorb) && tracer->calcvol_sed)
    {
      param->vtransp_por_sed[param->n_vtransp_por_sed] = n;
      param->n_vtransp_por_sed++;
    }

  /* establish at this time which tracers should have porosity in the water column calculated */
    if (!tracer->diagn && tracer->insed && tracer->partic && !(tracer->adsorb) && tracer->calcvol_wc)
    {
      param->vtransp_por_wc[param->n_vtransp_por_wc] = n;
      param->n_vtransp_por_wc++;
    }

    if (tracer->insed && tracer->partic && !tracer->adsorb)
    {
      param->vtransp_insed_partic[param->n_vtransp_insed_partic] = n;
      param->n_vtransp_insed_partic++;
    }

    if (!tracer->diagn && tracer->insed && tracer->partic && !tracer->adsorb)
    {
      param->vtransp_nd_insed_partic[param->n_vtransp_nd_insed_partic] = n;
      param->n_vtransp_nd_insed_partic++;
    }


    if (tracer->partic && !tracer->adsorb)
    {
      param->vtransp_partic[param->n_vtransp_partic] = n;
      param->n_vtransp_partic++;
    }

    if (!tracer->diagn && tracer->insed && tracer->adsorb)
    {
      param->vtransp_adsorb[param->n_vtransp_adsorb] = n;
      param->n_vtransp_adsorb++;      
    }

    if (!tracer->diagn && tracer->decay)
    {
      param->vtransp_decay[param->n_vtransp_decay] = n;
      param->n_vtransp_decay++;
    }

    if (tracer->partic && !(tracer->adsorb) && tracer->floc)
    {
      param->vtransp_floc[param->n_vtransp_floc] = n;
      param->n_vtransp_floc++;
    }

    if (!tracer->diagn && tracer->insed && tracer->dissol)
    {
      param->vtransp_dissolv[param->n_vtransp_dissolv] = n;
      param->n_vtransp_dissolv++;
      /* one wc layer */
    }

  }
}

/*-------------------------------------------------------------------*/
/* Routine to read the sediment parameters                           */
/*-------------------------------------------------------------------*/
FILE *get_sed_params(FILE *fp, sediment_t *sediment)
{
  int i;
  char fname[MAXSTRLEN];
  FILE *rp;
  struct {
    char *name;
    char *description;
    void (*init) (sediment_t *sediment);
  } param_list[] = {
    {"standard","Standard sediment values",  sed_params_std},
    {"estuary", "Estuarine sediment values", sed_params_est},
    {"shelf", "Shelf sediment values", sed_params_shf},
    {"basic", "Basic sediment values", sed_params_bsc},
    {NULL, NULL, NULL}
  };
  void (*init) (sediment_t *sediment)= NULL;

  if (prm_read_char(fp, "SEDFILE", fname)) {
    if ((rp = fopen(fname, "r"))) {
      return(rp);
    } else if (strcmp(fname, "auto") == 0) {
      sed_params_build(sediment);
      sed_params_auto(sediment);
      return(NULL);
    } else {
      for (i = 0; (init == NULL) && param_list[i].name; ++i) {
	if (strcasecmp(fname, param_list[i].name) == 0) {
	  init = param_list[i].init;
	}
      }
      if (init != NULL) {
	sed_params_build(sediment);
	init(sediment);
      }
      return(NULL);
    }
  } else
    return(fp);
}

/* END get_sed_params()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Allocates and initialises the sediment parameter structure        */
/*-------------------------------------------------------------------*/
void sed_params_build(sediment_t *sediment)
{
  void *hmodel = sediment->hmodel;

  /* Allocate memory for sed_params_t Structure */
  sed_params_t *param = (sed_params_t *)malloc(sizeof(sed_params_t));
  if (param == NULL) {
      sedtag(LFATAL,"sed:sed_init:sed_tracers_init"," Could not allocate memory for a sediment parameters.");
      exit(1);
  }
  memset(param, 0, sizeof(sed_params_t));
  sediment->msparam = param;

  param->t=sinterface_getmodeltime(hmodel);
  param->dt =  sinterface_getmodeltimestep(hmodel);
  param->nz = sinterface_getnumberwclayers(hmodel)+1;
  /* If there is only one water layer, mecosed will split it in two layers*/
  if(param->nz == 2)
      param->nz = 3;
  param->sednz = sinterface_getnumbersedlayers(hmodel)+1;
  /* Note that in mecosed param->sednz and param->nz are defined as
     number of interfaces (not sediment layers) */
  if(param->sednz < 3) {
      sedtag(LFATAL,"sed:sed_init:sed_tracers_init","At least two sediment layers must be specified");
      exit(1);
  }

  param->nstep=0;
  param->ncol =  sinterface_getnumbercolumns(hmodel);
  param->ntr = sinterface_getnumberoftracers(hmodel);

  sinterface_getmap2diagtracer_2d(hmodel,
				  &param->n_hripple,  &param->n_lripple,
				  &param->n_ustrcw_skin,
				  &param->n_depth_sedi,  &param->n_dzactive,
				  &param->n_erdepflux_total,
				  &param->n_erdepflux_total_ac);
  si_set_errfn_warn(hmodel);
}

/* END sed_params_build()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Default sediment parameter values. These mirror the values        */
/* assigned in sediments/sediments.c as defaults.                    */
/*-------------------------------------------------------------------*/
void sed_params_est(sediment_t *sediment)
{
  sed_params_t *param = sediment->msparam;

  param->verbose_sed = 0;
  param->quad_bfc = 0.00003;
  param->geomorph = 0;
  param->consolidate =  0;
  param->minpor_wc = 0.1;
  param->minpor_sed = 0.1;
  param->minseddz = 0.000001;
  param->finpor_sed = 0.4;
  param->consolrate = 10.;
  param->cssmode = 4;
  param->css = 0.2;
  param->css_er_val = d_alloc_1d(param->sednz);
  param->css_er_val[0] = 0.7;
  param->css_er_val[1] = 0.6;
  param->css_er_val[2] = 0.4;
  param->css_er_val[3] = 0.2;
  param->css_er_depth = d_alloc_1d(param->sednz);
  param->css_er_depth[0] = -0.080;
  param->css_er_depth[1] = -0.020;
  param->css_er_depth[2] = -0.005;
  param->css_er_depth[3] = -0.0;
  param->css_dep = 1.e+13;
  param->flocmode = 4;
  param->flocprm1 = 2.e-4;
  param->flocprm2 = 3.0;
  param->bbl_nonlinear = 1;
  param->calc_ripples = 0;
  param->physriph = 0.02;
  param->bioriph = 0.0;
  param->physripl = 0.5;
  param->bioripl = 0.5;
  param->biodens= 30.0;
  param->maxbiodepth= 0.2;
  param->bi_dissol_kz= 1.e-9;
  param->bt_partic_kz= 1.e-10;
  param->bi_dissol_kz_i= 1.e-13;
  param->bt_partic_kz_i= 2.e-13;
  param->biosedprofile = 'p';
  param->max_thick_sed = 10;
}

/* END sed_params_est()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Standard sediment parameter set.                                  */
/*-------------------------------------------------------------------*/
void sed_params_std(sediment_t *sediment)
{
  sed_params_t *param = sediment->msparam;

  param->verbose_sed = 0;
  param->quad_bfc = 0.00003;
  param->geomorph = 0;
  param->consolidate =  0;
  param->minpor_wc = 0.1;
  param->minpor_sed = 0.1;
  param->minseddz = 0.000001;
  param->finpor_sed = 0.4;
  param->consolrate =  10.0;
  param->cssmode = 4;
  param->css = 0.2;
  param->css_er_val = d_alloc_1d(param->sednz);
  param->css_er_val[0] = 0.7;
  param->css_er_val[1] = 0.6;
  param->css_er_val[2] = 0.4;
  param->css_er_val[3] = 0.2;
  param->css_er_depth = d_alloc_1d(param->sednz);
  param->css_er_depth[0] = -0.080;
  param->css_er_depth[1] = -0.020;
  param->css_er_depth[2] = -0.005;
  param->css_er_depth[3] = -0.0;
  param->css_dep = 1.e+13;
  param->flocmode = 4;
  param->flocprm1 = 2.e-4;
  param->flocprm2 = 3.0;
  param->bbl_nonlinear = 1;
  param->calc_ripples = 0;
  param->physriph = 0.01;
  param->bioriph = 0.0;
  param->physripl = 0.5;
  param->bioripl = 0.03;
  param->biodens= 30.0;
  param->maxbiodepth= 0.2;
  param->bi_dissol_kz= 1.0e-9;
  param->bt_partic_kz= 1.0e-10;
  param->bi_dissol_kz_i= 1.0e-13;
  param->bt_partic_kz_i= 2.0e-13;
  param->biosedprofile = 'p';
  param->max_thick_sed = 10;
}

/* END sed_params_std()                                              */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Shelf sediment parameter set.                                  */
/*-------------------------------------------------------------------*/
void sed_params_shf(sediment_t *sediment)
{
  sed_params_t *param = sediment->msparam;

  param->verbose_sed = 0;
  param->quad_bfc = 0.00003;
  param->geomorph = 0;
  param->consolidate =  0;
  param->minpor_wc = 0.1;
  param->minpor_sed = 0.1;
  param->minseddz = 0.000001;
  param->finpor_sed = 0.4;
  param->consolrate =  10.0;
  param->cssmode = 4;
  param->css = 0.2;
  param->css_er_val = d_alloc_1d(param->sednz);
  param->css_er_val[0] = 0.7;
  param->css_er_val[1] = 0.6;
  param->css_er_val[2] = 0.4;
  param->css_er_val[3] = 0.2;
  param->css_er_depth = d_alloc_1d(param->sednz);
  param->css_er_depth[0] = -0.080;
  param->css_er_depth[1] = -0.020;
  param->css_er_depth[2] = -0.005;
  param->css_er_depth[3] = -0.0;
  param->css_dep = 1.e+13;
  param->flocmode = 4;
  param->flocprm1 = 2.e-4;
  param->flocprm2 = 3.0;
  param->bbl_nonlinear = 1;
  param->calc_ripples = 0;
  param->physriph = 0.001;
  param->bioriph = 0.0;
  param->physripl = 0.5;
  param->bioripl = 0.03;
  param->biodens= 30.0;
  param->maxbiodepth= 0.2;
  param->bi_dissol_kz= 1.0e-9;
  param->bt_partic_kz= 1.0e-10;
  param->bi_dissol_kz_i= 1.0e-13;
  param->bt_partic_kz_i= 2.0e-13;
  param->biosedprofile = 'p';
  param->max_thick_sed = 10;
}

/* END sed_params_shf()                                              */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Basic sediment parameter set.                                  */
/*-------------------------------------------------------------------*/
void sed_params_bsc(sediment_t *sediment)
{
  sed_params_t *param = sediment->msparam;

  param->verbose_sed = 0;
  param->quad_bfc = 0.00003;
  param->geomorph = 0;
  param->consolidate =  0;
  param->minpor_wc = 0.1;
  param->minpor_sed = 0.1;
  param->minseddz = 0.000001;
  param->finpor_sed = 0.4;
  param->consolrate =  10.0;
  param->cssmode = 4;
  param->css = 0.2;
  param->css_er_val = d_alloc_1d(param->sednz);
  param->css_er_val[0] = 0.2;
  param->css_er_val[1] = 0.2;
  param->css_er_val[2] = 0.2;
  param->css_er_val[3] = 0.2;
  param->css_er_depth = d_alloc_1d(param->sednz);
  param->css_er_depth[0] = -0.080;
  param->css_er_depth[1] = -0.020;
  param->css_er_depth[2] = -0.005;
  param->css_er_depth[3] = -0.0;
  param->css_dep = 1.e+13;
  param->flocmode = 0;
  param->flocprm1 = 2.e-4;
  param->flocprm2 = 3.0;
  param->bbl_nonlinear = 1;
  param->calc_ripples = 0;
  param->physriph = 0.01;
  param->bioriph = 0.0;
  param->physripl = 0.5;
  param->bioripl = 0.03;
  param->biodens= 1.0;
  param->maxbiodepth= 0.2;
  param->bi_dissol_kz= 0.0;
  param->bt_partic_kz= 0.0;
  param->bi_dissol_kz_i= 0.0;
  param->bt_partic_kz_i= 0.0;
  param->biosedprofile = 'p';
  param->max_thick_sed = 100;
}

/* END sed_params_basic()                                              */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Dynamically prescribed sediment parameter specification           */
/*-------------------------------------------------------------------*/
void sed_params_auto(sediment_t *sediment)
{
  void *hmodel = sediment->hmodel;
  sed_params_t *param = sediment->msparam;

  /* Set defaults */
  sed_params_std(sediment);

  /* Prescribe values based on known system attributes */
}

/* END sed_params_auto()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Writes the sediment parameters to logfile                         */
/*-------------------------------------------------------------------*/
void sed_params_write(sediment_t *sediment)
{
  sed_params_t *param = sediment->msparam;
  FILE *flog;
  int k;

  if (param->verbose_sed) {
    if ((flog = fopen("sedlog.txt","w")) == NULL)
      return;

    fprintf(flog,"sed_init \n\n");
    fprintf(flog,"param->t=%f\n", param->t);
    fprintf(flog,"param->dt=%f\n", param->dt);
    fprintf(flog,"param->nz=%d\n", param->nz);
    fprintf(flog,"param->sednz=%d\n", param->sednz);
    fprintf(flog,"param->ncol=%d\n", param->ncol);
    fprintf(flog,"param->nstep=%d\n", param->nstep);
    fprintf(flog,"param->ntr=%d\n", param->ntr);
    fprintf(flog,"param->quad_bfc=%f\n", param->quad_bfc);
    fprintf(flog,"GEOMPRPH =%d\n", param->geomorph);
    fprintf(flog,"CONSOLIDATE =%d\n", param->consolidate);
    fprintf(flog,"MINPOR_SED =%f\n", param->minpor_sed);
    fprintf(flog,"MINPOR_WC =%f\n", param->minpor_wc);
    fprintf(flog,"FINPOR_SED =%f\n", param->finpor_sed);
    fprintf(flog,"CONSOLRATE =%f\n", param->consolrate);
    fprintf(flog,"CSSERMODE = %d \n", param->cssmode);
    fprintf(flog,"CSS =%f \n", param->css );
    fprintf(flog,"CSSDEP=%f \n", param->css_dep );
    fprintf(flog,"FLOC_MODE = %d \n", param->flocmode);
    fprintf(flog,"FLOCPRM1=%f \n", param->flocprm1 );
    fprintf(flog,"FLOCPRM2=%f \n", param->flocprm2 );
    fprintf(flog,"BBL_NONLINEAR=%d \n", param->bbl_nonlinear );
    fprintf(flog,"CALC_RIPPLES=%d \n", param->calc_ripples);
    if (param->physriph_spv)
      fprintf(flog,"PHYSRIPH spatially varying\n");
    else
      fprintf(flog,"PHYSRIPH=%f \n", param->physriph);
    if (param->physripl_spv)
      fprintf(flog,"PHYSRIPL=spatially varying\n");
    else
      fprintf(flog,"PHYSRIPL=%f \n", param->physripl);
    fprintf(flog,"BIORIPH=%f \n", param->bioriph );
    fprintf(flog,"BIORIPL=%f \n", param->bioripl );
    if (param->biodens_spv)
      fprintf(flog,"BIODENS spatially varying\n");
    else
      fprintf(flog,"BIODENS =%f\n", param->biodens);
    if (param->maxbiodepth_spv)
      fprintf(flog,"MAXBIODEPTH spatially varying\n");
    else
      fprintf(flog,"MAXBIODEPTH =%f\n", param->maxbiodepth);
    if (param->bi_dissol_kz_spv)
      fprintf(flog,"BI_DISSOL_KZ spatially varying\n");
    else
      fprintf(flog,"BI_DISSOL_KZ =%e\n", param->bi_dissol_kz);
    if (param->bi_dissol_kz_i_spv)
      fprintf(flog,"BI_DISSOL_I_KZ spatially varying\n");
    else
      fprintf(flog,"BI_DISSOL_I_KZ =%e\n", param->bi_dissol_kz_i);
    if (param->bt_partic_kz_spv)
      fprintf(flog,"BT_PARTIC_KZ spatially varying\n");
    else
      fprintf(flog,"BT_PARTIC_KZ =%e\n", param->bt_partic_kz);
    if (param->bt_partic_kz_i_spv)
      fprintf(flog,"BT_PARTIC_KZ_I spatially varying\n");
    else
      fprintf(flog,"BT_PARTIC_KZ_I =%e\n", param->bt_partic_kz_i);
    fprintf(flog,"BIOSEDPROFILE =%c\n", param->biosedprofile );
    if (param->css_scale_spv)
      fprintf(flog,"CSS_SCALE spatially varying\n");
    else
      fprintf(flog,"CSS_SCALE =%e\n", param->css_scale);

    if (param->cssmode == 4) {
      for(k=param->sednz-1;k>=0;k--)
	fprintf(flog, "er_depth=%f, er_val=%f \n", param->css_er_depth[k],param->css_er_val[k]);
    }
    fclose(flog);
  }
}

/* END sed_params_write()                                            */
/*-------------------------------------------------------------------*/

#if AS_EMSLIB_MODULE

/* prototypes - defined in emslogger.c */
int errfn_log_default (int level, char *format, ...);
int errfn_tag_default (int level, char *tag, char *format, ...);

sedlogfn  sedlog = errfn_log_default;
sedlogtag sedtag = errfn_tag_default;

#else
/* log redirection */

char* sedlevel2nice(int t)
{
  if(t == LERROR)
    return "ERROR ";
  if(t == LFATAL)
    return "FATAL ";
  if(t == LPANIC)
    return "PANIC ";
  if(t == LDEBUG)
    return "DEBUG ";
  if(t == LMETRIC)
    return "METRIC";
  if(t == LTRACE)
    return "TRACE ";
  if(t == LMAIN)
    return "MAIN  ";
  if(t == LINFO)
    return "INFO  ";

  return   "  NA  ";
}


int sedwritelog(int level, const char *tag, const char *s,va_list args)
{
  /* Check that the log flag is set.
   */

    char buffer[MAXSTRLEN];
    time_t now = time(NULL);
    struct tm *t = localtime(&now);

    vsprintf(buffer, s, args);
    if (strlen(buffer) > 0) {
      printf("[%d/%2.2d/%2.2d %2.2d:%2.2d:%2.2d][%s](%s) %s\n",
              1900 + t->tm_year, t->tm_mon, t->tm_mday,
              t->tm_hour, t->tm_min, t->tm_sec,sedlevel2nice(level),tag, buffer);
    }
    return 1;
}


int sedlogimpl(int level, char *format, ...)
{
  int i = 0;
  va_list args;
  va_start(args, format);
  i =sedwritelog(level,"",format,args);
  va_end(args);
  return i;
}

int sedtagimpl(int level, char *tag, char *format, ...)
{
  int i = 0;
  va_list args;
  va_start(args, format);
  i =sedwritelog(level,tag,format,args);
  va_end(args);
  return i;
}


sedlogfn  sedlog = sedlogimpl;
sedlogtag sedtag = sedtagimpl;

#endif


#endif













/**
        * Ansi C "itoa" based on Kernighan & Ritchie's "Ansi C":
	 */
	void strreverse(char* begin, char* end) {
		char aux;
		while(end>begin)
			aux=*end, *end--=*begin, *begin++=aux;
	}
	char* _itoa(int value, char* str, int base) {
		static char num[] = "0123456789abcdefghijklmnopqrstuvwxyz";
		char* wstr=str;
		int sign;
	
		// Validate base
		if (base<2 || base>35){ *wstr='\0'; return 0; }
	
		// Take care of sign
		if ((sign=value) < 0) value = -value;
	
		// Conversion. Number is reversed.
		do *wstr++ = num[value%base]; while(value/=base);
		if(sign<0) *wstr++='-';
		*wstr='\0';
	
		// Reverse string
		strreverse(str,wstr-1);
		return(str);
	}
	







#ifdef __cplusplus
}
#endif







