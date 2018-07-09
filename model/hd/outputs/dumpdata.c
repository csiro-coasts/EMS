/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/outputs/dumpdata.c
 *  
 *  Description:
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: dumpdata.c 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

#include <stdlib.h>
#include <string.h>
#include "hd.h"
#include "tracer.h"

void sigma_vmap(master_t *master, dump_data_t *dumpdata);
static int contains_string(char *p, char *tag);

/*-------------------------------------------------------------------*/
/* Routine to allocate memory for the parameter data structure       */
/*-------------------------------------------------------------------*/
dump_data_t *create_dumpdata(void)
{
  dump_data_t *dumpdata = (dump_data_t *)malloc(sizeof(dump_data_t));
  memset(dumpdata, 0, sizeof(dump_data_t));
  return dumpdata;
}

/* END create_dumpdata()                                             */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Routine to populate the dump structure with 2d data               */
/*-------------------------------------------------------------------*/
dump_data_t *dumpdata_build(parameters_t *params, /* Input parameter data
                                                     structure */
                            geometry_t *geom, /* Sparse global geometry
                                                 structure */
                            master_t *master, /* Master data structure */
                            unsigned long ***flag,  /* Flags array */
                            double *layers, /* Vertical grid spacing */
                            double **bathy  /* Bathymetry */
  )
{
  dump_data_t *dumpdata;
  int nce1 = geom->nce1;
  int nce2 = geom->nce2;
  int nfe1 = geom->nfe1;
  int nfe2 = geom->nfe2;
  int nz = geom->nz;
  int i, j, k;                  /* Cartesian indexes */
  int sednz = geom->sednz;
  char buf[MAXSTRLEN];

  dumpdata = create_dumpdata();
  master->dumpdata = dumpdata;
  dumpdata->momgrid = NULL;
  dumpdata->romsgrid = NULL;
  dumpdata->large = 1;
  strcpy(buf, params->oname);
  if (contains_string(buf, "short"))
    dumpdata->large = 0;

  dumpdata->nz = nz;
  dumpdata->sednz = sednz;
  dumpdata->nce1 = nce1;
  dumpdata->nce2 = nce2;
  dumpdata->nfe1 = nfe1;
  dumpdata->nfe2 = nfe2;
  dumpdata->ntr = min(master->ntr, 20);
  dumpdata->ntr = master->ntr;
  dumpdata->ntrS = master->ntrS;
  dumpdata->nsed = master->nsed;
  dumpdata->trinfo_3d = master->trinfo_3d;
  dumpdata->trinfo_2d = master->trinfo_2d;
  dumpdata->trinfo_sed = master->trinfo_sed;
  dumpdata->ns2 = params->ns2;
  dumpdata->ns3 = params->ns3;
  dumpdata->tmode = params->tmode;
  dumpdata->togn = master->togn;
  strcpy(dumpdata->rev, params->rev);
  dumpdata->runno = params->runno;

  strcpy(dumpdata->lenunit, params->lenunit);
  strcpy(dumpdata->prmname, params->prmname);
  strcpy(dumpdata->timeunit, params->timeunit);
  strcpy(dumpdata->output_tunit, params->output_tunit);
  strcpy(dumpdata->idumpname, params->idumpname);
  strcpy(dumpdata->gridtype, params->gridtype);
  if (strcmp(params->gridtype, "RECTANGULAR") == 0) {
    dumpdata->rg = (grid_rect_t *)malloc(sizeof(grid_rect_t));
    dumpdata->rg->x0 = params->rg->x0;
    dumpdata->rg->y0 = params->rg->y0;
    dumpdata->rg->dx = params->rg->dx;
    dumpdata->rg->dy = params->rg->dy;
    dumpdata->rg->th = params->rg->th;
    dumpdata->rg->sinth = params->rg->sinth;
    dumpdata->rg->costh = params->rg->costh;
  }
  if (strcmp(params->gridtype, "POLAR") == 0) {
    dumpdata->pg = (grid_polar_t *)malloc(sizeof(grid_polar_t));
    dumpdata->pg->x0 = params->pg->x0;
    dumpdata->pg->y0 = params->pg->y0;
    dumpdata->pg->arc = params->pg->arc;
    dumpdata->pg->rmin = params->pg->rmin;
    dumpdata->pg->rotation = params->pg->rotation;
  }

  /* Allocate memory for all dump variables */
  /* Variables dumped when the dumpfile is created */
  dumpdata->gridz = d_alloc_1d(nz + 1);
  dumpdata->cellz = d_alloc_1d(nz);
  dumpdata->gridx = d_alloc_2d(nce1 + 1, nce2 + 1);
  dumpdata->gridy = d_alloc_2d(nce1 + 1, nce2 + 1);
  dumpdata->cellx = d_alloc_2d(nce1, nce2);
  dumpdata->celly = d_alloc_2d(nce1, nce2);
  if (sednz) {
    dumpdata->gridz_sed = d_alloc_3d(nce1, nce2, sednz + 1);
    dumpdata->cellz_sed = d_alloc_3d(nce1, nce2, sednz);
    dumpdata->cellzcsed = d_alloc_1d(sednz);
  }
  dumpdata->u1x = d_alloc_2d(nce1 + 1, nce2);
  dumpdata->u1y = d_alloc_2d(nce1 + 1, nce2);
  dumpdata->u2x = d_alloc_2d(nce1, nce2 + 1);
  dumpdata->u2y = d_alloc_2d(nce1, nce2 + 1);
  dumpdata->botz = d_alloc_2d(nce1, nce2);
  dumpdata->h1au1 = d_alloc_2d(nce1 + 1, nce2);
  dumpdata->h1au2 = d_alloc_2d(nce1, nce2 + 1);
  dumpdata->h1acell = d_alloc_2d(nce1, nce2);
  dumpdata->h1agrid = d_alloc_2d(nce1 + 1, nce2 + 1);
  dumpdata->thetau1 = d_alloc_2d(nce1 + 1, nce2);
  dumpdata->h2au1 = d_alloc_2d(nce1 + 1, nce2);
  dumpdata->h2au2 = d_alloc_2d(nce1, nce2 + 1);
  dumpdata->h2acell = d_alloc_2d(nce1, nce2);
  dumpdata->h2agrid = d_alloc_2d(nce1 + 1, nce2 + 1);
  dumpdata->thetau2 = d_alloc_2d(nce1, nce2 + 1);
  dumpdata->coriolis = d_alloc_2d(nce1, nce2);
  dumpdata->flag = (unsigned long ***)l_alloc_3d(nfe1, nfe2, nz);

  /* Time dependent dump variables */
  dumpdata->u1av = d_alloc_2d(nfe1, nce2);
  dumpdata->u2av = d_alloc_2d(nce1, nfe2);
  dumpdata->wind1 = d_alloc_2d(nfe1, nce2);
  dumpdata->wind2 = d_alloc_2d(nce1, nfe2);
  dumpdata->wtop = d_alloc_2d(nce1, nce2);
  dumpdata->topz = d_alloc_2d(nce1, nce2);
  dumpdata->eta = d_alloc_2d(nce1, nce2);
  dumpdata->patm = d_alloc_2d(nce1, nce2);
  dumpdata->u1 = d_alloc_3d(nfe1, nce2, nz);
  dumpdata->u2 = d_alloc_3d(nce1, nfe2, nz);
  dumpdata->w = d_alloc_3d(nce1, nce2, nz);
  /* 
   * Optional velocity means, volume fluxes and streamline origin and
   * courant numbers
   */
  if (master->u1m)
    dumpdata->u1m = d_alloc_3d(nfe1, nce2, nz);
  if (master->u2m)
    dumpdata->u2m = d_alloc_3d(nce1, nfe2, nz);
  if (master->u1vm)
    dumpdata->u1vm = d_alloc_3d(nfe1, nce2, nz);
  if (master->u2vm)
    dumpdata->u2vm = d_alloc_3d(nce1, nfe2, nz);
  if (master->origin)
    dumpdata->origin = d_alloc_3d(nce1, nce2, nz);
  if (master->pc)
    dumpdata->pc = d_alloc_3d(nce1, nce2, nz);
  if (master->qc)
    dumpdata->qc = d_alloc_3d(nce1, nce2, nz);
  if (master->rc)
    dumpdata->rc = d_alloc_3d(nce1, nce2, nz);
  
  dumpdata->dzu1 = d_alloc_3d(nfe1, nce2, nz);
  dumpdata->dzu2 = d_alloc_3d(nce1, nfe2, nz);
  dumpdata->dzcell = d_alloc_3d(nfe1, nfe2, nz);
  dumpdata->u1vh = d_alloc_2d(nce1, nce2);
  dumpdata->u2vh = d_alloc_2d(nce1, nce2);
  dumpdata->u1bot = d_alloc_2d(nfe1, nce2);
  dumpdata->u2bot = d_alloc_2d(nce1, nfe2);
  dumpdata->tr_wc = d_alloc_4d(nce1, nce2, nz, dumpdata->ntr);
  if (master->ntrS)
    dumpdata->tr_wcS = d_alloc_3d(nce1, nce2, master->ntrS);
  if (sednz && master->nsed)
    dumpdata->tr_sed = d_alloc_4d(nce1, nce2, geom->sednz, master->nsed);

  /* Time dependent diagnostic dump variables */
  dumpdata->dens = d_alloc_3d(nce1, nce2, nz);
  dumpdata->dens_0 = d_alloc_3d(nce1, nce2, nz);
  dumpdata->Kz = d_alloc_3d(nce1, nce2, nz);
  dumpdata->Vz = d_alloc_3d(nce1, nce2, nz);
  dumpdata->Cd = d_alloc_2d(nce1, nce2);

  /* Cell mapping */
  dumpdata->crci = s_alloc_1d(nce1);
  dumpdata->clci = s_alloc_1d(nce1);
  dumpdata->frci = s_alloc_1d(nce1);
  dumpdata->flci = s_alloc_1d(nce1);
  dumpdata->crfi = s_alloc_1d(nfe1);
  dumpdata->clfi = s_alloc_1d(nfe1);
  dumpdata->frfi = s_alloc_1d(nfe1);
  dumpdata->flfi = s_alloc_1d(nfe1);
  dumpdata->cfcj = s_alloc_1d(nce2);
  dumpdata->cbcj = s_alloc_1d(nce2);
  dumpdata->ffcj = s_alloc_1d(nce2);
  dumpdata->fbcj = s_alloc_1d(nce2);
  dumpdata->cffj = s_alloc_1d(nfe2);
  dumpdata->cbfj = s_alloc_1d(nfe2);
  dumpdata->fffj = s_alloc_1d(nfe2);
  dumpdata->fbfj = s_alloc_1d(nfe2);

  /* Dummies, maps */
  dumpdata->w1 = d_alloc_1d(geom->sgsiz);
  dumpdata->w2 = d_alloc_2d(nfe1, nfe2);
  dumpdata->w3 = d_alloc_3d(nfe1, nfe2, nz);
  dumpdata->i1 = geom->wsa;
  /* Extra sparse array used to pack arrays in df_sparse */
  dumpdata->w1s = d_alloc_1d(geom->sgsiz);

  /* Layer geometry and flags */
  for (k = 0; k < nz; k++) {
    dumpdata->gridz[k] = layers[k];
    dumpdata->cellz[k] = 0.5 * (layers[k] + layers[k + 1]);
    for (j = 0; j < nfe2; j++)
      for (i = 0; i < nfe1; i++)
        dumpdata->flag[k][j][i] = (params->sigma) ? flag[nz - 1][j][i] : 
	  flag[k][j][i];
  }
  dumpdata->gridz[nz] = MAXGRIDZ;

  /* Sediments */
  if (sednz) {
    for (j = 0; j < nce2; j++)
      for (i = 0; i < nce1; i++) {
        for (k = 0; k < sednz; k++) {
	  dumpdata->gridz_sed[k][j][i] = params->gridz_sed[k];
	  dumpdata->cellz_sed[k][j][i] = 0.5 * (params->gridz_sed[k] +
						  params->gridz_sed[k + 1]);
	}
	dumpdata->gridz_sed[sednz][j][i] = params->gridz_sed[sednz];
      }
    /* This is the nominal z centre in sediments */
    for (k = 0; k < sednz; k++)
      dumpdata->cellzcsed[k] = dumpdata->cellz_sed[k][0][0];
  }

  /* Bathymetry is overwritten if dump_read() is invoked */
  for (j = 0; j < nce2; j++)
    for (i = 0; i < nce1; i++) {
      if (flag[nz - 1][j][i] & OUTSIDE)
        dumpdata->botz[j][i] = NaN;
      else
        dumpdata->botz[j][i] = bathy[j][i];
    }

  /* Get the mom/roms grid structure if required */
  dumpdata->momgrid  = convert_mom_grid(params, dumpdata);
  dumpdata->romsgrid = convert_roms_grid(params, dumpdata);
  dumpdata->fmap_i = dumpdata->fmap_j = NULL;

  /* Get the vertical mapping function */
  dumpdata->vmap = NULL;
  if (params->sigma)
    sigma_vmap(master, dumpdata);    

  /* Copy the output path */
  strcpy(dumpdata->opath, params->opath);
  /* Copy the transport key */
  strcpy(dumpdata->trkey, params->trkey);

  dumpdata->master = master;
  return (dumpdata);
}

/* END dumpdata_build()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to populate the dump structure with 2d data               */
/*-------------------------------------------------------------------*/
void dumpdata_init(dump_data_t *dumpdata, /* Dump data structure */
                   geometry_t *geom,  /* Sparse global geometry structure */
                   master_t *master /* Master data structure */
  )
{
  int nce1 = dumpdata->nce1;
  int nce2 = dumpdata->nce2;
  int nfe1 = dumpdata->nfe1;
  int nfe2 = dumpdata->nfe2;
  dumpdata->t = master->t;

  /* Fill the initialisation variables with sparse data */
  s2c_2d(geom, geom->gridx, dumpdata->gridx, nfe1, nfe2);
  s2c_2d(geom, geom->gridy, dumpdata->gridy, nfe1, nfe2);
  s2c_2d(geom, geom->cellx, dumpdata->cellx, nce1, nce2);
  s2c_2d(geom, geom->celly, dumpdata->celly, nce1, nce2);
  s2c_2d(geom, geom->u1x, dumpdata->u1x, nfe1, nce2);
  s2c_2d(geom, geom->u1y, dumpdata->u1y, nfe1, nce2);
  s2c_2d(geom, geom->u2x, dumpdata->u2x, nce1, nfe2);
  s2c_2d(geom, geom->u2y, dumpdata->u2y, nce1, nfe2);
  s2c_2d(geom, geom->botz, dumpdata->botz, nce1, nce2);
  s2c_2d(geom, geom->h1au1, dumpdata->h1au1, nfe1, nce2);
  s2c_2d(geom, geom->h1au2, dumpdata->h1au2, nce1, nfe2);
  s2c_2d(geom, geom->h1acell, dumpdata->h1acell, nce1, nce2);
  /* s2c_2d(geom,geom->h1agrid,dumpdata->h1agrid,nfe1,nfe2); */
  s2c_2d(geom, geom->thetau1, dumpdata->thetau1, nfe1, nce2);
  s2c_2d(geom, geom->h2au1, dumpdata->h2au1, nfe1, nce2);
  s2c_2d(geom, geom->h2au2, dumpdata->h2au2, nce1, nfe2);
  s2c_2d(geom, geom->h2acell, dumpdata->h2acell, nce1, nce2);
  /* s2c_2d(geom,geom->h2agrid,dumpdata->h2agrid,nfe1,nfe2); */
  s2c_2d(geom, geom->thetau2, dumpdata->thetau2, nfe1, nce2);
  s2c_2d(geom, master->coriolis, dumpdata->coriolis, nce1, nce2);
  emstag(LDEBUG,"hd:dumpdata:dumpdata_init","Dump data initialised OK");
}

/* END dumpdata_init()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to free memory in the dump data structure                 */
/*-------------------------------------------------------------------*/
void dumpdata_cleanup(dump_data_t *dumpdata,  /* Dump data structure */
                      int mode  /* Code for variables to free */
  )
{
  if (mode & 1) {
    d_free_1d(dumpdata->gridz);
    d_free_1d(dumpdata->cellz);
    d_free_2d(dumpdata->gridx);
    d_free_2d(dumpdata->gridy);
    d_free_2d(dumpdata->cellx);
    d_free_2d(dumpdata->celly);
    d_free_2d(dumpdata->u1x);
    d_free_2d(dumpdata->u1y);
    d_free_2d(dumpdata->u2x);
    d_free_2d(dumpdata->u2y);
    d_free_2d(dumpdata->botz);
    d_free_2d(dumpdata->h1au1);
    d_free_2d(dumpdata->h1au2);
    d_free_2d(dumpdata->h1acell);
    d_free_2d(dumpdata->h1agrid);
    d_free_2d(dumpdata->thetau1);
    d_free_2d(dumpdata->h2au1);
    d_free_2d(dumpdata->h2au2);
    d_free_2d(dumpdata->h2acell);
    d_free_2d(dumpdata->h2agrid);
    d_free_2d(dumpdata->thetau2);
    d_free_2d(dumpdata->coriolis);
    l_free_3d((long ***)dumpdata->flag);
    s_free_1d(dumpdata->clci);
    s_free_1d(dumpdata->crci);
    s_free_1d(dumpdata->flci);
    s_free_1d(dumpdata->frci);
    s_free_1d(dumpdata->clfi);
    s_free_1d(dumpdata->crfi);
    s_free_1d(dumpdata->flfi);
    s_free_1d(dumpdata->frfi);
    s_free_1d(dumpdata->cbcj);
    s_free_1d(dumpdata->cfcj);
    s_free_1d(dumpdata->fbcj);
    s_free_1d(dumpdata->ffcj);
    s_free_1d(dumpdata->cbfj);
    s_free_1d(dumpdata->cffj);
    s_free_1d(dumpdata->fbfj);
    s_free_1d(dumpdata->fffj);
    if (dumpdata->vmap) i_free_3d(dumpdata->vmap);
  }
  if (mode & 2) {
    d_free_2d(dumpdata->u1av);
    d_free_2d(dumpdata->u2av);
    d_free_2d(dumpdata->wind1);
    d_free_2d(dumpdata->wind2);
    d_free_2d(dumpdata->wtop);
    d_free_2d(dumpdata->topz);
    d_free_2d(dumpdata->eta);
    d_free_2d(dumpdata->patm);
    d_free_3d(dumpdata->u1);
    d_free_3d(dumpdata->u2);
    d_free_3d(dumpdata->w);
    d_free_3d(dumpdata->dzu1);
    d_free_3d(dumpdata->dzu2);
    d_free_3d(dumpdata->dzcell);
    d_free_2d(dumpdata->u1vh);
    d_free_2d(dumpdata->u2vh);
    d_free_2d(dumpdata->u1bot);
    d_free_2d(dumpdata->u2bot);
    d_free_4d(dumpdata->tr_wc);
    if (dumpdata->ntrS)
      d_free_3d(dumpdata->tr_wcS);
    if (dumpdata->sednz && dumpdata->nsed)
      d_free_4d(dumpdata->tr_sed);

    /* Time dependent diagnostic dump variables */
    d_free_3d(dumpdata->dens);
    d_free_3d(dumpdata->dens_0);
    d_free_3d(dumpdata->Kz);
    d_free_3d(dumpdata->Vz);
    d_free_2d(dumpdata->Cd);

    d_free_1d(dumpdata->w1);
    d_free_2d(dumpdata->w2);
    d_free_3d(dumpdata->w3);
    if (dumpdata->fmap_i != NULL)
      i_free_3d(dumpdata->fmap_i);
    if (dumpdata->fmap_j != NULL)
      i_free_3d(dumpdata->fmap_j);
    if (dumpdata->momgrid != NULL)
      mom_grid_free(dumpdata->momgrid, ALL);
    if (dumpdata->romsgrid != NULL)
      roms_grid_free(dumpdata->romsgrid);
  }
}

/* END dumpdata_cleanup()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to populate the dump structure with 2d data               */
/*-------------------------------------------------------------------*/
void dumpdata_fill(geometry_t *geom,  /* Sparse global geometry structure */
                   master_t *master,  /* Master data structure */
                   dump_data_t *dumpdata  /* Dump data structure */
  )
{
  int nce1 = geom->nce1;
  int nce2 = geom->nce2;
  int nfe1 = geom->nfe1;
  int nfe2 = geom->nfe2;
  int nz = geom->nz;
  int n, k = 0;
  int sednz = geom->sednz;

  /* Only the master dumps */
  if (master->mpi_rank) return;

  dumpdata->t = master->t;
  dumpdata->runcode = master->crf;

#if defined(HAVE_OMP) || defined(HAVE_MPI)
#pragma omp parallel sections 
#endif
  {
#if defined(HAVE_OMP) || defined(HAVE_MPI)
#pragma omp section
#endif
    s2c_2d(geom, master->u1av, dumpdata->u1av, nfe1, nce2);
#if defined(HAVE_OMP) || defined(HAVE_MPI)
#pragma omp section
#endif
    s2c_2d(geom, master->u2av, dumpdata->u2av, nce1, nfe2);
#if defined(HAVE_OMP) || defined(HAVE_MPI)
#pragma omp section
#endif
    s2c_2d(geom, master->wind1, dumpdata->wind1, nfe1, nce2);
#if defined(HAVE_OMP) || defined(HAVE_MPI)
#pragma omp section
#endif
    s2c_2d(geom, master->wind2, dumpdata->wind2, nce1, nfe2);
#if defined(HAVE_OMP) || defined(HAVE_MPI)
#pragma omp section
#endif
    s2c_2d(geom, master->wtop, dumpdata->wtop, nce1, nce2);
#if defined(HAVE_OMP) || defined(HAVE_MPI)
#pragma omp section
#endif
    s2c_2d(geom, master->topz, dumpdata->topz, nce1, nce2);
#if defined(HAVE_OMP) || defined(HAVE_MPI)
#pragma omp section
#endif
    s2c_2d(geom, master->eta, dumpdata->eta, nce1, nce2);
#if defined(HAVE_OMP) || defined(HAVE_MPI)
#pragma omp section
#endif
    s2c_2d(geom, master->patm, dumpdata->patm, nce1, nce2);
    /*
      s2c_3d_e1(geom, master->u1, dumpdata->u1, nfe1, nce2, nz);
      s2c_3d_e2(geom, master->u2, dumpdata->u2, nce1, nfe2, nz);
    */
#if defined(HAVE_OMP) || defined(HAVE_MPI)
#pragma omp section
#endif
    s2c_3d(geom, master->u1, dumpdata->u1, nfe1, nce2, nz);
#if defined(HAVE_OMP) || defined(HAVE_MPI)
#pragma omp section
#endif
    s2c_3d(geom, master->u2, dumpdata->u2, nce1, nfe2, nz);
#if defined(HAVE_OMP) || defined(HAVE_MPI)
#pragma omp section
#endif
    s2c_3d(geom, master->w, dumpdata->w, nce1, nce2, nz);
#if defined(HAVE_OMP) || defined(HAVE_MPI)
#pragma omp section
#endif
    s2c_3d(geom, master->dzu1, dumpdata->dzu1, nfe1, nce2, nz);
#if defined(HAVE_OMP) || defined(HAVE_MPI)
#pragma omp section
#endif
    s2c_3d(geom, master->dzu2, dumpdata->dzu2, nce1, nfe2, nz);
#if defined(HAVE_OMP) || defined(HAVE_MPI)
#pragma omp section
#endif
    s2c_3d(geom, master->dz, dumpdata->dzcell, nfe1, nfe2, nz);
#if defined(HAVE_OMP) || defined(HAVE_MPI)
#pragma omp section
#endif
    s2c_2d(geom, master->u1vh, dumpdata->u1vh, nce1, nce2);
#if defined(HAVE_OMP) || defined(HAVE_MPI)
#pragma omp section
#endif
    s2c_2d(geom, master->u2vh, dumpdata->u2vh, nce1, nce2);
#if defined(HAVE_OMP) || defined(HAVE_MPI)
#pragma omp section
#endif
    s2c_2d(geom, master->u1bot, dumpdata->u1bot, nfe1, nce2);
#if defined(HAVE_OMP) || defined(HAVE_MPI)
#pragma omp section
#endif
    s2c_2d(geom, master->u2bot, dumpdata->u2bot, nce1, nfe2);
  /* 
   * Even though these means are tracers, they still need to be listed
   * explicitly as they are face values
   */
#if defined(HAVE_OMP) || defined(HAVE_MPI)
#pragma omp section
#endif
    {
      if (master->u1m)
	s2c_3d(geom, master->u1m, dumpdata->u1m, nfe1, nce2, nz);
    }
#if defined(HAVE_OMP) || defined(HAVE_MPI)
#pragma omp section
#endif
    {
      if (master->u2m)
	s2c_3d(geom, master->u2m, dumpdata->u2m, nce1, nfe2, nz);
    }
 
#if defined(HAVE_OMP) || defined(HAVE_MPI)
#pragma omp section
#endif
    {
      if (master->u1vm)
	s2c_3d(geom, master->u1vm, dumpdata->u1vm, nfe1, nce2, nz);
    }
#if defined(HAVE_OMP) || defined(HAVE_MPI)
#pragma omp section
#endif
    {
      if (master->u2vm)
	s2c_3d(geom, master->u2vm, dumpdata->u2vm, nce1, nfe2, nz);
    }
  }
#if defined(HAVE_OMP) || defined(HAVE_MPI)
#pragma omp parallel for
#endif
  for (n = 0; n < dumpdata->ntr; n++) {
    s2c_3d(geom, master->tr_wc[n], dumpdata->tr_wc[n], nce1, nce2, nz);
  }
#if defined(HAVE_OMP) || defined(HAVE_MPI)
#pragma omp parallel for
#endif
  for (n = 0; n < master->ntrS; n++) {
    s2c_2d(geom, master->tr_wcS[n], dumpdata->tr_wcS[n], nce1, nce2);
  }
  if (sednz) {
    /*
    s2c_sed(geom, geom->cellz_sed, dumpdata->cellz_sed, nce1, nce2, sednz);
    s2c_sed(geom, geom->gridz_sed, dumpdata->gridz_sed, nce1, nce2, sednz+1);
    for (n = 0; n < master->nsed; n++) {
      s2c_sed(geom, master->tr_sed[n], dumpdata->tr_sed[n], nce1, nce2, sednz);
    } 
    */
    for (k = 0; k < sednz; k++) {
      s2c_2d(geom, geom->cellz_sed[k], dumpdata->cellz_sed[k], nce1, nce2);
      s2c_2d(geom, geom->gridz_sed[k], dumpdata->gridz_sed[k], nce1, nce2);
      for (n = 0; n < master->nsed; n++) {
        s2c_2d(geom, master->tr_sed[n][k], dumpdata->tr_sed[n][k], nce1,
               nce2);
      }
    }
    s2c_2d(geom, geom->gridz_sed[sednz], dumpdata->gridz_sed[sednz], nce1,
           nce2);
  }
#if defined(HAVE_OMP) || defined(HAVE_MPI)
#pragma omp parallel sections 
#endif
  {
#if defined(HAVE_OMP) || defined(HAVE_MPI)
#pragma omp section
#endif
    s2c_3d(geom, master->dens, dumpdata->dens, nce1, nce2, nz);
#if defined(HAVE_OMP) || defined(HAVE_MPI)
#pragma omp section
#endif
    s2c_3d(geom, master->dens_0, dumpdata->dens_0, nce1, nce2, nz);

#if defined(HAVE_OMP) || defined(HAVE_MPI)
#pragma omp section
#endif
    s2c_3d(geom, master->Kz, dumpdata->Kz, nce1, nce2, nz);

#if defined(HAVE_OMP) || defined(HAVE_MPI)
#pragma omp section
#endif
    s2c_3d(geom, master->Vz, dumpdata->Vz, nce1, nce2, nz);

#if defined(HAVE_OMP) || defined(HAVE_MPI)
#pragma omp section
#endif
    s2c_2d(geom, master->Cd, dumpdata->Cd, nce1, nce2);
  }
}

/* END dumpdata_fill()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Creates a map from sigma coordinates to Cartesian                 */
/*-------------------------------------------------------------------*/
void sigma_vmap(master_t *master,       /* Global geometry structure */
		dump_data_t *dumpdata   /* Dump data structure       */
		)
{
  geometry_t *geom = master->geom;
  int i, j, k, kk, c, cs;
  int nz = geom->nz - 1;
  double d;

  dumpdata->vmap = i_alloc_3d(geom->nfe1, geom->nfe2, geom->nz);

  for (j = 0; j < geom->nce2; j++) {	
    for (i = 0; i < geom->nce1; i++) {
      cs = geom->map[nz][j][i];
      for (k = nz; k >= 0; k--) {
	if (dumpdata->flag[nz][j][i] & (SOLID|OUTSIDE)) {
	  dumpdata->vmap[k][j][i] = 0;
	  continue;
	}
	c = cs;
	kk = nz;
	while (c != geom->zm1[c] && (fabs(geom->layers[k]) > 
	       fabs(dumpdata->cellz[kk] * dumpdata->botz[j][i]))) {
	  c = geom->zm1[c];
	  kk--;
	}
	dumpdata->vmap[k][j][i] = geom->zp1[c];
	if (geom->layers[k] < dumpdata->botz[j][i])
	  dumpdata->vmap[k][j][i] = 0;
      }
    }
  }
}

/* END 
/*-------------------------------------------------------------------*/

static int contains_string(char *p, char *tag)
{
  char *args[256];
  int len = strlen(p);
  int i, n;

  for (i = 0; i < len; ++i) {
    if (p[i] == '_' || p[i] == '.')
      p[i] = ' ';
  }

  i = parseline(p, args, 256);
  for (n = 0; n < i; n++)
    if (strcmp(args[n], tag) == 0)
      return 1;
  return 0;
}



