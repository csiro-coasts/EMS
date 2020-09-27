/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/boundary/bdryeval.c
 *  
 *  Description:
 *  Routines to read (and write) lists of specified
 *  boundary conditions from a text file
 *  (usually the parameter file)
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: bdryeval.c 6384 2019-11-21 22:54:21Z her127 $
 *
 */

#include <stdio.h>
#include "hd.h"

void read_bdry_zone(master_t *master, open_bdrys_t *open, int cc, int mode);
double adjust_OBC_mass(geometry_t *window, window_t *windat, win_priv_t *wincon, open_bdrys_t *open, int tn);
double adjust_OBC_maccready(geometry_t *window, window_t *windat, win_priv_t *wincon, open_bdrys_t *open, int tn);

/*-------------------------------------------------------------------*/
/* Call the boundary initialisation routines for all custom routines */
/*-------------------------------------------------------------------*/
void bdry_init_m(master_t *master)
{
  geometry_t *geom = master->sgrid;  /* Global geometry               */
  open_bdrys_t **open = geom->open; /* Open boudary data             */
  int n, t;                         /* Counters                      */

  for (n = 0; n < geom->nobc; ++n) {

    /*---------------------------------------------------------------*/
    /* Initialize the boundary data for eta, u1 or u2.               */
    /*---------------------------------------------------------------*/
    /* eta boundary data                                             */
    if (open[n]->etadata.init_m != NULL) {
      open[n]->etadata.init_m(master, open[n], &open[n]->etadata);
    } else if (!open[n]->etadata.explct && open[n]->bcond_ele & FILEIN)
      hd_ts_multifile_check(open[n]->ntsfiles, open[n]->tsfiles,
                            open[n]->filenames, open[n]->etadata.name,
                            schedule->start_time, schedule->stop_time);

    /*---------------------------------------------------------------*/
    /* Boundary velocity                                             */
    /* u1 (3D) boundary data                                         */
    if (open[n]->datau1.explct && (open[n]->bcond_nor & FILEIN ||
				   open[n]->bcond_tan & FILEIN))
      hd_warn("WARNING : Use the CUSTOM specification for boundary %d\n", n);		

    /* Normal CUSTOM forcing on a boundary                           */
    if (open[n]->datau1.init_m != NULL)
      open[n]->datau1.init_m(master, open[n], &open[n]->datau1);
    else if (!open[n]->datau1.explct && open[n]->bcond_nor & FILEIN) {
      /* Normal component FILEIN forcing on a boundary               */
      hd_ts_multifile_check(open[n]->ntsfiles, open[n]->tsfiles,
			    open[n]->filenames, open[n]->datau1.name,
			    schedule->start_time, schedule->stop_time);
    }

    /* Tangential CUSTOM forcing on a boundary                       */
    if (open[n]->datau2.init_m != NULL)
      open[n]->datau2.init_m(master, open[n], &open[n]->datau2);
    else if (!open[n]->datau2.explct && open[n]->bcond_tan & FILEIN) {
      /* Tangential component FILEIN forcing on a boundary           */
      hd_ts_multifile_check(open[n]->ntsfiles, open[n]->tsfiles,
			    open[n]->filenames, open[n]->datau2.name,
			    schedule->start_time, schedule->stop_time);
    }

    /* u1av (2D) boundary data                                       */
    if (open[n]->datau1av.explct && (open[n]->bcond_nor2d & FILEIN ||
				     open[n]->bcond_tan2d & FILEIN))
      hd_warn("WARNING : Use the CUSTOM specification for boundary %d\n", n);
    
    /* 2D CUSTOM forcing on a boundary                               */
    if (open[n]->datau1av.init_m != NULL)
      open[n]->datau1av.init_m(master, open[n], &open[n]->datau1av);
    else if (!open[n]->datau1av.explct && open[n]->bcond_nor2d & FILEIN) {
      /* Normal 2D component FILEIN forcing on a boundary            */
      hd_ts_multifile_check(open[n]->ntsfiles, open[n]->tsfiles,
			    open[n]->filenames, open[n]->datau1av.name,
			    schedule->start_time, schedule->stop_time);
    }
    
    /* Tangential 2D CUSTOM forcing on a boundary                    */
    if (open[n]->datau2av.init_m != NULL)
      open[n]->datau2av.init_m(master, open[n], &open[n]->datau2av);
    else if (!open[n]->datau2av.explct && open[n]->bcond_tan2d & FILEIN) {
      /* Tangential 2D component FILEIN forcing on a boundary        */
      hd_ts_multifile_check(open[n]->ntsfiles, open[n]->tsfiles,
			    open[n]->filenames, open[n]->datau2av.name,
			    schedule->start_time, schedule->stop_time);
    }

    /*---------------------------------------------------------------*/
    /* Initialize boundary data for tracers                          */
    for (t = 0; t < open[n]->ntr; ++t) {
      if (open[n]->bdata_t[t].init_m != NULL) {
        open[n]->bdata_t[t].init_m(master, open[n], &open[n]->bdata_t[t]);
      } else if (!open[n]->bdata_t[t].explct &&
                 open[n]->bcond_tra[t] & FILEIN) {
        hd_ts_multifile_check(open[n]->ntsfiles, open[n]->tsfiles,
                              open[n]->filenames, open[n]->bdata_t[t].name,
                              schedule->start_time, schedule->stop_time);
      }
    }
  }
}

/* END bdry_init_m()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Call the boundary initialisation routines for all custom routines */
/* for the slave.                                                    */
/*-------------------------------------------------------------------*/
void bdry_init_w(geometry_t *window,  /* Window geometry             */
                 open_bdrys_t *open,  /* Window open boudary data    */
                 open_bdrys_t *gopen  /* Global OBC structure in this
                                         window */
  )
{
  int t;                              /* Counters                    */

  /*-----------------------------------------------------------------*/
  /* Initialize the boundary data for eta or u1. Note: reading       */
  /* forcing data from file is always done on the master so there    */
  /* is no need to chech timeseries files here.                      */
  if (open->etadata.init_w != NULL)
    open->etadata.init_w(window, open, &open->etadata, &gopen->etadata);
  if (open->datau1.init_w != NULL)
    open->datau1.init_w(window, open, &open->datau1, &gopen->datau1);
  if (open->datau2.init_w != NULL)
    open->datau2.init_w(window, open, &open->datau2, &gopen->datau2);
  if (open->datau1av.init_w != NULL)
    open->datau1av.init_w(window, open, &open->datau1av, &gopen->datau1av);
  if (open->datau2av.init_w != NULL)
    open->datau2av.init_w(window, open, &open->datau2av, &gopen->datau2av);

  /*-----------------------------------------------------------------*/
  /* Initialize boundary data for tracers                            */
  for (t = 0; t < open->ntr; ++t) {
    if (open->bdata_t[t].init_w != NULL) {
      open->bdata_t[t].init_w(window, open, &open->bdata_t[t],
                              &gopen->bdata_t[t]);
    }
  }
}

/* END bdry_init_w()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Evaluate this bdry_details_t for the specified location in sparse */
/* coordinates.                                                      */
/*-------------------------------------------------------------------*/
double bdry_value_w(geometry_t *window, /* Window geometry           */
                    window_t *windat,   /* Window data               */
                    win_priv_t *wincon, /* Window constants          */
                    open_bdrys_t *open, /* Open boundary data        */
                    bdry_details_t *data, /* Custom data             */
                    int c,                /* 3D Sparse location      */
                    int cc,               /* 2D Sparse location      */
                    double t              /* Time (s)                */
  )
{

  /*-----------------------------------------------------------------*/
  /* Evaluate at the position */
  if (data->explct) {
    if (data->custom_w != NULL) {
      return data->custom_w(window, windat, wincon, open, t, c, cc, data);
    }
  }
  return data->fill_value;
}

/* END bdry_value_w()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to invoke custom routines for u1 velocity on the master.  */
/* This routine is called when the master is filled with nu1         */
/* velocities from the slaves.                                       */
/*-------------------------------------------------------------------*/
void bdry_eval_u1_m(geometry_t *geom, /* Global geometry             */
                    master_t *master  /* Global data                 */
  )
{
  int n;                        /* Counters                          */
  int e, ee, e2;                /* Edge coordinates / counters       */
  int c, c2;                    /* Cell centre coordinate            */
  double x, y, z;               /* Geographic coordinates            */
  open_bdrys_t **open = geom->open;

  /*-----------------------------------------------------------------*/
  /* Loop through the open boundaries                                */
  for (n = 0; n < geom->nobc; n++) {
    /* Set the normal velocity to the boundary                       */
    if (open[n]->bcond_nor & (FILEIN | CUSTOM)) {
      bdry_details_t *data = &open[n]->datau1;
      /* Call the custom boundary function                           */
      if (data->explct) {
	if (data->custom_m != NULL) {
	  data->custom_m(geom, master, open[n], data);
	} else {
	  for (ee = 1; ee <= open[n]->no3_e1; ee++) {
	    e = open[n]->obc_e1[ee];
	    open[n]->transfer_u1[ee] = data->fill_value;
	  }
	}
      }
      /* Read the forcing data from file                             */
      else {
	double ramp = (master->rampf & FILEIN) ? master->rampval : 1.0;
	for (ee = 1; ee <= open[n]->no3_e1; ee++) {
	  e = open[n]->obc_e1[ee];
	  c = open[n]->obc_e2[ee];
	  e2 = geom->m2de[e];
	  c2 = geom->m2d[c];
	  x = geom->u1x[e2];
	  y = geom->u1y[e2];
	  z = geom->cellz[c] * master->Ds[c2];
	  open[n]->transfer_u1[ee] = ramp * master->Hn1[c2] * 
	    hd_ts_multifile_eval_xyz_by_name(open[n]->ntsfiles, 
					     open[n]->tsfiles,
					     open[n]->filenames, data->name,
					     master->t, x, y, z);
	  if (open[n]->relax_zone_nor)
	    read_bdry_zone(master, open[n], ee, U1BDRY|U1GEN);
	}
      }
    }

    /* Set the tangential velocity to the boundary                   */
    if (open[n]->bcond_tan & (FILEIN | CUSTOM)) {
      bdry_details_t *data = &open[n]->datau2;
      /* Call the custom boundary function                           */
      if (data->explct) {
	if (data->custom_m != NULL) {
	  data->custom_m(geom, master, open[n], data);
	}
      }
      /* Read the forcing data from file                             */
      else {
	double ramp = (master->rampf & FILEIN) ? master->rampval : 1.0;
	for (ee = open[n]->no3_e1 + 1; ee <= open[n]->to3_e1; ee++) {
	  e = open[n]->obc_e1[ee];
	  c = open[n]->obc_e2[ee];
	  e2 = geom->m2de[e];
	  c2 = geom->m2d[c];
	  x = geom->u1x[e2];
	  y = geom->u1y[e2];
	  z = geom->cellz[c] * master->Ds[c2];
	  open[n]->transfer_u2[ee] = ramp * master->Hn2[c2] *
	    hd_ts_multifile_eval_xyz_by_name(open[n]->ntsfiles, 
					     open[n]->tsfiles,
					     open[n]->filenames, 
					     data->name,
					     master->t, x, y, z);
	  if (open[n]->relax_zone_tan)
	    read_bdry_zone(master, open[n], ee, U1BDRY|U2GEN);
	}
      }
    }
  }
  if (master->regf == RS_OBCSET) 
    sprintf(master->runerror, "OBC %d u1\n",open[n]->id);
}

/* END bdry_eval_u1_m()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to transfer the custom u1 boundary data to the slaves     */
/*-------------------------------------------------------------------*/
void bdry_transfer_u1(master_t *master, geometry_t *window,
                      window_t *windat)
{
  geometry_t *geom = master->geom;    /* Global geometry             */
  open_bdrys_t **open_m = geom->open; /* Global open boundary data   */
  open_bdrys_t **open_w = window->open; /* Window open boundary data */
  int n;                                /* Counters                  */
  int e, ee;                            /* Sparse location / counter */

  /* Always do the custom transfers for 1 window; the window trasfer */
  /* vectors point to the master vectors for 1 window.               */
  for (n = 0; n < window->nobc; n++) {

    /* Transfer normal velocity to the boundary                      */
    if (open_w[n]->bcond_nor & (FILEIN | CUSTOM)) {
      bdry_details_t *data = &open_w[n]->datau1;
      /* Do the custom transfer if required                          */
      if (data->explct) {
	if(data->trans)
	  data->trans(master, open_w[n], data, window, windat);
	else {
	  for (ee = 1; ee <= open_w[n]->no3_e1; ee++)
	    open_w[n]->transfer_u1[ee] = data->fill_value;
	}
      }
    }

    /* Transfer tangential velocity to the boundary                  */
    if (open_w[n]->bcond_tan & (FILEIN | CUSTOM)) {
      bdry_details_t *data = &open_w[n]->datau2;
      /* Do the custom transfer if required                          */
      if (data->explct) {
	if(data->trans)
	  data->trans(master, open_w[n], data, window, windat);
	else {
	  for (ee = open_w[n]->no3_e1 + 1; ee <= open_w[n]->to3_e1; ee++)
	    open_w[n]->transfer_u2[ee] = data->fill_value;
	}
      }
    }
  }

  /* Transfer data read from file from the master to the slave if    */
  /* multiple windows exist.                                         */
  if (master->nwindows > 1) {
    for (n = 0; n < window->nobc; n++) {

      /* Transfer normal velocity to the boundary                    */
      if (open_w[n]->bcond_nor & (FILEIN | CUSTOM)) {
	int nz = (open_w[n]->relax_zone_nor) ? open_w[n]->relax_zone_nor : 1;
	int mwn = open_w[n]->mwn;
	for (ee = 1; ee <= nz * open_w[n]->no3_e1; ee++) {
	  e = open_w[n]->tmap_u1[ee];
	  open_w[n]->transfer_u1[ee] = open_m[mwn]->transfer_u1[e];
	}
      }
      /* Transfer tangential velocity to the boundary                */
      if (open_w[n]->bcond_tan & (FILEIN | CUSTOM)) {
	int mwn = open_w[n]->mwn;

	for (ee = open_w[n]->no3_e1 + 1; ee <= open_w[n]->to3_e1; ee++) {
	  e = open_w[n]->tmap_u2[ee];
	  open_w[n]->transfer_u2[ee] = open_m[mwn]->transfer_u2[e];
	}
      }
    }
  }
}

/* END bdry_transfer_u1()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to invoke custom routines for elevation on the master.    */
/* This routine is called when the master is filled with elevation   */
/* values from the slaves.                                           */
/*-------------------------------------------------------------------*/
void bdry_eval_eta_m(geometry_t *geom,  /* Global geometry */
                     master_t *master /* Global data */
  )
{
  int n;                        /* Counters */
  int c, cc;                    /* Sparse coordinates / counters */
  double x, y, z;               /* Geographic coordinates */
  open_bdrys_t **open = geom->open;

  /*-----------------------------------------------------------------*/
  /* Loop through the open boundaries */
  for (n = 0; n < geom->nobc; n++) {

    /*---------------------------------------------------------------*/
    /* Set the elevation */
    if (open[n]->bcond_ele & (FILEIN | CUSTOM)) {
      bdry_details_t *data = &open[n]->etadata;
      /* Call the custom boundary function */
      if (data->explct) {
        if (data->custom_m != NULL) {
          data->custom_m(geom, master, open[n], data);
        } else {
	  for (cc = 1; cc <= open[n]->no2_t; cc++) {
	    c = open[n]->obc_t[cc];
	    open[n]->transfer_eta[cc] = data->fill_value;
	  }
	}
      }
      /* Read the forcing data from file. Note the time is that of   */
      /* the 2D mode at the forward time, t+1.                       */
      else {
	double ramp = (master->rampf & FILEIN) ? master->rampval : 1.0;
	if (open[n]->file_dt && master->t < open[n]->file_next - DT_EPS)
	  continue;

        for (cc = 1; cc <= open[n]->no2_t; cc++) {
          c = open[n]->obc_t[cc];
          x = geom->cellx[c];
          y = geom->celly[c];
          z = geom->cellz[c] * master->Ds[c];
          open[n]->transfer_eta[cc] = ramp *
            hd_ts_multifile_eval_xyz_by_name(open[n]->ntsfiles, 
					     open[n]->tsfiles,
					     open[n]->filenames, 
					     data->name,
					     master->t3d,
					     /*master->t + master->dt2d, */
					     x, y, z);
        }
	open[n]->file_next += open[n]->file_dt;
      }
    }
    if (master->regf == RS_OBCSET) 
      sprintf(master->runerror, "OBC %d eta\n",open[n]->id);
  }
}

/* END bdry_eval_eta_m()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to transfer the custom eta boundary data to the slaves    */
/*-------------------------------------------------------------------*/
void bdry_transfer_eta(master_t *master, geometry_t *window,
                       window_t *windat)
{
  geometry_t *geom = master->geom;
  open_bdrys_t **open_m = geom->open;
  open_bdrys_t **open_w = window->open;
  int n;                        /* Counters */
  int c, cc;                    /* Sparse location / counter */

  /* Always do the custom transfers for 1 window; the window trasfer */
  /* vectors point to the master vectors for 1 window.  */
  for (n = 0; n < window->nobc; n++) {
    if (open_w[n]->bcond_ele & (FILEIN | CUSTOM)) {
      bdry_details_t *data = &open_w[n]->etadata;

      /* Do the custom transfer if required */
      if (data->explct) {
	if(data->trans)
	  data->trans(master, open_w[n], data, window, windat);
	else {
	  for (cc = 1; cc <= open_w[n]->no2_t; cc++)
	    open_w[n]->transfer_eta[cc] = data->fill_value;
	}
      }
    }
  }

  /* Transfer data read from file from the master to the slave if */
  /* multiple windows exist.  */
  if (master->nwindows > 1) {
    for (n = 0; n < window->nobc; n++) {
      if (open_w[n]->bcond_ele & (FILEIN | CUSTOM)) {
        int mwn = open_w[n]->mwn;
        for (cc = 1; cc <= open_w[n]->no2_t; cc++) {
          c = open_w[n]->tmap[cc];
          open_w[n]->transfer_eta[cc] = open_m[mwn]->transfer_eta[c];
        }
      }
    }
  }
}

/* END bdry_transfer_eta()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to invoke custom routines for tracers on the master.      */
/* This routine is called when the master is filled with tracer      */
/* values from the slaves.                                           */
/*-------------------------------------------------------------------*/
void bdry_eval_tr_m(geometry_t *geom, /* Global geometry             */
                    master_t *master  /* Global data                 */
  )
{
  int n, m, tn, tm;             /* Counters                          */
  int c, cc, c2;                /* Centre coordinates / counters     */
  int e, ee;                    /* Edge coordinate / counter         */
  double x, y, z;               /* Geographic coordinates            */
  open_bdrys_t **open = geom->open;

  /*-----------------------------------------------------------------*/
  /* Loop through the open boundaries                                */
  for (n = 0; n < geom->nobc; n++) {

    /* Set the maps if TRCONC is used on any tracer                  */
    if (open[n]->bgz) {
      for (ee = 1; ee <= open[n]->no3_e1; ee++) {
	c = open[n]->obc_e2[ee];
	open[n]->omape[ee][c] = open[n]->ogc_t[ee];
      }
    }

    for (tn = 0; tn < open[n]->ntr; tn++) {

      /*-------------------------------------------------------------*/
      /* Set the tracers                                             */
      if (open[n]->bcond_tra[tn] & (FILEIN | CUSTOM)) {
        bdry_details_t *data = &open[n]->bdata_t[tn];

        /* Call the custom boundary function                         */
        if (data->explct) {
          if (data->custom_m != NULL) {
            data->custom_m(geom, master, open[n], data);
          } else {
	    tm = open[n]->trm[tn];
	    for (cc = 1; cc <= open[n]->ttsz; cc++) {
	      open[n]->t_transfer[tm][cc] = data->fill_value;
	    }
	  }
        } else if (open[n]->bcond_tra[tn] & (TRCONC|TRCONF)) {
	  /* TRCONC : interpolate data into OBC ghost cells          */
	  /* Write into OBC cell centres if required.                */
	  if (open[n]->options & OP_OWRITE) {
	    for (cc = 1; cc <= open[n]->no3_e1; cc++) {
	      c = open[n]->obc_t[ee];
	      c2 = geom->m2d[c];
	      z = geom->cellz[c] * master->Ds[c2];
	      tm = open[n]->trm[tn];
	      x = geom->cellx[c2];
	      y = geom->celly[c2];
	      open[n]->t_transfer[tm][cc] =
		hd_ts_multifile_eval_xyz_by_name(open[n]->ntsfiles, 
						 open[n]->tsfiles,
						 open[n]->filenames, 
						 data->name,
						 master->t, x, y, z);
	    }
	  }

	  /* Write into OBC ghost centres                            */
	  for (ee = 1; ee <= open[n]->no3_e1; ee++) {
	    c = open[n]->obc_e2[ee];
	    e = open[n]->obc_e1[ee];
	    c2 = geom->m2d[c];
	    z = geom->cellz[c] * master->Ds[c2];
	    tm = open[n]->trm[tn];
	    for(m = 1; m <= open[n]->bgz; m++) {
	      double tr;
	      cc = open[n]->t_imap[ee][m];
	      c2 = open[n]->omape[ee][c2];
	      x = geom->cellx[c2];
	      y = geom->celly[c2];
	      tr = hd_ts_multifile_eval_xyz_by_name(open[n]->ntsfiles, 
						    open[n]->tsfiles,
						    open[n]->filenames, 
						    data->name,
						    master->t, x, y, z);
	      /* NaN check : if interpolating from a depth below the */
	      /* sea bed, and the interpolated value is NaN, then    */
	      /* set to the layer above.                             */
	      if (isnan(tr))
		hd_quit("bdry_eval_tr_m: NaN found in boundary data for tracer %d.\n", tn);

	      open[n]->t_transfer[tm][cc] = tr;
	    }
	  }
	} else {
	  /* Read the forcing data from file                         */
	  for (cc = 1; cc <= open[n]->no3_t; cc++) {
	    c = open[n]->obc_t[cc];
	    c2 = geom->m2d[c];
	    x = geom->cellx[c2];
	    y = geom->celly[c2];
	    z = geom->cellz[c] * master->Ds[c2];
	    tm = open[n]->trm[tn];
	    open[n]->t_transfer[tm][cc] =
	      hd_ts_multifile_eval_xyz_by_name(open[n]->ntsfiles, 
					       open[n]->tsfiles,
					       open[n]->filenames, 
					       data->name,
					       master->t, x, y, z);
	  }
        }
      }
    }
    /* Reset the maps to self mapping across OBCs                    */
    if (open[n]->bgz) {
      for (ee = 1; ee <= open[n]->no3_e1; ee++) {
	c = open[n]->obc_e2[ee];
	open[n]->omape[ee][c] = c;
      }
    }
  }
}

/* END bdry_eval_tr_m()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to transfer the custom tracer boundary data to the slaves */
/*-------------------------------------------------------------------*/
void bdry_transfer_tr(master_t *master, geometry_t *window,
                      window_t *windat)
{
  geometry_t *geom = master->geom;
  win_priv_t *wincon = window->wincon;
  open_bdrys_t **open_m = geom->open;
  open_bdrys_t **open_w = window->open;
  int n, nn, m, tn, tm;         /* Counters */
  int c, ee, cc, cb;            /* Sparse location / counter */

  /* Always do the custom transfers for 1 window; the window trasfer */
  /* vectors point to the master vectors for 1 window.  */
  for (n = 0; n < window->nobc; n++) {
    for (tn = 0; tn < open_w[n]->ntr; tn++) {
      if (open_w[n]->bcond_tra[tn] & (FILEIN | CUSTOM)) {
        bdry_details_t *data = &open_w[n]->bdata_t[tn];
        /* Do the custom transfer if required */
        if (data->explct) {
	  if(data->trans) {
	    data->trans(master, open_w[n], data, window, windat);
	  }
	  else {
	    for (cc = 1; cc <= open_w[n]->ttsz; cc++) {
	      tm = open_w[n]->trm[tn];
	      open_w[n]->t_transfer[tm][cc] = data->fill_value;
	    }
	  }
        }
      }
    }
  }

  /* Transfer data read from file from the master to the slave if */
  /* multiple windows exist.  */
  if (master->nwindows > 1) {
    for (n = 0; n < window->nobc; n++) {
      for (tn = 0; tn < open_w[n]->ntr; tn++) {
	int mwn = open_w[n]->mwn;
        if (open_w[n]->bcond_tra[tn] & (FILEIN | CUSTOM)) {
	  tm = open_w[n]->trm[tn];
	  if (open_w[n]->bcond_tra[tn] & (TRCONC|TRCONF)) {
	    /* Copy OBC ghost cells only */
	    for (cc = 1; cc <= open_w[n]->ttsz; cc++) {
	      c = open_w[n]->t_tmap[cc];
	      open_w[n]->t_transfer[tm][cc] = 
		open_m[mwn]->t_transfer[tm][c];
	    }
	  } else {
	    for (cc = 1; cc <= open_w[n]->no3_t; cc++) {
	      c = open_w[n]->t_tmap[cc];
	      open_w[n]->t_transfer[tm][cc] = open_m[mwn]->t_transfer[tm][c];
	    }
	  }
        }
      }
    }
  }

  /* Modify the salinity in the boundary cell resulting from a salt  */
  /* mass balance if required, Balance is:                           */
  /* salt = (mass_river + mass_landward)/(vol_river + vol_landward)  */
  for (n = 0; n < window->nobc; n++) {
    for (tn = 0; tn < open_w[n]->ntr; tn++) {
      if(tn == windat->sno && 
	 (open_w[n]->bcond_tra[tn] == (TRCONC|CUSTOM) ||
	  open_w[n]->bcond_tra[tn] == (TRCONF|CUSTOM))) {
	if (!(open_w[n]->options & OP_NOSALT)) {
	  double val;
	  tm = open_w[n]->trm[tn];
	  if (open_w[n]->options & OP_MACREADY)
	    val = adjust_OBC_maccready(window, windat, wincon, open_w[n], tn);
	  else
	    val = adjust_OBC_mass(window, windat, wincon, open_w[n], tn);
	  for (cc = 1; cc <= open_w[n]->ttsz; cc++) {
	    open_w[n]->t_transfer[tm][cc] = val;	      
	  }
	}
	/* Save the ghost cell salinity                                 */
	if (open_w[n]->options & OP_UPSTRM) {
	  tm = open_w[n]->trm[tn];
	  for (ee = 1; ee <= open_w[n]->no3_e1; ee++) {
	    c = open_w[n]->ogc_t[ee];
	    wincon->w10[c] = windat->tr_wc[tn][c];
	  }
	}
      }
    }
  }

  /* Copy the value directly into the tracer OBC ghost cells for TRCONC */
  for (n = 0; n < window->nobc; n++) {
    for (nn = 0; nn < wincon->ntbdy; nn++) {
      tn = wincon->tbdy[nn];
      if (open_w[n]->bgz) {
	if (open_w[n]->bcond_tra[tn] == (TRCONC|FILEIN) ||
	    open_w[n]->bcond_tra[tn] == (TRCONC|CUSTOM) ||
	    open_w[n]->bcond_tra[tn] == (TRCONF|FILEIN) ||
	    open_w[n]->bcond_tra[tn] == (TRCONF|CUSTOM)) {
	  tm = open_w[n]->trm[tn];
	  /* Copy the value into the tracer ghost cells                 */
	  if (open_w[n]->options & OP_OWRITE) {
	    for (cc = 1; cc <= open_w[n]->no3_t; cc++) {
	      cb = open_w[n]->obc_t[cc];
	      windat->tr_wc[tn][cb] = open_w[n]->t_transfer[tm][cc];
	    }
	  }
	  for (ee = 1; ee <= open_w[n]->no3_e1; ee++) {
	    c = open_w[n]->ogc_t[ee];
	    for (m = 1; m <= open_w[n]->bgz; m++) {
	      cc = open_w[n]->t_imap[ee][m];
	      windat->tr_wc[tn][c] = open_w[n]->t_transfer[tm][cc];
	      c = open_w[n]->omape[ee][c];
	    }
	  }


	  /* Riverflow salinity options.                                */
	  /* Modify ghost values with an UPSTRM condition if required.  */
	  if (tn == windat->sno) {
	    if (open_w[n]->options & OP_UPSTRM) {
	      double sf = 1.0;
	      /* Uncomment for use with FACE in upstrm() 
	      memset(wincon->w10, 0, window->sgsiz * sizeof(double));
	      for (cc = 1; cc <= open_w[n]->no3_t; cc++) {
		c = open_w[n]->obc_t[cc];
		wincon->w10[c] = open_w[n]->t_transfer[tm][cc];
	      }
	      */
	      upstrm(window, windat, wincon, open_w[n], wincon->w10, windat->tr_wc[tn], GHOST);
	      for (ee = 1; ee <= open_w[n]->no3_e1; ee++) {
		c = open_w[n]->ogc_t[ee];
		cb = open_w[n]->obc_e2[ee];
		windat->tr_wc[tn][c] = sf * wincon->w10[c] + (1.0 - sf) * windat->tr_wc[tn][cb];

		/* Mix laterally 
		   windat->tr_wc[tn][c] = wincon->u1kh[cb] * (windat->tr_wc[tn][cb] + windat->tr_wc[tn][open_w[n]->omap[c]] - 2.0 * windat->tr_wc[tn][c]) / 2.0 * window->h2acell[window->m2d[cb]]; */

		/* Set a no gradient in other boundary ghost cells */
		cb = c;
		c = open_w[n]->omape[ee][c];
		for (m = 1; m < open_w[n]->bgz; m++) {
		  windat->tr_wc[tn][c] = windat->tr_wc[tn][cb];
		  c = open_w[n]->omape[ee][c];
		}
	      }
	    }
	    if (windat->riversalt) {
	      for (cc = 1; cc <= open_w[n]->no3_t; cc++) {
		c = open_w[n]->ogc_t[cc];
		cb = open_w[n]->obc_t[cc];
		windat->riversalt[cb] = windat->tr_wc[tn][c];
	      }
	    }
	  }
	  /* Save the ghost cell value for TRCONF implementation     */
	  if (open_w[n]->bcond_tra[tn] & TRCONF) {
	    int mm;
	    /* Get the index in the tracer dummy array corresponding */
	    /* to tn.                                                */
	    for (mm = 0; mm < open_w[n]->ntflx; mm++)
	      if (open_w[n]->tflx[mm] == tn) break;
	    /* Save the contribution to total flux through this OBC  */
	    for (ee = 1; ee <= open_w[n]->no3_e1; ee++) {
	      c = open_w[n]->ogc_t[ee];
	      open_w[n]->dumtr[mm][ee] = windat->tr_wc[tn][c];
	    }
	    /* Set a no-gradient condition in OBC ghost cells          */
	    OBC_bgz_nogradb(window, open_w[n], windat->tr_wc[tn]);
	  }
	} else if (open_w[n]->bcond_tra[tn] == (TRCONC|NOTHIN) ||
		   open_w[n]->bcond_tra[tn] == (TRCONC|NOGRAD)) {
	  /* Set a no-gradient condition in OBC ghost cells          */
	  OBC_bgz_nogradb(window, open_w[n], windat->tr_wc[tn]);
	} else if (!(open_w[n]->bcond_tra[tn] & (TRCONC|TRCONF))) {
	  /* Set a no-gradient over ghost OBCs */
	  OBC_bgz_nogradb(window, open_w[n], windat->tr_wc[tn]);
	}
      }
    }
  }
}

/* END bdry_transfer_tr()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to invoke custom routines for u1av velocity on the        */
/* master. This routine is called when the master is filled with     */
/* nu1av velocities from the slaves.                                 */
/*-------------------------------------------------------------------*/
void bdry_eval_u1av_m(geometry_t *geom, /* Global geometry           */
                    master_t *master    /* Global data               */
  )
{
  int n;                        /* Counters                          */
  int e, ee, e2;                /* Edge coordinates / counters       */
  int c, c2;                    /* Cell centre coordinate            */
  double x, y, z;               /* Geographic coordinates            */
  open_bdrys_t **open = geom->open;

  /*-----------------------------------------------------------------*/
  /* Loop through the open boundaries                                */
  for (n = 0; n < geom->nobc; n++) {
    /* Set the normal depth averaged velocity to the boundary        */
    if (open[n]->bcond_nor2d & (FILEIN | CUSTOM)) {
      bdry_details_t *data = &open[n]->datau1av;
      /* Call the custom boundary function                           */
      if (data->explct) {
	if (data->custom_m != NULL) {
	  data->custom_m(geom, master, open[n], data);
	} else {
	  for (ee = 1; ee <= open[n]->no2_e1; ee++) {
	    e = open[n]->obc_e1[ee];
	    open[n]->transfer_u1av[ee] = data->fill_value;
	  }
	}
      }
      /* Read the forcing data from file                             */
      else {
	double ramp = (master->rampf & FILEIN) ? master->rampval : 1.0;
	for (ee = 1; ee <= open[n]->no2_e1; ee++) {
	  e = open[n]->obc_e1[ee];
	  c = open[n]->obc_e2[ee];
	  e2 = geom->m2de[e];
	  c2 = geom->m2d[c];
	  x = geom->u1x[e2];
	  y = geom->u1y[e2];
	  z = geom->cellz[c] * master->Ds[c2];
	  open[n]->transfer_u1av[ee] = ramp * master->Hn1[c2] *
	    hd_ts_multifile_eval_xyz_by_name(open[n]->ntsfiles, 
					     open[n]->tsfiles,
					     open[n]->filenames, 
					     data->name,
					     master->t3d, x, y, z);
	}
      }
    }

    /* Set the tangential depth averaged velocity to the boundary    */
    if (open[n]->bcond_tan2d & (FILEIN | CUSTOM)) {
      bdry_details_t *data = &open[n]->datau2av;
      /* Call the custom boundary function                           */
      if (data->explct) {
	if (data->custom_m != NULL) {
	  data->custom_m(geom, master, open[n], data);
	}
      }
      /* Read the forcing data from file                             */
      else {
	double ramp = (master->rampf & FILEIN) ? master->rampval : 1.0;
	for (ee = open[n]->no3_e1 + 1; ee <= open[n]->to2_e1; ee++) {
	  e = open[n]->obc_e1[ee];
	  c = open[n]->obc_e2[ee];
	  e2 = geom->m2de[e];
	  c2 = geom->m2d[c];
	  x = geom->u1x[e2];
	  y = geom->u1y[e2];
	  z = geom->cellz[c] * master->Ds[c2];
	  open[n]->transfer_u2av[ee] = ramp * master->Hn1[c2] *
	    hd_ts_multifile_eval_xyz_by_name(open[n]->ntsfiles, 
					     open[n]->tsfiles,
					     open[n]->filenames, 
					     data->name,
					     master->t3d, x, y, z);
	}
      }
    }
  }
}

/* END bdry_eval_u1av_m()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to transfer the custom u1av boundary data to the slaves   */
/*-------------------------------------------------------------------*/
void bdry_transfer_u1av(master_t *master, geometry_t *window,
			window_t *windat)
{
  geometry_t *geom = master->geom;
  open_bdrys_t **open_m = geom->open;
  open_bdrys_t **open_w = window->open;
  int n;                        /* Counters */
  int e, ee;                            /* Sparse location / counter */

  /* Always do the custom transfers for 1 window; the window trasfer */
  /* vectors point to the master vectors for 1 window.  */
  for (n = 0; n < window->nobc; n++) {
    if (open_w[n]->bcond_nor2d & (FILEIN | CUSTOM)) {
      bdry_details_t *data = &open_w[n]->datau1av;
      /* Do the custom transfer if required */
      if (data->explct) {
	if(data->trans)
	  data->trans(master, open_w[n], data, window, windat);
	else {
	  for (ee = 1; ee <= open_w[n]->no2_e1; ee++)
	    open_w[n]->transfer_u1av[ee] = data->fill_value;
	}
      }
    }

    if (open_w[n]->bcond_tan2d & (FILEIN | CUSTOM)) {
      bdry_details_t *data = &open_w[n]->datau2av;
      /* Do the custom transfer if required */
      if (data->explct) {
	if(data->trans)
	  data->trans(master, open_w[n], data, window, windat);
	else {
	  for (ee = open_w[n]->no3_e1 + 1; ee <= open_w[n]->to2_e1; ee++)
	    open_w[n]->transfer_u2av[ee] = data->fill_value;
	}
      }
    }
  }

  /* Transfer data read from file from the master to the slave if */
  /* multiple windows exist.  */
  if (master->nwindows > 1) {
    for (n = 0; n < window->nobc; n++) {
      if (open_w[n]->bcond_nor2d & (FILEIN | CUSTOM)) {
        int mwn = open_w[n]->mwn;
        for (ee = 1; ee <= open_w[n]->no2_e1; ee++) {
          e = open_w[n]->tmap_u1[ee];
          open_w[n]->transfer_u1av[ee] = open_m[mwn]->transfer_u1av[e];
        }
      }
      if (open_w[n]->bcond_tan2d & (FILEIN | CUSTOM)) {
        int mwn = open_w[n]->mwn;
        for (ee = open_w[n]->no3_e1 + 1; ee <= open_w[n]->to2_e1; ee++) {
          e = open_w[n]->tmap_u2[ee];
          open_w[n]->transfer_u2av[ee] = open_m[mwn]->transfer_u2av[e];
        }
      }
    }
  }
}

/* END bdry_transfer_u1av()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets a no-gradient condition over open boundary ghost cells       */
/*-------------------------------------------------------------------*/
void OBC_bgz_nograd(geometry_t *window)
{
  window_t *windat = window->windat;
  int ee, c, cb, tn, n, m;

  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    if (open->bgz) {
      for (tn = 0; tn < open->ntr; tn++) {
	for (ee = 1; ee <= open->no3_e1; ee++) {
	  c = open->ogc_t[ee];
	  cb = open->obc_e2[ee];
	  for (m = 0; m < open->bgz; m++) {
	    windat->tr_wc[tn][c] = windat->tr_wc[tn][cb];
	    windat->dens[c] = windat->dens[cb];
	    if (open->options & (OP_UPSTRM) && tn == windat->sno) {
	      windat->tr_wc[tn][c] = 0.0;
	      windat->tr_wc[tn][c] = windat->tr_wc[tn][cb];
	    }
	    /*c = open->omape[ee][c];*/
	  }
	}
      }
    }
  }
}

void OBC_bgz_nogradb(geometry_t *window, open_bdrys_t *open, double *tr)
{

  int ee, c, cb, m;

  for (ee = 1; ee <= open->no3_e1; ee++) {
    c = open->ogc_t[ee];
    cb = open->obc_e2[ee];
    for (m = 0; m < open->bgz; m++) {
      tr[c] = tr[cb];
      c = open->omape[ee][c];
    }
  }
}

/* END OBC_bgz_nograd()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads in 3D blocks of data for flow relaxation and nudging        */
/*-------------------------------------------------------------------*/
void read_bdry_zone(master_t *master, open_bdrys_t *open, int ee, int mode)
{
  geometry_t *geom = master->geom;
  bdry_details_t *data;
  double *tvec, x, y, z;
  double ramp = (master->rampf & FILEIN) ? master->rampval : 1.0;
  int *obc;
  int tinc;
  int e2, c, c2, i, dir;
  int zone;

  if (mode & U1BDRY && mode & U1GEN) {
    data = &open->datau1;
    tvec = open->transfer_u1;
    tinc = open->no3_e1;
    zone = open->relax_zone_nor;    
    obc = open->obc_e1;
  } else if (mode & U1BDRY && mode & U2GEN) {
    data = &open->datau2;
    tvec = open->transfer_u2;
    tinc = open->to3_e1;
    obc = open->obc_e1;
    zone = open->relax_zone_tan;
  }
  c = open->obc_e2[ee];
  dir = open->ini[ee];
  for (i = 1; i < zone; i++) {
    c = open->nmape[ee][c];
    c2 = geom->m2d[c];
    e2 = geom->c2e[dir][c2];
    x = geom->u1x[e2];
    y = geom->u1y[e2];
    z = geom->cellz[c] * master->Ds[c2];
    tvec[ee + i * tinc] = ramp * master->Hn1[c2] *
      hd_ts_multifile_eval_xyz_by_name(open->ntsfiles, 
				       open->tsfiles,
				       open->filenames, data->name,
				       master->t3d, x, y, z);
  }
}

/* END read_bdry_zone()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Updates external tracer values via a mass balance                 */
/*-------------------------------------------------------------------*/	
double adjust_OBC_mass(geometry_t *window, 
		       window_t *windat, 
		       win_priv_t *wincon,
		       open_bdrys_t *open,  /* Boundary              */
		       int tn               /* Tracer number         */
		       )
{
  int c, cc, ee, e, ci, c2, tm, vc, *obc;
  double *vel, *dz, *hat;
  double d1, sgn = 1.0;
  double volin, volout, massin, massout, val;
  double tmax = 0.0;

  if (open->bcond_tra[tn] & CUSTOM) {
    dz = windat->dzu1;
    hat = window->h2au1;
    vel = windat->u1;
    vc = open->no3_e1;
    obc = open->obc_e1;

    /* Integrate the flow over the face (barotropic component) */
    tm = open->trm[tn];
    volin = volout = massin = massout = 0.0;
    /* Add estuarine volume if length specified */
    if (open->rlen) {
      int tm = open->trm[tn];
      for (cc = 1; cc <= open->no3_t; cc++) {
	c = open->obc_t[cc];
	d1 = open->rlen * dz[c] * hat[window->m2d[c]];
	volin += d1;
	massin += d1 * 0.5 * (windat->tr_wc[tn][c] + open->t_transfer[tm][cc]);
      }
    }
    for (ee = 1; ee <= vc; ee++) {
      e = obc[ee];
      ci = open->obc_e2[ee];
      c2 = window->m2d[ci];
      cc = open->e2c_e1[ee];
      sgn = (open->ceni[ee]) ? -1.0 : 1.0;
      if (window->gridz[ci] > windat->eta[c2]) continue;
      val = open->t_transfer[tm][cc];
      /*
      val = 0.5 * (open->t_transfer[tm][cc] +
      open->t_transfer[tm][open->t_imap[ee][1]]);
      */

      d1 = vel[e] * dz[e] * hat[window->m2de[e]];
      volin += d1;
      massin +=  val * fabs(d1);
      if (sgn * vel[e] < 0.0) {
	volout += fabs(d1);
	massout += windat->tr_wc[tn][ci] * fabs(d1);
	tmax = max(tmax, windat->tr_wc[tn][ci]);
      }
    }
    volin = fabs(volin);
    d1 = volin + volout;
    val = (d1) ? (massin + massout) / d1 : val;
    val = min(val, tmax);
    return(val);
  }
  return(0.0);
}

/* END adjust_OBC_mass()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Updates external salinity using methodology outlined in:          */
/* MacCready, P., Geyer, W.P. (2010) Advances in estuarine physics.  */
/* Annual Review of Marine Science, 2, 35 - 58.                      */
/* In particular, Eq. 19, with Eq. 16.                               */
/*-------------------------------------------------------------------*/	
double adjust_OBC_maccready(geometry_t *window, 
			    window_t *windat, 
			    win_priv_t *wincon,
			    open_bdrys_t *open,  /* Boundary         */
			    int tn               /* Tracer number    */
			    )
{
  int cc, e, ee, cs, cb, ci, c2, vc, *obc;
  double *vel, *dz, *hat;
  double d1, sgn = 1.0, vd;
  double sr, flow, depth, rho;
  double rho0 = 1024.0;         /* Standard reference density */
  double rhof = 1000.0;         /* Freshwater density */
  
  if (open->bcond_tra[tn] & CUSTOM) {
    dz = windat->dzu1;
    hat = window->h2au1;
    vel = windat->u1av;
    vc = open->no2_e1;
    obc = open->obc_e1;
    
    sr = 0.0;
    for (ee = 1; ee <= vc; ee++) {
      e = obc[ee];
      ci = cs = open->obc_e2[ee];
      cc = open->e2c_e1[ee];
      cb = open->nmape[ee][open->bot_t[cc]];
      c2 = window->m2d[ci];
      sgn = (open->ceni[ee]) ? -1.0 : 1.0;
      depth = windat->eta[cs] - window->botz[cs];
      rho = windat->dens_0[cb];
      rhof = (open->options & OP_IFRESH) ? 1000.0 : windat->dens_0[cs];

      /*
      flow = area = 0.0;
      while (c != window->zm1[c]) {
	if (window->gridz[ci] > windat->eta[c2]) continue;
	d1 = dz[c] * hat[window->m2d[c]];
	area += d1;
	flow += vel[c] * d1;
	c = window->zm1[c];
      }
      flow = (area) ? flow / area : 0.0;
      */
      flow = vel[e];

      /* Internal wave speed                                         */
      vd = sqrt(window->wincon->g * depth * (rho - rhof) / rho0);
      /* Salinity ratio, Eq. 19                                      */
      d1 = (vd) ? min(7.0 * pow(flow / vd, 2.0 / 3.0), 1.0) : 1.0;
      /* Outflow salinity                                            */
      sr += windat->sal[cb] * (1.0 - d1);
    }
    sr /= (double)vc;
    return(max(sr, 0.0));
  }
  return(0.0);
}

/* END adjust_OBC_maccready()                                        */
/*-------------------------------------------------------------------*/
