/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/particles/pl.c
 *  
 *  Description:
 *  plastic particle tracking code for meco
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: pt.c 7241 2022-10-25 00:26:55Z her127 $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <netcdf.h>
#include "hd.h"

char *MACROPCLASS[] = {
  "Small_Pl",
  "Large_PL",
  "Bottle",
  "Bag",
  "Foam"
};

#define NUM_MACRO_P (int)(sizeof(MACROPCLASS)/sizeof(char *))

struct {
  char *name;
  char *ptype;
  double diam;   /* Mean diameter (m) */
  double dens;   /* Density (kg/m) */
  double srat;   /* Surface ratio above/below water */
  int type;
} macro_plastic[] = {
  {
    "small_pl",
    "PET:50 PE:50",
    0.01,
    960.0,
    0.0,
    P_SMALL},
  {
    "large_pl",
    "PET:50 PE:50",
    0.1,
    960.0,
    0.0,
    P_LARGE},
  {
    "bottle",
    "PET",
    0.2,
    920.0,
    1.0,
    P_BOT},
  {
    "bag",
    "PE",
    0.0067,
    920.0,
    0.0,
    P_BAG},
  {
    "foam",
    "PS",
    0.2,
    330.0,
    3.0,
    P_FOAM},
  {"", "", 0, 0, 0, 0}
};

struct {
  char *name;
  char *long_name;
  double cc;            /* Carbon content */
  double dens;          /* Density: from Stubbins et al., Science 373, 5155 (2021) */
  int type;
} plastic_type[] = {
  {
    "PE",
    "polyethylene",
    0.86,
    945.0,
    P_PE},
  {
    "PP",
    "polypropylene",
    0.86,
    895.0,
    P_PP},
  {
    "PET",
    "polyethylene terephthalate",
    0.83,
    1350.0,
    P_PET},
  {
    "PVC",
    "polyvinyl chloride",
    0.38,
    1450.0,
    P_PVC},
  {
    "PS",
    "polystyrene",
    0.92,
    330.0,
    P_PS},
  {"", "", 0, 0, 0}
};

/*
typedef struct {
  double e1;
  double e2;
  double e3;
  int c;
  short flag;
  int dumpf;
  double age;
  unsigned char out_age;
  double size;
  unsigned char out_size;
  double svel;
} particle_t;

typedef struct {
  char name[MAXSTRLEN];
  char ptype[MAXSTRLEN];
  int type;
  double diam;
  double dens;
  double srat;

} plastic_t;

*/



/*-------------------------------------------------------------------*/
/* Sets up macroplastics type structures                             */
/*-------------------------------------------------------------------*/
void set_macrop_class(master_t *master) {
  parameters_t *params = master->params;
  FILE *fp = master->prmfd;
  char key[MAXSTRLEN], buf[MAXSTRLEN], val[MAXSTRLEN];
  int n, i, ic, ip, tn, ipr, nt, found;
  int nmacrop = 0;
  int nptype = 0;
  char *fields[MAXSTRLEN * MAXNUMARGS];
  double d1;

  /* Count the number of plastic types                               */
  for (i = 0; strlen(plastic_type[i].name); ++i) nptype++;
  for (i = 0; strlen(macro_plastic[i].name); i++) nmacrop++;

  if (strlen(params->cmacrop) == 0) {
    hd_warn("Macro-plastic types not defined.\n");
    return;
  }
  master->nmacrop = parseline(params->cmacrop, fields, MAXNUMARGS);

  /* Set up and populate the master macro-plastic structure          */
  master->macrop = (plastic_t **)malloc(sizeof(plastic_t *) * master->nmacrop);
  for (n = 0; n < master->nmacrop; n++) {
    master->macrop[n] = (plastic_t *)malloc(sizeof(plastic_t));
    memset(master->macrop[n], 0, sizeof(plastic_t));
    master->macrop[n]->ptype = i_alloc_1d(nptype);
    master->macrop[n]->ptr = d_alloc_1d(nptype);
    master->macrop[n]->microtr = i_alloc_1d(nptype);
    strcpy(master->macrop[n]->name, fields[n]);

    /* Find the specified plastic type in supported types            */
    found = 0;
    for (ic = 0; ic < nmacrop; ic++) {
      if (strcmp(master->macrop[n]->name, macro_plastic[ic].name) == 0) {
	found = 1;
	break;
      }
    }
    if (!found) hd_quit("Can't find plastic type %s in supported types.\n", 
			master->macrop[n]->name);

    master->macrop[n]->diam = macro_plastic[ic].diam;
    master->macrop[n]->dens = macro_plastic[ic].dens;
    master->macrop[n]->srat = macro_plastic[ic].srat;
    master->macrop[n]->type = macro_plastic[ic].type;

    /* Plastic type */
    found = ipr = 0;
    for (i = 0; i < nptype; i++) {
      master->macrop[n]->ptype[i] = -1;
      master->macrop[n]->ptr[i] = 0.0;
      master->macrop[n]->microtr[i] = -1;
      if ((strcmp(macro_plastic[ic].ptype, plastic_type[i].name) == 0) ||
	  (nt = decode_tag(macro_plastic[ic].ptype, plastic_type[i].name, val))) {
	master->macrop[n]->ptype[ipr] = i;
	if (nt == 1) {
	  d1 = atof(val);
	  master->macrop[n]->ptr[ipr] = (d1 > 1.0) ? d1 / 100.0 : d1;
	} else
	  master->macrop[n]->ptr[ipr] = 1.0;
	/* Set sediment tracers for micro-plastics                   */
	if ((tn = tracer_find_index(plastic_type[i].name, master->ntr, 
				    master->trinfo_3d)) >= 0) {
	  master->macrop[n]->microtr[ipr] = tn;
	} else 
	  hd_quit("Can't find plastic type %s for macro-plastic %s.\n",
		  plastic_type[i].name, master->macrop[n]->name); 
	ipr++;
      }
    }
  }

  /* Allocate particle sources to macro-plastic types                */
  if (master->pt_nsource > 0) {
    master->pt_ptype = i_alloc_1d(master->pt_nsource);

    master->pt_sizelim = 0.0;
    for (n = 0; n < master->pt_nsource; ++n) {
      sprintf(key, "PT_Source%d.Plastic",n);
      if (prm_read_char(fp, key, buf)) {
	master->pt_stype[n] &= ~PT_STANDARD;
	master->pt_stype[n] |= (PT_PLASTIC|PT_SVEL);
	master->do_pt |= PT_PLASTIC;
	for (i = 0; i < master->nmacrop; ++i) {
	  if (strcmp(buf, master->macrop[i]->name) == 0) {
	    master->pt_ptype[n] = i;
	    break;
	  }
	}
	if (i == master->nmacrop) hd_quit("Can't find plastic type %s\n", buf);
	ip = i;

	sprintf(key, "PT_Source%d.Mass", n);
	if (prm_read_double(fp, key, &master->pt_size[n]))
	  master->pt_sizelim = max(master->pt_sizelim, master->pt_size[n]);
	else
	  hd_warn("pt_params_init: Macro-plastic source%d mass is zero\n", n);
	sprintf(key, "PT_Source%d.Loss", n);
	if (!prm_read_double(fp, key, &master->pt_decay[n]))
	  hd_warn("pt_params_init: Macro-plastic source%d loss is zero\n", n);

	master->svel_type[n] = PLSTOKES;
	master->pt_dens[n] = master->macrop[ip]->dens;
	/*master->pt_dens[n] = plastic_type[master->macrop[ip]->ptype].dens;*/
      }
    }
    master->pt_dumpf |= PT_SIZE;
  }
}

/* END set_macrop_class()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
double plas_init(master_t *master)
{
  int n, sm, class, type;

  for (n = 0; n < master->ptn; n++) {
    particle_t *p = &master->ptp[n];
    sm = master->pt_sm[n];          /* Particle source map          */
    class = master->pt_ptype[sm];
    type = master->macrop[class]->type;
  }
}

/*-------------------------------------------------------------------*/
