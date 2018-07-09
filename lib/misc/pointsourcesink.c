/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/misc/pointsourcesink.c
 *
 *  \brief Point source-sink code
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: pointsourcesink.c 5831 2018-06-26 23:48:06Z riz008 $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "pointsourcesink.h"

//int parse_location(pointsourcesink_t *p, char *line);
static void pss_locate(pss_t *p, double t);
static int parse_ss_location( pss_t *p, char *line);

//extern int verbose;

/** Reads a list of point sources and sinks from an ascii file.
  * The sources sinks are specified as shown below\n
  * pssnn means pss followed by the integer point source/sink number nn
  * (with the minimum number of digits needed).
  * Here, S is a string (which must not contain whitespace),
  * N is an integer and X, Y, Z and A are floating point
  * numbers.
  * 
  * This routine depends on the existence of a global variable:
  * verbose -   Sets level of messages printed
  * 
  * \# Number of point source/sinks
  * npointss     N
  *
  * \# Parameters for each point source/sink
  * \# Point source/sink name
  * pssnn.name    S
  *
  * \# Location ( x y z)
  * pssnn.location X Y Z
  *
  * \# Data - the next line is a time series definition
  * \# as used by my timeseries routines in sjwlib. The example
  * \# below assumes that the data is in an ascii or netCDF file.
  * pssnn.data  filename
  *
  * @param name Name of ascii file containing point source/sink list
  * @param t_units Time units to be used for time series
  * @param pss Returned pointer to point source/sink list
  * @param np Returned number of point sources/sinks
  * @param model_data Data for xyzijk and trI
  * @param xyzijk Pointer to function which converts (x,y,z) to (i,j,k)
  * @param trI Pointer to function which finds model index of tracer name
  */
void
pss_read(
char *name,
char *t_units,
pss_t **pss,
int *np,
void *model_data,
int (*xyzijk) (void *, double, double, double, int *, int *, int *),
int (*trI) (void *, char *)
)
{
    FILE *fp;
    int npss = 0;
    pss_t *p = NULL;
    int i, j;
    int n, t;
    int id_counter;
    
    /* Open the file */
    if( (fp=fopen(name,"r")) == NULL )
        quit("pss_read: Can't open %s\n",name);


    /* Get the number of point inputs */
    prm_read_int(fp,"npointss",&npss);
    
    /* Allocate memory for list of point inputs */
    if( npss > 0  && ((p=(pss_t *)malloc(npss*sizeof(pss_t))) == NULL) )
        quit("pss_read: Can't allocate memory for point source/sink list\n");
    
    /* Read each point input */
    for(i=0, j=0; i<npss; i++) {
        char key[MAXLINELEN];
	char buf[MAXLINELEN];
        
	/* Store index routine pointer and data needed by it */
	p[j].xyzijk = xyzijk;
	p[j].model_data = model_data;

        /* Name */
        sprintf(key,"pss%d.name",i);
        prm_read_char(fp,key,p[j].name);

        /* Location */
        sprintf(key,"pss%d.location",i);
	prm_read_char(fp,key,buf);
	if (!parse_ss_location(&p[j],buf))
	    continue;
	
        /* Data */
        sprintf(key,"pss%d.data",i);
	prm_read_char(fp,key,buf);
	p[j].ntsfiles = ts_multifile_read(buf, p[j].ts);

	/* Check data time units */
	for (n=0; n<p[j].ntsfiles; n++) {
	  if( strcmp(p[j].ts[n].t_units,t_units) != 0 ) {
	    ts_convert_time_units(&p[j].ts[n],t_units);
	  }
	}
	
	/* Find variable indices for each of the time series variables */
	p[j].nv  = ts_multifile_get_nv(p[j].ntsfiles, p[j].ts);
	p[j].tmap = i_alloc_1d(p[j].nv);
	p[j].watertsid = -1;
	id_counter = 0;
	for (t = 0; t < p[j].ntsfiles; t++) {
	  for(n=0; n<p[j].ts[t].nv; n++) {
	    if( strcasecmp(p[j].ts[t].varname[n],"water") == 0 ||
		strcasecmp(p[j].ts[t].varname[n],"flow")  == 0 ) {
	      p[j].watertsid = id_counter;
	    }
	    /* 
	     * Cache the tracer index to the master
	     */
	    p[j].tmap[id_counter] = (*trI) (model_data, p[j].ts[t].varname[n]);
	    if (
		p[j].tmap[id_counter] < 0 &&   /* invalid flux id */
		n != p[j].ts[t].ti        &&   /* and not time variable */
		id_counter != p[j].watertsid   /* and not the water/flow id */
		)
	      warn("pss_read: %s in %s in file %s isn't a model tracer!\n",
		   p[j].ts[t].varname[n], p[j].name, p[j].ts[t].name);
	    /* Increment the global (across tsfiles) id counter */
	    id_counter++;
	  }
	}
	j++;
    }
    
    npss = j;
    
    /* Close the file */
    fclose(fp);
    
    /* Store pointer to list of point source/sinks */
    *pss = p;
    *np = npss;
}

static void pss_locate(pss_t *p, double t)
{
    /* Check if the source/sink has a time varying location */
    if( p->loc == NULL )
        return;
    
    /* Get x, y and z values */
    p->x = ts_eval(p->loc, p->x_id, t);
    p->y = ts_eval(p->loc, p->y_id, t);
    p->zlow = ts_eval(p->loc, p->zl_id, t);
    p->zhigh = ts_eval(p->loc, p->zh_id, t);
    p->z = (p->zlow+p->zhigh)/2;

    /* Convert to model indices */
    if( (*(p->xyzijk))(p->model_data, p->x,p->y,p->z,
		       &p->e1[0],&p->e2[0],&p->e3[0]) < 0 )
	quit("pss_locate: %s location (%.10g,%.10g,%.10g) can't be converted to indices\n",p->name,p->x,p->y,p->z);
}

static int parse_ss_location( pss_t *p, char *line)
{
    int e1, e2, e3;
    FILE *fp;

    /* Clear location time series pointer */
    p->loc = NULL;

    /* Time independent, range of z values */
    if( sscanf(line,"%lf %lf %lf %lf",&p->x,&p->y,&p->zlow,&p->zhigh) == 4 ) {
        if( p->zlow >= p->zhigh )
	    quit("parse_ss_location: %s has bad z range (must be low then high)\n",p->name);
	p->z = (p->zlow+p->zhigh)/2;
	if( (*(p->xyzijk))(p->model_data, p->x,p->y,p->z,&e1,&e2,&e3) < 0 ) {
	    warn("parse_ss_location: %s location (%.10g,%.10g,%.10g) can't be converted to indices\n",p->name,p->x,p->y,p->z);
	    return 0;
	}
	p->vc = 1;
	p->e1 = i_alloc_1d(p->vc);
	p->e2 = i_alloc_1d(p->vc);
	p->e3 = i_alloc_1d(p->vc);
	p->e1[0] = e1;
	p->e2[0] = e2;
	p->e3[0] = e3;
    }
    /* Time independent, single  z value */
    else if( sscanf(line,"%lf %lf %lf",&p->x,&p->y,&p->z) == 3 ) {
		if( (*(p->xyzijk))(p->model_data, p->x,p->y,p->z,&e1,&e2,&e3) < 0 ) {
		    warn("readPointSourceSink: %s location (%.10g,%.10g,%.10g) can't be converted to indices\n",p->name,p->x,p->y,p->z);
		    return 0;
		}
		p->vc = 1;
		p->e1 = i_alloc_1d(p->vc);
		p->e2 = i_alloc_1d(p->vc);
		p->e3 = i_alloc_1d(p->vc);
		p->e1[0] = e1;
		p->e2[0] = e2;
		p->e3[0] = e3;
		p->zlow = p->z;
		p->zhigh = p->z;

    }
    /* Time dependent, all coords from a time series */
    else if( (fp = fopen("line", "+r")) != NULL ) {
		fclose(fp);
	
		if( (p->loc = (timeseries_t*) malloc(sizeof(timeseries_t))) == NULL )
		     quit("readPointSourceSink: No memory for location time series\n");
	        ts_read(line,p->loc);
		/* Get coordinate indices */
		if( (p->x_id = ts_get_index(p->loc, "x")) < 0 && (p->x_id = ts_get_index(p->loc, "X")) < 0 )
		    quit("readPointSourceSink: Location timeseries must have a variable named x (or X)\n");
		if( (p->y_id = ts_get_index(p->loc, "y")) < 0 && (p->y_id = ts_get_index(p->loc, "Y")) < 0 )
		    quit("readPointSourceSink: Location timeseries must have a variable named y (or Y)\n");
		if( (p->zl_id = ts_get_index(p->loc, "z_low")) < 0 || (p->zh_id = ts_get_index(p->loc, "z_high")) < 0 )
		    quit("readPointSourceSink: Location timeseries must have variables named z_low and z_high\n");
		
		/* Set location at some arbitrary initial time.
		 * Probably don't really need to do this here, as
		 * any routines using the structure should call
		 * the location routine themselves.
		 */
		pss_locate(p, 0.0);
    }
    else {
		quit("readPointSourceSink: Can't understand location definition for %s\n",p->name);
    }

    return 1;
}
