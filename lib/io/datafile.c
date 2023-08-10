/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/io/datafile.c
 *
 *  \brief Functions for reading, writing and interpolating datafiles
 *
 *  The data may have associated coordinate variable and be multi-dimensional
 *   
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: datafile.c 7292 2023-02-20 03:47:10Z riz008 $
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "netcdf.h"
#include "ems.h"

static char *lon_names[7] = {
  "degrees_east",
  "degree_east",
  "degree_E",
  "degrees_E",
  "degreeE",
  "degreesE",
  NULL
};

static char *lat_names[7] = {
  "degrees_north",
  "degree_north",
  "degree_N",
  "degrees_N",
  "degreeN",
  "degreesN",
  NULL
};

static char *v_names[4] = {
  "metres",
  "meters",
  "m",
  NULL
};

static char *sst_names[5] = {
  "analysed_sst",
  "analysis_error",
  "sea_surface_temperature",
  "quality_level",
  NULL
};

/* Prototypes for local routines */
void netcdf_read(int fid, datafile_t *df, int type);
void netcdf_read_records(datafile_t *df, int varid);
void netcdf_read_data(datafile_t *df, df_variable_t *v, int rec, int roffset);
void netcdf_free(datafile_t *df);
void netcdf_read_attrib(datafile_t *df, int varid, int attnum,
                        df_attribute_t *a);
void multi_netcdf_read(FILE *fp, datafile_t *df);
void multi_netcdf_read_records(datafile_t *df, int varid);
void multi_netcdf_read_data(datafile_t *df, df_variable_t *v, int rec);
void multi_netcdf_free(datafile_t *df);
void multi_netcdf_read_attrib(datafile_t *df, int varid, int attnum,
                        df_attribute_t *a);
void multi_netcdf_read_atts(datafile_t *df, df_multi_file_t *fd);
void ascii_read(FILE * fp, datafile_t *df);
void ascii_write(FILE * fp, datafile_t *df);
void ascii_free(datafile_t *df);

void alloc_data_records(datafile_t *df, df_variable_t *v, int nrecs);
void free_data_records(datafile_t *df, df_variable_t *v);
double *alloc_data(datafile_t *df, df_variable_t *v);
void df_free_data(datafile_t *df, df_variable_t *v, double *data);
/*void read_data(datafile_t *df, df_variable_t *v, int rec);*/
void left_shift_data(datafile_t *df, df_variable_t *v, int nrecs);
void right_shift_data(datafile_t *df, df_variable_t *v, int nrecs);

static void destroy_hq_within_ht(void *ht);
static int has_name(char *name, char *array[]);

/* Located in dfcoords.c */
extern void decode_coords(datafile_t *df, df_variable_t *v, char *coords);
extern void decode_coord_type(char *ctype, VariableType * type,
                              char **coord_domain);
void free_coord_system(datafile_t *df, df_coord_system_t *csystem);
void df_set_record(datafile_t *df, int varid);
int df_reset_record(datafile_t *df);

extern int df_default_geotype;
extern char *df_default_projection;

/* MH mempack */
void mempack_read(char *name, datafile_t *df, int type);
void mempack_read_records(datafile_t *df, int varid);
void mempack_read_data(datafile_t *df, df_variable_t *v, int r, int roffset);
int mempack_read_file(char *mpname, datafile_t *df);
void mempack_free(datafile_t *df);
int MPK_SVARS = 4;

/** Read and set up a datafile.
  *
  * There are several ways of defining the data file.
  * 
  * 1. Ascii file of values stored in a column array.
  * If the file is an ascii file, it must have the
  * following format:
  * 
\verbatim
# Comments
## COLUMNS n
## 
## COLUMN1.name XXXX
## COLUMN1.long_name XXXX
## COLUMN1.units XXXX
## COLUMN1.missing_value XXXX
## COLUMN1.fill_value XXXX
##
## COLUMN2.name XXXX
## COLUMN2.long_name XXXX
## COLUMN2.units XXXX
## COLUMN2.missing_value XXXX
## COLUMN2.fill_value XXXX
##
 .
 .
 .
##
v   v   v   v   ...
v   v   v   v   ...
v   v   v   v   ...
 .
 .
 .
\endverbatim
  * 
  * 2. netCDF file with 0, 1 or 2 spatial
  * dimensions. The file may support contain one
  * unlimited dimension. This dimension is handled
  * specially, and it is assumed that the values
  * associated with the dimension monotonically
  * increase/decrease.
  * 
  * @param name file containing the data file data.
  * @param df pointer to datafile structure (assumed
  * 	 not previously initialised).
  * 
  * @see quit() If anything goes wrong.
  */
void df_read(char *name, datafile_t *df)
{
  int fid;
  FILE *fp;

  /* Clear the structure */
  memset(df, 0, sizeof(datafile_t));

  /* Store the name */
  if ((df->name = (char *)malloc(strlen(name) + 1)) == NULL)
    quit("df_read: Can't allocate memory for name\n");
  strcpy(df->name, name);

  /* Try to open the file */
  if (nc_open(name, NC_NOWRITE, &fid) == NC_NOERR) {
    /* It is a netCDF file! */
    netcdf_read(fid, df, DFT_NETCDF);
  }

  else if ((fp = fopen(name, "r")) != NULL) {

    /* Read the first line. If the file is of the type
     * multi-netcdf then process, else treat as an ascii file. */
    prm_set_errfn(quiet);
    if (prm_skip_to_end_of_key(fp, "multi-netcdf-version")) {
       prm_flush_line(fp);
       multi_netcdf_read(fp, df);       
    } else if (prm_skip_to_end_of_key(fp, "mempack-version")) { /* MH mempack */
       prm_flush_line(fp);
       mempack_read(name, df, DFT_MEMPACK);
       return;
    } else
       ascii_read(fp, df);
  }
    
  else if (endswith(name, ".mpk")) {
    mempack_read(name, df, DFT_MEMPACK);
    return;
  }

  else
    quit("df_read: Can't open %s\n", name);

  if (df->records != NULL)
    df_check_records(df);

}


/** Allocates memory for the datafile structure
 */
datafile_t *df_alloc(void)
{
  datafile_t *df = (datafile_t*)malloc(sizeof(datafile_t));
  if (df == NULL)
    quit("df_alloc: Memory allocation error\n");

  return(df);
}


/** Set the record variable. This can be quite expensive as it
  * requires all variables that depend on the record dimension
  * to be adjusted.
  *
  * The record variable MUST be one dimensional, and the
  * record dimension must be the first dimension of ALL
  * variables in the datafile_t that use this dimension.
  *
  * @param df pointer to datafile.
  * @param varid variable index/identifier.
  */
void df_set_record(datafile_t *df, int varid)
{
  df_variable_t *rv = NULL;
  int recdimid = -1;
  int i;

  /* Check if this is a valid variable */
  if ((varid < 0) && (varid >= df->nv))
    quit("df_set_record: Invalid variable identifier.\n");
  else
    rv = &df->variables[varid];

  /* Confirm that the variable contains only one dimension. */
  if ((rv->nd == 1) && (df->dimensions[rv->dimids[0]].size > 0))
    recdimid = rv->dimids[0];
  else
    quit
      ("df_set_record: Record variable must contain single non-zero dimension.\n");

  /* Read in the record data and populate the record info. and record
     variable in the datafile_t structure. */
  if (df->type == DFT_NETCDF) {
    netcdf_read_records(df, varid);
  } else if (df->type == DFT_MULTI_NETCDF) {
    multi_netcdf_read_records(df, varid);
  } else if (df->type == DFT_MEMPACK) {
    mempack_read_records(df, varid);
  }

  rv->start_record = 0;
  rv->nrecords = df->dimensions[recdimid].size;


  /* Copy over the record information. */
  df->records = rv->data;
  df->nrecords = (int)df->dimensions[recdimid].size;
  df->rec_name = rv->name;
  df->rec_units = rv->units;
  df->rec_longname = rv->longname;
  df->ri = varid;
  df->geotype = df_default_geotype;
  df->projection = df_default_projection;
  df->t0 = df->records[0];
  df->r0 = (df->rec_modulus) ? df->nrecords - 1 : 0;
  df->r1 = 0;
  df->frac = 0.0;

  /* Sweep through all of the variables and reduce any that contain the
     recdimid as the first dimension. Nullify the data. */
  for (i = 0; i < df->nv; ++i) {
    df_variable_t *v = &df->variables[i];
    if (v->dimids && v->dimids[0] == recdimid) {
      int i;

      v->dim_as_record = 1;
      if (v->data != NULL) {
        v->start_record = 0;
        v->nrecords = df->nrecords;
      } else {
        v->start_record = 0;
        v->nrecords = 0;
      }

      --(v->nd);
      for (i = 0; i < v->nd; ++i)
        v->dimids[i] = v->dimids[i + 1];

    } else
      v->dim_as_record = 0;
  }

  /* Check for modulo */
  df->rec_modulus = 0;
  if (df->ri >= 0) {
    df_variable_t *tv = &df->variables[df->ri];
    for (i = 0; i < tv->na; ++i) {
      df_attribute_t *a = &tv->attributes[i];
      if ((strcasecmp(a->name, "modulo") == 0) && CHK_TYPE(a, AT_TEXT)) {
        df->rec_modulus = 1;
        df->rec_mod_scale = atof(ATT_TEXT(a));
        break;
      }
    }
  }
}


/** Re-set the time records.
  * The record variable MUST be one dimensional, and the
  * record dimension must be the first dimension of ALL
  * variables in the datafile_t that use this dimension.
  *
  * @param df pointer to datafile.
  */
int df_reset_record(datafile_t *df)
{
  df_variable_t *rv = &df->variables[df->ri];
  int recdimid = rv->dimids[0];
  size_t recsize;
  int doread = 0;
  int i, fid;

  /* Read in the record data and populate the record info. and record
     variable in the datafile_t structure. */
  if (df->type == DFT_NETCDF) {
    if (nc_open(df->name, NC_NOWRITE, &fid) == NC_NOERR) {
      nc_close(df->ncid);
      df->ncid = fid;
      nc_inq_dimlen(df->ncid, df->dimensions[recdimid].dimid, &recsize);
      if (recsize > df->nrecords) {
	df->dimensions[recdimid].size = recsize;
	netcdf_read_records(df, df->ri);
	doread = 1;
      }
    }
  } else if (df->type == DFT_MULTI_NETCDF) {
    /*
     * Not implemented yet. Need to:
     *  o) Redo multi_netcdf_read to see if more files have been added
     *  o) And then check the total records
     */
    doread = 0;
  }

  if (doread) {
    rv->start_record = 0;
    rv->nrecords = df->dimensions[recdimid].size;
    df->records = rv->data;
    df->nrecords = (int)df->dimensions[recdimid].size;
    if (df->rec_mod_scale != 1.0) {
      for (i = 0; i < df->nrecords; i++)
	df->records[i] *= df->rec_mod_scale;
    }
    warn("df_reset_record: re-read time records to %.1f %s\n", 
	 df->records[df->nrecords - 1], df->rec_units);
    return(0);
  }
  return(1);
}


/**
  * Find record indices which bracket a requested record value.
  * 
  * @param df pointer to datafile structure
  * @param r specified sample record
  * @param before index value just before record specified
  * @param after index value just after record specified
  * @param frac fraction of record interval
  * @return 1 for interpolated value, 0 for extrapolated
  */
int df_find_record(datafile_t *df, double r, int *before, int *after,
                   double *frac)
{
  int imid;
  int ilow = 0;
  int ihigh = df->nrecords - 1;
  long iclk, clk, clock;
  int f1, f2, n = 0;
  int wait_inc = 60;

  if (df->records == NULL) {
    *before = df->r0 = 0;
    *after = df->r1 = 0;
    *frac = df->frac = 0.0;
    return (0);
  }

  /* If the time has not changed since the last function call, then
     return the previously calculated bounds and fraction. */
  if (r == df->t0) {
    *before = df->r0;
    *after = df->r1;
    *frac = df->frac;
    df->t0 = r;
    return (0);
  }

  if (df->rec_modulus) {
    while (r < 0)
      r += df->rec_mod_scale;
    r = fmod(r, df->rec_mod_scale);
  }

  if (df->type == DFT_MEMPACK) {
    df_mempack_t *mp = (df_mempack_t *)df->private_data;
    if (mp->runcode == DF_FAIL) {
      *before = df->r0;
      *after = df->r1;
      *frac = df->frac;
      df->runcode = DF_FAIL;
      df->t0 = r;
    } else
      df->runcode = DF_RUN;
    mempack_read_records(df, df->ri);
    if (strcmp(mp->tunits, df->rec_units) != 0)
      quit("df_find_record: mempack units must be %s\n", df->rec_units);
    /*    
    if (r < mp->next)
      printf("using file %s @ t=%f days mp->time=%f next=%f\n", df->name, r / 86400.0, mp->time/86400, mp->next/86400);
    */
    if (r >= mp->next) {
      /*printf("waiting for data: file %s @ %f days\n", df->name, r / 86400.0);*/
      /*printf("waiting for data: file %s @ %f, time = %f next = %f\n", df->name, r/86400, mp->time/86400, mp->next/86400);*/
      /*warn("waiting for data: file %s @ %f days\n", df->name, r / 86400.0);*/
      f1 = f2 = 1;
      time(&iclk);
      while (r >= mp->next) {
	mempack_read_records(df, df->ri);
	mp->next = df->records[ihigh] + mp->dt;
	clock = time(&clk) - iclk;
	if (clock == 1 && f1) {
	  warn("df_find_record: waiting (%ds) for data: file %s @ %f days\n", clock, df->name, r / 86400.0);
	  f1 = 0;
	}
	if (clock > 1 && clock%wait_inc == 0 && f2) {
	  warn("df_find_record: still waiting (%ds) for data: file %s @ %f days\n", clock, df->name, r / 86400.0);
	  f2 = 0;
	}
	if (clock%wait_inc >= 1 && !f2) f2 = 1;
	n++;
      }
      /*printf("using wait file %s @ t=%f days mp->time=%f next=%f\n", df->name, r / 86400.0, mp->time/86400, mp->next/86400);*/
      /*printf("Read record from %s at %f : %f. Continue?", df->name, df->records[ihigh]/86400);
      scanf("%d", &imid);*/
    } else {
      /*printf("File %s\n", df->name);
	printf("Continue : model time %f <= record time %f\n",r/86400,(df->records[ihigh]+mp->dt)/86400);*/
    }
  }

  /* first check whether t is within the table range */
  if (r <= df->records[ilow]) {
    if (df->rec_modulus) {
      /* Here we are interpolating for the wrap around */
      double nl = df->records[ihigh] - df->rec_mod_scale;
      double dl = df->records[ilow] - nl;
      *before = df->r0 = ihigh;
      *after = df->r1 = ilow;
      *frac = 0.0;
      if (dl)
	df->frac = (r - nl) / dl;
      df->t0 = r;
      return 1;
    } else {
      *before = df->r0 = ilow;
      *after = df->r1 = ilow;
      *frac = df->frac = 0.0;
      df->t0 = r;
      return (0);
    }
  }
  if (r >= df->records[ihigh]) {
    if (df->rec_modulus) {
      double nh = df->records[ilow] + df->rec_mod_scale;
      *before = df->r0 = ihigh;
      *after = df->r1 = ilow;
      *frac = df->frac = (r - df->records[ihigh]) / (nh - df->records[ihigh]);
      df->t0 = r;
      return 1;
    } else {
      if (df_reset_record(df)) {
	*before = df->r0 = ihigh;
	*after = df->r1 = ihigh;
	*frac = df->frac = 0.0;
	df->t0 = r;
	if (df->type == DFT_MULTI_NETCDF)
	  df->t0 = -1;

	return (0);
      } else {
	/* df->nrecords have changed so update value */
	ihigh = df->nrecords - 1;
      }
    }
  }

  /* perform binary chop to determine values either side of t */
  while (ihigh - ilow > 1) {
    imid = (ilow + ihigh) / 2;
    if (r >= df->records[imid])
      ilow = imid;
    else
      ihigh = imid;
  }

  /* Store results and return */
  *before = df->r0 = ilow;
  *after = df->r1 = ihigh;
  *frac = df->frac =
    (r - df->records[ilow]) / (df->records[ihigh] - df->records[ilow]);

  df->t0 = r;
  return (1);
}


/**
  * Read a region of data records into a memory buffer.
  *
  * If the variable has not been read before, then memory
  * will be allocated.
  * 
  * If the variable data has already been read then the
  * function will not re-read the data.
  *
  * If a portion of the specified record region has already
  * been read, then the buffers will be adjusted, and only the
  * new data read in.
  * 
  * @param df pointer to datafile structure
  * @param v pointer to variable 
  * @param start_rec record number (negative if var has not record info).
  * @param nrecs record number (must be 1 if start_rec is negative).
  */
void df_read_records(datafile_t *df, df_variable_t *v, int start_rec,
                     int nrecs)
{
  int i, n;

  if (df->type == DFT_ASCII)
    return;

  /* Is there a coordinate system. If so, then for a read of all of the
     variables in the coordinate system. */
  if (v->csystem != NULL) {
    df_coord_system_t *cs = v->csystem;
    for (i = 0; i < cs->nc; ++i) {
      df_variable_t *cv = &df->variables[cs->coordids[i]];
      if (cv->dim_as_record)
        df_read_records(df, cv, start_rec, nrecs);
      else
        df_read_records(df, cv, 0, 1);
    }
  }

  if (df->type == DFT_MEMPACK) {
    df_mempack_t *mp = (df_mempack_t *)df->private_data;
    if (v->data == NULL) {
      alloc_data_records(df, v, df->dimensions[df->ri].size);
      mempack_read_data(df, v, 0, 0);
      mempack_read_data(df, v, 1, 0);
    }
    if (mp->status[v->varid] & DF_GO) {
      left_shift_data(df, v, 1);
      mempack_read_data(df, v, 1, 0);
      mp->status[v->varid] = DF_WAIT;
    }
    return;
  }

  /* Check whether data has been allocated. If one of the dimension is a
     record, check that the size is the same as requested. If not clear
     the data, it's too much effort to regig. */
  if ((v->data != NULL) && (v->dim_as_record)) {
    if (v->nrecords != nrecs) {
      free_data_records(df, v);
      v->data = NULL;
      v->start_record = 0;
      v->nrecords = 0;
    }
  }

  /* Check whether the record range is appropriate for the specified
     variable. Quit if incorrect values were passed through. */
  if (start_rec >= 0) {
    if (v->dim_as_record) {
      if (!(df->rec_modulus) && (start_rec >= df->nrecords))
        quit
          ("df_read_records: Start buffer record (%d) out of bounds for variable '%s'.\n",
           start_rec, v->name);

      if (!(df->rec_modulus) && ((start_rec + nrecs - 1) >= df->nrecords))
        quit
          ("df_read_records: End buffer record (%d) out of bounds for variable '%s'.\n",
           start_rec + nrecs - 1, v->name);
    } else if (nrecs != 1)
      quit("df_read_records: nrecs should be 1 for variable '%s' not %d\n",
           v->name, nrecs);
  }

  else
    quit("df_read_records: start_rec was negative.");

  /* Check to see whether the data has already been read for the range
     specified. If so, then ignore. */
  if (!((v->data != NULL) && (v->start_record == start_rec)
        && (v->nrecords == nrecs))) {
    int sr = start_rec;
    int er = start_rec + nrecs - 1;

    /* Define the range (sr and er) which need to be read from the file.
       By default read the whole region specified. */
    if ((v->data == NULL) || (!v->dim_as_record)) {

      /* df_variable_t with non-record dimension or NULL data. Don't need
         to worry about swapping. */
      if (v->data == NULL) {
        alloc_data_records(df, v, nrecs);
      }
    } else {
      /* df_variable_t with record dimension and previously * allocated
         memory. If new data is entirely * outside of the existing range,
         then re-read. If overlap * then swap buffers and read
         non-overlapping regions. */
      int csr = v->start_record;
      int cer = v->start_record + v->nrecords - 1;
      if (!((er < csr) || (sr > cer))) {
        if (sr < csr) {
          /* Overlap at start - Shift right */
          right_shift_data(df, v, csr - sr);
          er = csr - 1;
        }

        else {
          /* Overlap at end - Shift left */
          left_shift_data(df, v, sr - csr);
          sr = cer + 1;
        }
      }
    }

    v->start_record = start_rec;
    v->nrecords = nrecs;
    for (i = sr; i <= er; ++i) {

      if (df->type == DFT_NETCDF) 
         netcdf_read_data(df, v, i, 0);

      else if (df->type == DFT_MULTI_NETCDF)
         multi_netcdf_read_data(df, v, i);
    }
  }
}


/**
  * Free the memory associated with the specified variable.
  *
  * @param df pointer to data file structure.
  * @param varid The variable id.
  */
void df_flush_variable(datafile_t *df, int varid) {

  df_variable_t *v = df_get_variable(df, varid);

  if (v != NULL) {
      if (v->data != NULL) {
         free_data_records(df, v);
         v->data = NULL;
      }
  }
}

/**
  * Free memory associated with a data file.
  * 
  * @param df pointer to data file structure.
  * 
  * @see quit() If anything goes wrong.
  */
void df_free(datafile_t *df)
{
  if (df != NULL) {
    int i, j;

    /* Free any memory specifically associated with the file type */
    if (df->type == DFT_ASCII)
      ascii_free(df);
    else
      netcdf_free(df);

    /* Deallocate memory associated with variables */
    if (df->variables != NULL) {
     
      for (i = 0; i < df->nv; ++i) {
        df_variable_t *v = &df->variables[i]; 
        free(v->name);

        if (v->longname != NULL)
          free(v->longname);

        if (v->units != NULL)
          free(v->units);

        if (v->ncdimids != NULL) {
          free(v->ncdimids);
        }

        if (v->coord_domain != NULL)
          free(v->coord_domain);

        if (v->coordids != NULL)
          free(v->coordids);

        if (v->attributes != NULL) {        
          for (j = 0; j < v->na; ++j)
          {
            free(v->attributes[j].name);
          }        
          free(v->attributes);
        }

        if (v->csystem != NULL)
          free_coord_system(df, v->csystem);
         

        if (v->dimids != NULL)
          free(v->dimids);

	if (v->gs0_master_2d != NULL) {
	  GRID_SPECS *gs = v->gs0_master_2d;
	  grid_specs_destroy(gs);
	}
	if (v->gs1_master_2d != NULL) {
	  GRID_SPECS *gs = v->gs1_master_2d;
	  delaunay *d = gs->d;
	  point *p = d->points;
	  free((point *)p);
	  delaunay_destroy(d);
	  grid_specs_destroy(gs);
	}
	if (v->gs0_master_3d != NULL) {
	  GRID_SPECS **gs = v->gs0_master_3d;
	  int nz = v->nz, k;
	  for (k = 0; k < nz; k++) {
	    if (gs[k] != NULL) {
	      gs[k]->destroy_pbathy = 0;
	      grid_specs_destroy(gs[k]);
	    }
	  }
	  free((GRID_SPECS **)gs);
	}
	if (v->gs1_master_3d != NULL) {
	  GRID_SPECS **gs = v->gs1_master_3d;
	  int nz = v->nz, k;
	  for (k = 0; k < nz; k++) {
	    if (gs[k] != NULL) {
	      delaunay *d = gs[k]->d;
	      point *p = d->points;
	      for (i = 0; i < v->kn[k]; i++) d_free_1d(p[i].v);
	      free_1d((point *)p);
	      delaunay_destroy(d);
	      gs[k]->destroy_pbathy = 0;
	      grid_specs_destroy(gs[k]);
	    }
	  }
	  free((GRID_SPECS **)gs);
	  i_free_1d(v->kn);
	  i_free_2d(v->kmap);
	}

        if (v->data != NULL)
          free_data_records(df, v);      }
      free(df->variables);    }

    /* Deallocate memory associated with the dimensions */
    if (df->dimensions != NULL) {
      for (i = 0; i < df->nd; ++i)
      {
        if(df->dimensions[i].name != NULL)
        {
          free(df->dimensions[i].name);
        }
      }
      free(df->dimensions);
    }

    /* Deallocate memory associated with the attributes */
    if (df->attributes != NULL) {
      for (i = 0; i < df->na; ++i)
      {
        free(df->attributes[i].name);
        if(df->attributes[i].value != NULL)/*UR-FIXED */
          free(df->attributes[i].value);
      }
      free(df->attributes);
    }

    // Remove hash memory
    if (df->ht_dist != NULL) {
      // Free all the data
      ht_process(df->ht_dist, free);
      ht_destroy(df->ht_dist);
    }
    if (df->ht_k != NULL) {
      // Free all the data
      ht_process(df->ht_k, free);
      ht_destroy(df->ht_k);
    }
    if (df->ht_master_1d != NULL) {
      // Free the hash-queue's within each table first
      ht_process(df->ht_master_1d, destroy_hq_within_ht);
      ht_destroy(df->ht_master_1d);
    }
    if (df->ht_master_2d != NULL) {
      // Free the hash-queue's within each table first
      ht_process(df->ht_master_2d, destroy_hq_within_ht);
      ht_destroy(df->ht_master_2d);
    }
    emstag(LMETRIC,"lib:datafile:df_free","Destroyed datafile %s ",df->name);
    free(df);  
  }
}


/**
  * Read all data records into memory buffer.
  * 
  * @param df pointer to datafile structure
  * @param v  pointer to df_variable structure
  *
  * @see dfReadRecords 
  */
void df_read_all_records(datafile_t *df, df_variable_t *v)
{
  df_read_records(df, v, 0, df->nrecords);
}


/**
  * Read a data record into memory buffer.
  * This function is effectively an interface to dfReadRecords.
  * 
  * @param df pointer to datafile structure
  * @param v pointer to datafile variable structure
  * @param rec record number (Use 0 if variable is indep of record dim).
  *
  * @see dfReadRecords 
  */
void df_read_record(datafile_t *df, df_variable_t *v, int rec)
{
  df_read_records(df, v, rec, 1);
}


/** Find the index of a variable in a data file.
  * 
  * @param df pointer to data file structure
  * @param name name of variable to find
  * @return variable index. -1 if unsuccessful.
  */
int df_get_index(datafile_t *df, const char *name)
{
  int i;

  if (name == NULL)
    return (-1);

  for (i = 0; i < df->nv; i++)
  {
    if (strcasecmp(name, df->variables[i].name) == 0 ||
        ((df->variables[i].longname != NULL) &&
         (strcasecmp(name, df->variables[i].longname) == 0)
        ))
      return (i);
  }
  return (-1);
}

/** Get the number of dimensions for the specified variable.
  *
  * @param df pointer to datafile structure.
  * @param v pointer to datafile variable.
  * @return Number of dimensions.
  * @see dfGetDimIds
  */
int df_get_num_dims(datafile_t *df, df_variable_t *v)
{
  return v->nd;
}


/** Get the dimension ids for the specified variable.
  *
  * @param df pointer to datafile structure.
  * @param v pointer to datafile variable.
  * @return pointer to array of dimension ids.
  * @see dfGetNumDimensions
  */
int *df_get_dim_ids(datafile_t *df, df_variable_t *v)
{
  return v->dimids;
}


/** Convienience function to find the number of variables.
  *
  * @param df pointer to data file structure.
  * @return Number of variables in the datafile structure.
  */
int df_get_num_variables(datafile_t *df)
{
  return df->nv;
}


/** Convienience function to find the variable structure for the
  * specified variable id.
  *
  * @param df pointer to data file structure.
  * @param varid variable index.
  * @return pointer to variable. NULL if out of bounds.
  */
df_variable_t *df_get_variable(datafile_t *df, int varid)
{
  if ((varid < 0) || (varid >= df->nv))
    return NULL;

  return &df->variables[varid];
}


/** Convienience function to find the variable structure by name.
  *
  * @param df pointer to data file structure.
  * @param name variable name.
  * @return pointer to variable. NULL if out of bounds.
  */
df_variable_t *df_get_variable_by_name(datafile_t *df, const char *name)
{
  return df_get_variable(df, df_get_index(df, name));
}


/** Convienience function to find the variable name by variable id
  *
  * @param df pointer to data file structure.
  * @param vid variable id
  * @return pointer to name
  */
const char *df_get_variable_name_by_id(datafile_t *df, int vid)
{
  return (df_get_variable(df, vid))->name;
}


/** Convienience function to find the dimension structure for the
  * specified dimension id.
  *
  * @param df pointer to data file structure.
  * @param dimid dimension index.
  * @return pointer to dimension. NULL if out of bounds.
  */
df_dimension_t *df_get_dimension(datafile_t *df, int dimid)
{
  if ((dimid < 0) || (dimid >= df->nd))
    return NULL;

  return &df->dimensions[dimid];
}



/** Find the global attribute by name, and return the attribute
  * structure.
  *
  * @param df pointer to datafile structure.
  * @param name name of the global attribute variable.
  * @return pointer to the attribute structure. NULL if unsuccessful.
  */
df_attribute_t *df_get_global_attribute(datafile_t *df, const char *name)
{
  int i;

  for (i = 0; i < df->na; ++i) {
    if (strcasecmp(df->attributes[i].name, name) == 0)
      return &df->attributes[i];
  }

  return NULL;
}


/** Find the global attribute by name, and return the attribute
  * structure.
  *
  * @param df pointer to datafile structure.
  * @param v pointer to the variable structure.
  * @param name name of the global attribute variable.
  * @return pointer to the attribute structure. NULL if unsuccessful.
  */
df_attribute_t *df_get_attribute(datafile_t *df, df_variable_t *v,
                                 const char *name)
{
  int i;

  for (i = 0; i < v->na; ++i) {
    if (strcasecmp(v->attributes[i].name, name) == 0)
      return &v->attributes[i];
  }

  return NULL;
}

/** Check if this is a RECOM file
 *
 */
int df_is_recom(datafile_t *df)
{
  if (df_get_global_attribute(df, "RECOM") != NULL)
    return(1);
  return(0);  
}


/** See if this is a UGRID file
 *
 */
int df_is_ugrid(datafile_t *df)
{
  df_attribute_t *at = df_get_global_attribute(df, "Conventions");
  if (at != NULL)
    return(strncmp(ATT_TEXT(at),"UGRID",5) == 0);
  return(0);  
}

/** See if this is a UGRID3 file
 *
 */
int df_is_ugrid3(datafile_t *df)
{
  df_attribute_t *at = df_get_global_attribute(df, "Conventions");
  if (at != NULL) {
    if (strncmp(ATT_TEXT(at),"UGRID",5) == 0) {
      at = df_get_global_attribute(df, "UGRID");
      return(strncmp(ATT_TEXT(at),"3D",2) == 0);
    }
  }
  return(0);  
}

int df_is_irule(datafile_t *df)
{
  if (strlen(df->i_rule))
    return(1);
  else
    return(0);
}


/** Check if positive down
 *
 */
int df_var_z_is_depth(df_variable_t *v)
{
  return(v->z_is_depth);
}

/** Add a text attribute to the variable. This is quite an
  * expensive routine as it requires the table to be reallocated.
  *
  * @param df pointer to the datafile structure.
  * @param v pointer to the variable structure.
  * @param name attribute name.
  * @param text attribute text.
  */
void df_add_text_attribute(datafile_t *df, df_variable_t *v,
                           const char *name, const char *text)
{
  int i;
  int na = v->na + 1;
  df_attribute_t *atts = NULL;

  atts = (df_attribute_t *)malloc(na * sizeof(df_attribute_t));
  memset(atts, 0, na * sizeof(df_attribute_t));

  /* Copy over the attribute information to the new atts array */
  for (i = 0; i < v->na; ++i) {
    atts[i].name = v->attributes[i].name;
    atts[i].n = v->attributes[i].n;
    atts[i].type = v->attributes[i].type;
    atts[i].value = v->attributes[i].value;
  }

  /* Now add the new information to the last element in the table. */
  atts[v->na].name = (char *)malloc((strlen(name) + 1) * sizeof(char));
  strcpy(atts[v->na].name, name);
  atts[v->na].n = strlen(text);
  atts[v->na].type = AT_TEXT;
  atts[v->na].value = (void *)malloc(sizeof(char) * (atts[v->na].n + 1));
  strcpy((char *)atts[v->na].value, text);


  /* Swizzle the old for the new, and free the old. */
  if (v->attributes)
    free(v->attributes);
  v->na = na;
  v->attributes = atts;

}


/**
  * Check that record values in a data file are
  * monotonic increasing
  * 
  * @param df pointer to data file structure
  * 
  * @see quit() If anything goes wrong.
  */
void df_check_records(datafile_t *df)
{
  int i;

  for (i = 1; i < df->nrecords; i++) {
    if (df->records[i] <= df->records[i - 1])
      quit("tsCheckRecords: Records out of order, %s, record=%g\n",
	   df->name, df->records[i]);
  }
}

/**
  * Print summary information aboout a data file.
  * 
  * @param df pointer to data file structure
  * @param fp output FILE pointer
  * @param level verbosity level
  * 
  * @see quit() If anything goes wrong.
  */
void df_print_info(datafile_t *df, FILE * fp, int level)
{
  int last = df->nrecords - 1;
  int i, j, n;

  fprintf(fp, "Data file: %s\n", df->name);
  if ((df->nrecords > 0) && (df->records != NULL)) {
    fprintf(fp, "Units %s\n", df->rec_units);
    fprintf(fp, "Starts=%.12g - ", df->records[0]);
    fprintf(fp, "Ends=%.12g\n", df->records[last]);

    if (level > 2) {
      for (i = 0; i < df->na; ++i) {
        df_attribute_t *a = &df->attributes[i];
        fprintf(fp, "Global attribute %d: %s = ", i, a->name);
        if (a->type != AT_TEXT) {
          for (j = 0; j < (int)a->n; ++j) {
            switch (a->type) {
            case AT_BYTE:
              fprintf(fp, "%c", ATT_BYTE(a, j));
              break;

            case AT_FLOAT:
              fprintf(fp, "%g", ATT_FLOAT(a, j));
              break;

            case AT_DOUBLE:
              fprintf(fp, "%g", ATT_DOUBLE(a, j));
              break;

            case AT_SHORT:
              fprintf(fp, "%d", (int)ATT_SHORT(a, j));
              break;

            case AT_INT:
              fprintf(fp, "%d", ATT_INT(a, j));
              break;

            default:
              break;
            }
          }
        } else
          fprintf(fp, "%s",(char *)a->value);
        fprintf(fp, "\n");
      }
    }
  } else
    fprintf(fp, "No named record available.\n");

  if (level <= 0) {
    if (df->type == DFT_ASCII)
      fprintf(fp, "Column\tVariable\n");
    else
      fprintf(fp, "ID\tVariable\n");

    for (n = 0; n < df->nv; n++) {
      fprintf(fp, "%d\t%s\n", n + (df->type == DFT_ASCII),
              df->variables[n].name);
    }
  } else {
    for (n = 0; n < df->nv; n++) {
      df_variable_t *v = &df->variables[n];

      if (df->type == DFT_ASCII)
        fprintf(fp, "Column:       %d\n", n + 1);
      else
        fprintf(fp, "ID:           %d\n", n);

      fprintf(fp, "Variable:     %s\n", v->name);
      fprintf(fp, "Long name:    %s\n",
              (v->longname) ? v->longname : v->name);
      fprintf(fp, "Units:        %s\n",
              (v->units) ? v->units : "no units");
      fprintf(fp, "Missing:      %g\n", v->missing);
      fprintf(fp, "Fill value:   %g\n", v->fillvalue);
      if (level > 2) {
        for (i = 0; i < v->na; ++i) {
          df_attribute_t *a = &v->attributes[i];
          fprintf(fp, "Attrib %d:     %s = ", i, a->name);
          if (a->type != AT_TEXT) {
            for (j = 0; j < (int)a->n; ++j) {
              switch (a->type) {
              case AT_BYTE:
                fprintf(fp, "%c", ATT_BYTE(a, j));
                break;

              case AT_FLOAT:
                fprintf(fp, "%g", ATT_FLOAT(a, j));
                break;

              case AT_DOUBLE:
                fprintf(fp, "%g", ATT_DOUBLE(a, j));
                break;

              case AT_SHORT:
                fprintf(fp, "%d", (int)ATT_SHORT(a, j));
                break;

              case AT_INT:
                fprintf(fp, "%d", ATT_INT(a, j));
                break;

              default:
                break;
              }
            }
          } else
            fprintf(fp, "%s", (char *)a->value);
          fprintf(fp, "\n");
        }
      }
      if (level > 1) {
        fprintf(fp, "Dimensions:  ");
        if (v->dim_as_record) {
          fprintf(fp, " %s",
                  df->dimensions[df->variables[df->ri].dimids[0]].name);

        }
        for (i = 0; i < v->nd; ++i)
          fprintf(fp, " %s", df->dimensions[v->dimids[i]].name);
        fprintf(fp, "\n");

        fprintf(fp, "Coordinates: ");
        for (i = 0; i < v->nc; ++i) {
          fprintf(fp, " %s", df->variables[v->coordids[i]].name);
          if (level > 2) {
            switch (df->variables[v->coordids[i]].type) {
            case VT_X:
              fprintf(fp, "(X)");
              break;

            case VT_Y:
              fprintf(fp, "(Y)");
              break;

            case VT_Z:
              fprintf(fp, "(Z)");
              break;

            case VT_TIME:
              fprintf(fp, "(TIME)");
              break;

            case VT_LONGITUDE:
              fprintf(fp, "(LONGITUDE)");
              break;

            case VT_LATITUDE:
              fprintf(fp, "(LATITUDE)");
              break;

            default:
              break;
            }
          }
        }
        fprintf(fp, "\n");
      }
      fprintf(fp, "\n\n");
    }
  }

  fprintf(fp, "\n");
}


/** Get the data value for the specified record and indicies.
  *
  * @param df pointer to datafile structure.
  * @param v pointer to variable structure.
  * @param record record index (not relative record index).
  * @param is array of indice values.
  * @return Value at specified record and index.
  *
  * @see dfGetNumDimensions
  */
double df_get_data_value(datafile_t *df, df_variable_t *v, int record,
                         int *is)
{
  int ri = record - v->start_record;
  double val = 0.0;

  if (ri < 0)
    ri = 0;
  if (ri >= v->nrecords)
    ri = v->nrecords - 1;

  switch (v->nd) {
  case 0:
    val = VAR_0D(v)[ri];
    break;

  case 1:
    val = VAR_1D(v)[ri][is[0]];
    break;

  case 2:
    val = VAR_2D(v)[ri][is[0]][is[1]];
    break;

  case 3:
    val = VAR_3D(v)[ri][is[0]][is[1]][is[2]];
    break;

  case 4:
    val = VAR_4D(v)[ri][is[0]][is[1]][is[2]][is[3]];
    break;

  default:
    quit("dfgetDataValue: Bad number of dimensions\n");
  }

  return val;
}


/* PRIVATE or PROTECTED functions */


/*
Routine to read data file from a netCDF file. The netcdf read
overrides any record name already specified.

Arguments:
    fid     -    netCDF file id
    df      -    pointer to datafile structure (assumed
	         not previously initialised)

This routine calls quit() if anything goes wrong
*/
void netcdf_read(int fid, datafile_t *df, int type)
{
  int i, j;
  char line[MAXLINELEN];
  int unlimited = -1;
  int inferred = 0;

  /* Query the netCDF file for the total number of dimensions, allocate
     the array, and populate. Special care must be taken define which
     dimension is the special record (or unlimited dimension) */
  df->type = type;
  df->ncid = fid;
  df->ri = -1;
  nc_inq_ndims(df->ncid, &df->nd);
  if (df->nd > 0) {
    df->dimensions =
      (df_dimension_t *)malloc(df->nd * sizeof(df_dimension_t));
    memset(df->dimensions, 0, df->nd * sizeof(df_dimension_t));
    for (i = 0; i < df->nd; ++i) {
      df->dimensions[i].dimid = i;
      nc_inq_dim(df->ncid, i, line, &df->dimensions[i].size);
      df->dimensions[i].name =
        (char *)malloc(sizeof(char) * (strlen(line) + 1));
      strcpy(df->dimensions[i].name, line);
    }

    /* Is there an unlimited dimension ? */
    nc_inq_unlimdim(df->ncid, &unlimited);
  } else {
    warn("No dimensions were specified in this netCDF file.\n");
    return;
  }

  /* Query the netCDF file for the total number of variables, and allocate 
     the array */
  nc_inq_nvars(df->ncid, &df->nv);
  if (df->nv > 0) {
    df->variables =
      (df_variable_t *)malloc(df->nv * sizeof(df_variable_t));
    memset(df->variables, 0, df->nv * sizeof(df_variable_t));
    for (i = 0; i < df->nv; ++i) {
      df_variable_t *v = &df->variables[i];

      /* df_variable_t type and netCDF id - default to data */
      df->variables[i].type = VT_DATA;
      v->varid = i;
      v->dim_as_record = 0;

      /* df_variable_t name */
      nc_inq_varname(df->ncid, v->varid, line);
      v->name = (char *)malloc(sizeof(char *) * (strlen(line) + 1));
      strcpy(v->name, line);

      /* Extract the dimension information for each variable */
      nc_inq_varndims(df->ncid, v->varid, &v->nd);

      if (v->nd == 0) {
	warn("netcdf_read: 0 coords found for variable %s\n", v->name);
	continue;
      }
      v->dimids = (int *)malloc(v->nd * sizeof(int));
      nc_inq_vardimid(df->ncid, v->varid, v->dimids);

      /* MH. Allow special treatment for GHRSST netCDF files. These files */
      /* have a time dimension of 1 rather than UNLIMITED, no coordinate  */
      /* types for lat and lon, and only 2 coordinates for the data       */
      /* (analysed_sst). If v->type |= VT_COORD, then these issues are    */
      /* forced.                                                          */
      if (v->nd == 3) {
	j = 0;
	while (sst_names[j] != NULL) {
	  if (strcmp(v->name, sst_names[j]) == 0) {
	    v->type |= VT_COORD;
	    break;
	  }
	  j++;
	}
      }
      /*
      if (strcmp(v->name, "analysed_sst") == 0 && v->nd == 3) 
	v->type |= VT_COORD;
      */
      /* Loop through all of the attributes. Take special care to with the 
         long_name, units, missing, _FillValue, and * coord type.
         Coordinates are checked later. */
      nc_inq_varnatts(df->ncid, v->varid, &v->na);
      if (v->na > 0) {
        v->attributes =
          (df_attribute_t *)malloc(v->na * sizeof(df_attribute_t));
        memset(v->attributes, 0, v->na * sizeof(df_attribute_t));
        v->scale_factor = 1.0;
        v->add_offset = 0.0;
        for (j = 0; j < v->na; ++j) {
          df_attribute_t *a = &v->attributes[j];

          netcdf_read_attrib(df, v->varid, j, a);

          /* Check for known attribute types */
          if (strcasecmp(a->name, "long_name") == 0)
            v->longname = ATT_TEXT(a);

          else if (strcasecmp(a->name, "units") == 0)
            v->units = ATT_TEXT(a);

          else if (strcasecmp(a->name, "missing_value") == 0) {
	    if (CHK_TYPE(a, AT_BYTE))
	      v->missing = ATT_BYTE(a, 0);
	    else if (CHK_TYPE(a, AT_FLOAT))
	      v->missing = ATT_FLOAT(a, 0);
	    else if (CHK_TYPE(a, AT_SHORT))
	      v->missing = ATT_SHORT(a, 0);
	    else
	      v->missing = ATT_DOUBLE(a, 0);
	    v->type |= VT_MV;
	  }
          else if (strcasecmp(a->name, "_FillValue") == 0) {
	    if (CHK_TYPE(a, AT_BYTE))
	      v->fillvalue = ATT_BYTE(a, 0);
	    else if (CHK_TYPE(a, AT_SHORT))
	      v->fillvalue = ATT_SHORT(a, 0);
	    else if (CHK_TYPE(a, AT_FLOAT))
	      v->fillvalue = ATT_FLOAT(a, 0);
	    else
	      v->fillvalue = ATT_DOUBLE(a, 0);
	    v->type |= VT_FV;
	  }
          else if (strcasecmp(a->name, "scale_factor") == 0) {
	    if (CHK_TYPE(a, AT_DOUBLE))
	      v->scale_factor = ATT_DOUBLE(a, 0);
	    else
	      v->scale_factor = ATT_FLOAT(a, 0);
	  }
          else if (strcasecmp(a->name, "add_offset") == 0) {
	    if (CHK_TYPE(a, AT_DOUBLE))
	      v->add_offset = ATT_DOUBLE(a, 0);
	    else
	      v->add_offset = ATT_FLOAT(a, 0);
	  }
          else if (strcasecmp(a->name, "positive") == 0)
            v->z_is_depth = strcasecmp(ATT_TEXT(a), "down") == 0;

          else if (strcasecmp(a->name, "coordinate_type") == 0) {
            /* Decode the coordinate types */
            decode_coord_type(ATT_TEXT(a), &v->type, &v->coord_domain);
          }

          else if (strcasecmp(a->name, "land_flag") == 0) {
	    if (CHK_TYPE(a, AT_BYTE))
	      v->lflag = ATT_BYTE(a, 0);
	    else
	      v->lflag = ATT_DOUBLE(a, 0);
	  }

          else if (strcasecmp(a->name, "cloud_flag") == 0) {
	    if (CHK_TYPE(a, AT_BYTE))
	      v->cflag = ATT_BYTE(a, 0);
	    else
	      v->cflag = ATT_DOUBLE(a, 0);
	  }
        }

        v->missing = v->missing * v->scale_factor + v->add_offset;
        v->fillvalue = v->fillvalue * v->scale_factor + v->add_offset;
        v->lflag = v->lflag * v->scale_factor + v->add_offset;
        v->cflag = v->cflag * v->scale_factor + v->add_offset;
      }

      /* MH. If a variable has the same name as a coordinate, and */
      /* its units are degrees or metres, then if its type is     */
      /* VT_DATA, reset to the appropriate coordinate.            */ 
      if (v->type & VT_DATA) {
	for (j = 0; j < df->nd; ++j) {
	  if (strcmp(v->name, df->dimensions[j].name) == 0) {
	    if (has_name(v->units, lon_names)) v->type = VT_LONGITUDE;
	    if (has_name(v->units, lat_names)) v->type = VT_LATITUDE;
	    if (has_name(v->units, v_names)) v->type = VT_Z;
	    inferred = 1;
	  }
	}
      }
    }

    /* Read the global attributes */
    nc_inq_natts(df->ncid, &df->na);
    if (df->na > 0) {
      df->attributes =
        (df_attribute_t *)malloc(df->na * sizeof(df_attribute_t));
      memset(df->attributes, 0, df->na * sizeof(df_attribute_t));
      for (j = 0; j < df->na; ++j)
        netcdf_read_attrib(df, NC_GLOBAL, j, &df->attributes[j]);
    }

    /* Search for a variable that has the same name as the unlimited
       dimension name. This will become the coordinate dimension. */
    for (i = 0; i < df->nv; ++i) {
      df_variable_t *v = &df->variables[i];

      if ((v->nd == 1) && (v->dimids[0] == unlimited)) {
        if (strcasecmp(v->name, df->dimensions[v->dimids[0]].name) == 0) {
          df_set_record(df, i);
          break;
        }
      }
      /* MH. Find the time dimension for GHRSST netCDF and set to coordinate dimension. */
      if (v->type & VT_COORD) {
	if ((v->nd == 1) && strcmp(v->name, "time") == 0 && df->dimensions[v->dimids[0]].size == 1) {
	  if (strcasecmp(v->name, df->dimensions[v->dimids[0]].name) == 0) {
	    df_set_record(df, i);
	    break;
	  } else
	    v->type &= ~VT_COORD;
	}
      }
      /* END MH */
    }

    /* Sweep back trough the variables and decode the associated
       coordinates information if the attribute is present. This * has to
       be left until the end to ensure all variables were * read. */
    for (i = 0; i < df->nv; ++i) {
      df_variable_t *v = &df->variables[i];
      char *text = NULL;

      for (j = 0; j < v->na; ++j) {
        df_attribute_t *a = &v->attributes[j];
        if (strcasecmp(a->name, "coordinates") == 0) {
          text = ATT_TEXT(a);
        }
      }

      decode_coords(df, v, text);

      /* MH : set the coordinate types (for GHRSST data) if not already set     */
      if (strcmp(v->name, "time") == 0 && v->nc == 0 && v->type & VT_DATA)
	v->type = VT_TIME;
      if (strcmp(v->name, "lon") == 0 && v->nc == 0 && v->type & VT_DATA)
	v->type = VT_LONGITUDE;
      if (strcmp(v->name, "lat") == 0 && v->nc == 0 && v->type & VT_DATA)
	v->type = VT_LATITUDE;

      /* MH : If the variable is botz and no coordinates were found, assume the */
      /* coordinates are x_centre, y_centre and try again.                      */
      if (strcmp(v->name, "botz") == 0 && v->nc == 0) {
	char buf[MAXSTRLEN];
	strcpy(buf, "x_centre, y_centre");
	df_add_text_attribute(df, v, "coordinates", buf);
	for (j = 0; j < v->na; ++j) {
	  df_attribute_t *a = &v->attributes[j];
	  if (strcasecmp(a->name, "coordinates") == 0) {
	    text = ATT_TEXT(a);
	  }
	}
	decode_coords(df, v, text);
	v->type |= VT_BATHY;
      }
      /* MH : If the variable is height and no coordinates were found, assume   */
      /* the coordinates are lon, lat and try again.    */
      if (strcmp(v->name, "height") == 0 && v->nc == 0 && df->records == NULL) {
	int ii;
	char buf[MAXSTRLEN];
	for (ii = 0; ii < df->nv; ++ii) {
	  df_variable_t *vc = &df->variables[ii];
	  if (strcmp(vc->name, "lon") == 0) {
	    strcpy(vc->units, "degrees_E");
	    strcpy(buf, "geographic");
	    df_add_text_attribute(df, vc, "projection", buf);
	  }
	  if (strcmp(vc->name, "lat") == 0) {
	    strcpy(vc->units, "degrees_N");
	    strcpy(buf, "geographic");
	    df_add_text_attribute(df, vc, "projection", buf);
	  }
	}
	df->nrecords = 1;
	v->type |= VT_BATHY;
      }

      /* MH : Check if the coordinates can be constructed from existing variables. */
      /* Only for 2d and 3d spatial data in lat/long coordinates. */
      if (inferred) v->type |= VT_INFERRED;
      if (v->type & VT_DATA && v->nc == 0 && v->nd > 1) {
	int n;
	v->coordids = (int *)malloc(v->nd * sizeof(int));
	memset(v->coordids, 0, sizeof(v->nd * sizeof(int)));
	for (j = 0; j < v->nd; ++j) {
	  for (n = 0; n < df->nv; ++n) {
	    df_variable_t *nv = &df->variables[n];
	    if(strcmp(df->dimensions[v->dimids[j]].name, nv->name) == 0) {
	      v->coordids[v->nd-(j+1)] = n;
	      v->nc++;
	    }
	  }
	}

	if (v->nc == 0) {
	  free(v->coordids);
	  v->coordids = NULL;
	} else {
	  for (j = 0; j < v->nc; ++j) {
	    df_variable_t *cv = &df->variables[v->coordids[j]];
	    if (cv->type & VT_DATA) {
	      if (has_name(cv->units, lon_names)) cv->type = VT_LONGITUDE;
	      if (has_name(cv->units, lat_names)) cv->type = VT_LATITUDE;
	      if (has_name(cv->units, v_names)) cv->type = VT_Z;
	    }
	    warn("netcdf_read: Coordinates for %s not found; using coord%d = %s\n",
		 v->name, j, cv->name);
	  }
	  v->type |= VT_INFERRED;
	}
      }

      /* Set the coordinate types for ghrsst netCDF. */
      if (v->type & VT_COORD) {
	if (v->type & VT_DATA && v->nc == 0 && v->nd >= 1) {
	  int n;
	  v->coordids = (int *)malloc(v->nd * sizeof(int));
	  memset(v->coordids, 0, sizeof(v->nd * sizeof(int)));
	  for (j = 0; j < v->nd; ++j) {
	    for (n = 0; n < df->nv; ++n) {
	      df_variable_t *nv = &df->variables[n];
	      if(strcmp(df->dimensions[v->dimids[j]].name, nv->name) == 0) {
		v->coordids[v->nd-(j+1)] = n;
		v->nc++;
	      }
	    }
	  }
	  if (v->nc == 0) {
	    free(v->coordids);
	    v->coordids = NULL;
	  } else {
	    for (j = 0; j < v->nc; ++j) {
	      df_variable_t *cv = &df->variables[v->coordids[j]];
	      if (cv->type & VT_DATA) {
		if (has_name(cv->units, lon_names)) cv->type = VT_LONGITUDE;
		if (has_name(cv->units, lat_names)) cv->type = VT_LATITUDE;
		if (has_name(cv->units, v_names)) cv->type = VT_Z;
	      }
	      warn("netcdf_read: Coordinates for %s not found; using coord%d = %s\n",
		   v->name, j, cv->name);
	    }
	  }
	}
      }
      /* end MH */
    }
  } else {
    warn("No variables were specified in this netCDF file.\n");
    return;
  }

}


void netcdf_read_records(datafile_t *df, int varid) {
    df_variable_t *rv = &df->variables[varid];
    int recdimid = rv->dimids[0];

    if (rv->data != NULL)
      free(rv->data);
    rv->data = (double *)malloc(sizeof(double) * df->dimensions[recdimid].size);
    nc_get_var_double(df->ncid, varid, rv->data);
}


/* Read a data record into the specified variable from a netCDF
 * file. It is assumed that the memory has already been allocated.
 */
void netcdf_read_data(datafile_t *df, df_variable_t *v, int r, int roffset)
{
  size_t start[5];
  size_t count[5];
  int i = 0, j, k, l;
  int index = r - v->start_record;
  int fid = df->ncid;

  if (!(df->rec_modulus) && ((r < 0) || (r >= df->nrecords)))
    quit
      ("netcdf_read_data: Attempt to read an invalid record for variable '%s'.\n",
       v->name);

  if (v->dim_as_record) {
    if (df->rec_modulus)
      start[i] = (r % df->nrecords) - roffset;
    else
      start[i] = r - roffset;
    count[i++] = 1;
  }

  switch (v->nd) {
  case 0:
    nc_get_vara_double(fid, v->varid, start, count, &v->data[index]);
    v->data[index] = v->data[index] * v->scale_factor + v->add_offset;
    break;

  case 1:
    start[i] = 0;
    count[i++] = df->dimensions[v->dimids[0]].size;
    nc_get_vara_double(fid, v->varid, start, count, VAR_1D(v)[index]);
    for (i = 0; i < (int)df->dimensions[v->dimids[0]].size; ++i)
      VAR_1D(v)[index][i]
        = VAR_1D(v)[index][i] * v->scale_factor + v->add_offset;
    break;

  case 2:
    start[i] = 0;
    count[i++] = df->dimensions[v->dimids[0]].size;
    start[i] = 0;
    count[i++] = df->dimensions[v->dimids[1]].size;
    nc_get_vara_double(fid, v->varid, start, count, VAR_2D(v)[index][0]);
    for (i = 0; i < (int)df->dimensions[v->dimids[0]].size; ++i)
      for (j = 0; j < (int)df->dimensions[v->dimids[1]].size; ++j)
        VAR_2D(v)[index][i][j]
          = VAR_2D(v)[index][i][j] * v->scale_factor + v->add_offset;
    break;

  case 3:
    start[i] = 0;
    count[i++] = df->dimensions[v->dimids[0]].size;
    start[i] = 0;
    count[i++] = df->dimensions[v->dimids[1]].size;
    start[i] = 0;
    count[i++] = df->dimensions[v->dimids[2]].size;
    nc_get_vara_double(fid, v->varid, start, count,
                       VAR_3D(v)[index][0][0]);
    for (i = 0; i < (int)df->dimensions[v->dimids[0]].size; ++i)
      for (j = 0; j < (int)df->dimensions[v->dimids[1]].size; ++j)
        for (k = 0; k < (int)df->dimensions[v->dimids[2]].size; ++k)
          VAR_3D(v)[index][i][j][k]
            = VAR_3D(v)[index][i][j][k] * v->scale_factor + v->add_offset;
    break;

  case 4:
    start[i] = 0;
    count[i++] = df->dimensions[v->dimids[0]].size;
    start[i] = 0;
    count[i++] = df->dimensions[v->dimids[1]].size;
    start[i] = 0;
    count[i++] = df->dimensions[v->dimids[2]].size;
    start[i] = 0;
    count[i++] = df->dimensions[v->dimids[3]].size;
    nc_get_vara_double(fid, v->varid, start, count,
                       VAR_4D(v)[index][0][0][0]);
    for (i = 0; i < (int)df->dimensions[v->dimids[0]].size; ++i)
      for (j = 0; j < (int)df->dimensions[v->dimids[1]].size; ++j)
        for (k = 0; k < (int)df->dimensions[v->dimids[2]].size; ++k)
          for (l = 0; l < (int)df->dimensions[v->dimids[3]].size; ++l)
            VAR_4D(v)[index][i][j][k][l]
              = VAR_4D(v)[index][i][j][k][l] * v->scale_factor +
              v->add_offset;
    break;

  default:
    quit("read_data: Bad number of dimensions\n");
  }
}



/*
 * Free up the memory associated with a netcdf data file.
 */
void netcdf_free(datafile_t *df)
{
  nc_close(df->ncid);

}

/* Multi-netCDF reader. It is assumed that the files all contain
 * the same data but simply for different records. The dataset may
 * overlap in time
 */
void multi_netcdf_read(FILE *fp, datafile_t *df)
{
  int fid = 0;
  char buf[MAXLINELEN];
  datafile_t *df0 = malloc(sizeof(datafile_t));
  df_multi_t *fd = (df_multi_t *)malloc(sizeof(df_multi_t));

  memset(fd,  0, sizeof(df_multi_t));
  memset(df0, 0, sizeof(datafile_t));

  df->private_data = fd;

  prm_set_errfn(quit);

  prm_read_double(fp, "multi-netcdf-version", &fd->version);
  prm_read_int(fp, "nfiles", &fd->nfiles);

  if (fd->nfiles > 0) {
     int i=0;
     fd->files = (df_multi_file_t *)malloc(fd->nfiles
               * sizeof(df_multi_file_t));
     memset(fd->files, 0, fd->nfiles * sizeof(df_multi_file_t));
     
     for (i=0; i<fd->nfiles; ++i) {
       sprintf(buf, "file%d.filename", i);
       prm_read_char(fp, buf, fd->files[i].filename);
     }

     /* Check and open df0 */
     if (nc_open(fd->files[0].filename, NC_NOWRITE, &fid) == NC_NOERR)
       netcdf_read(fid, df0, DFT_NETCDF);
     else
       quit("Error opening file '%s'\n", fd->files[0].filename);

     /*
      * Loop over every file and check for conformity
      */
     for (i=1; i<fd->nfiles; ++i) {
       int d;
       datafile_t *df1 = malloc(sizeof(datafile_t));
       memset(df1, 0, sizeof(datafile_t));
       if (nc_open(fd->files[i].filename, NC_NOWRITE, &fid) == NC_NOERR) {
	 netcdf_read(fid, df1, DFT_NETCDF);
	 /* Check number of dimensions */
	 if (df0->nd != df1->nd)
	   quit("Error: multinecdf file '%s' number of dimensions (%d) does not match between %s and %s\n",
		df->name, df0->nd, fd->files[0].filename, fd->files[i].filename);

	 /* Check dimension order */
	 for (d=0; d<df0->nd; d++) {
	   if (strcmp(df0->dimensions[d].name, df1->dimensions[d].name))
	     quit("Error: multinecdf file '%s' dimension order (%d) mismatch (%s:%s) vs (%s:%s)\n",
		  df->name, d, 
		  df0->dimensions[d].name, fd->files[0].filename, 
		  df1->dimensions[d].name, fd->files[i].filename);
	 }

	 /* Check number of variables */
	 if (df0->nv != df1->nv)
	   quit("Error: multinecdf file '%s' number of variables does not match between %s(%d) and %s(%d)\n",
		df->name, fd->files[0].filename, df0->nv, fd->files[i].filename, df1->nv);

	 /* Check variable order */
	 for (d=0; d<df0->nv; d++) {
	   if (strcmp(df0->variables[d].name, df1->variables[d].name))
	     quit("Error: multinecdf file '%s' variable order (%d) mismatch (%s:%s) vs (%s:%s)\n",
		  df->name, d, 
		  df0->variables[d].name, fd->files[0].filename, 
		  df1->variables[d].name, fd->files[i].filename);
	 }
	 /* MH : Read the saved attributes for this file */
	 multi_netcdf_read_atts(df1, &fd->files[i]);
	 /* Clean up */
	 df_free(df1);
	 df1 = NULL;

       } else
	 quit("Error opening the netcdf file '%d:%s'\n", 
	      i, fd->files[i].filename);
     }

     /* Cleanup df0 */
     df_free(df0);
     
     /*
      * Open the first file and populate the datafile structure
      */
     nc_open(fd->files[0].filename, NC_NOWRITE, &fid);
     netcdf_read(fid, df, DFT_MULTI_NETCDF);
     /* MH : Read the saved attributes for this file */
     multi_netcdf_read_atts(df, &fd->files[0]);
     nc_close(df->ncid);
     df->ncid = -1;
  }

  fclose(fp);
}

/* MH. Routine to read and save selected attributes for multi_netcdf files */
void multi_netcdf_read_atts(datafile_t *df, df_multi_file_t *fd) {
  int n, nv = df->nv;
  fd->offset = d_alloc_1d(nv);
  fd->scale = d_alloc_1d(nv);
  fd->fill = d_alloc_1d(nv);

  for (n = 0; n < nv; n++) {
    df_variable_t *v = &df->variables[n];
    df_attribute_t *a = df_get_attribute(df, v, "add_offset");
    fd->scale[n] = 1.0;
    fd->offset[n] = 0.0;
    if (a != NULL) {
      if (CHK_TYPE(a, AT_DOUBLE))
	v->add_offset = ATT_DOUBLE(a, 0);
      else
	v->add_offset = ATT_FLOAT(a, 0);
      fd->offset[n] = (double)v->add_offset;
      a = df_get_attribute(df, v, "scale_factor");
      if (a != NULL) {
	if (CHK_TYPE(a, AT_BYTE))
	  v->scale_factor = ATT_BYTE(a, 0);
	else if (CHK_TYPE(a, AT_SHORT))
	  v->scale_factor = ATT_SHORT(a, 0);
	else if (CHK_TYPE(a, AT_FLOAT))
	  v->scale_factor = ATT_FLOAT(a, 0);
	else
	  v->scale_factor = ATT_DOUBLE(a, 0);
	if (v->scale_factor == 0.0) v->scale_factor = 1.0;
	fd->scale[n] = (double)v->scale_factor;
	a = df_get_attribute(df, v, "_FillValue");
	if (a != NULL) {
	  if (CHK_TYPE(a, AT_BYTE))
	    v->fillvalue = ATT_BYTE(a, 0);
	  else if (CHK_TYPE(a, AT_SHORT))
	    v->fillvalue = ATT_SHORT(a, 0);
	  else if (CHK_TYPE(a, AT_FLOAT))
	    v->fillvalue = ATT_FLOAT(a, 0);
	  else
	    v->fillvalue = ATT_DOUBLE(a, 0);
	  v->fillvalue = v->fillvalue * v->scale_factor + v->add_offset;
	  fd->fill[n] = (double)v->fillvalue;
	}
      }
    }
  }
}

void multi_netcdf_read_records(datafile_t *df, int varid) {
  df_multi_t *fd = (df_multi_t *)df->private_data;
  df_variable_t *rv = &df->variables[varid];
  int recdimid = rv->dimids[0];
  int i, j, k;
  int nr = 0;

  if (rv->data != NULL)
    free(rv->data);

  /* Read the records from the inidividual files. */
  for (i=0; i< fd->nfiles; ++i) {
    df_multi_file_t *f = &fd->files[i];
    if (f->records != NULL)
      free(f->records);

    /*
     * Count up the upper bound of records
     */
    if (nc_open(f->filename, NC_NOWRITE, &df->ncid) == NC_NOERR) {
       nc_inq_dimlen(df->ncid, df->dimensions[recdimid].dimid, &f->nrecords);
       f->records = d_alloc_1d(f->nrecords);
       nc_get_var_double(df->ncid, varid, f->records);
       nr += f->nrecords;
       nc_close(df->ncid);
       df->ncid = -1;
    }
    else
      quit("Error opening file '%s'\n", f->filename);
  }

  /*
   * Here we discard any overlapping regions
   */
  rv->data = d_alloc_1d(nr); // alloc upper bound
  k = 0;
  for (i=0; i< fd->nfiles; ++i) {
    df_multi_file_t *f = &fd->files[i];
    /* Figure out overlapping offset */
    j = 0;
    if (i>0) {
      df_multi_file_t *fp = &fd->files[i-1];
      double tLast = fp->records[fp->nrecords-1];
      /* Assume monotonically increasing */
      while (tLast >= f->records[j])
	j++;
    }
    f->rec_offset = j;
    /* Now go ahead and fill in */
    for (; j<f->nrecords; ++j, ++k) {
      rv->data[k] = f->records[j];
      if ((k > 0) && (rv->data[k] < rv->data[k-1]))
	quit("multi-netcdf record do not monotonically increase.\n"	\
	     "File %d out of order or corrupted? t1=%f, t2=%f in '%s'\n",
	     i, rv->data[k-1], rv->data[k], fd->files[i].filename);
    }
  }
  df->dimensions[recdimid].size = k;
}

/*
 * Handles overlapping region by utilising the record offset field of
 * each file
 */
void multi_netcdf_read_data(datafile_t *df, df_variable_t *v, int rec) {

  df_multi_t *fd = (df_multi_t *)df->private_data;
  if (fd != NULL) {
    int i, j, nr = 0, jv;
    
    /* Start counter where we left off last */
    for (i=0; i<fd->nfiles; ++i) {
      df_multi_file_t *f = &fd->files[i];
      int roff = f->rec_offset;
      
      if ((rec >= nr) && (rec < (nr+f->nrecords - roff))) {

	/* MH : Update selected attributes if saved */
	for (j = 0; j < df->nv; j++) {
	  if (f->scale[j] != 1.0 && f->offset[j] != 0.0) {
	    df_variable_t *rv = &df->variables[j];
	    if (strcmp(rv->name, v->name) == 0) {
	      v->add_offset = f->offset[j];
	      v->scale_factor = f->scale[j];
	      v->fillvalue = f->fill[j] * v->scale_factor + v->add_offset;
	      break;
	    }
	  }
	}

         if (nc_open(f->filename, NC_NOWRITE, &df->ncid) == NC_NOERR) {
	    netcdf_read_data(df, v, rec, nr-roff);
            nc_close(df->ncid);
            df->ncid = -1;
            return;
         } else
            quit("Error opening file '%s'\n", f->filename);
      }

      nr += f->nrecords - roff;
    }
  }
}

void multi_netcdf_free(datafile_t *df)
{
  df_multi_t *fd = (df_multi_t *)df->private_data;

  if (fd != NULL) {
    int i = 0;
    for (i=0; i<fd->nfiles; ++i) {
      if (fd->files[i].records != NULL)
         free(fd->files[i].records);
    }

    free(fd->files);
    free(fd);
    df->private_data = NULL;
  }
}


/* MH mempack */
/*
Routine to read data file from a transferred memory packet.

Arguments:
    df      -    pointer to datafile structure (assumed
	         not previously initialised)

This routine calls quit() if anything goes wrong
*/


void mempack_read(char *name, datafile_t *df, int type)
{
  int id;
  int i, j, m;
  char *vars[MAXLINELEN], *units[MAXLINELEN], *xyzunits[MAXLINELEN];
  char coords[MAXLINELEN];
  df_mempack_t *mp;

  /* Allocate memory for the memory packet */
  mp = (df_mempack_t *)malloc(sizeof(df_mempack_t));
  memset(mp, 0, sizeof(df_mempack_t));
  df->private_data = mp;

  /* Use the name to get the memory address */
  /* Farhan to do */
  strcpy(mp->id, name);
  mempack_read_file(name, df);
  mp->next = mp->time + mp->dt;
  /*memcpy(mp, id, sizeof(df_mempack_t));*/

  /* Memory packets always have 2 dimensions (npoints, nz). Allocate
     the array, and populate. */
  df->type = type;
  df->ncid = -1;
  df->ri = -1;
  df->nd = 3;
  df->dimensions =
    (df_dimension_t *)malloc(df->nd * sizeof(df_dimension_t));
  memset(df->dimensions, 0, df->nd * sizeof(df_dimension_t));
  /* Time */
  df->dimensions[0].dimid = 0;
  df->dimensions[0].size = 2;
  strcpy(coords, "time");
  df->dimensions[0].name = (char *)malloc(sizeof(char *)*MAXLINELEN);
  strcpy(df->dimensions[0].name, coords);
  /* Number of points */
  df->dimensions[1].dimid = 1;
  strcpy(coords, "npoints");
  df->dimensions[1].size = mp->npoints;
  df->dimensions[1].name = (char *)malloc(sizeof(char *)*MAXLINELEN);
  strcpy(df->dimensions[1].name, coords);
  /* Number of levels */
  df->dimensions[2].dimid = 2;
  strcpy(coords, "nz");
  df->dimensions[2].size = mp->nz;
  df->dimensions[2].name = (char *)malloc(sizeof(char *)*MAXLINELEN);
  strcpy(df->dimensions[2].name, coords);
  /* Query the structure for the total number of variables, and allocate 
     the array */
  i = parseline(mp->units, units, MAXLINELEN);
  i = parseline(mp->xyzunits, xyzunits, MAXLINELEN);
  df->nv = MPK_SVARS; /* Time, x, y and z */
  df->nv += parseline(mp->vars, vars, MAXLINELEN);
  if (df->nv > 0) {
    df_variable_t *v;
    df->variables =
      (df_variable_t *)malloc(df->nv * sizeof(df_variable_t));
    memset(df->variables, 0, df->nv * sizeof(df_variable_t));
    mp->status = i_alloc_1d(df->nv);

    /* First variable is always time */
    m = 0;
    df->variables[m].type = VT_DATA;
    v = &df->variables[m];
    v->varid = m;   /* id = 0 */
    v->dim_as_record = 0;
    strcpy(coords, "time");
    v->name = (char *)malloc(sizeof(char) * MAXLINELEN);
    v->longname = (char *)malloc(sizeof(char) * MAXLINELEN);
    v->units = (char *)malloc(sizeof(char) * MAXLINELEN);
    strcpy(v->name, coords);
    strcpy(v->longname, coords);
    strcpy(v->units, mp->tunits);
    v->nd = 1;
    v->dimids = (int *)malloc(v->nd * sizeof(int));
    v->dimids[0] = 0;
    decode_coord_type(coords, &v->type, &v->coord_domain);
    v->na = 0;
    v->scale_factor = 1.0;
    v->add_offset = 0.0;
    m++;

    /* z position */
    df->variables[m].type = VT_DATA;
    v = &df->variables[m];
    v->varid = m;   /* id = 1 */
    v->dim_as_record = 0;
    strcpy(coords, "z");
    v->name = (char *)malloc(sizeof(char) * MAXLINELEN);
    v->longname = (char *)malloc(sizeof(char) * MAXLINELEN);
    strcpy(v->name, coords);
    strcpy(v->longname, coords);
    v->nd = 1;
    v->dimids = (int *)malloc(v->nd * sizeof(int));
    v->dimids[0] = 2;
    if (strlen(xyzunits[2])) {
      v->units = (char *)malloc(sizeof(char) * MAXLINELEN);
      strcpy(v->units, xyzunits[2]);
    }
    v->na = 0;
    v->scale_factor = 1.0;
    v->add_offset = 0.0;
    decode_coord_type(coords, &v->type, &v->coord_domain);
    m++;

    /* x position */
    v = &df->variables[m];
    df->variables[m].type = VT_DATA;
    v->varid = m;   /* id = 2 */
    v->dim_as_record = 0;
    if ((strcmp(xyzunits[0], "degrees_east") == 0))
      strcpy(coords, "longitude");
    else
      strcpy(coords, "x");
    v->name = (char *)malloc(sizeof(char) * MAXLINELEN);
    v->longname = (char *)malloc(sizeof(char) * MAXLINELEN);
    strcpy(v->name, coords);
    strcpy(v->longname, coords);
    v->nd = 1;
    v->dimids = (int *)malloc(v->nd * sizeof(int));
    v->dimids[0] = 1;
    if (strlen(xyzunits[0])) {
      v->units = (char *)malloc(sizeof(char) * MAXLINELEN);
      strcpy(v->units, xyzunits[0]);
    }
    v->na = 0;
    v->scale_factor = 1.0;
    v->add_offset = 0.0;
    decode_coord_type(coords, &v->type, &v->coord_domain);
    m++;

    /* y position */
    v = &df->variables[m];
    df->variables[m].type = VT_DATA;
    v->varid = m;   /* id = 3 */
    v->dim_as_record = 0;
    if ((strcmp(xyzunits[1], "degrees_north") == 0))
      strcpy(coords, "latitude");
    else
      strcpy(coords, "y");
    v->name = (char *)malloc(sizeof(char) * MAXLINELEN);
    v->longname = (char *)malloc(sizeof(char) * MAXLINELEN);
    strcpy(v->name, coords);
    strcpy(v->longname, coords);
    v->nd = 1;
    v->dimids = (int *)malloc(v->nd * sizeof(int));
    v->dimids[0] = 1;
    if (strlen(xyzunits[1])) {
      v->units = (char *)malloc(sizeof(char) * MAXLINELEN);
      strcpy(v->units, xyzunits[1]);
    }
    v->na = 0;
    v->scale_factor = 1.0;
    v->add_offset = 0.0;
    decode_coord_type(coords, &v->type, &v->coord_domain);
    m++;

    /* Remaining variables */
    for (i = m; i < df->nv; ++i) {
      df_variable_t *v = &df->variables[i];
      j = i - MPK_SVARS;

      /* df_variable_t type and netCDF id - default to data */
      df->variables[i].type = VT_DATA;
      v->varid = i;
      v->dim_as_record = 0;

      /* df_variable_t name */
      v->name = (char *)malloc(sizeof(char) * MAXLINELEN);
      strcpy(v->name, vars[j]);
      v->longname = (char *)malloc(sizeof(char) * MAXLINELEN);
      strcpy(v->longname, v->name);
      /* Extract the dimension information for each variable */
      if (j < mp->n2d) {
	v->nd = 2;
	v->dimids = (int *)malloc(v->nd * sizeof(int));
	v->dimids[0] = 0;
	v->dimids[1] = 1;
      } else {
	v->nd = 3;
	v->dimids = (int *)malloc(v->nd * sizeof(int));
	v->dimids[0] = 0;
	v->dimids[1] = 2;
	v->dimids[2] = 1;
      }
      v->units = (char *)malloc(sizeof(char) * MAXLINELEN);
      strcpy(v->units, units[j]);
      v->na = 0;
      v->scale_factor = 1.0;
      v->add_offset = 0.0;
    }

    /* No global attributes */
    df->na = 0;

    /* Time is always set to variable 0.
       This will become the coordinate dimension. */
    df_set_record(df, 0);

    /* Sweep back trough the variables and decode the associated
       coordinates information. This has to be left until the end to 
       ensure all variables were read. */
    for (i = MPK_SVARS; i < df->nv; ++i) {
      df_variable_t *v = &df->variables[i];
      char *text = NULL;

      j = i - MPK_SVARS;
      if (j < mp->n2d) {
	if ((strcmp(xyzunits[0], "degrees_east") == 0) &&
	    (strcmp(xyzunits[1], "degrees_north") == 0))
	  strcpy(coords, "latitude, longitude");
	else
	  strcpy(coords, "y, x");
      } else {
	if ((strcmp(xyzunits[0], "degrees_east") == 0) &&
	    (strcmp(xyzunits[1], "degrees_north") == 0) &&
	    (strcmp(xyzunits[2], "m") == 0))
	  strcpy(coords, "z, latitude, longitude");
	else
	  strcpy(coords, "z, y, x");
      }
      text = coords;
      decode_coords(df, v, text);
    }
  } else {
    warn("No variables were specified in this netCDF file.\n");
    return;
  }

}


void mempack_read_records(datafile_t *df, int varid) {
  df_variable_t *rv = &df->variables[varid];
  df_mempack_t *mp = (df_mempack_t *)df->private_data;
  int i, recdimid = rv->dimids[0];

  FILE *fp;
  while ((fp = fopen(mp->id, "r")) == NULL);
  prm_read_double(fp, "time", &mp->time);
  fclose(fp);

  if (rv->data == NULL)
    rv->data = (double *)malloc(sizeof(double) * 
				df->dimensions[recdimid].size);

  if (mp->time > rv->data[1]) {
    rv->data[0] = rv->data[1];
    rv->data[1] = mp->time;
    for (i = 0; i < df->nv; i++)
      mp->status[i] = DF_GO;
  }
}


/* Read a data record into the specified variable from a memory
 * packet. It is assumed that the memory has already been allocated.
 */
void mempack_read_data(datafile_t *df, df_variable_t *v, int r, int roffset)
{
  int start[2];
  int count[2];
  int i = 0, j, n;
  int index = r - v->start_record;
  df_mempack_t *mp = (df_mempack_t *)df->private_data;

  mempack_read_file(mp->id, df);

  if (!(df->rec_modulus) && ((r < 0) || (r >= df->nrecords)))
    quit
      ("mempack_read_data: Attempt to read an invalid record for variable '%s'.\n",
       v->name);

  if (v->dim_as_record) {
    if (df->rec_modulus)
      start[i] = (r % df->nrecords) - roffset;
    else
      start[i] = r - roffset;
    count[i++] = 1;
  }

  count[0] = df->dimensions[v->dimids[0]].size;
  if (v->varid == 1) {      /* z position data */
    count[1] = 2 * mp->npoints;
    for (i = 0; i < count[0]; i++)
      VAR_1D(v)[index][i] = mp->xyz[i + count[1]];
  }
  else if (v->varid == 2) { /* x position data */
    count[1] = 0;
    for (i = 0; i < count[0]; i++)
      VAR_1D(v)[index][i] = mp->xyz[i + count[1]];
  }
  else if (v->varid == 3) { /* y position data */
    count[1] = mp->npoints;
    for (i = 0; i < count[0]; i++)
      VAR_1D(v)[index][i] = mp->xyz[i + count[1]];
  } else {                  /* Other variables */
    switch (v->nd) {
    case 1:                 /* 1D data */
      count[0] = df->dimensions[v->dimids[0]].size;
      n = (v->varid - MPK_SVARS) * count[0];
      for (i = 0; i < count[0]; i++)
	VAR_1D(v)[index][i] = mp->v2d[n++];
      break;

    case 2:                 /* 2D data */
      count[0] = df->dimensions[v->dimids[0]].size;
      count[1] = df->dimensions[v->dimids[1]].size;
      n = (v->varid - mp->n2d - MPK_SVARS) * count[0] * count[1];
      for (i = 0; i < count[1]; ++i)
	for (j = 0; j < count[0]; ++j)
	  VAR_2D(v)[index][j][i] = mp->v3d[n++];
      break;

    default:
      quit("mempack_read_data: Bad number of dimensions\n");
    }
  }
}

/*
 * Free up the memory associated with a memory packet. */
void mempack_free(datafile_t *df)
{
  df_mempack_t *mp = (df_mempack_t *)df->private_data;
  d_free_1d(mp->xyz);
  if (mp->status) i_free_1d(mp->status);
  if (mp->n2d) d_free_1d(mp->v2d);
  if (mp->n3d) d_free_1d(mp->v3d);
  free(mp);
}

int mempack_read_file(char *mpname, datafile_t *df)
{
  int m, n;
  FILE *fp;
  df_mempack_t *data = (df_mempack_t *)df->private_data;

  while ((fp = fopen(mpname, "r")) == NULL);
  while (!(prm_skip_to_end_of_key(fp, "end")));

  prm_set_errfn(quit);
  prm_read_double(fp, "multi-netcdf-version", &data->version);
  prm_read_int(fp, "runcode", &data->runcode);
  prm_read_int(fp, "npoints", &data->npoints);
  prm_read_int(fp, "nz", &data->nz);
  prm_read_int(fp, "n2d", &data->n2d);
  prm_read_int(fp, "n3d", &data->n3d);
  prm_read_double(fp, "time", &data->time);
  prm_read_double(fp, "dt", &data->dt);
  prm_read_char(fp, "tunits", data->tunits);
  prm_read_char(fp, "xyzunits", data->xyzunits);
  prm_read_char(fp, "vars", data->vars);
  prm_read_char(fp, "units", data->units);
  prm_read_darray(fp, "xyzdata", &data->xyz, &m);
  if (data->n2d)
    prm_read_darray(fp, "var2d", &data->v2d, &m);
  if (data->n3d)
    prm_read_darray(fp, "var3d", &data->v3d, &m);
  fclose(fp);
  return(0);
}


/*
Routine to read column data from an ascii file

Arguments:
    fp      -    FILE pointer for ascii input file
    df      -    pointer to data file structure (assumed
	         not previously initialised)

This routine calls quit() if anything goes wrong
*/
void ascii_read(FILE * fp, datafile_t *df)
{
  int n;
  int i;
  char line[MAXLINELEN];
  char *str[MAXLINELEN];
  char key[MAXLINELEN];

  df->type = DFT_ASCII;

  /* Set the dimensions there can only be 1 */
  df->dimensions = (df_dimension_t *)malloc(sizeof(df_dimension_t));
  memset(df->dimensions, 0, sizeof(df_dimension_t));
  df->dimensions[0].name = (char *)malloc(sizeof(char) * MAXLINELEN);
  strcpy(df->dimensions[0].name, "row");


  /* Get number of variables, allocate memory and zero the array */
  prm_set_errfn(quit);
  prm_skip_to_end_of_key(fp, "## COLUMNS");
  if (!fscanf(fp, "%d", &n))
    quit("ascii_read: number of columns in %s not found\n", df->name);

  df->variables = (df_variable_t *)malloc(n * sizeof(df_variable_t));
  memset(df->variables, 0, n * sizeof(df_variable_t));

  /* Populate the arrays with meta data */
  if (n < 1)
    quit("ascii_read: No variables! (reading file %s)\n", df->name);

  /* Read variable info */
  for (i = 0; i < n; i++) {

    /* df_variable_t type - default to data */
    df->variables[i].type = VT_DATA;

    /* df_variable_t name */
    sprintf(key, "## COLUMN%d.name", i + 1);
    prm_read_char(fp, key, line);
    if ((df->variables[i].name = (char *)malloc(strlen(line) + 1)) == NULL)
      quit("ascii_read: No memory for variable names\n");
    strcpy(df->variables[i].name, line);
    prm_set_errfn(warn);

    /* df_variable_t long name */
    sprintf(key, "## COLUMN%d.long_name", i + 1);
    prm_read_char(fp, key, line);
    if ((df->variables[i].longname =
         (char *)malloc(strlen(line) + 1)) == NULL)
      quit("ascii_read: No memory for long variable names\n");
    strcpy(df->variables[i].longname, line);

    /* df_variable_t units */
    sprintf(key, "## COLUMN%d.units", i + 1);
    prm_read_char(fp, key, line);
    if ((df->variables[i].units =
         (char *)malloc(strlen(line) + 1)) == NULL)
      quit("ascii_read: No memory for variable units when reading file %s\n", df->name);
    strcpy(df->variables[i].units, line);

    prm_set_errfn(quiet);
    /* Missing values */
    sprintf(key, "## COLUMN%d.missing_value", i + 1);
    prm_read_double(fp, key, &df->variables[i].missing);

    /* Fill value values */
    sprintf(key, "## COLUMN%d.fill_value", i + 1);
    prm_read_double(fp, key, &df->variables[i].fillvalue);

    /* Check for modulo */
    sprintf(key, "## COLUMN%d.modulo", i + 1);
    if (prm_read_char(fp, key, line)) {
       df_variable_t *v = &df->variables[i];
       v->na = 1;
       v->attributes = (df_attribute_t *)malloc(v->na * sizeof(df_attribute_t));
       memset(v->attributes, 0, v->na * sizeof(df_attribute_t));
       v->attributes[0].name = strdup("modulo");
       v->attributes[0].n = strlen(line);
       v->attributes[0].type = AT_TEXT;
       v->attributes[0].value = (void *)strdup(line);
    }

    prm_set_errfn(quit);
  }


  /* Store number of variables and pointer to record units */
  df->nv = n;
  df->ri = -1;

  /* Count number of data lines */
  fseek(fp, 0L, 0);
  for (n = 0; prm_next_line(line, MAXLINELEN, fp); n++) /* loop */
    ;
  if (n < 1)
    quit("ascii_read: No data in file %s\n", df->name);
  df->dimensions[0].size = n;
  
  emslog(LTRACE,"reading datarows: %d \n",n);

  /* Allocate data memory and store pointer to record value (if applic) */
  for (i = 0; i < df->nv; ++i) {
    df->variables[i].nd = 1;
    df->variables[i].dimids =
      (int *)malloc(sizeof(int) * df->variables[i].nd);
    df->variables[i].dimids[0] = 0;
    df->variables[i].start_record = 0;
    df->variables[i].nrecords = 1;
    df->variables[i].dim_as_record = 0;
    df->variables[i].scale_factor = 1.0;
    df->variables[i].add_offset = 0.0;
    df->variables[i].data =
      (double *)malloc(df->dimensions[0].size * sizeof(double));
  }


  /* Read data */
  fseek(fp, 0L, 0);
  for (n = 0; prm_next_line(line, MAXLINELEN, fp); n++) {
    /* Check for correct number of strings on each line */
    if (parseline(line, str, MAXLINELEN) != df->nv)
      quit("ascii_read: Wrong number of data values on line %d of %s\n", n, df->name);
    /* Loop over each value */
    for (i = 0; i < df->nv; i++) {
      if (sscanf(str[i], "%lf", &df->variables[i].data[n]) != 1)
        quit("ascii_read: Can't read data value for column %d, line %d of %s\n",
             i + 1, n, df->name);
    }
  }

  fclose(fp);
}


void ascii_write(FILE * fp, datafile_t *df)
{
  int i, n;

  /* Write variable info */
  fprintf(fp, "## COLUMNS %d\n##\n", df->nv);
  for (i = 0; i < df->nv; i++) {
    fprintf(fp, "## COLUMN%d.name  %s\n", i + 1, df->variables[i].name);
    fprintf(fp, "## COLUMN%d.long_name  %s\n", i + 1,
            df->variables[i].longname);
    fprintf(fp, "## COLUMN%d.units  %s\n", i + 1, df->variables[i].units);
    fprintf(fp, "## COLUMN%d.missing_value  %g\n", i + 1,
            df->variables[i].missing);
    fprintf(fp, "## COLUMN%d.fill_value  %g\n", i + 1,
            df->variables[i].fillvalue);
    fprintf(fp, "##\n");
  }

  /* Write values */
  for (i = 0; i < df->nrecords; i++) {
    for (n = 0; n < df->nv; n++) {
      fprintf(fp, "%.10g", df->variables[n].data[i]);
      if (n < df->nv - 1)
        fprintf(fp, " ");
    }
    fprintf(fp, "\n");
  }
}


/*
 * Free up the memory associated with an ascii column data file.
 */
void ascii_free(datafile_t *df)
{
}




/*
Routine to read an attribute from the netCDF file.

Arguments:
    df          -    pointer to data file structure
    varid       -    df_variable_t id or NC_GLOBAL.
    attnum      -    df_attribute_t number.
    att         -    df_attribute_t structure to be populated.

This routine calls quit() if anything goes wrong
*/
void netcdf_read_attrib(datafile_t *df, int varid, int attnum,
                        df_attribute_t *a)
{
  char line[MAXLINELEN];
  nc_type type;

  /* Read attribute name */
  nc_inq_attname(df->ncid, varid, attnum, line);
  a->name = (char *)malloc((sizeof(line) + 1) * sizeof(char));
  strcpy(a->name, line);

  /* Get the attribute type and length */
  nc_inq_att(df->ncid, varid, a->name, &type, &a->n);
  if (type == NC_CHAR) {
    a->value = (void *)malloc(sizeof(char) * (a->n + 1));
    memset(a->value, 0, sizeof(char) * (a->n + 1));
    nc_get_att_text(df->ncid, varid, a->name, (char *)a->value);
    a->type = AT_TEXT;
  }

  else if (a->n > 0) {
    switch (type) {
    case NC_BYTE:
      a->value = (void *)malloc(sizeof(char) * a->n);
      nc_get_att_schar(df->ncid, varid, a->name, (signed char *)a->value);
      a->type = AT_BYTE;
      break;

    case NC_FLOAT:
      a->value = (void *)malloc(sizeof(float) * a->n);
      nc_get_att_float(df->ncid, varid, a->name, (float *)a->value);
      a->type = AT_FLOAT;
      break;

    case NC_DOUBLE:
      a->value = (void *)malloc(sizeof(double) * a->n);
      nc_get_att_double(df->ncid, varid, a->name, (double *)a->value);
      a->type = AT_DOUBLE;
      break;

    case NC_SHORT:
      a->value = (void *)malloc(sizeof(short) * a->n);
      nc_get_att_short(df->ncid, varid, a->name, (short *)a->value);
      a->type = AT_SHORT;
      break;

    case NC_INT:
      a->value = (void *)malloc(sizeof(int) * a->n);
      nc_get_att_int(df->ncid, varid, a->name, (int *)a->value);
      a->type = AT_INT;
      break;

    default:
      break;
    }
  }
}


/* Allocate a data record. */
void alloc_data_record(datafile_t *df, df_variable_t *v, int r)
{
  int e1 = 0;
  int e2 = 0;
  int e3 = 0;
  int e4 = 0;

  switch (v->nd) {
  case 0:
    break;

  case 1:
    e1 = df->dimensions[v->dimids[0]].size;
    VAR_1D(v)[r] = d_alloc_1d(e1);
    break;

  case 2:
    e1 = df->dimensions[v->dimids[0]].size;
    e2 = df->dimensions[v->dimids[1]].size;
    VAR_2D(v)[r] = d_alloc_2d(e2, e1);
    break;

  case 3:
    e1 = df->dimensions[v->dimids[0]].size;
    e2 = df->dimensions[v->dimids[1]].size;
    e3 = df->dimensions[v->dimids[2]].size;
    VAR_3D(v)[r] = d_alloc_3d(e3, e2, e1);
    break;

  case 4:
    e1 = df->dimensions[v->dimids[0]].size;
    e2 = df->dimensions[v->dimids[1]].size;
    e3 = df->dimensions[v->dimids[2]].size;
    e4 = df->dimensions[v->dimids[3]].size;
    VAR_4D(v)[r] = d_alloc_4d(e4, e3, e2, e1);
    break;

  default:
    quit
      ("alloc_data_record: Unable to allocate memory for '%d' dimensional array.\n",
       v->nd);
  }
}


double *alloc_data(datafile_t *df, df_variable_t *v)
{
  double *data = 0;

  switch (v->nd) {
  case 0:
    data = (double *)malloc(sizeof(double) * v->nrecords);
    break;

  case 1:
  case 2:
  case 3:
  case 4:
    data = (double *)malloc(sizeof(double *) * v->nrecords);
    break;

  default:
    quit
      ("alloc_data: Unable to allocate memory for '%d' dimensional array.\n",
       v->nd);
  }

  return data;
}


void alloc_data_records(datafile_t *df, df_variable_t *v, int nrecs)
{
  int i = 0;

  v->nrecords = nrecs;
  v->data = alloc_data(df, v);
  for (i = 0; i < nrecs; ++i)
    alloc_data_record(df, v, i);
}



void free_data_record(datafile_t *df, df_variable_t *v, int r)
{

  if (v->data != NULL) {
    switch (v->nd) {
    case 0:
      break;

    case 1:
      if (VAR_1D(v)[r] != NULL)
        d_free_1d(VAR_1D(v)[r]);
      VAR_1D(v)[r] = NULL;
      break;

    case 2:
      if (VAR_2D(v)[r] != NULL)
        d_free_2d(VAR_2D(v)[r]);
      VAR_2D(v)[r] = NULL;
      break;

    case 3:
      if (VAR_3D(v)[r] != NULL)
        d_free_3d(VAR_3D(v)[r]);
      VAR_3D(v)[r] = NULL;
      break;

    case 4:
      if (VAR_4D(v)[r] != NULL)
        d_free_4d(VAR_4D(v)[r]);
      VAR_4D(v)[r] = NULL;
      break;
    }
  }
}


void free_data_records(datafile_t *df, df_variable_t *v)
{
  int i;

  if (v->data != NULL) {
    for (i = 0; i < v->nrecords; ++i)
      free_data_record(df, v, i);
    free(v->data);
  } 

  v->data = NULL;
  v->nrecords = 0;
  v->start_record = 0;
}



/* Read a data record into the specified variable from a netCDF
 * file. It is assumed that the memory has already been allocated.
 */
/*
void read_data(datafile_t *df, df_variable_t *v, int r)
{
  size_t start[5];
  size_t count[5];
  int i = 0, j, k, l;
  int fid = df->ncid;
  int index = r - v->start_record;

  if (!(df->rec_modulus) && ((r < 0) || (r >= df->nrecords)))
    quit
      ("read_data: Attempt to read an invalid record for variable '%s'.\n",
       v->name);

  if (v->dim_as_record) {
    if (df->rec_modulus)
      start[i] = r % df->nrecords;
    else
      start[i] = r;
    count[i++] = 1;
  }

  switch (v->nd) {
  case 0:
    nc_get_vara_double(fid, v->varid, start, count, &v->data[index]);
    v->data[index] = v->data[index] * v->scale_factor + v->add_offset;
    break;

  case 1:
    start[i] = 0;
    count[i++] = df->dimensions[v->dimids[0]].size;
    nc_get_vara_double(fid, v->varid, start, count, VAR_1D(v)[index]);
    for (i = 0; i < (int)df->dimensions[v->dimids[0]].size; ++i)
      VAR_1D(v)[index][i]
        = VAR_1D(v)[index][i] * v->scale_factor + v->add_offset;
    break;

  case 2:
    start[i] = 0;
    count[i++] = df->dimensions[v->dimids[0]].size;
    start[i] = 0;
    count[i++] = df->dimensions[v->dimids[1]].size;
    nc_get_vara_double(fid, v->varid, start, count, VAR_2D(v)[index][0]);
    for (i = 0; i < (int)df->dimensions[v->dimids[0]].size; ++i)
      for (j = 0; j < (int)df->dimensions[v->dimids[1]].size; ++j)
        VAR_2D(v)[index][i][j]
          = VAR_2D(v)[index][i][j] * v->scale_factor + v->add_offset;
    break;

  case 3:
    start[i] = 0;
    count[i++] = df->dimensions[v->dimids[0]].size;
    start[i] = 0;
    count[i++] = df->dimensions[v->dimids[1]].size;
    start[i] = 0;
    count[i++] = df->dimensions[v->dimids[2]].size;
    nc_get_vara_double(fid, v->varid, start, count,
                       VAR_3D(v)[index][0][0]);
    for (i = 0; i < (int)df->dimensions[v->dimids[0]].size; ++i)
      for (j = 0; j < (int)df->dimensions[v->dimids[1]].size; ++j)
        for (k = 0; k < (int)df->dimensions[v->dimids[2]].size; ++k)
          VAR_3D(v)[index][i][j][k]
            = VAR_3D(v)[index][i][j][k] * v->scale_factor + v->add_offset;
    break;

  case 4:
    start[i] = 0;
    count[i++] = df->dimensions[v->dimids[0]].size;
    start[i] = 0;
    count[i++] = df->dimensions[v->dimids[1]].size;
    start[i] = 0;
    count[i++] = df->dimensions[v->dimids[2]].size;
    start[i] = 0;
    count[i++] = df->dimensions[v->dimids[3]].size;
    nc_get_vara_double(fid, v->varid, start, count,
                       VAR_4D(v)[index][0][0][0]);
    for (i = 0; i < (int)df->dimensions[v->dimids[0]].size; ++i)
      for (j = 0; j < (int)df->dimensions[v->dimids[1]].size; ++j)
        for (k = 0; k < (int)df->dimensions[v->dimids[2]].size; ++k)
          for (l = 0; l < (int)df->dimensions[v->dimids[3]].size; ++l)
            VAR_4D(v)[index][i][j][k][l]
              = VAR_4D(v)[index][i][j][k][l] * v->scale_factor +
              v->add_offset;
    break;

  default:
    quit("read_data: Bad number of dimensions\n");
  }
}
*/

/* Shift the data to the left by nrec places */
void left_shift_data(datafile_t *df, df_variable_t *v, int nrecs)
{
  int i;
  double *newdata = alloc_data(df, v);

  for (i = 0; i < v->nrecords; ++i) {
    int j = (i - nrecs + v->nrecords) % v->nrecords;

    /* Copy the data */
    switch (v->nd) {
    case 0:
      newdata[j] = v->data[i];
      break;

    case 1:
      ((double **)newdata)[j] = VAR_1D(v)[i];
      break;

    case 2:
      ((double ***)newdata)[j] = VAR_2D(v)[i];
      break;

    case 3:
      ((double ****)newdata)[j] = VAR_3D(v)[i];
      break;

    case 4:
      ((double *****)newdata)[j] = VAR_4D(v)[i];
      break;
    }
  }

  free(v->data);
  v->data = newdata;
}


/* Shift the data to the right by nrec places */
void right_shift_data(datafile_t *df, df_variable_t *v, int nrecs)
{
  int i;
  double *newdata = alloc_data(df, v);

  for (i = 0; i < v->nrecords; ++i) {
    int j = (i + nrecs + v->nrecords) % v->nrecords;

    /* Copy the data */
    switch (v->nd) {
    case 0:
      newdata[j] = v->data[i];
      break;

    case 1:
      ((double **)newdata)[j] = VAR_1D(v)[i];
      break;

    case 2:
      ((double ***)newdata)[j] = VAR_2D(v)[i];
      break;

    case 3:
      ((double ****)newdata)[j] = VAR_3D(v)[i];
      break;

    case 4:
      ((double *****)newdata)[j] = VAR_4D(v)[i];
      break;
    }
  }

  free(v->data);
  v->data = newdata;
}


/* Shift the data record along. A negative number of records implies
 * a shift to the left, a positive to the right.
 */
void shift_data(datafile_t *df, df_variable_t *v, int nrecs)
{
  int i = 0;
  int j = 0;
  int src_start = 0;
  int src_end = 0;
  int dst_start = 0;
  int dir = (nrecs < 0) ? 1 : -1;
  int nstart = v->start_record;

  nrecs = abs(nrecs);
  /* Only scroll if necessary */
  if ((nrecs == 0) || (!v->dim_as_record) || (v->data == NULL))
    return;

  /* Adjust so we don't scroll past the ends. */
  nstart += dir * nrecs;
  if (nstart < 0) {
    nrecs += nstart;
    nstart = 0;
  } else if ((nstart + v->nrecords) > df->nrecords) {
    nrecs -= (df->nrecords - (nstart + v->nrecords));
    nstart = df->nrecords - v->nrecords - 1;
  }

  if (nrecs < v->nrecords) {
    if (dir > 0) {
      src_start = nrecs;
      src_end = v->nrecords - 1;
      dst_start = 0;
    } else {
      src_start = v->nrecords - nrecs - 1;
      src_end = 0;
      dst_start = v->nrecords - 1;
    }

    for (i = src_start, j = dst_start;
         (dir < 0) ? (i >= src_end) : (i <= src_end); i += dir, j += dir) {

      /* Swap the values or addresses */
      switch (v->nd) {
      case 0:
        {
          double tmp = v->data[j];
          v->data[j] = v->data[i];
          v->data[i] = tmp;
        }
        break;

      case 1:
        {
          double *tmp = VAR_1D(v)[j];
          VAR_1D(v)[j] = VAR_1D(v)[i];
          VAR_1D(v)[i] = tmp;
        }
        break;

      case 2:
        {
          double **tmp = VAR_2D(v)[j];
          VAR_2D(v)[j] = VAR_2D(v)[i];
          VAR_2D(v)[i] = tmp;
        }
        break;

      case 3:
        {
          double ***tmp = VAR_3D(v)[j];
          VAR_3D(v)[j] = VAR_3D(v)[i];
          VAR_3D(v)[i] = tmp;
        }
        break;

      case 4:
        {
          double ****tmp = VAR_4D(v)[j];
          VAR_4D(v)[j] = VAR_4D(v)[i];
          VAR_4D(v)[i] = tmp;
        }
        break;
      }
    }
  } else
    quit("shift_data: Attempt to shift more records that in buffer.\n");

  v->start_record = nstart;
}



/*
 * Local functions
 */
static void destroy_hq_within_ht(void *ht_void)
{
  // special case destroy function that destroys the hash-queue with
  // the given hash-table.
  hash_table_t *ht = (hash_table_t *)ht_void;
  ht_process(ht, hq_destroy);
}

static int has_name(char *name, char *array[])
{
  int j = 0;
  while (array[j] != NULL) {
    if (strcmp(name, array[j]) == 0) {
      return(1);
    }
    j++;
  }
  return(0);
}

// EOF


