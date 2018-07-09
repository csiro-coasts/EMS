/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/io/ptrack.c
 *
 *  \brief Routines for particle tracking
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: ptrack.c 5833 2018-06-27 00:21:35Z riz008 $
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "netcdf.h"
#include "ems.h"


/** Read an array of particles from a netCDF particle file
  * at a specified record.
  *
  * @param name netCDF particle filename.
  * @param rec record number.
  * @param np returned number of particles.
  * @param p returned pointer to a particle structure.
  * @param t returned time for specified record (NULL if not required).
  * @param t_units returned time units (NULL if not required).
  * @param ndump returned number of dumps in particle file (NULL if not required).
  */
void
pt_read(char *name, int rec, long int *np, particle_t **p,
        double *t, char *t_units, int *ndump)
{
  int fid;
  int ndims;                    /* Number of dimensions */
  int nvars;                    /* Number of variables */
  int natts;                    /* Number of attributes */
  int recdim;                   /* ID of unlimited dimension */
  nc_type dt;                   /* type of variable */
  int t_did;                    /* Time dimension id */
  int p_did = -1;               /* p dimension id */
  int dims[NC_MAX_DIMS];
  size_t start[2];
  size_t count[2];
  size_t n;
  int t_vid;
  int x_vid;
  int y_vid;
  int z_vid;
  int f_vid;
  double *b;
  short *f;
  size_t nrec;

  /* Open the netCDF file */
  if ((nc_open(name, NC_NOWRITE, &fid)) != NC_NOERR)
    quit("pt_read: Can't open %s\n", name);

  /* Inquire about this file */
  if (nc_inq(fid, &ndims, &nvars, &natts, &recdim) == -1)
    quit("pt_read: nc_inq failed\n");
  if (ndims != 2)
    quit("pt_read: not enough dimensions\n");
  if (nvars < 5)
    quit("pt_read: not enough variables\n");

  /* Check dimensions are as expected */
  if ((nc_inq_dimid(fid, "t", &t_did)) != NC_NOERR)
    quit("pt_read: no t dimension");
  if (t_did != recdim)
    quit("pt_read: t dimension not unlimited\n");
  if ((nc_inq_dimid(fid, "n", &p_did)) != NC_NOERR)
    quit("pt_read: no n dimension");
  /* Get number of particles */
  nc_inq_dimlen(fid, p_did, &n);

  /* Check that there is a variable called "t" */
  if ((t_vid = ncw_var_id(fid, "t")) < 0)
    quit("pt_read: No t variable\n");
  nc_inq_var(fid, t_vid, NULL, &dt, &ndims, dims, &natts);
  if (ndims != 1 || dims[0] != t_did)
    quit("pt_read: t variable has wrong dimensions\n");
  if (dt != NC_DOUBLE)
    quit("pt_read: t variable must have type NC_DOUBLE\n");
  /* Get number of records and check that requested record exists */
  nc_inq_dimlen(fid, t_did, &nrec);
  if (rec >= (int)nrec)
    quit("pt_read: Record %d not in %s (only %d records)\n", rec, name,
         nrec);
  if (ndump)
    *ndump = (int)nrec;

  /* Check that there is a variable called "x" */
  if ((x_vid = ncw_var_id(fid, "x")) < 0)
    quit("pt_read: No x variable\n");
  nc_inq_var(fid, x_vid, NULL, &dt, &ndims, dims, &natts);
  if (ndims != 2 || dims[0] != t_did || dims[1] != p_did)
    quit("pt_read: x variable has wrong dimensions\n");
  if (dt != NC_DOUBLE)
    quit("pt_read: x variable must have type NC_DOUBLE\n");

  /* Check that there is a variable called "y" */
  if ((y_vid = ncw_var_id(fid, "y")) < 0)
    quit("pt_read: No y variable\n");
  nc_inq_var(fid, y_vid, NULL, &dt, &ndims, dims, &natts);
  if (ndims != 2 || dims[0] != t_did || dims[1] != p_did)
    quit("pt_read: y variable has wrong dimensions\n");
  if (dt != NC_DOUBLE)
    quit("pt_read: y variable must have type NC_DOUBLE\n");

  /* Check that there is a variable called "z" */
  if ((z_vid = ncw_var_id(fid, "z")) < 0)
    quit("pt_read: No z variable\n");
  nc_inq_var(fid, z_vid, NULL, &dt, &ndims, dims, &natts);
  if (ndims != 2 || dims[0] != t_did || dims[1] != p_did)
    quit("pt_read: z variable has wrong dimensions\n");
  if (dt != NC_DOUBLE)
    quit("pt_read: z variable must have type NC_DOUBLE\n");

  /* Check that there is a variable called "flag" */
  if ((f_vid = ncw_var_id(fid, "flag")) < 0)
    quit("pt_read: No flag variable\n");
  nc_inq_var(fid, f_vid, NULL, &dt, &ndims, dims, &natts);
  if (ndims != 2 || dims[0] != t_did || dims[1] != p_did)
    quit("pt_read: flag variable has wrong dimensions\n");
  if (dt != NC_SHORT)
    quit("pt_read: flag variable must have type NC_SHORT\n");

  /* Allocate space, if not already done */
  if (*p == NULL) {
    *np = (long)n;
    if ((*p = (particle_t *)malloc((*np) * sizeof(particle_t))) == NULL)
      quit("pt_read: not enough memory for particles\n");
  } else if (*np != (int)n)
    quit
      ("pt_read: Number of particles doesn't match space already allocated\n");

  /* Read time information if required */
  start[0] = rec;
  count[0] = 1;
  if (t)
    nc_get_vara_double(fid, t_vid, start, count, t);
  if (t_units)
    nc_get_att_text(fid, t_vid, "units", t_units);

  /* Read the particle data */
  b = d_alloc_1d(*np);
  f = s_alloc_1d(*np);
  start[0] = rec;
  start[1] = 0;
  count[0] = 1;
  count[1] = *np;
  nc_get_vara_double(fid, x_vid, start, count, b);
  for (n = 0; n < (size_t) * np; n++)
    (*p)[n].e1 = b[n];
  nc_get_vara_double(fid, y_vid, start, count, b);
  for (n = 0; n < (size_t) * np; n++)
    (*p)[n].e2 = b[n];
  nc_get_vara_double(fid, z_vid, start, count, b);
  for (n = 0; n < (size_t) * np; n++)
    (*p)[n].e3 = b[n];
  nc_get_vara_short(fid, f_vid, start, count, f);
  for (n = 0; n < (size_t) * np; n++)
    (*p)[n].flag = f[n];
  d_free_1d(b);
  s_free_1d(f);

  /* Close the file */
  nc_close(fid);
}

/** Create particle tracking output file.
  *
  * @param name particle filename.
  * @param np number of particles per record dump.
  * @param t_units time units of each dump.
  * @param dumpf variables to include in dumpfile.
  * @return netCDF file descriptor.
  */
int pt_create(char *name, long int np, char *t_units, int dumpf)
{
  int ncid;                     /* netCDF id */
  /* dimension ids */
  int t_dim, n_dim;
  /* variable ids */
  int t_id, x_id, y_id, z_id, flag_id;
  int age_id, size_id;
  /* variable shapes */
  int dims[2];

  /* enter define mode */
  if (nc_create(name, NC_NOCLOBBER, &ncid) != NC_NOERR)
    quit("pt_create: Unable to create the particle file '%s'\n", name);

  /* define dimensions */
  nc_def_dim(ncid, "t", NC_UNLIMITED, &t_dim);
  nc_def_dim(ncid, "n", np, &n_dim);

  /* define variables */
  dims[0] = t_dim;
  nc_def_var(ncid, "t", NC_DOUBLE, 1, dims, &t_id);
  dims[1] = n_dim;
  nc_def_var(ncid, "x", NC_DOUBLE, 2, dims, &x_id);
  nc_def_var(ncid, "y", NC_DOUBLE, 2, dims, &y_id);
  nc_def_var(ncid, "z", NC_DOUBLE, 2, dims, &z_id);
  nc_def_var(ncid, "flag", NC_SHORT, 2, dims, &flag_id);
  if (dumpf & PT_AGE)
    nc_def_var(ncid, "age", NC_BYTE, 2, dims, &age_id);
  if (dumpf & PT_SIZE)
    nc_def_var(ncid, "size", NC_BYTE, 2, dims, &size_id);

  /* Time variable units attribute */
  nc_put_att_text(ncid, t_id, "units", strlen(t_units), t_units);

  /* leave define mode */
  nc_enddef(ncid);

  return (ncid);
}


/** Write particles at a known record into the particle file.
  *
  * @param fid file descriptor to open and writable netCDF
  *            particle file.
  * @param rec record number to write.
  * @param t time at specified record.
  * @param np number of particles to write.
  * @param p array of particles to write.
  */
void pt_write(int fid, int rec, double t, long int np, particle_t *p)
{
  double *b;
  short *f;
  unsigned char *c;
  size_t start[2];
  size_t count[2];
  long n;

  /* Write time value */
  start[0] = rec;
  count[0] = 1;
  nc_put_vara_double(fid, ncw_var_id(fid, "t"), start, count, &t);

  start[1] = 0;
  count[1] = np;

  /* Allocate buffers */
  b = d_alloc_1d(np);
  f = s_alloc_1d(np);
  c = uc_alloc_1d(np);

  /* Transfer and write data */
  for (n = 0; n < np; n++)
    b[n] = p[n].e1;
  nc_put_vara_double(fid, ncw_var_id(fid, "x"), start, count, b);

  for (n = 0; n < np; n++)
    b[n] = p[n].e2;
  nc_put_vara_double(fid, ncw_var_id(fid, "y"), start, count, b);

  for (n = 0; n < np; n++)
    b[n] = p[n].e3;
  nc_put_vara_double(fid, ncw_var_id(fid, "z"), start, count, b);

  for (n = 0; n < np; n++)
    f[n] = p[n].flag;
  nc_put_vara_short(fid, ncw_var_id(fid, "flag"), start, count, f);

  for (n = 0; n < np; n++)
    if (p[n].dumpf & PT_AGE)
      c[n] = p[n].out_age;
  nc_put_vara_uchar(fid, ncw_var_id(fid, "age"), start, count, c);

  for (n = 0; n < np; n++)
    if (p[n].dumpf & PT_SIZE)
      c[n] = p[n].out_size;
  nc_put_vara_uchar(fid, ncw_var_id(fid, "size"), start, count, c);

  d_free_1d(b);
  s_free_1d(f);
  free_1d(c);

  nc_sync(fid);
}
