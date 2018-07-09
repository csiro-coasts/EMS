/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/io/ncw.c
 *
 *  \brief netCDF wrappers
 *
 * Simple wrappers to netCDF library procedures for better error messaging.
 * Some straightforward extensions to netcdf library procedures.
 *
 * NOTE: these functions will quit the application on error,
 *       be aware that if an emergency dump should be allowed to write 
 *       inconsistent data, these functions should not be used! 
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: ncw.c 5833 2018-06-27 00:21:35Z riz008 $
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <errno.h>
#include "netcdf.h"
#include "ems.h"

const char ncw_version[] = "0.11";

/* This macro is substituted in error messages instead of the name of a
 * variable in cases when the name could not be found by the variable id.
 */
#define STR_UNKNOWN "???"

/* Used in ncw_find_vars()
 */
#define NALLOCATED_START 10

/** Locate a netCDF dimension name and return the dimension id.
  * A replacement function to mimic the NetCDF 2.X function ncdimid.
  * @param fid file descriptor of an open netcdf file.
  * @param name netcdf dimension name.
  * @return dimension identifier.
  */
int ncw_dim_id(int fid, const char *name)
{
  int id = -1;

  if (nc_inq_dimid(fid, name, &id) != NC_NOERR)
    id = -1;

  return id;
}


/** Locate a netCDF variable name and return the variable id.
  * A replacement function to mimic the NetCDF 2.X function ncvarid.
  * @param fid file descriptor of an open netcdf file.
  * @param name netcdf variable name.
  * @return variable identifier.
  */
int ncw_var_id(int fid, const char *name)
{
  int id = -1;

  if (nc_inq_varid(fid, name, &id) != NC_NOERR)
    id = -1;

  return id;
}


/** Finds which variables in a netCDF file have a certain
  * set of dimensions and a particular attribute with a particular
  * value. Arguments are as follows:
  * 
  * @param fid file descriptor of an open netcdf file.
  * @param nvdims number of dimensions for matching variables. Only
  *               checked if > 0.
  * @param vdims list of dimension ids for matching variables. Only
  *              checked if not NULL.
  * @param attr name of attribute. Only checked if not NULL.
  * @param attval value of attribute (must be a string).
  *               Only checked if not NULL.
  * @param list Output list of variable ids which have the above
  * 	        properties
  * @return number of matching variables found.
  */
int ncw_var_find(int fid, int nvdims, int *vdims, const char *attr,
             const char *attval, int *list)
{
  int ndims;
  int nvars;
  int ngatts;
  int recdim;
  int id;
  int i;
  int n;

  /* Find out file properties */
  nc_inq(fid, &ndims, &nvars, &ngatts, &recdim);

  if (ndims < nvdims)
    return (0);

  /* Loop over the variables */
  n = 0;
  for (id = 0; id < nvars; id++) {
    char name[MAXLINELEN];
    nc_type dt;
    int nd;
    int d[NC_MAX_VAR_DIMS];
    int na;
    char aval[MAXLINELEN];
    int attnum;
    int match;

    match = 1;
    memset(aval, 0, MAXLINELEN);
    /* Check number of dimensions */
    nc_inq_var(fid, id, name, &dt, &nd, d, &na);
    if ((nvdims > 0) && (nd != nvdims))
      match = 0;

    /* Check each dimension */
    for (i = 0; vdims != NULL && match && i < nd; i++)
      if ((vdims[i] >= 0) && (d[i] != vdims[i]))
        match = 0;

    /* Check attribute existence */
    if (match && (attr != NULL) &&
        (nc_inq_attid(fid, id, attr, &attnum) != NC_NOERR))
      match = 0;

    /* Check attribute value */
    if (match && (attr != NULL) && (attval != NULL)) {
      if (nc_get_att_text(fid, id, attr, aval) == NC_NOERR) {
        if (strcmp(aval, attval) != 0)
          match = 0;
      } else
        match = 0;
    }

    /* Store id if we have a match */
    if (match) {
      list[n] = id;
      n++;
    }
  }

  return (n);
}


/** Determine the number of bytes needed to store a
  * single value of a specified netCDF variable.
  * @param fid file descriptor of an open netcdf file.
  * @param vid id of the variable.
  * @return number of bytes.
  */
int ncw_var_size(int fid, int vid)
{
  nc_type datatype;
  size_t size = 0;

  nc_inq_vartype(fid, vid, &datatype);
  switch (datatype) {
  case NC_BYTE:
  case NC_CHAR:
    size = sizeof(char);
    break;

  case NC_SHORT:
    size = sizeof(short);
    break;

  case NC_INT:
    size = sizeof(int);
    break;

  case NC_FLOAT:
    size = sizeof(float);
    break;

  case NC_DOUBLE:
    size = sizeof(double);
    break;

  default:
    quit("ncw_var_size: Unsupported netCDF data type.");
    break;
  }

  return ((int)size);
}

/** Reads a hyperslab of a single netcdf variable, checking that
  * the variable type matches a specified size.
  *
  * This routine mimics the functionality available below netCDF 3.
  * It is best to avoid this function in any new code.
  *
  * @param fid file descriptor of an open netcdf file.
  * @param name netcdf variable name.
  * @param size variable size.
  * @param start pointer to array of hyperslab start positions. 
  * @param count pointer to array of hyperslab sizes.
  * @param buf pointer to the array in which the values will be read.
  */
void ncw_var_read(int fid, char *name, int size, long int *start,
                  long int *count, void *buf)
{
  int n;
  nc_type datatype;
  int vid = ncw_var_id(fid, name);

  if ((n = ncw_var_size(fid, vid)) != size)
    quit("ncw_var_read: %s has %d bytes per value in file, code expects %d\n",
       n, size);

  nc_inq_vartype(fid, vid, &datatype);
  switch (datatype) {
  case NC_BYTE:
  case NC_CHAR:
    nc_get_vara_schar(fid, vid, (size_t *) start, (size_t *) count,
                      (signed char *)buf);
    break;

  case NC_SHORT:
    nc_get_vara_short(fid, vid, (size_t *) start, (size_t *) count,
                      (short *)buf);
    break;

  case NC_INT:
    nc_get_vara_int(fid, vid, (size_t *) start, (size_t *) count,
                    (int *)buf);
    break;

  case NC_FLOAT:
    nc_get_vara_float(fid, vid, (size_t *) start, (size_t *) count,
                      (float *)buf);
    break;

  case NC_DOUBLE:
    nc_get_vara_double(fid, vid, (size_t *) start, (size_t *) count,
                       (double *)buf);
    break;

  default:
    quit("ncw_var_read: Unsupported netCDF data type.");
    break;
  }
}


static void _ncw_inq_varname(const char fname[], int ncid, int varid, char varname[])
{
    int status;

    if (varid == -1)
        strcpy(varname, "NC_GLOBAL");
    else if ((status = nc_inq_varname(ncid, varid, varname)) != NC_NOERR)
        quit("\"%s\": nc_inq_varname(): failed for varid = %d: %s\n", (fname != NULL?fname:EMPTY), varid, nc_strerror(status));
}

static void _ncw_inq_varid(const char fname[], int ncid, const char varname[], int* varid)
{
    int status;

    if (strcmp(varname, "NC_GLOBAL") == 0)
        *varid = -1;
    else if ((status = nc_inq_varid(ncid, varname, varid)) != NC_NOERR)
        quit("\"%s\": nc_inq_varid(): failed for varid = %d: %s\n", (fname != NULL?fname:EMPTY), varid, nc_strerror(status));
}

#define STRBUFSIZE 1024

/* Prints array of integers to a string. E.g., {1,2,5} will be printed as
 * "(1,2,5)".
 */
static char* int2str(int n, const int v[])
{
    int i;
    char all[STRBUFSIZE] = "(";
    char next[STRBUFSIZE];

    for (i = 0; i < n; ++i) {
        if (i < n - 1)
            sprintf(next, "%d,", v[i]);
        else
            sprintf(next, "%d", v[i]);
        assert(strlen(all) + strlen(next) + 2 < STRBUFSIZE);
        strcat(all, next);
    }
    assert(strlen(all) + 2 < STRBUFSIZE);
    strcat(all, ")");

    return strdup(all);
}


/*UR ADDED since some platforms implement size_t as long
 * 
 */
static char* usize_t2str(int n, const size_t  v[])
{
       int i;
    char all[STRBUFSIZE] = "(";
    char next[STRBUFSIZE];

    for (i = 0; i < n; ++i) {
        if (i < n - 1)
            sprintf(next, "%zi,", v[i]);
        else
            sprintf(next, "%zi", v[i]);
        assert(strlen(all) + strlen(next) + 2 < STRBUFSIZE);
        strcat(all, next);
    }
    assert(strlen(all) + 2 < STRBUFSIZE);
    strcat(all, ")");

    return strdup(all);
}


static char* uint2str(int n, const unsigned int v[])
{
    int i;
    char all[STRBUFSIZE] = "(";
    char next[STRBUFSIZE];

    for (i = 0; i < n; ++i) {
        if (i < n - 1)
            sprintf(next, "%d,", v[i]);
        else
            sprintf(next, "%d", v[i]);
        assert(strlen(all) + strlen(next) + 2 < STRBUFSIZE);
        strcat(all, next);
    }
    assert(strlen(all) + 2 < STRBUFSIZE);
    strcat(all, ")");

    return strdup(all);
}

static char* double2str(int n, const double v[])
{
    int i;
    char all[STRBUFSIZE] = "(";
    char next[STRBUFSIZE];

    for (i = 0; i < n && i < 3; ++i) {
        if (i < n - 1)
            sprintf(next, "%.4g,", v[i]);
        else
            sprintf(next, "%.4g", v[i]);
        assert(strlen(all) + strlen(next) + 2 < STRBUFSIZE);
        strcat(all, next);
    }
    if (n > 3) {
        assert(strlen(all) + strlen(",...") + 1 < STRBUFSIZE);
        strcat(all, ",...");
    }
    assert(strlen(all) + strlen(")") + 1 < STRBUFSIZE);
    strcat(all, ")");

    return strdup(all);
}

static char* float2str(int n, const float v[])
{
    int i;
    char all[STRBUFSIZE] = "(";
    char next[STRBUFSIZE];

    for (i = 0; i < n && i < 3; ++i) {
        if (i < n - 1)
            sprintf(next, "%.4g,", v[i]);
        else
            sprintf(next, "%.4g", v[i]);
        assert(strlen(all) + strlen(next) + 2 < STRBUFSIZE);
        strcat(all, next);
    }
    if (n > 3) {
        assert(strlen(all) + strlen(",...") + 1 < STRBUFSIZE);
        strcat(all, ",...");
    }
    assert(strlen(all) + strlen(")") + 1 < STRBUFSIZE);
    strcat(all, ")");

    return strdup(all);
}

/* *INDENT-OFF* */
static struct nctype2str {
    nc_type type;
    char* str;
} nctypes2str[] = {
    {-1, "UNKNOWN"},
    {NC_BYTE, "NC_BYTE"},
    {NC_CHAR, "NC_CHAR"},
    {NC_SHORT, "NC_SHORT"},
    {NC_INT, "NC_INT"},
    {NC_FLOAT, "NC_FLOAT"},
    {NC_DOUBLE, "NC_DOUBLE"},
};
/* *INDENT-ON* */

int ncw_create(const char fname[], int mode, int* ncid)
{
  int status = nc_create((fname != NULL?fname:EMPTY), mode, ncid);

  if (status != NC_NOERR)
    quit("\"%s\": nc_create(): failed: %s\n", (fname != NULL?fname:EMPTY), nc_strerror(status));
  return status;
}

void ncw_open(const char fname[], int mode, int* ncid)
{
    int status = nc_open((fname != NULL?fname:EMPTY), mode, ncid);

    if (status != NC_NOERR)
        quit("\"%s\": nc_open(): failed: %s\n", (fname != NULL?fname:EMPTY), nc_strerror(status));
}

void ncw_close(const char fname[], int ncid)
{
    int status = nc_close(ncid);

    if (status != NC_NOERR)
        quit("\"%s\": nc_close(): failed: %s\n", (fname != NULL?fname:EMPTY), nc_strerror(status));
}


void ncw_enddef(const char fname[], int ncid)
{
  int status = nc_enddef(ncid);
  if (status != NC_NOERR)
    quit("\"%s\": nc_enddef(): failed: %s\n", 
	                   (fname != NULL?fname:EMPTY), nc_strerror(status));
}


void ncw_sync(const char fname[], int ncid)
{
    int status;

    status = nc_sync(ncid);
    if (status != NC_NOERR)
      quit("\"%s\": ncsync(): failed: %s\n", (fname != NULL?fname:EMPTY), nc_strerror(status));
}

void ncw_copy_att(const char fname_src[], int ncid_src, int varid_src, const char attname[], const char fname_dst[], int ncid_dst, int varid_dst)
{
    int status = nc_copy_att(ncid_src, varid_src, attname, ncid_dst, varid_dst);

    if (status != NC_NOERR) {
        char varname_src[NC_MAX_NAME] = STR_UNKNOWN;
        char varname_dst[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname_src, ncid_src, varid_src, varname_src);
        _ncw_inq_varname(fname_dst, ncid_dst, varid_dst, varname_dst);

        quit("\"%s\": -> %s:  nc_copy_att(): failed for varid_src = %d (varname = \"%s\"), attname = \"%s\", varid_dst = %d, (varname = \"%s\"): %s\n", fname_src, fname_dst, varid_src, varname_src, attname, varid_dst, varname_dst, nc_strerror(status));
    }
}

void ncw_del_att(const char fname[], int ncid, int varid, const char attname[])
{
    int status = nc_del_att(ncid, varid, attname);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_del_att(): failed for varid = %d (varname = \"%s\"): %s\n", (fname != NULL?fname:EMPTY), varid, varname, nc_strerror(status));
    }
}

void ncw_inq_attname(const char fname[], int ncid, int varid, int attrid, char attname[])
{
    int status = nc_inq_attname(ncid, varid, attrid, attname);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_inq_attname(): failed for varid = %d (varname = \"%s\"): %s\n", (fname != NULL?fname:EMPTY), varid, varname, nc_strerror(status));
    }
}

void ncw_inq_attlen(const char fname[], int ncid, int varid, const char attname[], size_t* len)
{
    int status = nc_inq_attlen(ncid, varid, attname, len);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_inq_attlen(): failed for varid = %d (varname = \"%s\"), attname = \"%s\": %s\n", (fname != NULL?fname:EMPTY), varid, varname, attname, nc_strerror(status));
    }
}

void ncw_inq_ndims(const char fname[], int ncid, int* ndims)
{
    int status = nc_inq_ndims(ncid, ndims);

    if (status != NC_NOERR)
        quit("\"%s\": nc_inq_ndims(): failed: %s\n", (fname != NULL?fname:EMPTY), nc_strerror(status));
}

void ncw_inq_dim(const char fname[], int ncid, int dimid, char dimname[], size_t* len)
{
    int status = nc_inq_dim(ncid, dimid, dimname, len);

    if (status != NC_NOERR)
        quit("\"%s\": nc_inq_dim(): failed for dimid = %d: %s\n", (fname != NULL?fname:EMPTY), dimid, nc_strerror(status));
}

void ncw_inq_dimid(const char fname[], int ncid, const char dimname[], int* dimid)
{
    int status = nc_inq_dimid(ncid, dimname, dimid);

    if (status != NC_NOERR)
        quit("\"%s\": nc_inq_dimid(): failed for dimname = \"%s\": %s\n", (fname != NULL?fname:EMPTY), dimname, nc_strerror(status));
}

void ncw_inq_dimname(const char fname[], int ncid, int dimid, char dimname[])
{
    int status = nc_inq_dimname(ncid, dimid, dimname);

    if (status != NC_NOERR)
        quit("\"%s\": nc_inq_dimname(): failed for dimid = %d: %s\n", (fname != NULL?fname:EMPTY), dimid, nc_strerror(status));
}

void ncw_inq_dimlen(const char fname[], int ncid, int dimid, size_t* len)
{
    int status = nc_inq_dimlen(ncid, dimid, len);

    if (status != NC_NOERR) {
        char dimname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_dimname(fname, ncid, dimid, dimname);
        quit("\"%s\": nc_inq_dimlen(): failed for dimid = %d (dimname = \"%s\"): %s\n", (fname != NULL?fname:EMPTY), dimid, dimname, nc_strerror(status));
    }
}

void ncw_inq_nvars(const char fname[], int ncid, int* nvars)
{
    int status = nc_inq_nvars(ncid, nvars);

    if (status != NC_NOERR)
        quit("\"%s\": nc_inq_nvars(): failed: %s\n", (fname != NULL?fname:EMPTY), nc_strerror(status));
}

void ncw_inq_var(const char fname[], int ncid, int varid, char varname[], nc_type* xtype, int* ndims, int dimids[], int* natts)
{
    int status = nc_inq_var(ncid, varid, varname, xtype, ndims, dimids, natts);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_inq_var(): failed for varid = %d (varname = \"%s\"): %s\n", (fname != NULL?fname:EMPTY), varid, varname, nc_strerror(status));
    }
}

void ncw_inq_varid(const char* fname, int ncid, const char varname[], int* varid)
{
    int status = nc_inq_varid(ncid, varname, varid);

    if (status != NC_NOERR)
        quit("\"%s\": nc_inq_varid(): failed for varname = \"%s\": %s\n", (fname != NULL?fname:EMPTY), varname, nc_strerror(status));
}

void ncw_inq_varname(const char fname[], int ncid, int varid, char varname[])
{
    int status = nc_inq_varname(ncid, varid, varname);

    if (status != NC_NOERR)
        quit("\"%s\": nc_inq_varname(): failed for varid = %d: %s\n", (fname != NULL?fname:EMPTY), varid, nc_strerror(status));
}

void ncw_inq_vartype(const char fname[], int ncid, int varid, nc_type* xtype)
{
    int status = nc_inq_vartype(ncid, varid, xtype);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_inq_vartype(): failed for varid = %d (varname = %s): %s\n", (fname != NULL?fname:EMPTY), varid, varname, nc_strerror(status));
    }
}

void ncw_get_var_double(const char fname[], int ncid, int varid, double v[])
{
    int status = nc_get_var_double(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_get_var_double(): failed for varid = %d (varname = \"%s\"): %s\n", (fname != NULL?fname:EMPTY), varid, varname, nc_strerror(status));
    }
}

void ncw_get_vara_double(const char fname[], int ncid, int varid, const size_t start[], const size_t count[], double v[])
{
    int status = nc_get_vara_double(ncid, varid, start, count, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_get_vara_double(): failed for varid = %d (varname = \"%s\"): %s\n", (fname != NULL?fname:EMPTY), varid, varname, nc_strerror(status));
    }
}

void ncw_get_var1_double(const char fname[], int ncid, int varid, const size_t len[], double* in)
{
    int status = nc_get_var1_double(ncid, varid, len, in);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_get_var1_double(): failed for varid = %d (varname = \"%s\"): %s\n", (fname != NULL?fname:EMPTY), varid, varname, nc_strerror(status));
    }
}

void ncw_get_var_int(const char fname[], int ncid, int varid, int v[])
{
    int status = nc_get_var_int(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_get_var_int(): failed for varid = %d (varname = \"%s\"): %s\n", (fname != NULL?fname:EMPTY), varid, varname, nc_strerror(status));
    }
}

void ncw_get_vara_int(const char fname[], int ncid, int varid, const size_t start[], const size_t count[], int v[])
{
    int status = nc_get_vara_int(ncid, varid, start, count, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_get_vara_int(): failed for varid = %d (varname = \"%s\"): %s\n", (fname != NULL?fname:EMPTY), varid, varname, nc_strerror(status));
    }
}

void ncw_get_att_text(const char fname[], int ncid, int varid, const char attname[], char v[])
{
    nc_type xtype;
    size_t len;
    int status;

    ncw_inq_att(fname, ncid, varid, attname, &xtype, &len);

    if (xtype != NC_CHAR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname((fname != NULL?fname:EMPTY), ncid, varid, varname);
        quit("ncw_get_att_text(): failed for varid = %d (varname = \"%s\"), attname = \"%s\": the attribute is of %s type\n", (fname != NULL?fname:EMPTY), varid, varname, attname, ncw_nctype2str(xtype));
    }

    status = nc_get_att_text(ncid, varid, attname, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_get_att_text(): failed for varid = %d (varname = \"%s\"), attname = \"%s\": %s\n", (fname != NULL?fname:EMPTY), varid, varname, attname, nc_strerror(status));
    }
}

void ncw_get_att_int(const char fname[], int ncid, int varid, const char attname[], int v[])
{
    nc_type xtype;
    size_t len;
    int status;

    ncw_inq_att(fname, ncid, varid, attname, &xtype, &len);

    if (xtype != NC_INT && xtype != NC_BYTE && xtype != NC_SHORT) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("ncw_get_att_int(): failed for varid = %d (varname = \"%s\"), attname = \"%s\": the attribute is of %s type\n", (fname != NULL?fname:EMPTY), varid, varname, attname, ncw_nctype2str(xtype));
    }

    status = nc_get_att_int(ncid, varid, attname, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_get_att_int(): failed for varid = %d (varname = \"%s\"), attname = \"%s\": %s\n", (fname != NULL?fname:EMPTY), varid, varname, attname, nc_strerror(status));
    }
}

void ncw_get_att_double(const char fname[], int ncid, int varid, const char attname[], double v[])
{
    nc_type xtype;
    size_t len;
    int status;

    ncw_inq_att(fname, ncid, varid, attname, &xtype, &len);

    if (xtype != NC_DOUBLE && xtype != NC_FLOAT) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("ncw_get_att_double(): failed for varid = %d (varname = \"%s\"), attname = \"%s\": the attribute is of %s type\n", (fname != NULL?fname:EMPTY), varid, varname, attname, ncw_nctype2str(xtype));
    }

    status = nc_get_att_double(ncid, varid, attname, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": ncw_get_att_double(): failed for varid = %d (varname = \"%s\"), attname = \"%s\": %s\n", (fname != NULL?fname:EMPTY), varid, varname, attname, nc_strerror(status));
    }
}

void ncw_def_dim(const char fname[], int ncid, const char dimname[], size_t len, int* dimid)
{
    int status = nc_def_dim(ncid, dimname, len, dimid);

    if (status != NC_NOERR)
        quit("\"%s\": nc_def_dim(): failed for dimname = \"%s\", dimlen = %d: %s\n", (fname != NULL?fname:EMPTY), dimname, len, nc_strerror(status));
}


void ncw_def_var(const char fname[], int ncid, const char varname[], nc_type xtype, int ndims, const int dimids[], int* varid)
{
    int status = nc_def_var(ncid, varname, xtype, ndims, dimids, varid);

    if (status != NC_NOERR)
        quit("\"%s\": nc_def_var(): failed for varname = \"%s\", vartype = %s, ndims = %d, dimids = %s: %s\n", (fname != NULL?fname:EMPTY), varname, ncw_nctype2str(xtype), ndims, uint2str(ndims, (const unsigned int*) dimids), nc_strerror(status));
}

/* Same as above except with an extra option */
void ncw_def_var2(const char fname[], int ncid, const char varname[], nc_type xtype, int ndims, const int dimids[], int* varid, int compress)
{
  int status;
  
  /* Call above function first */
  ncw_def_var(fname, ncid, varname, xtype, ndims, dimids, varid);
  
  /*
   * Shuffle shuffles the byte order of the data in the hope of
   * gathering more higher order bits of similar type to boost
   * compression
   */
#ifdef NC_NETCDF4 /* The ifdef is just for compilation purposes */
  if (compress) {
    /* I got a deflate value of 2 from some research on the net (FR) */
    status = nc_def_var_deflate(ncid, varid[0], NC_SHUFFLE, 1, 2);
    if (status != NC_NOERR)
      quit("nc_def_var_deflate error : %s\n", nc_strerror(status));
  }
#endif
}

void ncw_put_att_text(const char fname[], int ncid, int varid, const char attname[], const char v[])
{
    int status = nc_put_att_text(ncid, varid, attname, strlen(v), v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_put_att_text(): failed for varid = %d (varname = \"%s\"), attname = \"%s\", attvalue = \"%s\": %s\n", (fname != NULL?fname:EMPTY), varid, varname, attname, v, nc_strerror(status));
    }
}

void ncw_put_att_double(const char fname[], int ncid, int varid, const char attname[], size_t len, const double v[])
{
    int status = nc_put_att_double(ncid, varid, attname, NC_DOUBLE, len, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_put_att_double(): failed for varid = %d (varname = \"%s\"), attname = \"%s\", attvalue(s) = \"%s\": %s\n", (fname != NULL?fname:EMPTY), varid, varname, attname, double2str(len, v), nc_strerror(status));
    }

    if (strcmp(attname, "_FillValue") == 0 && varid != NC_GLOBAL) {
        nc_type xtype;

        ncw_inq_vartype(fname, ncid, varid, &xtype);
        if (xtype != NC_DOUBLE) {
            char varname[NC_MAX_NAME] = "STR_UNKNOWN";

            ncw_inq_varname(fname, ncid, varid, varname);

            quit("\"%s\": ncw_put_att_float(): fatal inconsistency for varid = %d (varname = \"%s\"): attype = NC_DOUBLE differs from vartype = %s for attname = \"_FillValue\"\n", (fname != NULL?fname:EMPTY), varid, varname, ncw_nctype2str(xtype));
        }
    }
}

void ncw_put_att_float(const char fname[], int ncid, int varid, const char attname[], size_t len, const float v[])
{
    int status = nc_put_att_float(ncid, varid, attname, NC_FLOAT, len, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = "STR_UNKNOWN";

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_put_att_float(): failed for varid = %d (varname = \"%s\"), attname = \"%s\", attvalue = \"%s\": %s\n", (fname != NULL?fname:EMPTY), varid, varname, attname, float2str(len, v), nc_strerror(status));
    }

    if (strcmp(attname, "_FillValue") == 0 && varid != NC_GLOBAL) {
        nc_type xtype;

        ncw_inq_vartype(fname, ncid, varid, &xtype);
        if (xtype != NC_FLOAT) {
            char varname[NC_MAX_NAME] = "STR_UNKNOWN";

            ncw_inq_varname(fname, ncid, varid, varname);

            quit("\"%s\": ncw_put_att_float(): failed for varid = %d (varname = \"%s\"): attype = NC_FLOAT differs from vartype = %s for attname = \"_FillValue\"\n", (fname != NULL?fname:EMPTY), varid, varname, ncw_nctype2str(xtype));
        }
    }
}

void ncw_put_att_int(const char fname[], int ncid, int varid, const char attname[], size_t len, const int v[])
{
    int status = nc_put_att_int(ncid, varid, attname, NC_INT, len, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_put_att_int(): failed for varid = %d (varname = \"%s\"), attname = \"%s\", attvalue(s) = \"%s\": %s\n", (fname != NULL?fname:EMPTY), varid, varname, attname, int2str(len, v), nc_strerror(status));
    }
}

void ncw_put_var_double(const char fname[], int ncid, int varid, const double v[])
{
    int status = nc_put_var_double(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_put_var_double(): failed for varid = %d (varname = \"%s\"): %s\n", (fname != NULL?fname:EMPTY), varid, varname, nc_strerror(status));
    }
}

void ncw_put_vara_double(const char fname[], int ncid, int varid, const size_t start[], const size_t count[], const double v[])
{
    int status = nc_put_vara_double(ncid, varid, start, count, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;
        int ndims = 0;

        ncw_inq_varname(fname, ncid, varid, varname);
        ncw_inq_varndims(fname, ncid, varid, &ndims);
        quit("\"%s\": nc_put_vara_double(): failed for varid = %d (varname = \"%s\"), start = %s, count = %s: %s\n", (fname != NULL?fname:EMPTY), varid, varname, usize_t2str(ndims, start), usize_t2str(ndims, count), nc_strerror(status));
    }
}

void ncw_put_vara_int(const char fname[], int ncid, int varid, const size_t start[], const size_t count[], const int v[])
{
    int status = nc_put_vara_int(ncid, varid, start, count, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;
        int ndims = 0;

        ncw_inq_varname(fname, ncid, varid, varname);
        ncw_inq_varndims(fname, ncid, varid, &ndims);
        quit("\"%s\": ncw_put_vara_int(): failed for varid = %d (varname = \"%s\"), start = %s, count = %s: %s\n", (fname != NULL?fname:EMPTY), varid, varname, usize_t2str(ndims, start), usize_t2str(ndims, count), nc_strerror(status));
    }
}


void ncw_put_vara_long(const char fname[], int ncid, int varid, const size_t start[], const size_t count[], const long v[])
{
    int status = nc_put_vara_long(ncid, varid, start, count, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;
        int ndims = 0;

        ncw_inq_varname(fname, ncid, varid, varname);
        ncw_inq_varndims(fname, ncid, varid, &ndims);
        quit("\"%s\": ncw_put_vara_int(): failed for varid = %d (varname = \"%s\"), start = %s, count = %s: %s\n", (fname != NULL?fname:EMPTY), varid, varname, usize_t2str(ndims, start), usize_t2str(ndims, count), nc_strerror(status));
    }
}


void ncw_inq_vardimid(const char fname[], int ncid, int varid, int dimids[])
{
    int status = nc_inq_vardimid(ncid, varid, dimids);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": ncw_inq_vardimid(): failed for varid = %d (varname = \"%s\"): %s\n", (fname != NULL?fname:EMPTY), varid, varname, nc_strerror(status));
    }
}


void ncw_inq_varndims(const char fname[], int ncid, int varid, int* ndims)
{
    int status = nc_inq_varndims(ncid, varid, ndims);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = "STR_UNKNOWN";

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_inq_varndims(): failed for varid = %d (varname = \"%s\"): %s\n", (fname != NULL?fname:EMPTY), varid, varname, nc_strerror(status));
    }
}

void ncw_inq_varnatts(const char fname[], int ncid, int varid, int* natts)
{
    int status = nc_inq_varnatts(ncid, varid, natts);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = "STR_UNKNOWN";

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_inq_varnatts(): failed for varid = %d (varnatts = \"%s\"): %s\n", (fname != NULL?fname:EMPTY), varid, varname, nc_strerror(status));
    }
}

void ncw_rename_dim(const char fname[], int ncid, const char oldname[], const char newname[])
{
    int dimid;
    int status;

    ncw_inq_dimid(fname, ncid, oldname, &dimid);
    status = nc_rename_dim(ncid, dimid, newname);

    if (status != NC_NOERR)
        quit("\"%s\": nc_rename_dim(): failed for dimid = %d (oldname = \"%s\", newname = \"%s\"): %s\n", dimid, oldname, newname, nc_strerror(status));
}

void ncw_rename_var(const char fname[], int ncid, const char oldname[], const char newname[])
{
    int varid;
    int status;

    ncw_inq_varid(fname, ncid, oldname, &varid);
    status = nc_rename_var(ncid, varid, newname);

    if (status != NC_NOERR)
        quit("\"%s\": nc_rename_var(): failed for varid = %d (oldname = \"%s\", newname = \"%s\"): %s\n", varid, oldname, newname, nc_strerror(status));
}

void ncw_rename_att(const char fname[], int ncid, const char varname[], const char oldname[], const char newname[])
{
    int varid;
    int status;

    _ncw_inq_varid(fname, ncid, varname, &varid);
    status = nc_rename_att(ncid, varid, oldname, newname);

    if (status != NC_NOERR)
        quit("\"%s\": nc_rename_att(): failed for varid = %d, oldname = \"%s\", newname = \"%s\": %s\n", varid, oldname, newname, nc_strerror(status));
}

/*
 * The following procedures do not have direct analogues in the netcdf library.
 */

const char* ncw_nctype2str(nc_type type)
{
    int i;

    for (i = 1; i < sizeof(nctypes2str) / sizeof(struct nctype2str); ++i) {
        if (type == nctypes2str[i].type)
            return nctypes2str[i].str;
    }

    return nctypes2str[0].str;
}

size_t ncw_sizeof(nc_type type)
{
    switch (type) {
    case NC_BYTE:
    case NC_CHAR:
        return sizeof(char);
    case NC_SHORT:
        return sizeof(short);
    case NC_INT:
        return sizeof(int);
    case NC_FLOAT:
        return sizeof(float);
    case NC_DOUBLE:
        return sizeof(double);
    default:
        quit("ncw_sizeof(): unknown type\n");
    }

    return UINT_MAX;
}

/** Copies all dimensions from one NetCDF file to another.
 *
 * @param fname_src Source file name
 * @param ncid_src Source file id
 * @param fname_dst Destination file name
 * @param ncid_dst Destination file id
 */
void ncw_copy_dims(const char* fname_src, int ncid_src, const char* fname_dst, int ncid_dst)
{
    int ndims;
    int unlimdimid = -1;
    int i;

    ncw_inq_ndims(fname_src, ncid_src, &ndims);
    nc_inq_unlimdim(ncid_src, &unlimdimid);

    for (i = 0; i < ndims; ++i) {
        char dimname[NC_MAX_NAME] = STR_UNKNOWN;
        size_t size;
        int dimid;

        memset(dimname, 0, NC_MAX_NAME);

        ncw_inq_dim(fname_src, ncid_src, i, dimname, &size);
        if (i == unlimdimid)
            size = NC_UNLIMITED;
        ncw_def_dim(fname_dst, ncid_dst, dimname, size, &dimid);
    }
}

/** Copies a specified variable from one NetCDF file to another.
 *
 * @param fname_src Source file name
 * @param ncid_src Source file id
 * @param vid_src Variable id
 * @param fname_dst Destination file name
 * @param ncid_dst Destination file id
 */
void ncw_copy_var(const char* fname_src, int ncid_src, int vid_src, const char* fname_dst, int ncid_dst)
{
    char varname[NC_MAX_NAME] = STR_UNKNOWN;
    int vid_dst;
    nc_type type;
    int ndims;
    int dimids_src[NC_MAX_DIMS], dimids_dst[NC_MAX_DIMS];
    int natts;
    int i;

    ncw_inq_varname(fname_src, ncid_src, vid_src, varname);
    ncw_inq_var(fname_src, ncid_src, vid_src, NULL, &type, &ndims, dimids_src, &natts);

    for (i = 0; i < ndims; ++i) {
        char dimname[NC_MAX_NAME];
        size_t len;

        ncw_inq_dim(fname_src, ncid_src, dimids_src[i], dimname, &len);
        ncw_inq_dimid(fname_dst, ncid_dst, dimname, &dimids_dst[i]);
    }

    ncw_def_var(fname_dst, ncid_dst, varname, type, ndims, dimids_dst, &vid_dst);

    ncw_copy_atts(fname_src, ncid_src, vid_src, fname_dst, ncid_dst, vid_dst);
}

void ncw_inq(const char fname[], int ncid, int* ndims, int* nvars, int* natts, int* unlimdimid)
{
    int status = nc_inq(ncid, ndims, nvars, natts, unlimdimid);

    if (status != NC_NOERR)
        quit("\"%s\": nc_inq(): failed: %s\n", (fname != NULL?fname:EMPTY), nc_strerror(status));
}

void ncw_inq_unlimdimid(const char fname[], int ncid, int* unlimdimid)
{
    int status = nc_inq_unlimdim(ncid, unlimdimid);

    if (status != NC_NOERR)
        quit("\"%s\": nc_inq_unlimdim()(): failed: %s\n", (fname != NULL?fname:EMPTY), nc_strerror(status));
}

int ncw_inq_nrecords(const char fname[], int ncid)
{
    int unlimdimid;
    size_t nrecords;

    ncw_inq_unlimdimid(fname, ncid, &unlimdimid);
    ncw_inq_dimlen(fname, ncid, unlimdimid, &nrecords);

    return nrecords;
}

void ncw_inq_att(const char fname[], int ncid, int varid, const char attname[], nc_type* xtype, size_t* len)
{
    int status = nc_inq_att(ncid, varid, attname, xtype, len);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_inq_att(): failed for varid = %d (varname = \"%s\"), attname = \"%s\", atttype = %d: %s\n", (fname != NULL?fname:EMPTY), varid, varname, attname, ncw_nctype2str(*xtype), nc_strerror(status));
    }
}

/** Gets the id for the first dimension found to be present in a NetCDF file
 * out of two dimensions specified by names. Useful for handling both new and
 * old data formats.
 *
 * @param fname NetCDF file name
 * @param ncid NetCDF file id
 * @param dimname1 The first dimension name to be tried
 * @param dimname2 The second dimension name to be tried
 * @param dimid Dimension id (output)
 */
void ncw_inq_dimid2(const char fname[], int ncid, const char dimname1[], const char dimname2[], int* dimid)
{
    int status1 = nc_inq_dimid(ncid, dimname1, dimid);

    if (status1 != NC_NOERR) {
        int status2 = nc_inq_dimid(ncid, dimname2, dimid);

        if (status2 != NC_NOERR) {
            quit("\"%s\": nc_inq_dimid(): failed for dimname = \"%s\": %s, and for dimname = \"%s\": %s\n", (fname != NULL?fname:EMPTY), dimname1, nc_strerror(status1), dimname2, nc_strerror(status2));
        }
    }
}

/** Gets the value(s) (converted to int type) of the first attribute found to
 * be present in a NetCDF file out of two specified by attribute names. Useful
 * for handling both new and old data formats.
 *
 * @param fname NetCDF file name
 * @param ncid NetCDF file id
 * @param varid Name of the variable the attrubute belongs to (use NC_GLOBAL
 *              for global attributes)
 * @param attname1 The first attribute name to be tried
 * @param attname2 The second attribute name to be tried
 * @param v Attribute value(s) (output)
 */
void ncw_get_att_int2(const char fname[], int ncid, int varid, const char attname1[], const char attname2[], int v[])
{
    nc_type xtype;
    size_t len;
    int status1 = nc_get_att_int(ncid, varid, attname1, v);

    if (status1 != NC_NOERR) {
        int status2 = nc_get_att_int(ncid, varid, attname2, v);

        if (status2 != NC_NOERR) {
            char varname[NC_MAX_NAME] = STR_UNKNOWN;

            _ncw_inq_varname(fname, ncid, varid, varname);
            quit("\"%s\": nc_get_att_int(): failed for varid = %d (varname = \"%s\"), attname = \"%s\": %s, and attname = \"%s\": %s\n", (fname != NULL?fname:EMPTY), varid, varname, attname1, nc_strerror(status1), attname2, nc_strerror(status2));
        }
    }

    ncw_inq_att(fname, ncid, varid, (status1 == NC_NOERR) ? attname1 : attname2, &xtype, &len);

    if (xtype != NC_INT && xtype != NC_BYTE && xtype != NC_SHORT) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("ncw_get_att_int2(): failed for varid = %d (varname = \"%s\"), attname = \"%s\": the attribute is of %s type\n", (fname != NULL?fname:EMPTY), varid, varname, (status1 == NC_NOERR) ? attname1 : attname2, ncw_nctype2str(xtype));
    }
}

/** Finds all variables with specified dimensions and specified attribute.
 *
 * @param fname NetCDF file name
 * @param ncid NetCDF file id
 * @param ndims Number of dimensions (0 if any)
 * @param dims Array of dimension ids
 * @param attname Attribute name that is be set (NULL if any)
 * @param attval Attribute value (NULL if any)
 * @param nvars Return -- number of variables found
 * @param vids Return -- array of variable ids (needs to be deallocated when
 *             necessary)
 */
void ncw_find_vars(const char fname[], int ncid, int ndims, const int dims[], const char attname[], const void* attval, int* nvars, int** vids)
{
    int nv = -1;                /* total number of variables */
    int nallocated = 0;
    int vid;

    *nvars = 0;
    *vids = NULL;

    ncw_inq_nvars(fname, ncid, &nv);

    for (vid = 0; vid < nv; ++vid) {
        char varname[NC_MAX_NAME] = "STR_UNKNOWN";
        nc_type vtype = -1;
        int nd = -1;
        int dids[NC_MAX_DIMS];
        int natts = -1;
        int i;

        /*
         * check if the dimensions are right for the variable
         */
        if (ndims > 0) {
            ncw_inq_var(fname, ncid, vid, varname, &vtype, &nd, dids, &natts);

            if (ndims != nd)
                goto nextvar;   /* (not using "continue" to be consistent) */

            /*
             * check dimensions
             */
            for (i = 0; i < ndims; ++i) {
                if (dids[i] != dims[i])
                    goto nextvar;
            }
        }

        /*
         * (so far so good)
         */

        /*
         * check if the attribute is present for the variable
         */
        if (attname != NULL) {
            for (i = 0; i < natts; ++i) {
                char aname[NC_MAX_NAME] = STR_UNKNOWN;

                ncw_inq_attname(fname, ncid, vid, i, aname);

                if (strcmp(aname, attname) != 0)
                    continue;

                if (attval != NULL) {
                    size_t alen = UINT_MAX;
                    nc_type atype = -1;
                    char* aval = NULL;

                    ncw_inq_att(fname, ncid, vid, aname, &atype, &alen);
                    if (alen <= 0)
                        continue;
                    aval = calloc(alen, ncw_sizeof(atype));
                    if (aval == NULL)
                        quit("\"%s\": ncw_find_vars(): could not allocate memory for attribute = \"%s\", type = %s, length = %d, varid = %d (varname = \"%s\"): %s\n", (fname != NULL?fname:EMPTY), aname, ncw_nctype2str(atype), alen, vid, (vid == NC_GLOBAL) ? "NC_GLOBAL" : varname, strerror(errno));
                    ncw_get_att_text(fname, ncid, vid, attname, aval);
                    if (memcmp(attval, aval, alen * sizeof(atype)) == 0) {
                        free(aval);
                        break;
                    }
                    free(aval);
                } else
                    break;
            }
            if (i >= natts)
                goto nextvar;
        }

        /*
         * the variable matches the criteria; add the id to the id array
         */

        /*
         * expand alocated space for the id array if necessary
         */
        if (nallocated == 0) {
            nallocated = NALLOCATED_START;
            *vids = malloc(nallocated * sizeof(int));
        } else if (*nvars == nallocated) {
            nallocated *= 2;
            *vids = realloc(vids, nallocated * sizeof(int));
        }

        (*vids)[*nvars] = vid;
        (*nvars)++;

      nextvar:
        ;
    }
}

/** Checks if specified attribute exists in the file.
 *
 * @param ncid NetCDF file id
 * @param varid Variable id (NC_GLOBAL for a global attribute)
 * @param attname Attribute name
 * @return 1 if attribute exists, 0 if not
 */
int ncw_att_exists(int ncid, int varid, const char attname[])
{
    if (nc_inq_attid(ncid, varid, attname, NULL) == NC_NOERR)
        return 1;
    else
        return 0;
}

/** Checks if specified variable exists in the file.
 *
 * @param ncid NetCDF file id
 * @param varname Variable name
 * @return 1 if attribute exists, 0 if not
 */
int ncw_var_exists(int ncid, const char varname[])
{
    int varid;

    if (nc_inq_varid(ncid, varname, &varid) == NC_NOERR)
        return 1;
    else
        return 0;
}

/** Checks if specified dimension exists in the file.
 *
 * @param ncid NetCDF file id
 * @param dimname Dimension name
 * @return 1 if attribute exists, 0 if not
 */
int ncw_dim_exists(int ncid, const char dimname[])
{
    int dimid;

    if (nc_inq_dimid(ncid, dimname, &dimid) == NC_NOERR)
        return 1;
    else
        return 0;
}

/** Copies all attributes of a specified variable from one NetCDF file to
 * another.
 *
 * @param fname_src Source file name
 * @param ncid_src Source file id
 * @param varid_src Variable id
 * @param fname_dst Destination file name
 * @param ncid_dst Destination file id
 * @param varid_dst Destination variable id
 */
void ncw_copy_atts(const char fname_src[], int ncid_src, int varid_src, const char* fname_dst, int ncid_dst, int varid_dst)
{
    char varname_src[NC_MAX_NAME] = STR_UNKNOWN;
    char varname_dst[NC_MAX_NAME] = STR_UNKNOWN;
    int status;
    int natts;
    int i;

    _ncw_inq_varname(fname_src, ncid_src, varid_src, varname_src);
    _ncw_inq_varname(fname_dst, ncid_dst, varid_dst, varname_dst);
    ncw_inq_varnatts(fname_src, ncid_src, varid_src, &natts);

    for (i = 0; i < natts; ++i) {
        char attname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_attname(fname_src, ncid_src, varid_src, i, attname);
        if ((status = nc_copy_att(ncid_src, varid_src, attname, ncid_dst, varid_dst)) != NC_NOERR)
            quit("\"%s\": -> %s:  nc_copy_att(): failed for varid_in = %d (varname = \"%s\"), attname = \"%s\", varid_dst = %d, (varname = \"%s\"): %s\n", fname_src, fname_dst, varid_src, varname_src, attname, varid_dst, varname_dst, nc_strerror(status));
    }
}

/** Defines a new variable using an existing variable as a template.
 *
 * @param fname NetCDF file name
 * @param ncid NetCDF file id
 * @param oldvarname Name of the existing variable
 * @param newvarname Name of the new variable
 */
void ncw_def_var_as(const char fname[], int ncid, const char oldvarname[], const char newvarname[])
{
    int oldvarid, newvarid;
    nc_type type;
    int ndims;
    int dimids[NC_MAX_DIMS];
    int natts;
    int i;

    ncw_inq_varid(fname, ncid, oldvarname, &oldvarid);
    ncw_inq_var(fname, ncid, oldvarid, NULL, &type, &ndims, dimids, &natts);

    ncw_def_var(fname, ncid, newvarname, type, ndims, dimids, &newvarid);

    for (i = 0; i < natts; ++i) {
        char attname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_attname(fname, ncid, oldvarid, i, attname);
        ncw_copy_att(fname, ncid, oldvarid, attname, fname, ncid, newvarid);
    }
}

// EOF
