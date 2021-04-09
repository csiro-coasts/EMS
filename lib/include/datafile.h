/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/datafile.h
 *
 *  \brief Header for datafile structures
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: datafile.h 6597 2020-09-03 05:29:26Z riz008 $
 */


#ifndef _DATAFILE_H
#define _DATAFILE_H

/* Enumerations */
/*
typedef enum { DFT_ASCII, DFT_NETCDF, DFT_MULTI_NETCDF } DataFileType;
*/
/* MH mempack */
typedef enum { DFT_ASCII, DFT_NETCDF, DFT_MULTI_NETCDF, DFT_MEMPACK} DataFileType;
typedef enum { AT_TEXT, AT_BYTE, AT_FLOAT, AT_DOUBLE, AT_SHORT, AT_INT }
  AttributeType;

/* Defines the set of variables understood by the datafile package.
 */
#define VT_UNKNOWN_COORD 0x00000001 /* Unsupported coordinate type */
#define VT_DATA		 0x00000002   /* Data var. */
#define VT_TIME		 0x00000004   /* Time information */
#define VT_X 		 0x00000008     /* X cartesian coord var. */
#define VT_Y 		 0x00000010     /* Y cartesian coord var. */
#define VT_Z 		 0x00000020     /* Z cartesian coord var. */
#define VT_LATITUDE	 0x00000040 /* Lat. geographic coord var. */
#define VT_LONGITUDE	 0x00000080 /* Long. geographic coord var. */
#define VT_INFERRED      0x00000100 /* Data var. with inferred coordinates */
#define VT_BATHY         0x00000200 /* Bathymetry variable                 */
#define VT_COORD         0x00000400 /* Data var. with inferred coordinates */
typedef int VariableType;

/* Define the geographic types understood by datafile */
#define GT_NONE		 0x00000001   /* No-geographic info. */
#define GT_GEOGRAPHIC	 0x00000002 /* Lat/Long. */
#define GT_PROJECTION	 0x00000004 /* Projection */
typedef int GeoType;

#define GEOGRAPHIC_TAG "geographic"

#if !defined(MAXNUMDIMS)
#define MAXNUMDIMS 10
#endif

#if !defined(MAXNUMCOORDS)
#define MAXNUMCOORDS (MAXNUMDIMS) /* must be same as dims */
#endif

#if !defined(MAXNUMCMAPS)
#define MAXNUMCMAPS (MAXNUMCOORDS)  /* must be same as coords */
#endif



/* Structures for the supported analytic coordinate systems. */

/* Structure for rectangular grid */
typedef struct {
  double x0;                    /* X origin */
  double y0;                    /* Y origin */
  double dx;                    /* X increment */
  double dy;                    /* Y increment */
  double th;                    /* rotation (from X axis, in degrees) */
  double sinth;
  double costh;
} grid_rect_t;

/* Structure for polar grid */
typedef struct {
  double x0;                    /* X origin (m) */
  double y0;                    /* Y origin (m) */
  double arc;                   /* azimuthal span of grid (in degrees) */
  double rmin;                  /* inside radius of polar grid (m). */
  double rotation;              /* rotation (from -ve X axis, in degrees) */
} grid_polar_t;

typedef struct {
  double ioffset;               /* Grid offset in i direction */
  double joffset;               /* Grid offset in j direction */
  int ni;                       /* Number grid points in i direction */
  int nj;                       /* Number grid points in j direction */
  double x0;                    /* X origin */
  double y0;                    /* Y origin */
  double dx;                    /* X increment */
  double dy;                    /* Y increment */
  double th;                    /* rotation (from X axis, in degrees) */
  double sinth;
  double costh;
} gi_rect_t;


typedef struct {
  double ioffset;               /* Grid offset in i direction */
  double joffset;               /* Grid offset in j direction */
  int ni;                       /* Number grid points in i direction */
  int nj;                       /* Number grid points in j direction */
  double x0;                    /* X origin (m) */
  double y0;                    /* Y origin (m) */
  double arc;                   /* azimuthal span of grid (in degrees) */
  double rmin;                  /* inside radius of polar grid (m). */
  double rotation;              /* rotation (from +ve X axis, in degrees) */
  double fac;                   /* Scale factor */
  double dth;                   /* delta theta */
  double logrmin;               /* log of rmin */
} gi_polar_t;


/* Multi netcdf file structure.
 */
typedef struct {
  char filename[MAXFNAMELEN];   /* Name of the file */
  size_t nrecords;              /* Number of records */
  double *records;              /* Records */
  int rec_offset;               /* Overlapping offset from the file before */
  double *offset;               /* Offset value */
  double *scale;                /* Scale value */
  double *fill;                 /* Fill value */
} df_multi_file_t;

typedef struct {
  double version;               /* Multifile version number */
  int nfiles;                   /* Number of files */
  df_multi_file_t *files;       /* The files */
  int last_fidx;                /* Index of the last file - ONLY for
				   reading trans files */
} df_multi_t;


/* MH: mempack structure containing output data */
typedef struct {
  double version;               /* Version flag */    
  int *status;                  /* Status flag */
  char id[MAXSTRLEN];           /* identifier name */
  char vars[MAXSTRLEN];         /* Variable list */
  char units[MAXSTRLEN];        /* Variable units */
  int runcode;                  /* Runcode */
  int npoints;                  /* Number of points */
  int nz;                       /* Number of layers */
  int n2d;                      /* Number of 2d variables */
  int n3d;                      /* Number of 2d variables */
  double time;                  /* Time stamp */
  double next;                  /* Next time to read */
  double dt;                    /* Time step of dump */
  char tunits[MAXSTRLEN];       /* Time units */
  double *xyz;                  /* Position data */
  char xyzunits[MAXSTRLEN];     /* Position units */
  double *v2d;                  /* 2d variables */
  double *v3d;                  /* 3d variables */
} df_mempack_t;

typedef struct {
  int nz;
  int *kn, **kmap;
  int last_record;
} df_gs_t;

#define DF_FIRST  0x00002
#define DF_WAIT   0x00004
#define DF_GO     0x00008
#define DF_RUN    0x00010
#define DF_FAIL   0x00020

typedef struct df_coord_system df_coord_system_t;
typedef struct df_coord_transform df_coord_transform_t;
typedef struct df_coord_mapping df_coord_mapping_t;
typedef struct df_geo_transform df_geo_transform_t;
typedef struct df_attribute df_attribute_t;
typedef struct df_dimension df_dimension_t;
typedef struct df_variable df_variable_t;
typedef struct datafile datafile_t;


/* An attribute contains name information, type, and number.
 */
#define CHK_TYPE(a,t) ((a)->type == (t))
#define CHK_TYPE_RANGE(a, t, i) (((a)->type == (t)) && ((i) >= 0) && ((i) < (int)((a)->n)))
#define ATT_TEXT(a) ((char*)(CHK_TYPE(a, AT_TEXT) ? (a)->value: ""))
#define ATT_BYTE(a,i) (CHK_TYPE_RANGE(a,AT_BYTE,i) ? ((char*)((a)->value))[(i)] : 0)
#define ATT_FLOAT(a,i) (CHK_TYPE_RANGE(a,AT_FLOAT,i) ? ((float*)((a)->value))[(i)] : 0.0)
#define ATT_DOUBLE(a,i) (CHK_TYPE_RANGE(a,AT_DOUBLE,i) ? ((double*)((a)->value))[(i)] : 0.0)
#define ATT_INT(a,i) (CHK_TYPE_RANGE(a,AT_INT,i) ? ((int*)((a)->value))[(i)] : 0)
#define ATT_SHORT(a,i) (CHK_TYPE_RANGE(a,AT_SHORT,i) ? ((short*)((a)->value))[(i)] : 0)

struct df_attribute {
  char *name;                   /* df_attribute_t name */
  size_t n;                     /* Number of values for this attribute */
  AttributeType type;           /* df_attribute_t type */
  void *value;                  /* df_attribute_t value(s) */
};



/* A dimension contains a name and size.
 */
struct df_dimension {
  char *name;                   /* df_dimension_t name */
  size_t size;                  /* Size of dimension */

  int dimid;                    /* netCDF dim id */
};

/* A variable can either be a data variable or a coordinate variable.
 */
#define VAR_0D(v) ((double*)((v)->data))
#define VAR_1D(v) ((double**)((v)->data))
#define VAR_2D(v) ((double***)((v)->data))
#define VAR_3D(v) ((double****)((v)->data))
#define VAR_4D(v) ((double*****)((v)->data))
/*
#define VAR_5D(v) ((double******)((v)->data))
#define VAR_6D(v) ((double*******)((v)->data))
*/

struct df_variable {

  VariableType type;            /* Enumerated variable type */

  /* Primary information */
  char *name;                   /* df_variable_t name */
  char *longname;               /* Long name */
  char *units;                  /* Units */
  double missing;               /* Missing value */
  double fillvalue;             /* Fill value */
  double scale_factor;          /* df_variable_t scale factor */
  double add_offset;            /* Offset to add */
  double lflag;                 /* Land flag */
  double cflag;                 /* Cloud flag */
  int z_is_depth;               /* Z coordinate is depth */
  int na;                       /* Number of attributes */
  df_attribute_t *attributes;   /* Unsupported attributes */

  /* Dimensional information */
  int nd;                       /* Number of dimensions */
  int *dimids;                  /* Dimensions */
  int dim_as_record;            /* If non-zero then dims include record
                                   dim */

  /* Data buffer - Always double but maybe multi-dimensional, so casting
     maybe required. */
  int start_record;             /* Record number corresponding to start of
                                   * data array.  */
  int nrecords;                 /* Number of records in record buffer */
  double *data;                 /* Array of data values (records by *
                                   array dims) */

  /* Coordinate specific information */
  char *coord_domain;           /* Defines the coordinate domain */
  /* e.g. projection=UTM(55), spheroid=WGS84 */
  int nc;                       /* Number of coordinates */
  int *coordids;                /* List of coord ids associated with *
                                   variable. NULL if none */
  df_coord_system_t *csystem;   /* Coordinate systems associated with *
                                   this variable. NULL if none */
  int nncd;                     /* Number of non-coordinate dimensions */
  int *ncdimids;                /* The number of non-coordinate dimension
                                   ids */

  /* netCDF identifiers - used for speed ups */
  int varid;                    /* Identfiers for var of netCDF file */
  /* Interpolation method */
  double (*interp) (datafile_t *df, df_variable_t *v, int record,
		    double coords[]);

  GRID_SPECS *gs0_master_2d;
  GRID_SPECS *gs1_master_2d;
  GRID_SPECS **gs0_master_3d;
  GRID_SPECS **gs1_master_3d;
  int nz;
  int *kn, **kmap, **kflag;
  int rid;
  int rv[2];
  /*
   * Buffers for parallel reads
   */
  void *nc_bufs;
};


/* Defines data and methods for transforming a coordinate coordinate
 * variables.
 */
struct df_coord_transform {
  void *private_data;           /* Transformation data. */
  void (*forward) (void *data, int n, double from[], double to[]);
  void (*inverse) (void *data, int n, double from[], double to[]);
  void (*free) (void *data);
};

/** Geographic transformation data.
  *
  * Transforms that data from the internal 'file' projection to
  * an external 'user' projection.
  */
struct df_geo_transform {
  int int_geotype;              /* The internal data type the file data
                                   should * be interpreted as. */
  map_proj_t *int_mp;           /* Projection (if required) of internal
                                   data */
  int ext_geotype;              /* The external data the data should be *
                                   tranformed into. */
  map_proj_t *ext_mp;           /* Projection (if required) of external
                                   data */
};


/* A coordinate map allows the definition of coordinate groupings
 * and the methods by which the mapping to and from indice
 * space occurs.
 */
struct df_coord_mapping {
  int nc;                       /* Number of coordinates in map */
  int *coordids;                /* Coordinate ids */
  VariableType *coordtypes;     /* Coordinate type */
  df_coord_transform_t *transform;  /* Transforming coordinates */

  /* This is private to the dfSetCoordSystem, do not set these. */
  void (*free) (datafile_t *df, /* Deallocate memory assoc with mapping */
                df_coord_mapping_t *cm);
  int (*coords_to_indices) (    /* Function to convert from coords to
                                   indices */
                             datafile_t *df, df_coord_mapping_t *cm, const int depindicies[], /* Dependent 
                                                                                                 indices 
                                                                                               */
                             const double coords[], /* Coordinate values */
                             double indices[]); /* Coordinate indicies */
  int (*indices_to_coords) (    /* Function to convert from indices to
                                   coords */
                             datafile_t *df, df_coord_mapping_t *cm, const int depindicies[], /* Dependent 
                                                                                                 indices 
                                                                                               */
                             const double indices[], double coords[]);  /* Coordinate 
                                                                           values 
                                                                         */
  void *special_data;           /* Special data which might be used by *
                                   the transformation functions */
  int nd;                       /* Number of dimensions these coords *
                                   define */
  int *dimids;                  /* df_dimension_t ids of defined coords */
  int ndepd;                    /* Number of dependent dimensions without
                                   * which the coordinate values cannot be
                                   * evaluated. */
  int *depdimids;               /* df_dimension_t ids of dependent
                                   dimensions */
};


/* A group of coordinates defines a coordinate system.
 * To convert between indice space and coordinate space
 * may be non-unique for a particular coordinate, and it
 * may be necessary to group coordinates (e.g. XY to IJ).
 * These coordinate mappings are defines as CoordinateMappings.
 *
 * Example:
 * If we have four variable T(time), X, Y, Z, and each is defined
 * using with the dimensions:
 *
 *   T(nrecords)
 *   X(ni, nj)
 *   Y(ni, nj)
 *   Z(ni, nj, nk)
 *
 * then we would have three coordinatate mappings.
 * 1 - Maps the time to record and record to time.
 * 2 - Maps XY to IJ and IJ to XY.
 * 3 - Maps Z to K knowing IJ from the XY to IJ calculation.
 *
 * In this case the simplist case is evaluated first (T),
 * followed by XY and then Z which relies on the ni and nj
 * (or X and Y) to compute Z (or K). This internal dependency
 * is critical to efficiently evaluating coordinates.
 */
struct df_coord_system {
  int nc;                       /* Number of coordinates */
  int *coordids;                /* Array of coordinate ids */
  int *coordtypes;              /* Array of coordinate types */
  int ncm;                      /* Number of coordinate mappings */
  df_coord_mapping_t *cmaps;    /* Coordinate maps. */
  int *cmorder;                 /* Defines the order in which the *
                                   coordinate should be evaluated. */
};

/* Data file data. The principle information is contained in the
 * record (if appropriate), but the data variables are stored separately.
 */
struct datafile {
  char *name;                   /* Data file name */
  DataFileType type;            /* Data file type */
  long nrecords;                /* Number of records */
  int ri;                       /* netCDF index of record variable */
  char *rec_name;               /* Record name */
  char *rec_units;              /* Record units string */
  char *rec_longname;           /* Record long name */
  int geotype;                  /* Geographic transform type */
  char *projection;             /* Projection information */
  double *records;              /* record values (if appropriate) */
  int ncid;                     /* netCDF file id, if appropriate */
  int nv;                       /* Number of variables */
  df_variable_t *variables;     /* Variables */
  int nd;                       /* Number of dimensions */
  df_dimension_t *dimensions;   /* Dimensions */
  int na;                       /* Number of attributes */
  df_attribute_t *attributes;   /* Unsupported attributes */
  int rec_modulus;              /* Non-zero if record should be moulo'd */
  double rec_mod_scale;         /* Scale of record modulo. */
  double t0;                    /* Time of previous df_eval */
  int r0, r1;                   /* Lower and upper bounds bracketing record no. */
  double frac;                  /* Fraction for record interpolation */
  int runcode;                  /* Runcode */
  void *private_data;           /* Private data for the reader */
  /*
   *  These hashtables facilitate the inverse weighting interpolation
   *  routines for performance
   */
  hash_table_t* ht_master_1d;
  hash_table_t* ht_master_2d;
  hash_table_t* ht_dist;
  hash_table_t* ht_k;

};


datafile_t *df_alloc(void);
void df_read(char *name, datafile_t *df);
void df_flush_variable(datafile_t *df, int varid);
void df_free(datafile_t *df);
void df_set_record(datafile_t *df, int varid);
double df_eval(datafile_t *df, df_variable_t *v, double r);
double df_eval_coords(datafile_t *df, df_variable_t *v, double r,
                      double coords[]);
int df_get_num_variables(datafile_t *df);
df_variable_t *df_get_variable(datafile_t *df, int varid);
df_variable_t *df_get_variable_by_name(datafile_t *df, const char *name);
const char *df_get_variable_name_by_id(datafile_t *df, int vid);
df_dimension_t *df_get_dimension(datafile_t *df, int dimid);
df_attribute_t *df_get_global_attribute(datafile_t *df, const char *name);
df_attribute_t *df_get_attribute(datafile_t *df, df_variable_t *v,
                                 const char *name);
void df_add_text_attribute(datafile_t *df, df_variable_t *v,
                           const char *name, const char *text);
int df_get_index(datafile_t *df, const char *varname);
void df_check_records(datafile_t *df);
void df_print_info(datafile_t *df, FILE * fp, int level);
int df_find_record(datafile_t *df, double r,
                   int *before, int *after, double *frac);
void df_read_all_records(datafile_t *df, df_variable_t *v);
void df_read_records(datafile_t *df, df_variable_t *v, int start_rec,
                     int nrecs);
void df_read_record(datafile_t *df, df_variable_t *v, int start_rec);
void df_fixup_variable_sign(datafile_t *df, int vid);
int df_eval_period(datafile_t *df, int vid, double r, double dr, double **ptr);
double df_get_data_value(datafile_t *df, df_variable_t *v, int record,
                         int *is);

int df_set_coord_system(datafile_t *df, df_variable_t *v,
                        int ncr, df_coord_mapping_t *requests);
int df_infer_coord_system(datafile_t *df, df_variable_t *v);
int df_get_num_dims(datafile_t *df, df_variable_t *v);
int *df_get_dim_ids(datafile_t *df, df_variable_t *v);
int df_get_num_coords(datafile_t *df, df_variable_t *v);
int *df_get_coord_ids(datafile_t *df, df_variable_t *v);
int *df_get_coord_types(datafile_t *df, df_variable_t *v);
void df_set_proj_type(datafile_t *df, char *type);
void df_set_default_proj_type(char *type);
void df_set_default_hashtable_size(int hsize);
int df_ctoi(datafile_t *df, df_variable_t *v,
            const double coords[], double indices[]);
int df_itoc(datafile_t *df, df_variable_t *v,
            const double indices[], double coords[]);
int df_is_recom(datafile_t *df);
int df_is_ugrid(datafile_t *df);
double interp_1d_inv_weight(datafile_t *df, df_variable_t *v, int record,
                            double coords[]);
double interp_2d_inv_weight(datafile_t *df, df_variable_t *v, int record,
                            double coords[]);
double interp_nearest_within_eps(datafile_t *df, df_variable_t *v, int record,
				 double coords[]);
double interp_us_2d(datafile_t *df, df_variable_t *v, int record,
		    double coords[]);
double interp_us_3d(datafile_t *df, df_variable_t *v, int record,
		    double coords[]);

/* Special threaded & buffered functions */
#ifdef __cplusplus
extern "C" {
#endif
void nc_buffered_get_vara_double(int ncid, df_variable_t *v, int varid, 
				 size_t start[], size_t count[], double *var);
void nc_buffered_clean_up_all_vars(datafile_t *df);
#ifdef __cplusplus
}
#endif
#endif                          /* _DATAFILE_H */
