/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/include/dumpfile.h
 *  
 *  Description:
 *  Structures and prototypes for dumpfile output.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: dumpfile.h 7525 2024-04-03 21:45:04Z her127 $
 *
 */

#if !defined(_DUMPFILE_H)
#define _DUMPFILE_H 1

#define CL_NONE   0x0000
#define CL_GRID   0x0001
#define CL_CENTRE 0x0002
#define CL_LEFT   0x0004
#define CL_BACK   0x0008
#define CL_SP2    0x0010
#define CL_SP3    0x0020
#define CL_B2     0x0040
#define CL_B3     0x0080
#define CL_B2T    0x0100
#define CL_B3T    0x0200
#define CL_FACE   0x0400
#define CL_EDGE   0x0800
#define CL_VERTEX 0x1000

#define M3PT   0x001
#define M5PT   0x002
#define W3PT   0x004
#define W5PT   0x008
#define HP3PT  0x010
#define SHU3PT 0x020

#define MAXNUMVARS 350

typedef struct dump_data dump_data_t;
typedef struct dump_file dump_file_t;

/* Land fill function mapping */
typedef void (*landfillfn_t) (dump_data_t *dumpdata);
typedef struct {
  char *name;                   /* name tag */
  landfillfn_t landfill;
} landfill_map_t;

typedef struct {
  char name[MAXSTRLEN]; /* Filter name */
  int type;             /* Filter type */
  int size;             /* Stencil size */
  double *k;            /* Stencil kernel */
  int *map;             /* Stencil map */
  int map3;             /* 1 ring stencil size */
  int map5;             /* 2 ring stencil size */
  double f;             /* Implementation factor */
  double v;             /* Damping factor (optional) */
} df_filter_t;

/* Structure to describe each dump file time dep variable */
typedef struct {
  void **v;                     /* Pointer to values */
  void **vt;                    /* Pointer to values */
  int ndims;                    /* Number of spatial dimensions */
  nc_type type;                 /* netCDF type of this variable */
  int xylocation;               /* Location on the horizontal mesh */
  int zlocation;                /* Location on the vertical mesh */
  int sediment;                 /* Sediment flag */
  int *hmap;                    /* Map to horizontal location */
  int *vmap;                    /* Map to layer */
  int *m2d;                     /* 3D to 2D map */
} df_ugrid_var_t;

/* Structure to describe each dump file */
typedef struct {
  int fid;                      /* Netcdf file id */
  int nextrec;                  /* Netcdf record number */
  df_ugrid_var_t *vars;         /* List of dump file variables */
} df_ugrid_data_t;

/*-------------------------------------------------------------------*/
/* MOM Grid structure                                                */
/*-------------------------------------------------------------------*/
typedef struct {
  FILE *prmfd;                  /* Parameter file pointer            */
  int cfid;                     /* netCDF file handle                */
  char momfile[MAXSTRLEN];      /* MOM conversion file name          */
  char ofile[MAXSTRLEN];        /* Output file name                  */
  int domom;                    /* Create MOM grid structure         */
  int nce1;                     /* Grid dimension in the e1 dir      */
  int nce2;                     /* Grid dimension in the e2 dir      */
  int nfe1;                     /* Grid face dimension in e1 dir     */
  int nfe2;                     /* Grid face dimension in e2 dir     */
  int grid_x_T;                 /* T cell gridsize x direction       */
  int grid_y_T;                 /* T cell gridsize y direction       */
  int grid_x_C;                 /* C cell gridsize x direction       */
  int grid_y_C;                 /* C cell gridsize y direction       */
  int positive;                 /* z direction                       */
  int nz;                       /* Number of layers                  */
  int vertex;                   /* Number of verticies               */
  grid_rect_t *rg;              /* Rectangular grid information      */
  grid_polar_t *pg;             /* Polar grid information            */
  double flon, flat;            /* False pole coordinates            */
  char projection[MAXSTRLEN];   /* Grid projection                   */
  char gridtype[MAXSTRLEN];     /* Grid name                         */
  char tunits[MAXSTRLEN];       /* Time units                        */
  int gridcode;                 /* Code to specify grid type         */
  double *layers;               /* Depths of the layer faces         */
  double *bathy;                /* Bathymetry array                  */
  double **topo;                /* Bathymetry array                  */
  double bmin;                  /* Minimum bathymetry (m)            */
  double bmax;                  /* Maximum bathymetry (m)            */
  double etamax;                /* Maximum elevation (m)             */
  char mct[MAXSTRLEN];          /* Minimum cell thickness            */
  double maxgrad;               /* Maximum bathymetry gradient       */
  int smooth;                   /* Bathymetry smoothing              */
  int edgef;                    /* Grid expansion                    */
  int nvals;                    /* Number of bathymetry values       */
  int ilower, jlower;           /* (i,j) start dump indices          */
  int iupper, jupper;           /* (i,j) end dump indices            */
  int cyoffs;                   /* Cyclic offset                     */
  int cytype;                   /* Cyclic type                       */  
  double **x;                   /* Grid e1 coordinate array          */
  double **y;                   /* Grid e2 coordinate array          */
  double **h1;                  /* Grid e1 metric array              */
  double **h2;                  /* Grid e2 metric array              */
  double **a1;                  /* Grid e1 angle array               */
  double **a2;                  /* Grid e2 angle array               */
  double **h1au1;               /* e1 centered e1 grid spacing (m)   */
  double **h2au2;               /* e2 centered e2 grid spacing (m)   */
  double *gridz;                /* Depth of layer faces              */
  double *cellz;                /* Depth of layer centers            */
  double **x_T;                 /* T cell centre x locations         */
  double **y_T;                 /* T cell centre y locations         */
  double **x_C;                 /* C cell centre x locations         */
  double **y_C;                 /* C cell centre y locations         */
  float *zt;                    /* Layer centre depths               */
  float *zb;                    /* Layer face depths                 */
  unsigned long ***flag;        /* Land mask                         */
  int nobc;                     /* Number of open boundaries         */
  int *obcis;                   /* Start i index of open boundaries  */
  int *obcjs;                   /* Start j index of open boundaries  */
  int *obcie;                   /* End i index of open boundaries    */
  int *obcje;                   /* End j index of open boundaries    */

  /* MOM metric quantities                                           */
  float *grid_x_Ta;
  float *grid_y_Ta;
  float *grid_x_Ca;
  float *grid_y_Ca;
  double **x_g;
  double **y_g;
  double **area;
  double **angle;
  double **ds_00_02;
  double **ds_20_22;
  double **ds_02_22;
  double **ds_00_20;
  double **ds_00_01;
  double **ds_01_02;
  double **ds_02_12;
  double **ds_12_22;
  double **ds_21_22;
  double **ds_20_21;
  double **ds_10_20;
  double **ds_00_10;
  double **ds_01_11;
  double **ds_11_12;
  double **ds_11_21;
  double **ds_10_11;
  double **ds_01_21;
  double **ds_10_12;
  double ***x_vert;
  double ***y_vert;
  double **depth;
  double **wet;
  double **num_levels;
  double **depth_c;
  double **wet_c;
  double **num_levels_c;
  double **T2, ***T3, ***T3G;
  double **C2, ***C3, ***C3G;

} momgrid_t;


/*-------------------------------------------------------------------*/
/* ROMS Grid structure                                               */
/*-------------------------------------------------------------------*/
typedef struct {
  FILE *prmfd;                  /* Parameter file pointer            */
  int cfid;                     /* netCDF file handle                */
  char romsfile[MAXSTRLEN];     /* MOM conversion file name          */
  char ofile[MAXSTRLEN];        /* Output file name                  */
  int doroms;                   /* Create MOM grid structure         */
  int nce1;                     /* Grid dimension in the e1 dir      */
  int nce2;                     /* Grid dimension in the e2 dir      */
  int nfe1;                     /* Grid face dimension in e1 dir     */
  int nfe2;                     /* Grid face dimension in e2 dir     */
  int xi_rho;                   /* x gridsize at RHO points          */
  int eta_rho;                  /* y gridsize at RHO points          */
  int xi_psi;                   /* x gridsize at PSI points          */
  int eta_psi;                  /* y gridsize at PSI points          */
  int xi_u;                     /* x gridsize at U points            */
  int eta_u;                    /* y gridsize at U points            */
  int xi_v;                     /* x gridsize at V points            */
  int eta_v;                    /* y gridsize at V points            */
  int positive;                 /* z direction                       */
  int nz;                       /* Number of z layers                */
  int nsigma;                   /* Number of sigma layers            */
  int vertex;                   /* Number of verticies               */
  grid_rect_t *rg;              /* Rectangular grid information      */
  grid_polar_t *pg;             /* Polar grid information            */
  double flon, flat;            /* False pole coordinates            */
  char projection[MAXSTRLEN];   /* Grid projection                   */
  char gridtype[MAXSTRLEN];     /* Grid name                         */
  char tunits[MAXSTRLEN];       /* Time units                        */
  int gridcode;                 /* Code to specify grid type         */
  double *layers;               /* Depths of the layer faces         */
  double *bathy;                /* Bathymetry array                  */
  double **topo;                /* Bathymetry array                  */
  double bmin;                  /* Minimum bathymetry (m)            */
  double bmax;                  /* Maximum bathymetry (m)            */
  double etamax;                /* Maximum elevation (m)             */
  char mct[MAXSTRLEN];          /* Minimum cell thickness            */
  double maxgrad;               /* Maximum bathymetry gradient       */
  double xl;                    /* Domain length in xi               */
  double el;                    /* Domain length in eta              */
  int smooth;                   /* Bathymetry smoothing              */
  int edgef;                    /* Grid expansion                    */
  int nvals;                    /* Number of bathymetry values       */
  double **x;                   /* Grid e1 coordinate array          */
  double **y;                   /* Grid e2 coordinate array          */
  double **h1;                  /* Grid e1 metric array              */
  double **h2;                  /* Grid e2 metric array              */
  double **a1;                  /* Grid e1 angle array               */
  double **a2;                  /* Grid e2 angle array               */
  double **h1au1;               /* e1 centered e1 grid spacing (m)   */
  double **h2au2;               /* e2 centered e2 grid spacing (m)   */
  double *gridz;                /* Depth of layer faces              */
  double *cellz;                /* Depth of layer centers            */
  double **thetax;              /* Grid angle to e1 axis             */
  double **thetay;              /* Grid angle to e2 axis             */
  double **lon_psi;             /* Longitude of PSI points           */
  double **lat_psi;             /* Latitude of PSI points            */
  double **lon_rho;             /* Longitude of PSI points           */
  double **lat_rho;             /* Latitude of PSI points            */
  double **lon_u;               /* Longitude of PSI points           */
  double **lat_u;               /* Latitude of PSI points            */
  double **lon_v;               /* Longitude of PSI points           */
  double **lat_v;               /* Latitude of PSI points            */
  double *zt;                   /* Layer centre depths               */
  double *zb;                   /* Layer face depths                 */
  unsigned long ***flag;        /* Land mask                         */
  int nobc;                     /* Number of open boundaries         */
  int *obcis;                   /* Start i index of open boundaries  */
  int *obcjs;                   /* Start j index of open boundaries  */
  int *obcie;                   /* End i index of open boundaries    */
  int *obcje;                   /* End j index of open boundaries    */

  float *sc_rho;                /* Sigma levels on rho points        */
  float *Cs_rho;                /* Stretching function on rho points */
  float *sc_w;                  /* Sigma levels on W points          */
  float *Cs_w;                  /* Stretching function on W points   */

  /* ROMS metric quantities                                           */
  double **pm;
  double **pn;
  double **dndx;
  double **dmde;
  double **theta;
  double ***hraw;
  double **h;
  double **mask_psi;
  double **mask_rho;
  double **mask_u;
  double **mask_v;
  double **num_levels;
  double **f;                   /* Coriolis parameter                */

} romsgrid_t;

/*-------------------------------------------------------------------*/
/* EMS output dump data structure                                    */
/*-------------------------------------------------------------------*/
struct dump_data {
  int us_type;
  int ndf;
  int nce1, nce2;
  int nfe1, nfe2;
  int nz;
  int ntr;
  int ntrS;
  int nsed;
  int sednz;
  int ns2;
  int ns3;
  int start_index;
  int face_dim;
  int edge_dim;
  int npe;
  int szmS;
  int nface2;
  int nedge2;
  int nvertex2;
  int nface3;
  int nedge3;
  int nvertex3;
  double runno;
  char runnoc[MAXSTRLEN];
  char rev[MAXSTRLEN];          /* Version number for parameter file */
  char reference[MAXSTRLEN];
  char runcode[MAXSTRLEN];
  char trl[MAXSTRLEN];
  char grid_name[MAXSTRLEN];

  double t;
  momgrid_t *momgrid;           /* MOM format grid structure */
  romsgrid_t *romsgrid;         /* ROMS format grid structure */
  int large;                    /* Dump large format file */
  tracer_info_t *trinfo_3d;     /* 3D Tracer info data structure */
  tracer_info_t *trinfo_2d;     /* 2D Tracer info data structure */
  tracer_info_t *trinfo_sed;    /* Sediment Tracer info data structure */
  char lenunit[MAXSTRLEN];      /* Length units */
  char timeunit[MAXSTRLEN];     /* Time units */
  char output_tunit[MAXSTRLEN]; /* Output time units */
  char prmname[MAXSTRLEN];      /* Parameter file name */
  char idumpname[MAXSTRLEN];    /* Input dump file name */
  char gridtype[MAXSTRLEN];     /* Grid name */
  char opath[MAXSTRLEN];        /* Output path for files */
  char trkey[MAXSTRLEN];        /* Transport keyname */
  grid_rect_t *rg;              /* Rectangular grid information */
  grid_polar_t *pg;             /* Polar grid information */
  double start_time;            /* Start time of simulation */
  double stop_time;             /* End time of simulation */
  int crf;                      /* Crash recovery flag */
  unsigned long ***flag;
  double *gridz;
  double *cellz;
  double *cellzcsed;
  double ***gridz_sed;
  double ***cellz_sed;

  double **gridx;
  double **gridy;
  double **h1agrid;
  double **h2agrid;
  double **cellx;
  double **celly;
  double **h1acell;
  double **h2acell;
  double **coriolis;
  double **botz;
  double **u1x;
  double **u1y;
  double **h1au1;
  double **h2au1;
  double **thetau1;
  double **u2x;
  double **u2y;
  double **h1au2;
  double **h2au2;
  double **thetau2;

  double **u1av;
  double **u2av;
  double **wind1;
  double **wind2;
  double **wtop;
  double **topz;
  double **eta;
  double **patm;
  double ***u1;
  double ***u2;
  double ***u;
  double ***v;
  double ***w;
  double ***u1m;
  double ***u2m;
  double ***wm;
  double ***u1vm;
  double ***u2vm;
  double ***dzu1;
  double ***dzu2;
  double ***dzcell;
  double **u1vh;
  double **u2vh;
  double **u1bot;
  double **u2bot;
  double ****tr_wc;
  double ***tr_wcS;
  double ****tr_sed;

  double ***dens;
  double ***dens_0;
  double ***Kz;
  double ***Vz;
  double **Cd;
  
  double ***origin;
  double ***pc;
  double ***qc;
  double ***rc;

  double *w1;
  double *w1s;
  double **w2s;
  double **wc;
  double **we;
  double **wv;
  double **w2;
  double ***w3;
  int *i1;     
  int ***fmap_i;
  int ***fmap_j;
  int ***vmap;
  int *i1s;
  int **i2s;
  int **i3s;
  int **ij2c;

  short *crci;
  short *clci;
  short *frci;
  short *flci;
  short *crfi;
  short *clfi;
  short *frfi;
  short *flfi;
  short *cfcj;
  short *cbcj;
  short *ffcj;
  short *fbcj;
  short *cffj;
  short *cbfj;
  short *fffj;
  short *fbfj;

  GRID_SPECS **gs;              /* Grid spec for parray interpolation */
  delaunay **d;                 /* Delaunay data structure for interpolation */
  int **c2p;                    /* Sparse to delaunay index map */
  dump_file_t *dumplist;        /* Table of dump files */
  int df_diagn_set;             /* Flag if the `diagn' flag has been set */
  int tmode;                    /* Transport mode */
  int togn;                     /* Origin for tri-linear interpolation */
  /* for one of the output files */

  char **vars;                  /* Percentile variable names */
  int nvars;                    /* Number of percentile variables */

#if defined(HAVE_SEDIMENT_MODULE)
  int do_bm;
  int do_bmphys;
#endif

#if defined(HAVE_ECOLOGY_MODULE)
  double ecodt;                 /* Ecology timestep */
  int eco_timestep_ratio;       /* Ecological: physical timestep */
  int do_eco;
  /*epi_info *einfo;*/
#endif

  master_t *master;
};

/* Data file output handle */
struct dump_file {
  char name[MAXSTRLEN];         /* Output filename */
  char tunit[MAXSTRLEN];        /* Output time unit */
  char origname[MAXSTRLEN];     /* Original output filename from prmfile */
  char type[MAXSTRLEN];         /* file type */
  char modulo[MAXSTRLEN];       /* Time modulo type */
  char irule[MAXSTRLEN];        /* Interpolation rule */
  int bpv;                      /* Bytes per value hint */
  int nvars;                    /* Number of output variables */
  char *vars[MAXNUMVARS];       /* Output variable names */
  int global;                   /* Is the file def global ? */
  int ilower;                   /* Lower i index for output */
  int jlower;                   /* Lower j index for output */
  int klower;                   /* Lower k index for output */
  int nce1;                     /* Number of i cells to output */
  int nfe1;                     /* Number of i faces to output */
  int nce2;                     /* Number of j cells to output */
  int nfe2;                     /* Number of j faces to output */
  int nz;                       /* Number of vertical cells to output */
  int nz_sed;                   /* Number of vertical cells in sediment to
                                   output */
  int ns2;                      /* 3D sparse size */
  int ns3;                      /* 2D sparse size */
  int nface2;                   /* Unstructured 2D face size */
  int nedge2;                   /* Unstructured 2D edge size */
  int nvertex2;                 /* Unstructured 2D vertex size */
  int nface3;                   /* Unstructured 3D face size */
  int nedge3;                   /* Unstructured 3D edge size */
  int nvertex3;                 /* Unstructured 3D vertex size */
  int finished;                 /* Non-zero if output has finished. */
  int append;                   /* Append */
  int osl;                      /* Interpolation method */
  double tout;                  /* Next output time */
  double tinc;                  /* Output interval */
  double sinc;                  /* 2-way nest time interval */
  double tstop;                 /* Stop time for this file */
  double tstart;                /* Start time for this file */
  int incf;                     /* Increment flag */
  void *(*create) (dump_data_t *dumpdata, dump_file_t *df);
  void (*write) (dump_data_t *dumpdata, dump_file_t *df, double t);
  void (*reset) (dump_data_t *dumpdata, dump_file_t *df, double t);
  void (*close) (dump_data_t *dumpdata, dump_file_t *df);
  landfillfn_t landfill;
  int curr_ymd[3];              /* Current year, month and day, respectively */
  int npoints;                  /* Number output points */
  double *x;                    /* X coordinates */
  double *y;                    /* Y coordinates */
  char **label;                 /* Label per point */
  int diagn;                    /* Flag if diagnostic tracers are to be *
                                   initialized after each write. */
  int long0_360;                /* Longitude in range 0 to 360 */
  int da_cycle;                 /* Data assimilation cycle flag */
  int compress;                 /* Activate compression */
  int chunk;                    /* Output in chunks */
  int flag;                     /* General purpose flag */
  int obcid;                    /* Id for obc dumps */
  double bathymask;             /* Bathymetry mask */
  df_filter_t *filter;          /* Filtering function */
  void *private_data;           /* Private data. */
};




/* Prototypes */
int dumpfiles(dump_data_t *dumpdata, double t, FILE * prmfp);
void dumpfile_init(dump_data_t *dumpdata, double t, FILE * fp, int *n,
                   dump_file_t **);
void dumpfile_close(dump_data_t *dumpdata, double t, FILE * prmfp);
int dump_create(dump_data_t *dumpdata, char *name);
int opendump_o(dump_data_t *dumpdata, char *name, int check);
int dump_choose(dump_data_t *dumpdata, int fid);
int choosedump_by_time_o(dump_data_t *dumpdata, int fid, double t);
void dump_write(dump_data_t *dumpdata, int cdfid, int ti);
void write_dump_attributes(dump_data_t *dumpdata, int cdfid,
                           nc_type fptype, int ilower, int nce1, int nfe1,
                           int jlower, int nce2, int nfe2, char *modulo,
			   char *tunit);
void write_analytic_att(dump_data_t *dumpdata, int cdfid, int vid,
                        int ilower, int ne1, int jlower, int ne2,
                        int xylocation);
void read_grid_atts_o(dump_data_t *dumpdata, int cdfid);
int check_sparse_dumpfile(geometry_t *geom, int ntsfiles, 
			  timeseries_t **tsfiles,
			  cstring * names);
void set_longitude(dump_data_t *dumpdata, dump_file_t *df, int mode);
void set_chunk_name(dump_data_t *dumpdata, dump_file_t *df, double t);
void df_parse_vars(dump_data_t *dumpdata, dump_file_t *df, char* excludevars, char* all_vars);
void dump_eta_snapshot(master_t *master, geometry_t **window,
		       window_t **windat, win_priv_t **wincon);

#endif
