/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/include/boundary.h
 *  
 *  Description:
 *  Include file for structure definitions for
 *  boundary lists
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: boundary.h 7161 2022-07-07 02:32:47Z her127 $
 *
 */

#ifndef _BOUNDARY_H
#define _BOUNDARY_H

#define NBSTD 2

typedef struct bdry_details bdry_details_t;
typedef struct scale_details scale_details_t;
typedef struct bdryfn_map bdryfn_map_t;
typedef struct bdry_old_values bdry_old_values_t;
typedef struct tidal_memory tidal_memory_t;
typedef struct tide_details tide_details_t;
typedef struct tidal_consts tidal_consts_t;
typedef struct open_bdrys open_bdrys_t;

typedef void (*bdrycustom_init_m) (master_t *master, open_bdrys_t *open,
                                   bdry_details_t *d);
typedef void (*bdrycustom_init_w) (geometry_t *window, open_bdrys_t *open,
                                   bdry_details_t *d, bdry_details_t *din);
typedef void (*bdrycustom_m) (geometry_t *geom, master_t *master,
                              open_bdrys_t *open, bdry_details_t *bdata);
typedef double (*bdrycustom_w) (geometry_t *window, window_t *windat,
                                win_priv_t *wincon, open_bdrys_t *open,
                                double t, int c, int cc,
                                bdry_details_t *data);
typedef void (*bdrytransfer) (master_t *master, open_bdrys_t *open,
                              bdry_details_t *data, geometry_t *window,
                              window_t *windat);
typedef void (*bdryfree_m) (master_t *master, bdry_details_t *data);
typedef void (*bdryfree_w) (geometry_t *window, bdry_details_t *data);

/*-------------------------------------------------------------------*/
/* Custom boundary function map                                      */
struct bdryfn_map {
  char *name;                   /* Name tag */
  bdrycustom_init_m init_m;     /* Initilisation function for master */
  bdrycustom_init_w init_w;     /* Initilisation function for slave */
  bdrycustom_m func_m;          /* Evaluation function on the master */
  bdrycustom_w func_w;          /* Evaluation function on the slave */
  bdrytransfer trans;           /* Transfer of data from master to slave */
  bdryfree_m free_m;            /* Free custom master data */
  bdryfree_w free_w;            /* Free custom window data */
};


/*-------------------------------------------------------------------*/
/* Structure to hold old boundary values                             */
struct bdry_old_values {
  double ti;                    /* Value on the boundary at time t */
  double ti1;                   /* Value at the inside cell at time t */
  double ti2;                   /* Value at the inside inside cell at time 
                                   t */
  double t1i;                   /* Value on the boundary at time t-1 */
  double t1i1;                  /* Value at the inside cell at time t-1 */
  double t1i2;                  /* Value at inside inside cell at time t-1 
                                 */
  double tjp;                   /* Tangential boundary value at b+1 time t 
                                 */
  double tjm;                   /* Tangential boundary value at b-1 time t 
                                 */
  double tj1p;                  /* Tangential boundary-1 value at b+1 time 
                                   t */
  double tj1m;                  /* Tangential boundary-1 value at b-1 time 
                                   t */
};

/*-------------------------------------------------------------------*/
/* Structure for tidal flow memory                                   */
struct tidal_memory {
  double f_ebb;                 /* First tracer value on the ebb tide */
  double l_ebb;                 /* Last tracer value on the ebb tide */
  double f_fld;                 /* First tracer value on the flood tide */
  double l_fld;                 /* Last tracer value on the flood tide */
  double avt;                   /* Average tracer value on the ebb tide */
  double avc;                   /* Average current for the ebb tide */
  double ttime;                 /* Time for the tide to turn */
  double avtc;                  /* Counter for avt */
  double avcc;                  /* Counter for avc */
  double timec;                 /* Time counter */
  int tdir;                     /* Tidal direction : 0=ebb, 1=flood */
};

/*-------------------------------------------------------------------*/
/* Structure for tidal forcing                                       */
struct tide_details {
  int ntides;                   /* Number of tidal constituents */
  char tname[MAXSTRLEN];        /* Name of tidal constituent */
  double ic;                    /* x co-ordinates of tidal amplitude */
  double jc;                    /* y co-ordinates of tidal amplitude */
  double amp;                   /* Tidal amplitude at cell (ic,jc) (m) */
  double per;                   /* Tidal period (hours) */
  double mta;                   /* Rate of modulation of tidal amplitude
                                   (cm/km) */
  double dta;                   /* Direction of progression of tide
                                   amplitude */
  double mtp;                   /* Rate of modulation of tidal phase
                                   (degrees/km) */
  double dtp;                   /* Direction of progression of tide phase */
};


/*-------------------------------------------------------------------*/
/* Structure for tidal harmonic synthesis                            */
struct tidal_consts {
  int nt;                       /* Number of tidal constituents */
  char **tname;                 /* Name of tidal constituents */
  double *sigma;                /* Frequency of tidal constituents
                                   (deg/hr) */
  double **amp;                 /* Tidal amplitude (m) */
  double **pha;                 /* Tidal phase (deg) for each constituent */
  double *yr;                   /* Years of nodal corrections */
  double **j;                   /* Nodal correction term for each
                                   constituent */
  double **v;                   /* Nodal correction term for each
                                   constituent */
  double **vpv;                 /* Nodal correction term for each
                                   constituent */
  double z0;                    /* Mean sea level above datum */
  int k;                        /* Year for nodal interpolations */
  int type;                     /* eta or velocity components */
  int *map;                     /* Map to tidal constituents */
};

typedef struct {
  int rlx;                 /* Relaxation flag */
  int tctype;              /* Time constant type */
  double rate;             /* Relaxation rate used */
  double dv0, dv1;         /* Adaptive relaxation endpoints */
  double tc0, tc1;         /* Adaptive relaxation endpoints */
  double slope;            /* Adaptive relaxation slope */
} fa_info_t;

/*-------------------------------------------------------------------*/
/* Structure containing the boundary value/custom routines for each  */
/* variable that wishes to explicity define how the boundary value   */
/* is computed, rather than taking them from a file.                 */
struct bdry_details {
  char name[MAXSTRLEN];         /* Name of boundary variable/tracer */
  int explct;                   /* Was data info. explicitly specified */
  double fill_value;            /* Default fill value */
  int nargs;                    /* Number of custom arguments */
  cstring *args;                /* String arguments */
  char custom_tag[MAXSTRLEN];   /* Custom function name */
  char i_rule[MAXSTRLEN];       /* Interpolation rule */
  int type;                     /* Type flag (e.g. normal or tangential) */
  bdrycustom_init_m init_m;     /* Master custom initialisation */
  bdrycustom_init_w init_w;     /* Slave custom initialisation */
  bdrycustom_m custom_m;        /* Custom function for the master */
  bdrycustom_w custom_w;        /* Custom function for the slave */
  bdrytransfer trans;           /* Transfer data from master to slave */
  bdryfree_m free_m;            /* Free custom master data */
  bdryfree_w free_w;            /* Free custom window data */
  void *custdata;               /* Data to be used by the custom func. */
};

/*-------------------------------------------------------------------*/
/* Structure containing the boundary scaling information for each    */
/* tracer variable that wishes to post-scale its calculated value.   */
/* Scaling can be accomplished by:                                   */
/* 1. Using a supplied scalar value                                  */
/* 2. Using a the values from a specified tracer variable in the     */
/*    tracer list. In this case the name of the tracer is supplied.  */
/*    The scaling tracer values may be updated with the reset        */
/*    function.                                                      */
/* 3. Supplying a custom routine for the scaling (not implemented).  */
struct scale_details {
  char name[MAXSTRLEN];         /* Name of tracer                    */
  int type;                     /* Flag for scale method             */
  int ntr;                      /* Scaling tracer number             */
  double fact;                  /* Scaling factor                    */
  double val;                   /* Endpoint value                    */
  int flag;                     /* Scaling flags                     */
  int nargs;                    /* Number of custom arguments */
  cstring *args;                /* String arguments */
  char custom_tag[MAXSTRLEN];   /* Custom function name */
  bdrycustom_w custom_w;        /* Custom function for the slave */
  void *custdata;               /* Data to be used by the custom func. */
};


/*------------------------------------------------------------------*/
/* Open boundary data structure                                     */
/*------------------------------------------------------------------*/
struct open_bdrys {
  char name[MAXSTRLEN];         /* Name of the open boundary */
  int type;                     /* U1BDRY, U2BDRY or VELOCITY */
  int id;                       /* Number identifier */
  int ocodex;                   /* Orientation of edge in x direction */
  int ocodey;                   /* Orientation of edge in y direction */
  int ocodec;                   /* Orientation of centre from edge */
  int *dir;                     /* Direction of edge vector (1=in, -1=out) */
  int *outi;                    /* Index pointing into the domain */
  int *ini;                     /* Index pointing out of the domain */
  int *ceni;                    /* Orientation of centre from edge */
  int *inc;                     /* Index pointing in for centres */
  int **bec;                    /* Active boundary edges for each centre */
  int **bcc;                    /* Active boundary centres for each centre */
  int *nmap;                    /* Centre map in an interior direction normal to bdry */
  int *omap;                    /* Centre map in an exterior direction normal to bdry */
  int **nmape;                  /* Edge map in an interior direction normal to bdry */
  int **omape;                  /* Edge map in an exterior direction normal to bdry */
  int *tmpp;                    /* Map in a + tangential dir. to bdry */
  int *tmpm;                    /* Map in a - tangential dir. to bdry */
  int no3_t;                    /* Number of 3D cells for tracers */
  int no2_t;                    /* Number of 2D cells for elevation */
  int no2_a;                    /* Number of 2D auxiliary cells */
  int no3_a;                    /* Number of 3D auxiliary cells */
  int no2_ta;                   /* Number of 2D tangential auxiliary cells */
  int *obc_t;                   /* Sparse OBC cells for tracers */
  int *oi1_t;                   /* Sparse 1 interior cell for tracers */
  int *oi2_t;                   /* Sparse 2 interior cells for tracers */
  int *cyc_t;                   /* Sparse cyclic OBC cells for tracers */
  int *ogc_t;                   /* Sparse ghost OBC cells for tracers */
  int *obc_a;                   /* Boundary auxiliary cells */
  int *obc_ta;                  /* Boundary tangential auxiliary cells */
  int *bot_t;                   /* Bottom cell centered sparse coordinate */
  int no3_e1;                   /* Number of 3D cells for u1 velocity */
  int no2_e1;                   /* Number of 2D cells for u1 velocity */
  int to3_e1;                   /* No. 3D cells for u1 tangential vel.  */
  int to2_e1;                   /* No. 2D cells for u1 tangential vel.  */
  int *obc_e1;                  /* Sparse OBC cells for u1 velocity */
  int *oi1_e1;                  /* Sparse 1 interior cell for u1 */
  int *oi2_e1;                  /* Sparse 2 interior cells for u1 */
  int *cyc_e1;                  /* Sparse cyclic OBC cells for u1 */
  int no3_e2;                   /* Number of 3D cells for u2 velocity */
  int no2_e2;                   /* Number of 2D cells for u2 velocity */
  int to3_e2;                   /* No. 3D cells for u2 tangential vel.  */
  int to2_e2;                   /* No. 2D cells for u2 tangential vel.  */
  int *obc_e2;                  /* Sparse OBC cells for u2 velocity */
  int *oi1_e2;                  /* Sparse 1 interior cell for u2 */
  int *oi2_e2;                  /* Sparse 2 interior cells for u2 */
  int *cyc_e2;                  /* Sparse cyclic OBC cells for u1 */
  int *e2c_e1;                  /* Edge to centre map */
  int ntflx;                    /* Number of tracers with flux OBCs */
  int *tflx;                    /* Tracers with flux OBCs */
  double *transfer_u1;          /* Data for master - slave transfers */
  double *transfer_u2;          /* Data for master - slave transfers */
  double *transfer_u1av;        /* Data for master - slave transfers */
  double *transfer_u2av;        /* Data for master - slave transfers */
  double *transfer_eta;         /* Data for master - slave transfers */
  int *tmap;                    /* Map from master-slave OBC vectors */
  int *tmap_u1;                 /* Map from master-slave OBC vectors */
  int *tmap_u2;                 /* Map from master-slave OBC vectors */
  int ntt;                      /* No. of tracers requiring transfers */
  double **t_transfer;          /* Data for tracer m - s transfers */
  int bgz;                      /* OBC ghost zone */
  int ttsz;                     /* Tracer transfer vector size */
  int *t_tmap;                  /* Map from tracer m - s OBC vectors */
  int **t_imap;                 /* Map from boundary - transfer index */
  int *trm;                     /* Tracers requiring m - s transfers */
  int mwn;                      /* Window number in the master */
  int inverse_barometer;        /* 1 if enabled, else 0 */
  int upmeth;                   /* Velocities to use in UPSTRM */
  int bcond_nor;                /* Normal boundary condition to invoke */
  int bcond_nor2d;              /* Normal boundary condition to invoke */
  int bcond_tan;                /* Tangential boundary condition */
  int bcond_tan2d;              /* Tangential boundary condition */
  int bcond_ele;                /* Elevation boundary condition */
  int bcond_w;                  /* Vertical velocity boundary cond */
  int bcond_Vz;                 /* Vertical viscosity boundary cond */
  int bcond_Kz;                 /* Vertical diffusivity boundary cond */
  int *bcond_tra;               /* Tracer boundary condition to invoke */
  double *clampv;               /* Value for clamped OBC's */
  int ntr;                      /* Number of tracers to apply OBCs */
  int atr;                      /* Number of automatically defined tracers */
  double relax_time;            /* Forcing relaxation time scale */
  double relax_timei;           /* Incoming waves relaxation time scale */
  double *trpc;                 /* Pycnocline depth for multi-BCs */
  int linear_zone_nor;          /* Linear normal velocity zone */
  int linear_zone_tan;          /* Linear tangential velocity zone */
  int relax_zone_nor;           /* Normal velocity flow relaxation zone */
  int relax_zone_tan;           /* Tangential velocity relaxation zone */
  int relax_zone_ele;           /* Elevation flow relaxation zone */
  int *relax_zone_tra;          /* Tracer flow relaxation zone */
  char nzone[MAXSTRLEN];        /* Nudging zone for T and S */
  double *rtra_b;               /* Relaxation timescale on boundary */
  double *rtra_i;               /* Relaxation timescale in interior */
  int relax_ele;                /* Elevation relaxation zone */
  double rele_b;                /* Relaxation timescale on boundary */
  double rele_i;                /* Relaxation timescale in interior */
  double rnor_b;                /* Relaxation timescale on boundary */
  double rnor_i;                /* Relaxation timescale in interior */
  double rtan_b;                /* Relaxation timescale on boundary */
  double rtan_i;                /* Relaxation timescale in interior */
  int sponge_zone;              /* Sponge zone for bottom friction */
  double sponge_zone_h;         /* Sponge zone for horizontal friction */
  double sponge_f;              /* Multiplication factor for sponges */
  int nspc, nspe1;              /* Number of cells / edges in the zone */
  int *spc, *spe1;              /* Cells / edges in the zone */
  double *swc, *swe1;           /* Weights for linear interpolation */
  int *snc, *smc;               /* Closest cell to boundary / outer perimeter */
  int *sne1, *sme1;             /* Closest edge to boundary / outer perimeter */
  double ncells;                /* Number of cells in this boundary */
  double *nepc;                 /* 1 / Number of edges per cell */
  double area;                  /* Area of boundary */
  double mindep, maxdep, meandep; /* Minimum, maximum, mean depth */
  double length;                /* Length of boundary */
  int cdeep;                    /* Index for deepest boundary location */
  int bathycon;                 /* Constant boundary depth flag */
  int npts;                     /* Number of cells for this boundary */
  int stagger;                  /* Normal velocity stagger arrangement */
  int *iloc;                    /* x location of the open boundary */
  int *jloc;                    /* y location of the open boundary */
  int ncyc;                     /* Number of CYCLED cells for this boundary */
  int *ilocc;                   /* x location of CYCLED cells */
  int *jlocc;                   /* y location of CYCLED cells */
  int *locu;                    /* Unstructured cell location */
  double **posx;                /* Unstructured x edge coordinates */
  double **posy;                /* Unstructured y edge coordinates */
  int bout;                     /* Boundary zone to set to OUTSIDE */
  bdry_old_values_t *etabold;   /* Inside cells and backward time eta */
  bdry_old_values_t *u1_on;     /* Inside cells and backward time u1 */
  bdry_old_values_t *u2_on;     /* Inside cells and backward time u2 */
  bdry_old_values_t *u1_ot;     /* Inside cells and backward time u1 */
  bdry_old_values_t *u2_ot;     /* Inside cells and backward time u2 */
  bdry_old_values_t *u1av_on;   /* Inside cells and backward time u1av */
  bdry_old_values_t *u2av_on;   /* Inside cells and backward time u2av */
  bdry_old_values_t *u1av_ot;   /* Inside cells and backward time u1av */
  bdry_old_values_t *u2av_ot;   /* Inside cells and backward time u2av */
  double tidemem_depth;         /* Depth below which to invoke TIDALF */
  int ntide;                    /* Number of tidal constituents */
  int options;                  /* Additional options */
  tide_details_t *tideforce;    /* Tidal forcing data structure */
  tidal_memory_t *tide;         /* Tidal memory data structure */
  tidal_consts_t tc;            /* Tidal harmonic structure */
  tidal_consts_t tun;           /* Tidal harmonic normal u velocity structure */
  tidal_consts_t tvn;           /* Tidal harmonic normal v velocity structure */
  tidal_consts_t tut;           /* Tidal harmonic tangential u velocity structure */
  tidal_consts_t tvt;           /* Tidal harmonic tangential v velocity structure */
  char tide_con[MAXSTRLEN];     /* Custom tidal constituent list */
  char *scale_s[MAXNUMVARS];    /* Scaling function (sum) for tracers */
  char *scale_p[MAXNUMVARS];    /* Scaling function (product) for tracers */
  char *scale_d[MAXNUMVARS];    /* Scaling function (density) for tracers */
  char scale_e[MAXSTRLEN];      /* Scaling function for eta */
  scale_details_t *sdata_t;     /* Scaling data for tracers */
  scale_details_t *sdata_e;     /* Scaling data for eta */
  double adjust_flux;           /* Flag to adjust normal flow */
  double adjust_flux_s;         /* Flag to adjust normal flow (short time-scale) */
  double afr;                   /* Adjust flux ratio relative to 2D timestep */
  fa_info_t *fas;               /* Flux adjustment relaxation structure */
  double bflux_2d;              /* 2D normal integrated flux */
  double bflux_3d;              /* 3D normal integrated flux */
  double meanc;                 /* Time counter for means */
  double *sphase;               /* Temporally smoothed phase */
  double spf;                   /* Phase smoothing factor */
  char bflow[MAXSTRLEN];        /* Boundary flow data filename */
  double bhc;                   /* Boundary flow halocline depth */
  double rlen;                  /* Boundary river length */
  double *dsv;                  /* Bottom density scaling value */
  double file_dt;               /* Time increment for FILEIN input */
  double file_next;             /* Time of next FILEIN input */
  int flag;                     /* General purpose flag */
  int bstdf;                    /* Standard OBC in use */
  int sbcond;                   /* Standard boundary code */
  int nbstd;                    /* Number of standard OBCs */
  int smooth_z;                 /* Boundary smoothing zone */
  int smooth_n;                 /* Boundary smoothing passes */
  char **bstd;                  /* Standard OBC definition */
  int *olap;                    /* Number of wet faces in the cell */
  double *u1d, *u2d;            /* u1 and u2 velocity dummies */
  double *flow;                 /* River flow buffer */
  int intype;                   /* Type of input specification */
  double slon, slat;            /* Start lon / lat for boundary segment */
  double elon, elat;            /* End lon / lat for boundary segment */
  double mlon, mlat;            /* Mid lon / lat for boundary segment */
  double *dum, *dum1;           /* Dummy array */
  double **dumn, **dumtr;       /* Dummy array */
  double *d1, *d2, *d3;
  double v1, v2, v3;
  int *i1;

  double maxlat, minlat;
  double maxlon, minlon;
  int nedges, *edges;
  char i_rule[MAXSTRLEN];       /* Interpolation rule          */

  /* Time series input data */
  char tsfn[MAXSTRLEN];         /* Time seris Filenames */
  int ntsfiles;                 /* Number of 'global' timeseries files */
  timeseries_t **tsfiles;       /* 'global' timeseries_t files */
  cstring *filenames;           /* 'global' filenames */

  /* Custom data */
  int custype;                  /* Type of custom routine (u1 or u2) */
  char cusname_u1[MAXSTRLEN];   /* Name of u1 custom function */
  char cusname_u2[MAXSTRLEN];   /* Name of u2 custom function */
  char cusname_u1av[MAXSTRLEN]; /* Name of u1av custom function */
  char cusname_u2av[MAXSTRLEN]; /* Name of u2av custom function */
  char cusname_eta[MAXSTRLEN];  /* Name of eta custom function */
  char *cusname_t[MAXNUMVARS];  /* Name of tracer custom function */
  bdry_details_t datau1;        /* Boundary data for u1 */
  bdry_details_t datau2;        /* Boundary data for u2 */
  bdry_details_t datau1av;      /* Boundary data for u1av */
  bdry_details_t datau2av;      /* Boundary data for u2av */
  bdry_details_t etadata;       /* Boundary data for eta */
  bdry_details_t *bdata_t;      /* Boundary data for tracers */
};


/*-------------------------------------------------------------------*/
/* Custom data associated with the u2 flow boundary                  */
typedef struct {
  double *v_river;              /* River flow velocity profile */
  timeseries_t *tsflow;         /* Time series river flow */
  int flowid;                   /* River flow id */
  double flow;                  /* Flow rate at time t */
  double hc;                    /* Depth of surface layer */
  double rlen;                  /* River length */
  double ncells;                /* Number of OBC cells */
  int init;                     /* Non-zero if initialisation was
                                   completed */
  int ocode;                    /* Orientation code */
  int options;                  /* Riverflow options */
} flow_data_t;


/*-------------------------------------------------------------------*/
/* Custom data for tracer function reconstructions                   */
typedef struct {
  char name[MAXSTRLEN];
  int mid;
  int ncells;
  int *cells;
  double *x;
  double dist;
  double tstart;
  int flag;
} tra_data_t;

/*-------------------------------------------------------------------*/
/* MIKE : Include boundary routines                                  */

double bc_orlanski(double F, double Fi1, double Fi2, double FFi1,
                   double FB, double FBi1, double *ps, double sf);
double bc_camerlengo(double F, double Fi1, double Fi2, double FFi1,
                     double FB, double FBi1, double *ps, double sf);
double bc_miller(double F, double Fi1, double Fi2, double FFi1, double FB,
                 double FBi1, double FBi2, double *ps, double sf);
double bc_raymond(double F, double Fi1, double Fi2, double FFi1,
                  double FFi2, double Fjp, double Fjm, double Fj1p,
                  double Fj1m, double *ps, double sf);
double bc_gravity(geometry_t *window, window_t *windat, win_priv_t *wincon,
                  double hat, double dt, int c2, double F, double FF, 
		  double *ps, double sf);
double bc_leastsq(double y1, double y2, double y3, double y4);
double bc_polint(double y1, double y2, double y3);
double bc_tidal(geometry_t *window, window_t *windat,
                win_priv_t *wincon, tide_details_t *tide, int s);
void bc_nograd_2d(double *vel, double *depth, double *md, double *h,
		  double hmin, int sb, int eb, int *obc, int *oi1);
void u1av_local(geometry_t *window, window_t *windat, 
		win_priv_t *wincon, open_bdrys_t *open, double *vel,
		int sb, int eb, int dir);
void u2av_local(geometry_t *window, window_t *windat, 
		win_priv_t *wincon, open_bdrys_t *open, double *vel,
		int sb, int eb, int dir);

/* END MIKE                                                          */
/*-------------------------------------------------------------------*/

#endif
