/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/include/hd_params.h
 *  
 *  Description:
 *  Include file for defining various
 *  program parameters for
 *  meco model
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: hd_params.h 6012 2018-10-31 22:43:43Z riz008 $
 *
 */

/* sizes */

#define MAXNUMTSFILES 20
#define MAXNUMARGS 60

/* z co-ordinate at top of top level */
#define	MAXGRIDZ	1e20

/* Value for pi */
#ifndef PI
#define PI (3.14159265358979323846)
#endif

/* Von Karmen's constant */
#define VON_KAR 0.408

/* Tolerence for timesteps */ 
#define DT_EPS 0.1

/* Fraction of hmin that defines a dry cell */
#define DRY_FRAC 0.05

#define NOTVALID -9999
#define LANDCELL 99

/* face flags */
#define U1SOLID		0x00000001    /* u1 point is a solid wall */
#define U2SOLID		0x00000002    /* u2 point is a solid wall */
#define U1OUTSIDE	0x00000004    /* u1 point outside model */
#define U2OUTSIDE	0x00000008    /* u2 point outside model */
#define U1BDRY		0x00000040    /* u1 point open boundary */
#define L_EDGE		0x00000080    /* u1 bdry point on left side of interior
                                   cell */
#define R_EDGE		0x00000100    /* u1 bdry point on right side of interior 
                                   cell */
#define U2BDRY		0x00000200    /* u2 point open boundary */
#define B_EDGE		0x00000400    /* u2 bdry point on back side of interior
                                   cell */
#define F_EDGE		0x00000800    /* u2 bdry point on forward side of
                                   interior cell */
#define U1AVZERO	0x00001000    /* columns either side of u1av point dry */
#define U2AVZERO	0x00002000    /* columns either side of u2av point dry */

/* cell flags */
#define ETASPEC		0x00004000    /* specified elevation */
#define DRY		0x00008000        /* contains no water */
#define	SOLID		0x00010000      /* solid */
#define	OUTSIDE		0x00020000    /* outside model volume */
#define	TRACERSPEC	0x00040000  /* tracer specified */

/* Grid corner flags */
#define	ALLWATER	0x00080000    /* There is water in all 4 cells
                                   surrounding this corner */
#define	MID_B_EDGE	0x00400000  /* On back boundary with water in front */
#define	MID_F_EDGE	0x00800000  /* On front boundary with water behind */
#define	MID_L_EDGE	0x01000000  /* On left boundary with water to right */
#define	MID_R_EDGE	0x02000000  /* On right boundary with water to left */

/* specified tracer face flux flags */
#define U1TFLUX		0x00100000    /* tracer flux defined at u1 face */
#define U2TFLUX		0x00200000    /* tracer flux defined at u2 face */

#define U1GEN           0x04000000  /* General purpose u1 flag */
#define U2GEN           0x08000000  /* General purpose u2 flag */

#define LO_EDGE         0x10000000  /* u1 boundary is left overwrite cell */
#define BO_EDGE         0x20000000  /* u2 boundary is back overwrite cell */

#define PRESERVE_B_FACE 0x00000001
#define PRESERVE_F_FACE 0x00000002
#define PRESERVE_L_FACE 0x00000004
#define PRESERVE_R_FACE 0x00000008

/* all the u1 flags */
#define U1FLAGS (U1SOLID|U1OUTSIDE|U1BDRY|U1TFLUX|L_EDGE|R_EDGE)
/* all the u2 flags */
#define U2FLAGS (U2SOLID|U2OUTSIDE|U2BDRY|U2TFLUX|B_EDGE|F_EDGE)

/* flags for surfacevelocity() routine */
#define	CALCWTOP	0
#define CALCDETADT	1

#define NONE            0x80000000  /* No action taken */

/* Number of auxiliary cells to include laterally */
#define laux 2           
       
/* Grid refinement codes */                                   
#define NOZOOM    0
#define DOZOOM    1
#define ZF       -1             /* Code for unused cells in a zoom zone */
#define ZN        1
#define ZC        2
#define ZE1       4
#define ZE2       8
#define ZCB       16
#define ZE1B      32
#define ZE2B      64
#define ZE1NB     128
#define ZE2NB     256
#define ZE1TB     512
#define ZE2TB     1024
#define PRECOND   2048
#define ZZ        4096   

/* Explicit mapping codes */
#define E1_MAP    1
#define E2_MAP    2
#define E1_INNER  4
#define E2_INNER  8

/* Boundary condition flags                                                  */
#define bcnum    32             /* Number of boundary conditions supported   */
#define maxbc    0x80000000     /* Maximum boundary condition number */
#define FILEIN   0x000001       /* Boundary data input from file (clamped
                                   B.C.) */
#define CUSTOM   0x000002       /* Custom OBC */
#define NOGRAD   0x000004       /* No gradient boundary condition */
#define LINEXT   0x000008       /* Least squares linear extrapolation */
#define POLEXT   0x000010       /* Polynomial extrapolation (2nd order) */
#define GRAVTY   0x000020       /* Gravity wave radiation boundary
                                   condition */
#define ORLANS   0x000040       /* Orlanski radiation boundary condition */
#define CAMOBR   0x000080       /* Camerlengo & O'Brien boundary condition 
                                 */
#define MILLER   0x000100       /* Miller & Thorpe boundary condition */
#define VERTIN   0x000200       /* Vertical integral of 3D velocity */
#define CYCLIC   0x000400       /* Cyclic boundary condition */
#define NOTHIN   0x000800       /* No boundary condition */
#define RAYMND   0x001000
#define TIDEBC   0x002000       /* Tidal elevation prescription */

#define CLAMPD   0x004000       /* Clamped to zero boundary condition */
#define STATIS   0x008000       /* Statistical prescription */
#define UPSTRM   0x010000       /* Upstream advection boundary condition */
#define CYCLED   0x020000       /* Cyclic via spatial maps */
#define TIDALH   0x040000       /* Tidal harmonics boundary condition
                                   (eta) */
#define TIDALC   0x080000       /* Custom tidal harmonics boundary condition
                                   (eta) */
#define PROFIL   0x100000       /* Tracer profile */
#define DEPROF   0x200000       /* Density matched tracer profile */
#define DESCAL   0x400000       /* Density profile tracer scaling */
#define FLATHR   0x800000       /* Flather boundary condition */
#define TRFLUX  0x1000000       /* Tracer flux */
#define TRCONC  0x2000000       /* Tracer concentration */
#define TRCONF  0x4000000       /* Tracer concentration * volume flux */
#define LINEAR  0x8000000       /* Linear boundary condition */
#define LOCALN 0x10000000       /* Local normal velocity boundary condition */
#define LOCALT 0x20000000       /* Local tamgential velocity boundary condition */
#define LOCALE 0x40000000       /* Local elevation boundary condition */
#define FLATHE 0x80000000       /* Elevation Flather component */
/*#define TIDALM 0x008000*/         /* Tidal memory boundary condition */
#define TIDALM 0x100000000         /* Tidal memory boundary condition */

/* Custom boundary flags */
#define NOR             0x001
#define TAN             0x002
#define INFACE          0x004
#define OUTFACE         0x008
#define INFACET         0x010
#define NEST1WAY        0x020
#define NEST2WAY        0x040
#define NEST_STD        0x080
#define NEST_ETA        0x100
#define RIVER           0x200
#define NEST_RAD        0x400
/* NOTHIN = 0x800 */
#define NEST_FLA        0x1000
#define NEST_BARO       0x2000
#define NEST_CPD        0x4000

/* Optional boundary flags */
#define OP_UPSTRM       0x0001
#define OP_GEOSTR       0x0002
#define OP_RLEN         0x0004
#define OP_DYNAHC       0x0008
#define OP_NOHDIF       0x0010
#define OP_YANKOVSKY    0x0020
#define OP_NOSALT       0x0040
#define OP_FDEPTH       0x0080
#define OP_IFRESH       0x0100
#define OP_TRUNCL       0x0200
#define OP_MULTF        0x0400
#define OP_ETAFIL       0x0800
#define OP_OBCCM        0x1000
#define OP_ISPNG        0x2000
#define OP_PARFLOW      0x4000
#define OP_NESTBARO     0x8000
#define OP_MACREADY     0x01000
#define OP_OWRITE       0x02000

/* Tidal boundary flags */
#define TD_ETA          0x0002
#define TD_NU           0x0004
#define TD_NV           0x0008
#define TD_TU           0x0010
#define TD_TV           0x0020

/* Stability adjustment method flags                                         */
#define ORIGINAL         1       /* Original stability adjustment */
#define SUB_STEP         2       /* Sub-time stepping stability adjustment */
#define SUB_STEP_NOSURF  4       /* SUB-STEP with no u1/u2 Courant
                                   sub-timestep */
                             /* calculation in the surface layer.  */
#define SUB_STEP_TRACER  16      /* SUB-STEP for tracers only */
#define MULTI_DT         32      /* Multi time-steps for windows */

/* Mixed layer depth flags                                                   */
#define DENS_MIX        1
#define TKE_MIX         2
#define TEMP_MIX        4

/* Heat flux options                                                         */
#define ADVANCED        2
#define INVERSE         4
#define NET_HEAT        16
#define SURF_RELAX      32
#define COMP_HEAT       64
#define AVHRR           128
#define COMP_HEAT_MOM   256
#define COMP_HEAT_NONE  512

/* Salt flux options                                                         */
#define BULK            4

/* Mean diagnostic options                                                   */
#define FLUX            0x0001
#define VEL3D           0x0002
#define WIND            0x0004
#define TENDENCY        0x0008
#define VEL2D           0x0010
#define TIDAL           0x0020
#define TRANSPORT       0x0040
#define RESET           0x0080
#define ETA_M           0x0100
#define TS              0x0200
#define KZ_M            0x0400
#define VOLFLUX         0x0800
#define PSSFLUX         0x1000
#define MMM             0x2000
#define MTRA3D          0x4000
#define MTRA2D          0x8000

/* Tracer increment flags */
#define U1VEL           1
#define U2VEL           2
#define TEMP            4
#define SALT            16

/* Window transfer flags                                                     */
#define NVELOCITY       1
#define VELOCITY        2
#define MIXING          4
#define WVEL            8
#define TMASS           16
#define DEPTH           32
#define U2FLUX          64
#define CFL             128
#define TRACERS         256

/* Master data memory free flags                                             */
#define ALL             1
#define UNUSED          2

/* Start options                                                             */
#define AUTO            1
#define MANUAL          2
#define DUMP            4
#define EXIT            8
#define ROAM            16
#define RE_ROAM         32
#define RE_NOTS         64
#define TRANS           128
#define PRE_MARVL       256

/* Auto config options */
#define A_ROAM_CPD1     2
#define A_ROAM_CPD2     4
#define A_ROAM_FLA      8
#define A_ROAM_R1       16
#define A_ROAM_R2       32
#define A_RECOM_R1      64
#define A_RECOM_R2      128
#define A_ROAM_R3       256

/* Advection scheme flags                                                    */
#define ORDER1        0x000002
#define ORDER2        0x000004
#define QUICKEST_AD   0x000008
#define ORDER4        0x000010
#define QUICKEST      0x000020
#define QUICKEST_CO   0x000040
#define VANLEER       0x000080
#define LAGRANGE      0x000100
#define HIORDER       254
#define ANGULAR       0x000200
#define ANGULAR3D     0x000400
#define EXPLICIT      0x000800
#define WIMPLICIT     0x001000
#define ORDER2_UW     0x002000
#define ADVECT_FORM   0x004000
#define WTOP_O2       0x008000
#define WTOP_O4       0x010000
#define ZERO_DRYK     0x020000
#define SHAPIRO       0x040000
#define FFSL          0x080000

/* Momentum ommission flags                                                  */
#define ADVECT        1
#define HDIFF         2
#define VDIFF         4
#define PRESS_BT      16
#define PRESS_BC      32
#define CORIOLIS      64

/* Debugging flags */
#define D_ADVECT      0x000001
#define D_HDIFF       0x000002
#define D_POST        0x000004
#define D_BDRY        0x000008
#define D_INIT        0x000010
#define D_UA          0x000020
#define D_VA          0x000040
#define D_ETA         0x000080
#define D_TS          0x000100
#define D_U           0x000200
#define D_V           0x000400
#define D_PRESSURE    0x000800
#define D_VZ          0x001000
#define D_APPEND      0x002000
#define D_AFTER       0x004000
#define D_OPEN        0x008000
#define D_STRML       0x010000
#define D_STEP        0x100000

/* CFL options                                                               */
#define PASSIVE       1
#define ACTIVE        2
#define ACTIVE3D      4
#define CFL_DIAG      8
#define CFL_ALERT     16

/* Vorticity options                                                         */
#define RELATIVE      1
#define ABSOLUTE      2
#define POTENTIAL     4

/* Tracer type options */
#define SEDIM         0x0001
#define INTER         0x0002
#define WATER         0x0004
#define HYDRO         0x0008
#define SEDIMENT      0x0010
#define ECOLOGY       0x0020
#define WAVE          0x0040
#define TRACERSTAT    0x0080
#define PROGNOSTIC    0x0100
#define DIAGNOSTIC    0x0200
#define PARAMETER     0x0400
#define FORCING       0x0800
#define E1VAR         0x1000
#define E2VAR         0x2000
#define CLOSURE       0x4000

/* Velocity faces / centers to use in UPSTRM condition */
#define CENTER        1
#define FACE          2
#define INTERIOR      4
#define ADAPTIVE      8
#define GHOST         16

/* Alert levels */
#define LEVEL1_U      1
#define LEVEL1_UA     4
#define LEVEL1_V      8
#define LEVEL1_VA     16
#define LEVEL2        32
#define LEVEL2w       64
#define ETA_A         128
#define LEVEL1_W      256
#define LEVEL1_U1VH   512
#define LEVEL1_U2VH   1024
/* Size of nalert[] */
#define nna           8

/* Numbers */
#define BRUNT         0x000001
#define INT_WAVE      0x000002
#define RICHARD_GR    0x000004
#define RICHARD_FL    0x000008
#define REYNOLDS      0x000010
#define FROUDE        0x000020
#define ROSSBY_IN     0x000040
#define ROSSBY_EX     0x000080
#define SOUND         0x000100
#define SHEAR_V       0x000200
#define BUOY_PROD     0x000400
#define SHEAR_PROD    0x000800
#define SPEED_3D      0x001000
#define SPEED_2D      0x002000
#define WIND_CD       0x004000
#define OBC_PHASE     0x008000
#define SPEED_SQ      0x010000
#define SIGMA_T       0x020000
#define ENERGY        0x040000
#define WET_CELLS     0x080000
#define SURF_LAYER    0x100000
#define SLOPE         0x200000
#define DUMMIES       0x400000
#define KINETIC       0x800000
#define BOTSTRESS    0x1000000
#define UNIT         0x2000000
#define EKPUMP       0x4000000
#define GLIDER       0x8000000
#define PASS         0x10000000

/* Wind input */
#define SPEED         2
#define STRESS        4
#define LARGEPOND     8
#define BUNKER        16
#define KITIAG        32
#define KONDO         64
#define MASAG         128

/* Horizontal mixing flags */
#define U1_SP          0x000001  /* u1vh points to sdc */
#define U2_SP          0x000002  /* u2vh points to sdc */
#define U1_SA          0x000004  /* Smagorinsky, u1vh allocated, sponges */
#define U2_SA          0x000008  /* Smagorinsky, u2vh allocated, sponges */
#define U1_A           0x000010  /* No Smagorinsky, u1vh allocated */
#define U2_A           0x000020  /* No Smagorinsky, u2vh allocated */
#define U1_SPK         0x000040  /* u1kh points to sdc */
#define U2_SPK         0x000080  /* u2kh points to sdc */
#define U1_SAK         0x000100  /* Smagorinsky, u1kh allocated */
#define U2_SAK         0x000200  /* Smagorinsky, u2kh allocated */
#define U1_AK          0x000400  /* No Smagorinsky, u1kh allocated */
#define U2_AK          0x000800  /* No Smagorinsky, u2kh allocated */
#define NONLIN         0x002000  /* Nonlinear diffusion scaling */
#define SMAG           0x004000  /* Smagorinsky diffusion */

/* Save input forcing flags */
#define OTEMP          1
#define OSALT          2
#define OETA           4
#define OVELU          128
#define OVELV          256

/* Eta relaxation flags */
#define RELAX          1
#define ALERT          2
#define BOUNDARY       4
#define INCREMENT      8
#define ETA_RELAX      16
#define FLUX_ADJUST    32
#define ETA_ADPT       64

/* Horizontal diffusion/viscosity methods */
#define LAPLACIAN      1
#define SIMPLE         2
#define PRE794         4
#define STANDARD       8
#define ROAM_WIN       16

/* General velocity flags */
#define U1_F1          1
#define U2_F1          2
#define U1_F2          4
#define U2_F2          8
#define U1_F3          16
#define U2_F3          32
#define U1_F4          64
#define U2_F4          128
#define NANF           256

/* Active ALERT tide removal falgs */
#define CSR_R  1     /* Tide removal using csr tide model */
#define MEAN_R 2     /* Tide removal using mean elevation */

/* Tracer boundary scaling flags */
#define TRSC_SUM       1
#define TRSC_PCT       2
#define TRSC_NUM       4
#define TRSC_TRA       8
#define TRSC_DEN       32
#define TRSC_CPY       64
#define TRSC_TRN       128

/* Grid code : TYPE_INCREMENT_ORIENTATION */
#define RECT_DXY_ROT      1
#define GRECT_DXY_ROT     2
#define GRECT_DLATLON_ROT 4
#define GRECT_DLATLON_FP  8
#define POLAR_DARC_ROT    16
#define ELLIPTIC_DTAU_ROT 32
#define NUMERICAL         64

/* Transport modes */
#define SP_EXACT          0x0001
#define SP_TINT           0x0002
#define XYZ_TINT          0x0004
#define GLOBAL            0x0008
#define MONOTONIC         0x0010
#define EXACT             0x0020
#define SUBSET            0x0040
#define INEXACT           0x0080
#define EQUALZ            0x0100
#define UNEQUALZ          0x0200
#define SP_ORIGIN         0x0400
#define SET_BDRY          0x0800
#define SET_BDRY_ETA      0x1000
#define OBC_ADJUST        0x2000
#define WEIGHTED          0x4000
#define DIAGNOSE          0x8000
#define SP_CHECK          0x10000
#define BOTMLEFT          0x100000
#define TOPRIGHT          0x200000
#define DIAGNOSE_BGC      0x400000
#define LOCAL             0x800000
#define TR_CHECK          0x1000000
#define SET_AIJ           0x2000000
#define MONGLOB           0x4000000
#define LOCALER           0x8000000
#define SP_FFSL           0x10000000
#define DO_SWR            0x20000000
#define DO_OBC            0x40000000

#define REINIT            0x0001
#define NOINIT            0x0002
#define CONS_W            0x0004
#define CONS_ETA          0x0008
#define CONS_WS           0x0010
#define CONS_PSS          0x0020
#define CONS_NOF          0x0040
#define CONS_SUB          0x0080
#define CONS_NGR          0x0100
#define CONS_MRG          0x0200

/* Waves */
#define BOT_STR           0x0002
#define TAN_RAD           0x0004
#define ORBITAL           0x0008
#define VERTMIX           0x0010
#define WAVE_FOR          0x0020
#define JONES             0x0040
#define WOM               0x0080
#define BVM               0x0100
#define STOKES            0x0200
#define STOKES_MIX        0x0400
#define STOKES_DRIFT      0x0800
#define STOKES_WIND       0x1000
#define SPECTRAL          0x2000

/* Specific humidity data */
#define WETBULB           2
#define DEWPOINT          4
#define RELHUM            8

/* Water types */
#define TYPE_I    1
#define TYPE_IA   2
#define TYPE_IB   3
#define TYPE_II   4
#define TYPE_III  5

/* MOM conversion modes */
#define CDL   1
#define CDF   2
#define REV   4
#define WDATA 8
#define METRIC 16
#define TEMSAL 32
#define RDGRID 64

/* ROMS grid options */
#define ROMS_GRID_2D    0x01
#define ROMS_GRID_BATHY 0x02
#define ROMS_GRID_MASK  0x04

/* Backwards compatibility */
#define V794    0x000001
#define V1246   0x000002
#define V1283   0x000004
#define V1562   0x000008
#define V1598   0x000010
#define V1652   0x000020
#define V1670   0x000040
#define V1957   0x000080
#define V4201   0x000100
#define V5342   0x000200
#define V5895   0x000400

/* Seasons */
#define DAILY    -1
#define MONTHLY  -2
#define SEASONAL -4
#define YEARLY   -8
#define SUMMER 2
#define AUTUMN 4
#define WINTER 8
#define SPRING 16

/* Directional flags */
#define XDIR      2
#define YDIR      4
#define ZDIR      8

/* Tracer tags */
#define SWR_INVERSE       0x0001
#define DO_SWR_INVERSE    0x0002
#define NO_SWR_INVERSE    0x0004
#define DA_ANL            0x0008
#define DA_OBS_VAL        0x0010
#define DA_OBS_ERR        0x0020
#define SC_SC             0x0040
#define SC_ST             0x0080
#define SC_PC             0x0100
#define SC_PT             0x0200
#define DE_TR2            0x0400
#define DE_TR3            0x0800
#define V_HP              0x1000

/* Dump tags */
#define DF_ETA            0x0001
#define DF_V2D            0x0002
#define DF_MPK            0x0004
#define DF_BARO           0x0008
#define DF_WRITESED       0x0010

/* DA flags */
#define NO_DA             1
#define DO_DA             2
#define FCST_DA           4

/* Blending */
#define U12D    2
#define U13D    4
#define U22D    8
#define U23D    16
#define U1_VH   32
#define U2_VH   64
#define BLM1    128
#define BLM2    256
#define BLR1    512
#define BLR2    1028

/* Relaxation time constant type */
#define RLX_CONS  0x0002
#define RLX_FILE  0x0004
#define RLX_ADPT  0x0008
#define RLX_LINR  0x0010
#define RLX_EXP   0x0020
#define RLX_TIME  0x0040
#define RLX_DEP   0x0080
#define RLX_EDEP  0x0100
#define RLX_CDEP  0x0200
#define RLX_GRD   0x0400

/* Regions */
#define RG_NONE    0
#define RG_START   0x0001
#define RG_WRITE   0x0002
#define RG_TRANS   0x0004
#define RG_COUPD   0x0008
#define RG_BDRY    0x0010
#define RG_AREA    0x0020
#define RG_BEGIN   0x0040
#define RG_VERR    0x0080
#define RG_GLOB    0x0100
#define RG_PSS     0x0200
#define RG_SED     0x0400
#define RG_ECO     0x0800
#define RG_ALLT    0x1000

/* Crash recovery resets */
#define RS_ALL     0x0001
#define RS_VH      0x0002
#define RS_RESTART 0x0004
#define RS_RESUME  0x0008
#define RS_RESET   0x0010
#define RS_WINSET  0x0020
#define RS_FAIL    0x0040
#define RS_OBCSET  0x0080
#define RS_DFSET   0x0100
#define RS_VHSET   0x0200
#define RS_TSSET   0x0400
#define RS_PSSSET  0x0800

/* Process exclusion */
#define EX_TRAN    0x0001
#define EX_BGC     0x0002
#define EX_SED     0x0004
#define EX_WAVE    0x0008
#define EX_TRST    0x0010

/* Partitioning types */
#define STRIPE_E1  0x0001
#define STRIPE_E2  0x0002
#define BLOCK_E1   0x0004
#define BLOCK_E2   0x0008
#define WIN_EXP    0x0010

/* Point source/sinks */
# define PSS_AW    0x0001
# define PSS_VW    0x0002

/* Tracer filtering */
# define TRF_FILL    0x0001
# define TRF_SMOO    0x0002
# define TRF_SHAP    0x0004
# define TRF_MEDI    0x0008
# define TRF_SHUM    0x0010
# define TRF_FILL2D  0x0020
# define TRF_FILL3D  0x0040
# define TRF_FILLSED 0x0080

/* Decorrelation length scale */
# define DEC_ETA    0x0001
# define DEC_U1     0x0002
# define DEC_U2     0x0004
# define DEC_TRA    0x0008

/* Short wave radiation attenuation */
# define SWR_2D     0x0001
# define SWR_3D     0x0002
# define SWR_SPLIT  0x0004

/* Generic interface */
# define GINT_OK   0x0001
# define GINT_QUIT 0x0002
# define GINT_WARN 0x0004

/* Private tracer data */
# define  PTR_SED  0x001
# define  PTR_BGC  0x002
# define  PTR_WAV  0x004
# define  PTR_TRS  0x008

/* Sediments */
#define   LIB_DO    0x001
#define   LIB_WRITE 0x002

/* Velocity initialisation */
#define VINIT_FLOW 1
#define VINIT_GEO  2

/* Time series comparison metrics */
#define TS_NONE    0x0000
#define TS_DIFF    0x0001
#define TS_MEAN    0x0002
#define TS_RMSE    0x0004
#define TS_CAT     0x0008
#define TS_HK      0x0010
#define TS_TS      0x0020
#define TS_H       0x0040
#define TS_F       0x0080
#define TS_CLOS    0x0100
#define TS_PRED    0x0200
#define TS_GLIDER  0x0400

/* Degree heating diagnostic */
#define DHW_NOAA   1
#define DHW_RT     2

/* Misc */
#define INV_BARO 8
#define GHRSST   0x0100

typedef char cstring[MAXSTRLEN];
