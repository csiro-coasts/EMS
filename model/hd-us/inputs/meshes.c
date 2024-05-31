/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/inputs/meshes.c
 *  
 *  Description:
 *  Routine to create grid meshes.
 *  from a file.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: meshes.c 7546 2024-05-10 05:24:35Z her127 $
 *
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hd.h"
/* JIGSAW grid generation */
#ifdef HAVE_JIGSAWLIB
#include "lib_jigsaw.h"
#endif

#define DEG2RAD(d) ((d)*M_PI/180.0)
#define RAD2DEG(r) ((r)*180.0/M_PI)

#define RADIUS 6370997.0
#define ECC 0.0

#define GM_W    0x0001
#define GM_E    0x0002
#define GM_N    0x0004
#define GM_S    0x0008
#define GM_NW   0x0010
#define GM_NE   0x0020
#define GM_B    0x0040
#define GM_F    0x0080
#define GM_L    0x0100
#define GM_R    0x0200
#define GM_U    0x0400
#define GM_D    0x0800

/* Stereogrphic projection */
typedef enum {
  ST3PROJ_FWD,
  ST3PROJ_INV
} ST3PROJ;

/*------------------------------------------------------------------*/
/* Valid bathymetry netCDF dimension names                          */
static char *bathy_dims[4][6] = {
  {"botz", "i_centre", "j_centre", "x_centre", "y_centre", "standard"},
  {"height", "lon", "lat", "lon", "lat", "nc_bathy"},
  {"height", "longitude", "latitude", "longitue", "latitude", "nc_bathy"},
  {NULL, NULL, NULL, NULL, NULL, NULL}
};

int iswetc(unsigned long flag);
int isdryc(unsigned long flag);
void addpoint(parameters_t *params, int i, int j, int *n, point *pin, int **mask);
double triarea(double ax, double ay, double bx, double by, double cx, double cy);
void reorder_pin(int np, point *pin, int ns, int *sin, int nh, int *hin);
int prm_skip_to_end_of_tok(FILE * fp, char *key, char *token, char *ret);
int SortCornersClockwise(double x0, double y0, double x1, double y1, double rx, double ry);
int find_mesh_vertex(int c, double x, double y, double **xloc, double **yloc, int *mask, int ns2, int *npe);
void order_c(double *a, double *b, int *c1, int *c2);
void sort_circle(delaunay *d, int *vedge, int nedge, int dir);
void sort_circle_g(double *x, double *y, int nedge, int dir);
int find_mindex(double slon, double slat, double *lon, double *lat, int nsl, double *d);
void mesh_init_OBC(parameters_t *params, mesh_t *mesh);
void tria_ball_2d(double *bb, double *p1, double *p2, double *p3);
void tria_edge_2d(double *bb, double *p1, double *p2);
void tria_ortho_2d(double *bb, double *p1, double *p2, double *p3);
void tria_ortho_edge_2d(double *bb, double *p1, double *p2);
void tria_com(double *bb, double *p1, double *p2, double *p3);
void remove_duplicates(int ns2, double **x, double **y, mesh_t *mesh);
void mesh_expand(parameters_t *params, double *bathy, double **xc, double **yc);
int mesh_expand_do(geometry_t *window, double *u, int *vec, int *ni, int *filla);
void xy_to_d(delaunay *d, int np, double *x, double *y);
void circen(double *p1, double *p2, double *p3);
double is_obtuse(double *p0, double *p1, double *p2);
void init_J_jig(jigsaw_jig_t *J_jig);
double coast_dist(msh_t *msh, double xloc, double yloc);
double point_dist(int npoints, point *p, double xloc, double yloc);
double poly_dist(int npoly, poly_t **pl, double hmin, double xloc, double yloc, int *mask);
double bathyset(double b, double bmin, double bmax, double hmin,
		double hmax, double expf);
static void st_transform(jigsaw_msh_t *J_msh, ST3PROJ kind,
			 double xmid, double ymid);
static void stereo3_fwd(double xpos, double ypos, double xmid, double ymid,
			double *xnew, double *ynew, double *scal);
static void stereo3_inv(double xpos, double ypos, double xmid, double ymid,
			double *xnew, double *ynew, double *scal);
static void convert_jigsaw_grid_to_mesh(jigsaw_msh_t *imsh,
					jigsaw_msh_t *omsh);
void inv_2x2(int la, double *aa, int lx, double *xx, double *da);
void inv_3x3(int la, double *aa, int lx, double *xx, double *da);
double det_2x2(int la, double *aa);
double det_3x3(int la, double *aa);
void tria_ball_2s(double *bb, double *p1, double *p2, double *p3, double  rr);
void tria_edge_2s(double *bb, double *p1, double *p2, double  rr);
void tria_ortho_2s(double *bb, double *p1, double *p2, double *p3, double  rr);
void tria_ortho_edge_2s(double *bb, double *p1, double *p2, double  rr);
void tria_norm_vect_3d(double *nv, double *pa, double *pb, double *pc);
void edge_circ_ball_3d (double *bb, double *p1,	double *p2);
void tria_circ_ball_3d(double *bb, double *p1, double *p2, double *p3);
void edge_orth_ball_3d(double *bb, double *p1, double *p2);
void tria_orth_ball_3d(double *bb, double *p1, double *p2, double *p3);
void lonlat_to_xyz(double *pp, double *ee, double rr);
void xyz_to_lonlat(double *ee, double *pp, double rr);
void tri_cen(int centref, int geogf, double *p1, double *p2, double *p3, double *po);
void edge_cen(int centref, int geogf, double *v1, double *v2, double *po);
int in_tri(double *po, double *p0, double *p1, double *p2);
void add_quad_grido(parameters_t *params, char *iname);
void add_quad_grid(parameters_t *params, char *iname, int mode);
void reorder_mesh(parameters_t *params, mesh_t *mesh, double *bathy);
int make_ellipsen(double as, double bs, double *xe, double *ye,
		 double xc, double yc, double x0, double y0,
		 double rot, int mode);
void create_bounded_mesh(int npts,     /* Number of perimeter points */
			 double *x,    /* Perimeter x coordinates    */
			 double *y,    /* Perimeter y coordinates    */
			 double rmin,  /* Min. resolution in m       */
			 double rmax,  /* Max. resolution in m       */
			 jigsaw_msh_t *J_mesh,
			 int filef);

/*-------------------------------------------------------------------*/
/* Compute the geographic metrics on the sphere using a false pole   */
/*-------------------------------------------------------------------*/
void geog_false_pole_coord(double **x,  /* where to store grid x values */
                           double **y,  /* where to store grid y values */
                           double **h1, /* where to store h1 metric values
                                         */
                           double **h2, /* where to store h2 metric values
                                         */
                           double **a1, /* where to store a1 angle values */
                           double **a2, /* where to store a2 angle values */
                           long int nce1, /* number of cells in e1
                                             direction */
                           long int nce2, /* number of cells in e2
                                             direction */
                           double x00,  /* x origin offset */
                           double y00,  /* y origin offset */
                           double flon, /* False longitude */
                           double flat, /* False latitude */
                           double xinc, /* cell size in x direction */
                           double yinc  /* cell size in y direction */
  )
{
  long i, j;
  double fx00, fy00, xval, yval;
  double rflon = DEG2RAD(flon);
  double rflat = DEG2RAD(flat);

  geod_fwd_spherical_rot(DEG2RAD(x00), DEG2RAD(y00), rflon, rflat, &fx00,
                         &fy00);
  xinc = DEG2RAD(xinc);
  yinc = DEG2RAD(yinc);

  for (j = 0; j < nce2 + 1; j++) {
    yval = j * yinc;
    for (i = 0; i < nce1 + 1; i++) {
      xval = i * xinc;
      geod_inv_spherical_rot(fx00 + xval, fy00 + yval, rflon, rflat,
                             &x[j][i], &y[j][i]);
      x[j][i] = RAD2DEG(x[j][i]);
      y[j][i] = RAD2DEG(y[j][i]);
    }
  }

  /* Calculate h1 and h2 numerically */
  grid_get_geog_metrics(x, y, nce1, nce2, h1, h2);

  /* calculate a1, a2 numerically */
  grid_get_geog_angle(x, y, nce1, nce2, a1, a2);

  /* Check ranges; must be 0 to 360 */
  for (j = 0; j < nce2 + 1; j++) {
    for (i = 0; i < nce1 + 1; i++) {
      while (x[j][i] < 0.0 || x[j][i] > 360.0) {
	if (x[j][i] < 0.0) x[j][i] += 360.0;
	if (x[j][i] > 360.0) x[j][i] -= 360.0;
      }
    }
  }
}

/* END geog_false_pole_coord()                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Compute the geographic metrics on the sphere by computing the     */
/* latitude and longitude by dead reckoning. We use the Sodano's     */
/* direct formulation to compute the latitude/longitude given a      */
/* range and bearing.                                                */
/*-------------------------------------------------------------------*/
void geog_dreckon_coord(double **x, /* where to store grid x values */
                        double **y, /* where to store grid y values */
                        double **h1,  /* where to store h1 metric values */
                        double **h2,  /* where to store h2 metric values */
                        double **a1,  /* where to store a1 angle values */
                        double **a2,  /* where to store a2 angle values */
                        long int nce1,  /* number of cells in e1 direction
                                         */
                        long int nce2,  /* number of cells in e2 direction
                                         */
                        double x00, /* x origin offset */
                        double y00, /* y origin offset */
                        double rotn,  /* Rotation */
                        double xinc,  /* cell size in x direction */
                        double yinc /* cell size in y direction */
  )
{
  long i, j;
  double rx00 = DEG2RAD(x00);
  double ry00 = DEG2RAD(y00);
  double xaz = M_PI / 2 - DEG2RAD(rotn);
  double yaz = xaz - M_PI / 2;

#if 0
  for (j = 0; j < nce2 + 1; j++) {
    double slat;                /* Start latitude of line */
    double slon;                /* Start longitude of line */
    geod_fwd_sodanos(rx00, ry00, yaz, j * yinc, RADIUS, ECC, &slon, &slat);
    for (i = 0; i < nce1 + 1; i++) {
      geod_fwd_sodanos(slon, slat, xaz, i * xinc, RADIUS, ECC, &x[j][i],
                       &y[j][i]);
      x[j][i] = RAD2DEG(x[j][i]);
      y[j][i] = RAD2DEG(y[j][i]);
    }
  }
#else
  x[0][0] = rx00;
  y[0][0] = ry00;
  for (j = 0; j < nce2 + 1; j++) {
    if (j > 0)
      geod_fwd_sodanos(x[j - 1][0], y[j - 1][0], yaz, yinc, RADIUS, ECC,
                       &x[j][0], &y[j][0]);
    for (i = 1; i < nce1 + 1; i++)
      geod_fwd_sodanos(x[j][i - 1], y[j][i - 1], xaz, xinc, RADIUS, ECC,
                       &x[j][i], &y[j][i]);
  }

  for (j = 0; j < nce2 + 1; j++) {
    for (i = 0; i < nce1 + 1; i++) {
      x[j][i] = RAD2DEG(x[j][i]);
      y[j][i] = RAD2DEG(y[j][i]);
    }
  }
#endif

  /* Calculate h1 and h2 numerically */
  grid_get_geog_metrics(x, y, nce1, nce2, h1, h2);

  /* calculate a1, a2 numerically */
  grid_get_geog_angle(x, y, nce1, nce2, a1, a2);
}

/* END geog_dreckon_coord()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Compute the metrics on a plane for a delaunay grid.               */
/*-------------------------------------------------------------------*/
int delaunay_rect_coord(double *x,   /* where to store grid x values */
			double *y,   /* where to store grid y values */
			long int nce1,  /* number of cells in e1 direction
					 */
			long int nce2,  /* number of cells in e2 direction
					 */
			double x00,  /* x origin offset */
			double y00,  /* y origin offset */
			double rotn, /* Rotation */
			double xinc  /* cell size in x direction */
			)
{
  long i, j, n;
  double xval, yval;
  double sinth;
  double costh;
  double **ax, **ay;
  double yinc = (0.5 * xinc) * tan(DEG2RAD(60.0));
  sinth = sin(DEG2RAD(rotn));
  costh = cos(DEG2RAD(rotn));
  ax = d_alloc_2d(nce1+1, nce2+1);
  ay = d_alloc_2d(nce1+1, nce2+1);
  int iend;

  ax[0][0] = x00;
  ay[0][0] = y00;
  for (j = 0; j < nce2 + 1; j++) {
    double oset;
    if (j%2 == 1) {
      iend = nce1;
      oset = 0.5 * xinc;
    } else {
      iend = nce1 + 1;
      oset = 0.0;
    }
    yval = j * yinc;
    for (i = 0; i < iend; i++) {
      xval = i * xinc;
      ax[j][i] = x00 + oset + xval * costh - yval * sinth;
      ay[j][i] = y00 + xval * sinth + yval * costh;
    }
  }

  /* Count the cells */
  n = 0;
  for (j = 0; j < nce2 + 1; j++) {
    for (i = 0; i < nce1 + 1; i++) {
      if (j%2 == 0) {
	n++;	
      } else {
	if (i < nce1) {
	  n++;	
	}
      }
    }
  }
  n = 0;
  for (j = 0; j < nce2 + 1; j++) {
    for (i = 0; i < nce1 + 1; i++) {
      if (j%2 == 0) {
	x[n] = ax[j][i];
	y[n] = ay[j][i];
	n++;	
      } else {
	if (i < nce1) {
	  x[n] = ax[j][i];
	  y[n] = ay[j][i];
	  n++;	
	}
      }
    }
  }
  d_free_2d(ax);
  d_free_2d(ay);
  return(n);
}

/* END delaunay_rect_coord()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Compute the geographic metrics on the sphere for a delaunay grid  */
/* by computing the latitude and longitude by dead reckoning. We use */
/* the Sodano's direct formulation to compute the latitude/longitude */
/* given a range and bearing.                                        */
/*-------------------------------------------------------------------*/
int delaunay_dreckon_coord(double *x,   /* where to store grid x values */
			   double *y,   /* where to store grid y values */
			   long int nce1,  /* number of cells in e1 direction
					    */
			   long int nce2,  /* number of cells in e2 direction
					    */
			   double x00,  /* x origin offset */
			   double y00,  /* y origin offset */
			   double rotn, /* Rotation */
			   double xinc  /* cell size in x direction */
			   )
{
  long i, j, n;
  double rx00 = DEG2RAD(x00);
  double ry00 = DEG2RAD(y00);
  double xaz = M_PI / 2 - DEG2RAD(rotn);
  double yaz = xaz - M_PI / 2;
  double yinc = (0.5 * xinc) * tan(DEG2RAD(60.0));
  int iend;
  double **ax, **ay;
  int filef = 1;

  ax = d_alloc_2d(nce1+1, nce2+1);
  ay = d_alloc_2d(nce1+1, nce2+1);

  ax[0][0] = rx00;
  ay[0][0] = ry00;
  for (j = 0; j < nce2 + 1; j++) {
    double oset;
    if (j == 0) {
      iend = nce1 + 1;
      oset = 0;
    } else if (j%2 == 1) {
      iend = nce1;
      oset = 0.5 * (ax[j - 1][1] - ax[j - 1][0]);
    } else {
      iend = nce1 + 1;
      oset = -oset;
    }
    if (j > 0) {
      geod_fwd_sodanos(oset + ax[j - 1][0], ay[j - 1][0], yaz, yinc, RADIUS, ECC,
                       &ax[j][0], &ay[j][0]);
    }
    for (i = 1; i < iend; i++)
      geod_fwd_sodanos(ax[j][i - 1], ay[j][i - 1], xaz, xinc, RADIUS, ECC,
                       &ax[j][i], &ay[j][i]);
  }

  /* Count the cells */
  n = 0;
  for (j = 0; j < nce2 + 1; j++) {
    for (i = 0; i < nce1 + 1; i++) {
      if (j%2 == 0) {
	n++;	
      } else {
	if (i < nce1) {
	  n++;	
	}
      }
    }
  }

  n = 0;
  for (j = 0; j < nce2 + 1; j++) {
    for (i = 0; i < nce1 + 1; i++) {
      if (j%2 == 0) {
	x[n] = RAD2DEG(ax[j][i]);
	y[n] = RAD2DEG(ay[j][i]);
	n++;	
      } else {
	if (i < nce1) {
	  x[n] = RAD2DEG(ax[j][i]);
	  y[n] = RAD2DEG(ay[j][i]);
	  n++;	
	}
      }
    }
  }

  if (filef) {
    FILE *op;
    op = fopen("ddg.txt", "w");
    for (i = 0; i < n; i++)
      fprintf(op, "%f %f\n", x[i], y[i]);
    fclose(op);
  }

  d_free_2d(ax);
  d_free_2d(ay);

  return(n);
}

/* END void delaunay_dreckon_coord()                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Finds a false pole                                                */
/*-------------------------------------------------------------------*/
void f_pole(double ilat, double ilon, double ang, double *olat,
            double *olon)
{
  double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0, az = 0.0;
  double qdist = (2.0 * M_PI * RADIUS) / 4.0;

  x1 = DEG2RAD(ilon);
  y1 = DEG2RAD(ilat);
  az = DEG2RAD(ang);

  geod_fwd_sodanos(x1, y1, az, qdist, RADIUS, ECC, &x2, &y2);
  *olat = RAD2DEG(x2);
  *olon = RAD2DEG(y2);
}

/* END f_pole()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Converts a curvilinear mesh to a triangulation using grid corners */
/* for even j and the u2 location for odd j.                         */
/*-------------------------------------------------------------------*/
void convert_quad_mesh(parameters_t *params)
{
  FILE *fp;
  int c, i, j, n, m;
  int **mask, **map;
  double **bathy;
  int nce1 = params->nce1;
  int nce2 = params->nce2;
  int nfe1 = params->nfe1;
  int nfe2 = params->nfe2;
  int nz = params->nz;
  unsigned long ***flg, **flag;
  int nbath;
  double *x, *y, *b;
  double **xin, **yin;
  GRID_SPECS *gs = NULL;
  char i_rule[MAXSTRLEN];
  point *pin;
  int *sin = NULL;
  int np, ns = 0;
  int verbose = 0;  /* Print grid info                               */
  int bverbose = 0; /* Print bathymetry info                         */
  int doseg = 1;    /* Create segments around the grid perimeter for the triangulation */
  int landf = 0;    /* Include land cells in the conversion          */
  int filef = 1;    /* Print grid to file                            */
  int limf = 1;     /* Limit bathymetry to bmin and bmax             */
  int dopoints = 1;
  int ***tri;
  double area, na;
  int ip, im, jp, jm;
  double bmean;

  filef = (params->meshinfo) ? 1 : 0;

  /*-----------------------------------------------------------------*/
  /* Interpolation method. Options:                                  */
  /* linear, cubic, nn_sibson, nn_non_sibson, average                */
  if (params->runmode & ROAM) {
    strcpy(i_rule, "linear");
    landf = 1;
  } else 
    strcpy(i_rule, "linear");
  /*strcpy(i_rule, "nn_non_sibson");*/

  /*-----------------------------------------------------------------*/
  /* Set the (i,j) indexed bathymetry array                          */
  bathy = d_alloc_2d(nce1, nce2);
  n = 0;
  for (j = 0; j < nce2; j++)
    for (i = 0; i < nce1; i++) {
      bathy[j][i] = -params->bathy[n++];
    }
  if (params->runmode & DUMP) {
    for (n = 1; n <= params->nland; n++) {
      i = params->lande1[n];
      j = params->lande2[n];
      bathy[j][i] = LANDCELL;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Make the flags array                                            */
  flg = (unsigned long ***)l_alloc_3d(nce1 + 1, nce2 + 1, nz);
  params->flag = (unsigned long **)l_alloc_2d(nce1 + 1, nce2 + 1);
  make_flags(params, flg, bathy, params->layers, nce1, nce2, nz);	     
  u1_flags(flg, nce1, nce2, nz);
  u2_flags(flg, nce1, nce2, nz);
  if (params->sigma)
    sigma_flags(nfe1, nfe2, nz-1, flg);
  for (j = 0; j < nfe2; j++)
    for (i = 0; i < nfe1; i++) {
      params->flag[j][i] = flg[nz-1][j][i];
      if (landf && params->flag[j][i] & SOLID) params->flag[j][i] &= ~SOLID;
    }
  for (i = 0; i < nfe1; i++) params->flag[nce2][i] |= OUTSIDE;
  for (j = 0; j < nfe2; j++) params->flag[j][nce1] |= OUTSIDE;
  l_free_3d((long ***)flg);
  flag = params->flag;

  /*-----------------------------------------------------------------*/
  /* Set the wet bathymetry vector (to interpolate from)             */
  nbath = n = 0;
  for (j = 0; j < nce2; j++) {
    for (i = 0; i < nce1; i++) {
      if (iswetc(params->flag[j][i])) nbath++;
    }
  }
  x = d_alloc_1d(nbath);
  y = d_alloc_1d(nbath);
  b = d_alloc_1d(nbath);
  bmean = 0.0;
  for (j = 0; j < nce2; j++) {
    for (i = 0; i < nce1; i++) {
      if (iswetc(params->flag[j][i])) {
      x[n] = params->x[j*2][i*2];
      y[n] = params->y[j*2][i*2];
      b[n] = bathy[j][i];
      bmean += b[n];
      n++;
      }
    }
  }
  if (n) bmean /= (double)n;

  /*-----------------------------------------------------------------*/
  /* Double the resolution for hex conversions                       */
  if (params->uscf & US_HEX) {
    nce1 *= 2;
    nce2 *= 2;
    nfe1 = nce1 + 1;
    nfe2 = nce2 + 1;
    xin = d_alloc_2d(2*nce1+1, 2*nce2+1);
    yin = d_alloc_2d(2*nce1+1, 2*nce2+1);

    /* Create a flag array at double the resolution                  */
    flg = (unsigned long ***)l_alloc_3d(nfe1, nfe2, 1);
    for (j = 0; j < nfe2; j++)
      for (i = 0; i < nfe1; i++) {
	flg[0][j][i] = params->flag[j/2][i/2];
      }
    for (i = 0; i < nce1; i++) flg[0][nce2][i] |= OUTSIDE;
    for (j = 0; j < nce2; j++) flg[0][j][nce1] |= OUTSIDE;
    flag = flg[0];

    /* Create grid metric arrays at double the resolution            */
    for (j = 0; j < 2*nce2+1; j++) {
      for (i = 0; i < 2*nce1+1; i++) {
	xin[j][i] = NOTVALID;
	yin[j][i] = NOTVALID;
      }
    }
    for (j = 0; j < nfe2; j++) {
      for (i = 0; i < nfe1; i++) {
	xin[2*j][2*i] = params->x[j][i];
	yin[2*j][2*i] = params->y[j][i];
	xin[2*j][2*i] = params->x[j][i];
	yin[2*j][2*i] = params->y[j][i];

      }
    }
    for (j = 0; j < 2*nce2+1; j++) {
      for (i = 0; i < 2*nce1+1; i++) {
	if (j%2==0) {
	  if (xin[j][i] == NOTVALID) xin[j][i] = 0.5 * (xin[j][i-1] + xin[j][i+1]);
	  if (yin[j][i] == NOTVALID) yin[j][i] = 0.5 * (yin[j][i-1] + yin[j][i+1]);
	} else {
	  if (xin[j][i] == NOTVALID) xin[j][i] = 0.5 * (xin[j-1][i] + xin[j+1][i]);
	  if (yin[j][i] == NOTVALID) yin[j][i] = 0.5 * (yin[j-1][i] + yin[j+1][i]);
	}
      }
    }
  } else {
    xin = params->x;
    yin = params->y;
  }

  /*-----------------------------------------------------------------*/
  /* Set the vertices to the boundaries of the grid.                 */ 
  mask = i_alloc_2d(nfe1, nfe2);
  map = i_alloc_2d(nfe1, nfe2);
  for (j = 0; j < nfe2; j++)
    for (i = 0; i < nfe1; i++) {
      mask[j][i] = 0;
    }

  /* Count the boundary vertices.                                    */
  np = 0;
  for (j = 0; j < nce2; j++) {
    for (i = 0; i < nce1; i++) {
      if (iswetc(flag[j][i])) {
	if (j%2 == 0) {
	  if (!(mask[j][i] & GM_W)) {
	    mask[j][i] |= GM_W;
	    np++;
	  }
	  if (!(mask[j][i+1] & GM_W)) {
	    mask[j][i+1] |= GM_W;
	    mask[j][i] |= GM_E;
	    np++;
	  }
	  if (!(mask[j+1][i] & GM_S)) {
	    mask[j+1][i] |= GM_S;
	    mask[j][i] |= GM_N;
	    np++;
	  }
	} else {
	  if (!(mask[j][i] & GM_S)) {
	    mask[j][i] |= GM_S;
	    np++;
	  }
	  if (!(mask[j+1][i] & GM_W)) {
	    mask[j+1][i] |= GM_W;
	    mask[j][i] |= GM_W;
	    np++;
	  }
	  if (!(mask[j+1][i+1] & GM_W)) {
	    mask[j+1][i+1] |= GM_W;
	    mask[j+1][i] |= GM_E;
	    np++;
	  }
	}
      }
    }
  }
  for (j = 0; j < nfe2; j++)
    for (i = 0; i < nfe1; i++) {
      mask[j][i] = 0;
      map[j][i] = 0;
    }

  /*-----------------------------------------------------------------*/
  /* Allocate the points                                             */
  n = m = 0;
  pin = malloc(np * sizeof(point));
  for (j = 0; j < nce2; j++) {
    for (i = 0; i < nce1; i++) {
      if (iswetc(flag[j][i])) {
	if (j%2 == 0) {
	  if (!(mask[j][i] & GM_W)) {
	    mask[j][i] |= GM_W;
	    pin[n].x = xin[j*2][i*2];
	    pin[n].y = yin[j*2][i*2];
	    if (verbose) printf("%d e(%d %d) : %f %f\n",n, i, j, pin[n].x, pin[n].y);
	    map[j][i] = n;
	    n++;
	  }
	  if (!(mask[j][i+1] & GM_W)) {
	    mask[j][i+1] |= GM_W;
	    mask[j][i] |= GM_E;
	    pin[n].x = xin[j*2][(i+1)*2];
	    pin[n].y = yin[j*2][(i+1)*2];
	    if (verbose) printf("%d e+(%d %d): %f %f\n",n, i+1, j, pin[n].x, pin[n].y);
	    map[j][i+1] = n;
	    n++;
	  }
	  if (!(mask[j+1][i] & GM_S)) {
	    mask[j+1][i] |= GM_S;
	    mask[j][i] |= GM_N;
	    pin[n].x = xin[(j+1)*2][i*2+1];
	    pin[n].y = yin[(j+1)*2][i*2+1];
	    if (verbose) printf("%d o+(%d %d): %f %f\n",n, i, j+1, pin[n].x, pin[n].y);
	    map[j+1][i] = n;
	    n++;
	  }
	} else {
	  if (!(mask[j][i] & GM_S)) {
	    mask[j][i] |= GM_S;
	    pin[n].x = xin[j*2][i*2+1];
	    pin[n].y = yin[j*2][i*2+1];
	    if (verbose) printf("%d o(%d %d) : %f %f\n",n, i, j, pin[n].x, pin[n].y);
	    map[j][i] = n;
	    n++;
	  }
	  if (!(mask[j+1][i] & GM_W)) {
	    mask[j+1][i] |= GM_W;
	    mask[j][i] |= GM_NW;
	    pin[n].x = xin[(j+1)*2][i*2];
	    pin[n].y = yin[(j+1)*2][i*2];
	    if (verbose) printf("%d e+(%d %d): %f %f\n",n, i, j+1, pin[n].x, pin[n].y);
	    map[j+1][i] = n;
	    n++;
	  }
	  if (!(mask[j+1][i+1] & GM_W)) {
	    mask[j+1][i+1] |= GM_W;
	    mask[j+1][i] |= GM_E;
	    mask[j][i] |= GM_NE;
	    pin[n].x = xin[(j+1)*2][(i+1)*2];
	    pin[n].y = yin[(j+1)*2][(i+1)*2];
	    if (verbose) printf("%d e++(%d %d): %f %f\n",n, i+1, j+1, pin[n].x, pin[n].y);
	    map[j+1][i+1] = n;
	    n++;
	  }
	}
      }
    }
  }

  /* Check for duplicates                                            */
  for (i = 0; i < np; i++) {
    double x = pin[i].x;
    double y = pin[i].y;
    for (n = 0; n < np; n++) {
      if (verbose && n != i && x == pin[n].x && y == pin[n].y)
	printf("Duplicate n=%d, i=%d\n", n, i);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Count the segments                                              */
  if (doseg) {
    ns = 0;
    area = na = 0.0;
    tri = i_alloc_3d(nce1, nce2, 3);
    for (j = 0; j < nce2; j++) {
      for (i = 0; i < nce1; i++) {
	if (iswetc(flag[j][i])) {
	  jp = j + 1;
	  jm = (j == 0) ? nce2 : j - 1;
	  ip = i + 1;
	  im = (i == 0) ? nce1 : i - 1;
	  if (j%2 == 0) {
	    tri[0][j][i] = map[j][i];
	    tri[1][j][i] = map[j][ip];
	    tri[2][j][i] = map[jp][i];
	    area += triarea(pin[map[j][i]].x, pin[map[j][i]].y,
			    pin[map[j][ip]].x, pin[map[j][ip]].y,
			    pin[map[jp][i]].x, pin[map[jp][i]].y);
	    na += 1.0;
	    if (isdryc(flag[jm][i])) {
	      mask[j][i] |= (GM_B|GM_U);
	      ns++;
	    }
	    if (isdryc(flag[jp][i]) && iswetc(flag[j][ip])) {
	      mask[j][i] |= (GM_F|GM_U);
	      ns++;
	    }
	  } else {
	    tri[2][j][i] = map[j][i];
	    tri[0][j][i] = map[jp][i];
	    tri[1][j][i] = map[jp][ip];
	    area += triarea(pin[map[jp][i]].x, pin[map[jp][i]].y,
			    pin[map[jp][ip]].x, pin[map[jp][ip]].y,
			    pin[map[j][i]].x, pin[map[j][i]].y);
	    na += 1.0;
	    if (isdryc(flag[jp][i])) {
	      mask[j][i] |= (GM_F|GM_D);
	      ns++;
	    }
	    if (isdryc(flag[jm][i]) && iswetc(flag[j][ip])) {
	      mask[j][i] |= (GM_B|GM_D);
	      ns++;
	    }
	  }
	  if (isdryc(flag[j][im])) {
	    mask[j][i] |= GM_L;
	    if (j%2 == 0)
	      mask[j][i] |= GM_U;
	    else
	      mask[j][i] |= GM_D;
	    ns++;
	  } 
	  if (isdryc(flag[j][ip])) {
	    mask[j][i] |= GM_R;
	    if (j%2 == 0) 
	      mask[j][i] |= GM_U;
	    else
	      mask[j][i] |= GM_D;
	    ns++;
	  }
	}
      }
    }
    if (verbose) {
      printf("Segments = %d\n",ns);
      for (j = 0; j < nce2; j++) {
	for (i = 0; i < nce1; i++) {
	  if (mask[j][i]) {
	    printf("(%d %d) %x : ", i,j,mask[j][i]);
	    if (mask[j][i] & GM_U) printf("up ");
	    if (mask[j][i] & GM_D) printf("down ");
	    if (mask[j][i] & GM_L) printf("left ");
	    if (mask[j][i] & GM_R) printf("right ");
	    if (mask[j][i] & GM_F) printf("front ");
	    if (mask[j][i] & GM_B) printf("back ");
	    printf("\n");
	    /*printf(" : (%d %d %d)\n",tri[0][j][i], tri[1][j][i], tri[2][j][i]);*/
	  }
	}
      }
    }
    sin = i_alloc_1d(2*ns);
    n = 0;
    for (j = 0; j < nce2; j++) {
      for (i = 0; i < nce1; i++) {
	if (iswetc(flag[j][i])) {
	  if (mask[j][i] & GM_F) {
	    if (mask[j][i] & GM_D) {
	      sin[n++] = tri[0][j][i];
	      sin[n++] = tri[1][j][i];
	      if (verbose) printf("segment %d = %d %d (FD %d,%d)\n",n-2, sin[n-2], sin[n-1], i, j);
	    }
	    if (mask[j][i] & GM_U) {
	      sin[n++] = tri[2][j][i];
	      sin[n++] = tri[2][j][i+1];
	      if (verbose) printf("segment %d = %d %d (FU %d,%d)\n",n-2, sin[n-2], sin[n-1], i, j);
	    }
	  }
	  if (mask[j][i] & GM_B) {
	    if (mask[j][i] & GM_U) {
	      sin[n++] = tri[0][j][i];
	      sin[n++] = tri[1][j][i];
	      if (verbose) printf("segment %d = %d %d (BU %d,%d)\n",n-2, sin[n-2], sin[n-1], i, j);
	    }
	    if (mask[j][i] & GM_D) {
	      sin[n++] = tri[2][j][i];
	      sin[n++] = tri[2][j][i+1];
	      if (verbose) printf("segment %d = %d %d (BD %d,%d)\n",n-2, sin[n-2], sin[n-1], i, j);
	    }
	  }
	  if (mask[j][i] & GM_L) {
	    sin[n++] = tri[0][j][i];
	    sin[n++] = tri[2][j][i];
	    if (verbose) printf("segment %d = %d %d (LUD %d,%d)\n",n-2, sin[n-2], sin[n-1], i, j);
	  }
	  if (mask[j][i] & GM_R) {
	    sin[n++] = tri[1][j][i];
	    sin[n++] = tri[2][j][i];
	    if (verbose) printf("segment %d = %d %d (RUD %d,%d)\n",n-2, sin[n-2], sin[n-1], i, j);
	  }
	}
      }
    }

    if (n != 2*ns) hd_quit("quad_mesh: inconsistent segment numbers (%d != %d)\n",n, ns);
    area /= na;
  }

  /*-----------------------------------------------------------------*/
  /* Write to file                                                   */
  if (filef) {
    if ((fp = fopen("points.txt", "w")) != NULL) {
      for (i = 0; i < np; i++)
	fprintf(fp, "%f %f\n",pin[i].x, pin[i].y);
      fclose(fp);
    }
    if (doseg) {
      if ((fp = fopen("segment.txt", "w")) != NULL) {
	for (n = 0; n < 2*ns; n+=2) {
	  fprintf(fp,"%f %f\n",pin[sin[n]].x, pin[sin[n]].y); 
	  fprintf(fp,"%f %f\n",pin[sin[n+1]].x, pin[sin[n+1]].y); 
	  fprintf(fp, "NaN NaN\n");
	}
	fclose(fp);
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Create the triangulation                                        */
  reorder_pin(np, pin, ns, sin, 0, NULL);
  if (dopoints) 
    area = 0.0;
  else
    np = ns;
  params->d = delaunay_voronoi_build(np, pin, ns, sin, 0, NULL, area, NULL);
  
  if (params->uscf & US_TRI)
    convert_tri_mesh(params, params->d);
  if (params->uscf & US_HEX)
    convert_hex_mesh(params, params->d, 1);

  /*-----------------------------------------------------------------*/
  /* Reset flags                                                     */
  strcpy(params->gridtype, "UNSTRUCTURED");
  params->gridcode = UNSTRUCTURED;
  params->us_type |= US_G;
  if (params->us_type & US_IJ)
    params->us_type &= ~US_IJ;
  if (params->us_type & US_RS)
    params->us_type &= ~US_RS;
  if (params->us_type & US_WS)
    params->us_type &= ~US_RS;
  params->us_type |= (US_RUS|US_WUS);

  /*-----------------------------------------------------------------*/
  /* Interpolate the bathymetry                                      */
  gs = grid_interp_init(x, y, b, nbath, i_rule);
  d_free_1d(params->bathy);
  params->bathy = d_alloc_1d(params->ns2+1);
  for (c = 1; c <= params->ns2; c++) {
    params->bathy[c] = grid_interp_on_point(gs, params->x[c][0], params->y[c][0]);
    if (limf) {
      if (params->bmax && fabs(params->bathy[c]) > params->bmax) 
	params->bathy[c] = -params->bmax;
      if (params->bmin && fabs(params->bathy[c]) < params->bmin) 
	params->bathy[c] = -params->bmin;
      if (params->bmin && params->bathy[c] > params->bmin)
	params->bathy[c] = LANDCELL;
    }
    if (isnan(params->bathy[c])) params->bathy[c] = bmean;
    if (bverbose) printf("%d %f : %f %f\n",c, params->bathy[c], params->x[c][0], params->y[c][0]);
  }
  grid_specs_destroy(gs);
  params->nvals = params->ns2-1;

  d_free_1d(x);
  d_free_1d(y);
  d_free_1d(b);
  free((point *)pin);
  d_free_2d(bathy);
  i_free_2d(mask);
  i_free_2d(map);
  if (doseg) {
    i_free_1d(sin);
    i_free_3d(tri);
  }
}

/* END convert_quad_mesh()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Re-orders the points list so that the segments occupy the first   */
/* ns points, followed by nh holes points, followed by points in the */
/* interior of the domain. This allows the segments and holes to be  */
/* retained when refining the mesh points in the interior (i.e.      */
/* chaning the points list from (ns + nh) onwards).                  */
/*-------------------------------------------------------------------*/
void reorder_pin(int np, point *pin, int ns, int *sin, int nh, int *hin)
{
  int i, n, m, mp;
  point *npin= malloc(np * sizeof(point));
  int *mask = i_alloc_1d(np);
  int *map = i_alloc_1d(np);
  int verbose = 0;

  /* Reorder the segments and make a map of old to new points.       */
  mp = 0;
  memset(mask, 0, np * sizeof(int));
  if (ns) {
    for (n = 0; n < 2*ns; n++) {
      m = sin[n];
      if (!mask[m]) {
	map[m] = mp;
	npin[mp].x = pin[m].x;
	npin[mp++].y = pin[m].y;
	for (i = n; i < 2*ns; i++) {	
	  if (sin[i] == sin[n]) {
	    mask[sin[i]] = 1;
	  }
	}
      }
    }

    /* Reorder the segments using the map                            */
    for (n = 0; n < 2*ns; n++) {
      sin[n] = map[sin[n]];
    }
  }

  /* Reorder the holes and make a map of old to new points.          */
  if (nh) {
    for (n = 0; n < 2*nh; n++) {
      m = hin[n];
      if (!mask[m]) {
	npin[mp].x = pin[m].x;
	npin[mp++].y = pin[m].y;
	for (i = n; i < 2*nh; i++) {
	  if (hin[i] == hin[n]) {
	    mask[hin[i]] = 1;
	  }
	}
      }
    }
    for (n = 0; n < 2*nh; n++) {
      hin[n] = map[hin[n]];
    }
  }

  /* Fill the remaining points                                       */
  for (n = 0; n < np; n++) {
    if (!mask[n]) {
      npin[mp].x = pin[n].x;
      npin[mp++].y = pin[n].y;
    }
  }
  for (n = 0; n < np; n++) {
    pin[n].x = npin[n].x;
    pin[n].y = npin[n].y;
  }

  /* Print if required                                               */
  if (verbose) {
    for (n = 0; n < np; n++) 
      printf("point %d %f %f\n",n, npin[n].x, npin[n].y);
    for (n = 0; n < 2*ns; n+=2) 
      printf("segment %d : %d %d\n",n,sin[n],sin[n+1]);
    for (n = 0; n < 2*nh; n+=2) 
      printf("hole %d : %d %d\n",n,hin[n],hin[n+1]);
  }
  free((point *)npin);
  i_free_1d(mask);
  i_free_1d(map);
}

/* END reorder_pin()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to create a Delaunay and Voronoi mesh given a set of      */
/* points.                                                           */
/*-------------------------------------------------------------------*/
delaunay *create_tri_mesh(int np, point *pin, int ns, int *sin, int nh, double *hin, char *code)
{
  FILE *fp;
  delaunay *d;
  int i, j, n, m;
  point *npin;
  double area, areaf = 0.5;
  int verbose = 0;
  int iterations = 0;
  int filef = 0;

  /*-----------------------------------------------------------------*/
  /* Write to file                                                   */
  if (filef) {
    if ((fp = fopen("trigrid.txt", "w")) == NULL) filef = 0;
    for (i = 0; i < np; i++) {
      /*
      n = (int)(pin[i].x * 1e4);
      pin[i].x = (double)n / 1e4;
      m = (int)(pin[i].y * 1e4);
      pin[i].y = (double)m / 1e4;
      */
      fprintf(fp, "%f %f\n",pin[i].x, pin[i].y);
    }
    fclose(fp);
  }

  /*-----------------------------------------------------------------*/
  /* Do the triangulation                                            */
  if (sin == NULL)
    d = delaunay_voronoi_build(np, pin, 0, NULL, 0, NULL, 0, code);
  else
    d = delaunay_voronoi_build(np, pin, ns, sin, nh, hin, 0, code);

  if (verbose) {
    printf("Triangulation has %d points\n", d->npoints);
    printf("Voronoi has %d points\n", d->nvoints);
  }

  npin = malloc(d->npoints * sizeof(point));

  /*-----------------------------------------------------------------*/
  /* Iterate to improve the triangulation                            */
  for (m = 0; m < iterations; m++) {
    /* Get the new set of points                                     */
    npin = malloc(d->npoints * sizeof(point));
    for (i = 0; i < d->npoints; ++i) {
      point* p = &d->points[i];
      npin[i].x = p->x;
      npin[i].y = p->y;
      if (verbose) printf("Delaunay points %d: %f %f\n", i, p->x, p->y);
    }

    /* Get the mean triangle area. A fraction of this is used to     */
    /* refine the next iteration of the mesh.                        */
    area = 0;
    for (n = 0; n < d->ntriangles; n++) {
      triangle* t = &d->triangles[n];
      triangle_neighbours* nb = &d->neighbours[n];
      circle* c = &d->circles[n];
      area += triarea(d->points[t->vids[0]].x, d->points[t->vids[0]].y,
		      d->points[t->vids[1]].x, d->points[t->vids[1]].y,
		      d->points[t->vids[2]].x, d->points[t->vids[2]].y);
      if (verbose) printf("triangle %d %d %d %d\n",n, t->vids[0], t->vids[1], t->vids[2]);
    }
    area /= (double)d->ntriangles;

    for (n = 0, j = 0; n < d->nedges-1; n++) {
      if (verbose) printf("edge %d %d %d\n",n, d->edges[j], d->edges[j+1]);
      j += 2;
    }
    area *= areaf;
    if (verbose) printf("Triangle area = %f\n",area);

    /* Print Voronoi points and edges if required                    */
    for (i = 0; i < d->nvoints; ++i) {
      point* p = &d->voints[i];
      if (verbose) printf("Voronoi points %d: %f %f\n", i, p->x, p->y);
    }

    for (n = 0, j = 0; n < d->nvdges; n++) {
      if (verbose) printf("Voronoi edge %d %d %d\n",n, d->vdges[j], d->vdges[j+1]);
      j += 2;
    }

    /* Create the next iteration of the mesh                         */
    free((delaunay *)d);

    if (sin == NULL)
      d = delaunay_voronoi_build(np, pin, 0, NULL, 0, NULL, area, code);
    else 
      d = delaunay_voronoi_build(np, pin, ns, sin, nh, hin, area, code);
    if (verbose) {
      printf("Triangulation has %d points\n", d->npoints);
      printf("Voronoi has %d points\n\n", d->nvoints);
    }
  }
  free((point *)npin);
  return(d);
}

/* END create_tri_mesh()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to create a Delaunay and Voronoi mesh with ordered        */
/* indexing using a delaunay datastructure as input. The output mesh */
/* is stored in params->x and params->y.                             */
/*-------------------------------------------------------------------*/
void convert_hex_mesh(parameters_t *params, delaunay *d, int mode)
{
  FILE *fp, *cf, *ef, *vf, *bf;
  char oname[MAXSTRLEN], key[MAXSTRLEN], buf[MAXSTRLEN];
  int **mask;
  int *vmask, **vedge, *nvedge, *bmask;
  point *pin;
  int npe;            /* Voronoi tesselation */
  int i, j, n, m, ne, nc;
  int donpe = 0;     /* 1 = only include cells with six sides */
  int dodel = 0;
  int verbose = 0;
  int filef = 1;
  int nnh = 0;
  double dist, dmin;
  double xn, yn, xs, ys;
  double mlon, mlat;
  int nskip = 0;
  int sortdir = 1; /* Sort verices; 1=clockwise, -1=anticlockwise    */
  int debug = -1;  /* Debugging information for points index         */
  int isdedge = 0; /* Debugging location is an edge (not a point)    */
  int edqu = 0;    /* 1 = exit if closed cells can't be made         */
  int regadj = 0;  /* Adjustment for regular hex meshes              */
  int nn;
  int vd;

  /*-----------------------------------------------------------------*/
  /* Allocate                                                        */
  filef = (params->meshinfo) ? 1 : 0;
  nvedge = i_alloc_1d(d->npoints);
  bmask = i_alloc_1d(d->npoints);
  memset(bmask, 0, d->npoints * sizeof(int));
  /*vmask = i_alloc_1d(d->npoints);*/
  if (params->gridcode & (GTRI_DXY_ROT|TRI_DXY_ROT)) donpe = 0;
  donpe = 1;
  if(mode==0) {
    donpe = 0;
    mode = 1;
  }
  if (params->us_type & US_HEX) regadj = 1;

  /*-----------------------------------------------------------------*/
  /* Get the centre and edges of each Voronoi cell. The centre is a  */
  /* point in the list of points. First store all the number of all  */
  /* d->edges emanting from points in vedge[]. The corresponding     */
  /* edge in d->vdges is the edge of the Voronoi cell whose centre   */
  /* is the point.                                                   */

  for (n = 0; n < d->npoints; n++) {
    nvedge[n] = 0;
    for (m = 0, j = 0;  m < d->nedges; m++) {
      if (n == d->edges[j] || n == d->edges[j+1]) {
	nvedge[n]++;
      }
      j += 2;
    }
    if (n == debug) printf("centre %d[%f %f] = %d vedges\n",n,
			   d->points[n].x, d->points[n].y, nvedge[n]);
  }
  npe = 6;

  if (donpe) {
    npe = 0;
    for (n = 0; n < d->npoints; n++) {
      if (nvedge[n] > npe) npe = nvedge[n];
    }

    /* If the maximum npe lies on a boundary that requires an extra  */
    /* edge, then nvedge[] is incremented by 1, so add 1 to npe for  */
    /* safety to account for this case.                              */
    npe++;
  }

  /*-----------------------------------------------------------------*/

  /* Convert debug edges to a point                                  */
  if (debug >= 0 && isdedge) {
    debug = d->edges[2*debug];
    printf("DEBUG = %d[%f %f]\n", debug, d->points[debug].x, d->points[debug].y);
  }

  /* Save the edges i (0:npe) emanating from point n in vedge[n][i]  */
  vedge = i_alloc_2d(npe, d->npoints);
  for (n = 0; n < d->npoints; n++) {
    nvedge[n] = 0;
    for (m = 0, j = 0;  m < d->nedges; m++) {
      if (n == d->edges[j] || n == d->edges[j+1]) {
	vedge[n][nvedge[n]] = m;
	if (n == debug) {
	  printf("vedge%d of edge%d = %d[%f %f] to %d[%f %f]\n", nvedge[n], m, d->vdges[2*m],
		 d->voints[d->vdges[2*m]].x, 
		 d->voints[d->vdges[2*m]].y, d->vdges[2*m+1],
		 d->voints[d->vdges[2*m+1]].x, 
		 d->voints[d->vdges[2*m+1]].y);
	}
	nvedge[n]++;
      }
      j += 2;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Adjust regular meshes by removing triangles where the           */
  /* triangulation has produced 'flipped' triangles with valence     */
  /* less than 6. These can occur on boundary edges using HEX_DX.    */
  /* For non-uniform meshes turn this off.                           */
  if (regadj) {
    double md;
    double area = 0;
    for (n = 0; n < d->ntriangles; n++) {
      triangle* t = &d->triangles[n];
      area += triarea(d->points[t->vids[0]].x, d->points[t->vids[0]].y,
		      d->points[t->vids[1]].x, d->points[t->vids[1]].y,
		      d->points[t->vids[2]].x, d->points[t->vids[2]].y);
    }
    area /= (double)d->ntriangles;
    for (n = 0; n < d->ntriangles; n++) {
      triangle* t = &d->triangles[n];
      md = triarea(d->points[t->vids[0]].x, d->points[t->vids[0]].y,
		   d->points[t->vids[1]].x, d->points[t->vids[1]].y,
		   d->points[t->vids[2]].x, d->points[t->vids[2]].y);
      
      if (md < area) {
	for (i = 0; i < 3; i++) {
	  if (d->n_point_triangles[t->vids[i]] <= 6) {
	    nvedge[t->vids[i]] = 0;
	  }
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Rearrange the edges to be continuous. Store the Voronoi edges   */
  /* in mask[].                                                      */
  mask = i_alloc_2d(npe+1, d->npoints);
  for (n = 0; n < d->npoints; n++) {
    int nn, ce, cs, c1, c2;
    int is_closed = 0;
    int cep, csp, cap[npe], cap2[npe], ip, ip2;
    char p1[MAXSTRLEN], p2[MAXSTRLEN];

    if (!nvedge[n]) continue;
    vmask = i_alloc_1d(nvedge[n]+1);
    memset(vmask, 0, (nvedge[n]+1) * sizeof(int));

    /* First edge                                                    */
    i = 0;
    vmask[i] = 1;
    j = 2 * vedge[n][0];
    mask[n][i++] = cs = d->vdges[j];  /* Save the start index        */
    j += 1;
    mask[n][i++] = ce = d->vdges[j];  /* Index to carry forward      */
    if (n == debug) {
      printf("Start index0 = %d [%f %f]\n", cs, d->voints[cs].x, d->voints[cs].y);
      printf("Carry index1 = %d [%f %f]\n", ce, d->voints[ce].x, d->voints[ce].y);
    }
    if (edqu) {
      csp = cs;
      cep = ce;
    }
    /* Remaining edges                                               */
    for (nn = 0; nn < nvedge[n]; nn++) {
      int found = 0;
      /* Loop over all edges and find an index equal to the one      */
      /* carried forward. If found, make this index the new one to   */
      /* carry forward. This may reside on either end of the         */
      /* Voronoi edge n.                                             */
      for (m = 0; m < nvedge[n]; m++) {
	if (found || vmask[m]) continue;
	j = 2 * vedge[n][m];
	if(n==debug)printf("%d %d %d %d\n",nn,m,d->vdges[j],d->vdges[j+1]);
	if (ce == d->vdges[j]) {
	  if (edqu) {
	    strcpy(p1, "Pass1(a)");
	    ip = i;
	  }
	  vmask[m] = 1;
	  mask[n][i++] = ce  = d->vdges[j + 1];
	  if (n == debug) printf("Pass1(a) index%d = %d[%f %f]\n",
				 i-1, ce, d->voints[ce].x, d->voints[ce].y);
	  found = 1;
	} else if (ce == d->vdges[j + 1]) {
	  if (edqu) {
	    strcpy(p1, "Pass1(b)");
	    ip = i;
	  }
	  vmask[m] = 1;
	  mask[n][i++] = ce = d->vdges[j];
	  if (n == debug) printf("Pass1(b) index%d = %d[%f %f]\n",
				 i-1, ce, d->voints[ce].x, d->voints[ce].y);
	  found = 1;
	}

	/* If the edge to carry forward is equal to the start index, */
	/* (i.e. a dead end is encountered) then a closed cell has   */
	/* been formed. Note that last value of mask contains the    */
	/* same index as that one before the last value (i.e.        */
	/* redundant information).                                   */
	if (ce == cs) is_closed = 1;
      }
      /*if (i == nvedge[n]) break;*/
    }
    if (n == debug && is_closed) printf("Cell %d is closed\n", n);

    /* If the loop is not closed we most likely have a boundary cell */
    /* that we require to close. Increment the number of points in   */
    /* the Voronoi cell (nvedge[n]) by one, and add indices in the   */
    /* opposite direction to that above, using the start index as    */
    /* that to carry forward, until a dead end is found. Note that   */
    /* the indexing in this direction starts from nvedge[n] and      */
    /* decrements.                                                   */
    /* The Voronoi edge to close the loop is that joining the two    */
    /* dead ends.                                                    */
    if (!is_closed) {
      int nve = nvedge[n] - 1;
      int ii = nvedge[n];
      int ni = i;
      nvedge[n]++;
      c1 = ce;
      ce = cs;
      /*
      printf("Cell not closed at %d: (nedges=%d) : %f %f\n", 
	     n, ii, d->points[n].x, d->points[n].y);
      */
      /* Do not include triangular boundary cells                    */
      if (nvedge[n] == 3) nvedge[n] = 0;
      /* Do not include quadrilateral boundary cells                 */
      /*if (nvedge[n] == 4) nvedge[n] = 0;*/

      /*
      if(nvedge[n])printf("Cell not closed at %d: (nedges=%d) : %f %f\n", 
			 n, ii, d->points[n].x, d->points[n].y);
      */

      /* Search for edges in the reverse direction */
      for (nn = nvedge[n]-1; nn >= 0 ; nn--) {
	int found = 0;
	for (m = nve; m >= 0; m--) {
	  if (found || vmask[m]) continue;
	  j = 2 * vedge[n][m];
	  if (ce == d->vdges[j]) {
	    if (edqu) {
	      strcpy(p2, "Pass2(a)");
	      ip2 = ii;
	    }
	    vmask[m] = 1;
	    mask[n][ii--] = ce  = d->vdges[j + 1];
	    if (n == debug) printf("Pass2(a) index%d = %d[%f %f]\n",
				   ni, ce, d->voints[ce].x, d->voints[ce].y);
	    ni++;
	    found = 1;
	  } else if (ce == d->vdges[j + 1]) {
	    if (edqu) {
	      strcpy(p2, "Pass2(b)");
	      ip2 = ii;
	    }
	    vmask[m] = 1;
	    mask[n][ii--] = ce = d->vdges[j];
	    if (n == debug) printf("Pass2(b) index%d = %d[%f %f]\n",
				   ni, ce, d->voints[ce].x, d->voints[ce].y);
	    ni++;
	    found = 1;
	  }
	}
      }
      c2 = ce;
      if(nvedge[n] && nvedge[n] != ni) {
	hd_warn("convert_hex_mesh: Can't create closed cell at %d [%f %f]. Removing cell.\n",n, d->points[n].x, d->points[n].y);
	if (edqu && nskip == 1) {
	  printf("nvedge[%d] = %d\n",n, nvedge[n]);
	  printf("ni = %d\n", ni);
	  for (nn = 0; nn < nvedge[n]; nn++) {
	    j = 2 * vedge[n][nn];
	    ce = d->vdges[j];
	    printf("vertex%d = %d[%f %f]\n",nn, ce, d->voints[ce].x, d->voints[ce].y);
	  }
	  cs = mask[n][0];
	  ce = mask[n][1];
	  printf("Start index0 = %d [%f %f]\n", cs, d->voints[cs].x, d->voints[cs].y);
	  printf("Carry index1 = %d [%f %f]\n", ce, d->voints[ce].x, d->voints[ce].y);
	  printf("Forward #points = 0:%d\n", ip);
	  for (nn = 0; nn < ip; nn++) {
	    ce = mask[n][nn];
	    printf("%s index%d = %d[%f %f]\n",p1, nn, ce, d->voints[ce].x, d->voints[ce].y);
	  }
	  printf("Backward #points = %d:%d\n", nve, ip2);
	  for (nn = nve; nn >= ip2; nn--) {
	    ce = mask[n][nn];
	    printf("%s index%d = %d[%f %f]\n",p2, nn, ce, d->voints[ce].x, d->voints[ce].y);
	  }
	  hd_quit("convert_hex_mesh: Can't create closed cell at %d [%f %f]. Removing cell.\n",n, d->points[n].x, d->points[n].y);
	}
	nvedge[n] = 0;
	nskip++;
      }
      bmask[n] = 1;
    }
    i_free_1d(vmask);
  }

  /*-----------------------------------------------------------------*/
  /* Remove edges to infinity if required                            */
  if (mode) {
    for (n = 0; n < d->npoints; n++) {
      /* 
      if (nvedge[n] != 6) {
	nvedge[n] = 0;
	break;
      }
      */
      for (m = 0; m < nvedge[n]; m++) {

	if (mask[n][m] == -1) {
	  nvedge[n] = 0;
	  break;
	}

	if (params->gridcode & (GTRI_DXY_ROT|GTRI_DLATLON_FP)) {
	  if (fabs(d->voints[mask[n][m]].x) > 360.0 || 
	      fabs(d->voints[mask[n][m]].y) > 360.0) {
	    nvedge[n] = 0;
	    break;
	  }
	  if (isnan(d->voints[mask[n][m]].x) || 
	      isnan(d->voints[mask[n][m]].y)) {
	    nvedge[n] = 0;
	    break;
	  }
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Reorder clockwise starting from most SW vertex. The vector      */
  /* containing the Voronoi edges (vedge[]) is overwritten with the  */
  /* edge numbers stored in mask[].                                  */
  if (mode) {
  for (n = 0; n < d->npoints; n++) {
    int nn, dir;
    double xsw = 1e10;
    double ysw = 1e10;
    double xn;
    double eps, epsm;
    double xr, yr;
    double x0, y0, x1, y1;
    int dof = 0;
    if (!nvedge[n]) continue;

    vmask = i_alloc_1d(nvedge[n]);
    memset(vmask, 0, nvedge[n] * sizeof(int));

    /* Sort in clockwise (or anti-clockwise) order. NOTE: All cells  */
    /* must have vertices ordered in the same sense (clockwise or    */
    /* anti-clockwise) or the preprocessor will fail (e2v will be in */
    /* the wrong direction if one cell adacent to an edge is in one  */
    /* sense and the other cell is in the opposite sense).           */
    sort_circle(d, mask[n], nvedge[n], sortdir);

    /* Get the most westerly                                         */
    /* Get the tolerance based on grid size                          */
    eps = 1e-10;
    epsm = HUGE;
    for (m = 0; m < nvedge[n]; m++) {
      double xn;
      xn = d->voints[mask[n][m]].x;
      for (nn = 0; nn < nvedge[n]; nn++) {
	double x;
	x = d->voints[mask[n][nn]].x;
	if (fabs(x-xn) > eps) eps = fabs(x-xn);
	if (fabs(x-xn) < epsm) epsm = fabs(x-xn);
      }
    }
    eps = max(0.25 * eps, epsm);

    /* Mask all vertices with distance differences > tolerance       */
    for (m = 0; m < nvedge[n]; m++) {
      double xn;
      /*if (mask[n][m] == -1) continue;*/
      xn = d->voints[mask[n][m]].x;
      for (nn = 0; nn < nvedge[n]; nn++) {
	double x;
	/*if (mask[n][nn] == -1) continue;*/
	x = d->voints[mask[n][nn]].x;
	if (x < xn && fabs(x-xn) > eps) vmask[m] = 1;
      }
    }

    /* Get the most southerly of un-masked vertices                  */
    for (m = 0; m < nvedge[n]; m++) {
      double y;
      /*if (mask[n][m] == -1) continue;*/
      y = d->voints[mask[n][m]].y;
      if (!vmask[m] && y < ysw) {
	ysw = y;
	i = m;
      }
    }

    /* Re-order                                                      */
    for (m = 0; m < nvedge[n]; m++) {
      vedge[n][m] = mask[n][i];
      if (n == debug) printf("Re-ordered: %d(%d)=>%d(%d)\n",i, mask[n][i], m, vedge[n][m]);
      i = (i == nvedge[n] - 1) ? 0 : i + 1;
    }

    if (dof) {
    xr = yr = 0.0;
    for (m = 0; m < nvedge[n]; m++) {
      xr += d->voints[vedge[n][m]].x;
      yr += d->voints[vedge[n][m]].y;
    }
    xr /= (double)nvedge[n];
    yr /= (double)nvedge[n];
    x0 = d->voints[vedge[n][0]].x;
    y0 = d->voints[vedge[n][0]].y;
    x1 = d->voints[vedge[n][1]].x;
    y1 = d->voints[vedge[n][1]].y;

    /* Set clockwise (in order of increasing latitude from vertex 0  */
    /* to 1). If latitudes are equal, then vertex 1 must be to the   */
    /* east of vertex 0 (vertex 0 is the SW corner) and order should */
    /* be reversed.                                                  */
    dir = SortCornersClockwise(x0, y0, x1, y1, xr, yr);

    if (dir == -1) {
      for (m = 0; m < nvedge[n]; m++)
	mask[n][m] = vedge[n][m];
      i = nvedge[n] - 1;
      for (m = 1; m < nvedge[n]; m++) {
	vedge[n][m] = mask[n][i--];
      }
    }
    }
    i_free_1d(vmask);
  }
  } else {
    /* Ideally truncate the point to the nearest boundary */
    for (n = 0; n < d->npoints; n++) {
      for (m = 0; m < nvedge[n]; m++) {
	if (mask[n][m] == -1) {
	  i = m;
	  while (mask[n][i] == -1)
	    i = (i == nvedge[n] - 1) ? 0 : i + 1;
	  vedge[n][m] = mask[n][i];
	} else
	  vedge[n][m] = mask[n][m];
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Write to file                                                   */
  if (filef) {
    mesh_ofile_name(params, key);
    sprintf(buf,"%s_us.us", key);
    if ((fp = fopen(buf, "w")) == NULL) filef = 0;
    sprintf(buf,"%s_c.txt", key);
    if ((cf = fopen(buf, "w")) == NULL) filef = 0;
    sprintf(buf,"%s_e.txt", key);
    if ((ef = fopen(buf, "w")) == NULL) filef = 0;
    sprintf(buf,"%s_v.txt", key);
    if ((vf = fopen(buf, "w")) == NULL) filef = 0;
    sprintf(buf,"%s_b.txt", key);
    if ((bf = fopen(buf, "w")) == NULL) filef = 0;
  }

  m = 0;
  for (n = 0; n < d->npoints; n++)
    if (nvedge[n]) m++;
  if (filef) {
    fprintf(fp, "Mesh2 unstructured v1.0\n");
    fprintf(fp, "nMaxMesh2_face_nodes %d\n", npe);
    fprintf(fp, "nMesh2_face          %d\n", m);
    fprintf(fp, "Mesh2_topology\n");
  }
  params->npe = npe;
  params->ns2 = 0;
  if (!donpe) {
    for (n = 0, i = 0; n < d->npoints; n++)
      if (nvedge[n] == npe)
	params->ns2++;
  } else {
    for (n = 0, i = 0; n < d->npoints; n++)
      if (nvedge[n])
	params->ns2++;
  }
  if (params->x != NULL) d_free_2d(params->x);
  if (params->y != NULL) d_free_2d(params->y);
  if (params->npe2 != NULL) i_free_1d(params->npe2);
  params->npe2 = i_alloc_1d(params->ns2+1);
  params->x = d_alloc_2d(params->npe+1, params->ns2+1);
  params->y = d_alloc_2d(params->npe+1, params->ns2+1);
  nc = max(d->npoints, d->nvoints);

  mask =  i_alloc_2d(nc, 2);
  memset(mask[0], 0, nc * sizeof(int));
  memset(mask[1], 0, nc * sizeof(int));
  nc = ne = 0;
  for (n = 0, i = 1; n < d->npoints; n++) {
    if (nvedge[n]) {
      double xm = 0.0, ym = 0.0;
      if (nvedge[n] != 6) nnh++;
      if (!donpe && nvedge[n] != npe) continue;
      params->npe2[i] = nvedge[n];
      if (verbose) printf("Voronoi cell %d %f %f\n", i, d->points[n].x, d->points[n].y);
      if (filef) {
	fprintf(fp, "%d %f %f\n", i, d->points[n].x, d->points[n].y);
      }
      if (debug == n) {
	printf("centre Voronoi cell %d[%f %f]\n", i, d->points[n].x, d->points[n].y);
	vd = i;
      }
      params->x[i][0] = d->points[n].x;
      params->y[i][0] = d->points[n].y;

      mask[0][n] = 1;
      nc++;
      for (m = 0; m < nvedge[n]; m++) {
	if (filef) {
	  fprintf(fp, "%f %f\n", d->voints[vedge[n][m]].x, d->voints[vedge[n][m]].y);
	  fprintf(ef, "%f %f\n", d->voints[vedge[n][m]].x, d->voints[vedge[n][m]].y);
	}
	if (n == debug) printf("vedge=%d[%f %f]\n", vedge[n][m], d->voints[vedge[n][m]].x, d->voints[vedge[n][m]].y);
	if (verbose) printf("  vertex %d=%d : %f %f\n", m, vedge[n][m],
			    d->voints[vedge[n][m]].x, d->voints[vedge[n][m]].y);
	params->x[i][m+1] = d->voints[vedge[n][m]].x;
	params->y[i][m+1] = d->voints[vedge[n][m]].y;
	xm += params->x[i][m+1];
	ym += params->y[i][m+1];
	if (mask[1][vedge[n][m]] == 0) {
	  mask[1][vedge[n][m]] = 1;
	  ne++;
	}
      }
      /* Replace the centre location with the centre of mass for     */
      /* perimeter cells.                                            */
      if (bmask[n]) {
	params->x[i][0] = xm / (double)nvedge[n];
	params->y[i][0] = ym / (double)nvedge[n];
      }

      if (filef) {
	fprintf(cf, "%f %f\n", params->x[i][0], params->y[i][0]);
	fprintf(ef, "%f %f\n", d->voints[vedge[n][0]].x, d->voints[vedge[n][0]].y);
	fprintf(ef, "NaN NaN\n");
	if (bmask[n])
	  fprintf(bf, "%f %f\n", d->points[n].x, d->points[n].y);
      }
      if (n == debug) {
	printf("vedge=%d[%f %f]\n", vedge[n][0], d->voints[vedge[n][0]].x, d->voints[vedge[n][0]].y);
      }
      i++;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set up a delaunay triangulation for xytoi interpolation using   */
  /* delaunay_xytoi().                                               */
  if (dodel) {
    pin = malloc((nc + ne) * sizeof(point));
  
    for (i = 0, j = 0; i < d->npoints; ++i) {
      if (mask[0][i]) {
	pin[j].x = d->points[i].x;
	pin[j++].y = d->points[i].y;
      }
    }
    for (i = 0; i < d->nvoints; ++i) {
      if (mask[1][i]) {
	pin[j].x = d->voints[i].x;
	pin[j++].y = d->voints[i].y;
      }
    }
    i_free_2d(mask);
    free((delaunay *)d);
    /*d = delaunay_voronoi_build(j, pin, 0, NULL, 0, NULL, 0.0, NULL);*/
    d = delaunay_build(j, pin, 0, NULL, 0, NULL);
    if (filef) {
      for (n = 0; n < d->ntriangles; n++) {
	triangle* t = &d->triangles[n];
	fprintf(vf, "%f %f\n", d->points[t->vids[0]].x, d->points[t->vids[0]].y);
	fprintf(vf, "%f %f\n", d->points[t->vids[1]].x, d->points[t->vids[1]].y);
	fprintf(vf, "%f %f\n", d->points[t->vids[2]].x, d->points[t->vids[2]].y);
	fprintf(vf, "%f %f\n", d->points[t->vids[0]].x, d->points[t->vids[0]].y);
	fprintf(vf, "NaN NaN\n");
      }
    }
  }
  params->us_type |= US_CUS;

  if (filef) {
    fclose(fp);
    fclose(cf);
    fclose(ef);
    fclose(vf);
    fclose(bf);
  }
  i_free_2d(mask);
  i_free_1d(nvedge);
  i_free_2d(vedge);
  i_free_1d(bmask);
}

/* END convert_hex_mesh()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Orders a list of coordinates in a clockwise direction, starting   */
/* from the coordinate with the smallest angle to the x-plane for    */
/* anti clockwise sorting (dir=-1) and largest angle for clockwise   */
/* sorting (dir=1).                                                  */
/*-------------------------------------------------------------------*/
void sort_circle(delaunay *d, int *vedge, int nedge, int dir)
{
  double x[nedge];
  double y[nedge];
  double a[nedge];
  double xr = 0.0, yr = 0.0;
  int m, i, j;

  /* Save the locations and index and find the central reference     */
  for(m = 0; m < nedge; m++) {
    x[m] = d->voints[vedge[m]].x;
    y[m] = d->voints[vedge[m]].y;
    xr += x[m];
    yr += y[m];
  }
  xr /= (double)nedge;
  yr /= (double)nedge;

  /* Get the atan2 values                                            */
  for(m = 0; m < nedge; m++) {
    a[m] = atan2(y[m] - yr, x[m] - xr);
    if (a[m] < 0.0) a[m] += 2.0 * PI;
  }

  /* Order the angles                                                */
  for (i = 0; i < nedge; i++)
    for (j = nedge-1; i < j; --j)
      order_c(&a[j-1], &a[j], &vedge[j-1], &vedge[j]);

  if (dir = 1) {
    int index[nedge];
    memcpy(index, vedge, nedge * sizeof(int));
    i = 0;
    for(m = nedge-1; m >= 0; m--) {
      vedge[i++] = index[m];
    }
  }
}


void sort_circle_g(double *x, double *y, int nedge, int dir)
{
  double a[nedge+1];
  int vedge[nedge+1];
  int vmask[nedge+1];
  double xr = x[0], yr = y[0];
  double xc[nedge+1];
  double yc[nedge+1];
  int m, i, j;
  double eps = 1e-10;
  double epsm = HUGE;
  double ysw = 1e10;

  /* Get the atan2 values                                            */
  for(m = 1; m <= nedge; m++) {
    xc[m] = x[m];
    yc[m] = y[m];
    vedge[m] = m;
    vmask[m] = 0;
    a[m] = atan2(y[m] - yr, x[m] - xr);
    if (a[m] < 0.0) a[m] += 2.0 * PI;
  }

  /* Order the angles                                                */
  for (i = 1; i <= nedge; i++)
    for (j = nedge; i < j; --j)
      order_c(&a[j-1], &a[j], &vedge[j-1], &vedge[j]);

  for (m = 1; m <= nedge; m++) {
    x[m] = xc[vedge[m]];
    y[m] = yc[vedge[m]];
  }

  if (dir = 1) {
    for (m = 1; m <= nedge; m++) {
      xc[m] = x[m];
      yc[m] = y[m];
    }
    i = 1;
    for(m = nedge; m >= 1; m--) {
      x[i] = xc[m];
      y[i++] = yc[m];
    }
  }


  /* Get the most westerly                                           */
  /* Get the tolerance based on grid size                            */
  for (m = 1; m <= nedge; m++) {
    double xm = x[m];
    for (j = 1; j <= nedge; j++) {
      double xn = x[j];
      if (fabs(xm-xn) > eps) eps = fabs(xm-xn);
      if (fabs(xm-xn) < epsm) epsm = fabs(xm-xn);
    }
    xc[m] = x[m];
    yc[m] = y[m];
  }
  eps = max(0.25 * eps, epsm);

  /* Mask all vertices with distance differences > tolerance         */
  for (m = 1; m <= nedge; m++) {
    double xm = x[m];
    for (j = 0; j < nedge; j++) {
      double xn = x[j];
      if (xn < xm && fabs(xn-xm) > eps) vmask[m] = 1;
    }
  }

  /* Get the most southerly of un-masked vertices                    */
  for (m = 1; m <= nedge; m++) {
    double ym = y[m];
    if (!vmask[m] && ym < ysw) {
      ysw = ym;
      i = m;
    }
  }

  /* Re-order                                                        */
  for (m = 1; m <= nedge; m++) {
    x[m] = xc[i];
    y[m] = yc[i];
    i = (i == nedge) ? 1 : i + 1;
  }
}

/* END sort_circle()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Orders a list of double values and thier associated index.        */
/*-------------------------------------------------------------------*/
void order_c(double *a, double *b, int *c1, int *c2)
{
  double val;
  int v;
  if (*a > *b) {
    val = *a;
    *a = *b;
    *b = val;
    v = *c1;
    *c1 = *c2;
    *c2 = v;
  }
}

/* END order_c()                                                     */
/*-------------------------------------------------------------------*/

int SortCornersClockwise(double x0, double y0, double x1, double y1, double rx, double ry)
{
  double aTanA, aTanB;

  if (x0 == x1 && y0 < y1) return 1;
  if (x0 == x1 && y0 > y1) return -1;
  if (y0 == y1 && x0 < x1) return -1;
  if (y0 == y1 && x0 > x1) return 1;

  aTanA = atan2(y0 - ry, x0 - rx);
  if (aTanA < 0.0) aTanA += 2.0 * PI;
  aTanB = atan2(y1 - ry, x1 - rx);
  if (aTanB < 0.0) aTanB += 2.0 * PI;

  if (aTanA < aTanB) return -1;
  else if (aTanB < aTanA) return 1;
  return 0;
}

/*-------------------------------------------------------------------*/
/* Routine to create a Delaunay and Voronoi mesh given a set of      */
/* points.                                                           */
/*-------------------------------------------------------------------*/
void convert_tri_mesh(parameters_t *params, delaunay *d)
{
  FILE *fp, *cf, *ef;
  char oname[MAXSTRLEN], key[MAXSTRLEN], buf[MAXSTRLEN];
  point *pin;
  int npe = 3;            /* Delaunay tesselation */
  int i, j, n, m;
  int verbose = 0;
  int filef = 1;

  /*-----------------------------------------------------------------*/
  /* Write to file                                                   */
  filef = (params->meshinfo) ? 1 : 0;
  if (filef) {
    mesh_ofile_name(params, key);
    sprintf(buf,"%s_us.us", key);
    if ((fp = fopen(buf, "w")) == NULL) filef = 0;
    sprintf(buf,"%s_c.txt", key);
    if ((cf = fopen(buf, "w")) == NULL) filef = 0;
    sprintf(buf,"%s_e.txt", key);
    if ((ef = fopen(buf, "w")) == NULL) filef = 0;
  }
  if (filef) {
    fprintf(fp, "Mesh2 unstructured v1.0\n");
    fprintf(fp, "nMaxMesh2_face_nodes %d\n", npe);
    fprintf(fp, "nMesh2_face          %d\n", d->ntriangles);
    fprintf(fp, "Mesh2_topology\n");
  }
  params->npe = npe;
  params->ns2 = d->ntriangles;
  if (params->x != NULL) d_free_2d(params->x);
  if (params->y != NULL) d_free_2d(params->y);
  if (params->npe2 != NULL) i_free_1d(params->npe2);
  params->npe2 = i_alloc_1d(params->ns2+1);
  for (n = 1; n <= params->ns2; n++)
    params->npe2[n] = npe;
  params->x = d_alloc_2d(params->npe+1, params->ns2+1);
  params->y = d_alloc_2d(params->npe+1, params->ns2+1);

  for (n = 0, i = 1; n < d->ntriangles; n++) {
    triangle* t = &d->triangles[n];
    double x[3], y[3], c1, c2, area, ec;
    x[0] = d->points[t->vids[0]].x;
    y[0] = d->points[t->vids[0]].y;
    x[1] = d->points[t->vids[1]].x;
    y[1] = d->points[t->vids[1]].y;
    x[2] = d->points[t->vids[2]].x;
    y[2] = d->points[t->vids[2]].y;
    c1 = (x[0] + x[1] + x[2]) / 3.0;
    c2 = (y[0] + y[1] + y[2]) / 3.0;

    /* Find the eastern-most coordinate                              */
    m = 0;
    for (j = 0; j < 3; j++) {
      if (x[j] < x[m]) m = j;
    }
    /* Reset coordinates clockwise from easternmost (note: triangle  */
    /* coordinates are listed in counter-clockwise sense.            */
    for (j = 0; j < 3; j++) {
      x[j] = d->points[t->vids[m]].x;
      y[j] = d->points[t->vids[m]].y;
      m = (m == 0) ? 2 : m - 1;
    }

    area = triarea(x[0], y[0], x[1], y[1], x[2], y[2]);

    if (verbose) printf("Delaunay cell %d %f %f\n", i, c1, c2);
    if (filef) {
      fprintf(fp, "%d %f %f\n", i, c1, c2);
      fprintf(cf, "%f %f\n", c1, c2);
    }
    params->x[i][0] = c1;
    params->y[i][0] = c2;
    for (m = 0; m < 3; m++) {
      if (filef) {
	fprintf(fp, "%f %f\n", x[m], y[m]);
	fprintf(ef, "%f %f\n", x[m], y[m]);
      }
      if (verbose) printf("  vertex %d : %f %f\n", m, x[m], y[m]);
      params->x[i][m+1] = x[m];
      params->y[i][m+1] = y[m];
    }
    if (filef) {
      fprintf(ef, "%f %f\n", x[0], y[0]);
      fprintf(ef, "NaN NaN\n");
    }
    i++;
  }

  params->us_type |= US_CUS;

  if (filef) {
    fclose(fp);
    fclose(cf);
    fclose(ef);
  }
}

/* END convert_tri_mesh()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Area of a triangle given vertex coordinates using Heron's formula */
/* and Pythagoras.                                                   */
/*-------------------------------------------------------------------*/
double triarea(double ax, double ay, double bx, double by, double cx, double cy)
{
  double a = sqrt((ax - bx) * (ax - bx) + (ay - by) * (ay - by));
  double b = sqrt((ax - cx) * (ax - cx) + (ay - cy) * (ay - cy));
  double c = sqrt((bx - cx) * (bx - cx) + (by - cy) * (by - cy));
  double s = 0.5 * (a + b + c);
  /*double area = 0.5*ax*(by - cy) + bx*(cy - ay) + cx*(ay - by);*/
  double area = sqrt(s * (s-a) * (s-b) * (s-c));
  return(area);
}

/* END triarea()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Converts a mesh given params->x[cc][1:npe] to a list of vertex    */
/* and centre locations (double) and a list of integers pointing to  */
/* these locations to define the polygon edges and centres.          */
/*-------------------------------------------------------------------*/
void convert_mesh_input(parameters_t *params,
			mesh_t *mesh,  /* Mesh structure             */
			double **xc,   /* Cell x coordinates         */
			double **yc,   /* Cell y coordinates         */
			double *bathy, /* Bathymetry                 */
			int wf         /* 1 = write to file          */
			)
{
  FILE *op;
  char key[MAXSTRLEN], buf[MAXSTRLEN];
  int ns2;               /* Number of cells                          */
  int npe;               /* Maximum number of vertices               */
  int *maskv;            /* Vertex mask                              */
  int cc, c, cci, j, jj, i, n, cco;
  double x0, y0, x1, y1, x2, y2;
  double *xv, *yv;
  int iv;
  int oldcode = 0;

  if (!params->meshinfo) wf = 0;
  if (DEBUG("init_m"))
    dlog("init_m", "\nStart mesh conversion\n");

  ns2 = mesh->ns2;
  npe = mesh->mnpe;

  /*-----------------------------------------------------------------*/
  /* Save all unique locations of vertices and make the mapping from */
  /* indices to coordinates.                                         */
  remove_duplicates(ns2, xc, yc, mesh);
  if (DEBUG("init_m"))
    dlog("init_m", "\nIndex - coordinate mappings created OK\n");

  if (oldcode) {
  maskv = i_alloc_1d(ns2 + 1);
  memset(maskv, 0, (ns2 + 1) * sizeof(int));
  mesh->ns = 1;

  for (cc = 1; cc <= ns2; cc++) {
    for (j = 1; j <= mesh->npe[cc]; j++) {
      x1 = xc[cc][j];
      y1 = yc[cc][j];
      i = find_mesh_vertex(cc, x1, y1, xc, yc, maskv, ns2, mesh->npe);
      if (!i) {
	mesh->xloc[mesh->ns] = x1;
	mesh->yloc[mesh->ns++] = y1;
      }
    }
    maskv[cc] = 1;
  }
  i_free_1d(maskv);
  if (DEBUG("init_m"))
    dlog("init_m", "\nUnique vertices found OK\n");

  /*-----------------------------------------------------------------*/
  /* Save the location of the cell centres                           */
  for (cc = 1; cc <= ns2; cc++) {
    x1 = xc[cc][0];
    y1 = yc[cc][0];
    mesh->xloc[mesh->ns] = x1;
    mesh->yloc[mesh->ns++] = y1;
  }
  mesh->ns--;
  if (DEBUG("init_m"))
    dlog("init_m", "\nCell centres found OK\n");

  /*-----------------------------------------------------------------*/
  /* Map the edges for each cell to the vertex location numbers, eg. */
  /* cell cc, edge j has the x vertices xloc[eloc[0][cc][j]] and     */
  /* xloc[eloc[1][cc][j]].                                           */
  for (cc = 1; cc <= ns2; cc++) {
    for (j = 1; j <= mesh->npe[cc]; j++) {
      x1 = xc[cc][j];
      y1 = yc[cc][j];
      jj = (j == mesh->npe[cc]) ? 1 : j + 1;
      x2 = xc[cc][jj];
      y2 = yc[cc][jj];
      for (n = 1; n <= mesh->ns; n++) {
	if (x1 == mesh->xloc[n] && y1 == mesh->yloc[n])
	  mesh->eloc[0][cc][j] = n;
	if (x2 == mesh->xloc[n] && y2 == mesh->yloc[n])
	  mesh->eloc[1][cc][j] = n;
      }
    }
    for (n = 1; n <= mesh->ns; n++) {
      x1 = xc[cc][0];
      y1 = yc[cc][0];
      if (x1 == mesh->xloc[n] && y1 == mesh->yloc[n])
	mesh->eloc[0][cc][0] = mesh->eloc[1][cc][0] = n;
    }
  }
  if (DEBUG("init_m"))
    dlog("init_m", "\nCell vertices mapped OK\n");
  }

  /*-----------------------------------------------------------------*/
  /* Get the neighbour mappings                                      */
  if (DEBUG("init_m"))
    dlog("init_m", "\nMesh neighbours\n");
  neighbour_finder(mesh);

  /*-----------------------------------------------------------------*/
  /* Expand the mesh if required                                     */
  if (DEBUG("init_m"))
    dlog("init_m", "\nMesh expansion\n");
  mesh_expand(params, params->bathy, xc, yc);
  if (params->runmode & (AUTO | DUMP) && !(params->us_type & US_IJ))
    mesh_reduce(params, params->bathy, xc, yc);

  /*-----------------------------------------------------------------*/
  /* Open boundaries. Remap the obc edge number to the vertex        */
  /* location indices.                                               */
  if (DEBUG("init_m"))
    dlog("init_m", "\nMesh OBCs\n");
  if (mesh->nobc) {
    for (n = 0; n < mesh->nobc; n++) {
      /* If mesh indices were read from the parameter file and       */
      /* copied directly to the mesh structure in meshstruct_us()    */
      /* then continue.                                              */ 
      if (mesh->loc[n][0] == NOTVALID) continue;
      for (cc = 1; cc <= mesh->npts[n]; cc++) {
	cco = mesh->loc[n][cc];
	x0 = xc[cco][0];
	y0 = yc[cco][0];
	j = mesh->obc[n][cc][0];
	x1 = xc[cco][j];
	y1 = yc[cco][j];
	j = mesh->obc[n][cc][1];
	x2 = xc[cco][j];
	y2 = yc[cco][j];

	for (cci = 1; cci <= mesh->ns2; cci++) {
	  if (x0 == mesh->xloc[mesh->eloc[0][cci][0]] && y0 == mesh->yloc[mesh->eloc[0][cci][0]])
	    mesh->loc[n][cc] = cci;
	}
	for (cci = 1; cci <= mesh->ns; cci++) {
	  if (x1 == mesh->xloc[cci] && y1 == mesh->yloc[cci])
	    mesh->obc[n][cc][0] = cci;
	  if (x2 == mesh->xloc[cci] && y2 == mesh->yloc[cci])
	    mesh->obc[n][cc][1] = cci;
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Write to file if required. Output file is <INPUT_FILE>_us.txt   */
  if (wf) {
    mesh_ofile_name(params, key);
    sprintf(buf,"%s_us.txt", key);
    if ((op = fopen(buf, "w")) == NULL)
      return;
    fprintf(op, "Mesh2 unstructured   v1.0\n");
    fprintf(op, "nMaxMesh2_face_nodes %d\n", mesh->mnpe);
    fprintf(op, "nMesh2_face_indices  %d\n", mesh->ns);
    fprintf(op, "nMesh2_face          %d\n",mesh->ns2);
    if (mesh->nce1)
      fprintf(op, "NCE1                 %d\n",mesh->nce1);
    if (mesh->nce2)
      fprintf(op, "NCE2                 %d\n",mesh->nce2);
    fprintf(op, "Mesh2_topology\n");

    write_mesh_us(params, bathy, op, 1);

    fclose(op);
  }

  if (DEBUG("init_m"))
    dlog("init_m", "\nMesh parameters\n");
  set_params_mesh(params, mesh, xc, yc);

  if (DEBUG("init_m"))
    dlog("init_m", "\nMesh conversion complete\n");
}

/* END convert_mesh_input()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Initializes the params flags based on information in the mesh     */
/* structure.                                                        */
/*-------------------------------------------------------------------*/
void set_params_mesh(parameters_t *params, mesh_t *mesh, double **x, double **y)
{
  int cc;

  /* Set the grid code                                               */
  if (mesh->mnpe == 3) 
    params->us_type |= US_TRI;
  else if (mesh->mnpe == 4) 
    params->us_type |= US_QUAD;
  else if (mesh->mnpe == 6) 
    params->us_type |= US_HEX;
  else
    params->us_type |= US_MIX;

  /* Get the grid metrics                                            */
  if (params->h1 != NULL) d_free_2d(params->h1);
  if (params->h2 != NULL) d_free_2d(params->h2);
  if (params->a1 != NULL) d_free_2d(params->a1);
  if (params->a2 != NULL) d_free_2d(params->a2);
  params->h1 = d_alloc_2d(mesh->mnpe+1, mesh->ns2+1);
  params->h2 = d_alloc_2d(mesh->mnpe+1, mesh->ns2+1);
  params->a1 = d_alloc_2d(mesh->mnpe+1, mesh->ns2+1);
  params->a2 = d_alloc_2d(mesh->mnpe+1, mesh->ns2+1);

  grid_get_metrics_us(mesh->npe, x, y, mesh->ns2,
		      params->h1, params->h2);
  grid_get_angle_us(mesh->npe, x, y, mesh->ns2,
		    params->a1, params->a2);
  /*params->gridcode = NUMERICAL;*/

  /* Set the time-series file projection                             */
  if (strlen(params->projection) > 0) {
    strcpy(projection, params->projection);
    ts_set_default_proj_type(projection);
  }

  /* If the the projection is specified and is geographic, then      */
  /* compute the metrics and angles for the sphere.                  */
  if (strncasecmp(params->projection,
                  GEOGRAPHIC_TAG, strlen(GEOGRAPHIC_TAG)) == 0) {
    grid_get_geog_metrics_us(mesh->npe, x, y, mesh->ns2,
			     params->h1, params->h2);
    grid_get_geog_angle_us(mesh->npe, x, y, mesh->ns2,
			   params->a1, params->a2);
  }
}

/* END set_params_mesh()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Build a mesh struture for structured orthogonal curvilinear grids */
/*-------------------------------------------------------------------*/
void meshstruct_s(parameters_t *params, geometry_t *geom) 
{
  int cc, c, i, j, cco, n;
  mesh_t *m;
  double **x, **y;
  int mo2;              /* Max. # of cells in each open boundary */
  int filef = 1;

  /* Allocate                                                        */
  filef = (params->meshinfo) ? 1 : 0;
  free_mesh(params->mesh);
  params->mesh = mesh_init(params, geom->b2_t, 4);
  m = params->mesh;
  strcpy(params->gridtype, "STRUCTURED");

  /*-----------------------------------------------------------------*/
  /* Cell centre and vertex locations                                */
  x = d_alloc_2d(m->mnpe+1, m->ns2+1);
  y = d_alloc_2d(m->mnpe+1, m->ns2+1);
  m->nce1 = geom->nce1;
  m->nce2 = geom->nce2;
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    i = geom->s2i[c];
    j = geom->s2j[c];
    x[cc][0] = params->x[j*2+1][i*2+1];
    y[cc][0] = params->y[j*2+1][i*2+1];
    if (params->us_type & US_IJ) {
      m->iloc[cc] = i;
      m->jloc[cc] = j;
    }
    x[cc][1] = params->x[j*2][i*2];
    y[cc][1] = params->y[j*2][i*2];
    x[cc][2] = params->x[(j+1)*2][i*2];
    y[cc][2] = params->y[(j+1)*2][i*2];
    x[cc][3] = params->x[(j+1)*2][(i+1)*2];
    y[cc][3] = params->y[(j+1)*2][(i+1)*2];
    x[cc][4] = params->x[j*2][(i+1)*2];
    y[cc][4] = params->y[j*2][(i+1)*2];
  }

  /*-----------------------------------------------------------------*/
  /* Open boundaries                                                 */
  if (params->nobc) {
    for (n = 0; n < params->nobc; n++) {
      open_bdrys_t *open = geom->open[n];
      convert_obc_list(params, open, n, geom, NULL);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Bathymetry                                                      */
  if (params->bathy) d_free_1d(params->bathy);
  params->bathy = d_alloc_1d(m->ns2+1);
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    params->bathy[cc] = geom->botz[c];
  }

  /*-----------------------------------------------------------------*/
  /* Write to file                                                   */
  /*-----------------------------------------------------------------*/
  /*
  if (filef) {
    FILE *ef;
    char key[MAXSTRLEN], buf[MAXSTRLEN];
    if (endswith(params->oname, ".nc")) {
      n = strlen(params->oname);
      for (i = 0; i < n-3; i++)
	key[i] = params->oname[i];
      key[i] = '\0';
    } else
      strcpy(key, params->oname);
    sprintf(buf,"%s_e.txt", key);
    if ((ef = fopen(buf, "w")) == NULL) filef = 0;
    if (filef) {
      for (cc = 1; cc <= geom->b2_t; cc++) {
	for (n = 1; n <= 4; n++) {
	  fprintf(ef, "%f %f\n",x[cc][n],y[cc][n]);
	}
	fprintf(ef, "%f %f\n",x[cc][1],y[cc][1]);
	fprintf(ef, "NaN Nan\n");
      }
    }
  }
  */

  /*-----------------------------------------------------------------*/
  /* Convert the coordinate information to mesh indices              */
  convert_mesh_input(params, m, x, y, params->bathy, 1);

  /* Write to file from the mesh structure                           */
  if (filef) {
    FILE *ef;
    char key[MAXSTRLEN], buf[MAXSTRLEN];
    mesh_ofile_name(params, key);
    sprintf(buf,"%s_e.txt", key);
    if ((ef = fopen(buf, "w")) == NULL) filef = 0;
    if (filef) {
      for (cc = 1; cc <= m->ns2; cc++) {
	for (n = 1; n <= m->npe[cc]; n++) {
	  fprintf(ef, "%f %f\n",m->xloc[m->eloc[0][cc][n]], m->yloc[m->eloc[0][cc][n]]);
	}
	fprintf(ef, "%f %f\n",m->xloc[m->eloc[0][cc][1]], m->yloc[m->eloc[0][cc][1]]);
	fprintf(ef, "NaN Nan\n");
      }
      fclose(ef);
    }
    sprintf(buf,"%s_c.txt", key);
    if ((ef = fopen(buf, "w")) == NULL) filef = 0;
    if (filef) {
      for (cc = 1; cc <= m->ns2; cc++)
	fprintf(ef, "%f %f\n",m->xloc[m->eloc[0][cc][0]], m->yloc[m->eloc[0][cc][0]]);
      fclose(ef);
    }
    if (m->nobc) {
      sprintf(buf,"%s_b.txt", key);
      if ((ef = fopen(buf, "w")) == NULL) filef = 0;
      if (filef) {
	for (n = 0; n < m->nobc; n++) {
	  for (cc = 1; cc <= m->npts[n]; cc++) {
	    fprintf(ef, "%f %f\n", m->xloc[m->obc[n][cc][0]], m->yloc[m->obc[n][cc][0]]);
	    fprintf(ef, "%f %f\n", m->xloc[m->obc[n][cc][1]], m->yloc[m->obc[n][cc][1]]);
	    fprintf(ef, "NaN NaN\n");
	  }
	}
	fclose(ef);
      }
    }
  }

  d_free_2d(x);
  d_free_2d(y);
}

/* END meshstruct_s()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Build an input mesh struture for unstructured grids               */
/*-------------------------------------------------------------------*/
void meshstruct_us(parameters_t *params)
{
  char key[MAXSTRLEN], buf[MAXSTRLEN];
  int cc, c, i, j, jj, cco, n;
  mesh_t *m;
  int mo2;              /* Max. # of cells in each open boundary */
  double bw, bo;
  double **x, **y;
  int *cmap;
  int verbose = 0;      /* Verbose output                            */
  int filef = 2;        /* Output mesh (1 from params, 2 from mesh)  */
  int testp = 0;        /* Test a closed perimeter exists            */
  int debug = -1;       /* Debug for params->x[debug][n]             */

  if (DEBUG("init_m"))
    dlog("init_m", "\nBegin building mesh structure\n");

  /*-----------------------------------------------------------------*/
  /* If an unstructured mesh configuration is input from file, then  */
  /* the mesh structure is populated in read_mesh_us(), so return.   */
  if (!params->meshinfo) filef = 0;
  if (params->us_type & US_IUS) {
    mesh_init_OBC(params, params->mesh);
    m = params->mesh;
    cmap = i_alloc_1d(params->ns2+1);
    for (cc = 1; cc <= params->ns2; cc++) {
      cmap[cc] = cc;
      if (m->map) cmap[cc] = m->map[cc];
    }
    /* Note: open boundary indices have not been adjusted for mesh   */
    /* reductions at this stage, so use the mapping cmap to account  */
    /* for this when the mesh structure is populated with            */
    /* params->open information.                                     */
    for (n = 0; n < params->nobc; n++) {
      convert_obc_list(params, params->open[n], n, NULL, cmap);
    }
    if (params->runmode & (AUTO | DUMP)) set_bathy_us(params);
    i_free_1d(cmap);

    /* Print any added curvilinear qrids if required                 */
    if (strlen(params->addquad)) {
      char *fields[MAXSTRLEN * MAXNUMARGS];
      n = parseline(params->addquad, fields, MAXNUMARGS);
      if (n == 1 && atoi(fields[0]) > 0) {
	for (i = 0; i < atoi(fields[0]); i++) {
	  sprintf(key, "QUAD%d", i);
	  prm_read_char(params->prmfd, key, buf);
	  add_quad_grid(params, buf, 0);
	  hd_warn("ADD_QUAD %s tested.\n", buf);
	}
      } else {
	for (i = 0; i < n; i++) {
	  add_quad_grid(params, fields[i], 0);
	  hd_warn("ADD_QUAD %s tested.\n", fields[i]);
	}
      }
    }
    return;
  }

  /*-----------------------------------------------------------------*/
  /* Merge the mesh with a curvilinear qrid if required              */
  if (strlen(params->addquad)) {
    char *fields[MAXSTRLEN * MAXNUMARGS];
    n = parseline(params->addquad, fields, MAXNUMARGS);
    if (n == 1 && atoi(fields[0]) > 0) {
      for (i = 0; i < atoi(fields[0]); i++) {
	sprintf(key, "QUAD%d", i);
	prm_read_char(params->prmfd, key, buf);
	add_quad_grid(params, buf, 1);
      }
    } else {
      for (i = 0; i < n; i++) {
	add_quad_grid(params, fields[i], 1);
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Strip out any land cells                                        */
  c = 1;
  cmap = i_alloc_1d(params->ns2+1);
  x = d_alloc_2d(params->npe+1, params->ns2+1);
  y = d_alloc_2d(params->npe+1, params->ns2+1);
  for (cc = 1; cc <= params->ns2; cc++) {
    i = (params->npe2[cc]) ? params->npe2[cc] : params->npe;
    /* Make a copy of the original mesh to get the OBCs below        */
    for (n = 0; n <= i; n++) {
      x[cc][n] = params->x[cc][n];
      y[cc][n] = params->y[cc][n];
    }
    if (params->bathyfill & B_REMOVE) {
      if (params->bathy[cc] != LANDCELL && params->bathy[cc] != NOTVALID) {
	for (n = 0; n <= i; n++) {
	  params->x[c][n] = params->x[cc][n];
	  params->y[c][n] = params->y[cc][n];
	}
	params->bathy[c] = params->bathy[cc];
	params->npe2[c] = params->npe2[cc];
	cmap[cc] = c;
	if (cc == debug) printf("New params->x[i] coordinate after land stripping = %d\n", c);
	c++;
      }
    } else if (params->bathyfill & B_MINIMUM) {
      if (params->bathy[cc] == LANDCELL || params->bathy[cc] == NOTVALID) 
	params->bathy[c] = -params->bmin;
      cmap[cc] = cc;
    }
    /* Bathymetry filling using surrounding average (B_AVERAGE)      */
    /* computed below after neighbour maps have been created.        */
  }
  if (params->bathyfill & B_REMOVE) params->ns2 = c - 1;

  /*-----------------------------------------------------------------*/
  /* Allocate                                                        */
  free_mesh(params->mesh);
  params->mesh = mesh_init(params, params->ns2, 0);
  m = params->mesh;
  if (!(params->us_type & US_IJ))
    strcpy(params->gridtype, "UNSTRUCTURED");

  if (DEBUG("init_m"))
    dlog("init_m", "\nMesh structure allocated OK\n");

  /*-----------------------------------------------------------------*/
  /* Open boundaries                                                 */
  if (m->nobc) {
    for (n = 0; n < params->nobc; n++) {
      open_bdrys_t *open = params->open[n];
      /* For mesh reductions, OBCs cell centres are referenced to    */
      /* the indices and need to be remapped. The edge locations are */
      /* referenced directly to the coordinates and do not.          */
      if (open->intype & (O_UPI|O_UPC)) {
	for (cc = 0; cc < open->npts; cc++) {
	  c = open->locu[cc];	
	  for (j = 1; j <= params->nland; j++) {
	    cco = params->lande1[j];
	    if (c >= cco) {
	      open->locu[cc]--;
	      hd_warn("Boundary cell remapping due to NLAND: %s %d->%d\n",open->name,c,open->locu[cc]);
	    }
	  }
	}
      }
      /* Set the open boundaries in the mesh structure               */
      convert_obc_list(params, open, n, NULL, cmap);
    }
  }

  if (verbose) {
    for (cc = 1; cc <= m->ns2; cc++) {
      printf("%d %d %f %f\n",cc,params->npe2[cc],params->x[cc][0],params->y[cc][0]);
      for (n = 1; n <= params->npe2[cc]; n++) {
	printf(" %d %f %f\n",n,params->x[cc][n],params->y[cc][n]);
      }
    }
  }

  if (DEBUG("init_m"))
    dlog("init_m", "\nMesh structure OBCs set OK\n");

  /*-----------------------------------------------------------------*/
  /* Write to file from the parameter structure                      */
  if (filef) {
    mesh_ofile_name(params, key);
  }
  if (filef == 1) {
    FILE *ef;
    sprintf(buf,"%s_e.txt", key);
    if ((ef = fopen(buf, "w")) == NULL) filef = 0;
    if (filef) {
      for (cc = 1; cc <= m->ns2; cc++) {
	for (n = 1; n <= params->npe2[cc]; n++) {
	  fprintf(ef, "%f %f\n",params->x[cc][n],params->y[cc][n]);
	}
	fprintf(ef, "%f %f\n",params->x[cc][1],params->y[cc][1]);
	fprintf(ef, "NaN Nan\n");
      }
      fclose(ef);
    }
    sprintf(buf,"%s_c.txt", key);
    if ((ef = fopen(buf, "w")) == NULL) filef = 0;
    if (filef) {
      for (cc = 1; cc <= m->ns2; cc++)
	fprintf(ef, "%f %f\n",params->x[cc][0],params->y[cc][0]);
      fclose(ef);
    }
    if (params->nobc) {
      sprintf(buf,"%s_b.txt", key);
      if ((ef = fopen(buf, "w")) == NULL) filef = 0;
      if (filef) {
	for (n = 0; n < params->nobc; n++) {
	  open_bdrys_t *open = params->open[n];
	  for (cc = 0; cc < open->npts; cc++) {
	    fprintf(ef, "%f %f\n", open->posx[cc][0], open->posy[cc][0]);
	    fprintf(ef, "%f %f\n", open->posx[cc][1], open->posy[cc][1]);
	    fprintf(ef, "NaN NaN\n");
	  }
	}
	fclose(ef);
      }
    }
  }

  convert_mesh_input(params, m, params->x, params->y, NULL, 1);
  if (DEBUG("init_m"))
    dlog("init_m", "\nMesh converted OK\n");

  /*-----------------------------------------------------------------*/
  /* Overwrite bathymetry values if required                         */
  if (params->runmode & (AUTO | DUMP)) set_bathy_us(params);

  /* Now performed in set_bathy_us()
  if (params->bathyfill & B_AVERAGE) {
    for (cc = 1; cc <= m->ns2; cc++) {
      if (params->bathy[cc] == LANDCELL || params->bathy[cc] == NOTVALID) {
	double bath = 0.0;
	double nb= 0.0;
	for (j = 1; j <= m->npe[cc]; j++) {
	  if ((c = m->neic[j][cc])) {
	    if (params->bathy[c] != LANDCELL && params->bathy[c] != NOTVALID) {
	      bath += params->bathy[c];
	      nb += 1.0;
	    }
	  }
	}
	params->bathy[cc] = (nb) ? bath / nb : -params->bmin;
      }
    }
  }
  */

  /*-----------------------------------------------------------------*/
  /* Write to file from the mesh structure                           */
  if (filef == 2) {
    FILE *ef;

    sprintf(buf,"%s_e.txt", key);
    if ((ef = fopen(buf, "w")) == NULL) filef = 0;
    if (filef) {
      for (cc = 1; cc <= m->ns2; cc++) {
	for (n = 1; n <= m->npe[cc]; n++) {
	  fprintf(ef, "%f %f\n",m->xloc[m->eloc[0][cc][n]], m->yloc[m->eloc[0][cc][n]]);
	}
	fprintf(ef, "%f %f\n",m->xloc[m->eloc[0][cc][1]], m->yloc[m->eloc[0][cc][1]]);
	fprintf(ef, "NaN Nan\n");
      }
      fclose(ef);
    }
    sprintf(buf,"%s_c.txt", key);
    if ((ef = fopen(buf, "w")) == NULL) filef = 0;
    if (filef) {
      for (cc = 1; cc <= m->ns2; cc++)
	fprintf(ef, "%f %f\n",m->xloc[m->eloc[0][cc][0]], m->yloc[m->eloc[0][cc][0]]);
      fclose(ef);
    }
    if (m->nobc) {
      sprintf(buf,"%s_b.txt", key);
      if ((ef = fopen(buf, "w")) == NULL) filef = 0;
      if (filef) {
	for (n = 0; n < m->nobc; n++) {
	  if (params->open[n]->intype & O_COR) {
	    fprintf(ef, "Boundary information written to file 'boundary.txt' and 'obc_spec.txt'.\n");
	    break;
	  }
	  for (cc = 1; cc <= m->npts[n]; cc++) {
	    fprintf(ef, "%f %f\n", m->xloc[m->obc[n][cc][0]], m->yloc[m->obc[n][cc][0]]);
	    fprintf(ef, "%f %f\n", m->xloc[m->obc[n][cc][1]], m->yloc[m->obc[n][cc][1]]);
	    fprintf(ef, "NaN NaN\n");
	  }
	}
	fclose(ef);
      }
    }
  }

  /* Decide whether to take the negative of the bathymetry value     */
  i = j = jj = 0;
  for (cc = 1; cc <= m->ns2; cc++) {
    double b = params->bathy[cc];
    if (b < 0) i++;
    if (b == NOTVALID) j++;
    if (b == LANDCELL) jj++;
  }
  if (!i && !j && !jj) {
    for (cc = 1; cc <= m->ns2; cc++) {
      params->bathy[cc] *= -1.0;
    }
  }

  /* Test if a closed perimeter can be made for open boundaries      */
  if (testp)
    perimeter_mask(params, m->neic);

  i_free_1d(cmap);
  d_free_2d(x);
  d_free_2d(y);
  if (DEBUG("init_m"))
    dlog("init_m", "\nMesh structure built OK\n");
}

/* END meshstruct_us()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Initialize a mesh structure.                                      */
/*-------------------------------------------------------------------*/
mesh_t *mesh_init(parameters_t *params, /* Input parameters          */
		  int ns2,              /* Mesh size                 */
		  int npe               /* Constant nodes per cell   */
		  )
{
  int cc, n;

  mesh_t *mesh = (mesh_t *)malloc(sizeof(mesh_t));
  memset(mesh, 0, sizeof(mesh_t));

  mesh->ns2 = ns2;

  /*-----------------------------------------------------------------*/
  /* Vertices per cell                                               */
  mesh->npe = i_alloc_1d(mesh->ns2+1);
  mesh->mnpe = 0;
  for (cc = 1; cc <= mesh->ns2; cc++) {
    if (npe)
      mesh->npe[cc] = npe;
    else
      mesh->npe[cc] = params->npe2[cc];
    if (mesh->npe[cc] > mesh->mnpe)
      mesh->mnpe = mesh->npe[cc];
  }

  /*-----------------------------------------------------------------*/
  /* Mesh coordinates                                                */
  mesh->yloc = d_alloc_1d((mesh->ns2 * (mesh->mnpe+1)) + 1);
  mesh->xloc = d_alloc_1d((mesh->ns2 * (mesh->mnpe+1)) + 1);
  mesh->eloc = i_alloc_3d(mesh->mnpe+1, mesh->ns2+1, 2);

  if (params->us_type & US_IJ) {
    mesh->nce1 = params->nce1;
    mesh->nce2 = params->nce2;
    mesh->iloc = i_alloc_1d(mesh->ns2+1);
    mesh->jloc = i_alloc_1d(mesh->ns2+1);
  }

  /*-----------------------------------------------------------------*/
  /* Open boundaries                                                 */
  mesh_init_OBC(params, mesh);

  return mesh;
}

/* END mesh_init()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Allocate the open boundary arrays in the mesh structure           */
/*-------------------------------------------------------------------*/
void mesh_init_OBC(parameters_t *params, /* Input parameters         */
		   mesh_t *mesh          /* Mesh info                */
		   )
{
  int n, mo2;
  int uf = -1, sf = -1;

  if (params->nobc) {

    /* Check that START_LOC and other OBC specifications are not     */
    /* mixed.                                                        */
    for (n = 0; n < params->nobc; n++) {
      open_bdrys_t *open = params->open[n];
      if (!(open->intype & O_COR)) uf = n;
      if (open->intype & O_COR) sf = n;
    }
    if (uf >=0 && sf >=0)
      hd_quit("read_grid_us: Can't mix OBC specification UPOINTS and START_LOC; boundaries %d and %d.\n", uf, sf);

    mesh->nobc = params->nobc;
    mesh->npts = i_alloc_1d(mesh->nobc);

    mo2 = 0;
    for (n = 0; n < params->nobc; n++) {
      open_bdrys_t *open = params->open[n];
      mesh->npts[n] = open->npts;
      if (open->npts > mo2) mo2 = open->npts;
    }
    mo2++;
    mesh->loc = i_alloc_2d(mo2, mesh->nobc);
    mesh->obc =  i_alloc_3d(2, mo2, mesh->nobc);
  }
}

/* END mesh_init_OBC()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Frees an input mesh structure                                     */
/*-------------------------------------------------------------------*/
void free_mesh(mesh_t *mesh)
{
  if (mesh != NULL) {
    i_free_1d(mesh->npe);
    d_free_1d(mesh->xloc);
    d_free_1d(mesh->yloc);
    i_free_3d(mesh->eloc);
    if (mesh->nobc) {
      i_free_1d(mesh->npts);
      i_free_2d(mesh->loc);
      i_free_3d(mesh->obc);
    }
    if (mesh->iloc) i_free_1d(mesh->iloc);
    if (mesh->jloc) i_free_1d(mesh->jloc);
    if (mesh->map) i_free_1d(mesh->map);
    free((mesh_t *)mesh);
  }
}

/* END free_mesh()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Creates a list of segments locations at the grid perimeter        */
/* corresponding to OBC locations.                                   */
/* WARNING: If an OBC is specified on a different segment to the     */
/* major segment (e.g. a river boundary on an island within the main */
/* domain) then the OBCs will not be specified properly. Ideally we  */
/* need a perimeter computed for every OBC using it's START_LOC as   */
/* the first point.                                                  */
/* NOTE: OBCs specified without a mid location (start location only  */
/* or start + end location only) are processed regardless of which   */
/* perimeter they reside on.                                         */
/*-------------------------------------------------------------------*/
int get_mesh_obc(parameters_t *params, 
		 int **neic
		 )
{
  FILE *fp, *sp, *op;
  int n, m, cc, c, cn, cs, i, j, jj, jn, js, jss;
  int np;               /* Size of the perimeter array               */
  int *perm;            /* Mesh indices of the perimeter             */
  double *lon, *lat;    /* Coordinates of the perimeter              */
  int filef = 1;        /* Print perimeter in 'perimeter.txt'        */
  int filep = 1;        /* Print OBCs in 'obc_spec.txt'              */
  int verbose = 0;      /* Print to screen                           */
  int dof = 1;          /* Do this routine when dof = 1              */
  int *si, *ei, mi;     /* Start, end and mid indices in perimeter   */
  int *dir;             /* Direction to traverse the perimeter       */
  int mobc = 0;         /* Maximum points in the OBCs                */
  int *mask;            /* Mask for perimeter cells visited          */
  int *donef;           /* Flag for completed obcs                   */
  int isclosed = 0;     /* Set for closed perimeters                 */
  mesh_t *mesh = params->mesh;
  double x, y, d, dist = HUGE;
  int typef = 1;        /* START_LOC only supplied. typef = 0; set   */
                        /* the obc to all unconnected edges.         */
                        /* typef = 1; set obc to the middle edge.    */
  if (!params->nobc) return(0);

  for (n = 0; n < params->nobc; n++) {
    open_bdrys_t *open = params->open[n];
    if (open->slat == NOTVALID || open->slon == NOTVALID) {
      /*open->elat == NOTVALID || open->elon == NOTVALID)*/
      dof = 0;
    }
  }
  if (!dof) return(0);

  mesh->nobc = params->nobc;
  mesh->npts = i_alloc_1d(mesh->nobc);
  donef = i_alloc_1d(mesh->nobc);
  memset(donef, 0, mesh->nobc * sizeof(int));

  /*-----------------------------------------------------------------*/
  /* Make a path of cells at the perimeter of the mesh               */
  /* First perimeter cell                                            */
  for (cc = 1; cc <= mesh->ns2; cc++) {
    if ((jss = has_neighbour(cc, mesh->npe[cc], neic)))
      break;
  }

  for (n = 0; n < params->nobc; n++) {
    open_bdrys_t *open = params->open[n];
    if (open->slon != NOTVALID && open->slat != NOTVALID) {
      for (c = 1; c <= mesh->ns2; c++) {
	x = mesh->xloc[mesh->eloc[0][c][0]];
	y = mesh->yloc[mesh->eloc[0][c][0]];
	d = sqrt((x - open->slon) * (x - open->slon) + 
		 (y - open->slat) * (y - open->slat));
	if (d < dist && (jss = has_neighbour(c, mesh->npe[c], neic))) {
	  dist = d;
	  cc = c;
	}
      }
      break;
    }
  }
  cs = cn = cc;
  n = 1;
  jn = jss;
  mask = i_alloc_1d(mesh->ns2+1);
  memset(mask, 0, (mesh->ns2+1) * sizeof(int));
  /*mask[cs] = 1;*/

  /* Count the number of perimeter cells                             */
  for (cc = 1; cc <= mesh->ns2; cc++) {
    int found = 0;
    /* Cycle edges clockwise from the perimeter edge                 */
    j = jn;
    for (jj = 1; jj <= mesh->npe[cn]; jj++) {
      c = neic[j][cn];
      j = (j >= mesh->npe[cn]) ? 1 : j + 1;
      if (c == 0) continue;
      if ((js = has_neighbour(c, mesh->npe[c], neic)) && !mask[c]) {
	cn = c;
	jn = js;
	mask[cn] = 1;
	found = 1;
	n++;
	break;
      }
    }
    /* Dead end : try to back out                                    */
    if (!found) {
      int je, cp = cn;
      j = jn;
      for (jj = 1; jj <= mesh->npe[cn]; jj++) {
	c = neic[j][cn];
	j = (j >= mesh->npe[cn]) ? 1 : j + 1;
	if (c == 0) continue;
	if ((js = has_neighbour(c, mesh->npe[c], neic))) {
	  cn = c;
	  for (je = 1; je <= mesh->npe[cn]; je++) {
	    if (neic[je][cn] == cp)
	      jn = (je + 1 > mesh->npe[cn]) ? 1 : je + 1;
	  }
	  break;
	}
      }
    }
    if (cn == cs) {
      isclosed = 1;
      break;
    }
  }
  if (!isclosed) {
    hd_warn("get_mesh_obc: Can't make a closed perimeter for mesh (%d != %d).\n", cs, cn);
    hd_warn("              Proceeding with partial perimeter: not all OBC's may be accounted for.\n");
  }
  np = n;

  /* Allocate                                                        */
  lon = d_alloc_1d(np+1);
  lat = d_alloc_1d(np+1);
  perm = i_alloc_1d(np+1);
  if (filep)
    fp = fopen("perimeter.txt", "w");

  /* Get the perimeter path                                          */
  memset(mask, 0, (mesh->ns2+1) * sizeof(int));
  n = 1;
  cn = cs;
  jn = jss;
  for (cc = 1; cc <= mesh->ns2; cc++) {
    int found = 0;
    j = jn;
    for (jj = 1; jj <= mesh->npe[cn]; jj++) {
      c = neic[j][cn];
      j = (j >= mesh->npe[cn]) ? 1 : j + 1;
      if (c == 0) continue;
      if ((js = has_neighbour(c, mesh->npe[c], neic)) && !mask[c]) {
	cn = c;
	jn = js;
	mask[cn] = 1;
	found = 1;
	lon[n] = mesh->xloc[mesh->eloc[0][cn][0]];
	lat[n] = mesh->yloc[mesh->eloc[0][cn][0]];
	if (verbose) printf("%f %f\n",lon[n],lat[n]);
	if (filep) fprintf(fp, "%f %f\n",lon[n],lat[n]);
	perm[n] = cn;
	n++;
	break;
      }
    }
    /* Dead end : try to back out                                    */
    if (!found) {
      int je, cp = cn;
      j = jn;
      for (jj = 1; jj <= mesh->npe[cn]; jj++) {
	c = neic[j][cn];
	j = (j >= mesh->npe[cn]) ? 1 : j + 1;
	if (c == 0) continue;
	if ((js = has_neighbour(c, mesh->npe[c], neic))) {
	  cn = c;
	  for (je = 1; je <= mesh->npe[cn]; je++) {
	    if (neic[je][cn] == cp)
	      jn = (je + 1 > mesh->npe[cn]) ? 1 : je + 1;
	  }
	  break;
	}
      }
    }
    if (cn == cs) break;
  }
  i_free_1d(mask);
  if (filep) fclose(fp);

  /*-----------------------------------------------------------------*/
  /* Count the boundary edges within each segment                    */
  dir = i_alloc_1d(mesh->nobc);
  si = i_alloc_1d(mesh->nobc);
  ei = i_alloc_1d(mesh->nobc);
  for (m = 0; m < mesh->nobc; m++) {
    open_bdrys_t *open = params->open[m];

    if (open->slat != NOTVALID && open->slon != NOTVALID) {
      if (open->mlat == NOTVALID && open->mlon == NOTVALID) {
	donef[m] = 1;
	if (open->elat != NOTVALID && open->elon != NOTVALID) donef[m] = 2;
      }
    }
    /*
    if (open->slat == NOTVALID || open->slon == NOTVALID ||
	open->elat == NOTVALID || open->elon == NOTVALID)
      continue;
    */
    if (donef[m]) continue;
    /* Get the indices of start, med and end locations               */
    mesh->npts[m] = 0;
    si[m] = find_mindex(open->slon, open->slat, lon, lat, np, NULL);
    ei[m] = find_mindex(open->elon, open->elat, lon, lat, np, NULL);
    /* Swap if required so that si < ei                              */
    if (ei[m] < si[m]) {
      n = si[m];
      si[m] = ei[m];
      ei[m] = n;
    }
    if (si[m] == ei[m]) {
      /* si and ei are the same location                             */
      dir[m] = 1;
    } else {
      mi = find_mindex(open->mlon, open->mlat, lon, lat, np, NULL);
      /* si is next to ei (mi==si or mi==ei)                         */
      if (mi == si[m] || mi == ei[m]) {
	dir[m] = 0;
	if (ei[m] == si[m] + 1)
	  dir[m] = 1;
	else
	  dir[m] = -1;
      } else {
	/* Get the direction to cycle along the path                 */
	dir[m] = 0;
	n = si[m];
	for (i = 1; i <= np; i++) {
	  if (n == mi) {
	    dir[m] = 1;
	    break;
	  }
	  n = (n == np) ? 1 : n+1;
	  if (n == ei[m]) break;
	}
	if (!dir[m]) {
	  n = si[m];
	  for (i = 1; i <= np; i++) {
	    if (n == mi) {
	      dir[m] = -1;
	      break;
	    }
	    n = (n == 1) ? np : n-1;
	    if (n == ei[m]) break;
	  }
	}
      }
    }
    if (!dir[m])
      hd_quit("get_mesh_obc: Can't find OBC mid-point (%f %f) in boundary %d\n", 
	      open->mlon, open->mlat, m);
    /* Count the boundary edges along the segment                    */
    i = 0;
    n = si[m];
    while (i <= np) {
      /* This can fail for closed OBCs depending on start and end
         indices.
	 for (i = si[m]; i <= ei[m]; i++) { */
      cc = perm[n];
      for (j = 1; j <= mesh->npe[cc]; j++) {
	if (!neic[j][cc]) mesh->npts[m]++;
      }
      if (n == ei[m]) break;
      n += dir[m];
      if (n <= 0 && dir[m] == -1) n = np;
      if (n > np && dir[m] == 1) n = 1;
      i++;
    }
    if (mesh->npts[m] > mobc) mobc = mesh->npts[m];
  }

  /*-----------------------------------------------------------------*/
  /* Allocate memory                                                 */
  if (mobc) {
    mesh->loc = i_alloc_2d(mobc+1, mesh->nobc);
    mesh->obc =  i_alloc_3d(2, mobc+1, mesh->nobc);
  }
  fp = fopen("boundary.txt", "w");
  sp = fopen("boundary.site", "w");
  op = fopen("obc_spec.txt", "w");

  /*-----------------------------------------------------------------*/
  /* Set boundaries without a mid point. These are included          */
  /* regardless of on which perimeter they lie.                      */
  for (m = 0; m < mesh->nobc; m++) {
    open_bdrys_t *open = params->open[m];
    if (donef[m]) {
      /* Get the information for the start location                  */
      dist = HUGE;
      mesh->npts[m] = 0;
      for (c = 1; c <= mesh->ns2; c++) {
	x = mesh->xloc[mesh->eloc[0][c][0]];
	y = mesh->yloc[mesh->eloc[0][c][0]];
	d = sqrt((x - open->slon) * (x - open->slon) + 
		 (y - open->slat) * (y - open->slat));
	if (d < dist) {
	  dist = d;
	  cs = c;
	}
      }
    }
    /* Get the information for the end location                      */
    if (donef[m] == 2) {
      dist = HUGE;
      for (c = 1; c <= mesh->ns2; c++) {
	x = mesh->xloc[mesh->eloc[0][c][0]];
	y = mesh->yloc[mesh->eloc[0][c][0]];
	d = sqrt((x - open->elon) * (x - open->elon) + 
		 (y - open->elat) * (y - open->elat));
	if (d < dist) {
	  dist = d;
	  cn = c;
	}
      }
      /* Find which vertex is common to the two edges and set the    */
      /* obc to the two common edges.                                */
      dof = 1;
      for (js = 1; js <= mesh->npe[cs] && dof; js++) {
	if (neic[js][cs]) continue;
	for (jn = 1; jn <= mesh->npe[cn] && dof; jn++) {
	  if (neic[jn][cn]) continue;
	  if (mesh->eloc[0][cs][js] == mesh->eloc[1][cn][jn] ||
	      mesh->eloc[1][cs][js] == mesh->eloc[0][cn][jn]) {
	    dof = 0;	    
	    break;
	  }
	}
      }
      if (dof) hd_warn("get_mesh_obc(): Can't find boundary cell for OBC%d %s\n", m, open->name);
      js--;
      mesh->npts[m] = 2;
      mesh->loc[m][1] = cs;
      mesh->obc[m][1][0] = mesh->eloc[0][cs][js];
      mesh->obc[m][1][1] = mesh->eloc[1][cs][js];
      mesh->loc[m][2] = cn;
      mesh->obc[m][2][0] = mesh->eloc[0][cn][jn];
      mesh->obc[m][2][1] = mesh->eloc[1][cn][jn];
    }
    /* The start coordinate only was supplied.                         */
    if (donef[m] == 1) {
      /* Set the obc to all edges that map outside the domain          */
      if (typef == 0) {
	mesh->npts[m] = 1;
	for (j = 1; j <= mesh->npe[cs]; j++) {
	  if (!neic[j][cs]) {
	    mesh->loc[m][mesh->npts[m]] = cs;
	    mesh->obc[m][mesh->npts[m]][0] = mesh->eloc[0][cs][j];
	    mesh->obc[m][mesh->npts[m]][1] = mesh->eloc[1][cs][j];
	    mesh->npts[m]++;
	  }
	}
	mesh->npts[m]--;
      } else {
	/* Set the obc to the middle (central) edge of all edges     */
	/* that map outside the domain (e.g. for a quad with three   */
	/* edges mapping outside, the obc is the second edge that    */
	/* maps out.                                                 */
	n = 0;
	for (j = 1; j <= mesh->npe[cs]; j++)
	  if (!neic[j][cs]) n++;
	jj = (n + 1) / 2;
	n = 0;
	mesh->npts[m] = 1;
	for (j = 1; j <= mesh->npe[cs]; j++) {
	  if (!neic[j][cs]) {
	    n++;
	    if (n == jj) {
	      mesh->loc[m][mesh->npts[m]] = cs;
	      mesh->obc[m][mesh->npts[m]][0] = mesh->eloc[0][cs][j];
	      mesh->obc[m][mesh->npts[m]][1] = mesh->eloc[1][cs][j];
	    }
	  }
	}
      }
      if (!mesh->npts[m]) hd_warn("get_mesh_obc(): Can't find boundary cell for OBC%d %s\n", m, open->name);
    }
    if (donef[m]) {
      /* File output if required                                     */
      if (filef) {
	if (fp != NULL) {
	  for (cc = 1; cc <= mesh->npts[m]; cc++) {
	    fprintf(fp, "%f %f\n", mesh->xloc[mesh->obc[m][cc][0]], mesh->yloc[mesh->obc[m][cc][0]]);
	    fprintf(fp, "%f %f\n", mesh->xloc[mesh->obc[m][cc][1]], mesh->yloc[mesh->obc[m][cc][1]]);
	  }
	  fprintf(fp, "NaN NaN\n");
	}
      }
      if (op  != NULL) {
	fprintf(op, "\nBOUNDARY%d.UPOINTS     %d\n", m, mesh->npts[m]);
	for (cc = 1; cc <= mesh->npts[m]; cc++) {
	  fprintf(op, "%d (%lf,%lf)-(%lf,%lf)\n", mesh->loc[m][cc], 
		  mesh->xloc[mesh->obc[m][cc][0]], mesh->yloc[mesh->obc[m][cc][0]],
		  mesh->xloc[mesh->obc[m][cc][1]], mesh->yloc[mesh->obc[m][cc][1]]);
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Save the OBC locations to the mesh structure                    */
  for (m = 0; m < mesh->nobc; m++) {
    open_bdrys_t *open = params->open[m];
    /*
    if (open->slat == NOTVALID || open->slon == NOTVALID ||
	open->elat == NOTVALID || open->elon == NOTVALID)
      continue;
    */
    if (donef[m]) continue;

    /* Save to mesh                                                  */
    i = 0;
    n = si[m];
    mesh->npts[m] = 1;
    while (i <= np) {
      cc = perm[n];
      for (j = 1; j <= mesh->npe[cc]; j++) {
	if (!neic[j][cc]) {
	  mesh->loc[m][mesh->npts[m]] = cc;
	  mesh->obc[m][mesh->npts[m]][0] = mesh->eloc[0][cc][j];
	  mesh->obc[m][mesh->npts[m]][1] = mesh->eloc[1][cc][j];
	  mesh->npts[m]++;
	}
      }
      if (n == ei[m]) break;
      n += dir[m];
      if (n <= 0 && dir[m] == -1) n = np;
      if (n > np && dir[m] == 1) n = 1;
      i++;
    }
    mesh->npts[m]--;

    /* Verbose output if required                                    */
    if (verbose) {
      printf("\nBOUNDARY%d     %d\n", m, mesh->npts[m]);
      for (cc = 1; cc <= mesh->npts[m]; cc++) {
	printf("%d (%f,%f)-(%f,%f)\n", mesh->loc[m][cc], mesh->xloc[mesh->obc[m][cc][0]],
	       mesh->yloc[mesh->obc[m][cc][0]], mesh->xloc[mesh->obc[m][cc][1]],
	       mesh->yloc[mesh->obc[m][cc][1]]);
      }
    }

    /* File output if required                                       */
    if (filef) {
      if (fp != NULL) {
	for (cc = 1; cc <= mesh->npts[m]; cc++) {
	  fprintf(fp, "%f %f\n", mesh->xloc[mesh->obc[m][cc][0]], mesh->yloc[mesh->obc[m][cc][0]]);
	  fprintf(fp, "%f %f\n", mesh->xloc[mesh->obc[m][cc][1]], mesh->yloc[mesh->obc[m][cc][1]]);
	}
	fprintf(fp, "NaN NaN\n");
      }
      if (sp != NULL) {
	double x = 0.5 * (mesh->xloc[mesh->obc[m][1][0]] + 
			  mesh->xloc[mesh->obc[m][1][1]]);
	double y = 0.5 * (mesh->yloc[mesh->obc[m][1][0]] + 
			  mesh->yloc[mesh->obc[m][1][1]]);
	fprintf(sp, "%f %f %s_s\n", x, y, open->name);
	if (donef[m] != 1) {
	  double x = 0.5 * (mesh->xloc[mesh->obc[m][mesh->npts[m]][0]] + 
			    mesh->xloc[mesh->obc[m][mesh->npts[m]][1]]);
	  double y = 0.5 * (mesh->yloc[mesh->obc[m][mesh->npts[m]][0]] + 
			    mesh->yloc[mesh->obc[m][mesh->npts[m]][1]]);
	  fprintf(sp, "%f %f %s_e\n", x, y, open->name);
	}
	if (donef[m] == 0) {
	  fprintf(sp, "%f %f %s_m\n", open->mlon, open->mlat, open->name);
	}
      }
    }
    if (op  != NULL) {
      fprintf(op, "\nBOUNDARY%d.UPOINTS     %d\n", m, mesh->npts[m]);
      for (cc = 1; cc <= mesh->npts[m]; cc++) {
	/*
	fprintf(op, "%d (%d %d)\n", mesh->loc[m][cc], mesh->obc[m][cc][0],
		mesh->obc[m][cc][1]);
	*/
	fprintf(op, "%d (%lf,%lf)-(%lf,%lf)\n", mesh->loc[m][cc], 
		mesh->xloc[mesh->obc[m][cc][0]], mesh->yloc[mesh->obc[m][cc][0]],
		mesh->xloc[mesh->obc[m][cc][1]], mesh->yloc[mesh->obc[m][cc][1]]);
      }
    }
  }
  fclose(fp);
  fclose(sp);
  fclose(op);
  i_free_1d(si);
  i_free_1d(ei);
  i_free_1d(dir);
  i_free_1d(perm);
  i_free_1d(donef);
  d_free_1d(lon);
  d_free_1d(lat);
  return(1);
}

/* END get_mesh_obc()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Gets the perimeter cell centres of a mesh                         */
/*-------------------------------------------------------------------*/
int perimeter_mask(parameters_t *params, 
		   int **neic
		   )
{
  FILE *fp;
  mesh_t *mesh = params->mesh;
  char buf[MAXSTRLEN], keyword[MAXSTRLEN];
  int n, m, cc, c, cn, cs, i, j, jj, jn, js, np;
  int *mask;
  int verbose = 0;
  int isclosed = 0;
  int found = 0;

  fp = fopen("perimeter.txt", "w");

  /*-----------------------------------------------------------------*/
  /* Make a path of cells at the perimeter of the mesh               */
  /* First perimeter cell                                            */
  for (n = 0; n < params->nobc; n++) {
    open_bdrys_t *open = params->open[n];
    double d, xd, yd, dmin = HUGE;
    if (open->slon != NOTVALID && open->slat != NOTVALID) {
      for (cc = 1; cc <= mesh->ns2; cc++) {
	xd = mesh->xloc[mesh->eloc[0][cc][0]] - open->slon;
	yd = mesh->yloc[mesh->eloc[0][cc][0]] - open->slat;
	d = sqrt(xd * xd + yd * yd);
	if (d < dmin) {
	  cn = cc;
	  dmin = d;
	}
      }
      cc = cn;
      js = has_neighbour(cc, mesh->npe[cc], neic);
      found = 1;
      break;
    }
  }
  if (!found) {
    for (cc = 1; cc <= mesh->ns2; cc++) {
      if ((js = has_neighbour(cc, mesh->npe[cc], neic)))
	break;
    }
  }

  /* Set the mask for all cells surrounded by other cells            */
  mask = i_alloc_1d(mesh->ns2+1);
  memset(mask, 0, (mesh->ns2+1) * sizeof(int));
  /*
  for (cc = 1; cc <= mesh->ns2; cc++) {
    if (has_neighbour(cc, mesh->npe[cc], neic) == 0 && !mask[cc])
      mask[cc] = 1;
  }
  */

  cs = cn = cc;
  n = 1;
  jn = js;
  /*mask[cs] = 1;*/
  printf("Start coordinate = %d[%f %f]\n",cs, mesh->xloc[mesh->eloc[0][cs][0]], mesh->yloc[mesh->eloc[0][cs][0]]);
  printf("Edge direction = %d\n", jn);
  /* Count the number of perimeter cells                             */
  for (cc = 1; cc <= mesh->ns2; cc++) {
    int found = 0;
    /* Cycle edges clockwise from the perimeter edge                 */
    j = jn;
    for (jj = 1; jj <= mesh->npe[cn]; jj++) {
      c = neic[j][cn];
      j = (j >= mesh->npe[cn]) ? 1 : j + 1;
      /*
      fprintf(fp, "cell %d : j=%d cn=%d[%f %f]\n",cn, j, c, mesh->xloc[mesh->eloc[0][c][0]], mesh->yloc[mesh->eloc[0][c][0]]);
      */
      if (c == 0) continue;
      /* If the cell hasn't been encountered, add it to the          */
      /* perimeter.                                                  */
      if ((js = has_neighbour(c, mesh->npe[c], neic)) && !mask[c]) {
	cn = c;
	jn = js;
	mask[cn] = 1;
	found = 1;
	fprintf(fp, "%f %f %d %d\n",mesh->xloc[mesh->eloc[0][cn][0]], mesh->yloc[mesh->eloc[0][cn][0]], cn, jn);
	n++;
	break;
      }
    }
    /* Dead end : try to back out                                    */
    if (!found) {
      int je, cp = cn;
      j = jn;
      for (jj = 1; jj <= mesh->npe[cn]; jj++) {
	c = neic[j][cn];
	j = (j >= mesh->npe[cn]) ? 1 : j + 1;
	if (c == 0) continue;
	if ((js = has_neighbour(c, mesh->npe[c], neic))) {
	  cn = c;
	  for (je = 1; je <= mesh->npe[cn]; je++) {
	    if (neic[je][cn] == cp)
	      jn = (je + 1 > mesh->npe[cn]) ? 1 : je + 1;
	  }
	  /*fprintf(fp, "dead end %d %f %f %d %d\n",cc, mesh->xloc[mesh->eloc[0][cn][0]], mesh->yloc[mesh->eloc[0][cn][0]], cn, jn);*/
	  break;
	}
      }
    }

    if (cn == cs) {
      isclosed = 1;
      break;
    }
  }
  printf("End coordinate = %d[%f %f]\n",cn, mesh->xloc[mesh->eloc[0][cn][0]], mesh->yloc[mesh->eloc[0][cn][0]]);
  printf("Number perimeter cells = %d\n", n);
  i_free_1d(mask);
  fclose(fp);
  if (!isclosed) hd_quit("perimeter_mask: Can't make a closed perimeter for mesh (%d != %d).\n", cs, cn);
  return(n);
}

/* END perimeter_mask()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Writes grid infromation in a mesh structure to a nominated file   */
/*-------------------------------------------------------------------*/
void write_mesh_us(parameters_t *params, 
		   double *bathy,
		   FILE *op,
		   int mode     /* 1 = _us.txt output with OBCs      */
		                /* 2 = _us.txt output without OBCs   */
		                /* 4 = grid output from params       */
		                /* 8 = grid output from mesh         */
		   )
{
  int n, cc, i, j, m;
  mesh_t *mesh = params->mesh;
  int ns2 = mesh->ns2;

  /*-----------------------------------------------------------------*/
  /* Write the mesh info to a file with pointer op                   */
  if (mode == 1 || mode == 2) {
    fprintf(op, "\nCoordinates\n");
    for (n = 1; n <= mesh->ns; n++) {
      fprintf(op, "%d %f %f\n", n, mesh->xloc[n], mesh->yloc[n]);
    }

    fprintf(op, "\nIndices\n");
    for (cc = 1; cc <= ns2; cc++) {
      if (mesh->iloc && mesh->jloc)
	fprintf(op, "%d %d %d : %d %d\n",cc, mesh->npe[cc], mesh->eloc[0][cc][0],
		mesh->iloc[cc], mesh->jloc[cc]);
      else
	fprintf(op, "%d %d %d\n",cc, mesh->npe[cc], mesh->eloc[0][cc][0]);
      for (j = 1; j <= mesh->npe[cc]; j++)
	fprintf(op, "%d %d %d\n", j, mesh->eloc[0][cc][j], mesh->eloc[1][cc][j]);
    }
    
    if (mode == 1) {
      if (mesh->nobc) {
	fprintf(op, "\nNBOUNDARIES    %d\n", mesh->nobc);
	for (n = 0; n < mesh->nobc; n++) {
	  fprintf(op, "BOUNDARY%1d.NPOINTS  %d\n", n, mesh->npts[n]);
	  for (cc = 1; cc <= mesh->npts[n]; cc++) {
	    fprintf(op, "%d (%d %d)\n", mesh->loc[n][cc], mesh->obc[n][cc][0], mesh->obc[n][cc][1]);
	  }
	}
      } else
	fprintf(op, "\nNBOUNDARIES    0\n");
    }
    if (params->bathy != NULL) {
      fprintf(op, "\nBATHY   %d\n", ns2);
      for (cc = 1; cc <= ns2; cc++)
	fprintf(op, "%f\n", params->bathy[cc]);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Write the mesh info to files suitable to plot; edge _e.txt,     */
  /* centre _c.txt and boundary _b.txt. Mesh information is          */
  /* extracted from the params structure.                            */
  if (mode == 4) {
    FILE *ef;
    char key[MAXSTRLEN], buf[MAXSTRLEN];
    int filef = (params->meshinfo) ? 1 : 0;

    if (!filef) return;
    mesh_ofile_name(params, key);
    sprintf(buf,"%s_e.txt", key);
    if ((ef = fopen(buf, "w")) == NULL) filef = 0;
    if (filef) {
      for (cc = 1; cc <= mesh->ns2; cc++) {
	for (n = 1; n <= mesh->npe[cc]; n++) {
	  fprintf(ef, "%f %f\n",mesh->xloc[mesh->eloc[0][cc][n]], mesh->yloc[mesh->eloc[0][cc][n]]);
	}
	fprintf(ef, "%f %f\n",mesh->xloc[mesh->eloc[0][cc][1]], mesh->yloc[mesh->eloc[0][cc][1]]);
	fprintf(ef, "NaN Nan\n");
      }
      fclose(ef);
    }
    sprintf(buf,"%s_c.txt", key);
    if ((ef = fopen(buf, "w")) == NULL) filef = 0;
    if (filef) {
      for (cc = 1; cc <= mesh->ns2; cc++)
	fprintf(ef, "%f %f\n",mesh->xloc[mesh->eloc[0][cc][0]], mesh->yloc[mesh->eloc[0][cc][0]]);
      fclose(ef);
    }
    if (mesh->nobc) {
      sprintf(buf,"%s_b.txt", key);
      if ((ef = fopen(buf, "w")) == NULL) filef = 0;
      if (filef) {
	for (n = 0; n < mesh->nobc; n++) {
	  for (cc = 1; cc <= mesh->npts[n]; cc++) {
	    fprintf(ef, "%f %f\n", mesh->xloc[mesh->obc[n][cc][0]], mesh->yloc[mesh->obc[n][cc][0]]);
	    fprintf(ef, "%f %f\n", mesh->xloc[mesh->obc[n][cc][1]], mesh->yloc[mesh->obc[n][cc][1]]);
	    fprintf(ef, "NaN NaN\n");
	  }
	}
	fclose(ef);
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Write the mesh info to files suitable to plot; edge _e.txt,     */
  /* centre _c.txt and boundary _b.txt. Mesh information is          */
  /* extracted from the mesh structure.                              */
  if (mode == 8) {
    FILE *ef;
    char key[MAXSTRLEN], buf[MAXSTRLEN];
    int filef = (params->meshinfo) ? 1 : 0;

    if (!filef) return;
    if (mesh == NULL) return;
    mesh_ofile_name(params, key);
    sprintf(buf,"%s_e.txt", key);
    if ((ef = fopen(buf, "w")) == NULL) filef = 0;
    if (filef) {
      for (cc = 1; cc <= mesh->ns2; cc++) {
	for (n = 1; n <= params->npe2[cc]; n++) {
	  fprintf(ef, "%f %f\n",params->x[cc][n],params->y[cc][n]);
	}
	fprintf(ef, "%f %f\n",params->x[cc][1],params->y[cc][1]);
	fprintf(ef, "NaN Nan\n");
      }
      fclose(ef);
    }
    sprintf(buf,"%s_c.txt", key);
    if ((ef = fopen(buf, "w")) == NULL) filef = 0;
    if (filef) {
      for (cc = 1; cc <= mesh->ns2; cc++)
	fprintf(ef, "%f %f\n",params->x[cc][0],params->y[cc][0]);
      fclose(ef);
    }
    if (params->nobc) {
      sprintf(buf,"%s_b.txt", key);
      if ((ef = fopen(buf, "w")) == NULL) filef = 0;
      if (filef) {
	for (n = 0; n < params->nobc; n++) {
	  open_bdrys_t *open = params->open[n];
	  for (cc = 0; cc < open->npts; cc++) {
	    fprintf(ef, "%f %f\n", open->posx[cc][0], open->posy[cc][0]);
	    fprintf(ef, "%f %f\n", open->posx[cc][1], open->posy[cc][1]);
	    fprintf(ef, "NaN NaN\n");
	  }
	}
	fclose(ef);
      }
    }
  }
}

/* END write_mesh_us()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads and populates the input mesh structure                      */
/*-------------------------------------------------------------------*/
int read_mesh_us(parameters_t *params)
{
  mesh_t *mesh;
  char buf[MAXSTRLEN], oname[MAXSTRLEN];
  char *fields[MAXSTRLEN * MAXNUMARGS];
  FILE *op;
  int i, n, c, cc, np, m, ns, j;
  int ns2, npe;
  int verbose = 0;
  int rbf = 0;
  int *mask;

  /* Read the mesh dimensions                                        */
  op = params->prmfd;
  memset(projection, 0, sizeof(projection));
  prm_read_char(op, "PROJECTION", params->projection);

  /*
  if (prm_read_char(op, "JIGSAW_FILE", buf))
    return(-1);
  prm_read_int(op, "nMaxMesh2_face_nodes", &npe);
  prm_read_int(op, "nMesh2_face_indices", &ns);
  prm_read_int(op, "nMesh2_face", &ns2);
  */
  npe = params->npe;
  ns = params->ns;
  ns2 = params->ns2;
  params->nce1 = params->nce2 = 0;
  if (params->nce1 && params->nce2) params->us_type |= US_IJ;
  params->us_type |= US_IUS;

  params->mesh = mesh_init(params, ns2, npe);
  mesh = params->mesh;
  mesh->ns = ns;

  /* Set the grid code                                               */
  params->us_type |= US_MIX;
  if (params->npe == 3) params->us_type |= US_TRI;
  if (params->npe == 4) params->us_type |= US_QUAD;
  if (params->npe == 6) params->us_type |= US_HEX;

  if (verbose) {
    printf("nMaxMesh2_face_nodes  %d\n", mesh->mnpe);
    printf("nMesh2_face           %d\n", mesh->ns2);
    printf("nMesh2_face_indices   %d\n", mesh->ns);
  }

  /* The grid infromation is read from netCDF in dumpdata_read_us()  */
  /* Check the number of entries                                     */
  prm_skip_to_end_of_key(op, "Indices");
  prm_read_char(op, "1", buf);
  n = parseline(buf, fields, MAXNUMARGS);
  if (!(params->us_type & US_IJ) && n > 3)
    hd_quit("read_mesh_us: Too many entries (%d) for grid topology (indicies not required).\n", n);
  if (params->us_type & US_IJ && n <= 3)
    hd_quit("read_mesh_us: Too few entries (%d) for grid topology (indicies required).\n", n);
    
  /* Read the entries                                                */
  /*
  fprintf(op, "\nCoordinates\n");
  for (n = 1; n <= mesh->ns; n++) {
    fprintf(op, "%d %f %f\n", n, mesh->xloc[n], mesh->yloc[n]);
  }
  */

  prm_skip_to_end_of_key(op, "Coordinates");
  if (verbose) printf("\nCoordinates\n");
  for (cc = 1; cc <= mesh->ns; cc++) {
    i = fscanf(op, "%d %lf %lf", &m, &mesh->xloc[cc], &mesh->yloc[cc]);
    if (verbose) printf("%d %f %f\n",cc, mesh->xloc[cc], mesh->yloc[cc]);
  }
  prm_skip_to_end_of_key(op, "Indices");
  if (verbose) printf("\nIndices\n");
  cc = 1;
  if (ns2) mask = i_alloc_1d(ns2 + 1);
  for (c = 1; c <= ns2; c++) {
    if (params->us_type & US_IJ) {
      fscanf(op, "%d %d %d : %d %d",&i, &mesh->npe[cc], &mesh->eloc[0][cc][0],
	     &mesh->iloc[cc], &mesh->jloc[cc]);
      if (verbose) printf("%d %d %d : %d %d\n",cc, mesh->npe[cc], mesh->eloc[0][cc][0], mesh->iloc[cc], mesh->jloc[cc]);
    } else {
      fscanf(op, "%d %d %d",&i, &mesh->npe[cc], &mesh->eloc[0][cc][0]);
      if (verbose) printf("%d %d %d\n",cc, mesh->npe[cc], mesh->eloc[0][cc][0]);
    }
    mesh->eloc[1][cc][0] = mesh->eloc[0][cc][0];
    for (j = 1; j <= mesh->npe[cc]; j++) {
      fscanf(op, "%d %d %d\n",&i, &mesh->eloc[0][cc][j], &mesh->eloc[1][cc][j]);
      if (verbose) printf(" %d %d %d\n", j, mesh->eloc[0][cc][j], mesh->eloc[1][cc][j]);
    }
    if (mesh->npe[cc] == 0) {
      mask[c] = 0;
      mesh->ns2--;
    } else {
      mask[c] = cc;
      cc++;
    }
  }

  /* Open boundaries                                                 */
  if (rbf) {
  prm_read_int(op, "NBOUNDARIES", &mesh->nobc);
  if (mesh->nobc) {
    int uf = -1, sf = -1;
    mesh->npts = i_alloc_1d(mesh->nobc);
    i = 0;
    for (n = 0; n < mesh->nobc; n++) {
      sprintf(buf, "BOUNDARY%1d.UPOINTS", n);
      if (prm_read_int(op, buf, &mesh->npts[n])) {
	if (mesh->npts[n] > i) i = mesh->npts[n];
	uf = n;
      }
      sprintf(buf, "BOUNDARY%1d.START_LOC", n);
      if (prm_skip_to_end_of_key(op, buf)) sf = n;
    }
    if (uf >=0 && sf >=0)
      hd_quit("read_grid_us: Can't mix OBC specification UPOINTS and START_LOC; boundaries %d and %d.\n", uf, sf);
    if (i) {
      mesh->loc = i_alloc_2d(i+1, mesh->nobc);
      mesh->obc =  i_alloc_3d(2, i+1, mesh->nobc);
      for (n = 0; n < mesh->nobc; n++) {
	sprintf(buf, "BOUNDARY%1d.UPOINTS", n);
	prm_skip_to_end_of_key(op, buf);
	prm_flush_line(op);
	for (cc = 1; cc <= mesh->npts[n]; cc++) {
	  if (fscanf(op, "%d (%d %d)", &c, &mesh->obc[n][cc][0], &mesh->obc[n][cc][1]) != 3)
	    hd_warn("read_mesh_us: Can't read edge coordinates in boundary points list. Format '1 (2 3)'\n");
	  /* Flush the remainder of the line */
	  prm_flush_line(op);
	  if (!mask[c]) hd_warn("read_mesh_us: OBC%d cell %d(%d) is an empty cell.\n", n, cc, c);
	  mesh->loc[n][cc] = mask[c];
	}
      }
    }
  } else {
    mesh->npts = NULL;
    mesh->loc = NULL;
    mesh->obc = NULL;
  }
  }

  /* Bathymetry                                                      */
  prm_read_int(op, "BATHY", &np);
  if (verbose) printf("\nBATHY\n");
  if (params->bathy != NULL) d_free_1d(params->bathy);
  params->bathy = d_alloc_1d(mesh->ns2+1);
  for (c = 1; c <= ns2; c++) {
    double d1;
    fscanf(op, "%lf", &d1);
    if ((cc = mask[c])) {
      params->bathy[cc] = d1;
      if (verbose) printf("%f\n", params->bathy[cc]);
    }
  }

  /* Reorder the input coordinates if required                       */
  if (params->mrf & MR_READ)
    reorder_mesh(params, mesh, params->bathy);

  /* Reconstruct the parameter coordinate arrays                     */
  if (params->x == NULL && params->y == NULL) {
    params->x = d_alloc_2d(mesh->mnpe+1, mesh->ns2+1);
    params->y = d_alloc_2d(mesh->mnpe+1, mesh->ns2+1);
    for (c = 1; c <= ns2; c++) {
      if ((cc = mask[c])) {
	params->x[cc][0] = mesh->xloc[mesh->eloc[0][cc][0]];
	params->y[cc][0] = mesh->yloc[mesh->eloc[0][cc][0]];
	for (j = 1; j <= mesh->npe[cc]; j++) {
	  params->x[cc][j] = mesh->xloc[mesh->eloc[0][cc][j]];
	  params->y[cc][j] = mesh->yloc[mesh->eloc[0][cc][j]];
	}
      }
    }
  }

  neighbour_finder(mesh);

  /*-----------------------------------------------------------------*/
  /* Expand the mesh if required                                     */
  mesh_expand(params, params->bathy, params->x, params->y);
  mesh_reduce(params, params->bathy, params->x, params->y);

  set_params_mesh(params, mesh, params->x, params->y);

  write_mesh_us(params, params->bathy, 0, 4);

  if (mask) i_free_1d(mask);
  return(1);
}

/* END read_mesh_us()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Re-orders the input mesh indices to those supplied by file.       */
/*-------------------------------------------------------------------*/
void reorder_mesh(parameters_t *params, mesh_t *mesh, double *bathy)
{
  FILE *fp;
  mesh_t *nmesh;
  double *nbathy;
  int *order;
  int mnpe, npe = params->npe;
  int ns;
  int ns2;
  int c, cc, j;

  /* Read the re-mapping                                             */
  if (strlen(params->mesh_reorder)) {
    if ((fp = fopen(params->mesh_reorder, "r")) == NULL) {
      hd_warn("Can't open reorder read file %s\n", params->mesh_reorder);
      return;
    }
    prm_read_int(fp, "nMaxMesh2_face_nodes", &mnpe);
    prm_read_int(fp, "nMesh2_face", &ns2);
    prm_read_int(fp, "nMesh2_face_indices", &ns);
    if (mnpe != mesh->mnpe || ns2 != mesh->ns2 || ns != mesh->ns) {
      hd_warn("Index reorder file has incompatible dimensions\n");
      hd_warn("mnpe = %d & %d\n", mnpe, mesh->mnpe);
      hd_warn("ns2 = %d & %d\n", ns2, mesh->ns2);
      hd_warn("ns2 = %d & %d\n", ns, mesh->ns);
      return;
    }
    order = i_alloc_1d(ns2 + 1);
    for (c = 1; c <= ns2; c++) {
      fscanf(fp, "%d", &order[c]);
    }
    fclose(fp);
  }

  /* Initialise a mesh copy                                          */
  nmesh = mesh_init(params, ns2, npe);
  nbathy = d_alloc_1d(ns2 + 1);

  /* Copy mesh information                                           */
  for (c = 1; c <= ns2; c++) {
    nmesh->npe[c] = mesh->npe[c];
    nmesh->eloc[0][c][0] = mesh->eloc[0][c][0];
    if (params->us_type & US_IJ) {
      nmesh->iloc[c] = mesh->iloc[c];
      nmesh->jloc[c] = mesh->jloc[c];
    }
    nmesh->eloc[1][c][0] = mesh->eloc[0][c][0];
    for (j = 1; j <= mesh->npe[c]; j++) {
      nmesh->eloc[0][c][j] = mesh->eloc[0][c][j];
      nmesh->eloc[1][c][j] = mesh->eloc[1][c][j];
    }
    nbathy[c] = bathy[c];
  }

  /* Reorder                                                         */
  for (cc = 1; cc <= ns2; cc++) {
    c = order[cc];
    mesh->npe[cc] = nmesh->npe[c];
    mesh->eloc[0][cc][0] = nmesh->eloc[0][c][0];
    if (params->us_type & US_IJ) {
      mesh->iloc[cc] = nmesh->iloc[c];
      mesh->jloc[cc] = nmesh->jloc[c];
    }
    mesh->eloc[1][cc][0] = nmesh->eloc[0][c][0];
    for (j = 1; j <= nmesh->npe[c]; j++) {
      mesh->eloc[0][cc][j] = nmesh->eloc[0][c][j];
      mesh->eloc[1][cc][j] = nmesh->eloc[1][c][j];
    }
    bathy[cc] = nbathy[c];
  }

  /* Free mesh                                                       */
  d_free_1d(nmesh->yloc);
  d_free_1d(nmesh->xloc);
  i_free_3d(nmesh->eloc);
  if (params->us_type & US_IJ) {
    i_free_1d(nmesh->iloc);
    i_free_1d(nmesh->jloc);
  }
  if(params->nobc) {
    i_free_1d(nmesh->npts);
    i_free_2d(nmesh->loc);
    i_free_3d(nmesh->obc);
  }
  free((mesh_t *)nmesh);
  d_free_1d(nbathy);
  i_free_1d(order);
}

/* END reorder_mesh()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Re-orders the input mesh indices to those supplied by file.       */
/*-------------------------------------------------------------------*/
void reorder_bathy(parameters_t *params)
{
  FILE *fp;
  mesh_t *mesh = params->mesh;
  double *bathy = params->bathy;
  double *nbathy;
  int *order;
  int mnpe, npe = params->npe;
  int ns;
  int ns2;
  int c, cc;

  /* Read the re-mapping                                             */
  if (strlen(params->mesh_reorder)) {
    if ((fp = fopen(params->mesh_reorder, "r")) == NULL) {
      hd_warn("Can't open reorder read file %s\n", params->mesh_reorder);
      return;
    }
    prm_read_int(fp, "nMaxMesh2_face_nodes", &mnpe);
    prm_read_int(fp, "nMesh2_face", &ns2);
    prm_read_int(fp, "nMesh2_face_indices", &ns);
    if (mnpe != mesh->mnpe || ns2 != mesh->ns2 || ns != mesh->ns) {
      hd_warn("Index reorder file has incompatible dimensions\n");
      hd_warn("mnpe = %d & %d\n", mnpe, mesh->mnpe);
      hd_warn("ns2 = %d & %d\n", ns2, mesh->ns2);
      hd_warn("ns2 = %d & %d\n", ns, mesh->ns);
      return;
    }
    order = i_alloc_1d(ns2 + 1);
    for (c = 1; c <= ns2; c++) {
      fscanf(fp, "%d", &order[c]);
    }
    fclose(fp);
  }

  /* Initialise a mesh copy                                          */
  nbathy = d_alloc_1d(ns2 + 1);

  /* Copy mesh information                                           */
  for (c = 1; c <= ns2; c++) {
    nbathy[c] = bathy[c];
  }

  /* Reorder                                                         */
  for (cc = 1; cc <= ns2; cc++) {
    c = order[cc];
    bathy[cc] = nbathy[c];
  }

  /* Free mesh                                                       */
  d_free_1d(nbathy);
  i_free_1d(order);
}

/* END reorder_bathy()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Creates a circular coastline msh and resolution hfun file         */
/* for JIGSAW. This used to be done in lat/lon space but is now done */
/* in stereographic projection with units of metres. The inverse     */
/* mapping back to geog space then ensures that the grid is evenly   */
/* spaced on the sphere.                                             */
/*-------------------------------------------------------------------*/
#ifdef HAVE_JIGSAWLIB
void create_hex_radius(double crad,   /* Radius in metres     */
		       double x00,    /* Centre longitude     */
		       double y00,    /* Centre latitude      */
		       double rmin,   /* Min. resolution in m */
		       double rmax,   /* Max. resolution in m */
		       double gscale, /* Gaussian func scaler */
		       jigsaw_msh_t *J_mesh,
		       int powf       /* Power mesh flag      */
		       )
{
  /* Constants */
  int cst_res  = 250;    /* resolution for coastline in m       */
  int nHfun    = 100;    /* (Nhfun x Nhfun) resolution matrix   */
  double deg2m = 60.0 * 1852.0;
  double x00r, y00r;
  
  /* local variables */
  double cang = crad / deg2m; /* Radius in degress */
  int npts, n, i, j;
  int jret = 0;
  
  /* Jigsaw variables */
  jigsaw_jig_t J_jig;
  jigsaw_msh_t J_geom;
  jigsaw_msh_t J_hfun;
  jigsaw_msh_t J_hfun_new;

  /* Convenient pointers */
  jigsaw_VERT2_array_t *jpoints = &J_geom._vert2;
  jigsaw_EDGE2_array_t *jedges  = &J_geom._edge2;

  jigsaw_REALS_array_t *hfx    = &J_hfun._xgrid;
  jigsaw_REALS_array_t *hfy    = &J_hfun._ygrid;
  jigsaw_REALS_array_t *hfvals = &J_hfun._value;

  /* Jigsaw initialisation routines */
  jigsaw_init_jig_t (&J_jig);
  jigsaw_init_msh_t (&J_geom);
  jigsaw_init_msh_t (&J_hfun);
  if (powf) J_jig._optm_dual = 1;

  /* Flag input mesh type */
  J_geom._flags = JIGSAW_EUCLIDEAN_MESH;

  /* Number of points needed for the coastline */
  npts = 2*PI*crad / cst_res;
  
  /* Allocate JIGSAW geom (input mesh) arrays */
  jigsaw_alloc_vert2(jpoints, npts);
  jigsaw_alloc_edge2(jedges,  npts);

  /* Origin in radians */
  x00r = DEG2RAD(x00);
  y00r = DEG2RAD(y00);
  
  /* Iterate using polar coordinates */
  for (n=0; n<npts; n++) {
    double theta = (2*PI*n)/npts;
    /* Points */
    jpoints->_data[n]._ppos[0] = crad * cos(theta);
    jpoints->_data[n]._ppos[1] = crad * sin(theta);
    jpoints->_data[n]._itag = 0;
    /* Connecting edges */
    jedges->_data[n]._itag = 0;
    jedges->_data[n]._node[0] = n;
    if (n<npts-1)
      jedges->_data[n]._node[1] = n+1;
    else
      jedges->_data[n]._node[1] = 0; // close polygon
  }
  /*
   * Create the gaussian resolution matrix
   *   Uniformly spaced grid
   */
  jigsaw_alloc_reals(hfx, nHfun+1);
  jigsaw_alloc_reals(hfy, nHfun+1);
  jigsaw_alloc_reals(hfvals, (nHfun+1)*(nHfun+1));

  /* Fill in the x & y grids */
  for (n=0; n<=nHfun; n++) {
    hfx->_data[n] = -crad + 2*n*crad/nHfun;
    hfy->_data[n] = -crad + 2*n*crad/nHfun;
  }

  /*
   * Fill in the gaussian values
   *
   * Note: this needs to be in COLUMN major format but it doesn't really
   *       matter in this case due to the symmetry in stereographic projection
   *       space.
   */
  /* wgt = (exp(-(x*x + y*y) / (4 * scale * scale))); */
  gscale *= crad; /* relative spread */
  n = 0;
  for (i=0; i<=nHfun; i++)
    for (j=0; j<=nHfun; j++) {
      double gval = wgt_gaussian_2d(hfx->_data[i], /* x grid value */
				    hfy->_data[j], /* y grid value */
				    gscale);       /* scale        */
      hfvals->_data[n++] = ((rmax * (1.-gval)) + rmin);
    }


  /* Set up some flags/parameters */
  J_hfun._flags = JIGSAW_EUCLIDEAN_GRID;
  J_jig._hfun_scal = JIGSAW_HFUN_ABSOLUTE ;
  J_jig._hfun_hmin = rmin;
  J_jig._hfun_hmax = rmax;

  // Call Jigsaw
  J_jig._verbosity = 0;
  jret = jigsaw_make_mesh(&J_jig, &J_geom, NULL, &J_hfun, J_mesh);
  if (jret) hd_quit("Error calling JIGSAW\n");

  /* Perform inverse stereographic projection of output data         */
  st_transform(J_mesh, ST3PROJ_INV, x00r, y00r);

  /* Cleanup */
  jigsaw_free_msh_t(&J_geom);
  jigsaw_free_msh_t(&J_hfun);
  
}

/* END create_hex_radius()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Creates elliptic mesh with an inner ellipse having a constant     */
/* high reolution, and transitioning to lower resolution on an outer */
/* ellipse.                                                          */
/*-------------------------------------------------------------------*/
void create_ellipse(parameters_t *params,
		    long int nce1, /* Number cells in e1 direction   */
		    long int nce2, /* Number cells in e2 direction   */
		    double x00,    /* x origin offset                */
		    double y00,    /* y origin offset                */
		    double flon,   /* False longitude                */
		    double flat,   /* False latitude                 */
		    double xinc,   /* x resolution                   */
		    double yinc,   /* y resolution                   */
		    double elf,    /* Stretch; large=circular        */   
		    double ores,   /* Outser resolution              */
		    jigsaw_msh_t *J_mesh,  /* JIGSAW mesh            */
		    jigsaw_msh_t *J_hfun   /* JIGSAW weighting       */
		    )
{
  char polyname[MAXSTRLEN], obcname[MAXSTRLEN];
  char buf[MAXSTRLEN];
  long i, j, m;
  double fx00, fy00, xval, yval;
  double rflon = DEG2RAD(flon);
  double rflat = DEG2RAD(flat);
  double x0, y0, xp, yp, dx, dy, rot, x1, y1, x2, y2, dist, xc, yc;
  double as, bs, as2, bs2, inc, d1, d2;
  double *xw, *yw, *xb, *yb;
  double mnlon;
  double mxlon;
  double mnlat;
  double mxlat;
  double xmid = 0.0, ymid = 0.0;
  FILE *fp, *op, *hf;
  int npw, npo;
  double hmin, hmax, bmin, bmax, res, s, *b;
  int mce1, mce2, mhfun, *pmask;
  double deg2m = 60.0 * 1852.0;
  int jret;
  int filef = 0;         /* Print to ellipses file                   */
  int powf = 1;          /* Power mesh flag                          */
  int stproj = 1;        /* Stereographic projection flag            */
  int expf = 0;          /* Resolution stretching flag               */
  int verbose = 0;       /* Verbosity                                */

  /* Jigsaw variables                                                */
  jigsaw_jig_t J_jig;
  jigsaw_msh_t J_geom;
  jigsaw_msh_t J_hfun_new;

  /* Convenient pointers                                             */
  jigsaw_VERT2_array_t *jpoints = &J_geom._vert2;
  jigsaw_EDGE2_array_t *jedges  = &J_geom._edge2;

  jigsaw_REALS_array_t *hfx    = &J_hfun->_xgrid;
  jigsaw_REALS_array_t *hfy    = &J_hfun->_ygrid;
  jigsaw_REALS_array_t *hfvals = &J_hfun->_value;
  poly_t **pl;
  double *polyres;

  /* Set defaults                                                   */
  strcpy(polyname, "poly.xy");
  strcpy(obcname, "obc.xy");

  /* Make the ellipses                                              */
  ellipse_mesh(nce1, nce2, x00, y00, flon, flat, xinc, yinc, elf, ores,
	       &mnlon, &mxlon, &mnlat, &mxlat,
	       &xw, &yw, &xb, &yb, &npw, &npo);

  /* Print the inner                                                */
  fp = fopen(polyname, "w");
  d1 = 0.0;
  for (j = 0; j < npw; j++) {
    fprintf(fp, "%f %f\n", xw[j], yw[j]);
    x1 = xw[j] - x00;
    y1 = yw[j] - y00;
    d1 += sqrt(x1 * x1 + y1 * y1);
  }
  d1 /= (double)npw;
  fclose(fp);

  /* Save to poly                                                   */
  pl = (poly_t **)malloc(npw * sizeof(poly_t));
  polyres = d_alloc_1d(1);
  polyres[0] = 0.5 * (xinc + yinc);
  pl[0] = poly_create();
  for (j = 0; j < npw; j++) {
    poly_add_point(pl[0], xw[j], yw[j]);
  }

  /* Print the bounding ellipse                                     */
  fp = fopen(obcname, "w");
  d2 = 0.0;
  for (j = 0; j < npo; j++) {
    fprintf(fp, "%f %f\n", xb[j], yb[j]);
    x2 = xb[j] - x00;
    y2 = yb[j] - y00;
    d2 += sqrt(x2 * x2 + y2 * y2);
  }
  d2 /= (double)npo;
  dist = d2 - d1;
  fclose(fp);
  if (verbose) {
    printf("dist %f %f %f\n",d1*deg2m,d2*deg2m,dist*deg2m);
  }

  /* Weighting: hmin is resolution, bmin is distance                */
  hmin = polyres[0];
  hmax = ores;
  bmin = 200.0 / deg2m;
  bmax = dist;
  if (verbose) printf("ores %f %f\n",ores,ores*deg2m);

  /* Jigsaw initialisation routines                                 */
  jigsaw_init_jig_t (&J_jig);
  jigsaw_init_msh_t (&J_geom);

  init_J_jig(&J_jig);

  J_jig._optm_dual = 1;

  /* Flag input mesh type                                           */
  J_geom._flags = JIGSAW_EUCLIDEAN_MESH;

  /* Allocate JIGSAW geom (input mesh) arrays                       */
  jigsaw_alloc_vert2(jpoints, npo);
  jigsaw_alloc_edge2(jedges,  npo);
  /* Iterate using polar coordinates                                */
  for (i=0; i < npo; i++) {
    /* Points                                                       */
    jpoints->_data[i]._ppos[0] = xb[i];
    jpoints->_data[i]._ppos[1] = yb[i];
    jpoints->_data[i]._itag = 0;
    /* Connecting edges                                             */
    jedges->_data[i]._itag = 0;
    jedges->_data[i]._node[0] = i;
    if (i<npo-1)
      jedges->_data[i]._node[1] = i+1;
    else
      jedges->_data[i]._node[1] = 0; // close polygon
  }

  /* Weighting grid                                                */
  res = 0.25 * (hmin + hmax);
  res = 2.0 * hmax;
  if (verbose) {
    printf("b %f %f\n",bmin*deg2m,bmax*deg2m);
    printf("h %f %f\n",hmin*deg2m,hmax*deg2m);
  }
  mce1 = (int)(fabs(mxlon - mnlon) / res);
  mce2 = (int)(fabs(mxlat - mnlat) / res);
  mce1 = mce2 = 100;
  mhfun = mce1 * mce2;
  if (verbose) printf("%d %d %d\n",mce1,mce2,mhfun);
  jigsaw_alloc_reals(hfx, mce1);
  jigsaw_alloc_reals(hfy, mce2);

  /* Fill in the x & y grids                                       */
  s = (mxlon - mnlon) / (double)(mce1 - 1);
  for (i = 0; i < mce1; i++)
    hfx->_data[i] = s * (double)i + mnlon;
  s = (mxlat - mnlat) / (double)(mce2 - 1);
  for (j = 0; j < mce2; j++)
    hfy->_data[j] = s * (double)j + mnlat;

  /* Make a gridded distance to poly array                         */
  jigsaw_alloc_reals(hfvals, mhfun);
  b = d_alloc_1d(mhfun);
  pmask = i_alloc_1d(mhfun);
  m = 0;
  for (i = 0; i < mce1; i++) {
    double xloc = hfx->_data[i];
    for (j = 0; j < mce2; j++) {
      double yloc = hfy->_data[j];
      b[m] = HUGE;
      b[m] = min(b[m], poly_dist(1, pl, bmin, xloc, yloc, &pmask[m]));
      m++;
    }
  }

  /* Get the minimum and maximum distances                         */
  /*
  bmin = HUGE;
  bmax = 0.0;
  for (m = 0; m < mhfun; m++) {
    bmin = min(bmin, b[m]);
    bmax = max(bmax, b[m]);
  }
  */

  /* Convert to the hfun function                                  */
  s = (bmin == bmax) ? 0.0 : (hmin - hmax) / (bmin - bmax);
  for (m = 0; m < mhfun; m++) {
    b[m] = max(min(bmax, b[m]), bmin);

    hfvals->_data[m] = bathyset(b[m], bmin, bmax, hmin, hmax, expf);
  }

  /* Set the perimeter to maximum resolution                       */
  for (m = 0; m < J_hfun->_edge2._size; m++) {
    i = J_hfun->_edge2._data[m]._node[0];
    hfvals->_data[i] = hmin;
    i = J_hfun->_edge2._data[m]._node[1];
    hfvals->_data[i] = hmin;
  }

  /* Set resolution inside polygons if required                    */
  for (i = 0; i < mhfun; i++) {
    if ((m = pmask[i]) >= 0) {
      hfvals->_data[i] = polyres[m];
    }
  }

  /* Set up some flags/parameters                                  */
  J_hfun->_flags = JIGSAW_EUCLIDEAN_GRID;
  J_jig._hfun_scal = JIGSAW_HFUN_ABSOLUTE ;
  if (stproj) {
    /* Leave in metres                                             */
    J_jig._hfun_hmin = hmin * deg2m;
    J_jig._hfun_hmax = hmax * deg2m;
  } else {
    J_jig._hfun_hmin = hmin;
    J_jig._hfun_hmax = hmax;
  }


  /* Call Jigsaw                                                   */
  J_jig._verbosity = 1;

  /* Stereographic projection                                      */
  if (stproj) {
    
    /* Calculate mid-points                                        */
    for (m = 0; m < npo; m++) {
      /* Figure out the centroid */
      xmid += jpoints->_data[m]._ppos[0];
      ymid += jpoints->_data[m]._ppos[1];
    }
    xmid /= npo; xmid = DEG2RAD(xmid);
    ymid /= npo; ymid = DEG2RAD(ymid);

    /* Perform forward stereographic projection of input data      */
    st_transform(&J_geom, ST3PROJ_FWD, xmid, ymid);

    /* Convert hfun grid to mesh, if needed                        */
    if (J_hfun->_flags == JIGSAW_EUCLIDEAN_GRID) {
      jigsaw_msh_t J_hfun_new;
      jigsaw_init_msh_t (&J_hfun_new);
      convert_jigsaw_grid_to_mesh(J_hfun, &J_hfun_new);
      /* Swap out                                                  */
      J_hfun = &J_hfun_new;
    }
    /* Projection of hfun                                          */
    st_transform(J_hfun, ST3PROJ_FWD, xmid, ymid);
  }

  jret = jigsaw_make_mesh(&J_jig, &J_geom, NULL, J_hfun, J_mesh);
  if (jret) hd_quit("Error calling JIGSAW\n");

  /* Perform inverse stereographic projection of output data       */
  if (stproj) st_transform(J_mesh, ST3PROJ_INV, xmid, ymid);

  /* Write to file                                                 */
  if (filef) {
    char key[MAXSTRLEN];
    FILE *ef, *hf;
    int n;
    mesh_ofile_name(params, key);
    sprintf(buf,"%s_h.txt", key);
    if ((ef = fopen(buf, "w")) != NULL) {
      /*
      fprintf(ef, "# Minimum distance to coast = %f m\n", bmin * deg2m);
      fprintf(ef, "# Maximum distance to coast = %f m\n", bmax * deg2m);
      */
      n = 0;
      for (i = 0; i < mce1; i++) {
	for (j = 0; j < mce2; j++) {
	  fprintf(ef, "%f %f %f %f\n",hfx->_data[i], hfy->_data[j],
		  hfvals->_data[n], b[n]);
	  n++;
	}
      }
      fclose(ef);
    }
    sprintf(buf,"%s_hfun.msh", key);
    if ((hf = fopen(buf, "w")) != NULL) {
      fprintf(hf, "# .msh geometry file\n");
      fprintf(hf, "MSHID=2;EUCLIDEAN-GRID\n");
      fprintf(hf, "NDIMS=2\n");
      fprintf(hf, "coord=1;%d\n",mce1);
      for (i = 0; i < mce1; i++)
	fprintf(hf, "%f\n",hfx->_data[i]);
      fprintf(hf, "coord=2;%d\n",mce2);
      for (j = 0; j < nce2; j++)
	fprintf(hf, "%f\n",hfy->_data[j]);
      
      fprintf(hf, "value=%d;1\n",mce1*mce2);
      n = 0;
      for (i = 0; i < mce1; i++)
	for (j = 0; j < mce2; j++) {
	  fprintf(hf, "%f\n",hfvals->_data[n++]);
	}
      fclose(hf);
    }
  }

  /* Interior coordinate                                             */
  params->mlon = x00;
  params->mlat = y00;

  /* Cleanup                                                         */
  jigsaw_free_msh_t(&J_geom);

  d_free_1d(xw);
  d_free_1d(yw);
  d_free_1d(xb);
  d_free_1d(yb);
  d_free_1d(b); 
  d_free_1d(polyres);
  i_free_1d(pmask);
  poly_destroy(pl[0]);
  free((poly_t **)pl);
}

/* END create_ellipse()                                              */
/*-------------------------------------------------------------------*/


#define I_NC    1
#define I_BTY   2
#define I_USR   4
#define I_MSH   8
#define I_MESH  16
#define I_GRID  32

/*-------------------------------------------------------------------*/
/* Reads a bathymetry file and creates a hfun file                   */
/*-------------------------------------------------------------------*/
void hfun_from_bathy(parameters_t *params, char *fname, coamsh_t *cm, jigsaw_msh_t *J_hfun, int mode)
{
  FILE *fp = params->prmfd, *ef, *bf, *hf;
  char buf[MAXSTRLEN], key[MAXSTRLEN];
  char i_rule[MAXSTRLEN];
  char vname[MAXSTRLEN];
  GRID_SPECS *gs = NULL;
  int nbath;
  int **nei;
  double *x, *y, *b, *r;
  double *hfun, xloc, yloc;
  int n, m, i, j, nhfun, intype, bmf = 1;
  int imeth = 0;
  int verbose = 0;
  int filef = 1;
  int gridf = 0;
  int orf = 0;
  int sn, smooth = 0;
  int nce1, nce2;
  int fid;
  int ncerr;
  int dims;
  int ugrid = 0;
  double stime;
  double bmin, bmax, bnin, bnax, hmin = 0.0, hmax = 0.0, s;
  double *exlon, *exlat, *exrad, *exbth;
  double deg2m = 60.0 * 1852.0;
  double expf = 0.0;
  double cres = 0.0;
  delaunay *d;
  jigsaw_REALS_array_t *hfvals = &J_hfun->_value;
  jigsaw_REALS_array_t *hfx    = &J_hfun->_xgrid;
  jigsaw_REALS_array_t *hfy    = &J_hfun->_ygrid;
  int perimf = 1; /* Set hfun perimeter to maximum resolution        */
  int npass = 0;
  double *hmx, *hmn, *bmx, *bmn, *exf;
  double gws = 9.81;
  int nord;

  filef = (params->meshinfo) ? 1 : 0;
  /*-----------------------------------------------------------------*/
  /* Get the interpolation method                                    */
  /*if(endswith(fname,".nc")) {*/
  if(mode & H_NC) {
    size_t ncx, ncy;
    /* Open the dump file for reading                                */
    if ((ncerr = nc_open(fname, NC_NOWRITE, &fid)) != NC_NOERR) {
      printf("Can't find input bathymetry file %s\n", fname);
      hd_quit((char *)nc_strerror(ncerr));
    }

    /* Get dimensions                                                */
    sprintf(vname, "%c", '\0');
    prm_read_char(fp, "HFUN_VAR", vname);
    stime = 0.0;
    if (prm_read_char(fp, "HFUN_TIME", buf)) {
      char *fields[MAXSTRLEN * MAXNUMARGS];
      i = parseline(buf, fields, MAXNUMARGS);
      stime = atoi(fields[0]);
    }
    i = 0;
    while (bathy_dims[i][0] != NULL) {
      if (nc_inq_dimlen(fid, ncw_dim_id(fid, bathy_dims[i][2]), &ncy) == 0 &&
	  nc_inq_dimlen(fid, ncw_dim_id(fid, bathy_dims[i][1]), &ncx) == 0) {
	intype = i;
	break;
      }
      i++;
    }
    if (intype < 0)      
      hd_quit("hfun_from_bathy: Can't find batyhmetry attributes in file %s.\n", fname);

    if (strcmp(bathy_dims[intype][5], "nc_bathy") == 0)
      imeth = (I_NC|I_GRID);
    else
      imeth = (I_NC|I_MESH);
  } else if (mode & H_BTY) {
    /*else if (endswith(fname,".bty")) {*/
    imeth = (I_BTY|I_MESH);
  } else if (mode & H_MSH) {
    /*else if (endswith(fname,".msh")) {*/
    imeth = I_MSH;
  } else {
    double d1, d2, d3;
    imeth = (I_USR|I_MESH);
    prm_skip_to_end_of_key(fp, "HFUN_BATHY_FILE");
    i = atoi(fname);    
    if (i == 0 || fscanf(fp, "%lf %lf %lf", &d1, &d2, &d3) != 3)
      cres = atof(fname) / deg2m;
  }

  /*-----------------------------------------------------------------*/
  /* Define and read parameters                                      */
  /* Note: bmin and bmax are negative for wet cells.                 */
  if (prm_read_int(fp, "NHFUN", &npass)) {
    bnin = bmax = hmax = -HUGE;
    bnax = bmin = hmin = HUGE;
    bmn = d_alloc_1d(npass);
    bmx = d_alloc_1d(npass);
    hmn = d_alloc_1d(npass);
    hmx = d_alloc_1d(npass);
    exf = d_alloc_1d(npass);
    for (n = 0; n < npass; n++) {
      sprintf(buf, "HMIN%d", n);
      prm_read_double(fp, buf, &hmn[n]);
      hmn[n] /= deg2m;
      sprintf(buf, "HMAX%d", n);
      prm_read_double(fp, buf, &hmx[n]);
      hmx[n] /= deg2m;
      sprintf(buf, "BMIN%d", n);
      prm_read_double(fp, buf, &bmn[n]);
      sprintf(buf, "BMAX%d", n);
      prm_read_double(fp, buf, &bmx[n]);
      sprintf(buf, "TYPE%d", n);
      prm_read_double(fp, buf, &exf[n]);
      /* If bathymetries from file are < 0, then bmax is the deepest */
      /* depth and bmin the shallowest such that bmax < bmin. In     */
      /* case we get the bmax from min(bmax, bmn[n]).                */
      if (bmn[n] < 0) {
	bnin = max(bnin, bmn[n]);
	bmin = bnin;
      } else
	bmin = min(bmin, bmn[n]);
      if (bmx[n] < 0) {
	bnax = min(bnax, bmx[n]);
	bmax = bnax;
      } else
	bmax = max(bmax, bmx[n]);
      hmax = max(hmax, hmx[n]);
      hmin = min(hmin, hmn[n]);
    }
    cm->hfun_min = hmin * deg2m;
    cm->hfun_max = hmax * deg2m;

  } else {
    hmin = cm->hfun_min / deg2m;
    hmax = cm->hfun_max / deg2m;
    bmin = bmax = 0.0;
    prm_read_double(fp, "HFUN_BMIN", &bmin);
    prm_read_double(fp, "HFUN_TYPE", &expf);
    if (prm_read_double(fp, "HFUN_BMAX", &bmax)) bmf = 1;
    /*
    if (strlen(vname) == 0 && bmin > 0.0 && bmax > 0.0 && bmin < bmax) {
      bmin *= -1.0;
      bmax *= -1.0;
    }
    */
  }
  prm_read_int(fp, "HFUN_SMOOTH", &smooth);

  /* Get the interpolation rule                                      */
  strcpy(i_rule, "linear");
  prm_read_char(params->prmfd, "HFUN_INTERP_RULE", i_rule);
 
  /* Set on a regular grid if required                               */
  if (prm_read_char(params->prmfd, "HFUN_GRID", buf)) {
    if (imeth & (I_NC|I_BTY|I_USR) && is_true(buf)) {
      if (imeth & I_MESH) imeth &= ~I_MESH;
      imeth |= I_GRID;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set up the grid to interpolate onto.                            */
  if (imeth & I_GRID) {
    /* Make a regular grid in the bounding box at the maximum        */
    /* resolution.                                                   */
    double res = (cres) ? cres : 0.25 * (hmin + hmax);
    /*double res = 0.05;*/
    nce1 = (int)(fabs(cm->belon - cm->bslon) / res);
    nce2 = (int)(fabs(cm->belat - cm->bslat) / res);
    nhfun = nce1 * nce2;
    jigsaw_alloc_reals(hfx, nce1);
    jigsaw_alloc_reals(hfy, nce2);
    /* Fill in the x & y grids                                       */
    s = (cm->belon - cm->bslon) / (double)(nce1 - 1);
    for (i = 0; i < nce1; i++)
      hfx->_data[i] = s * (double)i + cm->bslon;
    s = (cm->belat - cm->bslat) / (double)(nce2 - 1);
    for (j = 0; j < nce2; j++)
      hfy->_data[j] = s * (double)j + cm->bslat;
    jigsaw_alloc_reals(hfvals, nhfun);
  } else if (imeth & I_MESH) {
    /* Get a triangulation based on the mesh perimeter               */
    /*xy_to_d(d, cm->np, cm->x, cm->y);*/
    if (imeth & I_BTY)
      hd_quit("Can't create a mesh with .bty input. Use HFUN_GRID.\n");
    if (imeth & I_USR)
      hd_quit("Can't create a mesh with user input. Use HFUN_GRID.\n");

    create_bounded_mesh(cm->np, cm->x, cm->y, hmin, hmax, J_hfun, filef);
    nhfun = J_hfun->_vert2._size;
    jigsaw_alloc_reals(hfvals, nhfun);
  } else {
    FILE *hp;
    if ((hp = fopen(fname, "r")) == NULL)
      hd_quit("Can't open hfun file %s\n", fname);
    fgets(buf, MAXSTRLEN, hp);  /* Header comment */
    fgets(buf, MAXSTRLEN, hp);  /* MESHID */
    if (sscanf(buf, "MSHID=%d;EUCLIDEAN-GRID", &i) == 1) {
      fgets(buf, MAXSTRLEN, hp);  /* NDIMS  */
      if (sscanf(buf, "NDIMS=%d", &i) != 1)
	hd_quit("Can't read NDIMS in file %s\n", fname);

      if (!(prm_skip_to_end_of_tok(hp, "COORD=1", ";", buf)))
	hd_quit("Can't read COORD1 in file %s\n", fname);
      nce1 = atoi(buf);
      jigsaw_alloc_reals(hfx, nce1);
      for (i = 0; i < nce1; i++) {
	fscanf(hp, "%lf", &hfx->_data[i]);
      }

      if (!(prm_skip_to_end_of_tok(hp, "COORD=2", ";", buf)))
	hd_quit("Can't read COORD2 in file %s\n", fname);
      nce2 = atoi(buf);
      jigsaw_alloc_reals(hfy, nce2);
      for (i = 0; i < nce2; i++) {
	fscanf(hp, "%lf", &hfy->_data[i]);
      }

      if (!(prm_skip_to_end_of_tok(hp, "VALUE", "=", buf)))
	hd_quit("Can't read VALUE in file %s\n", fname);
      nhfun = atoi(buf);
      jigsaw_alloc_reals(hfvals, nhfun);
      for (i = 0; i < nhfun; i++) {
	fscanf(hp, "%lf", &hfvals->_data[i]);
      }
      J_hfun->_flags = JIGSAW_EUCLIDEAN_GRID;
    } else if (sscanf(buf, "MSHID=%d;EUCLIDEAN-MESH", &i) == 1) {
      jigsaw_VERT2_array_t *jpoints = &J_hfun->_vert2;
      fgets(buf, MAXSTRLEN, hp);  /* NDIMS  */
      if (sscanf(buf, "NDIMS=%d", &i) != 1)
	hd_quit("Can't read NDIMS in file %s\n", fname);
      fgets(buf, MAXSTRLEN, hp);  /* POINT  */
      if (sscanf(buf, "POINT=%d", &nhfun) != 1)
	hd_quit("Can't read POINT in file %s\n", fname);
      jigsaw_alloc_vert2(jpoints, nhfun);
      for (i = 0; i < nhfun; i++) {
	fscanf(hp, "%lf;%lf;0\n",&jpoints->_data[n]._ppos[0],
	       &jpoints->_data[n]._ppos[1]);
      }
      J_hfun->_flags = JIGSAW_EUCLIDEAN_MESH;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Read bathymetry over-ride values                                */
  if (prm_read_int(fp, "HFUN_OVERRIDE", &nord)) {
    exlon = d_alloc_1d(nord);
    exlat = d_alloc_1d(nord);
    exrad = d_alloc_1d(nord);
    exbth = d_alloc_1d(nord);
    for (n = 0; n < nord; n++) {
      if (fscanf(fp, "%lf %lf %lf %lf", &exlon[n], &exlat[n], &exrad[n], &exbth[n]) != 4)
	hd_quit("hfun_from_bathy: Format for HFUN_OVERRIDE is 'lon lat radius value'.\n");
      if (exbth[n] > 0.0) exbth[n] /= deg2m;
    }
    orf = 1;
  }

  /*-----------------------------------------------------------------*/
  /* Make a gridded bathymetry                                       */
  if(imeth & I_NC) {
    /* Interpolation from structured input using ts_read().          */
    size_t start[4];
    size_t count[4];
    timeseries_t *ts = NULL;
    int idb;
    int dimf;

    /* Initialize the bathymetry file                                */
    ts = (timeseries_t *)malloc(sizeof(timeseries_t));
    if (ts == NULL)
      hd_quit("hfun_from_bathy: No memory available.\n");
    memset(ts, 0, sizeof(timeseries_t));
    
    /* Read the bathymetry file                                      */
    ts_read(fname, ts);
    ugrid = df_is_ugrid(ts->df);

    if (strlen(vname)) {
      strcpy(key, vname);
      nc_inq_varndims(fid, ncw_var_id(fid, vname), &dims);
    } else
      strcpy(key, fv_get_varname(fname, bathy_dims[intype][0], buf));

    if ((idb = ts_get_index(ts, key)) == -1)
      hd_quit("hfun_from_bathy: Can't find variable %s in file %s\n", 
	      key, fname);

    b = d_alloc_1d(nhfun);
    if (ugrid)
      dimf = (dims == 3) ? 3 : 2;
    else
      dimf = (dims == 4) ? 3 : 2;

    if (imeth & I_GRID) {
      /* Interpolate bathymetry onto the grid                        */
      n = 0;
      for (i = 0; i < nce1; i++) {
	for (j = 0; j < nce2; j++) {
	  xloc = hfx->_data[i];
	  yloc = hfy->_data[j];
	  if (dimf == 3)
	    b[n] = ts_eval_xyz(ts, idb, stime, xloc, yloc, 0.0);
	  else
	    b[n] = ts_eval_xy(ts, idb, stime, xloc, yloc);
	  if (!bmf && b[n] < 0.0 && b[n] < bmax) bmax = b[n];
	  if (!bmf && b[n] > 0.0 && b[n] > bmax) bmax = b[n];
	  if (verbose) printf("%d %f %f : %f %f\n",n, xloc, yloc, b[n], bmax);
	  n++;
	}
      }
      J_hfun->_flags = JIGSAW_EUCLIDEAN_GRID;
    } else {
      /* Interpolate the bathymetry values onto the triangulation and  */
      /* find the maximum depth if required.                           */
      for (n = 0; n < nhfun; n++) {
	xloc = J_hfun->_vert2._data[n]._ppos[0];
	yloc = J_hfun->_vert2._data[n]._ppos[1];
	if (dims == 4)
	  b[n] = ts_eval_xyz(ts, idb, stime, xloc, yloc, 0.0);
	else
	  b[n] = ts_eval_xy(ts, idb, stime, xloc, yloc);
	if (!bmf && b[n] < bmax) bmax = b[n];
	if (verbose) printf("%d %f %f : %f %f\n",n, xloc, yloc, b[n], bmax);
      }
      J_hfun->_flags = JIGSAW_EUCLIDEAN_MESH;
    }
    ts_free((timeseries_t*)ts);
    free(ts);
  } else if(imeth & I_BTY) {

    /* Read the bathymetry values                                    */
    if ((bf = fopen(fname, "r")) == NULL)
      hd_quit("hfun_from_bathy: Can't open bathymetry file %s.\n", fname);
    nbath = 0;
    while (fgets(buf, MAXSTRLEN, bf) != NULL) {
      nbath++;
    }
    rewind(bf);
    x = d_alloc_1d(nbath);
    y = d_alloc_1d(nbath);
    r = d_alloc_1d(nbath);
    n = 0;
    while (fgets(buf, MAXSTRLEN, bf) != NULL) {
      char *fields[MAXSTRLEN * MAXNUMARGS];
      i = parseline(buf, fields, MAXNUMARGS);
      x[n] = atof(fields[0]);
      y[n] = atof(fields[1]);
      r[n] = atof(fields[2]);
      if (r[n] > 0) r[n] *= -1.0;
      if (!bmf && r[n] < bmax) bmax = r[n];
      if (verbose) printf("%d %f %f : %f %f\n",n, x[n], y[n], r[n], bmax);
      n++;
    }
    fclose(bf);

    /* Set up the interpolation triangulation                        */
    gs = grid_interp_init(x, y, r, nbath, i_rule);

    /* Interpolte bathymetry onto the mesh                           */
    b = d_alloc_1d(nhfun+1);
    if (imeth & I_GRID) {
      /* Interpolate bathymetry onto the grid                        */
      n = 0;
      for (i = 0; i < nce1; i++) {
	for (j = 0; j < nce2; j++) {
	  xloc = hfx->_data[i];
	  yloc = hfy->_data[j];
	  b[n] = grid_interp_on_point(gs, xloc, yloc);
	  if (!bmf && b[n] < bmax) bmax = b[n];
	  n++;
	}
      }
      grid_specs_destroy(gs);
      J_hfun->_flags = JIGSAW_EUCLIDEAN_GRID;
    } else {
      for (n = 0; n < nhfun; n++) {
	xloc = J_hfun->_vert2._data[n]._ppos[0];
	yloc = J_hfun->_vert2._data[n]._ppos[1];
	b[n] = grid_interp_on_point(gs, xloc, yloc);
	b[n] = min(max(bmax, b[n]), bmin);
      }
      J_hfun->_flags = JIGSAW_EUCLIDEAN_MESH;
      grid_specs_destroy(gs);
      d_free_1d(x);
      d_free_1d(y);
    }
    d_free_1d(r);
  } else {
    int constf = 0;
    double d1, d2, d3;
    if (!(imeth & I_MSH)) {
      prm_skip_to_end_of_key(fp, "HFUN_BATHY_FILE");
      nbath = atoi(fname);    
      if (nbath == 0 || fscanf(fp, "%lf %lf %lf", &d1, &d2, &d3) != 3) {
	constf = 1;
	d1 = atof(fname) / deg2m;
	if (imeth & I_GRID) {
	  n = 0;
	  for (i = 0; i < nce1; i++) {
	    for (j = 0; j < nce2; j++) {
	      hfvals->_data[n++] = d1;
	    }
	  }
	} else {
	  for (n = 0; n < nhfun; n++) {
	    hfvals->_data[n] = d1;
	  }
	}
	J_hfun->_flags = JIGSAW_EUCLIDEAN_GRID;
	if (orf) {
	  for (n = 0; n < nord; n++) exbth[n] = fabs(exbth[n]);
	}
      } else {
	prm_skip_to_end_of_key(fp, "HFUN_BATHY_FILE");
	x = d_alloc_1d(nbath);
	y = d_alloc_1d(nbath);
	b = d_alloc_1d(nbath);
	for (i = 0; i < nbath; i++) {
	  if ((fscanf(fp, "%lf %lf %lf", &x[i], &y[i], &b[i])) != 3)
	    hd_quit("hfun_from_bathy: Incorrect resolution specification.\n");
	}
	
	/* Interpolate from a triangulation                              */
	gs = grid_interp_init(x, y, b, nbath, i_rule);

	n = 0;
	for (i = 0; i < nce1; i++) {
	  for (j = 0; j < nce2; j++) {
	    xloc = hfx->_data[i];
	    yloc = hfy->_data[j];
	    hfvals->_data[n] = grid_interp_on_point(gs, xloc, yloc);
	    if (verbose) printf("%d %f : %f %f\n",n, hfvals->_data[n], xloc, yloc);
	    n++;
	  }
	}
	/*
	for (n = 0; n < nhfun; n++) {
	  xloc = J_hfun->_vert2._data[n]._ppos[0];
	  yloc = J_hfun->_vert2._data[n]._ppos[1];
	  hfvals->_data[n] = grid_interp_on_point(gs, xloc, yloc);
	  if (verbose) printf("%d %f : %f %f\n",n, hfvals->_data[n], xloc, yloc);
	}
	*/
	J_hfun->_flags = JIGSAW_EUCLIDEAN_GRID;
	grid_specs_destroy(gs);
	d_free_1d(x);
	d_free_1d(y);
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Gravity wave speed conversion if required.                      */
  if (mode & H_GWS && imeth & (I_NC|I_BTY)) {
    for (n = 0; n < nhfun; n++) {
      b[n] = sqrt(gws * fabs(b[n]));
    }
  }

  /*-----------------------------------------------------------------*/
  /* Gradient conversion if required.                                */
  if (mode & H_GRAD && imeth & (I_NC|I_BTY)) {
    double **bg;

    if (!(imeth & I_GRID))
      hd_quit("Can't create a bathy gradient: Use HFUN_GRID.\n");

    n = 0;
    bg = d_alloc_2d(nce1, nce2);
    /* Map the bathy to Cartesian coordinates                        */
    for (i = 0; i < nce1; i++) {
      for (j = 0; j < nce2; j++) {
	bg[j][i] = b[n++];
      }
    }
    /* Get the gradient                                              */
    n = 0;
    for (i = 0; i < nce1; i++) {
      for (j = 0; j < nce2; j++) {
	double dx, dy;
	int im = (i == 0) ? i : i-1;
	int ip = (i == nce1 - 1) ? i : i+1;
	int jm = (j == 0) ? j : j-1;
	int jp = (j == nce2 - 1) ? j : j+1;
	dx = (bg[ip] - bg[im]) / (deg2m * (hfx->_data[ip] - hfx->_data[im]));
	dy = (bg[jp] - bg[jm]) / (deg2m * (hfy->_data[jp] - hfy->_data[jm]));
	b[n++] = 0.5 * fabs(dx + dy);
      }
    }
    d_free_2d(bg);
  }

  /*-----------------------------------------------------------------*/
  /* Over-ride bathymetry values if required                         */
  if (orf) {

    if (imeth & I_GRID) {
      double *val = (imeth & I_USR) ? hfvals->_data : b;
      double smin[nord];
      int imin[nord], jmin[nord], nmin[nord];
      n = 0;
      for (m = 0; m < nord; m++) {
	smin[m] = HUGE;
	nmin[m] = -1;
      }
      for (i = 0; i < nce1; i++) {
	for (j = 0; j < nce2; j++) {
	  for (m = 0; m < nord; m++) {
	    if (exbth[m] > 0.0) continue;
	    xloc = exlon[m] - hfx->_data[i];
	    yloc = exlat[m] - hfy->_data[j];
	    s = sqrt(xloc * xloc + yloc * yloc) * deg2m;
	    if (s < exrad[m]) {
	      val[n] = (val[n] - exbth[m]) * s / exrad[m] + exbth[m];
	      if (s < smin[m]) {
		smin[m] = s;
		imin[m] = i;
		jmin[m] = j;
		nmin[m] = n;
	      }
	    }
	  }
	  n++;
	}
      }
      for (m = 0; m < nord; m++) {
	if (nmin[m] >= 0 && exbth[m] <= 0.0) {
	  val[nmin[m]] = exbth[m];
	  hd_warn("HFUN_OVERRIDE%d minimum dist @ %f %f\n", m,
		  hfx->_data[imin[m]],hfy->_data[jmin[m]]);
	}
      }
    } else {
      for (n = 0; n < nhfun; n++) {
	for (m = 0; m < nord; m++) {
	  if (exbth[m] > 0.0) continue;
	  xloc = exlon[m] - J_hfun->_vert2._data[n]._ppos[0];
	  yloc = exlat[m] - J_hfun->_vert2._data[n]._ppos[1];
	  s = sqrt(xloc * xloc + yloc * yloc) * deg2m;
	  if (s < exrad[m]) b[n] = (b[n] - exbth[m]) * s / exrad[m] + exbth[m];
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Convert to the hfun function                                    */
  if (imeth & (I_NC|I_BTY)) {
    int dof = 0;

    x = d_alloc_1d(nhfun);
    y = d_alloc_1d(nhfun);
    n = 0;
    for (i = 0; i < nce1; i++) {
      for (j = 0; j < nce2; j++) {
	x[n] = hfx->_data[i];
	y[n++] = hfy->_data[j];
      }
    }

    for (n = 0; n < nhfun; n++) {

      /* Replace NaNs with nearest valid value                       */
      bnax = HUGE;
      if (isnan(b[n])) {
	for (m = 0; m < nhfun; m++) {
	  if (!isnan(b[m])) {
	    s = sqrt((x[n]-x[m])*(x[n]-x[m])+(y[n]-y[m])*(y[n]-y[m]));
	    if (s < bnax) {
	      b[n] = b[m];
	      bnax = s;
	    }
	  }
	}
      }

      /* Limit the bathymetry. Nore for negative bathymetries the    */
      /* maximum is the deepest (i.e. most negative).                */
      if (b[n] < 0.0)
	b[n] = min(max(bmax, b[n]), bmin);
      else
	b[n] = max(min(bmax, b[n]), bmin);

      if (npass) {
	hfvals->_data[n] = 0.0;
	for (m = 0; m < npass; m++) {
	  if ((b[n] < 0.0 && (b[n] >= bmx[m] && b[n] <= bmn[m])) ||
	      (b[n] > 0.0 && (b[n] <= bmx[m] && b[n] >= bmn[m])))
	    hfvals->_data[n] = bathyset(b[n], bmn[m], bmx[m], 
					hmn[m], hmx[m], exf[m]);
	}
	if (!hfvals->_data[n]) {
	  if (b[n] < 0.0) {
	    if (b[n] < bmax) hfvals->_data[n] = hmax;
	    if (b[n] > bmin) hfvals->_data[n] = hmin;
	  } else {
	    if (b[n] > bmax) hfvals->_data[n] = hmax;
	    if (b[n] < bmin) hfvals->_data[n] = hmin;
	  }
	}
      } else 
	hfvals->_data[n] = bathyset(b[n], bmin, bmax, hmin, hmax, expf);

      /* Sanity check: JIGSAW only acceps values > 0                 */
      if (hfvals->_data[n] <= 0.0) hfvals->_data[n] = hmin;

      /*printf("%d %f %f min=%f max=%f\n",n, hfvals->_data[n]*deg2m, b[n],bmin,bmax);*/
      if (dof) {
      if (expf > 0.0) {
	/* Exponential                                               */
	hfvals->_data[n] = (hmin - hmax) * exp(b[n] / expf) + 
	  (hmax - hmin * exp(-fabs(bmax) / expf));
      } else if (expf == 0) {
	/* Linear                                                    */
	hfvals->_data[n] = ((b[n] - bmax) * s + hmax);
      } else {
	double bn = fabs(bmin);
	double bx = fabs(bmax);
	double ba = fabs(b[n]);
	double dd = bx - bn;
	if (ba < bn)
	  hfvals->_data[n] = hmin;
	else if (ba > bx)
	  hfvals->_data[n] = hmax;
	else
	  hfvals->_data[n] = 0.5 * ((hmin - hmax) * cos(ba * PI / dd - bn * PI / dd) + (hmin + hmax));
      }
      if (verbose) printf("%d %f %f\n",n, hfvals->_data[n], b[n]);
    }
    }
    d_free_1d(x);
    d_free_1d(y);
  }

  /*-----------------------------------------------------------------*/
  /* Set the perimeter to maximum resolution                         */
  if (perimf && imeth & I_MESH) {
    for (n = 0; n < J_hfun->_edge2._size; n++) {
      i = J_hfun->_edge2._data[n]._node[0];
      hfvals->_data[i] = hmin;
      i = J_hfun->_edge2._data[n]._node[1];
      hfvals->_data[i] = hmin;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Over-ride resolution values if required                         */
  if (orf) {
    double *val = hfvals->_data;
    if (imeth & I_GRID) {
      double smin[nord];
      int imin[nord], jmin[nord], nmin[nord];
      n = 0;
      for (m = 0; m < nord; m++) {
	smin[m] = HUGE;
	nmin[m] = -1;
      }
      for (i = 0; i < nce1; i++) {
	for (j = 0; j < nce2; j++) {
	  for (m = 0; m < nord; m++) {
	    if (exbth[m] <= 0.0) continue;
	    xloc = exlon[m] - hfx->_data[i];
	    yloc = exlat[m] - hfy->_data[j];
	    s = sqrt(xloc * xloc + yloc * yloc) * deg2m;
	    if (s < exrad[m]) {
	      val[n] = (val[n] - exbth[m]) * s / exrad[m] + exbth[m];
	      if (s < smin[m]) {
		smin[m] = s;
		imin[m] = i;
		jmin[m] = j;
		nmin[m] = n;
	      }
	    }
	  }
	  n++;
	}
      }
      for (m = 0; m < nord; m++) {
	if (nmin[m] >= 0 && exbth[m] > 0.0) {
	  val[nmin[m]] = exbth[m];
	  hmin = min(hmin, exbth[m]);
	  hmax = max(hmax, exbth[m]);
	  cm->hfun_min = hmin * deg2m;
	  cm->hfun_max = hmax * deg2m;
	  hd_warn("HFUN_OVERRIDER%d minimum dist @ %f %f\n", m,
		  hfx->_data[imin[m]],hfy->_data[jmin[m]]);
	}
      }
    } else {
      for (n = 0; n < nhfun; n++) {
	for (m = 0; m < nord; m++) {
	  if (exbth[m] <= 0.0) continue;
	  xloc = exlon[m] - J_hfun->_vert2._data[n]._ppos[0];
	  yloc = exlat[m] - J_hfun->_vert2._data[n]._ppos[1];
	  s = sqrt(xloc * xloc + yloc * yloc) * deg2m;
	  if (s < exrad[m]) val[n] = (val[n] - exbth[m]) * s / exrad[m] + exbth[m];
	}
      }
    }
    if (exlon) d_free_1d(exlon);
    if (exlat) d_free_1d(exlat);
    if (exrad) d_free_1d(exrad);
    if (exbth) d_free_1d(exbth);
  }

  /*-----------------------------------------------------------------*/
  /* Smooth if required                                              */
  if (smooth) {
    if (imeth & I_MESH) {
      int nv2 = J_hfun->_vert2._size;
      int *nv, **eiv;
      /*neighbour_finder_j(J_hfun, &nei);*/
      nv = i_alloc_1d(nv2);
      memset(nv, 0, nv2 * sizeof(int));
      for (n = 0; n < J_hfun->_tria3._size; n++) {
	for (j = 0; j < 3; j++) {
	  i = J_hfun->_tria3._data[n]._node[j];
	  nv[i]++;
	}
      }
      m = 0;
      for (n = 0; n < nv2; n++)
	if (nv[n] > m) m = nv[n];
      eiv = i_alloc_2d(nv2, m);
      memset(nv, 0, nv2 * sizeof(int));
      for (n = 0; n < J_hfun->_tria3._size; n++) {
	for (j = 0; j < 3; j++) {
	  i = J_hfun->_tria3._data[n]._node[j];
	  eiv[nv[i]][i] = n;
	  nv[i]++;
	}
      }
      r = d_alloc_1d(nhfun);
      for (sn = 0; sn < smooth; sn++) {
	/* Loop over all vertices                                    */
	for (m = 0; m < nv2; m++) {
	  r[m] = hfvals->_data[m];
	  s = 1.0;
	  /* Loop over all neighbouring triangles                    */
	  for (j = 0; j < nv[m]; j++) {
	    n = eiv[j][m];
	    for (i = 0; i < 3; i++) {
	      int nn = J_hfun->_tria3._data[n]._node[i];
	      if (nn != m && nn < nhfun) {
		r[m] += hfvals->_data[nn];
	      s += 1.0;
	      }
	    }
	  }
	  r[m] /= s;
	}
	memcpy(hfvals->_data, r, nhfun * sizeof(double));
      }
      d_free_1d(r);
      i_free_1d(nv);
      i_free_2d(eiv);
    }
    if (imeth & I_GRID) {
      int ii, jj, **ij2n;

      /* Get a mapping from ij to hfvals indices                     */
      ij2n = i_alloc_2d(nce1, nce2);
      n = 0;
      for (i = 0; i < nce1; i++) {
	for (j = 0; j < nce2; j++) {
	  ij2n[j][i] = n;
	  n++;
	}
      }

      /* Smooth                                                      */
      for (sn = 0; sn < smooth; sn++) {
	/* Loop over all vertices                                    */
	r = d_alloc_1d(nhfun);
	for (i = 0; i < nce1; i++) {
	  int i1 = max(i-1, 0), i2 = min(i+1, nce1-1);
	  for (j = 0; j < nce2; j++) {
	    int j1 = max(j-1, 0), j2 = min(j+1, nce2-1);
	    s = 0.0;
	    n = ij2n[j][i];
	    r[n] = 0.0;
	    for (ii = i1; ii <= i2; ii++) {
	      for (jj = j1; jj <= j2; jj++) {
		ii = max(0, min(nce1-1, ii));
		jj = max(0, min(nce2-1, jj));
		m = ij2n[jj][ii];
		r[n] += hfvals->_data[m];
		s += 1.0;
	      }
	    }
	    r[n] /= s;
	  }
	}
	memcpy(hfvals->_data, r, nhfun * sizeof(double));
      }
      d_free_1d(r);
      i_free_2d(ij2n);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Write to file                                                   */
  if (filef) {
    mesh_ofile_name(params, key);
    sprintf(buf,"%s_h.txt", key);
    if ((ef = fopen(buf, "w")) != NULL) {
      if (imeth & I_MESH) { 
	for (n = 0; n < nhfun; n++) {
	  if (imeth & I_USR)
	    fprintf(ef, "%f %f %f\n",J_hfun->_vert2._data[n]._ppos[0], 
		    J_hfun->_vert2._data[n]._ppos[1], 
		    hfvals->_data[n]);
	  else
	    fprintf(ef, "%f %f %f %f\n",J_hfun->_vert2._data[n]._ppos[0], 
		    J_hfun->_vert2._data[n]._ppos[1], 
		    hfvals->_data[n], b[n]);
	}
      } else {
	n = 0;
	for (i = 0; i < nce1; i++) {
	  for (j = 0; j < nce2; j++) {
	    if (imeth & (I_USR|I_MSH))
	      fprintf(ef, "%f %f %f\n",hfx->_data[i], hfy->_data[j],
		      hfvals->_data[n]);
	    else
	      fprintf(ef, "%f %f %f %f\n",hfx->_data[i], hfy->_data[j],
		      hfvals->_data[n], b[n]);
	    n++;
	  }
	}
      }
      fclose(ef);
    }
    sprintf(buf,"%s_hfun.msh", key);
    if ((hf = fopen(buf, "w")) != NULL) {
      fprintf(hf, "# .msh geometry file\n");
      if (imeth & I_MESH) {
	/* EUCLIDEAN-MESH */
	fprintf(hf, "MSHID=2;EUCLIDEAN-MESH\n");
      } else {
	/* EUCLIDEAN-GRID */
	fprintf(hf, "MSHID=2;EUCLIDEAN-GRID\n");
      }
      fprintf(hf, "NDIMS=2\n");
      if (imeth & I_MESH) {
	fprintf(hf, "POINT=%d\n",nhfun);
	for (n = 0; n < nhfun; n++)
	  fprintf(hf, "%f;%f;0\n",J_hfun->_vert2._data[n]._ppos[0],
		  J_hfun->_vert2._data[n]._ppos[1]);
      } else {
	fprintf(hf, "coord=1;%d\n",nce1);
	for (i = 0; i < nce1; i++)
	  fprintf(hf, "%f\n",hfx->_data[i]);
	fprintf(hf, "coord=2;%d\n",nce2);
	for (j = 0; j < nce2; j++)
	  fprintf(hf, "%f\n",hfy->_data[j]);
      }
      if (imeth & I_MESH) {
	fprintf(hf, "value=%d;1\n",nhfun);
	for (n = 0; n < nhfun; n++)
	  fprintf(hf, "%f\n",hfvals->_data[n]);
      } else {
	fprintf(hf, "value=%d;1\n",nce1*nce2);
	n = 0;
	for (i = 0; i < nce1; i++)
	  for (j = 0; j < nce2; j++)
	    fprintf(hf, "%f\n",hfvals->_data[n++]);
      }
      fclose(hf);
    }
  }
  if (b) d_free_1d(b); 
  if (DEBUG("init_m"))
    dlog("init_m", "Bathymetric weighting function computed OK\n");
}

/* END hfun_from_bathy()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computes the weighting function from bathymetry                   */
/*-------------------------------------------------------------------*/
double bathyset(double b, 
		double bmin,
		double bmax,
		double hmin,
		double hmax,
		double expf
		)
{
  double ret;
  double s = (bmin == bmax) ? 0.0 : (hmin - hmax) / (bmin - bmax);

  if (expf > 0.0) {
    /* Exponential                                                   */
    ret = (hmin - hmax) * exp(-fabs(b) / expf) + 
      (hmax - hmin * exp(-fabs(bmax) / expf));
  } else if (expf == 0) {
    /* Linear                                                        */
    ret = ((b - bmax) * s + hmax);
  } else {
    double bn = fabs(bmin);
    double bx = fabs(bmax);
    double ba = fabs(b);
    double dd = bx - bn;
    if (ba < bn)
      ret = hmin;
    else if (ba > bx)
      ret = hmax;
    else
      ret = 0.5 * ((hmin - hmax) * cos(ba * PI / dd - bn * PI / dd) + (hmin + hmax));
  }
  return(ret);
}

/* END bathyset()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads a coastline file and creates a hfun file                    */
/*-------------------------------------------------------------------*/
void hfun_from_coast(parameters_t *params, coamsh_t *cm, jigsaw_msh_t *J_hfun, int mode)
{
  FILE *fp = params->prmfd, *ef, *bf, *hf;
  char buf[MAXSTRLEN], key[MAXSTRLEN];
  int nbath;
  double *x, *y, *b, *r;
  double *hfun, xloc, yloc;
  int n, m, i, j, nhfun, intype, bmf = 0;
  int imeth = 0;
  int verbose = 0;
  int filef = 1;
  int sn, smooth = 0;
  int nce1, nce2;
  int nexp;
  double *exlat, *exlon, *exrad;
  double bmin, bmax, hmin, hmax, s;
  double deg2m = 60.0 * 1852.0;
  double expf = 0.0;
  jigsaw_REALS_array_t *hfvals = &J_hfun->_value;
  jigsaw_REALS_array_t *hfx    = &J_hfun->_xgrid;
  jigsaw_REALS_array_t *hfy    = &J_hfun->_ygrid;
  int perimf = 0; /* Set hfun perimeter to maximum resolution        */
  int npass = 0;
  int npoints;
  int npoly = 0;
  point *p;
  poly_t **pl;
  double *polyres;
  int *pmask;
  double *hmx, *hmn, *bmx, *bmn, *exf;
  msh_t *msh = cm->msh;
  int nord, orf = 0;
  double *exclat, *exclon, *excrad, *exccst;

  filef = (params->meshinfo) ? 1 : 0;
  /*-----------------------------------------------------------------*/
  /* Get the seed points if required                                 */
  if (mode & H_POINT) {
    prm_read_int(fp, "HFUN_POINTS", &npoints);
    p = malloc(npoints * sizeof(point));
    for(n = 0; n < npoints; n++) {
      if (fscanf(fp, "%lf %lf", &p[n].x, &p[n].y) != 2)
	hd_quit("hfun_from_coast: Format for HFUN_POINTS is 'lon lat'.\n");
    }
  }
  if (mode & H_POLY) {
    if (prm_read_char(fp, "HFUN_POLY", buf)) {
      char *fields[MAXSTRLEN * MAXNUMARGS];
      npoly = parseline(buf, fields, MAXNUMARGS);
      pl = (poly_t **)malloc(npoly * sizeof(poly_t));
      polyres = d_alloc_1d(npoly);
      for (n = 0; n < npoly; n++) {
	char *tok, *pname;
	strcpy(buf, fields[n]);
	pname = strtok(buf, ":");
	if ((tok = strtok(NULL, ":")) != NULL)
	  polyres[n] = atof(tok);
	else 
	  polyres[n] = 0.0;
	pl[n] = poly_create();
	if ((hf = fopen(pname, "r")) != NULL) {
	  m = poly_read(pl[n], hf);
	  fclose(hf);
	} else
	  hd_quit("hfun_from_coast: Can't open polygon file %s\n",pname);
      }
    }
    if (prm_read_int(fp, "NHFUN_POLY", &npoly)) {
      pl = (poly_t **)malloc(npoly * sizeof(poly_t));
      polyres = d_alloc_1d(npoly);
      for (n = 0; n < npoly; n++) {
	char *tok, *pname;
	sprintf(key, "HFUN_POLY%d", n);
	prm_read_char(fp, key, buf);
	pname = strtok(buf, ":");
	if ((tok = strtok(NULL, ":")) != NULL)
	  polyres[n] = atof(tok);
	else 
	  polyres[n] = 0.0;
	pl[n] = poly_create();
	if ((hf = fopen(pname, "r")) != NULL) {
	  m = poly_read(pl[n], hf);
	  fclose(hf);
	} else
	  hd_quit("hfun_from_coast: Can't open polygon file %s\n",pname);
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Define and read parameters                                      */
  /* Note: bmin and bmax are negative for wet cells.                 */
  /*-----------------------------------------------------------------*/
  /* Define and read parameters                                      */
  if (prm_read_int(fp, "NHFUN", &npass)) {
    bmax = hmax = -HUGE;
    bmin = hmin = HUGE;
    bmn = d_alloc_1d(npass);
    bmx = d_alloc_1d(npass);
    hmn = d_alloc_1d(npass);
    hmx = d_alloc_1d(npass);
    exf = d_alloc_1d(npass);
    for (n = 0; n < npass; n++) {
      sprintf(buf, "HMIN%d", n);
      prm_read_double(fp, buf, &hmn[n]);
      hmn[n] /= deg2m;
      sprintf(buf, "HMAX%d", n);
      prm_read_double(fp, buf, &hmx[n]);
      hmx[n] /= deg2m;
      sprintf(buf, "BMIN%d", n);
      prm_read_double(fp, buf, &bmn[n]);
      bmn[n] /= deg2m;
      sprintf(buf, "BMAX%d", n);
      prm_read_double(fp, buf, &bmx[n]);
      bmx[n] /= deg2m;
      sprintf(buf, "TYPE%d", n);
      prm_read_double(fp, buf, &exf[n]);
      bmin = min(bmin, bmn[n]);
      bmax = max(bmax, bmx[n]);
      hmax = max(hmax, hmx[n]);
      hmin = min(hmin, hmn[n]);
    }
    bmf = 1;
    cm->hfun_min = hmin * deg2m;
    cm->hfun_max = hmax * deg2m;
  } else {
    hmin = cm->hfun_min / deg2m;
    hmax = cm->hfun_max / deg2m;
    bmin = bmax = 0.0;
    if (prm_read_double(fp, "HFUN_BMIN", &bmin) && 
	prm_read_double(fp, "HFUN_BMAX", &bmax)) {
      bmf = 1;
      bmin /= deg2m;
      bmax /= deg2m;
    } else
      hd_quit("Can't find NHFUN or HFUN_BMIN & HFUN_BMAX.\n");
    prm_read_double(fp, "HFUN_TYPE", &expf);
  }
  prm_read_int(fp, "HFUN_SMOOTH", &smooth);

  nexp = 0;
  if (prm_read_int(fp, "HFUN_EXCLUDE", &nexp)) {
    exlon = d_alloc_1d(nexp);
    exlat = d_alloc_1d(nexp);
    exrad = d_alloc_1d(nexp);
    for (n = 0; n < nexp; n++) {
      if (fscanf(fp, "%lf %lf %lf", &exlon[n], &exlat[n], &exrad[n]) != 3)
	hd_quit("hfun_from_coast: Format for HFUN_EXCLUDE is 'lon lat radius'.\n");
    }
    for (i = 0; i < msh->npoint; i++) {
      if (msh->flag[i] & (S_OBC|S_LINK)) continue;
      for (n = 0; n < nexp; n++) {
	xloc = exlon[n] - msh->coords[0][i];
	yloc = exlat[n] - msh->coords[1][i];
	if (sqrt(xloc * xloc + yloc * yloc) * deg2m < exrad[n]) {
	  msh->flag[i] |= S_LINK;
	}
      }
    }
  }

  /* Set on a regular grid if required                               */
  imeth = I_MESH;
  if (prm_read_char(params->prmfd, "HFUN_GRID", buf)) {
    if (is_true(buf)) {
      if (imeth & I_MESH) imeth &= ~I_MESH;
      imeth |= I_GRID;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set up the grid to interpolate onto.                            */
  if (imeth & I_GRID) {
    /* Make a regular grid in the bounding box at the maximum        */
    /* resolution.                                                   */
    double res = 0.25 * (hmin + hmax);
    /*double res = 0.05;*/
    nce1 = (int)(fabs(cm->belon - cm->bslon) / res);
    nce2 = (int)(fabs(cm->belat - cm->bslat) / res);

    nhfun = nce1 * nce2;
    jigsaw_alloc_reals(hfx, nce1);
    jigsaw_alloc_reals(hfy, nce2);

    /* Fill in the x & y grids                                       */
    s = (cm->belon - cm->bslon) / (double)(nce1 - 1);
    for (i = 0; i < nce1; i++)
      hfx->_data[i] = s * (double)i + cm->bslon;
    s = (cm->belat - cm->bslat) / (double)(nce2 - 1);
    for (j = 0; j < nce2; j++)
      hfy->_data[j] = s * (double)j + cm->bslat;

  } else {
    /* Get a triangulation based on the mesh perimeter               */
    /*xy_to_d(d, cm->np, cm->x, cm->y);*/
    create_bounded_mesh(cm->np, cm->x, cm->y, hmin, hmax, J_hfun, filef);
    nhfun = J_hfun->_vert2._size;
  }
  jigsaw_alloc_reals(hfvals, nhfun);

  /*-----------------------------------------------------------------*/
  /* Read coastal over-ride values                                   */
  if (prm_read_int(fp, "HFUN_OVERRIDE", &nord)) {
    exclon = d_alloc_1d(nord);
    exclat = d_alloc_1d(nord);
    excrad = d_alloc_1d(nord);
    exccst = d_alloc_1d(nord);
    for (n = 0; n < nord; n++) {
      if (fscanf(fp, "%lf %lf %lf %lf", &exclon[n], &exclat[n], &excrad[n], &exccst[n]) != 4)
	hd_quit("hfun_from_coast: Format for HFUN_OVERRIDE is 'lon lat radius value'.\n");
      exccst[n] /= deg2m;
    }
    orf = 1;
  }

  /*-----------------------------------------------------------------*/
  /* Make a gridded distance to coast array                          */
  b = d_alloc_1d(nhfun);
  if (mode & H_POLY) pmask = i_alloc_1d(nhfun);
  if (imeth & I_GRID) {
    /* Get the minimum distance to the coast                         */
    n = 0;
    for (i = 0; i < nce1; i++) {
      xloc = hfx->_data[i];
      for (j = 0; j < nce2; j++) {
	yloc = hfy->_data[j];
	b[n] = HUGE;
	if (mode & H_CST)
	  b[n] = min(b[n], coast_dist(msh, xloc, yloc));
	if (mode & H_POINT)
	  b[n] = min(b[n], point_dist(npoints, p, xloc, yloc));
	if (mode & H_POLY)
	  b[n] = min(b[n], poly_dist(npoly, pl, bmin, xloc, yloc, &pmask[n]));
	if (verbose) printf("%d %f %f : %f\n",n, xloc, yloc, b[n]);
	n++;
      }
    }
    J_hfun->_flags = JIGSAW_EUCLIDEAN_GRID;
  } else {
    for (n = 0; n < nhfun; n++) {
      xloc = J_hfun->_vert2._data[n]._ppos[0];
      yloc = J_hfun->_vert2._data[n]._ppos[1];
      b[n] = HUGE;
      if (mode & H_CST) 
	b[n] = min(b[n], coast_dist(msh, xloc, yloc));
      if (mode & H_POINT)
	b[n] = min(b[n], point_dist(npoints, p, xloc, yloc));
      if (mode & H_POLY)
	b[n] = min(b[n], poly_dist(npoly, pl, bmin, xloc, yloc, &pmask[n]));
      if (verbose) printf("%d %f %f : %f\n",n, xloc, yloc, b[n]);
    }
    J_hfun->_flags = JIGSAW_EUCLIDEAN_MESH;
  }

  /*-----------------------------------------------------------------*/
  /* Over-ride distances values if required                          */
  if (orf) {

    if (imeth & I_GRID) {
      double *val = (imeth & I_USR) ? hfvals->_data : b;
      double smin[nord];
      int imin[nord], jmin[nord], nmin[nord];
      n = 0;
      for (m = 0; m < nord; m++) {
	smin[m] = HUGE;
	nmin[m] = -1;
      }
      for (i = 0; i < nce1; i++) {
	for (j = 0; j < nce2; j++) {
	  for (m = 0; m < nord; m++) {
	    if (exccst[m] > 0.0) continue;
	    xloc = exclon[m] - hfx->_data[i];
	    yloc = exclat[m] - hfy->_data[j];
	    s = sqrt(xloc * xloc + yloc * yloc) * deg2m;
	    if (s < excrad[m]) {
	      val[n] = (val[n] - fabs(exccst[m])) * s / excrad[m] + fabs(exccst[m]);
	      printf("%f\n",val[n]);
	      if (s < smin[m]) {
		smin[m] = s;
		imin[m] = i;
		jmin[m] = j;
		nmin[m] = n;
	      }
	    }
	  }
	  n++;
	}
      }
      for (m = 0; m < nord; m++) {
	if (nmin[m] >= 0 && exccst[m] <= 0.0) {
	  val[nmin[m]] = fabs(exccst[m]);
	  hd_warn("HFUN_OVERRIDE%d minimum dist @ %f %f\n", m,
		  hfx->_data[imin[m]],hfy->_data[jmin[m]]);
	}
      }
    } else {
      for (n = 0; n < nhfun; n++) {
	for (m = 0; m < nord; m++) {
	  if (exccst[m] > 0.0) continue;
	  xloc = exclon[m] - J_hfun->_vert2._data[n]._ppos[0];
	  yloc = exclat[m] - J_hfun->_vert2._data[n]._ppos[1];
	  s = sqrt(xloc * xloc + yloc * yloc) * deg2m;
	  if (s < excrad[m]) b[n] = (b[n] - fabs(exccst[m])) * s / excrad[m] + fabs(exccst[m]);
	}
      }
    }
  }

  /* Get the minimum and maximum distances                           */
  if (!bmf) {
    bmin = HUGE;
    bmax = 0.0;
    for (n = 0; n < nhfun; n++) {
      bmin = min(bmin, b[n]);
      bmax = max(bmax, b[n]);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Convert to the hfun function                                    */
  s = (bmin == bmax) ? 0.0 : (hmin - hmax) / (bmin - bmax);
  for (n = 0; n < nhfun; n++) {
    int dof = 0;
    b[n] = max(min(bmax, b[n]), bmin);

    if (npass) {
      hfvals->_data[n] = 0.0;
      for (m = 0; m < npass; m++) {
	if (b[n] >= bmn[m] && b[n] <= bmx[m])
	  hfvals->_data[n] = bathyset(b[n], bmn[m], bmx[m], 
				      hmn[m], hmx[m], exf[m]);
      }
      if (!hfvals->_data[n]) {
	if (b[n] >= bmax) hfvals->_data[n] = hmax;
	if (b[n] <= bmin) hfvals->_data[n] = hmin;
      }
    } else 
      hfvals->_data[n] = bathyset(b[n], bmin, bmax, hmin, hmax, expf);

    if (dof) {
    if (expf > 0.0) {

      /* Exponential                                               */
      hfvals->_data[n] = (hmin - hmax) * exp(b[n] / expf) + 
	(hmax - hmin * exp(-fabs(bmax) / expf));
    } else if (expf == 0) {
      /* Linear                                                    */
      hfvals->_data[n] = ((b[n] - bmax) * s + hmax);
    } else {
      double bn = fabs(bmin);
      double bx = fabs(bmax);
      double ba = fabs(b[n]);
      double dd = bx - bn;
      if (ba < bn)
	hfvals->_data[n] = hmin;
      else if (ba > bx)
	hfvals->_data[n] = hmax;
      else
	hfvals->_data[n] = 0.5 * ((hmin - hmax) * cos(ba * PI / dd - s) + (hmin + hmax));
    }
    if (verbose) printf("%d %f %f\n",n, hfvals->_data[n], b[n]);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the perimeter to maximum resolution                         */
  if (perimf && imeth & I_MESH) {
    for (n = 0; n < J_hfun->_edge2._size; n++) {
      i = J_hfun->_edge2._data[n]._node[0];
      hfvals->_data[i] = hmin;
      i = J_hfun->_edge2._data[n]._node[1];
      hfvals->_data[i] = hmin;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set resolution inside polygons if required                      */
  if (mode & H_POLY) {
    for (n = 0; n < npoly; n++) {
      if (polyres[n] == 0.0) 
	polyres[n] = hmin;
      else
	polyres[n] /= deg2m;
      cm->hfun_min = min(cm->hfun_min, polyres[n]);
      cm->hfun_max = max(cm->hfun_max, polyres[n]);
    }
    for (n = 0; n < nhfun; n++) {
      if ((m = pmask[n]) >= 0) {
	hfvals->_data[n] = polyres[m];
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Over-ride resolution values if required                         */
  if (orf) {
    double *val = hfvals->_data;
    if (imeth & I_GRID) {
      double smin[nord];
      int imin[nord], jmin[nord], nmin[nord];
      n = 0;
      for (m = 0; m < nord; m++) {
	smin[m] = HUGE;
	nmin[m] = -1;
      }
      for (i = 0; i < nce1; i++) {
	for (j = 0; j < nce2; j++) {
	  for (m = 0; m < nord; m++) {
	    if (exccst[m] <= 0.0) continue;
	    xloc = exclon[m] - hfx->_data[i];
	    yloc = exclat[m] - hfy->_data[j];
	    s = sqrt(xloc * xloc + yloc * yloc) * deg2m;
	    if (s < excrad[m]) {
	      val[n] = (val[n] - exccst[m]) * s / excrad[m] + exccst[m];
	      if (s < smin[m]) {
		smin[m] = s;
		imin[m] = i;
		jmin[m] = j;
		nmin[m] = n;
	      }
	    }
	  }
	  n++;
	}
      }
      for (m = 0; m < nord; m++) {
	if (nmin[m] >= 0 && exccst[m] > 0.0) {
	  val[nmin[m]] = exccst[m];
	  hmin = min(hmin, exccst[m]);
	  hmax = max(hmax, exccst[m]);
	  cm->hfun_min = hmin * deg2m;
	  cm->hfun_max = hmax * deg2m;
	  hd_warn("HFUN_OVERRIDER%d minimum dist @ %f %f\n", m,
		  hfx->_data[imin[m]],hfy->_data[jmin[m]]);
	}
      }
    } else {
      for (n = 0; n < nhfun; n++) {
	for (m = 0; m < nord; m++) {
	  if (exccst[m] <= 0.0) continue;
	  xloc = exclon[m] - J_hfun->_vert2._data[n]._ppos[0];
	  yloc = exclat[m] - J_hfun->_vert2._data[n]._ppos[1];
	  s = sqrt(xloc * xloc + yloc * yloc) * deg2m;
	  if (s < excrad[m]) val[n] = (val[n] - exccst[m]) * s / excrad[m] + exccst[m];
	}
      }
    }
    if (exclon) d_free_1d(exclon);
    if (exclat) d_free_1d(exclat);
    if (excrad) d_free_1d(excrad);
    if (exccst) d_free_1d(exccst);
  }

  /*-----------------------------------------------------------------*/
  /* Smooth if required                                              */
  if (smooth) {
    if (imeth & I_MESH) {
      int nv2 = J_hfun->_vert2._size;
      int *nv, **eiv;
      /*neighbour_finder_j(J_hfun, &nei);*/
      nv = i_alloc_1d(nv2);
      memset(nv, 0, nv2 * sizeof(int));
      for (n = 0; n < J_hfun->_tria3._size; n++) {
	for (j = 0; j < 3; j++) {
	  i = J_hfun->_tria3._data[n]._node[j];
	  nv[i]++;
	}
      }
      m = 0;
      for (n = 0; n < nv2; n++)
	if (nv[n] > m) m = nv[n];
      eiv = i_alloc_2d(nv2, m);
      memset(nv, 0, nv2 * sizeof(int));
      for (n = 0; n < J_hfun->_tria3._size; n++) {
	for (j = 0; j < 3; j++) {
	  i = J_hfun->_tria3._data[n]._node[j];
	  eiv[nv[i]][i] = n;
	  nv[i]++;
	}
      }
      r = d_alloc_1d(nhfun);
      for (sn = 0; sn < smooth; sn++) {
	/* Loop over all vertices                                      */
	for (m = 0; m < nv2; m++) {
	  r[m] = hfvals->_data[m];
	  s = 1.0;
	  /* Loop over all neighbouring triangles                       */
	  for (j = 0; j < nv[m]; j++) {
	    n = eiv[j][m];
	    for (i = 0; i < 3; i++) {
	      int nn = J_hfun->_tria3._data[n]._node[i];
	      if (nn != m && nn < nhfun) {
		r[m] += hfvals->_data[nn];
	      s += 1.0;
	      }
	    }
	  }
	  r[m] /= s;
	}
	memcpy(hfvals->_data, r, nhfun * sizeof(double));
      }
      d_free_1d(r);
      i_free_1d(nv);
      i_free_2d(eiv);
    }
    if (imeth & I_GRID) {
      int ii, jj, **ij2n;

      /* Get a mapping from ij to hfvals indices                     */
      ij2n = i_alloc_2d(nce1, nce2);
      n = 0;
      for (i = 0; i < nce1; i++) {
	for (j = 0; j < nce2; j++) {
	  ij2n[j][i] = n;
	  n++;
	}
      }

      /* Smooth                                                      */
      r = d_alloc_1d(nhfun);
      for (sn = 0; sn < smooth; sn++) {
	/* Loop over all vertices                                      */
	/*memcpy(r, hfvals->_data, nhfun * sizeof(double));*/
	for (i = 0; i < nce1; i++) {
	  int i1 = max(i-1, 0), i2 = min(i+1, nce1-1);
	  for (j = 0; j < nce2; j++) {
	    int j1 = max(j-1, 0), j2 = min(j+1, nce2-1);
	    s = 0.0;
	    n = ij2n[j][i];
	    r[n] = 0.0;
	    for (ii = i1; ii <= i2; ii++) {
	      for (jj = j1; jj <= j2; jj++) {
		ii = max(0, min(nce1-1, ii));
		jj = max(0, min(nce2-1, jj));
		m = ij2n[jj][ii];
		r[n] += hfvals->_data[m];
		s += 1.0;
	      }
	    }
	    if (s) r[n] /= s;
	  }
	}
	memcpy(hfvals->_data, r, nhfun * sizeof(double));
      }

      d_free_1d(r);
      i_free_2d(ij2n);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Write to file                                                   */
  if (filef) {
    mesh_ofile_name(params, key);
    sprintf(buf,"%s_h.txt", key);
    if ((ef = fopen(buf, "w")) != NULL) {
      fprintf(ef, "# Minimum distance to coast = %f m\n", bmin * deg2m);
      fprintf(ef, "# Maximum distance to coast = %f m\n", bmax * deg2m);
      if (imeth & I_MESH) { 
	for (n = 0; n < nhfun; n++) {
	  fprintf(ef, "%f %f %f %f\n",J_hfun->_vert2._data[n]._ppos[0], 
		  J_hfun->_vert2._data[n]._ppos[1], 
		  hfvals->_data[n], b[n]);
	}
      } else {
	n = 0;
	for (i = 0; i < nce1; i++) {
	  for (j = 0; j < nce2; j++) {
	    fprintf(ef, "%f %f %f %f\n",hfx->_data[i], hfy->_data[j],
		    hfvals->_data[n], b[n]);
	    n++;
	  }
	}
      }
      fclose(ef);
    }
    sprintf(buf,"%s_hfun.msh", key);
    if ((hf = fopen(buf, "w")) != NULL) {
      fprintf(hf, "# .msh geometry file\n");
      if (imeth & I_MESH) {
	/* EUCLIDEAN-MESH */
	fprintf(hf, "MSHID=2;EUCLIDEAN-MESH\n");
      } else {
	/* EUCLIDEAN-GRID */
	fprintf(hf, "MSHID=2;EUCLIDEAN-GRID\n");
      }
      fprintf(hf, "NDIMS=2\n");
      if (imeth & I_MESH) {
	fprintf(hf, "POINT=%d\n",nhfun);
	for (n = 0; n < nhfun; n++)
	  fprintf(hf, "%f;%f;0\n",J_hfun->_vert2._data[n]._ppos[0],
		  J_hfun->_vert2._data[n]._ppos[1]);
      } else {
	fprintf(hf, "coord=1;%d\n",nce1);
	for (i = 0; i < nce1; i++)
	  fprintf(hf, "%f\n",hfx->_data[i]);
	fprintf(hf, "coord=2;%d\n",nce2);
	for (j = 0; j < nce2; j++)
	  fprintf(hf, "%f\n",hfy->_data[j]);
      }
      if (imeth & I_MESH) {
	fprintf(hf, "value=%d;1\n",nhfun);
	for (n = 0; n < nhfun; n++)
	  fprintf(hf, "%f\n",hfvals->_data[n]);
      } else {
	fprintf(hf, "value=%d;1\n",nce1*nce2);
	n = 0;
	for (i = 0; i < nce1; i++)
	  for (j = 0; j < nce2; j++) {
	    fprintf(hf, "%f\n",hfvals->_data[n++]);
	  }
      }
      fclose(hf);
    }
  }
  /*for (n = 0; n < nhfun; n++) hfvals->_data[n] = 1.0;*/
  d_free_1d(b); 
  if(polyres) d_free_1d(polyres);
  if(pmask) i_free_1d(pmask);
  if (mode & H_POLY) {
    for (n = 0; n < npoly; n++)
      poly_destroy(pl[n]);
  free((poly_t **)pl);
  }
  if (DEBUG("init_m"))
    dlog("init_m", "Coastline weighting function computed OK\n");
}

/* END hfun_from_coast()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computs the minimum distance to a coastline                       */
/*-------------------------------------------------------------------*/
double coast_dist(msh_t *msh, double xloc, double yloc)
{
  int n;
  double d, x, y, dx, dy;
  double dist = HUGE;

  for (n = 0; n < msh->npoint; n++) {
    /* Coastline coordinates                                         */
    if (msh->flag[n] & (S_OBC|S_LINK)) continue;
    x = xloc - msh->coords[0][n];
    y = yloc - msh->coords[1][n];
    d = sqrt(x * x + y * y);
    dist = min(dist, d);
  }
  return(dist);
}

/* END coast_dist()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computs the minimum distance to a coastline                       */
/*-------------------------------------------------------------------*/
double point_dist(int npoints, point *p, double xloc, double yloc)
{
  int n;
  double d, x, y, dx, dy;
  double dist = HUGE;

  for (n = 0; n < npoints; n++) {
    x = xloc - p[n].x;
    y = yloc - p[n].y;
    d = sqrt(x * x + y * y);
    dist = min(dist, d);
  }
  return(dist);
}

/* END point_dist()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computs the minimum distance to a coastline                       */
/*-------------------------------------------------------------------*/
 double poly_dist(int npoly, poly_t **pl, double hmin, double xloc, double yloc, int *mask)
{
  int n, i;
  double d, x, y;
  double dist = HUGE;

  *mask = -1;
  for (n = 0; n < npoly; n++) {
    if (poly_contains_point(pl[n], xloc, yloc)) {
      *mask = n;
      return(hmin);
    } else {
      for (i = 0; i < pl[n]->n; i++) {
	x = xloc - pl[n]->x[i];
	y = yloc - pl[n]->y[i];
	d = sqrt(x * x + y * y);
	dist = min(dist, d);
      }
    }
  }
  return(dist);
}

/* END poly_dist()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Creates a constant triangulation given a bounding perimeter       */
/*-------------------------------------------------------------------*/
void create_bounded_mesh(int npts,     /* Number of perimeter points */
			 double *x,    /* Perimeter x coordinates    */
			 double *y,    /* Perimeter y coordinates    */
			 double rmin,  /* Min. resolution in m       */
			 double rmax,  /* Max. resolution in m       */
			 jigsaw_msh_t *J_mesh,
			 int filef
			 )
{
  FILE *ef;
  char buf[MAXSTRLEN], key[MAXSTRLEN];
  int np = npts - 1;

  /* local variables                                                 */
  int n, i, j;
  int jret = 0;
  
  /* Jigsaw variables                                                */
  jigsaw_jig_t J_jig;
  jigsaw_msh_t J_geom;

  /* Convenient pointers                                             */
  jigsaw_VERT2_array_t *jpoints = &J_geom._vert2;
  jigsaw_EDGE2_array_t *jedges  = &J_geom._edge2;

  /* Jigsaw initialisation routines                                  */
  jigsaw_init_jig_t (&J_jig);
  jigsaw_init_msh_t (&J_geom);

  /* Flag input mesh type                                            */
  J_geom._flags = JIGSAW_EUCLIDEAN_MESH;

  /* Allocate JIGSAW geom (input mesh) arrays                        */
  jigsaw_alloc_vert2(jpoints, npts);
  jigsaw_alloc_edge2(jedges,  npts);
  
  /* Set up the coordinates and points of the perimeter              */
  for (n = 0; n < np; n++) {
    /* Coordinates                                                   */
    jpoints->_data[n]._ppos[0] = x[n];
    jpoints->_data[n]._ppos[1] = y[n];
    jpoints->_data[n]._itag = 0;
    /* Connecting edges                                              */
    jedges->_data[n]._node[0] = n;
    if (n<np-1)
      jedges->_data[n]._node[1] = n+1;
    else
      jedges->_data[n]._node[1] = 0; // close polygon
    jedges->_data[n]._itag = 0;    
  }
  /*
  for (n = 0; n < np; n++) {
    printf("%f %f\n",jpoints->_data[jedges->_data[n]._node[0]]._ppos[0],
	   jpoints->_data[jedges->_data[n]._node[0]]._ppos[1]);
    printf("%f %f\n",jpoints->_data[jedges->_data[n]._node[1]]._ppos[0],
	   jpoints->_data[jedges->_data[n]._node[1]]._ppos[1]);
  }
  */
  /* Set up some flags/parameters                                    */
  J_jig._hfun_scal = JIGSAW_HFUN_ABSOLUTE;
  J_jig._hfun_hmin = 0.25 * (rmin + rmax);
  J_jig._hfun_hmax = 0.25 * (rmin + rmax);

  /*
  J_jig._hfun_scal = JIGSAW_HFUN_RELATIVE;
  J_jig._hfun_hmin = 0.0;
  J_jig._hfun_hmax = 0.001;
  */
  /* Call Jigsaw                                                     */
  J_jig._verbosity = 0;
  jret = jigsaw_make_mesh(&J_jig, &J_geom, NULL, NULL, J_mesh);
  if (jret) hd_quit("create_bounded_mesh: Error calling JIGSAW\n");

  /* Write to file                                                   */
  if (filef) {
    mesh_ofile_name(params, key);
    sprintf(buf,"%s_e.txt", key);
    if ((ef = fopen(buf, "w")) != NULL) {
      /*
      for (n = 0; n < J_mesh->_vert2._size; n++) {
	fprintf(ef, "%f %f\n", J_mesh->_vert2._data[n]._ppos[0],
		J_mesh->_vert2._data[n]._ppos[1]);
      }
      */
      for (n = 0; n < J_mesh->_tria3._size; n++) {
	fprintf(ef, "%f %f\n", 
		J_mesh->_vert2._data[J_mesh->_tria3._data[n]._node[0]]._ppos[0],
		J_mesh->_vert2._data[J_mesh->_tria3._data[n]._node[0]]._ppos[1]);
	fprintf(ef, "%f %f\n", 
		J_mesh->_vert2._data[J_mesh->_tria3._data[n]._node[1]]._ppos[0],
		J_mesh->_vert2._data[J_mesh->_tria3._data[n]._node[1]]._ppos[1]);
	fprintf(ef, "%f %f\n", 
		J_mesh->_vert2._data[J_mesh->_tria3._data[n]._node[2]]._ppos[0],
		J_mesh->_vert2._data[J_mesh->_tria3._data[n]._node[2]]._ppos[1]);
	fprintf(ef, "NaN NaN\n");
      }
      fclose(ef);
    }
  }
  /* Cleanup                                                         */
  jigsaw_free_msh_t(&J_geom);
}

/* END create_bounded_mesh()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void init_J_jig(jigsaw_jig_t *J_jig)
{
  J_jig->_verbosity = 0;
  J_jig->_geom_seed = 8;
  J_jig->_geom_feat = 0;
  J_jig->_geom_eta1 = 45.0;
  J_jig->_geom_eta2 = 45.0;
  J_jig->_hfun_scal = JIGSAW_HFUN_RELATIVE;
  J_jig->_hfun_hmax = 0.02;
  J_jig->_hfun_hmin = 0.0;
  J_jig->_mesh_dims = 3;
  /* JIGSAW_KERN_DELAUNAY or JIGSAW_KERN_DELFRONT */
  J_jig->_mesh_kern = JIGSAW_KERN_DELFRONT;
  /*J_jig->_mesh_iter = Inf;*/
  J_jig->_mesh_top1 = 0;
  J_jig->_mesh_top2 = 0;
  J_jig->_mesh_rad2 = 1.05;
  J_jig->_mesh_rad3 = 2.05;
  J_jig->_mesh_off2 = 0.9;
  J_jig->_mesh_off3 = 1.1;
  J_jig->_mesh_snk2 = 0.2;
  J_jig->_mesh_snk3 = 0.33;
  J_jig->_mesh_eps1 = 0.33;
  J_jig->_mesh_eps2 = 0.33;
  J_jig->_mesh_vol3 = 0.0;
  J_jig->_optm_iter = 32;
  J_jig->_optm_qtol = 1.E-04;
  J_jig->_optm_qlim = 0.9250;

}

/* END init_J_jig()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Creates a constant triangulation given a bounding perimeter       */
/*-------------------------------------------------------------------*/
void create_jigsaw_mesh(coamsh_t *cm, 
			jigsaw_msh_t *J_mesh,
			jigsaw_msh_t *J_hfun,
			int powf, int stproj,
			int filef
			)
{
  FILE *ef;
  char buf[MAXSTRLEN], key[MAXSTRLEN];
  double xmid = 0.0, ymid = 0.0;
  
  /* Constants                                                       */
  double deg2m = 60.0 * 1852.0;
  msh_t *msh = cm->msh;
  jig_t *jig = cm->jig;

  /* local variables                                                 */
  int n, i, j;
  int jret = 0;
  
  /* Jigsaw variables                                                */
  jigsaw_jig_t J_jig;
  jigsaw_msh_t J_geom;

  /* Convenient pointers                                             */
  jigsaw_VERT2_array_t *jpoints = &J_geom._vert2;
  jigsaw_EDGE2_array_t *jedges  = &J_geom._edge2;

  /* Jigsaw initialisation routines                                  */
  jigsaw_init_jig_t (&J_jig);
  jigsaw_init_msh_t (&J_geom);
  init_J_jig(&J_jig);
  /* Power mesh                                                      */
  if (powf) J_jig._optm_dual = 1;

  /* Flag input mesh type                                            */
  J_geom._flags = JIGSAW_EUCLIDEAN_MESH;

  /* Allocate JIGSAW geom (input mesh) arrays                        */
  jigsaw_alloc_vert2(jpoints, msh->npoint);
  jigsaw_alloc_edge2(jedges,  msh->npoint);
  
  /* Set up the coordinates and points of the perimeter              */
  for (n = 0; n < msh->npoint; n++) {
    /* Coordinates                                                   */
    jpoints->_data[n]._ppos[0] = msh->coords[0][n];
    jpoints->_data[n]._ppos[1] = msh->coords[1][n];
    jpoints->_data[n]._itag = 0;
  }
  for (n = 0; n < msh->edge2; n++) {
    /* Connecting edges                                              */
    jedges->_data[n]._node[0] = msh->edges[0][n];
    jedges->_data[n]._node[1] = msh->edges[1][n];
    jedges->_data[n]._itag = 0;
  }

  /* Set up some flags/parameters                                    */
  J_jig._hfun_scal = JIGSAW_HFUN_ABSOLUTE;
  if (stproj) {
    /* Leave in metres */
    J_jig._hfun_hmin = cm->hfun_min;
    J_jig._hfun_hmax = cm->hfun_max;
  } else {
    J_jig._hfun_hmin = cm->hfun_min / deg2m;
    J_jig._hfun_hmax = cm->hfun_max / deg2m;
  }

  J_jig._verbosity = 1;
  
  /*J_jig._mesh_eps1 = 10.0;*/
  /*
  J_jig._hfun_hmin = 0.01;
  J_jig._hfun_hmax = 0.02;
  */

  /* Stereographic projection */
  if (stproj) {
    /* Calculate mid-points */
    for (n = 0; n < msh->npoint; n++) {
      /* Figure out the centroid */
      xmid += jpoints->_data[n]._ppos[0];
      ymid += jpoints->_data[n]._ppos[1];
    }
    xmid /= msh->npoint; xmid = DEG2RAD(xmid);
    ymid /= msh->npoint; ymid = DEG2RAD(ymid);

    /* Perform forward stereographic projection of input data          */
    st_transform(&J_geom, ST3PROJ_FWD, xmid, ymid);

    /* Convert hfun grid to mesh, if needed */
    if (J_hfun->_flags == JIGSAW_EUCLIDEAN_GRID) {
      jigsaw_msh_t J_hfun_new;
      jigsaw_init_msh_t (&J_hfun_new);
      convert_jigsaw_grid_to_mesh(J_hfun, &J_hfun_new);
      /* Swap out */
      J_hfun = &J_hfun_new;
    }
    /* Projection of hfun */
    st_transform(J_hfun, ST3PROJ_FWD, xmid, ymid);
  }

  /* Call Jigsaw                                                     */
  jret = jigsaw_make_mesh(&J_jig, &J_geom, NULL, J_hfun, J_mesh);
  if (jret) hd_quit("create_jigsaw_mesh: Error calling JIGSAW\n");

  /* Perform inverse stereographic projection of output data         */
  if (stproj) st_transform(J_mesh, ST3PROJ_INV, xmid, ymid);

  /* Write to file                                                   */
  if (filef) {
    mesh_ofile_name(params, key);
    sprintf(buf,"%s_t.txt", key);
    if ((ef = fopen(buf, "w")) != NULL) {
      for (n = 0; n < J_mesh->_tria3._size; n++) {
	fprintf(ef, "%f %f\n", 
		J_mesh->_vert2._data[J_mesh->_tria3._data[n]._node[0]]._ppos[0],
		J_mesh->_vert2._data[J_mesh->_tria3._data[n]._node[0]]._ppos[1]);
	fprintf(ef, "%f %f\n", 
		J_mesh->_vert2._data[J_mesh->_tria3._data[n]._node[1]]._ppos[0],
		J_mesh->_vert2._data[J_mesh->_tria3._data[n]._node[1]]._ppos[1]);
	fprintf(ef, "%f %f\n", 
		J_mesh->_vert2._data[J_mesh->_tria3._data[n]._node[2]]._ppos[0],
		J_mesh->_vert2._data[J_mesh->_tria3._data[n]._node[2]]._ppos[1]);
	fprintf(ef, "NaN NaN\n");
      }
      fclose(ef);
    }
  }

  /* Cleanup                                                         */
  jigsaw_free_msh_t(&J_geom);
  jigsaw_free_msh_t(J_hfun);
  if (DEBUG("init_m"))
    dlog("init_m", "JIGSAW mesh created OK\n");
}

/* END create_jigsaw_mesh()                                          */
/*-------------------------------------------------------------------*/
#endif

/*-------------------------------------------------------------------*/
/* Reads a JIGSAW .msh triangulation file and converts to a Voronoi  */
/* mesh.                                                             */
/*-------------------------------------------------------------------*/
void convert_jigsaw_msh(parameters_t *params, char *infile,
			void *jmsh)
{
  FILE *fp, *ef, *vf, *tf, *cf, *jf, *cef;
  int meshid;
  int ndims;
  int npoints;
  int n, nn, i, j, m, mm, id;
  int e1, e2, t1, t2;
  char buf[MAXSTRLEN], code[MAXSTRLEN], key[MAXSTRLEN];
  char mshtype[MAXSTRLEN];
  double x1, y1, x2, y2, d1, d2;
  double nvoints, nvedges;
  double radius;
  point *voints;
  delaunay* d;
  triangle* t, tn;
  int ns, *sin, *edges, *mask, **nei, **e2t;
  int intype;
  int has_proj = (strlen(params->projection) > 0);
  int is_geog = has_proj && (strcasecmp(params->projection, GEOGRAPHIC_TAG) == 0);
  int verbose = 0;
  int filef = 1;
  int filefj = 0;         /* Print JIGSAW .msh output file           */
  int centref = 1;        /* 0 : centre of mass for triangle centres */
                          /* 1 : circumcentres for triangle centres  */
                          /* 2 : orthocentres for triangle centres   */
  int newcode = 1;        /* Optimised neighbour finding             */
  int obtusef = 1;        /* Use centroid for obtuse triangles       */
  int dtri = -1;          /* Debugging to find triangles             */
  int dedge = -1;         /* Debugging for edges                     */

#ifdef HAVE_JIGSAWLIB
  jigsaw_msh_t *msh = (jigsaw_msh_t*)jmsh;
#endif
  filef = (params->meshinfo) ? 1 : 0;
  params->d = d = malloc(sizeof(delaunay));
  params->us_type |= US_JUS;
  if (params->us_type & US_POW) centref = 2;
  if (!centref && obtusef) obtusef = 0;
  if (centref == 2 && obtusef) obtusef = 0;

  is_geog = 0;
  if (params->us_type & US_STER) is_geog = 1;

  /* centref = 1 with is_geog = 1 appears to generate some odd       */
  /* looking cells.                                                  */
  if (centref == 1 && is_geog) is_geog = 0;

  if (jmsh == NULL) {
    if ((fp = fopen(infile, "r")) == NULL)
      hd_quit("convert_msh: Can't open JIGSAW file %s\n", infile);

    prm_skip_to_end_of_tok(fp, "mshid", "=", buf);
    meshid = atoi(buf);
    sprintf(key, "MSHID=%d", meshid);
    prm_skip_to_end_of_tok(fp, key, ";", mshtype);
    prm_skip_to_end_of_tok(fp, "ndims", "=", buf);
    ndims = atoi(buf);
    prm_skip_to_end_of_tok(fp, "point", "=", buf);
    d->npoints = atoi(buf);
    intype = 1;
    /*
    if (strcmp(mshtype, "ELLIPSOID-MESH") == 0) {
      if (prm_skip_to_end_of_tok(fp, "RADII", "=", buf))
	radius = atof(buf);
      else
	hd_quit("Jigsaw file must specify RADII for ELLIPSOID-MESH files.\n");
      if (!is_geog)
	hd_quit("Must use PROJECTION GEOGRAPHIC with Jigsaw ELLIPSOID-MESH files.\n");
    }
    */
  } else {
#ifdef HAVE_JIGSAWLIB
    /* Always the case for now */
    ndims = 2;
    d->npoints = msh->_vert2._size;
    intype = 0;
#endif
  }

  /*
  prm_read_int(fp, "mshid=", &meshid);
  prm_read_int(fp, "ndims=", &ndims);
  prm_read_int(fp, "point=", &d->npoints);
  */

  /*-----------------------------------------------------------------*/
  /* Write to file                                                   */
  if (filef) {
    mesh_ofile_name(params, key);
    sprintf(buf,"%s_jv.txt", key);
    if ((vf = fopen(buf, "w")) == NULL) filef = 0;
    sprintf(buf,"%s_je.txt", key);
    if ((ef = fopen(buf, "w")) == NULL) filef = 0;
    sprintf(buf,"%s_jc.txt", key);
    if ((cf = fopen(buf, "w")) == NULL) filef = 0;
    sprintf(buf,"%s_jce.txt", key);
    if ((cef = fopen(buf, "w")) == NULL) filef = 0;
    sprintf(buf,"%s_jt.txt", key);
    if ((tf = fopen(buf, "w")) == NULL) filef = 0;
    sprintf(buf,"%s_jig.msh", key);
    if ((jf = fopen(buf, "w")) == NULL) 
      filefj = 0;
    else {
      fprintf(jf, "# .msh geometry file\n");
      fprintf(jf, "mshid=1\n");
      fprintf(jf, "ndims=2\n");
      fprintf(jf, "point=%d\n",d->npoints);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Read the triangulation specification.                           */
  /* Points in the triangulation.                                    */
  d->points = malloc(d->npoints * sizeof(point));
  for (i = 0; i < d->npoints; i++) {
    point* p = &d->points[i];
    if (intype) {
      if (fscanf(fp, "%lf;%lf;%d", &p->x, &p->y, &id) != 3)
	hd_warn("convert_jigsaw_msh: Incorrect number of point entries in %s, line %d\n", infile, i);
    } else {
#ifdef HAVE_JIGSAWLIB
	p->x = msh->_vert2._data[i]._ppos[0];
	p->y = msh->_vert2._data[i]._ppos[1];
	p->z = 0.0;
	if (centref == 2)
	  p->z = msh->_power._data[i];
#endif
    }
    if (filef) fprintf(cf, "%f %f\n", p->x, p->y);
    if (filefj) fprintf(jf, "%f;%f;0\n", p->x, p->y);
  }
  if (intype && centref == 2) {
    prm_skip_to_end_of_tok(fp, "power", "=", buf);
    n = atoi(buf);
    for (i = 0; i < n; i++) {
      point* p = &d->points[i];
      if (fscanf(fp, "%lf", &p->z) != 1)
	hd_warn("convert_jigsaw_msh: Incorrect number of power entries in %s, line %d\n", infile, i);
    }
  }

  if (filef) fclose(cf);

  /*-----------------------------------------------------------------*/
  /* Boundary edges in the triangulation.                            */
  /*prm_read_int(fp, "edge2=", &ns);*/
  ns = 0;
  if (intype) {
    if (prm_skip_to_end_of_tok(fp, "edge2", "=", buf)) {
      ns = atoi(buf);
    }
  } else {
#ifdef HAVE_JIGSAWLIB
    ns = msh->_edge2._size;
    /*
    if (ns) {
      sin = i_alloc_1d(2 * ns);
      j = 0;
      for (i = 0; i < ns; i++) {
	if (intype) {
	  if (fscanf(fp, "%d;%d;%d", &e1, &e2, &id) != 3)
	    hd_warn("convert_msh: Incorrect number of edge entries in %s, line %d\n", infile, i);
	} else {
	  e1 = msh->_edge2._data[i]._node[0];
	  e2 = msh->_edge2._data[i]._node[1];
	}
	sin[j++] = e1;
	sin[j++] = e2;
	if (filef) {
	  fprintf(ef, "%f %f\n", d->points[e1].x, d->points[e1].y);
	  fprintf(ef, "%f %f\n", d->points[e2].x, d->points[e2].y);
	  fprintf(ef, "NaN NaN\n");
	}
      }
    }
    */
#endif
  }
  if (filefj) fprintf(jf, "edge2=%d\n",ns);
  if (ns) {
    sin = i_alloc_1d(2 * ns);
    j = 0;
    for (i = 0; i < ns; i++) {
      if (intype) {
	if (fscanf(fp, "%d;%d;%d", &e1, &e2, &id) != 3)
	  hd_warn("convert_msh: Incorrect number of edge entries in %s, line %d\n", infile, i);
      } else {
#ifdef HAVE_JIGSAWLIB
	e1 = msh->_edge2._data[i]._node[0];
	e2 = msh->_edge2._data[i]._node[1];
#endif
      }
      sin[j++] = e1;
      sin[j++] = e2;
      if (filef) {
	fprintf(ef, "%f %f\n", d->points[e1].x, d->points[e1].y);
	fprintf(ef, "%f %f\n", d->points[e2].x, d->points[e2].y);
	fprintf(ef, "NaN NaN\n");
      }
      if (filefj) fprintf(jf, "%d;%d\n",e1, e2);
    }
  }
  if (filef) fclose(ef);

  /*-----------------------------------------------------------------*/
  /* Triangles in the triangulation.                                 */
  if (intype) {
    prm_skip_to_end_of_tok(fp, "tria3", "=", buf);
    d->ntriangles = atoi(buf);
  } else {
#ifdef HAVE_JIGSAWLIB
    d->ntriangles = msh->_tria3._size;
#endif
  }
  d->triangles = malloc(d->ntriangles * sizeof(triangle));
  if (filefj) fprintf(jf, "tria3=%d\n", d->ntriangles);
  for (i = 0; i < d->ntriangles; i++) {
    t = &d->triangles[i];
    if (intype) {
      if (fscanf(fp, "%d;%d;%d;%d", &t->vids[0], &t->vids[1], &t->vids[2], &id) != 4)
	hd_warn("convert_msh: Incorrect number of triangle entries in %s, line %d\n", infile, i);
    } else {
#ifdef HAVE_JIGSAWLIB
      t->vids[0] = msh->_tria3._data[i]._node[0];
      t->vids[1] = msh->_tria3._data[i]._node[1];
      t->vids[2] = msh->_tria3._data[i]._node[2];
#endif
    }
    if (filefj) fprintf(jf, "%d;%d;%d;0\n",t->vids[0], t->vids[1], t->vids[2]);
    /* Print the centre of mass if required                          */
    if (filef) {
      double po[2];
      double p1[3][2];
      for (n = 0; n < 3; n++) {
	j = t->vids[n];
	p1[n][0] = d->points[j].x;
	p1[n][1] = d->points[j].y;
      }
      tria_com(po, p1[0], p1[1], p1[2]);
      fprintf(cef, "%f %f\n", po[0], po[1]);

      /* Find a triange with location limits                         */
      if (dtri >= 0) {
	if(po[0] > 130.766 && po[0] < 130.768 && 
	   po[1] > -8.3175 && po[1] < -8.316) {
	  printf("TRI %d\n",i);
	  dtri = i;
	  for (n = 0; n < 3; n++) {
	    j = t->vids[n];
	    printf("Vertex%d = [%f %f]\n",n, d->points[j].x, d->points[j].y);
	  }
	}
      }
    }
  }
  if (filefj) fclose(jf);
  if (filef) fclose(cef);

  /*-----------------------------------------------------------------*/
  /* Extract the edges from the triangles and populate d->edges.      */
  /* Count the number of edges.                                      */
  d->nedges = 0;
  edges = i_alloc_1d(7 * d->ntriangles);
  e2t = i_alloc_2d(2, 7 * d->ntriangles);
  for (j = 0; j < d->ntriangles; j++) {
    t = &d->triangles[j];
    for (n = 0; n < 3; n++) {
      int found = 0;
      nn = (n == 2) ? 0 : n + 1;
      e1 = t->vids[n];
      e2 = t->vids[nn];
      for (i = 0; i < d->nedges; i++) {
	if ((edges[2*i] == e1 && edges[2*i+1] == e2) ||
	    (edges[2*i] == e2 && edges[2*i+1] == e1)) {
	  found = 1;
	  break;
	}
      }
      if (e1 < 0 || e1 >= d->npoints) hd_quit("convert_jigsaw_msh: found invalid coordinate index in triangle %d (%d)\n", j, e1);
      if (e2 < 0 || e2 >= d->npoints) hd_quit("convert_jigsaw_msh: found invalid coordinate index in triangle %d (%d)\n", j, e2);
      if (!found) {
	edges[2*d->nedges] = e1;
	edges[2*d->nedges+1] = e2;
	e2t[d->nedges][0] = j;
	e2t[d->nedges][1] = n;
	d->nedges++;
      }
    }
  }
  d->edges = i_alloc_1d(2*d->nedges);
  memcpy(d->edges, edges, 2 * d->nedges * sizeof(int));
  i_free_1d(edges);
  if (verbose) printf("npoints=%d nedges=%d ntriangles=%d\n", d->npoints, d->nedges, d->ntriangles);
  if (filef) {
    for (i = 0; i < d->nedges; i++) {
      fprintf(tf, "%f %f\n",d->points[d->edges[2*i]].x,d->points[d->edges[2*i]].y);
      fprintf(tf, "%f %f\n",d->points[d->edges[2*i+1]].x,d->points[d->edges[2*i+1]].y);
      fprintf(tf, "NaN NaN\n");
    }
  fclose(tf);
  }

  /*-----------------------------------------------------------------*/
  /* Find the triangle neighbours                                    */
  if (newcode) 
    neighbour_finder_b(d, &nei);
  /*
  else {
    nei = i_alloc_2d(3, d->ntriangles);
    for (j = 0; j < d->ntriangles; j++)
      for (n = 0; n < 2; n++)
	nei[j][n] = -1;
    for (j = 0; j < d->ntriangles; j++) {
      t = &d->triangles[j];
      for (n = 0; n < 3; n++) {
	nn = (n == 2) ? 0 : n + 1;
	t1 = t->vids[n];
	t2 = t->vids[nn];
	for (i = 0; i < d->ntriangles, i != j; i++) {
	  tn = &d->triangles[i];
	  for (m = 0; m < 3; m++) {
	    mm = (m == 2) ? 0 : m + 1;
	    t1 = tn->vids[m];
	    t2 = tn->vids[mm];
	  }
	}
      }
    }
  }
  */

  /*-----------------------------------------------------------------*/
  /* Find the triangles that share each edge and save the centroids. */
  /* Note: triangle centroids are the vertices of Voronoi cells      */
  /* (i.e. the endpoints of Voronoi edges; vdges[]) and are common   */
  /* to multiple Voronoi edges (3 for hexagons). mask[] ensures that */
  /* these vertices are only accounted for once.                     */
  d->nvoints = 0;
  d->nvdges = d->nedges;
  d->voints = malloc(3 * d->ntriangles * sizeof(point));
  d->vdges = malloc(d->nvdges * 2 * sizeof(int));
  mask = i_alloc_1d(d->ntriangles);
  for (j = 0; j < d->ntriangles; j++) 
    mask[j] = -1;

  for (i = 0; i < d->nedges; i++) {
    int tri[2];
    double po[2];

    /* Coordinate indices (points[e1], points[e2]) for the edge      */
    e1 = d->edges[2*i];
    e2 = d->edges[2*i+1];
    m = 0;

    /* Find the two triangles that share the edge. For boundaries     */
    /* there is only one trianlge.                                   */                            
    if (newcode) {
      m = 1;
      tri[0] = e2t[i][0];
      j = e2t[i][1];
      if ((tri[1] = nei[j][tri[0]]) >= 0) m++;
    } else {
    for (j = 0; j < d->ntriangles; j++) {
      t = &d->triangles[j];
      for (n = 0; n < 3; n++) {
	nn = (n == 2) ? 0 : n + 1;
	t1 = t->vids[n];
	t2 = t->vids[nn];
	if ((e1 == t1 && e2 == t2 ) ||
	    (e2 == t1 && e1 == t2 )) {
	  tri[m++] = j;
	  break;
	}
      }
      if (m == 2) break;
    }
    }

    /* Debugging for edges                                           */
    if (i == dedge)
      printf("edge %d shares %d triangles (%d & %d) %d %d\n",
	     i, m, tri[0], tri[1], e2t[i][0], e2t[i][1]);

    /* Debugging for triangles                                        */
    if (dtri >= 0) {
    if (tri[0] == dtri|| tri[1] == dtri)
      printf("tri0=%d, tri1=%d, edge=%d[%f %f]-[%f %f], m=%d\n",
	     tri[0], tri[1], i, d->points[e1].x, d->points[e1].y,
	     d->points[e2].x, d->points[e2].y, m);
    }

    if (verbose) {
      printf("%f %f\n",d->points[d->edges[2*i]].x,d->points[d->edges[2*i]].y);
      printf("%f %f\n",d->points[d->edges[2*i+1]].x,d->points[d->edges[2*i+1]].y);
      printf("NaN NaN\n");
    }

    /* Get the centres of the triangles. These define the vertices of */
    /* Voronoi cells.                                                */
    if (m == 1) {
      /* If only one triangle is a neighbour to the edge, then use   */
      /* the centroid of that triangle and the midpoint of the edge  */
      /* as the Voronoi edge.                                        */
      double p1[3][3], v1[3], v2[3];
      int o1;
      x1 = y1 = x2 = y2 = d1 = 0.0;
      for (n = 0; n < 3; n++) {
	t = &d->triangles[tri[0]];
	j = t->vids[n];
	p1[n][0] = d->points[j].x;
	p1[n][1] = d->points[j].y;
	p1[n][2] = d->points[j].z;

	d1 += 1.0;

	/* Set the vertex coordinates                                */
	v1[0] = d->points[e1].x;
	v1[1] = d->points[e1].y;
	v1[2] = d->points[e1].z;
	v2[0] = d->points[e2].x;
	v2[1] = d->points[e2].y;
	v2[2] = d->points[e2].z;

	/* Compute the mid-point                                     */
	edge_cen(centref, is_geog, v1, v2, po);

	/* Set the second triangle 'centre' to the mid-point         */
	x2 = po[0];
	y2 = po[1];
      }
      o1 = (centref) ? 0 : 1;
      if (obtusef) {
	if ((d2 = is_obtuse(p1[0],p1[1],p1[2]))) {
	  if (DEBUG("init_m"))
	    dlog("init_m", "convert_jigsaw_msh: Obtuse triangle at edge %d, tri %d (%5.2f)\n", i, tri[0], d2);
	  hd_warn("convert_jigsaw_msh: Obtuse triangle at edge %d, tri %d (%5.2f)\n", i, tri[0], d2);
	  o1 = 1;
	}
      }
      if (o1) {
	tria_com(po, p1[0], p1[1], p1[2]);
      	x1 = po[0]; y1 = po[1];
      } else {
	/* Compute the centre of the first triangle                  */
	tri_cen(centref, is_geog, p1[0], p1[1], p1[2], po);
	x1 = po[0]; y1 = po[1];
      }
      if(i == dedge) 
	printf("TRI%d centres [%f %f] and [%f %f] : %f\n",tri[0],x1,y1,x2,y2,d1);

    } else if (m == 2) {
      /* If two triangles are neighbours to the edge, then use the   */
      /* centroids of each triangle as the Voronoi edge.             */
      double p1[3][3], p2[3][3];
      int o1, o2;
      x1 = y1 = x2 = y2 = d1 = 0.0;
      for (n = 0; n < 3; n++) {
	t = &d->triangles[tri[0]];
	j = t->vids[n];
	p1[n][0] = d->points[j].x;
	p1[n][1] = d->points[j].y;
	p1[n][2] = d->points[j].z;

	t = &d->triangles[tri[1]];
	j = t->vids[n];
	p2[n][0] = d->points[j].x;
	p2[n][1] = d->points[j].y;
	p2[n][2] = d->points[j].z;

	d1 += 1.0;
      }
      o1 = o2 = (centref) ? 0 : 1;
      if (obtusef) {
	if ((d2 = is_obtuse(p1[0],p1[1],p1[2]))) {
	  if (DEBUG("init_m"))
	    dlog("init_m", "convert_jigsaw_msh: Obtuse triangle at edge %d, tri %d (%5.2f@[%f %f])\n", 
		 i, tri[0], d2, p1[0][0], p1[0][1]);
	  /*
	  hd_warn("convert_jigsaw_msh: Obtuse triangle at edge %d, tri %d (%5.2f@[%f %f])\n", 
		  i, tri[0], d2, p1[0][0], p1[0][1]);
	  */
	  o1 = 1;
	}
	if ((d2 = is_obtuse(p2[0],p2[1],p2[2]))) {
	  if (DEBUG("init_m"))
	    dlog("init_m", "convert_jigsaw_msh: Obtuse triangle at edge %d, tri %d (%5.2f@[%f %f])\n", 
		 i, tri[1], d2, p2[0][0], p2[0][1]);
	  /*
	  hd_warn("convert_jigsaw_msh: Obtuse triangle at edge %d, tri %d (%5.2f@[%f %f])\n", 
		  i, tri[1], d2, p2[0][0], p2[0][1]);
	  */
	  o2 = 1;
	}			      
      }
      if (o1) {
	tria_com(po, p1[0], p1[1], p1[2]);
      	x1 = po[0]; y1 = po[1];
      } else {
	tri_cen(centref, is_geog, p1[0], p1[1], p1[2], po);
	x1 = po[0]; y1 = po[1];
      }
      if (o2) {
	tria_com(po, p2[0], p2[1], p2[2]);
      	x2 = po[0]; y2 = po[1];
      } else {
	tri_cen(centref, is_geog, p2[0], p2[1], p2[2], po);
	x2 = po[0]; y2 = po[1];
      }
      if (obtusef && o1) {
	if (DEBUG("init_m"))
	  dlog("init_m", "  tri%d [%f %f], [%f %f], [%f %f]: centre [%f %f]\n",
	       tri[0], p1[0][0], p1[0][1],
	       p1[1][0], p1[1][1], p1[2][0], p1[2][1], x1, y1);
	/*
	hd_warn("  tri%d [%f %f], [%f %f], [%f %f]: centre [%f %f]\n",
		tri[0], p1[0][0], p1[0][1],
		p1[1][0], p1[1][1], p1[2][0], p1[2][1], x1, y1);
	*/
      }
      if (obtusef && o2) {
	if (DEBUG("init_m"))
	  dlog("init_m", "  tri%d [%f %f], [%f %f], [%f %f]: centre [%f %f]\n",
		tri[1], p2[0][0], p2[0][1],
		p2[1][0], p2[1][1], p2[2][0], p2[2][1], x2, y2);
	/*
	hd_warn("  tri%d [%f %f], [%f %f], [%f %f]: centre [%f %f]\n",
		tri[1], p2[0][0], p2[0][1],
		p2[1][0], p2[1][1], p2[2][0], p2[2][1], x2, y2);
	*/
      }
    }

    /* Set a new Voronoi point coordinate (vpoints[]) if it doesn't  */
    /* already exist and set the Voronoi edge endpoints (vdges[]) to */
    /* the corresponding coordinate indices. The triangle edges and  */
    /* Voronoi edges having the same index should be perpendicular.  */
    if (d1) {
      point *p;
      /* First Voronoi edge point. Should always be the centre of a  */
      /* triangle.                                                   */
      if (mask[tri[0]] < 0) {
	p = &d->voints[d->nvoints];
	p->x = x1;
	p->y = y1;
	mask[tri[0]] = d->nvoints;
	d->nvoints++;
      }
      d->vdges[2*i] = mask[tri[0]];
      if( i == dedge) 
	printf("first Voronoi edge = [%f %f]\n", 
	       d->voints[d->vdges[2*i]].x, d->voints[d->vdges[2*i]].y);

      /* Second Voronoi edge point. May be a triangle centroid or    */
      /* the midpoint of a triangle edge.                            */
      if (m == 1 || (m == 2 && mask[tri[1]] < 0)) {
	p = &d->voints[d->nvoints];
	p->x = x2;
	p->y = y2;
	if (m == 2) mask[tri[1]] = d->nvoints;
	mm = d->nvoints;
	d->nvoints++;
      } else
	mm = mask[tri[1]];
      d->vdges[2*i+1] = mm;
      if(i == dedge) 
	printf("second Voronoi edge = [%f %f]\n", 
	       d->voints[d->vdges[2*i+1]].x, d->voints[d->vdges[2*i+1]].y);

      if (filef) {
	fprintf(vf, "%f %f\n",x1, y1);
	fprintf(vf, "%f %f\n",x2, y2);
	fprintf(vf, "NaN NaN\n");
      }
    } else
      printf("Can't find Voronoi edge at %d\n", i);
  }
  if (filef) fclose(vf);

  /*
  if (params->uscf & US_TRI)
    convert_tri_mesh(params, params->d);
  if (params->uscf & US_HEX)
    convert_hex_mesh(params, params->d, 1);
  */
  if (DEBUG("init_m"))
    dlog("init_m", "JIGSAW mesh converted OK\n");
  convert_hex_mesh(params, params->d, 1);
  if (DEBUG("init_m"))
    dlog("init_m", "Voronoi mesh created OK\n");
}

/* END convert_jigsaw_msh()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computes the 'centre of a triangle                                */
/*-------------------------------------------------------------------*/
void tri_cen(int centref,     /* Type of triange centre              */
	     int geogf,       /* Geographic projection flag          */
	     double *p1,      /* [x,y,z] of vertex 1                 */
	     double *p2,      /* [x,y,z] of vertex 2                 */
	     double *p3,      /* [x,y,z] of vertex 3                 */
	     double *po       /* Coordinates of centre               */
	     )
{

  if (centref == 1) { /* Voronoi meshes                              */
    if (geogf) {
      /* Stereographic projections                                   */
      tria_ball_2s(po, p1, p2, p3, RADIUS);
    } else {
      /* Lat / long plane                                            */
      tria_ball_2d(po, p1, p2, p3); 
    }
  } else {            /* Power meshes                                */
    if (geogf) {
      /* Stereographic projections                                   */
      tria_ortho_2s(po, p1, p2, p3, RADIUS);
    } else {
      /* Lat / long plane                                            */
      tria_ortho_2d(po, p1, p2, p3);
    }
    /* If the centre lies outside the triangle then use the centre   */
    /* of mass.                                                      */
    if (!in_tri(po, p1, p2, p3))
      tria_com(po, p1, p2, p3);
  }
}

/* END tri_cen()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computes the 'centre of an edge of a triangle                     */
/*-------------------------------------------------------------------*/
void edge_cen(int centref,     /* Type of triange centre             */
	      int geogf,       /* Geographic projection flag         */
	      double *v1,      /* [x,y,z] of vertex 1                */
	      double *v2,      /* [x,y,z] of vertex 2                */
	      double *po       /* Coordinates of centre              */
	      )
{
  if (centref == 1) { /* Voronoi meshes                              */
    /* Voronoi meshes                                          */
    if (geogf) {
      /* Stereographic projections                                   */
      tria_edge_2s(po, v1, v2, RADIUS);
    } else {
      /* Lat / long plane                                            */
      tria_edge_2d(po, v1, v2);
    }
  } else {            /* Power meshes                                */
    if (geogf) {
      /* Stereographic projections                                   */
      tria_ortho_edge_2s(po, v1, v2, RADIUS);
    } else {
      /* Lat / long plane                                            */
      tria_ortho_edge_2d(po, v1, v2);
    }
  }
}

/* END edge_cen()                                                    */
/*-------------------------------------------------------------------*/


double is_obtuse(double *p0, double *p1, double *p2)
{
  double a, b, c, r = 0.0;
  double x1 = p0[0] - p1[0];
  double y1 = p0[1] - p1[1];
  a = sqrt(x1 * x1 + y1 * y1);
  x1 = p0[0] - p2[0];
  y1 = p0[1] - p2[1];
  b = sqrt(x1 * x1 + y1 * y1);
  x1 = p1[0] - p2[0];
  y1 = p1[1] - p2[1];
  c = sqrt(x1 * x1 + y1 * y1);

  if (a > b && a > c)
    r = (b * b + c * c - a * a) / (2.0 * b * c);
  else if (b > a && b > c)
    r = (a * a + c * c - b * b) / (2.0 * a * c);
  else
    r = (a * a + b * b - c * c) / (2.0 * a * b);
  r = acos(r);
  r = (r > PI/2.0) ? r * 180.0 / PI : 0.0;
  return(r);
}

      
int iswetc(unsigned long flag)
{
  if (!(flag & (SOLID|OUTSIDE)))
    return(1);
  else
    return(0);
}

int isdryc(unsigned long flag)
{
  if (flag & (SOLID|OUTSIDE))
    return(1);
  else
    return(0);
}

int find_mindex(double slon, double slat, double *lon, double *lat, int nsl, double *d)
{
  double dist, d1, d2;
  double dmin = HUGE;
  int i;
  int ret = -1;

  for (i = 1; i <= nsl; i++) {
    d1 = slon - lon[i];
    d2 = slat - lat[i];
    dist = sqrt(d1 * d1 + d2 * d2);
    if (dist < dmin) {
      dmin = dist;
      if (d) *d = dist;
      ret = i;
    }
  }
  return(ret);
}

void addpoint(parameters_t *params, 
	      int i, 
	      int j, 
	      int *n, 
	      point *pin,
	      int **mask)
{
  if (mask[j][i] < 0) {
    pin[*n].x = params->x[j*2][i*2];
    pin[*n].y = params->y[j*2][i*2];
    mask[j][i] = *n;
    *n += 1;
  }
}

int find_mesh_vertex(int c, double x, double y, double **xloc, double **yloc, int *mask, int ns2, int *npe)
{
  int cc, j;
  int found = 0;

  for (cc = 1; cc <= ns2; cc++) {
    if(c == cc) continue;
    if (mask[cc]) {
      for (j = 1; j <= npe[cc]; j++) {
	if (x == xloc[cc][j] && y == yloc[cc][j]) {
	  found = 1;
	}
      }
    }
  } 
 return(found);
}

/** Skip forward from the current file position to
  * the next line beginning with key, positioned at the character
  * immediately after the key.
  *
  * @param fp pointer to stdio FILE structure.
  * @param key keyname to locate in file.
  * @return non-zero if successful.
  */
int prm_skip_to_end_of_tok(FILE * fp, char *key, char *token, char *ret)
{
  char buf[MAXLINELEN], *tok;
  int len = strlen(key);
  long fpos;
  char *s;
  char *r;
  int rewound = 0;
  /*UR 10/06/2005
   * reduced failure of finding a key to 
   * 'TRACE'
   * since it is not vital and allowed. 
   */
  do {
    fpos = ftell(fp);
    s = fgets(buf, MAXLINELEN, fp);
    if (s == NULL) {
      if (!rewound) {
        fpos = 0L;
        if (fseek(fp, fpos, 0))
          quit("prm_skip_to_end_of_tok: %s\n", strerror(errno));
        s = fgets(buf, MAXLINELEN, fp);
        rewound = 1;
      } else {
        emstag(LTRACE,"lib:prmfile:prm_skip_to_end_of_tok"," key %s not found\n",
                         key);
        /*
        (*keyprm_errfn) ("prm_skip_to_end_of_key: key %s not found\n",
                         key);*/
        return (0);
      }
    }

    if (s == NULL)
      break;
    tok = strtok(s, token);

    /* Truncate the string at the first space after the key length. */
    if (strlen(s) > len && is_blank(s[len]))/*UR added length check to prevent invalid write*/
      s[len] = '\000';
  } while (strcasecmp(key, s) != 0);

  if (s == NULL) {
    emstag(LTRACE,"lib:prmfile:prm_skip_to_end_of_key","key %s not found\n", key);
    /*(*keyprm_errfn) ("prm_skip_to_end_of_key: key %s not found\n", key);*/
    return (0);
  }

  /* seek to character after key */
  if (fseek(fp, fpos + (s - buf) + len + 1, 0))
    quit("prm_skip_to_end_of_key: %s\n", strerror(errno));

  fgets(buf, MAXLINELEN, fp);

  /* Strip leading space */
  for (s = buf; *s && is_blank(*s); s++) /* loop */
    ;

  /* Copy out result */
  for (r = ret; *s && (*s != '\n'); *r++ = *s++)  /* loop */
    ;
  *r = 0;

  return (1);
}

void bisector(double x1a, 
	      double y1a, 
	      double x2a, 
	      double y2a,
	      double x1b, 
	      double y1b, 
	      double x2b, 
	      double y2b,
	      double *x,
	      double *y
	      )
{
  double mxa = 0.5 * (x1a + x2a);
  double mya = 0.5 * (x2a + y2a);
  double sa = -(x2a - x1a) / (y2a - y1a);
  double ia = mya - sa * mxa;
  double mxb = 0.5 * (x1b + x2b);
  double myb = 0.5 * (x2b + y2b);
  double sb = -(x2b - x1b) / (y2b - y1b);
  double ib = myb - sb * mxb;
  *x = (ia - ib) / (sb - sa);
  *y = sa * (*x) + ia;
}


/*------------------------------ unrolled matrix indexing */

#define uij(ir, ic, nr)  ((ir)+(ic)*(nr))

/*-------------------------------------------------------------------*/
/* Determinant of a 2x2 matrix                                       */
/* Contributed by Darren Engwirda (JIGSAW code).                     */
/*-------------------------------------------------------------------*/
double det_2x2(
	       int la ,    /* Leading dimension of aa                */
	       double *aa  /* Matrix                                 */
	       )
{   
  double ret;
  ret = aa[uij(0,0,la)] * aa[uij(1,1,la)] -
    aa[uij(0,1,la)] * aa[uij(1,0,la)] ;
  
  return(ret);
}

/* END det_2x2()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Determinant of a 3x3 matrix                                       */
/* Contributed by Darren Engwirda (JIGSAW code).                     */
/*-------------------------------------------------------------------*/
double det_3x3(
	       int la ,    /* Leading dimension of aa                */
	       double *aa  /* Matrix                                 */
	       )
{   
  double ret;
  ret = aa[uij(0,0,la)] * (aa[uij(1,1,la)] * 
			    aa[uij(2,2,la)] - 
			    aa[uij(1,2,la)] * 
			    aa[uij(2,1,la)] ) -

    aa[uij(0,1,la)] * (aa[uij(1,0,la)] * 
			aa[uij(2,2,la)] - 
			aa[uij(1,2,la)] * 
			aa[uij(2,0,la)] ) +

    aa[uij(0,2,la)] * (aa[uij(1,0,la)] * 
			aa[uij(2,1,la)] - 
			aa[uij(1,1,la)] * 
			aa[uij(2,0,la)] ) ;
  return(ret);
}

/* END det_3x3()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* 2-by-2 matrix inversion                                           */
/* Contributed by Darren Engwirda (JIGSAW code).                     */
/*-------------------------------------------------------------------*/
void inv_2x2(int la,     /* Leading dim. of aa                       */
	     double *aa, /* Matrix                                   */
	     int lx,     /* Leading dim. of xx                       */
	     double *xx, /* Matrix inversion * det                   */
	     double *da  /* matrix determinant                       */
	     )
{
  *da = det_2x2(la, aa);
    
  xx[uij(0,0,lx)] = aa[uij(1,1,la)];
  xx[uij(1,1,lx)] = aa[uij(0,0,la)];
  xx[uij(0,1,lx)] = -aa[uij(0,1,la)];
  xx[uij(1,0,lx)] = -aa[uij(1,0,la)];
}

/* END inv_2x2()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* 3-by-3 matrix inversion                                           */
/* Contributed by Darren Engwirda (JIGSAW code).                     */
/*-------------------------------------------------------------------*/
void inv_3x3(int la,     /* Leading dim. of aa                       */
	     double *aa, /* Matrix                                   */
	     int lx,     /* Leading dim. of xx                       */
	     double *xx, /* Matrix inversion * det                   */
	     double *da  /* matrix determinant                       */
	     )
{
  *da = det_3x3(la, aa);

  xx[uij(0,0,lx)] =
    aa[uij(2,2,la)] * 
    aa[uij(1,1,la)] - 
    aa[uij(2,1,la)] * 
    aa[uij(1,2,la)] ;
    
  xx[uij(0,1,lx)] =
    aa[uij(2,1,la)] * 
    aa[uij(0,2,la)] - 
    aa[uij(2,2,la)] * 
    aa[uij(0,1,la)] ;
  
  xx[uij(0,2,lx)] =
    aa[uij(1,2,la)] * 
    aa[uij(0,1,la)] - 
    aa[uij(1,1,la)] * 
    aa[uij(0,2,la)] ;
  
  xx[uij(1,0,lx)] =
    aa[uij(2,0,la)] * 
    aa[uij(1,2,la)] - 
    aa[uij(2,2,la)] * 
    aa[uij(1,0,la)] ;
    
  xx[uij(1,1,lx)] =
    aa[uij(2,2,la)] * 
    aa[uij(0,0,la)] - 
    aa[uij(2,0,la)] * 
    aa[uij(0,2,la)] ;
  
  xx[uij(1,2,lx)] =
    aa[uij(1,0,la)] * 
    aa[uij(0,2,la)] - 
    aa[uij(1,2,la)] * 
    aa[uij(0,0,la)] ;

  xx[uij(2,0,lx)] =
    aa[uij(2,1,la)] * 
    aa[uij(1,0,la)] - 
    aa[uij(2,0,la)] * 
    aa[uij(1,1,la)] ;

  xx[uij(2,1,lx)] =
    aa[uij(2,0,la)] * 
    aa[uij(0,1,la)] - 
    aa[uij(2,1,la)] * 
    aa[uij(0,0,la)] ;
    
  xx[uij(2,2,lx)] =
    aa[uij(1,1,la)] * 
    aa[uij(0,0,la)] - 
    aa[uij(1,0,la)] * 
    aa[uij(0,1,la)] ;

}

/* END inv_3x3()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Triangle circumcentre. Contributed by Darren Engwirda (JIGSAW     */
/* code).                                                            */
/*-------------------------------------------------------------------*/
void tria_ball_2d(double *bb,    /* centre: [x,y]                    */
		  double *p1,    /* vertex 1: [x,y]                  */
		  double *p2,    /* vertex 2: [x,y]                  */
		  double *p3     /* vertex 3: [x,y]                  */
		  )
{

  double xm[2*2];
  double xi[2*2];
  double xr[2*1];
  double dd;

  /* LHS matrix                                                      */
  xm[uij(0,0,2)] = p2[0]-p1[0] ;
  xm[uij(0,1,2)] = p2[1]-p1[1] ;
  xm[uij(1,0,2)] = p3[0]-p1[0] ;
  xm[uij(1,1,2)] = p3[1]-p1[1] ;

  /* RHS vector                                                      */
  xr[0] = (double)+0.5 * (xm[uij(0,0,2)] * xm[uij(0,0,2)] +
			 xm[uij(0,1,2)] * xm[uij(0,1,2)]);
        
  xr[1] = (double)+0.5 * (xm[uij(1,0,2)] * xm[uij(1,0,2)] +
			  xm[uij(1,1,2)] * xm[uij(1,1,2)]);

  /* Matrix inversion                                                */
  inv_2x2(2, xm, 2, xi, &dd) ;
    
  /* linear solver                                                   */
  bb[0] = (xi[uij(0,0,2)] * xr[0] + xi[uij(0,1,2)] * xr[1]);
  bb[1] = (xi[uij(1,0,2)] * xr[0] + xi[uij(1,1,2)] * xr[1]);
  bb[0] /= dd ;
  bb[1] /= dd ;
  
  /* Offset                                                          */    
  bb[0] += p1[0] ;
  bb[1] += p1[1] ;
}

/* END tri_ball_2d()                                                 */
/*-------------------------------------------------------------------*/	


/*-------------------------------------------------------------------*/	
/* Routines to compute the 'centres' for edges. This is the midpoint */
/* if circumcentres are used.                                        */
/* Contributed by Darren Engwirda (JIGSAW code).                     */
/*-------------------------------------------------------------------*/
void tria_edge_2d(double *bb,    /* centre: [x,y]                    */
		  double *p1,    /* vertex 1: [x,y]                  */
		  double *p2     /* vertex 2: [x,y]                  */
		  )
{
  bb[0] = 0.5 * (p1[0] + p2[0]);
  bb[1] = 0.5 * (p1[1] + p2[1]);
}

/* END tria_edge_2d()                                                */
/*-------------------------------------------------------------------*/	


/*-------------------------------------------------------------------*/
/* Triangle circumcentre on stereographic projections.               */
/* Contributed by Darren Engwirda (JIGSAW code).                     */
/*-------------------------------------------------------------------*/
void tria_ball_2s(double *bb,    /* centre: [x,y] (deg)              */
		  double *p1,    /* vert 1: [l,p] (rad)              */
		  double *p2,    /* vert 2: [l,p] (rad)              */
		  double *p3,    /* vert 3: [l,p] (rad)              */
		  double  rr     /* radius of sphere                 */
		  )
{
  double  e1[3] ;
  double  e2[3] ;
  double  e3[3] ;
  double  eb[3] ;
    
  lonlat_to_xyz(p1, e1, rr) ;
  lonlat_to_xyz(p2, e2, rr) ;
  lonlat_to_xyz(p3, e3, rr) ;
        
  tria_circ_ball_3d(eb, e1, e2, e3) ;

  xyz_to_lonlat(eb, bb, rr) ; 

}

/* END void tria_ball_2d()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/	
/* Routines to compute the 'centres' for edges on stereographic      */
/* projections. This is the midpoint if circumcentres are used.      */
/* Contributed by Darren Engwirda (JIGSAW code).                     */
/*-------------------------------------------------------------------*/
void tria_edge_2s(double *bb,    /* centre: [x,y] (deg)              */
		  double *p1,    /* vert 1: [l,p] (rad)              */
		  double *p2,    /* vert 2: [l,p] (rad)              */
		  double  rr     /* radius of sphere                 */
		  )
{
  double  e1[3] ;
  double  e2[3] ;
  double  eb[3] ;
    
  lonlat_to_xyz(p1, e1, rr) ;
  lonlat_to_xyz(p2, e2, rr) ;
        
  edge_circ_ball_3d(eb, e1, e2) ;
    
  xyz_to_lonlat(eb, bb, rr) ; 

}

/* END tria_edge_2s()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Triangle circumcentre in 3 dimensions.                            */
/* Contributed by Darren Engwirda (JIGSAW code).                     */
/*-------------------------------------------------------------------*/
void tria_circ_ball_3d(double *bb,    /* centre: [x,y,z]             */
		       double *p1,    /* vert 1: [x,y,z]             */
		       double *p2,    /* vert 2: [x,y,z]             */
		       double *p3     /* vert 3: [x,y,z]             */
		       )
{
  double  nv[3*1] ;
  double  xm[3*3] ;
  double  xi[3*3] ;
  double  xr[3*1] ;
  double  ee[3*1] ;
  double  db[3*1] ;
  double  dd ;
  int ii;
       
  /*------------------------------------ LHS matrix                  */
  xm[uij(0,0,3)] = p2[0]-p1[0] ;
  xm[uij(0,1,3)] = p2[1]-p1[1] ;
  xm[uij(0,2,3)] = p2[2]-p1[2] ;
        
  xm[uij(1,0,3)] = p3[0]-p1[0] ;
  xm[uij(1,1,3)] = p3[1]-p1[1] ;
  xm[uij(1,2,3)] = p3[2]-p1[2] ;
        
  tria_norm_vect_3d(nv, p1, p2, p3);
        
  xm[uij(2,0,3)] = nv[0] ;
  xm[uij(2,1,3)] = nv[1] ;
  xm[uij(2,2,3)] = nv[2] ;

  /*------------------------------------ RHS vector                  */
  xr[0] = (double)+.5 * (xm[uij(0,0,3)] * 
			 xm[uij(0,0,3)] +
			 xm[uij(0,1,3)] *
			 xm[uij(0,1,3)] +
			 xm[uij(0,2,3)] *
			 xm[uij(0,2,3)] ) ;

  xr[1] = (double)+.5 * (xm[uij(1,0,3)] * 
			 xm[uij(1,0,3)] +
			 xm[uij(1,1,3)] *
			 xm[uij(1,1,3)] +
			 xm[uij(1,2,3)] *
			 xm[uij(1,2,3)] ) ;
        
  xr[2] = (double)+.0 ;

  /*------------------------------------ matrix inv                  */
  inv_3x3(3, xm, 3, xi, &dd) ;
        
  /*------------------------------------ lin-solver                  */
  bb[0] = (xi[uij(0,0,3)] * xr[0] +
	   xi[uij(0,1,3)] * xr[1] +
	   xi[uij(0,2,3)] * xr[2] )  ;
        
  bb[1] = (xi[uij(1,0,3)] * xr[0] +
	   xi[uij(1,1,3)] * xr[1] +
	   xi[uij(1,2,3)] * xr[2] )  ;

  bb[2] = (xi[uij(2,0,3)] * xr[0] +
	   xi[uij(2,1,3)] * xr[1] +
	   xi[uij(2,2,3)] * xr[2] )  ;

  bb[0] /= dd ;
  bb[1] /= dd ;
  bb[2] /= dd ;

  /*------------------------------------ iter. ref.                  */
  for(ii = 1; ii-- != 0; ) {
    ee[0] = xr[0] - (xm[uij(0,0,3)] * bb[0] +
		     xm[uij(0,1,3)] * bb[1] +
		     xm[uij(0,2,3)] * bb[2] )  ;
        
    ee[1] = xr[1] - (xm[uij(1,0,3)] * bb[0] +
		     xm[uij(1,1,3)] * bb[1] +
		     xm[uij(1,2,3)] * bb[2] )  ;
        
    ee[2] = xr[2] - (xm[uij(2,0,3)] * bb[0] +
		     xm[uij(2,1,3)] * bb[1] +
		     xm[uij(2,2,3)] * bb[2] )  ;
        
    db[0] = (xi[uij(0,0,3)] * ee[0] +
	     xi[uij(0,1,3)] * ee[1] +
	     xi[uij(0,2,3)] * ee[2] )  ;
        
    db[1] = (xi[uij(1,0,3)] * ee[0] +
	     xi[uij(1,1,3)] * ee[1] +
	     xi[uij(1,2,3)] * ee[2] )  ;

    db[2] = (xi[uij(2,0,3)] * ee[0] +
	     xi[uij(2,1,3)] * ee[1] +
	     xi[uij(2,2,3)] * ee[2] )  ;
        
    bb[0] += db[0] / dd ;
    bb[1] += db[1] / dd ;
    bb[2] += db[2] / dd ;     
  }

  /*------------------------------------ offset mid                  */
  bb[0] += p1[0] ;
  bb[1] += p1[1] ;
  bb[2] += p1[2] ;
}

/* END tria_circ_ball_3d()                                           */
/*-------------------------------------------------------------------*/	


/*-------------------------------------------------------------------*/	
/* Routines to compute the 'centres' for edges i three dimensions.   */
/* This is the midpoint if circumcentres are used.                   */
/*-------------------------------------------------------------------*/    
void edge_circ_ball_3d (double *bb,    /* centre: [x,y,z]            */
			double *p1,    /* vert 1: [x,y,z]            */
			double *p2     /* vert 2: [x,y,z]            */
			)
{
  bb[0] = (double)+.5 * (p1[0] + p2[0]) ;
  bb[1] = (double)+.5 * (p1[1] + p2[1]) ;
  bb[2] = (double)+.5 * (p1[2] + p2[2]) ;
}

/* END edge_circ_ball_3d()                                           */
/*-------------------------------------------------------------------*/	


/*-------------------------------------------------------------------*/
/* Triangle ortho-centre. Contributed by Darren Engwirda (JIGSAW     */
/* code).                                                            */
/*-------------------------------------------------------------------*/
void tria_ortho_2d(double *bb,    /* centre: [x,y]                   */
		   double *p1,    /* vertex 1: [x,y,w]               */
		   double *p2,    /* vertex 2: [x,y,w]               */
		   double *p3     /* vertex 3: [x,y,w]               */
		   )
{
  int ii;
  double xm[2*2];
  double xi[2*2];
  double xr[2*1];
  double ee[2*1];
  double db[2*1];
  double dd;

  /* LHS matrix                                                      */
  xm[uij(0,0,2)] = p2[0]-p1[0] ;
  xm[uij(0,1,2)] = p2[1]-p1[1] ;
  xm[uij(1,0,2)] = p3[0]-p1[0] ;
  xm[uij(1,1,2)] = p3[1]-p1[1] ;

  /* RHS vector                                                      */
  xr[0] = (double)+0.5 * (xm[uij(0,0,2)] * xm[uij(0,0,2)] +
			 xm[uij(0,1,2)] * xm[uij(0,1,2)]);
        
  xr[1] = (double)+0.5 * (xm[uij(1,0,2)] * xm[uij(1,0,2)] +
			  xm[uij(1,1,2)] * xm[uij(1,1,2)]);


  xr[0] -= (double)+0.5 * (p2[2] - p1[2]);
  xr[1] -= (double)+0.5 * (p3[2] - p1[2]);

  /* Matrix inversion                                                */
  inv_2x2(2, xm, 2, xi, &dd) ;
    
  /* linear solver                                                   */
  bb[0] = (xi[uij(0,0,2)] * xr[0] + xi[uij(0,1,2)] * xr[1]);
  bb[1] = (xi[uij(1,0,2)] * xr[0] + xi[uij(1,1,2)] * xr[1]);
  bb[0] /= dd ;
  bb[1] /= dd ;

  /* Iteration                                                       */
  for (ii = 1; ii-- != 0;) {
    ee[0] = xr[0] - (xm[uij(0,0,2)] * bb[0] + xm[uij(0,1,2)] * bb[1]);
    ee[1] = xr[1] - (xm[uij(1,0,2)] * bb[0] + xm[uij(1,1,2)] * bb[1]);
        
    db[0] = (xi[uij(0,0,2)] * ee[0] + xi[uij(0,1,2)] * ee[1]);
    db[1] = (xi[uij(1,0,2)] * ee[0] + xi[uij(1,1,2)] * ee[1]);

    bb[0] += db[0] / dd;
    bb[1] += db[1] / dd;
  }

  /* Offset                                                          */    
  bb[0] += p1[0] ;
  bb[1] += p1[1] ;
  
  /*------------------------------------ iter. ref.                  */
  for(ii = 1; ii-- != 0; ) {
    ee[0] = xr[0] - (xm[uij(0,0,2)] * bb[0] +
		     xm[uij(0,1,2)] * bb[1] )  ;
        
    ee[1] = xr[1] - (xm[uij(1,0,2)] * bb[0] +
		     xm[uij(1,1,2)] * bb[1] )  ;
        
    db[0] = (xi[uij(0,0,2)] * ee[0] +
	     xi[uij(0,1,2)] * ee[1] )  ;
        
    db[1] = (xi[uij(1,0,2)] * ee[0] + xi[uij(1,1,2)] * ee[1] )  ;

    bb[0] += db[0] / dd ;
    bb[1] += db[1] / dd ;
  }
    
  /*------------------------------------ offset mid                  */
  bb[0] += p1[0] ;
  bb[1] += p1[1] ;
}

/* END tri_orth_2d()                                                 */
/*-------------------------------------------------------------------*/	


/*-------------------------------------------------------------------*/	
/* Routines to compute the 'centres' for edges. This is the midpoint */
/* if circumcentres are used.                                        */
/*-------------------------------------------------------------------*/
void tria_ortho_edge_2d(double *bb,    /* centre: [x,y]              */
			double *p1,    /* vertex 1: [x,y,w]          */
			double *p2     /* vertex 2: [x,y,w]          */
			)
{
  double dp, tt, dd[3*1] ;

  /*------------------------------------ DoF deltas */

  dd[0] = p1[0] - p2[0] ;
  dd[1] = p1[1] - p2[1] ;
  dd[2] = p1[2] - p2[2] ;
        
  dp = dd[0] * dd[0] + dd[1] * dd[1] ;

  /*------------------------------------ line param */
  tt = (double)+.5 * (dd[2] + dp) / dp ;
     
  /*------------------------------------ interp. 1d */           
  bb[0] = p1[0] - tt * dd[0] ;
  bb[1] = p1[1] - tt * dd[1] ; 
}

/* END tria_orth_edge_2d()                                           */
/*-------------------------------------------------------------------*/	


/*-------------------------------------------------------------------*/
/* Triangle ortho-centre for stereographic projections.              */
/* Contributed by Darren Engwirda (JIGSAW code).                     */
/*-------------------------------------------------------------------*/
void tria_ortho_2s(double *bb,    /* centre: [l,p] (rad)             */
		   double *p1,    /* vert 1: [l,p,w] (rad)           */
		   double *p2,    /* vert 2: [l,p,w] (rad)           */
		   double *p3,    /* vert 3: [l,p,w] (rad)           */
		   double  rr     /* radius of sphere                */
		   )   
{
  double  e1[4] ;
  double  e2[4] ;
  double  e3[4] ;
  double  eb[3] ;
    
  lonlat_to_xyz(p1, e1, rr) ;
  lonlat_to_xyz(p2, e2, rr) ;
  lonlat_to_xyz(p3, e3, rr) ;
   
  e1[3] = p1[2] ;
  e2[3] = p2[2] ;
  e3[3] = p3[2] ;
        
  tria_orth_ball_3d(eb, e1, e2, e3) ;
    
  xyz_to_lonlat(eb, bb, rr) ; 
}

/* END tria_ortho_2s()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/	
/* Routines to compute the 'centres' for edges in stereographic      */
/* projections. This is the midpoint if circumcentres are used.      */
/*-------------------------------------------------------------------*/
void tria_ortho_edge_2s(double *bb,    /* centre: [l,p] (rad)        */
			double *p1,    /* vert 1: [l,p,w] (rad)      */
			double *p2,    /* vert 2: [l,p,w] (rad)      */
			double  rr     /* radius of sphere           */
			)
{
  double  e1[4] ;
  double  e2[4] ;
  double  eb[3] ;
    
  lonlat_to_xyz(p1, e1, rr) ;
  lonlat_to_xyz(p2, e2, rr) ;
  
  e1[3] = p1[2] ;
  e2[3] = p2[2] ;
        
  edge_orth_ball_3d(eb, e1, e2) ;
    
  xyz_to_lonlat(eb, bb, rr) ; 
}

/* END tria_ortho_edge_2s()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/	
/* Routines to compute the 'centres' for edges in 3 dimensions. This */
/* is the midpoint if circumcentres are used.                        */
/*-------------------------------------------------------------------*/
void edge_orth_ball_3d(double *bb,    /* centre: [x,y,z]             */
		       double *p1,    /* vert 1: [x,y,z,w]           */
		       double *p2     /* vert 2: [x,y,z,w]           */
		       )   
{
  double dp, tt, dd[3*1] ;

  /*------------------------------------ DoF deltas                  */
  dd[0] = p1[0] - p2[0] ;
  dd[1] = p1[1] - p2[1] ;
  dd[2] = p1[2] - p2[2] ;
  dd[3] = p1[3] - p2[3] ;
        

  dp = dd[0] * dd[0] + dd[1] * dd[1] + dd[2] * dd[2] ;

  /*------------------------------------ line param                  */
  tt = (double)+.5 * ( dd[3] + dp ) / dp ;
     
  /*------------------------------------ interp. 1d                  */
  bb[0] = p1[0] - tt * dd[0] ;
  bb[1] = p1[1] - tt * dd[1] ;
  bb[2] = p1[2] - tt * dd[2] ; 
}
    
/* END edge_orth_ball_3d()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Triangle ortho-centre in three dimensions.                        */
/* Contributed by Darren Engwirda (JIGSAW code).                     */
/*-------------------------------------------------------------------*/
void tria_orth_ball_3d(double *bb,     /* centre: [x,y,z]            */
		       double *p1,     /* vert 1: [x,y,z,w]          */
		       double *p2,     /* vert 2: [x,y,z,w]          */
		       double *p3      /* vert 3: [x,y,z,w]          */
		       )
{
  double  nv[3*1] ;
  double  xm[3*3] ;
  double  xi[3*3] ;
  double  xr[3*1] ;
  double  ee[3*1] ;
  double  db[3*1] ;
  double  dd, w21, w31 ;
  int ii;
    
  /*------------------------------------ LHS matrix                  */
  xm[uij(0,0,3)] = p2[0] - p1[0] ;
  xm[uij(0,1,3)] = p2[1] - p1[1] ;
  xm[uij(0,2,3)] = p2[2] - p1[2] ;
        
  xm[uij(1,0,3)] = p3[0] - p1[0] ;
  xm[uij(1,1,3)] = p3[1] - p1[1] ;
  xm[uij(1,2,3)] = p3[2] - p1[2] ;
        
  tria_norm_vect_3d(nv, p1, p2, p3);
        
  xm[uij(2,0,3)] = nv[0] ;
  xm[uij(2,1,3)] = nv[1] ;
  xm[uij(2,2,3)] = nv[2] ;

  /*------------------------------------ RHS vector                  */
  xr[0] = (double)+.5 * (xm[uij(0,0,3)] * 
			 xm[uij(0,0,3)] +
			 xm[uij(0,1,3)] *
			 xm[uij(0,1,3)] +
			 xm[uij(0,2,3)] *
			 xm[uij(0,2,3)] ) ;

  xr[1] = (double)+.5 * (xm[uij(1,0,3)] * 
			 xm[uij(1,0,3)] +
			 xm[uij(1,1,3)] *
			 xm[uij(1,1,3)] +
			 xm[uij(1,2,3)] *
			 xm[uij(1,2,3)] ) ;
        
  w21 = p2[3] - p1[3] ;
  w31 = p3[3] - p1[3] ;
  
  xr[0]-= (double)+.5 * w21 ;
  xr[1]-= (double)+.5 * w31 ;
  xr[2] = (double)+.0 ;

  /*------------------------------------ matrix inv                  */
  inv_3x3(3, xm, 3, xi,&dd) ;
        
  /*------------------------------------ lin-solver                  */
  bb[0] = (xi[uij(0,0,3)] * xr[0] +
	   xi[uij(0,1,3)] * xr[1] +
	   xi[uij(0,2,3)] * xr[2] )  ;
        
  bb[1] = (xi[uij(1,0,3)] * xr[0] +
	   xi[uij(1,1,3)] * xr[1] +
	   xi[uij(1,2,3)] * xr[2] )  ;

  bb[2] = (xi[uij(2,0,3)] * xr[0] +
	   xi[uij(2,1,3)] * xr[1] +
	   xi[uij(2,2,3)] * xr[2] )  ;

  bb[0] /= dd ;
  bb[1] /= dd ;
  bb[2] /= dd ;

  /*------------------------------------ iter. ref.                  */
  for(ii = +1; ii-- != +0; ) {
    ee[0] = xr[0] - (xm[uij(0,0,3)] * bb[0] +
		     xm[uij(0,1,3)] * bb[1] +
		     xm[uij(0,2,3)] * bb[2] )  ;
        
    ee[1] = xr[1] - (xm[uij(1,0,3)] * bb[0] +
		     xm[uij(1,1,3)] * bb[1] +
		     xm[uij(1,2,3)] * bb[2] )  ;
        
    ee[2] = xr[2] - ( xm[uij(2,0,3)] * bb[0] +
		      xm[uij(2,1,3)] * bb[1] +
		      xm[uij(2,2,3)] * bb[2] )  ;
        
    db[0] = (xi[uij(0,0,3)] * ee[0] +
	     xi[uij(0,1,3)] * ee[1] +
	     xi[uij(0,2,3)] * ee[2] )  ;
        
    db[1] = (xi[uij(1,0,3)] * ee[0] +
	     xi[uij(1,1,3)] * ee[1] +
	     xi[uij(1,2,3)] * ee[2] )  ;

    db[2] = (xi[uij(2,0,3)] * ee[0] +
	     xi[uij(2,1,3)] * ee[1] +
	     xi[uij(2,2,3)] * ee[2] )  ;
    
    bb[0] += db[0] / dd ;
    bb[1] += db[1] / dd ;
    bb[2] += db[2] / dd ;
  }

  /*------------------------------------ offset mid                  */
  bb[0] += p1[0] ;
  bb[1] += p1[1] ;
  bb[2] += p1[2] ;
}

/* END tria_orth_ball_3d()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Coord. transformation (S^2 <=> E^3)                               */
/*-------------------------------------------------------------------*/
void lonlat_to_xyz(double *pp,     /* lonlat: [lat,lon] (deg)        */
		   double *ee,     /* euclid: [x,y,z]                */
		   double  rr      /* radius of sphere               */
		   )
{
  ee[0] = rr * cos(DEG2RAD(pp[0])) * cos(DEG2RAD(pp[1])) ;
  ee[1] = rr * sin(DEG2RAD(pp[0])) * cos(DEG2RAD(pp[1])) ;
  ee[2] = rr * sin(DEG2RAD(pp[1])) ;
}

/* END lonlat_to_xyz ()                                              */
/*-------------------------------------------------------------------*/
    
/*-------------------------------------------------------------------*/
/* Coord. transformation (S^2 <=> E^3)                               */
/*-------------------------------------------------------------------*/
void xyz_to_lonlat(double *ee,    /* euclid: [x,y,z]                 */
		   double *pp,    /* lonlat: [lon,lat] (deg)         */
		   double rr      /* radius of sphere                */
		   )
{
  pp[1] = RAD2DEG(asin(ee[2]/rr)) ;
  pp[0] = RAD2DEG(atan2(ee[1], ee[0])) ;
}

/* END xyz_to_lonlat()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computes the centre of mass of a triangle                         */
/*-------------------------------------------------------------------*/
void tria_com(double *bb,         /* centre: [x,y]                   */
	      double *p1,         /* vertex 1: [x,y]                 */
	      double *p2,         /* vertex 2: [x,y]                 */
	      double *p3          /* vertex 3: [x,y]                 */
	      )
{
  int n;

  bb[0] = 0.0;  bb[1] = 0.0;
  for (n = 0; n < 2; n++) {
    bb[n] = (p1[n] + p2[n] + p3[n]) / 3.0;
  }
}

/* END tri_com()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*  Calculate normal vector (TRIA, E^3)                              */
/*-------------------------------------------------------------------*/
void tria_norm_vect_3d(double *nv,    /* normal: [x,y,z]             */
		       double *pa,    /* vert 1: [x,y,z]             */
		       double *pb,    /* vert 2: [x,y,z]             */
		       double *pc     /* vert 3: [x,y,z]             */
		       )
{
  double  ab[3] ;
  double  ac[3] ;

  /*------------------------------------ edge vect.                  */
  ab[0] = pb[0] - pa[0] ;
  ab[1] = pb[1] - pa[1] ;
  ab[2] = pb[2] - pa[2] ;
  
  ac[0] = pc[0] - pa[0] ;
  ac[1] = pc[1] - pa[1] ;
  ac[2] = pc[2] - pa[2] ;

  /*------------------------------------ norm vect.                  */
  nv[0] = ab[1] * ac[2] - ab[2] * ac[1] ;
  nv[1] = ab[2] * ac[0] - ab[0] * ac[2] ;
  nv[2] = ab[0] * ac[1] - ab[1] * ac[0] ;
}

/* END tria_norm_vect_3d()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Returns 1 if a point lies within a triangle. Uses Baycentric      */
/* coordinates to find if this is so. The area is > 0 for vertices   */
/* ordered in an anti-clockwise sense, and negative for a clockwise  */
/* sense.                                                            */
/*-------------------------------------------------------------------*/
int in_tri(double *po, double *p0, double *p1, double *p2) {
  double p0x = p0[0];
  double p0y = p0[1];
  double p1x = p1[0];
  double p1y = p1[1];
  double p2x = p2[0];
  double p2y = p2[1];
  double area = 0.5 * (-p1y*p2x + p0y*(-p1x + p2x) + p0x*(p1y - p2y) + p1x*p2y);
  double ap = 2.0*fabs(area);
  double s = (p0y*p2x - p0x*p2y + (p2y - p0y)*po[0] + (p0x - p2x)*po[1])/(2.0*area);
  double t = (p0x*p1y - p0y*p1x + (p0y - p1y)*po[0] + (p1x - p0x)*po[1])/(2.0*area);
  double sp = s * ap;
  double tp = t * ap;
  int ret = (sp > 0.0 & tp > 0.0 & sp + tp < ap) ? 1 : 0;
  return(ret);
}

/* END intri()                                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Edge sort function for dual row / index sorting                   */
/* Written by Darren Engwirda, 9/2017                                */
/* return -1 if "edge-1" < "edge-2"                                  */
/* return +1 if "edge-1" > "edge-2"                                  */
/* return +0 if "edge-1" = "edge-2"                                  */
/*-------------------------------------------------------------------*/
int edge_sort_compare (void const *x1, void const *x2)
{
  int  const *e1 = (int *)x1 ;
  int  const *e2 = (int *)x2 ;
    
  if (e1[0]  < e2[0]) {         /* less on 1st col. */
    return -1 ;
  } else {
    if (e1[0]  > e2[0]) {       /* more on 1st col. */
      return +1 ;
    } else {
      if (e1[0] == e2[0]) {     /* test on 2nd col. */
	if (e1[1]  < e2[1]) {   /* less on 2nd col. */
	  return -1 ;
	} else {
	  if (e1[1]  > e2[1])   /* more on 2nd col. */
	    return +1 ;
	}
      }
    }
  }
  return +0 ;                   /* have exact match */
}

int edge_sort_double (void const *x1, void const *x2)
{
  double  const *e1 = (double *)x1 ;
  double  const *e2 = (double *)x2 ;
    
  if (e1[0]  < e2[0]) {         /* less on 1st col. */
    return -1 ;
  } else {
    if (e1[0]  > e2[0]) {       /* more on 1st col. */
      return +1 ;
    } else {
      if (e1[0] == e2[0]) {     /* test on 2nd col. */
	if (e1[1]  < e2[1]) {   /* less on 2nd col. */
	  return -1 ;
	} else {
	  if (e1[1]  > e2[1])   /* more on 2nd col. */
	    return +1 ;
	}
      }
    }
  }
  return +0 ;                   /* have exact match */
}

/* END edge_sort_compare()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to find the neighbors in the mesh structure.              */
/* Uses O(nlog(n)) operations.                                       */
/* Originally written by Darren Engwirda.                            */
/*-------------------------------------------------------------------*/
void neighbour_finder(mesh_t *mesh)
{
  int cc, c1, c2, j, j1, j2, nedge, next;
  int **edge;

  /* Allocate                                                        */
  nedge = 0;
  for (cc = 1; cc <= mesh->ns2; cc++) {
    nedge += mesh->npe[cc];
  }
  edge = i_alloc_2d(4, nedge);

  /* Populate the edge structure                                     */
  next = 0;
  for (cc = 1; cc <= mesh->ns2; cc++) {
    for (j = 1; j <= mesh->npe[cc]; j++) {
      j1 = mesh->eloc[0][cc][j];
      j2 = mesh->eloc[1][cc][j];
      if (j1 < j2) {
	edge[next][0] = j1;
	edge[next][1] = j2;
      } else {
	edge[next][0] = j2;
	edge[next][1] = j1;
      }
      edge[next][2] = cc;
      edge[next][3] = j;
      next++;
    }
  }

  /* Sort edges, duplicates are consecutive!                         */
  qsort(edge[0], nedge, sizeof(int)*4, edge_sort_compare);

  /* Create the mappings                                             */
  mesh->neic = i_alloc_2d(mesh->ns2+1, mesh->mnpe+1);
  mesh->neij = i_alloc_2d(mesh->ns2+1, mesh->mnpe+1);
  for (cc = 1; cc <= mesh->ns2; cc++) {
    for (j = 1; j <= mesh->npe[cc]; j++) {
      mesh->neic[j][cc] = 0;
      mesh->neij[j][cc] = j;
    }
  }

  /* Create the mappings                                             */
  for (cc = 0; cc < nedge-1; cc++) {
    if ((edge[cc][0] == edge[cc+1][0]) &&
	(edge[cc][1] == edge[cc+1][1])) {
      c1 = edge[cc][2];
      c2 = edge[cc+1][2];
      j1 = edge[cc][3];
      j2 = edge[cc+1][3];
      mesh->neic[j1][c1] = c2;
      mesh->neij[j1][c1] = j2;
      mesh->neic[j2][c2] = c1;
      mesh->neij[j2][c2] = j1;
    }
  }
  i_free_2d(edge);
}    

/* Neighbour maps from delaunay structure */
void neighbour_finder_b(delaunay* d, int ***neic)
{
  int cc, c1, c2, j, next, n, nn, j1, j2;
  int **edge, **eic;
  int mnpe = 3;
  int ns2 = d->ntriangles;
  int nedge = ns2 * 3;
  triangle *t;

  edge = i_alloc_2d(4, nedge);
  /* Populate the edge structure                                     */
  next = 0;
  for (cc = 0; cc < ns2; cc++) {
    t = &d->triangles[cc];
    for (n = 0; n < 3; n++) {
      nn = (n == 2) ? 0 : n + 1;
      j1 = t->vids[n];
      j2 = t->vids[nn];
      if (j1 < j2) {
	edge[next][0] = j1;
	edge[next][1] = j2;
      } else {
	edge[next][0] = j2;
	edge[next][1] = j1;
      }
      edge[next][2] = cc;
      edge[next][3] = n;
      next++;
    }
  }

  /* Sort edges, duplicates are consecutive!                         */
  qsort(edge[0], nedge, sizeof(int)*4, edge_sort_compare);

  /* Create the mappings                                             */
  eic = i_alloc_2d(ns2, mnpe);
  for (cc = 0; cc < ns2; cc++) {
    for (j = 0; j < mnpe; j++) {
      eic[j][cc] = -1;
    }
  }

  /* Create the mappings                                             */
  for (cc = 0; cc < nedge-1; cc++) {
    /*printf("%d %d %d\n",cc,edge[cc][0],edge[cc][1]);*/
    if ((edge[cc][0] == edge[cc+1][0]) &&
	(edge[cc][1] == edge[cc+1][1])) {
      c1 = edge[cc][2];
      c2 = edge[cc+1][2];
      j1 = edge[cc][3];
      j2 = edge[cc+1][3];
      eic[j1][c1] = c2;
      eic[j2][c2] = c1;
    }
  }
  *neic = eic;
  i_free_2d(edge);
}

#ifdef HAVE_JIGSAWLIB
/* Neighbour maps from JIGSAW .msh input file */
void neighbour_finder_j(jigsaw_msh_t *msh, int ***neic)
{
  int cc, c1, c2, j, next, n, nn, j1, j2;
  int **edge, **eic, *nv, **eiv, meiv;
  int mnpe = 3;
  int ns2 = msh->_tria3._size;
  int nv2 = msh->_vert2._size;
  int nedge = ns2 * 3;

  edge = i_alloc_2d(4, nedge);
  nv = i_alloc_1d(nv2);
  memset(nv, 0, nv2 * sizeof(int));
  /* Populate the edge structure                                     */
  next = 0;
  for (cc = 0; cc < ns2; cc++) {
    for (n = 0; n < 3; n++) {
      nn = (n == 2) ? 0 : n + 1;
      j1 = msh->_tria3._data[cc]._node[n];
      j2 = msh->_tria3._data[cc]._node[nn];
      nv[j1]++;
      if (j1 < j2) {
	edge[next][0] = j1;
	edge[next][1] = j2;
      } else {
	edge[next][0] = j2;
	edge[next][1] = j1;
      }
      edge[next][2] = cc;
      edge[next][3] = n;
      next++;
    }
  }
  meiv = 0;
  for (cc = 0; cc < nv2; cc++)
    if (nv[cc] > meiv) meiv = nv[cc];
  eiv = i_alloc_2d(nv2, meiv);

  /* Sort edges, duplicates are consecutive!                         */
  qsort(edge[0], nedge, sizeof(int)*4, edge_sort_compare);

  /* Create the mappings                                             */
  eic = i_alloc_2d(ns2, mnpe);
  for (cc = 0; cc < ns2; cc++) {
    for (j = 0; j < mnpe; j++) {
      eic[j][cc] = -1;
    }
  }

  /* Create the mappings                                             */
  for (cc = 0; cc < nedge-1; cc++) {
    /*printf("%d %d %d\n",cc,edge[cc][0],edge[cc][1]);*/
    if ((edge[cc][0] == edge[cc+1][0]) &&
	(edge[cc][1] == edge[cc+1][1])) {
      c1 = edge[cc][2];
      c2 = edge[cc+1][2];
      j1 = edge[cc][3];
      j2 = edge[cc+1][3];
      eic[j1][c1] = c2;
      eic[j2][c2] = c1;
    }
  }
  *neic = eic;
  i_free_2d(edge);
}
#endif

/* Neighbour maps from Voronoi edges set in convert_hex_mesh() */
void neighbour_finder_a(int ns2, int *npe, int **vedge, int **neic, int **neij)
{
  int cc, c1, c2, j, j1, j2, nedge, next;
  int **edge, **eic, **eij, mnpe;

  /* Allocate                                                        */
  nedge = mnpe = 0;
  for (cc = 0; cc < ns2; cc++) {
    nedge += npe[cc];
    mnpe = max(mnpe, npe[cc]);
  }
  edge = i_alloc_2d(4, nedge);

  /* Populate the edge structure                                     */
  next = 0;
  for (cc = 0; cc < ns2; cc++) {
    for (j = 0; j < npe[cc]; j++) {
      j1 = vedge[cc][j];
      j2 = (j1 == npe[cc]-1) ? 0: j1 + 1;
      if (j1 < j2) {
	edge[next][0] = j1;
	edge[next][1] = j2;
      } else {
	edge[next][0] = j2;
	edge[next][1] = j1;
      }
      edge[next][2] = cc;
      edge[next][3] = j;
      next++;
    }
  }

  /* Sort edges, duplicates are consecutive!                         */
  qsort(edge[0], nedge, sizeof(int)*4, edge_sort_compare);

  /* Create the mappings                                             */
  eic = i_alloc_2d(ns2+1, mnpe+1);
  eij = i_alloc_2d(ns2+1, mnpe+1);
  for (cc = 1; cc <= ns2; cc++) {
    for (j = 1; j <= npe[cc]; j++) {
      eic[j][cc] = -1;
      eij[j][cc] = -1;
    }
  }

  /* Create the mappings                                             */
  for (cc = 0; cc < nedge-1; cc++) {
    if ((edge[cc][0] == edge[cc+1][0]) &&
	(edge[cc][1] == edge[cc+1][1])) {
      c1 = edge[cc][2]+1;
      c2 = edge[cc+1][2]+1;
      j1 = edge[cc][3]+1;
      j2 = edge[cc+1][3]+1;
      eic[j1][c1] = c2;
      eij[j1][c1] = j2;
      eic[j2][c2] = c1;
      eij[j2][c2] = j1;
    }
  }
  neic = eic;
  neij = eij;
}

/* END neighbour_finder()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Removes duplicate locations in a list. O(nlog(n)) operations.     */
/*-------------------------------------------------------------------*/
void remove_duplicates(int ns2, double **x, double **y, mesh_t *mesh)
{
  int cc, j, jj, nn, n, next, nedge;
  double **edge;
  int **cc2n;
  int checkf = 0;

  /* Map the locations into a continuous vector                      */
  /* Get the vector size                                             */
  /* Allocate                                                        */
  nedge = ns2;
  for (cc = 1; cc <= ns2; cc++) {
    nedge += mesh->npe[cc];
  }
  edge = d_alloc_2d(4, nedge);

  if (mesh->xloc && mesh->yloc) {
    d_free_1d(mesh->xloc);
    d_free_1d(mesh->yloc);
    mesh->xloc = d_alloc_1d(nedge + ns2 + 1);
    mesh->yloc = d_alloc_1d(nedge + ns2 + 1);
  }
  cc2n = i_alloc_2d(mesh->mnpe+1, ns2+1);

  /* Make the vector                                                 */
  next = 0;
  for (cc = 1; cc <= ns2; cc++) {
    for (j = 1; j <= mesh->npe[cc]; j++) {
      edge[next][0] = x[cc][j];
      edge[next][1] = y[cc][j];
      edge[next][2] = (double)cc;
      edge[next][3] = (double)j;
      if (cc == checkf) printf("old%d = [%f %f]\n",j, x[cc][j], y[cc][j]);
      next++;
    }
  }
  for (cc = 1; cc <= ns2; cc++) {
      edge[next][0] = x[cc][0];
      edge[next][1] = y[cc][0];
      edge[next][2] = (double)cc;
      edge[next][3] = 0.0;
      if (cc == checkf) printf("old%d = [%f %f]\n",j, x[cc][0], y[cc][0]);
      next++;
  }

  /* Sort edges, duplicates are consecutive!                         */
  qsort(edge[0], nedge, sizeof(double)*4, edge_sort_double);

  if (checkf) {
    for (n = 0; n < nedge; n++) {
      cc = (int)edge[n][2];
      j = (int)edge[n][3];
      if (cc == checkf) printf("new%d = %d[%f %f]\n", j, n, edge[n][0], edge[n][1]);
    }
  }

  /* Remove the duplicates                                           */
  mesh->ns = 1;
  for (n = 0; n < nedge; n++) {
    nn = (n == nedge - 1) ? n - 1 : n + 1;
    if ((edge[n][0] == edge[nn][0]) &&
	(edge[n][1] == edge[nn][1])) {
      cc = (int)edge[n][2];
      j = (int)edge[n][3];
      cc2n[cc][j] = mesh->ns;
      /* If the last entries are duplicates, make sure these are     */
      /* also included in the mesh array.                            */
      if (n == nedge - 1) cc2n[cc][j]--;
      if (nn == nedge - 1) {
	mesh->xloc[mesh->ns] = edge[n][0];
	mesh->yloc[mesh->ns] = edge[n][1];
	mesh->ns++;
      }
    } else {
      mesh->xloc[mesh->ns] = edge[n][0];
      mesh->yloc[mesh->ns] = edge[n][1];
      cc = (int)edge[n][2];
      j = (int)edge[n][3];
      cc2n[cc][j] = mesh->ns;
      mesh->ns++;
    }
  }
  mesh->ns--;

  /* Make the mapping from indices to coordinates                    */
  for (cc = 1; cc <= ns2; cc++) {
     for (j = 1; j <= mesh->npe[cc]; j++) {
       jj = (j == mesh->npe[cc]) ? 1 : j + 1;
       mesh->eloc[0][cc][j] = cc2n[cc][j];
       mesh->eloc[1][cc][j] = cc2n[cc][jj];
     }
     mesh->eloc[0][cc][0] = mesh->eloc[1][cc][0] = cc2n[cc][0];
  }

  if (checkf) {
    for (cc = 1; cc <= ns2; cc++) {
      for (j = 0; j <= mesh->npe[cc]; j++) {
	double xc = x[cc][j];
	double yc = y[cc][j];
	int found = 0;
	for (n = 1; n <= mesh->ns; n++) {
	  if (xc == mesh->xloc[n] && yc == mesh->yloc[n]) {
	    found = 1;
	    break;
	  }
	}
	if (!found) printf("Can't find coordinate cc=%d, j=%d [%f %f]\n", cc, j, xc, yc);
	if (mesh->eloc[0][cc][j] < 1 || mesh->eloc[0][cc][j] > mesh->ns) 
	  printf("Invalid0 index cc=%d, j=%d [%f %f] : %d\n", cc, j, xc, yc, mesh->eloc[0][cc][j]);
	if (mesh->eloc[1][cc][j] < 1 || mesh->eloc[1][cc][j] > mesh->ns) 
	  printf("Invalid1 index cc=%d, j=%d [%f %f] : %d\n", cc, j, xc, yc, mesh->eloc[1][cc][j]);
	if (cc == checkf) {
	  printf("check: cc=%d j=%d %f %f\n",cc, j, mesh->xloc[mesh->eloc[0][cc][j]], mesh->yloc[mesh->eloc[0][cc][j]]);
	}

      }
    }
  }
  d_free_2d(edge);
}

/* END remove_duplicates()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Expand the mesh around a nominated point to isolate and remove    */
/* any 'lakes' in the mesh. The mesh indices are re-ordered to       */
/* reflect the removal of the lakes.                                 */
/*-------------------------------------------------------------------*/
void mesh_expand(parameters_t *params, double *bathy, double **xc, double **yc)
{
  mesh_t *mesh = params->mesh;
  int cc, c, ci, cn, j, n, ns, ns2i;
  double dist, dmin;
  int found, ni, *film, *filla;
  double mlat, mlon, x, y;
  int verbose = 0;
  int **neic;
  int *n2o, *o2n;

  /*-----------------------------------------------------------------*/
  /* Expand inwards from open boundaries to remove 'lakes'           */
  if (params->mlon != NOTVALID && params->mlat != NOTVALID) {
    /* Open file                                                       */
    FILE *fp;
    if (params->mrf & MR_WRITEX) {
      if ((fp = fopen(params->mesh_reorder, "w")) == NULL) {
	hd_warn("Can't open mesh reorder write file %s\n", params->mesh_reorder);
	params->mrf &= ~MR_WRITEX;
      } else {
	fprintf(fp, "nMaxMesh2_face_nodes %d\n", mesh->mnpe);
	fprintf(fp, "nMesh2_face          %d\n", mesh->ns2);
	fprintf(fp, "nMesh2_face_indices  %d\n", mesh->ns);
      }
    }

    /* Get the index of the closest cell to the interior coordinate  */
    dmin = HUGE;
    ns2i = mesh->ns2;
    for (cc = 1; cc <= mesh->ns2; cc++) {
      x = mesh->xloc[mesh->eloc[0][cc][0]];
      y = mesh->yloc[mesh->eloc[0][cc][0]];
      dist = sqrt((params->mlon - x) * (params->mlon - x) +
		  (params->mlat - y) * (params->mlat - y));
      if (dist < dmin) {
	dmin = dist;
	ci = cc;
      }
    }

    /* Set the mask fanning out from the interior coordinate         */
    ns = mesh->ns2 + 2;
    film = i_alloc_1d(ns);
    memset(film, 0, ns * sizeof(int));
    filla = i_alloc_1d(ns);
    memset(filla, 0, ns * sizeof(int));
    found = 1;
    ni = 1;
    if (params->mrf & MR_WRITEX) fprintf(fp, "%d\n", ci);
    film[ni++] = ci;
    filla[ci] = 1;
    while (found) {
      found = 0;
      for (cc = 1; cc < ni; cc++) {
	c = film[cc];
	if (filla[c] == 2) continue;
	for (j = 1; j <= mesh->npe[c]; j++) {
	  cn = mesh->neic[j][c];
	  if (!filla[cn]) {
	    filla[cn] = 1;
	    filla[c] = 2;
	    film[ni] = cn;
	    if (cn && params->mrf & MR_WRITEX) fprintf(fp, "%d\n", cn);
	    found = 1;
	    ni++;
	  }
	}
      }
    }

    /* Make a copy of the current map                                */
    neic = i_alloc_2d(mesh->ns2+1, mesh->mnpe+1);
    n2o = i_alloc_1d(mesh->ns2+1);
    o2n = i_alloc_1d(mesh->ns2+1);
    for (cc = 1; cc <= mesh->ns2; cc++) {
      for (j = 1; j <= mesh->npe[cc]; j++) {
	neic[j][cc] = mesh->neic[j][cc];
      }
    }

    /* Remove cells which are part of 'lakes'                        */
    /* Note; the coordinates remain unchanged, with some redundant   */
    /* information corresponding to the lakes.                       */
    cn = 1;
    memset(film, 0, ns * sizeof(int));
    for (cc = 1; cc <= mesh->ns2; cc++) {
      if (filla[cc]) {
	mesh->npe[cn] = mesh->npe[cc];
	bathy[cn] = bathy[cc];
	for (j = 0; j <= mesh->npe[cc]; j++) {
	  xc[cn][j] = xc[cc][j];
	  yc[cn][j] = yc[cc][j];
	  mesh->eloc[0][cn][j] = mesh->eloc[0][cc][j];
	  mesh->eloc[1][cn][j] = mesh->eloc[1][cc][j];
	  mesh->neij[j][cn] = mesh->neij[j][cc];
	  n2o[cn] = cc;
	  o2n[cc] = cn;
	}
	film[cn++] = cc;
      } else {
	if (verbose) printf("%d %f %f\n",cc, mesh->xloc[mesh->eloc[0][cc][0]],
			    mesh->yloc[mesh->eloc[0][cc][0]]);
      }
    }
    mesh->ns2 = cn - 1;

    /* Remap the mapping function                                    */
    for (cc = 1; cc <= mesh->ns2; cc++) {
      for (j = 0; j <= mesh->npe[cc]; j++) {
	c = n2o[cc];
	mesh->neic[j][cc] = o2n[neic[j][c]];
      }
    }
    i_free_2d(neic);
    i_free_1d(n2o);
    i_free_1d(o2n);

    hd_warn("mesh_expand: %d cells eliminated\n", ns2i - mesh->ns2);

    /* OBCs cell centres are referenced to the indices and need to   */
    /* be remapped. The edge locations are referenced directly to    */
    /* the coordinates and do not.                                   */
    for (n = 0; n < mesh->nobc; n++) {
      for (cc = 1; cc <= mesh->npts[n]; cc++) {
	c = mesh->loc[n][cc];
	mesh->loc[n][cc] = film[c];
      }
    }
    i_free_1d(film);
    i_free_1d(filla);
    if (params->mrf & MR_WRITEX) fclose(fp);
  }
}

/* END mesh_expand()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Expand the mesh around a nominated point to create a vector (vec) */
/* that is continuously connected; i.e. neighbours in the vector are */
/* geographic neighbours.                                            */
/*-------------------------------------------------------------------*/
int mesh_expand_w(geometry_t *window,  /* Window geometry            */
		  int *vec)            /* Return vector              */
{
  int cc, c, c2, ci, cn, j;
  int found, ni, *film, *filla;
  int verbose = 0;

  /* Set the mask fanning out from the interior coordinate           */
  ci = 1;
  memset(vec, 0, window->szcS * sizeof(int));
  filla = i_alloc_1d(window->szcS);
  memset(filla, 0, window->szcS * sizeof(int));
  found = 1;
  ni = 1;
  vec[ni++] = ci;
  filla[ci] = 1;

  while (found) {
    found = 0;
    for (cc = 1; cc < ni; cc++) {
      c = vec[cc];
      c2 = window->m2d[c];
      if (filla[c] == 2) continue;
      for (j = 1; j <= window->npe[c2]; j++) {
	cn = window->c2c[j][c];
	if (window->wgst[cn]) filla[cn] = 2;
	if (!filla[cn]) {
	  filla[cn] = 1;
	  filla[c] = 2;
	  vec[ni] = cn;
	  found = 1;
	  ni++;
	}
      }
    }
  }
  ni--;
  return(ni);
}

/* END mesh_expand_w()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Expand the mesh around a nominated point to create a vector (vec) */
/* that is continuously connected; i.e. neighbours in the vector are */
/* geographic neighbours. The mesh is connected in directions of     */
/* outgoing flow.                                                    */
/*-------------------------------------------------------------------*/
int mesh_expand_3d(geometry_t *window,  /* Window geometry           */
		   double *u,           /* Velocity array            */
		   int *vecin           /* Cells to process          */
		   )            
{
  win_priv_t *wincon = window->wincon;
  window_t *windat = window->windat;
  int cc, c, c2, ci, cn, j, cm, n;
  int e, ee;
  int found, ni, *film, *filla;
  int dir;
  double d1;
  int *vec = wincon->s6;
  int *c2cc = wincon->s7;
  int *ctp = wincon->s1;
  int nctp = wincon->vc;

  memset(vec, 0, window->szc * sizeof(int));
  memset(c2cc, 0, window->szc * sizeof(int));
  filla = i_alloc_1d(window->szc);
  memset(filla, 0, window->szc * sizeof(int));
  ni = 1;
  /* If vec is provided, then this contains non-monotonic cells,     */
  /* and only expand out from these.                                 */               
  if (vecin != NULL) {
    ctp = vecin;
    nctp = vecin[0];
  }

  if (vecin == NULL) {
    for (cc = 1; cc <=wincon->vc; cc++) {
      c = wincon->s1[cc];
      c2 = window->m2d[c];
      n = 0;
      for (ee = 1; ee <= window->npe[c2]; ee++) {
	e = window->c2e[ee][c];
	/* Direction of flow                                         */
	dir = (window->eSc[ee][c2] * u[e] >= 0.0) ? 1 : -1;
	/* Sum the edges with outgoing flow                          */
	if (dir == 1) n++;
      }
      /* All edges are outflow                                       */
      if (!filla[c] && n == window->npe[c2]) {
	/*printf("start %d %d %d(%d %d)\n",master->nstep,ni,c,window->s2i[c],window->s2j[c]);*/
	vec[ni++] = c;
	filla[c] = 1;
	mesh_expand_do(window, u, vec, &ni, filla);
      }
    }
  }
  /* Collect unassigned cells having at least one outflow edge       */
  for (cc = 1; cc <= nctp; cc++) {
    c = ctp[cc];
    c2 = window->m2d[c];
    for (ee = 1; ee <= window->npe[c2]; ee++) {
      e = window->c2e[ee][c];
      dir = (window->eSc[ee][c2] * u[e] >= 0.0) ? 1 : -1;
      if (!filla[c] && dir == 1) {
	vec[ni++] = c;
	filla[c] = 1;
	mesh_expand_do(window, u, vec, &ni, filla);
      }
    }
  }
  ni--;
  vec[0] = ni;

  for (cc = 1; cc <=wincon->vcs; cc++) {
    c = wincon->s1[cc];
    c2 = window->m2d[c];
    c2cc[c2] = cc;
  }

  return(ni);
}

/* END mesh_expand_3d()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Performs mesh expansion                                           */
/*-------------------------------------------------------------------*/
int mesh_expand_do(geometry_t *window,  /* Window geometry           */
		   double *u,           /* Velocity array            */
		   int *vec,            /* Return vector             */
		   int *ni,
		   int *filla
		   )
{
  int cc, c, c2, j;
  int e, dir, cn;
  int found = 1;

  while (found) {
    found = 0;
    for (cc = 1; cc < *ni; cc++) {
      c = vec[cc];
      c2 = window->m2d[c];
      if (filla[c] == 2) continue;
      /*printf("do cc=%d ni=%d c=%d\n",cc,*ni,c);*/
      for (j = 1; j <= window->npe[c2]; j++) {
	e = window->c2e[j][c];
	dir = (window->eSc[j][c2] * u[e] >= 0.0) ? 1 : -1;
	cn = window->c2c[j][c];
	if (window->wgst[cn]) continue;
	if (dir == 1) {
	  if (!filla[cn]) {
	    vec[*ni] = cn;
	    /*printf("found %d %d %d(%d %d)\n",j,*ni,cn,window->s2i[cn],window->s2j[cn]);*/
	    *ni+=1;
	    filla[cn] = 1;
	  }
	  found = 1;
	}
      }
      filla[c] = 2;
    }
  }
}

/* END mesh_expand_do()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Removes nominated cells from the mesh and re-orders the mesh      */
/* structure.                                                        */
/*-------------------------------------------------------------------*/
void mesh_reduce(parameters_t *params, double *bathy, double **xc, double **yc)
{
  mesh_t *mesh = params->mesh;
  int cc, c, cn, cs, j, n, ns, ns2i;
  int ni, *mask;
  int verbose = 0;
  int **neic;
  int *n2o, *o2n;
  int dored = 0;

  /* Set the mask for cells to include in the mesh                   */
  ns = mesh->ns2 + 2;
  mask = i_alloc_1d(ns);
  memset(mask, 0, ns * sizeof(int));
  for (cc = 1; cc <= params->nland; cc++) {
    if (params->polyland != NULL && strlen(params->polyland[cc])) {
      FILE *fp;
      poly_t *pl = poly_create();
      if ((fp = fopen(params->polyland[cc], "r")) != NULL) {
	j = poly_read(pl, fp);
	for (c = 1; c <= mesh->ns2; c++) {
	  if (poly_contains_point(pl, mesh->xloc[mesh->eloc[0][c][0]], 
				  mesh->yloc[mesh->eloc[0][c][0]])) {
	    mask[c] = c;
	    mesh->npe[c] = 0;
	    dored = 1;
	  }
	}
	fclose(fp);
	poly_destroy(pl);
      }
    } else {
      c = params->lande1[cc];
      mask[c] = c;
      mesh->npe[c] = 0;
      dored = 1;
    }
  }
  if (!dored) return;

  /* Make a copy of the current map                                  */
  if (verbose) printf("Mesh reduction on %d cells\n", params->nland);
  ns2i = mesh->ns2;
  neic = i_alloc_2d(mesh->ns2+1, mesh->mnpe+1);
  n2o = i_alloc_1d(mesh->ns2+1);
  o2n = i_alloc_1d(mesh->ns2+1);
  for (cc = 1; cc <= mesh->ns2; cc++) {
    for (j = 1; j <= mesh->npe[cc]; j++) {
      neic[j][cc] = mesh->neic[j][cc];
    }
  }

  /* Note; the coordinates remain unchanged, with some redundant   */
  /* information.                                                  */
  cn = 1;
  mesh->map = i_alloc_1d(ns);
  memset(mesh->map, 0, ns * sizeof(int));
  for (cc = 1; cc <= mesh->ns2; cc++) {
    if (!mask[cc]) {
      mesh->npe[cn] = mesh->npe[cc];
      bathy[cn] = bathy[cc];
      for (j = 0; j <= mesh->npe[cc]; j++) {
	xc[cn][j] = xc[cc][j];
	yc[cn][j] = yc[cc][j];
	mesh->eloc[0][cn][j] = mesh->eloc[0][cc][j];
	mesh->eloc[1][cn][j] = mesh->eloc[1][cc][j];
	cs = mesh->neic[j][cc];
	if (mask[cs]) mesh->neij[j][cc] = j;
	mesh->neij[j][cn] = mesh->neij[j][cc];
	n2o[cn] = cc;
	o2n[cc] = cn;
      }
      mesh->map[cn++] = cc;
    } else {
      if (verbose) printf("%d %f %f\n",cc, mesh->xloc[mesh->eloc[0][cc][0]],
			  mesh->yloc[mesh->eloc[0][cc][0]]);
    }
  }
  mesh->ns2 = params->ns2 = cn - 1;

  /* Remap the mapping function                                      */
  for (cc = 1; cc <= mesh->ns2; cc++) {
    for (j = 0; j <= mesh->npe[cc]; j++) {
      c = n2o[cc];
      mesh->neic[j][cc] = o2n[neic[j][c]];
    }
  }
  i_free_2d(neic);
  i_free_1d(n2o);
  i_free_1d(o2n);

  hd_warn("mesh_reduce: %d cells eliminated\n", ns2i - mesh->ns2);

  /* OBCs cell centres are referenced to the indices and need to be  */
  /* remapped. The edge locations are referenced directly to the     */
  /* coordinates and do not.                                         */
  if (mesh->nobc) {
    for (n = 0; n < mesh->nobc; n++) {
      for (cc = 1; cc <= mesh->npts[n]; cc++) {
	c = mesh->loc[n][cc];
	mesh->loc[n][cc] = mesh->map[c];
      }
    }
    i_free_1d(mesh->map);
  }
}

/* END mesh_reduce()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Creates a triangulation using triangle using a bounding           */
/* perimeter as input. Output is stored in the delaunay structure.   */
/*-------------------------------------------------------------------*/
void xy_to_d(delaunay *d, int np, double *x, double *y)
{
  int i, j;
  point *pin;
  int ns, *sin;

  pin = malloc(np * sizeof(point));
  for(i = 0; i < np; i++) {
    point* p = &params->d->points[i];
    pin[i].x = x[i];
    pin[i].y = y[i];
  }
  j = 0;
  ns = 2 * np;
  sin = i_alloc_1d(2 * np);
  for(i = 0; i < np; i++) {
    sin[j++] = 2 * i;
    sin[j++] = 2 * i + 1;
  }
  sin[np - 1] = 0;
  d = delaunay_build(np, pin, ns, sin, 0, NULL);
}

/* END xy_to_d()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Creates a name for mesh output files                              */
/*-------------------------------------------------------------------*/
void mesh_ofile_name(parameters_t *params, char *key)
{
  char buf[MAXSTRLEN];
  int n, i, j;

  if (endswith(params->oname, ".nc")) {
    j = 3;
    strcpy(buf, params->oname);
  } else if (endswith(params->prmname, ".tran")) {
    j = 5;
    strcpy(buf, params->prmname);
  } else if (endswith(params->prmname, ".prm")) {
    j = 4;
    strcpy(buf, params->prmname);
  }
  n = strlen(buf);
  for (i = 0; i < n-j; i++)
    key[i] = buf[i];
  key[i] = '\0';
}
/* END mesh_ofile_name()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Writes a summary of output files from unstructured grid           */
/* generation.                                                       */
/*-------------------------------------------------------------------*/
void write_mesh_desc(parameters_t *params, coamsh_t *cm, FILE *fp, int mode)
{
  FILE *op, *ef;
  char buf[MAXSTRLEN], key[MAXSTRLEN], key1[MAXSTRLEN];
  int n, i;

  mesh_ofile_name(params, key);

  if (fp != NULL)
    op = fp;
  else {
    sprintf(buf,"%s_meshfiles.txt", key);
    op = fopen(buf, "w");
  }

  fprintf(op, "\n---------------------------------------------\n");
  fprintf(op, "Output diagnostic files from mesh generation.\n");
  fprintf(op, "---------------------------------------------\n\n");
  if (cm != NULL) {
    if (strlen(cm->pname)) {
      fprintf(op, "---------------------------------------------\n");
      fprintf(op, "coastmesh files.\n");
      sprintf(buf, "%s_in.txt", cm->pname);
      fprintf(op, "FILE %s contains the input perimeter for coastmesh\n", buf);
      fprintf(op, "and some output diagnostics.\n\n");
      sprintf(buf, "%s_out.txt", cm->pname);
      fprintf(op, "FILE %s contains the output perimeter for coastmesh.\n\n", buf);
      sprintf(buf, "%s.site", cm->pname);
      fprintf(op, "FILE %s contains labelled segments in a '.site' file.\n\n", buf);
    }
    if (strlen(cm->ofile)) {
      fprintf(op, "FILE %s contains the input geom_file .msh file for JIGSAW.\n", cm->ofile);
      fprintf(op, "This may be used to run JIGSAW offline.\n\n");
      for (i = 0; i < strlen(cm->ofile); i++) {
	if (cm->ofile[i] == '.')
	  break;
	else
	  buf[i] = cm->ofile[i];
      }
      buf[i] = '\0';
      sprintf(key1, "%s.jig", buf);
      fprintf(op, "FILE %s contains the input parameter .jig file for JIGSAW.\n", key1);
      fprintf(op, "This may be used to run JIGSAW offline.\n\n");
    }
  }
  if (mode & M_BOUND_HFUN) {
    fprintf(op, "---------------------------------------------\n");
    fprintf(op, "Files produced in creating the hfun function.\n");
    sprintf(buf,"%s_e.txt", key);
    fprintf(op, "FILE %s contains edges of the hfun file. This is likely\n", buf);
    fprintf(op, "overwritten at a later date, unless the model fails\n");
    fprintf(op, "during the construction of the JIGSAW mesh. May be used for\n");
    fprintf(op, "plotting; e.g. 'plot(in_e(:,1),in_e(:,2),'g-')'.\n\n");

    sprintf(buf,"%s_h.txt", key);
    fprintf(op, "FILE %s contains the coordinates and value of the hfun function.\n", buf);
    fprintf(op, "The bathymetry used in the conversion is also included. May be used\n");
    fprintf(op, "for plotting; e.g. 'scatter(in_h(:,1),in_h(:,2),100,in_h(:,3),'filled')'.\n\n");

    sprintf(buf,"%s_hfun.msh", key);
    fprintf(op, "FILE %s contains the values of the hfun function in .msh format.\n", buf);
    fprintf(op, "This may be used when running JIGSAW offline.\n\n");
  }

  if (mode & M_CREATE_JIG) {
    fprintf(op, "---------------------------------------------\n");
    fprintf(op, "Files produced by invoking JIGSAW.\n");
    sprintf(buf,"%s_t.txt", key);
    fprintf(op, "FILE %s contains the triangulation of the JIGSAW output.\n", buf);
    fprintf(op, "May be used for plotting.\n\n");
  }

  if (mode & M_CONVERT_JIG) {
    fprintf(op, "---------------------------------------------\n");
    fprintf(op, "Files produced converting JIGSAW input to a Voronoi mesh.\n");
    sprintf(buf,"%s_jc.txt", key);
    fprintf(op, "FILE %s contains the vertices of the JIGSAW triangulation.\n", buf);
    sprintf(buf,"%s_jce.txt", key);
    fprintf(op, "FILE %s contains the centre of mass of the JIGSAW triangulation.\n", buf);
    fprintf(op, "May be used for plotting.\n\n");
    sprintf(buf,"%s_je.txt", key);
    fprintf(op, "FILE %s contains the perimeter edges used in JIGSAW.\n", buf);
    fprintf(op, "May be used for plotting.\n\n");
    sprintf(buf,"%s_jt.txt", key);
    fprintf(op, "FILE %s contains the Delaunay edges of the JIGSAW triangulation.\n", buf);
    fprintf(op, "May be used for plotting.\n\n");
    sprintf(buf,"%s_jv.txt", key);
    fprintf(op, "FILE %s contains the Voronoi edges.\n", buf);
    fprintf(op, "May be used for plotting.\n\n");
    sprintf(buf,"%s_jig.msh", key);
    fprintf(op, "FILE %s may be produced containing a .msh format of the JIGSAW output.\n\n", buf);
  }

  if (mode & M_MESHSTRUCT) {
    mesh_t *mesh = params->mesh;
    fprintf(op, "---------------------------------------------\n");
    fprintf(op, "Files produced converting to COMPAS mesh specification.\n");
    sprintf(buf,"%s_e.txt", key);
    fprintf(op, "FILE %s contains the Voronoi edges.\n", buf);
    fprintf(op, "May be used for plotting.\n\n");
    sprintf(buf,"%s_c.txt", key);
    fprintf(op, "FILE %s contains the Voronoi centres.\n", buf);
    fprintf(op, "May be used for plotting (e.g. using Matlab script 'plot_%s.m').\n\n", key);
    if (params->nobc) {
      sprintf(buf,"%s_b.txt", key);
      fprintf(op, "FILE %s contains the open boundary perimeter.\n", buf);
      fprintf(op, "May be used for plotting.\n\n");
    }
    sprintf(buf,"%s_us.txt", key);
    fprintf(op, "FILE %s contains the COMPAS specification of the mesh.\n", buf);
    fprintf(op, "May be used to define a grid in the input file.\n\n");

    sprintf(buf,"plot_%s.m", key);
    if ((ef = fopen(buf, "w")) != NULL) {
      fprintf(ef,"%% Matlab script to plot an unstructured mesh using files\n");
      fprintf(ef,"%% %s_e.txt and %s_c.txt output.\n\n", key, key);
      fprintf(ef,"close all;\n");
      fprintf(ef,"load('%s_e.txt');\n", key);
      fprintf(ef,"plot(%s_e(:,1), %s_e(:,2), 'b-');\n", key, key);
      fprintf(ef,"hold on;\n");
      fprintf(ef,"load('%s_c.txt');\n", key);
      fprintf(ef,"plot(%s_c(:,1), %s_c(:,2), 'r.');\n", key, key);
      if (params->nobc) {
	fprintf(ef,"load('boundary.txt');\n");
	fprintf(ef,"plot(boundary(:,1), boundary(:,2), 'g-');\n");
      }
      fprintf(ef,"\nprint -djpeg90 o1.jpg\n");
      fclose(ef);
    }
  }
  fclose(op);
}

/* END write_mesh_desc()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Makes a Delaunay dual from a Voronoi mesh, layer k                */
/*-------------------------------------------------------------------*/
delaunay*  make_dual(geometry_t *geom, int kin)
{
  delaunay* d = delaunay_create();
  int c, c2, ci, cn, cc, ee, e;
  int i, j, jj, k, n;
  int **nei;
  int *mask;
  int verbose = 1;

  mask = i_alloc_1d(geom->szc);
  memset(mask, 0, geom->szc * sizeof(int));
  geom->c2p = i_alloc_2d(geom->szcS, geom->nz+1);

  /*-----------------------------------------------------------------*/
  /* Get the points                                                  */
  d->npoints = 0;
  /* Count the points in layer k                                     */
  for (cc = 1; cc <= geom->n3_t; cc++) {
    c = geom->w3_t[cc];
    k = geom->s2k[c];
    if (k == kin) {
      d->npoints++;
    }
  }

  n = 0;
  d->points = malloc(d->npoints * sizeof(point));
  /* Wet cells                                                       */
  for (cc = 1; cc <= geom->b3_t; cc++) {
    c = geom->w3_t[cc];
    c2 = geom->m2d[c];
    k = geom->s2k[c];
    if (k == kin) {
      geom->c2p[k+1][c2] = n;
      d->points[n].x = geom->cellx[c2];
      d->points[n++].y = geom->celly[c2];
    }
  }
  /* Ghost cells                                                     */
  for (cc = geom->b3_t+1; cc <= geom->n3_t; cc++) {
    c = geom->w3_t[cc];
    ci = geom->wgst[c];
    c2 = geom->m2d[c];
    k = geom->s2k[ci];
    if (k == kin) {
      if (!mask[ci]) {
	geom->c2p[k][geom->m2d[c]] = n;
	d->points[n].x = geom->cellx[c2];
	d->points[n++].y = geom->celly[c2];
	mask[ci] = 1;
      }
      /*
      for (n = 1; n <= geom->npe[c2]; n++) {
	if (c == geom->c2c[n][ci])
	  break;
      }
      e = geom->m2de[geom->c2e[n][ci]];
      geom->c2p[k][geom->m2d[c]] = n;
      d->points[n]->x = geom->u1x[e];
      d->points[n++]->y = geom->u1y[e];
      */
    }
  }
  /* Open boundary ghosts                                            */
  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    for (ee = 1; ee <= open->no3_e1; ee++) {
      c = open->ogc_t[ee];
      ci = open->obc_e2[ee];
      c2 = geom->m2d[c];
      k = geom->s2k[ci];
      if (k == kin && geom->c2p[k][c2] == -1) {
	geom->c2p[k][c2] = n;
	d->points[n].x = geom->cellx[c2];
	d->points[n++].y = geom->celly[c2];
      }
    }
  }
  if (verbose) printf("npoints = %d\n", d->npoints);

  /*-----------------------------------------------------------------*/
  /* Get the edges                                                   */
  d->nedges = 0;
  memset(mask, 0, geom->szc * sizeof(int));
  for (cc = 1; cc <= geom->b3_t; cc++) {
    c = geom->w3_t[cc];
    c2 = geom->m2d[c];
    k = geom->s2k[c];
    if (k == kin) {
      for (j = 1; j <= geom->npe[c2]; j++) {
	ci = geom->c2c[j][c];
	if (!mask[c] && !mask[ci]) {
	  d->nedges++;
	  mask[c] = 1;
	  mask[ci] = 1;
	}
      }
    }
  }
  memset(mask, 0, geom->szc * sizeof(int));
  d->edges = i_alloc_1d(2 * d->nedges);
  n = 0;
  for (cc = 1; cc <= geom->b3_t; cc++) {
    c = geom->w3_t[cc];
    c2 = geom->m2d[c];
    k = geom->s2k[c];
    if (k == kin) {
      for (j = 1; j <= geom->npe[c2]; j++) {
	ci = geom->c2c[j][c];
	if (!mask[c] && !mask[ci]) {
	  d->edges[n++] = geom->c2p[k+1][c2];
	  d->edges[n++] = geom->c2p[k+1][geom->m2d[ci]];
	  mask[c] = 1;
	  mask[ci] = 1;
	}
      }
    }
  }
  if (verbose) printf("nedges = %d\n", d->nedges);

  /*-----------------------------------------------------------------*/
  /* Get the triangles. First count the triangles.                   */
  d->ntriangles = 0;
  memset(mask, 0, geom->szc * sizeof(int));
  for (cc = 1; cc <= geom->b3_t; cc++) {
    c = geom->w3_t[cc];
    c2 = geom->m2d[c];
    k = geom->s2k[c];
    if (k == kin) {
      for (j = 1; j <= geom->npe[c2]; j++) {
	jj = (j == geom->npe[c2]) ? 1 : j + 1;
	ci = geom->c2c[j][c];
	cn = geom->c2c[jj][c];
	if (!mask[c] || !mask[ci] || !mask[cn]) {
	  d->ntriangles++;
	  mask[c] = 1;
	  mask[ci] = 1;
	  mask[cn] = 1;
	  printf("%d : %d %d %d\n",d->ntriangles,c,ci,cn);
	}
      }
    }
  }
  if (verbose) printf("ntriangles = %d\n", d->ntriangles);
  /* Next set the triangles                                          */
  d->triangles = malloc(d->ntriangles * sizeof(triangle));
  n = 0;
  memset(mask, 0, geom->szc * sizeof(int));
  for (cc = 1; cc <= geom->b3_t; cc++) {
    c = geom->w3_t[cc];
    c2 = geom->m2d[c];
    k = geom->s2k[c];
    if (k == kin) {
      for (j = 1; j <= geom->npe[c2]; j++) {
	jj = (j == geom->npe[c2]) ? 1 : j + 1;
	ci = geom->c2c[j][c];
	cn = geom->c2c[jj][c];
	if (!mask[c] || !mask[ci] || !mask[cn]) {
	  triangle* t = &d->triangles[n++];
	  t->vids[0] = geom->c2p[k+1][c2];
	  t->vids[1] = geom->c2p[k+1][geom->m2d[ci]];
	  t->vids[2] = geom->c2p[k+1][geom->m2d[cn]];
	  mask[c] = 1;
	  mask[ci] = 1;
	  mask[cn] = 1;
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the neighbours                                              */
  neighbour_finder_b(d, &nei);
  d->neighbours = malloc(d->ntriangles * sizeof(triangle_neighbours));
  for (i = 0; i < d->ntriangles; ++i) {
    triangle_neighbours* ne = &d->neighbours[i];
    ne->tids[0] = nei[0][i];
    ne->tids[1] = nei[1][i];
    ne->tids[2] = nei[2][i];
  }
  i_free_2d(nei);
  i_free_1d(mask);
  return d;
}

/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Merges a quadrilateral grid with an unstructured mesh at a user   */
/* defined location.                                                 */
/*-------------------------------------------------------------------*/
void add_quad_grid(parameters_t *params, char *iname, int mode)
{
  FILE *fp, *dp, *op, *pp = NULL;
  mesh_t *mesh;
  char buf[MAXSTRLEN];
  char polyr[MAXSTRLEN];
  int i, j, n, nc, ns2, npe, ncells, nce1, nce2, ice1, ice2;
  int si, ei, sj, ej, prmf = 0;
  double *x, *y;
  double **lat, **lon;
  double **xc, **yc;
  double **xm, **ym;
  double *bathy, *bathyin = NULL, nbathy = NOTVALID, *b1 = NULL;
  double *x1 = NULL, *y1 = NULL;
  int *mpe;
  double xr1, yr1, xr2, yr2, d, dist1, dist2;
  int *npe2, ir, jr1, jr2, *mask;
  int verbose = 0;
  int dirf = 1;    /* Coordinates are ordered across the river       */
  int sortdir = 1; /* Sort verices; 1=clockwise, -1=anticlockwise    */

  if (mode) {
    ns2 = params->ns2;
    npe = params->npe;
    mpe = params->npe2;
  } else {
    mesh = params->mesh;
    ns2 = mesh->ns2;
    npe = mesh->mnpe;
    mpe = mesh->npe;
  }

  /*-----------------------------------------------------------------*/
  /* Open the file to merge and read the merge location              */
  if ((fp = fopen(iname, "r")) == NULL) {
    hd_warn("Can't open merge grid '%s'\n", iname);
    return;
  }

  if (prm_read_char(fp, "MERGE_LOC", buf))
    sscanf(buf, "(%lf %lf)-(%lf %lf)", &xr1, &yr1, &xr2, &yr2);
  else {
    hd_warn("Can't find MERGE_LOC for %s: no mesh merging.\n", iname);
    return;
  }

  /* Dimensions                                                      */
  prm_read_int(fp, "NCE1", &nce1);
  prm_read_int(fp, "NCE2", &nce2);
  if (nce1 != 1 && nce2 != 1) hd_warn("One of NCE1 or NCE2 must = 1: no mesh merging.\n");

  /* Cookie cut from structured .prm file                            */
  if (prm_read_char(fp, "PRM_FILE", buf)) {
    if ((pp = fopen(buf, "r")) == NULL)
      hd_quit("add_quad_grid: Can't open prm file %s\n", buf);
    prmf = 1;
    prm_read_int(pp, "NCE1", &ice1);
    prm_read_int(pp, "NCE2", &ice2);
    if (prm_read_char(fp, "SUBSECTION", buf))
    sscanf(buf, "(%d,%d)-(%d,%d)", &si, &sj, &ei, &ej);
    if (si > ei) {
      i = ei;
      ei = si;
      si = i;
    }
    if (sj > ej) {
      j = ej;
      ej = sj;
      sj = j;
    }
    nce1 = ei - si + 1;
    nce2 = ej - sj + 1;
    if (!prm_read_double(fp, "BATHYVAL", &nbathy)) {
      prm_read_darray(pp, "BATHY", &b1, &i);
      bathyin = d_alloc_1d(nce1 * nce2 + 1);
      nc = 0;
      for (j = sj; j <= ej; j++)
	for (i = si; i <= ei; i++) {
	  n = j * ice1 + (i);
	  bathyin[nc++] = b1[n];
	}
      d_free_1d(b1);
    }
  }
  ncells = nce1 * nce2;

  /* Bathymetry                                                      */
  if (prm_read_darray(fp, "BATHY", &bathyin, &i)) {
    if (i != nce1 && i != nce2)
      hd_quit("add_quad_grid: number of bathy values != NCE1 or NCE2 (%d!=%d,%d)\n",
	      i, nce1, nce2);
  } else {
    prm_read_double(fp, "BATHYVAL", &nbathy);
  }

  if (prm_read_char(fp, "MERGE_DIR", buf))
    dirf = is_true(buf);
  sprintf(polyr, "%c", '\0');
  prm_read_char(fp, "REMOVE", polyr);

  /* The merge location is set to the closest vertices in the mesh.  */
  /* Also make a copy of the current mesh.                           */
  dist1 = dist2 = HUGE;
  npe2 = i_alloc_1d(ns2+1);
  bathy = d_alloc_1d(ns2+1);
  xc = d_alloc_2d(npe+1, ns2+1);
  yc = d_alloc_2d(npe+1, ns2+1);
  mask = i_alloc_1d(ns2+1);
  memset(mask, 0, (ns2+1) * sizeof(int));
  for (i = 1; i <= ns2; i++) {
    npe2[i] = mpe[i];
    if (mode) bathy[i] = params->bathy[i];
    for (j = 0; j <= mpe[i]; j++) {
      if (mode) {
	xc[i][j] = params->x[i][j];
	yc[i][j] = params->y[i][j];
      } else {
	xc[i][j] = mesh->xloc[mesh->eloc[0][i][j]];
	yc[i][j] = mesh->yloc[mesh->eloc[0][i][j]];
      }
      d = sqrt((xc[i][j] - xr1) * (xc[i][j] - xr1) + (yc[i][j] - yr1) * (yc[i][j] - yr1));
      if (j > 0 && d < dist1) {
	dist1 = d;
	ir = i;
	jr1 = j;
      }
    }
  }
  xr1 = xc[ir][jr1];
  yr1 = yc[ir][jr1];
  for (i = 1; i <= ns2; i++) {
    for (j = 1; j <= mpe[i]; j++) {
      d = sqrt((xc[i][j] - xr2) * (xc[i][j] - xr2) + (yc[i][j] - yr2) * (yc[i][j] - yr2));
      if (j > 0 && d < dist2) {
	dist2 = d;
	ir = i;
	jr2 = j;
      }
    }
  }
  xr2 = xc[ir][jr2];
  yr2 = yc[ir][jr2];
  if (verbose) {
    printf("%f %f s r\n", xr1, yr1);
    printf("%f %f e r\n", xr2, yr2);
  }
  if (nbathy == NOTVALID) nbathy = bathy[ir];
  hd_warn("Adding QUAD grid %s into mesh at locations (%f %f)-(%f %f)\n", iname, xr1, yr1, xr2, yr2);
 
  /*-----------------------------------------------------------------*/
  /* Read and save the vertices and centre of the quadrilateral grid */
  if (prmf) {
    int mm = 2 * ice1 + 1;
    prm_read_darray(pp, "XCOORDS", &x1, &nc);
    prm_read_darray(pp, "YCOORDS", &y1, &nc);
    x = d_alloc_1d(nc);
    y = d_alloc_1d(nc);
    nc = 0;
    for (j = 2*sj; j <= 2*ej+2; j++)
      for (i = 2*si+1; i <= 2*ei+1; i++) {
	n = (j) * mm + i;
        x[nc] = x1[n-1];
        y[nc++] = y1[n-1];
      }
    d_free_1d(x1);
    d_free_1d(y1);
    fclose(pp);
  } else {
    prm_read_int(fp, "XCOORDS", &nc);
    if (nc != (2 * nce1 + 1) * (2 * nce2 + 1)) hd_quit("add_quad_grid: %s wrong grid size: (%d,%d) != %d\n", iname, i, j, nc);
    x = d_alloc_1d(nc);
    for (i = 0; i < nc; i++) fscanf(fp, "%lf", &x[i]);
    prm_read_int(fp, "YCOORDS", &nc);
    y = d_alloc_1d(nc);
    for (i = 0; i < nc; i++) fscanf(fp, "%lf", &y[i]);
  }
  lon = d_alloc_2d(5, ncells);
  lat = d_alloc_2d(5, ncells);
  n = 0;
  dist1 = dist2 = HUGE;
  if (verbose) printf("file=%s ncells=%d nc=%d dir=%d\n", iname, ncells, nc, dirf);

  /* If dirf = 1 then the orthogonal grid coordinates are ordered    */
  /* with values increasing sequentially across the channel. If      */
  /* dirf = 0 then the coordinates are ordered such that values      */
  /* sequentially increase along the channel. The user must discern  */
  /* which is the case and set MERGE_DIR accordingly.                */
  if (dirf) {
    for (i = 0; i < ncells; i++) {
      n = i * 6;
      lon[i][0] = x[n+4];
      lat[i][0] = y[n+4];
      lon[i][1] = x[n];
      lat[i][1] = y[n];
      lon[i][2] = x[n+2];
      lat[i][2] = y[n+2];
      lon[i][3] = x[n+6];
      lat[i][3] = y[n+6];
      lon[i][4] = x[n+8];
      lat[i][4] = y[n+8];
      if(verbose && i==0)for(n=0; n<=4; n++)printf("%f %f y%d\n",lon[i][n],lat[i][n],n);
    }
  } else {
    for (i = 0; i < ncells; i++) {
      n = i * 2;
      lon[i][0] = x[nc/3+n+1];
      lat[i][0] = y[nc/3+n+1];
      lon[i][1] = x[n];
      lat[i][1] = y[n];
      lon[i][2] = x[n+2];
      lat[i][2] = y[n+2];
      lon[i][3] = x[2*nc/3+n];
      lat[i][3] = y[2*nc/3+n];
      lon[i][4] = x[2*nc/3+n+2];
      lat[i][4] = y[2*nc/3+n+2];
      if(verbose && i==0)for(n=0; n<=4; n++)printf("%f %f a%d\n",lon[i][n],lat[i][n],n);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the merge location to the closest vertices in the grid      */
  if (verbose) printf("Old start merge: i=%d j=%d [%f %f]\n", ir, jr1, xr1, yr1);
  for (i = 0; i < ncells; i++) {
    sort_circle_g(lon[i], lat[i], 4, sortdir);
    for (j = 1; j <= 4; j++) {
      d = sqrt((lon[i][j] - xr1) * (lon[i][j] - xr1) + (lat[i][j] - yr1) * (lat[i][j] - yr1));
      if (d < dist1) {
	dist1 = d;
	ir = i;
	jr1 = j;
      }
    }
  }
  lon[ir][jr1] = xr1;
  lat[ir][jr1] = yr1;
  if (verbose) {
    printf("New start merge: i=%d j=%d\n", ir, jr1);
    printf("Old end merge: j=%d [%f %f]\n", jr2, xr2, yr2);
  }
  for (j = 0; j <= 4; j++) {
    d = sqrt((lon[ir][j] - xr2) * (lon[ir][j] - xr2) + (lat[ir][j] - yr2) * (lat[ir][j] - yr2));
    if (j > 0 && j != jr1 && d < dist2) {
      dist2 = d;
      jr2 = j;
    }
  }
  lon[ir][jr2] = xr2;
  lat[ir][jr2] = yr2;
  if(verbose) {
    printf("New end merge: j=%d\n", jr2);
    for(n=0; n<=4; n++)printf("%f %f b%d\n",lon[0][n],lat[0][n],n);
  }
  hd_warn("Curvilinear grid %s added at (%f %f)-(%f %f)\n", iname, xr1, yr1, xr2, yr2);

  /*-----------------------------------------------------------------*/
  /* Merge a second set of coordinates if required (e.g. for two     */
  /* ends of a channel.                                              */
  if (prm_read_char(fp, "MERGE_LOC2", buf)) {
    sscanf(buf, "(%lf %lf)-(%lf %lf)", &xr1, &yr1, &xr2, &yr2);

    /* The merge location is set to the closest vertices in the      */
    /* mesh.                                                         */
    dist1 = dist2 = HUGE;
    for (i = 1; i <= ns2; i++) {
      for (j = 0; j <= mpe[i]; j++) {
	d = sqrt((xc[i][j] - xr1) * (xc[i][j] - xr1) + (yc[i][j] - yr1) * (yc[i][j] - yr1));
	if (j > 0 && d < dist1) {
	  dist1 = d;
	  ir = i;
	  jr1 = j;
	}
      }
    }
    xr1 = xc[ir][jr1];
    yr1 = yc[ir][jr1];
    for (i = 1; i <= ns2; i++) {
      for (j = 0; j <= mpe[i]; j++) {
	d = sqrt((xc[i][j] - xr2) * (xc[i][j] - xr2) + (yc[i][j] - yr2) * (yc[i][j] - yr2));
	if (j > 0 && d < dist2) {
	dist2 = d;
	ir = i;
	jr2 = j;
	}
      }
    }
    xr2 = xc[ir][jr2];
    yr2 = yc[ir][jr2];
    hd_warn("  Second coordinates to merge are (%f %f)-(%f %f)\n", xr1, yr1, xr2, yr2);

    /* Set the merge location to the closest vertices in the grid    */
    dist1 = dist2 = HUGE;
    for (i = 0; i < ncells; i++) {
      /*sort_circle_g(lon[i], lat[i], 4, sortdir);*/
      for (j = 1; j <= 4; j++) {
	d = sqrt((lon[i][j] - xr1) * (lon[i][j] - xr1) + (lat[i][j] - yr1) * (lat[i][j] - yr1));
	if (d < dist1) {
	  dist1 = d;
	  ir = i;
	  jr1 = j;
	}
      }
    }
    lon[ir][jr1] = xr1;
    lat[ir][jr1] = yr1;
    for (j = 0; j <= 4; j++) {
      d = sqrt((lon[ir][j] - xr2) * (lon[ir][j] - xr2) + (lat[ir][j] - yr2) * (lat[ir][j] - yr2));
      if (j > 0 && j != jr1 && d < dist2) {
	dist2 = d;
	jr2 = j;
      }
    }
    lon[ir][jr2] = xr2;
    lat[ir][jr2] = yr2;
    if(verbose)for(n=0; n<=4; n++)printf("%f %f c%d\n",lon[0][n],lat[0][n],n);
    hd_warn("Curvilinear grid %s added at second location (%f %f)-(%f %f)\n", iname, xr1, yr1, xr2, yr2);
  }

  if (mode) {
    ns2 = params->ns2;
    /* Sometimes a few cells from the mesh should be removed so that */
    /* the quad grid can seamlessly interface with the mesh.         */
    /* Read polygons to remove some cells if required                */
    if (strlen(polyr)) {
      char *fields[MAXSTRLEN * MAXNUMARGS];
      int npoly = parseline(polyr, fields, MAXNUMARGS);
      FILE *pf;
      poly_t *pl;
      ns2 = params->ns2;
      for (n = 0; n < npoly; n++) {
	if ((pf = fopen(fields[n], "r")) == NULL) {
	  hd_warn("Can't find file %s to remove from quad addition %s\n", 
		  fields[0], iname);
	  continue;
	}
	pl = poly_create();
	j = poly_read(pl, pf);
	for (i = 1; i <= ns2; i++) {
	  if (poly_contains_point(pl, params->x[i][0], params->y[i][0])) {
	    mask[i] = i;
	    params->ns2--;
	  }
	}
	fclose(pf);
	poly_destroy(pl);
      }
    }
    
    /* Reallocate and reset the (x,y) arrays, including the points in  */
    /* the grid.                                                       */
    params->ns2 += ncells;
    i_free_1d(params->npe2);
    d_free_1d(params->bathy);
    d_free_2d(params->x);
    d_free_2d(params->y);
    params->npe2 = i_alloc_1d(params->ns2+1);
    params->bathy = d_alloc_1d(params->ns2+1);
    params->x = d_alloc_2d(params->npe+1, params->ns2+1);
    params->y = d_alloc_2d(params->npe+1, params->ns2+1);
    n = 1;
    for (i = 1; i <= ns2; i++) {
      if (mask[i]) continue;
      params->npe2[n] = npe2[i];
      params->bathy[n] = bathy[i];
      for (j = 0; j <= params->npe2[n]; j++) {
	params->x[n][j] = xc[i][j];
	params->y[n][j] = yc[i][j];
      }
      n++;
    }
    for (i = 0; i < ncells; i++) {
      params->npe2[n] = 4;
      if (bathyin)
	params->bathy[n] = -fabs(bathyin[i]);
      else
	params->bathy[n] = -fabs(nbathy);
      for (j = 0; j <= params->npe2[n]; j++) {
	params->x[n][j] = lon[i][j];
	params->y[n][j] = lat[i][j];
      }
      n++;
    }
  } else {
    /* Create ascii files with the added quad grid locations         */
    dp = fopen("quad_c.xy", "w");
    op = fopen("quad.site", "w");
    for (i = 0; i < ncells; i++) {
      fprintf(dp, "%f %f\n",lon[i][0], lat[i][0]);
      fprintf(op, "%f %f c%d r\n",lon[i][0], lat[i][0], i);
    }
    fclose(dp);
    dp = fopen("quad_e.xy", "w");
    n = 0;
    for (i = 0; i < ncells; i++) {
      for (j = 1; j <= 4; j++) {
	fprintf(dp, "%f %f\n",lon[i][j], lat[i][j]);
	fprintf(op, "%f %f e%d\n",lon[i][j], lat[i][j], n);
	n++;
      } 
      fprintf(dp, "%f %f\n",lon[i][1], lat[i][1]);
      fprintf(dp, "NaN NaN\n");
    }
    fclose(dp);
    fclose(op);
  }

  fclose(fp);
  d_free_1d(x);
  d_free_1d(y);
  i_free_1d(npe2);
  d_free_2d(lon);
  d_free_2d(lat);
  d_free_2d(xc);
  d_free_2d(yc);
  i_free_1d(mask);
  if (bathyin) d_free_1d(bathyin);
}

/* END add_quad_grid()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Merges a quadrilateral grid with an unstructured mesh at a user   */
/* defined location.                                                 */
/*-------------------------------------------------------------------*/
void add_quad_grido(parameters_t *params, char *iname)
{
  FILE *fp;
  char buf[MAXSTRLEN];
  char polyr[MAXSTRLEN];
  int i, j, n, nc, ns2, ncells, nce1, nce2;
  double *x, *y;
  double **lat, **lon;
  double **xc, **yc;
  double *bathy, nbathy = NOTVALID;
  double xr1, yr1, xr2, yr2, d, dist1, dist2;
  int *npe2, ir, jr1, jr2, *mask;
  int verbose = 0;
  int dirf = 1;    /* Coordinates are ordered across the river       */
  int sortdir = 1; /* Sort verices; 1=clockwise, -1=anticlockwise    */

  /*-----------------------------------------------------------------*/
  /* Open the file to merge and read the merge location              */
  if ((fp = fopen(iname, "r")) == NULL) {
    hd_warn("Can't open merge grid '%s'\n", iname);
    return;
  }
  if (prm_read_char(fp, "MERGE_LOC", buf))
    sscanf(buf, "(%lf %lf)-(%lf %lf)", &xr1, &yr1, &xr2, &yr2);
  else {
    hd_warn("Can't find MERGE_LOC for %s: no mesh merging.\n", iname);
    return;
  }
  prm_read_double(fp, "BATHYVAL", &nbathy);
  if (prm_read_char(fp, "MERGE_DIR", buf))
    dirf = is_true(buf);
  sprintf(polyr, "%c", '\0');
  prm_read_char(fp, "REMOVE", polyr);

  /* The merge location is set to the closest vertices in the mesh.  */
  /* Also make a copy of the current mesh.                           */
  dist1 = dist2 = HUGE;
  npe2 = i_alloc_1d(params->ns2+1);
  bathy = d_alloc_1d(params->ns2+1);
  xc = d_alloc_2d(params->npe+1, params->ns2+1);
  yc = d_alloc_2d(params->npe+1, params->ns2+1);
  mask = i_alloc_1d(params->ns2+1);
  memset(mask, 0, (params->ns2+1) * sizeof(int));
  for (i = 1; i <= params->ns2; i++) {
    npe2[i] = params->npe2[i];
    bathy[i] = params->bathy[i];
    for (j = 0; j <= params->npe2[i]; j++) {
      xc[i][j] = params->x[i][j];
      yc[i][j] = params->y[i][j];
      d = sqrt((xc[i][j] - xr1) * (xc[i][j] - xr1) + (yc[i][j] - yr1) * (yc[i][j] - yr1));
      if (j > 0 && d < dist1) {
	dist1 = d;
	ir = i;
	jr1 = j;
      }
    }
  }
  xr1 = xc[ir][jr1];
  yr1 = yc[ir][jr1];
  for (i = 1; i <= params->ns2; i++) {
    for (j = 1; j <= params->npe2[i]; j++) {
      d = sqrt((xc[i][j] - xr2) * (xc[i][j] - xr2) + (yc[i][j] - yr2) * (yc[i][j] - yr2));
      if (j > 0 && d < dist2) {
	dist2 = d;
	ir = i;
	jr2 = j;
      }
    }
  }
  xr2 = xc[ir][jr2];
  yr2 = yc[ir][jr2];
  if (nbathy == NOTVALID) nbathy = bathy[ir];
  hd_warn("Adding QUAD grid %s into mesh at locations (%f %f)-(%f %f)\n", iname, xr1, yr1, xr2, yr2);
 
  /*-----------------------------------------------------------------*/
  /* Read and save the vertices and centre of the quadrilateral grid */
  prm_read_int(fp, "NCE1", &nce1);
  prm_read_int(fp, "NCE2", &nce2);
  if (nce1 != 1 && nce2 != 1) hd_warn("One of NCE1 or NCE2 must = 1: no mesh merging.\n");
  ncells = nce1 * nce2;
  prm_read_int(fp, "XCOORDS", &nc);
  if (nc != (2 * nce1 + 1) * (2 * nce2 + 1)) hd_quit("add_quad_grid: %s wrong grid size: (%d,%d) != %d\n", iname, i, j, nc);
  x = d_alloc_1d(nc);
  for (i = 0; i < nc; i++) fscanf(fp, "%lf", &x[i]);
  prm_read_int(fp, "YCOORDS", &nc);
  y = d_alloc_1d(nc);
  for (i = 0; i < nc; i++) fscanf(fp, "%lf", &y[i]);
  lon = d_alloc_2d(5, ncells);
  lat = d_alloc_2d(5, ncells);
  n = 0;
  dist1 = dist2 = HUGE;
  if (verbose) printf("file=%s ncells=%d nc=%d dir=%d\n", iname, ncells, nc, dirf);

  /* If dirf = 1 then the orthogonal grid coordinates are ordered    */
  /* with values increasing sequentially across the channel. If      */
  /* dirf = 0 then the coordinates are ordered such that values      */
  /* sequentially increase along the channel. The user must discern  */
  /* which is the case and set MERGE_DIR accordingly.                */
  if (dirf) {
    for (i = 0; i < ncells; i++) {
      n = i * 6;
      lon[i][0] = x[n+4];
      lat[i][0] = y[n+4];
      lon[i][1] = x[n];
      lat[i][1] = y[n];
      lon[i][2] = x[n+2];
      lat[i][2] = y[n+2];
      lon[i][3] = x[n+6];
      lat[i][3] = y[n+6];
      lon[i][4] = x[n+8];
      lat[i][4] = y[n+8];
      if(verbose && i==0)for(n=0; n<=4; n++)printf("%f %f y%d\n",lon[i][n],lat[i][n],n);
    }
  } else {
    for (i = 0; i < ncells; i++) {
      n = i * 2;
      lon[i][0] = x[nc/3+n+1];
      lat[i][0] = y[nc/3+n+1];
      lon[i][1] = x[n];
      lat[i][1] = y[n];
      lon[i][2] = x[n+2];
      lat[i][2] = y[n+2];
      lon[i][3] = x[2*nc/3+n];
      lat[i][3] = y[2*nc/3+n];
      lon[i][4] = x[2*nc/3+n+2];
      lat[i][4] = y[2*nc/3+n+2];
      if(verbose && i==0)for(n=0; n<=4; n++)printf("%f %f n%d\n",lon[i][n],lat[i][n],n);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the merge location to the closest vertices in the grid      */
  for (i = 0; i < ncells; i++) {
    sort_circle_g(lon[i], lat[i], 4, sortdir);
    for (j = 1; j <= 4; j++) {
      d = sqrt((lon[i][j] - xr1) * (lon[i][j] - xr1) + (lat[i][j] - yr1) * (lat[i][j] - yr1));
      if (d < dist1) {
	dist1 = d;
	ir = i;
	jr1 = j;
      }
    }
  }
  lon[ir][jr1] = xr1;
  lat[ir][jr1] = yr1;
  if (verbose) {
    printf("New start merge: i=%d j=%d\n", ir, jr1);
    printf("Old end merge: j=%d [%f %f]\n", jr2, xr2, yr2);
  }
  /*
  xr1 = xc[ir][jr1];
  yr1 = yc[ir][jr1];
  */
  for (j = 0; j <= 4; j++) {
    d = sqrt((lon[ir][j] - xr2) * (lon[ir][j] - xr2) + (lat[ir][j] - yr2) * (lat[ir][j] - yr2));
    if (j > 0 && j != jr1 && d < dist2) {
      dist2 = d;
      jr2 = j;
    }
  }
  lon[ir][jr2] = xr2;
  lat[ir][jr2] = yr2;
  hd_warn("Curvilinear grid %s added at (%f %f)-(%f %f)\n", iname, xr1, yr1, xr2, yr2);

  /*-----------------------------------------------------------------*/
  /* Merge a second set of coordinates if required (e.g. for two     */
  /* ends of a channel.                                              */
  if (prm_read_char(fp, "MERGE_LOC2", buf)) {
    sscanf(buf, "(%lf %lf)-(%lf %lf)", &xr1, &yr1, &xr2, &yr2);

    /* The merge location is set to the closest vertices in the      */
    /* mesh.                                                         */
    dist1 = dist2 = HUGE;
    for (i = 1; i <= params->ns2; i++) {
      for (j = 0; j <= params->npe2[i]; j++) {
	d = sqrt((xc[i][j] - xr1) * (xc[i][j] - xr1) + (yc[i][j] - yr1) * (yc[i][j] - yr1));
	if (j > 0 && d < dist1) {
	  dist1 = d;
	  ir = i;
	  jr1 = j;
	}
      }
    }
    xr1 = xc[ir][jr1];
    yr1 = yc[ir][jr1];
    for (i = 1; i <= params->ns2; i++) {
      for (j = 0; j <= params->npe2[i]; j++) {
	d = sqrt((xc[i][j] - xr2) * (xc[i][j] - xr2) + (yc[i][j] - yr2) * (yc[i][j] - yr2));
	if (j > 0 && d < dist2) {
	dist2 = d;
	ir = i;
	jr2 = j;
	}
      }
    }
    xr2 = xc[ir][jr2];
    yr2 = yc[ir][jr2];
    hd_warn("  Second coordinates to merge are (%f %f)-(%f %f)\n", xr1, yr1, xr2, yr2);

    /* Set the merge location to the closest vertices in the grid    */
    dist1 = dist2 = HUGE;
    for (i = 0; i < ncells; i++) {
      /*sort_circle_g(lon[i], lat[i], 4, sortdir);*/
      for (j = 1; j <= 4; j++) {
	d = sqrt((lon[i][j] - xr1) * (lon[i][j] - xr1) + (lat[i][j] - yr1) * (lat[i][j] - yr1));
	if (d < dist1) {
	  dist1 = d;
	  ir = i;
	  jr1 = j;
	}
      }
    }
    lon[ir][jr1] = xr1;
    lat[ir][jr1] = yr1;
    /*
    xr1 = xc[ir][jr1];
    yr1 = yc[ir][jr1];
    */
    for (j = 0; j <= 4; j++) {
      d = sqrt((lon[ir][j] - xr2) * (lon[ir][j] - xr2) + (lat[ir][j] - yr2) * (lat[ir][j] - yr2));
      if (j > 0 && j != jr1 && d < dist2) {
	dist2 = d;
	jr2 = j;
      }
    }
    lon[ir][jr2] = xr2;
    lat[ir][jr2] = yr2;
    hd_warn("Curvilinear grid %s added at second location (%f %f)-(%f %f)\n", iname, xr1, yr1, xr2, yr2);
  }

  ns2 = params->ns2;
  /* Sometimes a few cells from the mesh should be removed so that   */
  /* the quad grid can seamlessly interface with the mesh.           */
  /* Read polygons to remove some cells if required                  */
  if (strlen(polyr)) {
    char *fields[MAXSTRLEN * MAXNUMARGS];
    int npoly = parseline(polyr, fields, MAXNUMARGS);
    FILE *fp;
    poly_t *pl = poly_create();
    ns2 = params->ns2;
    for (n = 0; n < npoly; n++) {
      if ((fp = fopen(fields[n], "r")) == NULL) {
	hd_warn("Can't find file %s to remove from quad addition %s\n", 
		fields[0], iname);
	continue;
      }
      j = poly_read(pl, fp);
      for (i = 1; i <= ns2; i++) {
	if (poly_contains_point(pl, params->x[i][0], params->y[i][0])) {
	  mask[i] = i;
	  params->ns2--;
	}
      }
      fclose(fp);
      poly_destroy(pl);
    }
  }

  /* Reallocate and reset the (x,y) arrays, including the points in  */
  /* the grid.                                                       */
  params->ns2 += ncells;
  i_free_1d(params->npe2);
  d_free_1d(params->bathy);
  d_free_2d(params->x);
  d_free_2d(params->y);
  params->npe2 = i_alloc_1d(params->ns2+1);
  params->bathy = d_alloc_1d(params->ns2+1);
  params->x = d_alloc_2d(params->npe+1, params->ns2+1);
  params->y = d_alloc_2d(params->npe+1, params->ns2+1);
  n = 1;
  for (i = 1; i <= ns2; i++) {
    if (mask[i]) continue;
    params->npe2[n] = npe2[i];
    params->bathy[n] = bathy[i];
    for (j = 0; j <= params->npe2[n]; j++) {
      params->x[n][j] = xc[i][j];
      params->y[n][j] = yc[i][j];
    }
    n++;
  }
  for (i = 0; i < ncells; i++) {
    params->npe2[n] = 4;
    params->bathy[n] = -fabs(nbathy);
    for (j = 0; j <= params->npe2[n]; j++) {
      params->x[n][j] = lon[i][j];
      params->y[n][j] = lat[i][j];
    }
    n++;
  }

  fclose(fp);
  d_free_1d(x);
  d_free_1d(y);
  d_free_1d(bathy);
  i_free_1d(npe2);
  d_free_2d(lon);
  d_free_2d(lat);
  d_free_2d(xc);
  d_free_2d(yc);
  i_free_1d(mask);
}

/* END add_quad_grido()                                              */
/*-------------------------------------------------------------------*/

/*
 * Forward stereographic projection of a jigsaw_msh structure
 *
 * Note : Only meshes are supported here
 */
static void st_transform(jigsaw_msh_t *J_msh, ST3PROJ kind,
			 double xmid, double ymid)
{
  int n;
  double deg2m = 60.0 * 1852.0;
  
  if (J_msh->_flags == JIGSAW_EUCLIDEAN_MESH) {
    jigsaw_VERT2_array_t *vert2 = &J_msh->_vert2;
    int sz = vert2->_size;
    int vsz = J_msh->_value._size;
    if (sz) {
      /* Do the transformation */
      for (n=0; n<sz; n++) {
	double x, y, scl;
	if (kind == ST3PROJ_FWD) {
	  /* Angles will be in lat/lon degrees so need to convert */
	  stereo3_fwd(DEG2RAD(vert2->_data[n]._ppos[0]),
		      DEG2RAD(vert2->_data[n]._ppos[1]),
		      xmid, ymid, &x, &y, &scl);

	  /* Apply value scaling if available */
	  if (vsz)
	    J_msh->_value._data[n] *= deg2m*scl;

	  /* Swap out new values */
	  vert2->_data[n]._ppos[0] = x;
	  vert2->_data[n]._ppos[1] = y;
	  
	} else if (kind == ST3PROJ_INV) {
	  /* Angles in the mesh will have been converted to radians */
	  stereo3_inv(vert2->_data[n]._ppos[0],
		      vert2->_data[n]._ppos[1],
		      xmid, ymid, &x, &y, NULL);

	  /* Swap out new values - what about wrapping? */
	  vert2->_data[n]._ppos[0] = RAD2DEG(x);
	  vert2->_data[n]._ppos[1] = RAD2DEG(y);
	  
	} else
	  hd_quit("st_transform: Unknown stereographic projection kind supplied");
      }
    } else
      hd_quit("st_fwd_transform: VERT2 is empty");
  } else
    hd_quit("st_fwd_transform: Only JIGSAW_EUCLIDEAN_MESH type is supported\n");
}


/* 
 * stereographic projection functions based on Darren Engwirda's 
 * jigsaw-geo-matlab/script/geom-util/stereo3.m
 */
static void stereo3_fwd(double xpos, double ypos, double xmid, double ymid,
			double *xnew, double *ynew, double *scal)
{
  double kden, kval;

  kden = +1. + sin(ymid)*sin(ypos) +
    cos(ymid)*cos(ypos)*cos(xpos-xmid) ;
  
  kval = +2. * RADIUS / kden;

  /* Calculate and assign new values */
  *xnew = kval * cos(ypos)*sin(xpos-xmid) ;
  *ynew = kval * (cos(ymid)*sin(ypos) -
		  sin(ymid)*cos(ypos)*cos(xpos-xmid));

  /* Calculate scaling, if asked */
  if (scal != NULL) {
    double dval, dkdy, dXdy, dYdy;
    
    dval = -2. * RADIUS / (kden * kden) ;
    dkdy = dval * (sin(ymid) * cos(ypos) -
		   cos(ymid) * sin(ypos)*cos(xpos-xmid));
    
    dXdy = dkdy * cos(ypos) * sin(xpos-xmid) -
      kval * sin(ypos)*sin(xpos-xmid) ;
  
    dYdy = dkdy * (cos(ymid) * sin(ypos) -
		   sin(ymid)*cos(ypos)*cos(xpos-xmid)) +
      kval*(cos(ymid) * cos(ypos) + sin(ymid)*sin(ypos)*cos(xpos-xmid));

    /* Assign output */
    *scal = +1./RADIUS * sqrt(dXdy*dXdy + dYdy*dYdy) ;
  }
}

/* Inverse stereographic */
static void stereo3_inv(double xpos, double ypos, double xmid, double ymid,
			double *xnew, double *ynew, double *scal)
{
  double rval, cval;
  
  rval = sqrt(xpos*xpos + ypos*ypos);
      
  cval = +2. * atan2(rval, +2. * RADIUS);

  /* Assign inverse transformations */
  *xnew = xmid + atan2(xpos * sin(cval),
		       rval * cos(ymid) * cos(cval) -
		       ypos * sin(ymid) * sin(cval));
  
  *ynew = asin(cos(cval)*sin(ymid) + 
	       (ypos*sin(cval)*cos(ymid)) / rval);


  /* Calculate scaling, if asked */
  if (scal != NULL) {
    double dval, kden, kval, dkdy, dXdy, dYdy;

    kden = +1. + sin(ymid)*sin(*xnew) + 
      cos(ymid)*cos(*ynew)*cos(*xnew-xmid) ;
    
    dval = -2. * RADIUS / (kden*kden) ;
    kval = +2. * RADIUS / kden ;
    
    dkdy = dval*(sin(ymid) * cos(*ynew) - 
		 cos(ymid) * sin(*ynew) * cos(*xnew-xmid));
    
    dXdy = dkdy * cos(*ynew) * sin(*xnew-xmid) -
      kval * sin(*ynew) * sin(*xnew-xmid);
  
    dYdy = dkdy * (cos(ymid) * sin(*ynew) - 
		   sin(ymid) * cos(*ynew)*cos(*xnew-xmid)) +
      kval*(cos(ymid) * cos(*ynew) + 
	    sin(ymid)*sin(*ynew)*cos(*xnew-xmid));
    
    *scal = RADIUS * +1/sqrt(dXdy*dXdy + dYdy*dYdy);
  }
}


/*
 * Converts JIGSAW x/y grid into tria3 mesh
 */
static void convert_jigsaw_grid_to_mesh(jigsaw_msh_t *imsh,
					jigsaw_msh_t *omsh)
{
  /* Input grid */
  jigsaw_REALS_array_t *xgrid = &imsh->_xgrid;
  jigsaw_REALS_array_t *ygrid = &imsh->_ygrid;
  jigsaw_REALS_array_t *ivals = &imsh->_value;

  /* Output mesh */
  jigsaw_TRIA3_array_t *tria3 = &omsh->_tria3;
  jigsaw_REALS_array_t *ovals = &omsh->_value;
  jigsaw_VERT2_array_t *vert2 = &omsh->_vert2;

  int nx = xgrid->_size;
  int ny = ygrid->_size;
  int i,j,n, ntria, nvert;

  /* Allocate memory for output mesh */
  ntria = 2*(nx-1)*(ny-1);
  jigsaw_alloc_tria3(tria3, ntria);

  nvert = nx*ny;
  jigsaw_alloc_vert2(vert2, nvert);
  jigsaw_alloc_reals(ovals, nvert);
  
  /* Must be x/y grid */
  if (nx == 0 || ny == 0)
    hd_quit("convert_jigsaw_grid_to_mesh: x/y grid empty!");
  
  /* Create vertices */
  n = 0;
  for (i=0; i<nx; i++) {
    for (j=0; j<ny; j++) {
      vert2->_data[n]._ppos[0] = xgrid->_data[i];
      vert2->_data[n]._ppos[1] = ygrid->_data[j];
      vert2->_data[n]._itag = 0;
      ovals->_data[n] = ivals->_data[n];
      n++;
    }
  }

  /* 
   * Traverse the grid and split each rectangular cell into two
   * triangles
   */
  n = 0;
  for (i=0; i<nx-1; i++) {
    for (j=0; j<ny-1; j++) {
      /* 
       * Start in the lower left hand corner, n0
       *
       * n1     n2
       *  o --- o
       *  |   / |   ^
       *  |  /  |   |
       *  | /   |   y
       *  o --- o    x->
       * n0     n3
       *
       */
      int xn0 = i,   yn0 = j;
      int xn1 = i,   yn1 = j+1;
      int xn2 = i+1, yn2 = j+1;
      int xn3 = i+1, yn3 = j;
      int n0 = ny*xn0 + yn0;
      int n1 = ny*xn1 + yn1;
      int n2 = ny*xn2 + yn2;
      int n3 = ny*xn3 + yn3;
      
      /* Clockwise node : n0->n1->n2 */
      tria3->_data[n]._node[0] = n0;
      tria3->_data[n]._node[1] = n1;
      tria3->_data[n]._node[2] = n2;
      n++;
      
      /* Anto-clockwise node : n0->n3->n2 */
      tria3->_data[n]._node[0] = n0;
      tria3->_data[n]._node[1] = n3;
      tria3->_data[n]._node[2] = n2;
      n++;
    }
  }

  /* Flag output mesh type */
  omsh->_flags = JIGSAW_EUCLIDEAN_MESH;

  /* Delete old mesh */
  jigsaw_free_msh_t(imsh);
  
}

point p0;
/*-------------------------------------------------------------------*/  
/* A utility function to swap two points                             */
/*-------------------------------------------------------------------*/  
void swap(point *p1, point *p2)
{
  point temp = *p1;
  *p1 = *p2;
  *p2 = temp;
}

/*-------------------------------------------------------------------*/  
/* A utility function to return square of distance between p1 and p2 */
/*-------------------------------------------------------------------*/  
double distSq(point p1, point p2)
{
  return((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y));
}

/*-------------------------------------------------------------------*/  
/* To find orientation of ordered triplet (p, q, r).                 */
/* The function returns following values:                            */
/* 0 --> p, q and r are colinear                                     */
/* 1 --> Clockwise                                                   */
/* 2 --> Counterclockwise                                            */
/*-------------------------------------------------------------------*/  
int orientation(point p, point q, point r)
{
  double val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
  if (val == 0) return 0;
  return (val > 0) ? 1 : 2;
}

/*-------------------------------------------------------------------*/
/* A function used by library function qsort() to sort an array of   */
/* points with respect to the first point.                           */
/*-------------------------------------------------------------------*/
int compare(const void *vp1, const void *vp2) 
{
  point *p1 = (point *)vp1;
  point *p2 = (point *)vp2;

  int o = orientation(p0, *p1, *p2);
  if (o == 0)
    return (distSq(p0, *p2) >= distSq(p0, *p1))? -1 : 1;

  return (o == 2) ? -1: 1;
}

/*-------------------------------------------------------------------*/
/* Finds convex hull of a set of n points using Graham Scan          */
/* algorithm.                                                        */
/* Adapted from:                                                     */
/* https://www.geeksforgeeks.org/convex-hull-set-2-graham-scan/      */
/* https://stackoverflow.com/questions/37635258/convex-hull-in-c     */
/*-------------------------------------------------------------------*/
point *convex_hull(point *v,  int *count)
{
  int n = *count;
  double ymin = v[0].y;
  int min = 0;
  int i, m;
  point *stack;

  /* Pick the bottom-most or chose the left most point in case of a  */
  /* tie.                                                            */
  for(i = 1; i < n; i++) {
    if((v[i].y < ymin) || ((v[i].y == ymin) && (v[i].x < v[min].x))) {
      ymin = v[i].y;
      min = i;
    }
  }

  /* Place the bottom-most point at first position                   */
  swap(&v[0], &v[min]);

  /* Sort n-1 points with respect to the first point. A point p1     */
  /* comes before p2 in sorted output if p2 has larger polar angle   */
  /* (in counterclockwise direction) than p1.                        */
  p0 = v[0];
  if(n > 1)
    qsort(&v[1], n - 1, sizeof(point), compare);

  /* If two or more points make same angle with p0, remove all but   */
  /* the one that is farthest from p0. Remember that, in above       */
  /* sorting, our criteria was to keep the farthest point at the end */
  /* when more than one points have same angle.                      */
  m = 1; /* Initialize size of modified array.                       */
  for(i = 1; i < n; i++) {
    while((i < n - 1) && orientation(v[0], v[i], v[i + 1]) == 0)
      i++;
    v[m++] = v[i];
  }
  *count = n = m;

  /* If modified array of points has less than 3 points, convex hull */
  /* is not possible.                                                */
  if(n < 3) return v;

  /* Allocate the hull points and push first three points to it.     */
  stack = (point *)malloc(n * sizeof(point));
  stack[0] = v[0];
  stack[1] = v[1];
  stack[2] = v[2];

  /* Process remaining n-3 points                                    */
  m = 2;
  for(i = 3; i < n; i++) {
    /* Keep removing top while the angle formed by points next-to-   */
    /* top, top, and points[i] makes a non-left turn.                */
    while(orientation(stack[m-1], stack[m], v[i]) != 2) {
      m--;
    }
    stack[++m] = v[i];
  }

  *count = n = ++m;

  return stack;
}

/*-------------------------------------------------------------------*/
/* Finds concave hull of a set of n points by finding the perimeter  */
/* of a Delaunay triangulation of those points.                      */
/*-------------------------------------------------------------------*/
point *concave_hull_c(point *v,  int *count)
{
  int n = *count;
  delaunay *d;
  int **nei, *pe;
  int i, npe = 0;
  int *index, *mask;
  point *stack;

  d = delaunay_build(n, v, 0, NULL, 0, NULL);
  neighbour_finder_b(d, &nei);

  /* Count the perimeter edges. These will always be less than the   */
  /* number of trianges.                                             */
  npe = 0;
  pe = i_alloc_1d(d->ntriangles); 
  for (i = 0; i < d->ntriangles; ++i) {
    for (n = 0; n < 3; n++) {
      if (nei[n][i] == -1) pe[npe++] = i;
    }
  }
  /* Get the indices of the perimeter edge points. A mapping == -1   */
  /* indicates that triange lies on the perimeter of the             */
  /* triangulation (an interior triangle should have 3 valid maps).  */
  index = i_alloc_1d(2 * npe);
  for (i = 0; i < npe; i++) {
    triangle* t = &d->triangles[pe[i]];
    for (n = 0; n < 3; n++) {
      if (nei[n][i] == -1) {
	int nn = (n == 2) ? 0 : n + 1;
	index[2*i] = t->vids[n];
	index[2*i+1] = t->vids[nn];
      }
    }
  }
  /* Remove duplicates and make a continuous polygon                 */
  mask = i_alloc_1d(2 * npe);
  memset(mask, 0, 2 * npe * sizeof(int));
  n = 0;
  pe[n] = -1;
  mask[n++] = index[0];
  mask[n++] = index[1];
  for (i = 0; i < npe; i++) {
    if (pe[i] != -1) {
      if (mask[n-1] == index[2*i]) {
	mask[n++] = index[2*i+1];
	pe[i] = -1;
	break;
      }
      if (mask[n-1] == index[2*i+1]) {
	mask[n++] = index[2*i];
	pe[i] = -1;
	break;
      }
    }
  }
  n--;

  /* Save the point locations                                        */
  stack = (point *)malloc((n + 1)* sizeof(point));
  for (i = 0; i < n; i++) {
    stack[i].x = d->points[mask[i]].x;
    stack[i].y = d->points[mask[i]].y;
  }
  n++;
  stack[n].x = stack[0].x;    /* Close the polygon                   */
  stack[n].y = stack[0].y;
  *count = n;
  i_free_1d(pe);
  i_free_1d(index);
  i_free_1d(mask);
  i_free_2d(nei);
  return stack;
}




