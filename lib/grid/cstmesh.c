/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/grid/cstmesh.c
 *
 *  \brief Coastal mesh routines as a preprocessor for JIGSAW
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: cstmesh.c 5860 2018-07-02 04:07:46Z her127 $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <sys/time.h>
#include <unistd.h>
#include "ems.h"

#define DEG2RAD(d) ((d)*M_PI/180.0)
#define RAD2DEG(r) ((r)*180.0/M_PI)

#define RADIUS 6370997.0
#define ECC 0.0
#define GEODESIC(x1, y1, x2, y2) (geod_inv_geod_fwd_sodanos(DEG2RAD(x1), \
							    DEG2RAD(y1), \
                                 DEG2RAD(x2), DEG2RAD(y2),\
                                 RADIUS, ECC))

static void order (int *a, int *b, int *c1, int *c2);
static int find_index(double slon, double slat, double *lon, double *lat, int nsl, double *d);
static void find_link_index(int ns, int *nsl, double **lat, double **lon, link_t *links, int dir);
static void connect_link(int ns, int nlink, int *nsl, int *nlm, double **lat, double **lon, int il, link_t *link);
static int find_link(int en, int eid, int *n, int *si, int *dir, int *size, int *nsl, int nlink,
		     link_t *link, double **lat, double **lon, double *nlat, double *nlon, int *flag);
static int find_dir(link_t *links, int *nsl, int n);
static msh_t *msh_alloc(int ndims, int npoint);
static void msh_free(msh_t *msh);

static jig_t *jig_alloc(char *infile, char *hfile, double hmin, double hmax);
static void jig_free(jig_t *jig);
// xxx What about jig_free?

/*
 * PUBLIC API's *
 */

/*---------------------------------------------------------------------*/
/* Allocates coastmesh memory                                          */
/*---------------------------------------------------------------------*/
coamsh_t *cm_alloc(void)
{
  coamsh_t *cm = (coamsh_t *)malloc(sizeof(coamsh_t));
  memset(cm, 0, sizeof(coamsh_t));
  return cm;
}

/* END cm_alloc()                                                      */
/*---------------------------------------------------------------------*/


/*---------------------------------------------------------------------*/
/* Frees coastmesh memory                                              */
/*---------------------------------------------------------------------*/
void cm_free(coamsh_t *cm)
{
  int m;

  if (cm->ss) i_free_2d(cm->ss);
  if (cm->sr) i_free_1d(cm->sr);
  if (cm->nobc) {
    for (m = 0; m < cm->nobc; m++) {
      cmobc_t *open = &cm->obc[m];
      d_free_1d(open->olon);
      d_free_1d(open->olat);
    }
    free(cm->obc);
  }
  if (cm->nlink) {
    for (m = 0; m < cm->nlink; m++) {
      link_t *links = &cm->link[m];
      i_free_1d(links->slink);
      i_free_1d(links->ilink);
      d_free_1d(links->llon);
      d_free_1d(links->llat);
      d_free_1d(links->alon);
      d_free_1d(links->alat);
      if (links->nseg) {
	d_free_1d(links->segx);
	d_free_1d(links->segy);
      }
    }
    free(cm->link);   
  }
  if (cm->msh) {
    msh_free(cm->msh);
  }
  free((coamsh_t *)cm);
}

/* END cm_free()                                                       */
/*---------------------------------------------------------------------*/

/*---------------------------------------------------------------------*/
/* Builds a msh specification from a coastmesh definition.             */
/*---------------------------------------------------------------------*/
int coastmesh(coamsh_t *cm, int mode)
{
  FILE *fp, *op = NULL, *bp = NULL, *pp = NULL, *opp = NULL;
  char buf[MAXSTRLEN];
  int n, m, i, j, nl, nlm, sm;
  int n1, n2, ii, i1, i2, i3, mi;
  int msl;           /* Length of the major segment                    */
  int sid, eid, mid; /* Start, end and mid index for the major segment */
  int obcf;          /* =0; sid closest to OBC start, =1; eid closest  */
  int obci;          /* Start segment for minor segments               */
  int obcs;          /* End index for major edges                      */
  int obcm;          /* Maximum obc size                               */
  double jd, sjd, njd;
  double **lat, **lon;
  double **rlat, **rlon;
  double *nlat, *nlon;
  double *tlon, *tlat;
  int *mask;          /* Mask for segments to include in the mesh      */
  int ns;             /* Number of segments in the coastfile           */
  int *nsl;           /* Number of points in each segment              */
  int *nsp;           /* Index of first coordinate in each segment     */
  int *nso;           /* Re-ordered segment list                       */
  int *nsr;           /* Reverse map of nso                            */
  int point;          /* Total number of points                        */
  int obcint = 0;     /* Integrate the OBC path into the major segment */
  int **flag;         /* Type of pont in the list (coast, link, obc)   */
  int *incf;
  double d1, d2, dist, dist1, dist2;
  int isclosed = 1;
  int checki = 0;
  int filef = 0;
  int cutoff;
  link_t *link = NULL;
  int dir;
  int sseg, ssid, ssize, sdir, sn, en, mn;
  msh_t *msh = NULL;
  jig_t *jig = NULL;
  int psize = 10;

  /*-------------------------------------------------------------------*/
  /* Open the ouput file                                               */
  if (mode & CM_FILE) {
    if((op = fopen(cm->ofile, "w")) == NULL) {
      warn("coastmesh: Can't open output file %s\n", cm->ofile);
    } else {
      fprintf(op, "# .msh geometry file\n");
      fprintf(op, "mshid=2\n");
      fprintf(op, "ndims=2\n");
      filef = 1;
    }
  }

  /*-------------------------------------------------------------------*/
  /* Read the input parameters                                         */
  if((fp = fopen(cm->cofile, "r")) == NULL) {
    quit("coastmesh: Can't open coastline file %s\n", cm->cofile);
  }

  pp = opp = NULL;
  if (strlen(cm->pname)) {
    sprintf(buf, "%s_in.txt", cm->pname);
    if((pp = fopen(buf, "w")) == NULL) {
      warn("coastmesh: Can't open plot file %s\n", buf);
    }
    sprintf(buf, "%s_out.txt", cm->pname);
    if((opp = fopen(buf, "w")) == NULL) {
      warn("coastmesh: Can't open plot file %s\n", buf);
    }
  }

  /*-------------------------------------------------------------------*/
  /* Read in the bounding box                                          */
  if (cm->bslon > cm->belon) {
    d1 = cm->bslon;
    cm->bslon = cm->belon;
    cm->belon = cm->bslon;
  }
  if (cm->bslat > cm->belat) {
    d1 = cm->bslat;
    cm->bslat = cm->belat;
    cm->belat = cm->bslat;
  }

  /*-------------------------------------------------------------------*/
  /* Open boundary paths                                               */
  obcm = 0;
  for (n = 0; n < cm->nobc; n++) {
    cmobc_t *open = &cm->obc[n];
    if (open->nobc > obcm) obcm = open->nobc;
  }

  /*-------------------------------------------------------------------*/
  /* Get the number of segments and maximum segment length             */
  m = ns = nl = nlm = n1 = 0;
  while (fgets(buf, 256, fp) != NULL) {
    char *fields[MAXSTRLEN * 100];
    if (strlen(buf) <= 1) continue;
    n = parseline(buf, fields, 100);
    if (m % cm->resample == 0) nl++;
    m++;
    if (strcmp(fields[0], "nan") == 0 ||
	strcmp(fields[0], "NaN") == 0) {
      if (nl > psize) {
	nlm = max(nlm, nl);
	ns ++;
      }
      nl = m = 0;
      n1++;
    }
  }

  /* Allocate                                                          */
  rlon = d_alloc_2d(nlm, ns);
  rlat = d_alloc_2d(nlm, ns);
  nsl = i_alloc_1d(ns);
  nsp = i_alloc_1d(ns);
  nso = i_alloc_1d(ns);
  mask = i_alloc_1d(ns);
  incf = i_alloc_1d(n1);
  flag = i_alloc_2d(nlm+obcm, ns);
  memset(incf, 0, n1 * sizeof(int));
  for (i = 0; i < ns; i++)
    mask[i] = 1;

  rewind(fp);
  nl = m = n1 = 0;
  while (fgets(buf, 256, fp) != NULL) {
    char *fields[MAXSTRLEN * 100];
    if (strlen(buf) <= 1) continue;
    n = parseline(buf, fields, 100);
    if (m % cm->resample == 0) nl++;
    m++;
    if (strcmp(fields[0], "nan") == 0 ||
	strcmp(fields[0], "NaN") == 0) {
      if (nl > psize) {
	incf[n1] = 1; 
      }
      nl = m = 0;
      n1++;
    }
  }

  /* Populate the coorinate arrays for each segment                    */
  rewind(fp);
  ns = nl = n1 = 0;
  while (fgets(buf, 256, fp) != NULL) {
    char *fields[MAXSTRLEN * 100];
    if (strlen(buf) <= 1) continue;
    n = parseline(buf, fields, 100);
    if (strcmp(fields[0], "nan") == 0 ||
	strcmp(fields[0], "NaN") == 0) {
      if (incf[n1]) {
	nsl[ns] = nl;
	ns++;
      }
      nl = m = 0;
      n1++;
    } else {
      if (incf[n1]) {
	if (m % cm->resample == 0) {
	  rlon[ns][nl] = atof(fields[0]);
	  rlat[ns][nl] = atof(fields[1]);
	  nl++;
	}
	m++;
      }
    }
  }

  /*-------------------------------------------------------------------*/
  /* Find the start, end , mid index and segment this occurs in.       */
  /* Start index.                                                      */
  dist1 = 1e10;
  sdir = 0;
  for (n = 0; n < ns; n++) {
    if ((i = find_index(cm->slon, cm->slat, rlon[n], rlat[n], nsl[n], &dist)) >= 0) {
      if (dist < dist1) {
	dist1 = dist;
	sid = i;
	sn = n;
      }
    }
  }
  if (sid >= 0) {
    /* End and mid indices                                             */
    dist1 = dist2 = 1e10;
    for (n = 0; n < ns; n++) {
      if ((i = find_index(cm->elon, cm->elat, rlon[n], rlat[n], nsl[n], &dist)) >= 0) {
	if (dist < dist1) {
	  dist1 = dist;
	  eid = i;
	  en = n;
	}
      }
    }
    if (sid >= 0 && eid >= 0) {
      obcint = 1;
      for (n = 0; n < ns; n++) {
	if ((i = find_index(cm->mlon, cm->mlat, rlon[n], rlat[n], nsl[n], &dist)) >= 0) {
	  if (dist < dist2) {
	    dist2 = dist;
	    mid = i;
	    mn = n;
	  }
	}
      }
    }
    if (sn != mn) fprintf(pp, "Error start segment %d is not the same as the mid segment %d\n", sn, mn);

    /* Find the direction to cycle the start index                     */
    sdir = -1;
    for( i = sid; i <= nsl[sn]; i++) {
      if (i == mid) {
	sdir = 1;
	break;
      }
    }
  }

  /*-------------------------------------------------------------------*/
  /* Find the link indices                                             */
  for (i = 0; i < cm->nlink; i++) {
    link_t *links = &cm->link[i];
    /* Find the indices of the links (segment specific)                */
    links->id = i;
    find_link_index(ns, nsl, rlat, rlon, links, sdir);
  }

  /*-------------------------------------------------------------------*/
  /* Create the links                                                  */
  if (obcint) {
    /* Count the number of points in the major segment                 */
    ssize = 0;
    ssid = sid;
    sseg = sn;
    dir = sdir;
    i = 0;
    while (ssid != eid || sseg != en) {
      m = find_link(en, eid, &sseg, &ssid, &dir, &ssize, nsl, cm->nlink, cm->link, rlat, rlon, NULL, NULL, NULL);
      /*printf("Pass%d sid = %d, seg = %d, size = %d, link = %d\n", i, ssid, sseg, ssize, m);*/
      if(i > cm->nlink) {
	printf("Can't find end index %d in segment %d\n", sid, sn);
	quit("coastmesh error\n");
      }
      i++;
    }

    /* Count the points from the last link to the end point. To do     */
    /* this we make a dummy link, where the end point and segment are  */
    /* placed in the last link, then call find_link() to capture these */
    /* last points. After counting, reset the last link to its         */
    /* original specification.                                         */
    /*
    n = cm->nlink-1;
    i1 = cm->link[n].ilink[0];
    i2 = cm->link[n].ilink[1];
    n1 = cm->link[n].slink[0];
    n2 = cm->link[n].slink[1];
    cm->link[n].slink[0] = cm->link[n].slink[1] = en;
    cm->link[n].ilink[0] = cm->link[n].ilink[1] = eid;
    m = find_link(en, eid, &sseg, &ssid, &dir, &ssize, nsl, cm->nlink, cm->link, rlat, rlon, NULL, NULL, NULL);
    cm->link[n].ilink[0] = i1;
    cm->link[n].ilink[1] = i2;
    cm->link[n].slink[0] = n1;
    cm->link[n].slink[1] = n2;
    */
    msl = ssize;

    /* Assign the points to the major segment                          */
    nlat = d_alloc_1d(msl);
    nlon = d_alloc_1d(msl);
    ssize = 0;
    ssid = sid;
    sseg = sn;
    dir = sdir;
    i = 0;
    while (ssid != eid || sseg != en) {
      m = find_link(en, eid, &sseg, &ssid, &dir, &ssize, nsl, cm->nlink, cm->link, rlat, rlon, nlat, nlon, flag[sn]);
      i++;
    }
    /*
    cm->link[n].slink[0] = cm->link[n].slink[1] = en;
    cm->link[n].ilink[0] = cm->link[n].ilink[1] = eid;
    m = find_link(en, eid, &sseg, &ssid, &dir, &ssize, nsl, cm->nlink, cm->link, rlat, rlon, nlat, nlon, flag[sn]);
    */
    /* Clean                                                           */
    for(i = 1; i < msl-1; i++) {
      if (nlon[i] == 0.0) nlon[i] = 0.5 * (nlon[i-1] + nlon[i+1]);
      if (nlat[i] == 0.0) nlat[i] = 0.5 * (nlat[i-1] + nlat[i+1]);
    }

    /* Copy the linked and un-linked segments to the lat/long arrays   */
    nlm = max(nlm, msl);
    lon = d_alloc_2d(nlm, ns);
    lat = d_alloc_2d(nlm, ns);
    for (n = 0; n < ns; n++) {
      if (n == sn) {
	for (i = 0; i < msl; i++) {
	  lat[n][i] = nlat[i];
	  lon[n][i] = nlon[i];
	}
	nsl[n] = msl;
      } else {
	for (i = 0; i < nsl[n]; i++) {
	  lat[n][i] = rlat[n][i];
	  lon[n][i] = rlon[n][i];
	  flag[n][i] = S_COAST;
	}
      }
    }

    /* Remove all segments in the link list other than the segment     */
    /* hosting the start index. This latter segment now contains the   */
    /* points in all those other segments.                             */
    for (n = 0; n < cm->nlink; n++) {
      link_t *links = &cm->link[n];
      m = links->slink[0];
      if (m != sn) nsl[m] = 0;
      m = links->slink[1];
      if (m != sn) nsl[m] = 0;
    }
    d_free_2d(rlat);
    d_free_2d(rlon);
  } else {
    /* Connect any links                                               */
    connect_link(ns, cm->nlink, nsl, &nlm, lat, lon, i, cm->link);
    for (n = 0; n < ns; n++) {
      for (i = 0; i < nsl[n]; i++) {
	lat[n][i] = rlat[n][i];
	lat[n][i] = rlat[n][i];
      }
    }
    d_free_2d(rlat);
    d_free_2d(rlon);
  }

  if (pp) {
    fprintf(pp, "Original coastline\n");
    for (n = 0; n < ns; n++) {
      ii = 0;
      for(i = 1; i < nsl[n]; i++) {
	if (lon[n][i] >= cm->bslon && lon[n][i] <= cm->belon && lat[n][i] >= cm->bslat && lat[n][i] <= cm->belat) {
	  fprintf(pp, "%f %f\n", lon[n][i], lat[n][i]);
	  ii = 1;
	}
      }
      if (ii) fprintf(pp, "NaN NaN\n");
    }
  }

  /*-------------------------------------------------------------------*/
  /* Pre-smooth if required                                            */
  for (m = 0; m < cm->ismooth; m++) {
    for (n = 0; n < ns; n++) {
      if (nsl[n] == 0) continue;
      tlon = d_alloc_1d(nsl[n]);
      tlat = d_alloc_1d(nsl[n]);
      memcpy(tlon, lon[n], nsl[n] * sizeof(double));
      memcpy(tlat, lat[n], nsl[n] * sizeof(double));
      for(i = 1; i < nsl[n]; i++) {
	double d3 = 1.0;
	d1 = tlon[i];
	d2 = tlat[i];
	for (j = 1; j <= cm->ismoothz; j++) {
	  ii = max(0, i-j);
	  d1 += tlon[ii];
	  d2 += tlat[ii];
	  ii = min(nsl[n]-1, i+j);
	  d1 += tlon[ii];
	  d2 += tlat[ii];
	  d3 += 2.0;
	}
	lon[n][i] = d1 / d3;
	lat[n][i] = d2 / d3;
      }
      d_free_1d(tlon);
      d_free_1d(tlat);
    }
  }

  /*-------------------------------------------------------------------*/
  /* Get the longest segment and order segment sizes                   */
  sm = 0;
  for (n = 0; n < ns; n++) {
    if (nsl[n] > sm) sm = nsl[n];
    nso[n] = n;
  }
  for (i = 0; i < ns; i++)
    for (j = ns-1; i < j; --j)
      order(&nsl[j-1], &nsl[j], &nso[j-1], &nso[j]);

  /* Map a map of reordered segments                                   */
  nsr = i_alloc_1d(ns);
  for (n = 0; n < ns; n++) {
    for (m = 0; m < ns; m++) {
      if (n == nso[m]) nsr[n] = m;
    }
  }
  /* Exclude the maor segment from the segment list                    */
  mask[nsr[sn]] = 0;

  /*-------------------------------------------------------------------*/
  /* Find the coordinate indices.                                      */
  /* Find the start coordinate.                                        */
  /*n = nso[ns-1];*/
  n = nsr[sn];
  if (sid >= 0)
    sid = find_index(cm->slon, cm->slat, lon[sn], lat[sn], nsl[n], NULL);

  /* Find the end coordinate                                           */
  if (eid >= 0)
    eid = find_index(cm->elon, cm->elat, lon[sn], lat[sn], nsl[n], NULL);

  /* Find the mid coordinate                                           */
  if (mid >= 0)
    mid = find_index(cm->mlon, cm->mlat, lon[sn], lat[sn], nsl[n], NULL);

  /*-------------------------------------------------------------------*/
  /* Pre-smooth individual segments if required                        */
  for (n1 = 0; n1 < cm->nss; n1++) {
    n = cm->ss[0][n1];
    n2 = nso[n];
    for (m = 0; m < cm->ss[1][n1]; m++) {
      if (nsl[n] == 0) continue;
      tlon = d_alloc_1d(nsl[n]);
      tlat = d_alloc_1d(nsl[n]);
      memcpy(tlon, lon[n2], nsl[n] * sizeof(double));
      memcpy(tlat, lat[n2], nsl[n] * sizeof(double));
      for(i = 1; i < nsl[n]; i++) {
	double d3 = 1.0;
	d1 = tlon[i];
	d2 = tlat[i];
	for (j = 1; j <= cm->ss[2][n1]; j++) {
	  ii = max(0, i-j);
	  d1 += tlon[ii];
	  d2 += tlat[ii];
	  ii = min(nsl[n]-1, i+j);
	  d1 += tlon[ii];
	  d2 += tlat[ii];
	  d3 += 2.0;
	}
	lon[n2][i] = d1 / d3;
	lat[n2][i] = d2 / d3;
      }
      d_free_1d(tlon);
      d_free_1d(tlat);
    }
  }

  /*-------------------------------------------------------------------*/
  /* Remove individual segments if required                            */
  for (n1 = 0; n1 < cm->nrs; n1++) {
    n = cm->sr[n1];
    nsl[n] = 0;
  }

  /*-------------------------------------------------------------------*/
  /* Print summary to plotfile                                         */
  if (pp) {
    fprintf(pp, "Number of segments = %d\n", ns);
    fprintf(pp, "Maximum segment size = %d\n", nlm);
    fprintf(pp, "Number of OBCs = %d\n", cm->nobc);
    /*n = nso[ns-1];*/
    n = sn;
    if (sid >= 0)
      fprintf(pp, "Start index = %d (%f %f)\n",sid, lon[n][sid], lat[n][sid]);
    if (eid >= 0)
      fprintf(pp, "End index = %d (%f %f)\n",eid, lon[n][eid], lat[n][eid]);
    if (mid >= 0)
      fprintf(pp, "Mid index = %d (%f %f)\n",mid, lon[n][mid], lat[n][mid]);
    fprintf(pp, "%d pre-smoothing passes with a width of %d\n", cm->ismooth, cm->ismoothz);
    fprintf(pp, "%d post-smoothing passes with a width of %d\n", cm->osmooth, cm->osmoothz);
    for (n = 0; n < cm->nss; n++) {
      fprintf(pp, "Segment %d has %d pre-smoothing passes with a width of %d\n", cm->ss[0][n], cm->ss[1][n], cm->ss[2][n]);
    }
    fprintf(pp, "Re-sampled every %d points\n", cm->resample);
    if (cm->nlink) fprintf(pp, "Links=%d:\n", cm->nlink);
    for (m = 0; m < cm->nlink; m++) {
      link_t *links = &cm->link[m];
      fprintf(pp, "link%d\n", m);
      if (strlen(links->text)) fprintf(pp, "  %s\n", links->text);
      i = links->ilink[0];
      /*n = nso[links->slink[0]];*/
      n = links->slink[0];
      fprintf(pp, "  Start index(%d) in segment %d: %d=(%4.8f %4.8f)\n", i, n, i, links->alon[0], links->alat[0]);
      i = links->ilink[1];
      /*n = nso[links->slink[1]];*/
      n = links->slink[1];
      fprintf(pp, "  End index(%d) in segment %d: %d=(%4.8f %4.8f)\n", i, n, i, links->alon[1], links->alat[1]);
      if (links->flag & (LINK_M|LINK_P)) {
	i = links->ilink[2];
	fprintf(pp, "  Mid index(%d): %d=(%4.8f %4.8f)\n", i, i, links->alon[2], links->alat[2]);
      }
      if (links->nseg)
	fprintf(pp, "  Path with %d points included in LINK%d\n", links->nseg, m);
    }
  }

  /*-------------------------------------------------------------------*/
  /* If no start/end index is supplied, the boundary forms a closed    */
  /* segment within which all othe segments lie.                       */
  if (sid == -1 && eid == -1) {
    cmobc_t *open = &cm->obc[0];
    if (!cm->nobc) {
      quit("coastmesh: Must supply a closed open boundary if no start/end coordinates are specified.\n");
    }
    obci = 1;
    obcs = 1;
    /* Include the OBC as the 1st segment                              */
    msl = open->nobc;
    nlon = d_alloc_1d(msl);
    nlat = d_alloc_1d(msl);
    /* Copy the locations of OBC0                                      */
    open = &cm->obc[0];
    for (i = 0; i < open->nobc; i++) {
      nlon[i] = open->olon[i];
      nlat[i] = open->olat[i];
    }
    cm->nobc = 0;
  } else {
    /* Copy the major segment back into nlon,nlat to account for any   */
    /* pre-processing of lon,lat (e,g, smoothing).                     */
    for (i = 0; i < msl; i++) {
      nlon[i] = lon[sn][i];
      nlat[i] = lat[sn][i];
    }

    /*-----------------------------------------------------------------*/
    /* Get the OBC index closest to the start index                    */
    obcf = 0;
    if (cm->nobc) {
      d1 = nlon[0] - cm->obc[0].olon[0];
      d2 = nlat[0] - cm->obc[0].olat[0];
      dist1 = sqrt(d1 * d1 + d2 * d2);
      d1 = nlon[msl-1] - cm->obc[0].olon[0];
      d2 = nlat[msl-1] - cm->obc[0].olat[0];
      dist2 = sqrt(d1 * d1 + d2 * d2);
      if (dist2 < dist1) obcf = 1;
    }
    obci = 2;
    obcs = 0;
  }

  /*-------------------------------------------------------------------*/
  /* Segments are written to the coastfile as closed loops, so we do   */
  /* not need to list the last coordinate entry in the segment; i.e.   */
  /* only loop to nsl[n] - 1. Decrement nsl by one to account for      */
  /* this.                                                             */
  if (isclosed) {
    for (n = 0; n < ns; n++) {
      j = nso[n];
      if(lon[j][0] == lon[j][nsl[n]-1] && lat[j][0] == lat[j][nsl[n]-1]) {
	nsl[n]--;
      }
    }
  }

  /*-------------------------------------------------------------------*/
  /* Eliminate segments outside the bounding box, less than the cutoff */
  /* and with mean radius greater than the minimum.                    */
  cutoff = ns * cm->cutoff / 100;
  if (cutoff < 0) cutoff = ns - 1;
  /* Exclude the 2nd segment in links                                  */
  /*
  for (m = 0; m < cm->nlink; m++) {
    link_t *links = &cm->link[m];
    n = nsr[links->slink[1]];
    if (nsr[links->slink[0] != n])
      mask[n] = 0;
  }
  */
  /* Exclude segments less than the cutoff size                        */
  for (n = 0; n < cutoff; n++)
    mask[n] = 0;
  /* Exclude segments with any coordinates outside the bounding box    */
  /*for (n = 0; n < ns-1; n++) {*/
  for (n = 0; n < ns; n++) {
    if (!mask[n]) continue;
    j = nso[n];
    for (i = 0; i < nsl[n]; i++) {
      if (lon[j][i] < cm->bslon || lon[j][i] > cm->belon || lat[j][i] < cm->bslat || lat[j][i] > cm->belat)
	mask[n] = 0;
    }
  }
  /* Exclude segments with a mean radius less than the minimum         */
  if (cm->radius < HUGE) {
    if (pp) fprintf(pp, "Minimum radius = %f (m)\n", cm->radius);
    /*for (n = 0; n < ns-2; n++) {*/
    for (n = 0; n < ns; n++) {
      double clon, clat;
      if (!mask[n]) continue;
      /* Get the centre of mass                                        */
      clon = clat = 0.0;
      for (i = 0; i < nsl[n]; i++) {
	j = nso[n];
	clon += lon[j][i];
	clat += lat[j][i];
      }
      clon /= (double)nsl[n];
      clat /= (double)nsl[n];
      /* Get the mean radius                                           */
      dist1 = 0.0;
      for (i = 0; i < nsl[n]; i++) {
	j = nso[n];
	dist1 += GEODESIC(lon[j][i], lat[j][i], clon, clat);
      }
      dist1 /= (double)nsl[n];
      if (pp && cm->radius < HUGE) fprintf(pp, "Segment %d radius = %f (m) from [%f %f]\n", n, dist1, clon, clat);
      if (dist1 < cm->radius) mask[n] = 0;
    }
  }
  /* Exclude segments with a maximum length less than the minimum      */
  if (cm->length < HUGE) {
    if (pp) fprintf(pp, "Minimum length = %f (m)\n", cm->length);
    /*for (n = 0; n < ns-2; n++) {*/
    for (n = 0; n < ns; n++) {
      if (!mask[n]) continue;
      dist1 = 0.0;
      j = nso[n];
      for (i = 0; i < nsl[n]; i++) {
	for (ii = i + 1; ii < nsl[n]; ii++) {
	  dist2 = GEODESIC(lon[j][ii], lat[j][ii], lon[j][i], lat[j][i]);
	  if (dist2 > dist1) dist1 = dist2;
	}
      }
      if (pp && cm->length < HUGE) fprintf(pp, "Segment %d length = %f (m)\n", n, dist1);
      if (dist1 < cm->length) mask[n] = 0;
    }
  }

  /*-------------------------------------------------------------------*/
  /* Get the number of points                                          */
  point = msl;
  /*nsp[ns-1] = 0;*/
  nsp[nsr[sn]] = 0;
  for (i = 0; i < cm->nobc; i++) {
    cmobc_t *open = &cm->obc[i];
    point += open->nobc;
  }
  /*for (n = ns-obci; n >= 0; n--) {*/
  for (n = ns-1; n >= 0; n--) {
    if (!mask[n]) continue;
    j = nso[n];
    nsp[n] = point;
    point += nsl[n];
  }
  if (filef) fprintf(op, "point=%d\n", point);

  /*-------------------------------------------------------------------*/
  /* Set the OBC endpoints to the start and end index (for continuity  */
  /* since the start/end coordinates of nlon/nlat may have changed due */
  /* to smoothing.                                                     */
  for (n2 = 0; n2 < cm->nobc; n2++) {
    cmobc_t *open = &cm->obc[n2];
    if (obcf == 0) {
      open->olon[0] = nlon[0];
      open->olat[0] = nlat[0];
      open->olon[open->nobc-1] = nlon[msl-1];
      open->olat[open->nobc-1] = nlat[msl-1];
    } else {
      open->olon[0] = nlon[msl-1];
      open->olat[0] = nlat[msl-1];
      open->olon[open->nobc-1] = nlon[0];
      open->olat[open->nobc-1] = nlat[0];
    }
  }

  /*-------------------------------------------------------------------*/
  /* Integrate the obc path into the major segment if required         */
  if (obcint) {
    tlon = d_alloc_1d(msl + obcm);
    tlat = d_alloc_1d(msl + obcm);
    for(i = 0; i < msl; i++) {
      tlon[i] = nlon[i];
      tlat[i] = nlat[i];
    }
    for (m = 0; m < cm->nobc; m++) {
      cmobc_t *open = &cm->obc[m];
      if (obcf == 0) {
	for (ii = open->nobc-1; ii >= 0; ii--) {
	  tlon[i] = open->olon[ii];
	  tlat[i] = open->olat[ii];
	  flag[sn][i] = S_OBC;
	  i++;
	}
      } else {
	for (ii = 0; ii < open->nobc; ii++) {
	  tlon[i] = open->olon[ii];
	  tlat[i] = open->olat[ii];
	  flag[sn][i] = S_OBC;
	  i++;
	}
      }
    }
    msl = i;

    d_free_1d(nlon);
    d_free_1d(nlat);
    nlon = tlon;
    nlat = tlat;
  }

  /*-------------------------------------------------------------------*/
  /* Post-smooth the major segment if required                         */
  for (m = 0; m < cm->osmooth; m++) {
    /*
    for (n = 0; n < ns; n++) {
      if (nsl[n] == 0) continue;
      if (!mask[n]) continue;
    */
    tlon = d_alloc_1d(msl);
    tlat = d_alloc_1d(msl);
    for(i = 0; i < msl; i++) {
      tlon[i] = nlon[i];
      tlat[i] = nlat[i];
    }
    for(i = 0; i < msl; i++) {
      double d3 = 1.0;
      d1 = tlon[i];
      d2 = tlat[i];
      i1 = i2 = i;
      for (j = 1; j <= cm->osmoothz; j++) {
	i1++;
	i2--;
	if (i1 == msl) i1 = 0;
	if (i2 < 0) i2 = msl - 1;
	/*
	i1 = max(0, i1);
	i2 = min(msl-1, i2);
	*/
	d1 += (tlon[i1] + tlon[i2]);
	d2 += (tlat[i1] + tlat[i2]);
	d3 += 2.0;
      }
      nlon[i] = d1 / d3;
      nlat[i] = d2 / d3;
    }
    d_free_1d(tlon);
    d_free_1d(tlat);
  }

  /*-------------------------------------------------------------------*/
  /* Print the coordinates in order of decreasing segment size above   */
  /* the cutoff fraction.                                              */
  /* Major segment coordinates                                         */
  /*if (opp) fprintf(opp, "Processed coastline\n");*/
  mi = 0;
  msh = msh_alloc(2, point);
  cm->np = (obcint) ? msl : msl + obcm;
  cm->x = d_alloc_1d(cm->np);
  cm->y = d_alloc_1d(cm->np);
  for(i = 0; i < msl; i++) {
    if (filef) {
      if (checki)
	fprintf(op, "%d %f %f\n", i, nlon[i], nlat[i]);
      else
	fprintf(op, "%f;%f;0\n", nlon[i], nlat[i]);
    }
    if (opp) fprintf(opp, "%f %f\n", nlon[i], nlat[i]);
    msh->coords[0][mi] = cm->x[mi] = nlon[i];
    msh->coords[1][mi] = cm->y[mi] = nlat[i];
    msh->flag[mi] = flag[sn][i];
    mi++;
  }

  /* Boundary path, from the location closest to the end coordinate to */
  /* that closest to the start coordinate.                             */
  if (!obcint) {
    m = msl;
    for (i = 0; i < cm->nobc; i++) {
      cmobc_t *open = &cm->obc[i];
      if (obcf == 0) {
	for (n = open->nobc-1; n >= 0; n--) {
	  if (filef) {
	    if (checki)
	      fprintf(op, "%d %f %f\n", m++, open->olon[n], open->olat[n]);
	    else
	      fprintf(op, "%f;%f;0\n", open->olon[n], open->olat[n]);
	  }
	  if (opp) fprintf(opp, "%f %f\n", open->olon[n], open->olat[n]);
	  msh->coords[0][mi] = cm->x[mi] = open->olon[n]; 
	  msh->coords[1][mi] = cm->y[mi] = open->olat[n];
	  msh->flag[mi] = S_OBC;
	  mi++;
	}
      } else {
	for (n = 0; n < open->nobc; n++) {
	  if (filef) {
	    if (checki)
	      fprintf(op, "%d %f %f\n", m++, open->olon[n], open->olat[n]);
	    else
	      fprintf(op, "%f;%f;0\n", open->olon[n], open->olat[n]);
	  }
	  if (opp) fprintf(opp, "%f %f\n", open->olon[n], open->olat[n]);
	  msh->coords[0][mi] = cm->x[mi] = open->olon[n]; 
	  msh->coords[1][mi] = cm->y[mi] = open->olat[n];
	  msh->flag[mi] = S_OBC;
	  mi++;
	}
      }
    }
  }
  if (opp){
    fprintf(opp, "%f %f\n", nlon[0], nlat[0]);
    fprintf(opp, "NaN NaN\n");
  }

  /* Other segments.                                                   */
  /*for (n = ns-obci; n >= cutoff; n--) {*/
  for (n = ns-1; n >= cutoff; n--) {
    if (!mask[n]) continue;
    j = nso[n];
    for (i = 0; i < nsl[n]; i++) {
      if (filef) {
	if (checki)
	  fprintf(op, "%d %f %f\n", m++, lon[j][i], lat[j][i]);
	else
	  fprintf(op, "%f;%f;0\n", lon[j][i], lat[j][i]);
      }
      if (opp) fprintf(opp, "%f %f\n", lon[j][i], lat[j][i]);
      msh->coords[0][mi] = lon[j][i]; 
      msh->coords[1][mi] = lat[j][i];
      mi++;
    }
    if (opp) {
      fprintf(opp, "%f %f\n", lon[j][0], lat[j][0]);
      fprintf(opp, "NaN NaN\n");
    }
  }

  /*-------------------------------------------------------------------*/
  /* Print the edges                                                   */
  msh->edge2 = point;
  if (filef) fprintf(op, "edge2=%d\n", point);
  /* Major segment                                                     */
  j = mi = 0;
  for(i = 0; i < msl-obcs-1; i++) {
    if (filef) fprintf(op, "%d;%d;0\n", j, j+1);
    msh->edges[0][mi] = j;
    msh->edges[1][mi] = j + 1;
    j++;
    mi++;
  }
  /* Close the major loop                                              */
  if (filef) fprintf(op, "%d;0;0\n", j);
  msh->edges[0][mi] = j;
  msh->edges[1][mi] = 0;
  j++;
  mi++;

  /* Boundary path. Note j+1 is the last coordinate in the path, and   */
  /* this must close to the start index, hence only loop to            */
  /* open->nobc-1.                                                     */
  if (!obcint) {
    for (i = 0; i < cm->nobc; i++) {
      cmobc_t *open = &cm->obc[i];
      for (n = 0; n < open->nobc-1; n++) {
	if (filef) fprintf(op, "%d;%d;0\n", j, j+1);
	msh->edges[0][mi] = j;
	msh->edges[1][mi] = j + 1;
	j++;
	mi++;
      }
    }
    /* Close the boundary path loop                                    */
    if (filef) fprintf(op, "%d;0;0\n", j);
    msh->edges[0][mi] = j;
    msh->edges[1][mi] = 0;
    j++;
    mi++;
  }

  /* Other segments                                                    */
  /*for (n = ns-obci; n >= cutoff; n--) {*/
  for (n = ns-1; n >= cutoff; n--) {
    if (!mask[n]) continue;
    m = nso[n];
    for (i = 0; i < nsl[n]-1; i++) {
      if (filef) {
	if (checki)
	  fprintf(op, "%d %d %d;%d;0\n", n, i, j, j+1);
	else
	  fprintf(op, "%d;%d;0\n", j, j+1);
      }
      msh->edges[0][mi] = j;
      msh->edges[1][mi] = j + 1;
      j++;
      mi++;
    }

    /* Explititly close the loop. Note; segments are written to the    */
    /* coastfile as closed loops, so we do not need to list the last   */
    /* coordinate entry in the segment; i.e. only loop to nsl[n] - 2   */
    /* since j++ takes care of nsl[n]-1.                               */
    if (filef) {
      if (checki)
	fprintf(op, "%d %d %d;%d;0\n",n, i, j, nsp[n]);
      else
	fprintf(op, "%d;%d;0\n",j, nsp[n]);
    }
    msh->edges[0][mi] = j;
    msh->edges[1][mi] = nsp[n];
    j++;
    mi++;
  }

  /*-------------------------------------------------------------------*/
  /* Populate the jig structure and write to file                      */
  /*-------------------------------------------------------------------*/
  if (filef) {
    jig = jig_alloc(cm->ofile, cm->hfunf, cm->hfun_min, cm->hfun_max);
    write_jigsaw(jig);
  }

  /* Close and exit                                                    */
  fclose(fp);
  if (op) fclose(op);
  if (pp) fclose(pp);
  if (opp) fclose(opp);
  i_free_1d(nsl);
  i_free_1d(nsp);
  i_free_1d(nso);
  i_free_1d(nsr);
  i_free_1d(mask);
  d_free_2d(lon);
  d_free_2d(lat);
  i_free_2d(flag);
  /*
  if (tlon) d_free_1d(tlon);
  if (tlat) d_free_1d(tlat);
  */
  if (nlon) d_free_1d(nlon);
  if (nlat) d_free_1d(nlat);
  /*if (flag) i_free_2d(flag);*/
  if (incf) i_free_1d(incf);
  cm->msh = msh;
  cm->jig = jig;
}

/* END coastmesh()                                                     */
/*---------------------------------------------------------------------*/

/*---------------------------------------------------------------------*/
/* Reads a coastmesh input specification from file ip                  */
/*---------------------------------------------------------------------*/
void cm_read(coamsh_t *cm, FILE *ip)
{
  FILE *bp;
  char fname[MAXSTRLEN], obcname[MAXSTRLEN];
  char buf[MAXSTRLEN], key[MAXSTRLEN], key1[MAXSTRLEN];
  int n, m, i, j;
  int sid, eid, mid; /* Start, end and mid index for the major segment */
  double bslon, bslat, belon, belat; /* Bounding box                   */


  /*-------------------------------------------------------------------*/
  /* Read the input parameters                                         */
  if (prm_read_char(ip, "JIG_GEOM_FILE", fname))
    strcpy(cm->ofile, fname);

  if (prm_read_char(ip, "COASTFILE", fname))
    strcpy(cm->cofile, fname);

  sprintf(cm->pname, "%c", '\0');
  if (prm_read_char(ip, "PLOTFILE", fname))
    strcpy(cm->pname, fname);

  sprintf(cm->hfunf, "%c", '\0');
  if (prm_read_char(ip, "HFUN_FILE", fname))
    strcpy(cm->hfunf, fname);

  /*-------------------------------------------------------------------*/
  /* Smoothing and resampling                                          */
  cm->ismooth = 0;
  cm->ismoothz = 1;
  cm->osmooth = 0;
  cm->osmoothz = 1;
  if (prm_read_int(ip, "SMOOTH_PRE", &cm->ismooth) || 
      prm_read_int(ip, "SMOOTH", &cm->ismooth)) {
    if (!(prm_read_int(ip, "SMOOTHZ_PRE", &cm->ismoothz)))
      prm_read_int(ip, "SMOOTHZ", &cm->ismoothz);
  }
  if (prm_read_int(ip, "SMOOTH_POST", &cm->osmooth)) {
    prm_read_int(ip, "SMOOTHZ_POST", &cm->osmoothz);
  }
  cm->nss = 0;
  if (prm_read_char(ip, "SMOOTH_SEGMENT", buf)) {
    char *fields[MAXSTRLEN * 100];
    cm->nss = parseline(buf, fields, 100);
    cm->ss = i_alloc_2d(cm->nss,3);
    for (n = 0; n < cm->nss; n++) {
      char *tok;
      i = 0;
      tok = strtok(fields[n], ":");
      while (tok != NULL) {
	cm->ss[i++][n] = atoi(tok);
	tok = strtok(NULL, ":");
	if (i == 3) break;
      }
    }
  }
  if (prm_read_char(ip, "REMOVE_SEGMENT", buf)) {
    char *fields[MAXSTRLEN * 100];
    cm->nrs = parseline(buf, fields, 100);
    cm->sr = i_alloc_1d(cm->nrs);
    for (n = 0; n < cm->nrs; n++) {
      cm->sr[n] = atoi(fields[n]);
    }
  }
  cm->resample = 1;
  prm_read_int(ip, "RESAMPLE", &cm->resample);
  cm->radius = HUGE;
  prm_read_double(ip, "RADIUS", &cm->radius);
  cm->length = HUGE;
  prm_read_double(ip, "LENGTH", &cm->length);
  cm->hfun_min = 0.01;
  prm_read_double(ip, "HFUN_HMIN", &cm->hfun_min);
  cm->hfun_max = 0.05;
  prm_read_double(ip, "HFUN_HMAX", &cm->hfun_max);

  /*-------------------------------------------------------------------*/
  /* Read in the bounding box                                          */
  cm->cutoff = 50;
  prm_read_int(ip, "CUTOFF", &cm->cutoff);
  prm_read_double(ip, "MINLON", &cm->bslon);
  prm_read_double(ip, "MINLAT", &cm->bslat);
  prm_read_double(ip, "MAXLON", &cm->belon);
  prm_read_double(ip, "MAXLAT", &cm->belat);

  /*-------------------------------------------------------------------*/
  /* Read the open boundary paths                                      */
  cm->nobc = 0;
  if (prm_read_int(ip, "NOBC", &cm->nobc)) {
    cm->obc = (cmobc_t *)malloc(sizeof(cmobc_t) * cm->nobc);
    /* Count the bounray path entries */
    for (i = 0; i < cm->nobc; i++) {
      cmobc_t *open = &cm->obc[i];
      sprintf(buf, "OBC%d", i);
      if (prm_read_char(ip, buf, obcname)) {
	if((bp = fopen(obcname, "r")) == NULL) {
	  printf("Can't open boundary path file %s\n", obcname);
	  exit(0);
	}
	open->nobc = 0;
	while (fgets(buf, 256, bp) != NULL) {
	  open->nobc++;
	}
	fclose(bp);
      } else {
	printf("Can't find keyword %s\n", buf);
	exit(0);
      }
      open->olon = d_alloc_1d(open->nobc);
      open->olat = d_alloc_1d(open->nobc);
    }

    /* Read the boundary path entries                                */
    for (i = 0; i < cm->nobc; i++) {
      cmobc_t *open = &cm->obc[i];
      sprintf(buf, "OBC%d", i);
      if (prm_read_char(ip, buf, obcname)) {
	bp = fopen(obcname, "r");
	for (n = 0; n < open->nobc; n++) {
	  if (!fscanf(bp, "%lf %lf", &open->olon[n], &open->olat[n]))
	    quit("cm_read: OBC read error\n");
	}
	fclose(bp);
      }
    }
  }

  /*-------------------------------------------------------------------*/
  /* Read the links                                                    */
  cm->nlink = 0;
  cm->link = NULL;
  if (prm_read_int(ip, "LINKS", &cm->nlink)) {
    if (cm->nlink) {
      cm->link = (link_t *)malloc(sizeof(link_t) * cm->nlink);
      for (i = 0; i < cm->nlink; i++) {
	char *fields[MAXSTRLEN * 100];
	link_t *links = &cm->link[i];
	links->llon = d_alloc_1d(3);
	links->llat = d_alloc_1d(3);
	links->alon = d_alloc_1d(3);
	links->alat = d_alloc_1d(3);
	links->slink = i_alloc_1d(2);
	links->ilink = i_alloc_1d(3);
	sprintf(key, "LINK%d.start", i);
	if (prm_read_char(ip, key, buf)) {
	  n = parseline(buf, fields, 100);
	  links->llon[0] = atof(fields[0]);
	  links->llat[0] = atof(fields[1]);
	} else {
	  printf("Can't find keyword %s\n", key);
	  exit(0);
	}
	sprintf(key, "LINK%d.end", i);
	if (prm_read_char(ip, key, buf)) {
	  n = parseline(buf, fields, 100);
	  links->llon[1] = atof(fields[0]);
	  links->llat[1] = atof(fields[1]);
	} else {
	  printf("Can't find keyword %s\n", key);
	  exit(0);
	}
	links->dir = 0;
	links->flag = 0;
	sprintf(key, "LINK%d.mid", i);
	sprintf(key1, "LINK%d.dir", i);
	if (prm_read_char(ip, key, buf)) {
	  n = parseline(buf, fields, 100);
	  links->llon[2] = atof(fields[0]);
	  links->llat[2] = atof(fields[1]);
	  links->flag = LINK_M;
	} else if (prm_read_char(ip, key1, buf)) {
	  if (strcmp(buf, "+") == 0)
	    links->dir = 1;
	  if (strcmp(buf, "-") == 0)
	    links->dir = -1;
	  links->flag = LINK_P;
	} else {
	  links->flag = LINK_N;
	  links->dir = 0;
	}
	sprintf(links->text, "%c", '\0');
	sprintf(key, "LINK%d.text", i);
	prm_read_char(ip, key, links->text);
	sprintf(key, "LINK%d.path", i);
	links->nseg = 0;
	if (prm_read_char(ip, key, buf)) {
	  FILE *sp;
	  if ((sp = fopen(buf, "r")) != NULL) {
	    links->nseg= 0;
	    while (fgets(buf, MAXSTRLEN, sp) != NULL) {
	      double d1, d2;
	      j = parseline(buf, fields, 100);
	      d1 = atof(fields[0]);
	      d2 = atof(fields[1]);
	      if (!isnan(d1) && !isnan(d2)) links->nseg++;
	    }
	    rewind(sp);
	    links->segx = d_alloc_1d(links->nseg);
	    links->segy = d_alloc_1d(links->nseg);
	    n = 0;
	    while (fgets(buf, MAXSTRLEN, sp) != NULL) {
	      j = parseline(buf, fields, 100);
	      links->segx[n] = atof(fields[0]);
	      links->segy[n] = atof(fields[1]);
	      if (!isnan(links->segx[n]) && !isnan(links->segy[n])) n++;
	    }
	    fclose(sp);
	  } else
	    printf("Can't open link segment file %s\n", buf);
	}
      }
    }
  }

  /*-------------------------------------------------------------------*/
  /* Read in the start, end and mid coordinates                        */
  if (prm_read_char(ip, "START_COORD", buf)) {
    char *fields[MAXSTRLEN * 100];
    n = parseline(buf, fields, 100);
    if (n ==2 ) {
      cm->slon = atof(fields[0]);
      cm->slat = atof(fields[1]);
    } else {
      printf("Start lon/lat incorrectly specified.\n");
      exit(0);
    }
    sid = 0;
  } else {
    sid = -1;
  }
  if (prm_read_char(ip, "END_COORD", buf)) {
    char *fields[MAXSTRLEN * 100];
    n = parseline(buf, fields, 100);
    if (n ==2 ) {
      cm->elon = atof(fields[0]);
      cm->elat = atof(fields[1]);
    } else {
      printf("End lon/lat incorrectly specified.\n");
      exit(0);
    }
    eid = 0;
  } else {
    eid = -1;
  }
  if (prm_read_char(ip, "MID_COORD", buf)) {
    char *fields[MAXSTRLEN * 100];
    n = parseline(buf, fields, 100);
    if (n ==2 ) {
      cm->mlon = atof(fields[0]);
      cm->mlat = atof(fields[1]);
    } else {
      printf("Mid lon/lat incorrectly specified.\n");
      exit(0);
    }
    mid = 0;
  } else {
    mid = -1;
  }
}

/* END cm_read()                                                       */
/*---------------------------------------------------------------------*/

/*---------------------------------------------------------------------*/
/* Writesd a JIGSAW compatible input file                              */
/*---------------------------------------------------------------------*/
void write_jigsaw(jig_t *jig)
{
  FILE *fp;
  int i;
  char buf[MAXSTRLEN], jigfile[MAXSTRLEN];
  double hfun_hmax;
  double hfun_hmin;
  double km2nm = 1852.0;
  int verbose = 0;

  hfun_hmin = jig->hfun_hmin / (60.0 * km2nm);
  hfun_hmax = jig->hfun_hmax / (60.0 * km2nm);
 
  for (i = 0; i < strlen(jig->geom_file); i++) {
    if (jig->geom_file[i] == '.')
      break;
    else
      buf[i] = jig->geom_file[i];
  }
  buf[i] = '\0';
  sprintf(jigfile, "%s.jig", buf);
  if ((fp = fopen(jigfile, "w")) == NULL) {
    printf("Can't open file %s\n", jigfile);
    exit(0);
  }
  fprintf(fp,"# JIGSAW input file for %s\n", buf);
  fprintf(fp, "geom_file=%s\n", jig->geom_file);
  fprintf(fp, "mesh_file=%s\n", jig->mesh_file);
  fprintf(fp, "\n# Number of seed vertices used to init-ialise mesh generation. Default = 8.\n");
  fprintf(fp, "geom_seed=%d\n", jig->geom_seed);
  fprintf(fp, "\n# Attempt to auto-detect sharp features in the input geometry.\n");
  fprintf(fp, "# Features can be adjacent to geometry 'edges' based on both\n");
  fprintf(fp, "# geometrical and/or topological constraints.\n");
  fprintf(fp, "# Geometrically, features are located between any neighbouring\n");
  fprintf(fp, "# entities that subtend angles less than GEOM_ETA1 degrees.\n");
  if (jig->geom_feat)
    fprintf(fp, "geom_feat=true\n");
  else
    fprintf(fp, "geom_feat=false\n");
  fprintf(fp, "\n# Feature-angle, features are located between any neighbouring\n");
  fprintf(fp, "# edges that subtend angles less than GEOM_ETA1 degrees. Default = 45deg.\n");
  fprintf(fp, "geom_eta1=%d\n", jig->geom_eta1);
  fprintf(fp, "geom_eta2=%d\n", jig->geom_eta2);
  /*
  fprintf(fp, "\n# Mesh-size kernal, choice between a constant size-function\n");
  fprintf(fp, "# (KERN='constant') and a Delaunay-based medial-axis method\n");
  fprintf(fp, "# (KERN='delaunay') that attempts to automatically generate geometry-adaptive\n");
  fprintf(fp, "# sizing data. Default = constant.\n");
  fprintf(fp, "hfun_kern=%s\n", jig->hfun_kern);
  */
  fprintf(fp, "\n# Scaling type for mesh-size fuction. SCAL='relative' interprets mesh-size\n");
  fprintf(fp, "# values as percentages of the (mean) length of the axis-aligned bounding-box\n");
  fprintf(fp, "# for the geometry. SCAL='absolute' interprets mesh-size values as absolute\n");
  fprintf(fp, "#  measures. Default = relative.\n");
  fprintf(fp, "hfun_scal=%s\n", jig->hfun_scal);
  fprintf(fp, "\n# Maximum and minimum mesh size function value. Default = 0.02 & 0.00.\n");
  fprintf(fp, "hfun_file=%s\n", jig->hfun_file);
  fprintf(fp, "hfun_hmax=%f\n", hfun_hmax);
  fprintf(fp, "hfun_hmin=%f\n", hfun_hmin);
  fprintf(fp, "\n# Maximum gradient in the mesh-size function. Default = 0.25.\n");
  fprintf(fp, "hfun_grad=%f\n", jig->hfun_grad);
  fprintf(fp, "\n# Number of topological dimensions to mesh. DIMS=K meshes K-dimensional\n");
  fprintf(fp, "# features, irrespective of the number of spatial dimensions of the problem\n");
  fprintf(fp, "# (i.e. if the geometry is 3-dimensional and DIMS=2 a surface mesh will be produced).\n");
  fprintf(fp, "# Default = 3.\n");
  fprintf(fp, "mesh_dims=%d\n", jig->mesh_dims);
  fprintf(fp, "\n# Meshing kernal, choice of the standard Delaunay-refinement algorithm \n");
  fprintf(fp, "# (KERN='delaunay') or the Frontal-Delaunay method (KERN='delfront').\n");
  fprintf(fp, "# default = delfront.\n");
  fprintf(fp, "mesh_kern=%s\n", jig->mesh_kern);
  fprintf(fp, "\n# Maximum number of mesh refinement iterations. Set ITER=N to see progress\n");
  fprintf(fp, "# after N iterations.\n");
  fprintf(fp, "# mesh_iter = %d\n", jig->mesh_iter);
  fprintf(fp, "\n# Enforce 1-dim. topological constraints. 1-dim. edges are refined until all\n");
  fprintf(fp, "# embedded nodes are 'locally 1-manifold', i.e. nodes are either centred at\n");
  fprintf(fp, "# topological 'features', or lie on 1-manifold complexes. Default = false.\n");
  if (jig->mesh_top1)
    fprintf(fp, "mesh_top1=true\n");
  else
    fprintf(fp, "mesh_top1=false\n");
  if (jig->mesh_top2)
    fprintf(fp, "mesh_top2=true\n");
  else
    fprintf(fp, "mesh_top2=false\n");
  fprintf(fp, "\n# Maximum radius-edge ratio for 2-tria elements. 2-trias are refined until the\n");
  fprintf(fp, "# ratio of the element circumradius to min. edge length is less-than MESH_RAD2.\n");
  fprintf(fp, "# Default = 1.05.\n");
  fprintf(fp, "mesh_rad2=%f\n", jig->mesh_rad2);
  fprintf(fp, "mesh_rad3=%f\n", jig->mesh_rad3);
  fprintf(fp, "\n# Maximum surface-discretisation error multiplier for 1-edge elements.\n");
  fprintf(fp, "# 1-edge elements are refined until the surface-disc. error is less-than\n");
  fprintf(fp, "# MESH_EPS1 * HFUN(X). Default = 0.33.\n");
  fprintf(fp, "mesh_eps1=%f\n", jig->mesh_eps1);
  fprintf(fp, "\n# Maximum surface-discretisation error multiplier for 2-tria elements.\n");
  fprintf(fp, "# 2-tria elements are refined until the surface-disc. error is less-than\n");
  fprintf(fp, "# MESH_EPS2 * HFUN(X). Default = 0.20\n");
  fprintf(fp, "mesh_eps2=%f\n", jig->mesh_eps2);
  fprintf(fp, "mesh_vol3=%f\n", jig->mesh_vol3);
  fprintf(fp, "\n# Verbosity of log-file output generated by JIGSAW. Set VERBOSITY>0 to print\n");
  fprintf(fp, "# additional output. Default = 0.\n");
  fprintf(fp, "verbosity=%d\n", jig->verbosity);
  fclose(fp);
  if (verbose) printf("\nRun 'jigsaw64r %s'\n\n", jigfile);
}

/* END write_jigsaw()                                                  */
/*---------------------------------------------------------------------*/

/*
 * PRIVATE API's *
 */

/*---------------------------------------------------------------------*/
/* Frees msh memory                                                    */
/*---------------------------------------------------------------------*/
static void msh_free(msh_t *msh)
{
  int m;

  d_free_2d(msh->coords);
  i_free_2d(msh->edges);
  i_free_1d(msh->flag);
  free((msh_t *)msh);
}

/* END msh_free()                                                      */
/*---------------------------------------------------------------------*/


/*---------------------------------------------------------------------*/
/* Allocates msh memory                                                */
/*---------------------------------------------------------------------*/
static msh_t *msh_alloc(int ndims, int npoint)
{
  msh_t *msh = (msh_t *)malloc(sizeof(msh_t));
  memset(msh, 0, sizeof(msh_t));

  msh->ndims = ndims;
  msh->npoint = npoint;
  msh->coords = d_alloc_2d(npoint, ndims);
  msh->edges = i_alloc_2d(npoint, ndims);
  msh->flag = i_alloc_1d(npoint);
  memset(msh->flag, 0, npoint * sizeof(int));
  return msh;
}

/* END msh_alloc()                                                     */
/*---------------------------------------------------------------------*/


/*---------------------------------------------------------------------*/
/* Allocates and initialises jig structure                             */
/*---------------------------------------------------------------------*/
static jig_t *jig_alloc(char *infile, char *hfile, double hmin, double hmax)
{
  int i;
  char buf[MAXSTRLEN];
  jig_t *jig = (jig_t *)malloc(sizeof(jig_t));
  memset(jig, 0, sizeof(jig_t));

  for (i = 0; i < strlen(infile); i++) {
    if (infile[i] == '.')
      break;
    else
      buf[i] = infile[i];
  }
  buf[i] = '\0';

  strcpy(jig->geom_file, infile);
  sprintf(jig->mesh_file, "%s_out.msh", buf);

  jig->geom_seed = 32;
  jig->geom_feat = 0;
  jig->geom_eta1 = 45;
  jig->geom_eta2 = 45;
  strcpy(jig->hfun_kern, "delaunay");
  sprintf(jig->hfun_file, "%c", '\0');
  if (strlen(hfile)) strcpy(jig->hfun_file,hfile);
  if (hmin && hmax) {
    strcpy(jig->hfun_scal, "absolute");
    jig->hfun_hmax = hmax;
    jig->hfun_hmin = hmin;
  } else {
    strcpy(jig->hfun_scal, "relative");
    jig->hfun_hmax = 0.02;
    jig->hfun_hmin = 0.0;
  }
  jig->hfun_grad = 0.25;
  jig->mesh_dims = 3;
  strcpy(jig->mesh_kern, "delfront");
  jig->mesh_iter = 1000;
  jig->mesh_top1 = 0;
  jig->mesh_top2 = 0;
  jig->mesh_rad2 = 1.05;
  jig->mesh_rad3 = 2.05;
  jig->mesh_eps1 = 0.33;
  jig->mesh_eps2 = 0.33;
  jig->mesh_vol3 = 0.0;
  jig->verbosity = 0;

  return jig;
}

static void jig_free(jig_t *jig)
{
  free((jig_t *)jig);
}

/* END jig_alloc()                                                     */
/*---------------------------------------------------------------------*/


static void order (int *a, int *b, int *c1, int *c2)
{
  int val;
  if (*a > *b) {
    val = *a;
    *a = *b;
    *b = val;
    val = *c1;
    *c1 = *c2;
    *c2 = val;
  }
}

static int find_index(double slon, double slat, double *lon, double *lat, int nsl, double *d)
{
  double dist, d1, d2;
  double dmin = HUGE;
  int i;
  int ret = -1;

  for (i = 0; i < nsl; i++) {
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

/*---------------------------------------------------------------------*/
/* Finds the index in a segment for a link                             */
/*---------------------------------------------------------------------*/
static
void find_link_index(int ns,       /* Number of segments               */
		     int *nsl,     /* Segment sizes                    */
		     double **lat, /* Coastline latitude               */
		     double **lon, /* Coastline longitide              */
		     link_t *links,/* Links                            */
		     int dir       /* Direction to cycle major segment */
		     )
{
  
  double dist;
  double dist1 = HUGE;
  double dist2 = HUGE;
  int j, n, m;

  for (n = 0; n < ns; n++) {      
    /* Get the segment number and index of the start link coordinate   */
    /* The minimum distance of all segments is used to specify this    */
    /* index.                                                          */
    j = find_index(links->llon[0], links->llat[0], lon[n], lat[n], nsl[n], &dist);
    if (j >= 0 && dist < dist1) {
      dist1 = dist;
      links->slink[0] = n;
      links->ilink[0] = j;
    }
    j = find_index(links->llon[1], links->llat[1], lon[n], lat[n], nsl[n], &dist);
    if (j >= 0 && dist < dist2) {
      dist2 = dist;
      links->slink[1] = n;
      links->ilink[1] = j;
      if (links->flag & (LINK_M|LINK_P))
	links->ilink[2] = find_index(links->llon[2], links->llat[2], lon[n], lat[n], nsl[n], NULL);
    }
  }
  /*
  printf("link%d start(s=%d i=%d) end(s=%d i=%d)\n",links->id, 
	 links->slink[0], links->ilink[0],
	 links->slink[1], links->ilink[1]);
  */

  if (links->flag & LINK_N && links->slink[0] != links->slink[1]) {
    printf("LINK%d requires a direction.\n", links->id);
    quit("coastmesh error.\n");
  }
  /* Swap if required                                                  */
  if (links->flag & LINK_N) {
    if (links->dir == 0) links->dir = dir;
    if (dir == 1 && links->ilink[0] > links->ilink[1]) {
      m = links->ilink[0];
      links->ilink[0] = links->ilink[1];
      links->ilink[1] = m;
    }
    if (dir == -1 && links->ilink[0] < links->ilink[1]) {
      m = links->ilink[0];
      links->ilink[0] = links->ilink[1];
      links->ilink[1] = m;
    }
    if (links->nseg) {
      double *x1, *y1;
      double s1 = links->llon[0] - links->segx[0];
      double s2 = links->llat[0] - links->segy[0];
      double e1 = links->llon[0] - links->segx[links->nseg - 1];
      double e2 = links->llat[0] - links->segy[links->nseg - 1];
      double d1 = sqrt(s1 * s1 + s2 * s2);
      double d2 = sqrt(e1 * e1 + e2 * e2);
      if (d2 < d1) { /* End of path closest to start coordinate - swap */
	double *x1 = d_alloc_1d(links->nseg);
	double *y1 = d_alloc_1d(links->nseg);
	memcpy(x1, links->segx, links->nseg * sizeof(double));
	memcpy(y1, links->segy, links->nseg * sizeof(double));
	for (m = 0; m < links->nseg; m++) {
	  links->segx[m] = x1[links->nseg - m - 1];
	  links->segy[m] = y1[links->nseg - m - 1];
	}
	d_free_1d(x1);
	d_free_1d(y1);
      }
    }
  }

  /* Populate with coordinates                                         */
  j = links->ilink[0];
  n = links->slink[0];
  links->alon[0] = lon[n][j];
  links->alat[0] = lat[n][j];
  j = links->ilink[1];
  n = links->slink[1];
  links->alon[1] = lon[n][j];
  links->alat[1] = lat[n][j];
  if (links->flag & (LINK_M|LINK_P)) {
    j = links->ilink[2];
    links->alon[2] = lon[n][j];
    links->alat[2] = lat[n][j];
  }
}

/* END find_link_index()                                               */
/*---------------------------------------------------------------------*/


/*---------------------------------------------------------------------*/
/* Creates or counts a list of points from a start index to the next   */
/* start index in the link list. Returns the end link index and        */
/* segment number. This routine can be called multiple times to create */
/* a continuous segment where different segments are linked, or        */
/* points in the same segment are linked.                              */
/*---------------------------------------------------------------------*/
static
int find_link(int en,          /* End segment number                   */
	      int eid,         /* End index number                     */
	      int *n,          /* Current segment                      */
	      int *si,         /* Current index                        */
	      int *dir,        /* Current direction                    */
	      int *size,       /* Size of major segment                */
	      int *nsl,        /* Segment sizes                        */
	      int nlink,       /* Number of links                      */
	      link_t *link,    /* Links                                */
	      double **lat,    /* Segment latitude                     */
	      double **lon,    /* Segment longitude                    */
	      double *nlat,    /* Major segment latitude               */
	      double *nlon,    /* Major segment longitude              */
	      int *flag        /* Status flag of the point             */
	      )
{
  int i, ii, m, ns;            /* Counters                             */
  int found = 0;               /* Flag to exit when next link found    */
  double *rlat, *rlon;         /* Updated coordinates of major segment */
  int psize = *size;           /* Size of major segment before update  */
  int psi = *si;               /* Start index before update            */
  int pn = *n;                 /* Segment before the update            */
  int pdir = *dir;             /* Traverse direction before update     */
  int verbose = 0;

  if (verbose)
    printf("\nstart i=%d n=%d dir=%d size=%d\n",*si, *n, *dir, *size);
  i = *si; ii = 0;
  while (!found) {
    for (m = 0; m < nlink; m++) {   /* Loop over all links             */
      if (*n == en && i == eid) {
	*n = en;
	*si = eid;
	found = 3;
	break;
      }
      if (*n == link[m].slink[0] && i == link[m].ilink[0]) {
	/* Find the direction to traverse in after the link            */
	if (link[m].slink[1] != *n)
	  *dir = find_dir(&link[m], nsl, link[m].slink[1]);
	if (verbose) printf("found l=%d i=%d s=%d ei=%d es=%d dir=%d\n",m,i,*n,link[m].ilink[1],link[m].slink[1],*dir);
	*n = link[m].slink[1];     /* Updated segment number          */
	*si = link[m].ilink[1];    /* Updated start index             */
	found = 2;
	break;
      }
    }
    i += *dir;
    if (*dir == 1 && i >= nsl[*n]) i = 0;
    if (*dir == -1 && i < 0) i = nsl[*n] - 1;
    ii++;
    if (!found && ii >= nsl[*n]) found = 1;
  }
  if (found == 1) {
    printf("Couldn't find any links in segment %d (%f %f)\n", *n, lon[*n][0], lat[*n][0]);
    quit("coastmesh error.\n");
  }
  *size += ii;

  if (verbose)
    printf("end i=%d n=%d dir=%d size=%d link=%d\n",*si, *n, *dir, *size, m);
  if (nlat == NULL && nlon == NULL) {
    if (found == 2 && link[m].nseg) *size += link[m].nseg;
    return(m);
  }

  /* Fill the coordinate arrays with updated locations                 */
  i = psi;
  for (ii = psize; ii < *size; ii++) {
    ns = (found == 2) ? link[m].slink[0] : en;
    nlat[ii] = lat[ns][i];
    nlon[ii] = lon[ns][i];
    flag[ii] = S_COAST;
    i += pdir;
    if (pdir == 1 && i >= nsl[pn]) i = 0;
    if (pdir == -1 && i < 0) i = nsl[pn] - 1;
  }
  if (link[m].nseg) {
    *size += link[m].nseg;
    for (i = 0; i < link[m].nseg; i++) {
      nlon[ii] = link[m].segx[i];
      nlat[ii] = link[m].segy[i];
      flag[ii] = S_LINK;
     ii++;
    }
  }

  return(m);
}

/* END find_link()                                                     */
/*---------------------------------------------------------------------*/


/*---------------------------------------------------------------------*/
/* Finds the direction to traverse a segment given link information    */
/*---------------------------------------------------------------------*/
static
int find_dir(link_t *links, int *nsl, int n)
{
  int ii, i1, i2, i3, dir;

  if (links->flag & LINK_M) {

    /* Find the direction to cycle the segment; dir=1 means we     */
    /* cycle in the +i direction to reach the mid index in the     */
    /* minimum number of steps, dir=-1 means cycle in the -i       */
    /* direction.                                                  */
    /* In the increasing i direction.                              */
    i1 = links->ilink[1];
    i2 = 0;
    for (ii = 0; ii < nsl[n]; ii++) {
      if (i1 == links->ilink[2]) break;
      i2++;
      i1 = (i1 == nsl[n]-1) ? 0 : i1+1;
    }
    /* In the decreasing i direction.                              */
    i1 = links->ilink[1];
    i3 = 0;
    for (ii = 0; ii < nsl[n]; ii++) {
      if (i1 == links->ilink[2]) break;
      i3++;
      i1 = (i1 == 0) ? nsl[n]-1 : i1-1;
    }
    dir = 1;
    if (i3 < i2) dir = -1;
  } else {
    dir = links->dir;
  }
  return(dir);
}

/* END find_dir()                                                      */
/*---------------------------------------------------------------------*/


/*---------------------------------------------------------------------*/
/* Connects a link to create a new segment                             */
/*---------------------------------------------------------------------*/
static
void connect_link(int ns,          /* Number of segments               */
		  int nlink,       /* Number of links                  */
		  int *nsl,        /* Segment sizes                    */
		  int *nlm,        /* Maximum segment size             */
		  double **lat,    /* Coastline latitude               */
		  double **lon,    /* Coastline longitide              */
		  int il,          /* Link to process                  */
		  link_t *link     /* Links                            */ 
		  )
{
  int n, n1, i, ii, i1, i2, i3, j;
  int dir;
  double **rlat, **rlon;
  link_t *links;

  /* Find the new maximum segment size                                 */
  for (i = 0; i < nlink; i++) {
    links = &link[i];
    if ((links->slink[0] != links->slink[1]) && nsl[links->slink[0]] + nsl[links->slink[1]] > *nlm)
      *nlm = nsl[links->slink[0]] + nsl[links->slink[1]] + links->nseg;
  }
  /* Allocate and copy                                                 */
  rlon = d_alloc_2d(*nlm, ns);
  rlat = d_alloc_2d(*nlm, ns);
  for (n = 0; n < ns; n++) {
    for (i = 0; i < nsl[n]; i++) {
      rlon[n][i] = lon[n][i];
      rlat[n][i] = lat[n][i];
    }
  }

  links = &link[il];
  n = links->slink[0];
  j = 0;
  for (i = 0; i < nsl[n]; i++) {
    rlon[n][j] = lon[n][i];
    rlat[n][j] = lat[n][i];
    j++;
    if (i == links->ilink[0]) {
      /* Add in the link path if supplied                              */
      for (ii = 0; ii < links->nseg; ii++) {
	rlon[n][j] = links->segx[ii];
	rlat[n][j] = links->segy[ii];
	j++;
      }
      n1 = links->slink[1];
      if (links->flag & LINK_M) {
	/* Find the direction to cycle the segment; dir=1 means we     */
	/* cycle in the +i direction to reach the mid index in the     */
	/* minimum number of steps, dir=-1 means cycle in the -i       */
	/* direction.                                                  */
	/* In the increasing i direction.                              */
	i1 = links->ilink[1];
	i2 = 0;
	for (ii = 0; ii < nsl[n1]; ii++) {
	  if (i1 == links->ilink[2]) break;
	  i2++;
	  i1 = (i1 == nsl[n1]-1) ? 0 : i1+1;
	}
	/* In the decreasing i direction.                              */
	i1 = links->ilink[1];
	i3 = 0;
	for (ii = 0; ii < nsl[n1]; ii++) {
	  if (i1 == links->ilink[2]) break;
	  i3++;
	  i1 = (i1 == 0) ? nsl[n1]-1 : i1-1;
	}
	dir = 1;
	if (i3 < i2) dir = -1;
      } else
	dir = links->dir;
      /* Cycle the end segment and add coordinates to rlat/rlon        */
      i1 = links->ilink[1];
      if (links->flag & LINK_N) i = i1;
      printf("start %f %f\n",lon[n1][i1],lat[n1][i1]);
      for (ii = 0; ii < nsl[n1]; ii++) {
	rlon[n][j] = lon[n1][i1];
	rlat[n][j] = lat[n1][i1];
	j++;
	if (i1 == links->ilink[1] - dir) break;
	if (dir == 1)
	  i1 = (i1 == nsl[n1]-1) ? 0 : i1+1;
	else
	  i1 = (i1 == 0) ? nsl[n1]-1 : i1-1;
      }
      if (n1 != n) nsl[n1] = 0;
    }
  }
  nsl[n] = j;

  for (n = 0; n < ns; n++) {
    for (i = 0; i < nsl[n]; i++) {
      lon[n][i] = rlon[n][i];
      lat[n][i] = rlat[n][i];
    }
  }
  d_free_2d(rlon);
  d_free_2d(rlat);
}

/* END connect_link()                                                  */
/*---------------------------------------------------------------------*/

