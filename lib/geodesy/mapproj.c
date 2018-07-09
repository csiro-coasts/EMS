/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/geodesy/mapproj.c
 *
 *  \brief Map projection routines.
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: mapproj.c 5831 2018-06-26 23:48:06Z riz008 $
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "ems.h"
#if defined(HAVE_USGS_PROJ)
#include "projects.h"
#endif

#define HALF_PI (M_PI/2.0)
#define QUARTER_PI (M_PI/4.0)
#define DEG2RAD(d) (M_PI*(d)/180.0)
#define RAD2DEG(r) (180.0*(r)/M_PI)

/* Common mandatory projections arguments. */
#define PROJECTION_NAME(n, a) get_argument_value("proj", n, a, 1)
#define FALSE_EASTING(n, a) atof(get_argument_value("x_0", n, a, 1))
#define FALSE_NORTHING(n, a) atof(get_argument_value("y_0", n, a, 1))
#define CENTRAL_MERIDIAN(n, a) atof(get_argument_value("lon_0", n, a, 1))
#define CENTRAL_LATITUDE(n, a) atof(get_argument_value("lat_0", n, a, 1))
#define FIRST_LATITUDE(n, a) atof(get_argument_value("lat_1", n, a, 1))
#define SECOND_LATITUDE(n, a) atof(get_argument_value("lat_2", n, a, 1))
#define SCALE_FACTOR(n, a) atof(get_argument_value("k_0", n, a, 1))
#define UTM_ZONE(n, a) atof(get_argument_value("zone", n, a, 1))

/* Prototypes */
static int find_argument(char *tag, int nargs, char *args[]);
static char *get_argument_value(char *tag, int nargs, char *arglist[],
                                int req);
static char **add_argument(char *arg, int nargs, char *args[],
                           int *new_nargs);
static void subst_argument_value(char *arg, char *value, int nargs,
                                 char *args[]);

static void *merc_init(map_proj_t *mp, int nargs, char *args[]);
static void merc_cleanup(void *data);
static void merc_forward(map_proj_t *mp, double lat, double lon,
                         double *east, double *north);
static void merc_inverse(map_proj_t *mp, double east, double north,
                         double *lat, double *lon);

static void *tcc_init(map_proj_t *mp, int nargs, char *args[]);
static void tcc_cleanup(void *data);
static void tcc_forward(map_proj_t *mp, double lat, double lon,
                        double *east, double *north);
static void tcc_inverse(map_proj_t *mp, double east, double north,
                        double *lat, double *lon);
static double tcc_footpoint_latitude(map_proj_t *mp, double M);

static void *lcc_init(map_proj_t *mp, int nargs, char *args[]);
static void lcc_cleanup(void *data);
static void lcc_forward(map_proj_t *mp, double lat, double lon,
                        double *east, double *north);
static void lcc_inverse(map_proj_t *mp, double east, double north,
                        double *lat, double *lon);

static void *serov_init(map_proj_t *mp, int nargs, char *args[]);
static void serov_cleanup(void *data);
static void serov_forward(map_proj_t *mp, double lat, double lon,
                          double *east, double *north);
static void serov_inverse(map_proj_t *mp, double east, double north,
                          double *lat, double *lon);

#if defined(HAVE_USGS_PROJ)
static void *usgs_init(map_proj_t *mp, int nargs, char *args[]);
static void usgs_cleanup(void *data);
static void usgs_forward(map_proj_t *mp, double lat, double lon,
                         double *east, double *north);
static void usgs_inverse(map_proj_t *mp, double east, double north,
                         double *lat, double *lon);
#endif

/* Known ellpsoids. */
typedef struct {
  char *id;                     /* ellipse keyword name */
  double major;                 /* major axis */
  double flat;                  /* flattening */
  char *name;                   /* comments */
} ellipsoid_t;

ellipsoid_t ellipsoids[] = {
  {"merit", 6378137.0, 1 / 298.257, "MERIT 1983"},
  {"grs80", 6378137.0, 1 / 298.257222, "GRS 1980(IUGG, 1980)"},
  {"iau76", 6378140.0, 1 / 298.257, "IAU 1976"},
  {"airy", 6377563.396, 1 / 299.324975315, "Airy 1830"},
  {"mod_airy", 6377340.189, 1 / 299.348780462, "Modified Airy"},
  {"aust_ntl", 6378160.0, 1 / 298.25, "Australian Natl, S. Amer., IAU 64"},
  {"grs67", 6378160.0, 1 / 247.247167, "GRS 67(IUGG 1967)"},
  {"bessel", 6377397.155, 1 / 299.1528128, "Bessel 1841"},
  {"bess_nam", 6377483.865, 1 / 299.1528128, "Bessel 1841 (Namibia)"},
  {"clrk66", 6378206.4, 1 / 294.98, "Clarke 1866"},
  {"clark66", 6378206.4, 1 / 294.98, "Clarke 1866"},
  {"clrk80", 6378249.145, 1 / 293.4663, "Clarke 1880 mod."},
  {"everest", 6377276.3452, 1 / 300.80, "Everest 1830"},
  {"hough", 6378270.0, 1 / 297.0, "Hough"},
  {"intl", 6378388.0, 1 / 297.0, "International 1909 (Hayford)"},
  {"krass", 6378245.0, 1 / 298.3, "Krassovsky, 1942"},
  {"mercury", 6378166.0, 1 / 298.3, "Mercury 1960"},
  {"mod_ever", 6377304.063, 1 / 300.8, "Modified Everest"},
  {"mod_merc", 6378150.0, 1 / 298.3, "Modified Merc 1968"},
  {"new_intl", 6378157.5, 1 / 298.24961539, "New International 1967"},
  {"seasia", 6378155.0, 1 / 298.3, "Southeast Asia"},
  {"walbeck", 6376896.0, 1 / 302.78, "Walbeck"},
  {"wgs66", 6378145.0, 1 / 298.25, "WGS 66"},
  {"wgs72", 6378135.0, 1 / 298.26, "WGS 72"},
  {"wgs84", 6378137.0, 1 / 298.257223563, "WGS 84"},
  {"agd66", 6378160.0, 1 / 298.25, "Same as aust_ntl"},
  {"agd84", 6378160.0, 1 / 298.25, "Same as aust_ntl"},
  {"gda94", 6378137.0, 1 / 298.25722101, "New Aust. ellip"},
  {"sphere", 6370997.0, 0.0, "Sphere of 6370997 m"},
  {(char *)0, 0, 0, (char *)0}
};


/** Create a map projection from a string of arguments.
  * All arguments are string valued and conform to the syntax
  * param[=value].
  *
  * param - A tag name associated with a parameter.
  * value - Optional parameter value.
  *
  * The projection is defined by the mandatory parameter
  * tag 'proj=\<prjname\>'.
  *
  * Global parameter (applicable to all projections) are as follows:
\verbatim
ellps  - Ellipsoid name.
es     - Eccentricity.
a      - Major ellipsoid axis radius.
b      - Minor ellipsoid axis radius.
rf     - Reverse flattening.
x_0    - explicit false easting. In UTM it is implicit.
y_0    - explicit false northing.
\end{verbatim}
  * Projections implicitly supported (and mandatory arguments) are:
\begin{verbatim}
amg     - Australian Map Grid
  zone  - UTM zone (1-60).

lcc     - Lambert Conformal Conic
  lon_0 - Central meridian.
  lat_0 - Central latitude.
  lat_1 - First standard parallel latitude.
  lat_2 - Second standard parallel latitude.

merc    - Mercator
  lon_0 - Central meridian.

mga     - Map Grid of Australia
  zone  - UTM zone (1-60).

serov   - CMR Remote sensing Serov projection.

tcc     - Transverse Central Cylindrical
  lon_0 - Central meridian.
  k_0   - Scale factor.

utm     - Universal Transverse Mercator
  zone  - UTM zone (1-60).
  south - Enabled if for southern hemisphere.
  north - Enabled if for northern hemisphere (default).
\endverbatim
  *
  * Projections supported by the USGS PROJ 4 library (if enabled) are:
\verbatim
aea     - Albers Egual Area
aeqd    - Azimuthal equidistant
alsk    - Alaska Mod.-Stereographics
apian   - Apian Globular
bipc    - Bipolar Conic
bonne   - Bonne
cass    - Cassini
cc      - Central Cylindrical
cea     - Cylindrical Equal Area
collg   - Collignon
eck1    - Eckert I
eck2    - Eckert II
eck3    - Eckert III
eck4    - Eckert IV
eck5    - Eckert V
eck6    - Eckert VI
eqc     - Equidistant Cylindrical
eqdc    - Equidistant Conic
gall    - Gall (Stereographic)
gnom    - Gnomonic
gs50    - 50 State U.S. Mod.-Stereographic
gs48    - 48 State U.S. Mod.-Stereographic
hataea  - Hatano Asymmetrical Equal Area
labrd   - Laborde
laea    - Lambert Azimuthal Equal Area
leac    - Lambert Equal Area Conic
lee_os  - Lee Oblate Stereographics Pacific
loxim   - Loximuthal
lsat    - LANDSAT Space Oblique Mercator
mbtfpp  - McBryde-Thomas Flat-Polar Parabolic
mbtfps  - McBryde-Thomas Flat-Polar Sinusoidal
mbtfpq  - McBryde-Thomas Flat-Polar Quartic
mill    - Miller
mill_os - Miller Oblate Stereographics Eur-Africa
moll    - Mollweides
nicol   - Nicolosi Globular
nsper   - General Vertical Persepective
nzmg    - New Zealand Map Grid
ocea    - Oblique Cylindrical Equal Area
omerc   - Oblique Mercator
ortho   - Orthographic
parab   - Caster Parabolic
poly_t    - Polyconic (American)
putp2   - Putnins P2'
putp5   - Putnins P5
quau    - Quartic Authalic
robin   - Robinson
sinu    - Sinusoidal
stere   - Stereographic
tcea    - Transverse Cylindrical Equal Area
tpers   - Tilted perspective
ups     - Universal Polar Stereographic
vandg   - Van der Grinten
wink1   - Winkel I
\endverbatim
  *
  * Please consult the USGS manuals for a full description of
  * relavent arguments for each projection.
  *
  * Supported ellipsoid names are:
\verbatim
merit     a=6378137.0 rf=298.257        - MERIT 1983.
grs80     a=6378137.0 rf=298.257222     - GRS 1980(IUGG, 1980).
iau76     a=6378140.0 rf=298.257        - IAU 1976.
airy      a=6377563.396 b=6356256.910   - Airy 1830.
mod_airy  a=6377340.189 b=6356036.143   - Modified Airy.
aust_ntl  a=6378160.0 rf=298.25         - Australian Natl, S. Amer., IAU 64.
grs67     a=6378160.0 rf=247.247167     - GRS 67(IUGG 1967).
bessel    a=6377397.155 rf=299.1528128  - Bessel 1841.
bess_nam  a=6377483.865 rf=299.1528128  - Bessel 1841 (Namibia).
clrk66    a=6378206.4 b=6356583.8       - Clarke 1866.
clark66   a=6378206.4 b=6356583.8       - Clarke 1866.
clrk80    a=6378249.145 rf=293.4663     - Clarke 1880 mod..
everest   a=6377276.3452 b=6356075.4133 - Everest 1830.
hough     a=6378270.0 b=6356794.343479  - Hough.
intl      a=6378388.0 rf=297.           - International 1909 (Hayford).
krass     a=6378245.0 rf=298.3          - Krassovsky, 1942.
mercury   a=6378166.0 b=6356784.283666  - Mercury 1960.
mod_ever  a=6377304.063 b=6356103.039   - Modified Everest.
mod_merc  a=6378150.0 b=6356768.337303  - Modified Merc 1968.
new_intl  a=6378157.5 b=6356772.2       - New International 1967.
SEasia    a=6378155.0 b=6356773.3205    - Southeast Asia.
walbeck   a=6376896.0 b=6355834.8467    - Walbeck.
wgs66     a=6378145.0 b=6356759.769356  - WGS 66.
wgs72     a=6378135.0 b=6356750.519915  - WGS 72.
wgs84     a=6378137.0 rf=298.257223563  - WGS 84.
agd66     a=6378160.0 rf=298.25         - Sames as aust_ntl.
agd84     a=6378160.0 rf=298.25         - Sames as aust_ntl.
gda94     a=6378137.0 rf=298.25722101   - New Aust. ellip.
sphere    a=6370997.0 es=0.0            - Sphere of 6370997 m.
\endverbatim
  *
  * warn, quit maybe called if invalid arguments are passed.
  * 
  * @param nargs Number of parameter arguements.
  * @param args Array of parameter arguments as strings. The
  * arguments conform to the proj conventions.
  * @return pointer to map projection information.
  */
map_proj_t *mp_init(int nargs, char *args[])
{
  int new_nargs;
  char **new_args;
  map_proj_t *mp = NULL;
  char *projname = PROJECTION_NAME(nargs, args);

/* Check for those projections that can be handled outside of 'proj'. */

  /* Australian Map Grid */
  if (strcasecmp(projname, "amg") == 0) {
    new_args = add_argument("ellps=agd66", nargs, args, &new_nargs);
    new_args = add_argument("south", new_nargs, new_args, &new_nargs);
    subst_argument_value("proj", "utm", new_nargs, new_args);
    mp = mp_init(new_nargs, new_args);
  }

  /* Map Grid of Australia */
  else if (strcasecmp(projname, "mga") == 0) {
    new_args = add_argument("ellps=gda94", nargs, args, &new_nargs);
    new_args = add_argument("south", new_nargs, new_args, &new_nargs);
    subst_argument_value("proj", "utm", new_nargs, new_args);
    mp = mp_init(new_nargs, new_args);
  }

  /* Universal Transverse Mercator */
  else if (strcasecmp(projname, "utm") == 0) {
    char buf[256];
    int is_north = (find_argument("south", nargs, args) < 0);

    new_args = add_argument("x_0=500000.0", nargs, args, &new_nargs);
    if (is_north)
      new_args = add_argument("y_0=0.0", new_nargs, new_args, &new_nargs);
    else
      new_args =
        add_argument("y_0=10000000.0", new_nargs, new_args, &new_nargs);

    sprintf(buf, "lon_0=%f", UTM_ZONE(new_nargs, new_args) * 6 - 183);
    new_args = add_argument(buf, new_nargs, new_args, &new_nargs);

    new_args = add_argument("k_0=0.9996", new_nargs, new_args, &new_nargs);

    subst_argument_value("proj", "tcc", new_nargs, new_args);
    mp = mp_init(new_nargs, new_args);
  }

  /* Mercator, Transverse Mercator, Lambert Conformal Conic and CMR remote 
     sensing Serov projection. */
  else if ((strcasecmp(projname, "merc") == 0)
           || (strcasecmp(projname, "tcc") == 0)
           || (strcasecmp(projname, "lcc") == 0)
           || (strcasecmp(projname, "serov") == 0)) {
    int got_ellip_data = 0;
    char *ename;
    char *falsestr;
    mp = (map_proj_t *)malloc(sizeof(map_proj_t));

    /* Set the default ellisoid and false eastings/northings, etc. */
    mp->ellip_major = 6378206.4;  /* Default Clark 1866. */
    mp->ellip_flat = 1 / 294.98;
    mp->falsex = 0.0;
    mp->falsey = 0.0;

    /* Extract the ellipsoid data. */
    ename = get_argument_value("ellps", nargs, args, 0);
    if (ename != NULL) {
      ellipsoid_t *e = &ellipsoids[0];
      while (e->id != NULL) {
        if (strcasecmp(e->id, ename) == 0) {
          mp->ellip_major = e->major;
          mp->ellip_flat = e->flat;
          got_ellip_data = 1;
          break;
        }
        *e++;
      }
      if (!got_ellip_data)
        quit("Ellipsoid '%s' is unknown.\n", ename);
    }

    if (!got_ellip_data) {
      char *major = get_argument_value("a", nargs, args, 0);
      char *minor = get_argument_value("b", nargs, args, 0);
      char *rf = get_argument_value("rf", nargs, args, 0);
      char *es = get_argument_value("es", nargs, args, 0);

      /* Get major axis */
      if (major != NULL)
        mp->ellip_major = atof(major);

      /* Get minor axis */
      if (minor != NULL) {
        double b = atof(minor);
        mp->ellip_flat = 1 - b / mp->ellip_major;
      }

      /* Get reverse flattening. */
      if (rf != NULL)
        mp->ellip_flat = 1 / atof(rf);

      /* Get eccentricity */
      if (es != NULL) {
        double ecc = atof(es);
        mp->ellip_flat = 1 / (1 - sqrt(1 - ecc * ecc));
      }
    }

    /* Extract the false eastings/northings */
    falsestr = get_argument_value("x_0", nargs, args, 0);
    if (falsestr != NULL)
      mp->falsex = atof(falsestr);
    falsestr = get_argument_value("y_0", nargs, args, 0);
    if (falsestr != NULL)
      mp->falsey = atof(falsestr);

    /* Setup the Mercator projection. */
    if (strcasecmp(projname, "merc") == 0) {
      mp->init = merc_init;
      mp->free = merc_cleanup;
      mp->forward = merc_forward;
      mp->inverse = merc_inverse;
    }

    /* Setup the Transverse Mercator projection. */
    else if (strcasecmp(projname, "tcc") == 0) {
      mp->init = tcc_init;
      mp->free = tcc_cleanup;
      mp->forward = tcc_forward;
      mp->inverse = tcc_inverse;
    }

    /* Setup the Lambert Conformal Conic projection. */
    else if (strcasecmp(projname, "lcc") == 0) {
      mp->init = lcc_init;
      mp->free = lcc_cleanup;
      mp->forward = lcc_forward;
      mp->inverse = lcc_inverse;
    }

    /* Setup the CMR remote sensing Serov projection. */
    else if (strcasecmp(projname, "serov") == 0) {
      mp->init = serov_init;
      mp->free = serov_cleanup;
      mp->forward = serov_forward;
      mp->inverse = serov_inverse;
    }

    mp->private_data = mp->init(mp, nargs, args);
  }
#if defined(HAVE_USGS_PROJ)

  else {
    mp = (map_proj_t *)malloc(sizeof(map_proj_t));
    mp->init = proj_init;
    mp->free = proj_cleanup;
    mp->forward = proj_forward;
    mp->inverse = proj_inverse;
    mp->private_data = mp->init(mp, nargs, args);
  }

#else

  else
    quit("The projection '%s' is not supported.", projname);

#endif

  return mp;
}


/** Deallocate any memory associated with a map projection.
  *
  * @param mp pointer to map projection.
  */
void mp_cleanup(map_proj_t *mp)
{
  if (mp != NULL) {
    if (mp->free != NULL)
      mp->free(mp->private_data);
    free(mp);
  }
}


/** Convert from latitude/longitude to projected X/Y for
  * the projection specified.
  *
  * @param mp Pointer to map projection.
  * @param lat Latitude in degrees (-90 to 90).
  * @param lon Longitude in degrees (-180 to 180).
  * @param east Pointer to easting value.
  * @param north Pointer to north value.
  */
void mp_forward(map_proj_t *mp, double lat, double lon,
                double *east, double *north)
{
  mp->forward(mp, lat, lon, east, north);
}


/** Convert from projected X/Y to latitude/longitude for
  * the projection specified.
  *
  * @param mp Pointer to map projection.
  * @param east Easting value.
  * @param north North value.
  * @param lat Pointer to latitude in degrees (-90 to 90).
  * @param lon Pointer to longitude in degrees (-180 to 180).
  */
void mp_inverse(map_proj_t *mp, double east, double north,
                double *lat, double *lon)
{
  mp->inverse(mp, east, north, lat, lon);
}


/* Search for the specified argument within the list, and return
 * the value as a string. NULL if it cannot be found.
 */
static char *get_argument_value(char *tag, int nargs, char *args[],
                                int req)
{
  int index = find_argument(tag, nargs, args);

  if (index >= 0) {
    char *arg = args[index];
    while (*arg) {
      if (*arg++ == '=')
        return arg;
    }
  }

  if (req)
    quit("Unable to locate the Map Projection argument '%s'.\n", tag);

  return NULL;
}

/* Locate an argument in the list by name.
 * @return index of tag in list or -1.
 */
static int find_argument(char *tag, int nargs, char *args[])
{
  int i;

  for (i = 0; i < nargs; ++i) {
    if (strncasecmp(tag, args[i], strlen(tag)) == 0)
      return i;
  }

  return -1;
}

/* Add a new argument to the list.
 */
static char **add_argument(char *arg, int nargs, char *args[],
                           int *new_nargs)
{
  char **new_args;
  int i;

  *new_nargs = nargs + 1;
  new_args = (char **)malloc(sizeof(char *) * (*new_nargs));
  new_args[0] = (char *)malloc(sizeof(char *) * strlen(arg) + 1);
  strcpy(new_args[0], arg);

  for (i = 0; i < nargs; ++i)
    new_args[i + 1] = args[i];

  return new_args;
}


/* Subsitute an argument.
 */
static void subst_argument_value(char *arg, char *value,
                                 int nargs, char *args[])
{
  int index = find_argument(arg, nargs, args);

  if (index >= 0)
    sprintf(args[index], "%s=%s", arg, value);
}


/*
 * Mercator Map Projection.
 *
 * The Mercator map projection is a cylindrical and conformal map
 * projection. It has the properties that all meridians are equally
 * spaced straight lines, parallels are unequally spaced (closer at
 * the equator), and Rhumb lines are down as strait lines.
 *
 * This implementation was taken from the USGS Map Projections
 * manual pg. 38-47.
 */
typedef struct {
  double cm;                    /* Centeral meridian in radians. */
  double e2;                    /* Eccentricity squared. */
  double e4;                    /* Eccentricity ** 4. */
  double e6;                    /* Eccentricity ** 6. */
} mercator_t;

/* Mercator map projection constructor.
 *
 * @param mp Pointer to map projection information.
 * @param nargs Number of parameter arguements.
 * @param args Array of parameter arguments as strings. The
 * arguments conform to the proj conventions.
 * @return Pointer to data private to the mercator projection.
 */
static void *merc_init(map_proj_t *mp, int nargs, char *args[])
{
  mercator_t *data = (mercator_t *)malloc(sizeof(mercator_t));
  data->cm = DEG2RAD(CENTRAL_MERIDIAN(nargs, args));
  data->e2 = (2.0 * mp->ellip_flat) - (mp->ellip_flat * mp->ellip_flat);
  data->e4 = data->e2 * data->e2;
  data->e6 = data->e2 * data->e2 * data->e2;

  return data;
}

/* Free the private memory associated with a mercator transformation.
 */
static void merc_cleanup(void *data)
{
  if (data != NULL)
    free(data);
}

/* Method to convert latitudes and longitudes
 * to a projected Mercator eastings and northings.
 *
 * @param mp    Map projection information.
 * @param lat   Latitude in degrees.
 * @param lon   Longitude in degrees.
 * @param east  Pointer to Easting in metres
 * @param north Pointer to Northing in metres
 */
static void merc_forward(map_proj_t *mp, double lat, double lon,
                         double *east, double *north)
{
  mercator_t *data = (mercator_t *)mp->private_data;
  double f1 = -(data->e2 + data->e4 / 4.0 + data->e6 / 8.0);
  double f2 = data->e4 / 12.0 + data->e6 / 16.0;
  double rlat = DEG2RAD(lat);
  double rlon = DEG2RAD(lon);
  double x = mp->ellip_major * (rlon - data->cm);
  double y = log(tan(QUARTER_PI + rlat / 2.0));
  y += f1 * sin(rlat) + f2 * sin(3.0 * rlat);
  y *= mp->ellip_major;
  *east = x + mp->falsex;
  *north = y + mp->falsey;
}

/* Method to convert projected Mercator eastings and northings
 * to latitudes and longitudes.
 *
 * @param mp    Map projection information.
 * @param east  Eastings in metres
 * @param north Northings in metres
 * @param lat   Pointer to Latitude in degrees.
 * @param lon   Pointer to Longitude in degrees.
 */
static void merc_inverse(map_proj_t *mp, double east, double north,
                         double *lat, double *lon)
{
  mercator_t *data = (mercator_t *)mp->private_data;
  double f1 = data->e2 / 2.0 + (5.0 / 24.0) * data->e4 + data->e6 / 12.0;
  double f2 = (7.0 / 48.0) * data->e4 + (29.0 / 240.0) * data->e6;
  double rlon = (east - mp->falsex) / mp->ellip_major + data->cm;
  double phi =
    HALF_PI - 2.0 * atan(exp(-(north - mp->falsey) / mp->ellip_major));
  double rlat = phi + f1 * sin(2.0 * phi) + f2 * sin(4.0 * phi);

  *lat = RAD2DEG(rlat);
  *lon = RAD2DEG(rlon);
}


/* Transverse Mercator grid projection.
 */
typedef struct {
  double cm;                    /* Centeral meridian in radians. */
  double k0;                    /* Central meridian scale factor. */
  double e2;
  double e4;
  double e6;
  double edash2;
  double mdF0;                  /* Meridianal distance factors. */
  double mdF1;
  double mdF2;
  double mdF3;
  double mdF4;
  double fpF0;                  /* Footpoint latitude factors. */
  double fpF1;
  double fpF2;
  double fpF3;
  double fpF4;
} trans_mercator_t;

/* Transverse Mercator map projection constructor.
 *
 * The transvere mercator is a cylindrical map projection (transverse).
 * Both the meridians and parallels are complex curves.
 *
 * @param mp Pointer to map projection information.
 * @param nargs Number of parameter arguements.
 * @param args Array of parameter arguments as strings. The
 * arguments conform to the proj conventions.
 * @return Pointer to data private to the trans. mercator projection.
 */
static void *tcc_init(map_proj_t *mp, int nargs, char *args[])
{
  double e2sq, e1, e1_2, e1_3, e1_4;
  trans_mercator_t *data =
    (trans_mercator_t *)malloc(sizeof(trans_mercator_t));

  data->cm = DEG2RAD(CENTRAL_MERIDIAN(nargs, args));
  data->k0 = SCALE_FACTOR(nargs, args);

  data->e2 = (2.0 * mp->ellip_flat) - (mp->ellip_flat * mp->ellip_flat);
  data->e4 = data->e2 * data->e2;
  data->e6 = data->e2 * data->e2 * data->e2;
  data->edash2 = data->e2 / (1.0 - data->e2);

  /* Compute meridian distance factors. */
  data->mdF0 = mp->ellip_major * (1.0 - data->e2 / 4.0
                                  - (3.0 / 64.0) * data->e4 -
                                  (5.0 / 256.0) * data->e6);
  data->mdF1 =
    1.0 - data->e2 / 4.0 - (3.0 / 64.0) * data->e4 -
    (5.0 / 256.0) * data->e6;
  data->mdF2 =
    (3.0 / 8.0) * data->e2 + (3.0 / 32.0) * data->e4 +
    (45.0 / 1024.0) * data->e6;
  data->mdF3 = (15.0 / 256.0) * data->e4 + (45.0 / 1024.0) * data->e6;
  data->mdF4 = (35.0 / 3072.0) * data->e6;

  /* Compute footpoint latitude factors. */
  e2sq = sqrt(1 - data->e2);
  e1 = (1 - e2sq) / (1 + e2sq);
  e1_2 = e1 * e1;
  e1_3 = e1_2 * e1;
  e1_4 = e1_3 * e1;

  data->fpF0 = (3.0 / 2.0) * e1 - (27.0 / 32.0) * e1_3;
  data->fpF1 = (21.0 / 16.0) * e1_2 - (55.0 / 32.0) * e1_4;
  data->fpF2 = (151.0 / 96.0) * e1_3;
  data->fpF3 = (1097.0 / 512.0) * e1_4;

  return data;
}

/* Free the private memory associated with a trans. mercator transformation.
 */
static void tcc_cleanup(void *data)
{
  if (data != NULL)
    free(data);
}


/* Method to convert latitudes and longitudes
 * to a projected Transverse Mercator eastings and northings.
 *
 * @param mp    Map projection information.
 * @param lat   Latitude in degrees.
 * @param lon   Longitude in degrees.
 * @param east  Pointer to Easting in metres
 * @param north Pointer to Northing in metres
 */
static void tcc_forward(map_proj_t *mp, double lat, double lon,
                        double *east, double *north)
{
  trans_mercator_t *data = (trans_mercator_t *)mp->private_data;
  double m, sincoslat, sin2lat, sin4lat, sin6lat, M;
  double x = 0;
  double y = 0;
  double rlat = DEG2RAD(lat);
  double rlon = DEG2RAD(lon);
  double M0 = 0.0;              /* merid dist at lat 0. */
  double coslat = cos(rlat);
  double sinlat = sin(rlat);
  double sinlat2 = sinlat * sinlat;
  double coslat2 = coslat * coslat;

/* If not at the poles then do the calculation, elese this
 * is a special case.
 */
  if (rlat != HALF_PI && rlat != -HALF_PI) {
    double tanlat = sinlat / coslat;
    double N = mp->ellip_major / sqrt(1 - data->e2 * sinlat2);
    double T = tanlat * tanlat;
    double T2 = T * T;
    double C = data->edash2 * coslat2;
    double C2 = C * C;
    double A = (rlon - data->cm) * coslat;
    double A2 = A * A;
    double A3 = A2 * A;
    double A4 = A3 * A;
    double A5 = A4 * A;
    double A6 = A5 * A;

    x = A + (1.0 - T + C) * A3 / 6.0
      + (5.0 - 18.0 * T + T2 + 72.0 * C -
         58.0 * data->edash2) * A5 / 120.0;
    x *= data->k0 * N;
    y = N * tanlat * (A2 / 2.0 + (5.0 - T + 9.0 * C + 4.0 * C2) * A4 / 24.0
                      + (61.0 - 58.0 * T + T2 + 600.0 * C -
                         330 * data->edash2) * A6 / 720.0);
  }

  /* Compute meridian distance. */
  m = data->mdF1 * rlat;
  sincoslat = sinlat * coslat;
  sin2lat = 2.0 * sincoslat;
  sin4lat = sin2lat * (4.0 * coslat2 - 2);
  sin6lat = 3 * sin2lat - 4 * sin2lat * sin2lat * sin2lat;
  m -= data->mdF2 * sin2lat;
  m += data->mdF3 * sin4lat;
  m -= data->mdF4 * sin6lat;
  M = mp->ellip_major * m;

  y = data->k0 * (M - M0 + y);

  /* reduce to false origins */
  *east = x + mp->falsex;
  *north = y + mp->falsey;
}

/* Method to convert projected Mercator eastings and northings
 * to latitudes and longitudes.
 *
 * @param mp    Map projection information.
 * @param east  Eastings in metres
 * @param north Northings in metres
 * @param lat   Pointer to Latitude in degrees.
 * @param lon   Pointer to Longitude in degrees.
 */
static void tcc_inverse(map_proj_t *mp, double east, double north,
                        double *lat, double *lon)
{
  trans_mercator_t *data = (trans_mercator_t *)mp->private_data;
  double M0 = 0.0;
  double M = M0 + (north - mp->falsey) / data->k0;
  double phi1 = tcc_footpoint_latitude(mp, M);
  double coslat = cos(phi1);
  double sinlat = sin(phi1);
  double sinlat2 = sinlat * sinlat;
  double coslat2 = coslat * coslat;
  double tanlat = sinlat / coslat;
  double N = mp->ellip_major / sqrt(1 - data->e2 * sinlat2);
  double R = mp->ellip_major * (1 - data->e2)
    / pow(1 - data->e2 * sinlat2, 3.0 / 2.0);
  double T = tanlat * tanlat;
  double T2 = T * T;
  double C = data->edash2 * coslat2;
  double C2 = C * C;
  double D = (east - mp->falsex) / (N * data->k0);
  double D2 = D * D;
  double D3 = D2 * D;
  double D4 = D3 * D;
  double D5 = D4 * D;
  double D6 = D5 * D;
  double tlat = (N * tanlat / R) * ((D2 / 2.0 -
                                     (5.0 + 3 * T + 10 * C - 4 * C2 -
                                      9 * data->edash2) * D4 / 24.0 + (61 +
                                                                       90 *
                                                                       T +
                                                                       298
                                                                       *
                                                                       C +
                                                                       45 *
                                                                       T2 -
                                                                       252
                                                                       *
                                                                       data->
                                                                       edash2
                                                                       -
                                                                       3 *
                                                                       C2)
                                     * D6 / 720.0));
  double tlon =
    (D - (1. + 2 * T + C) * D3 / 6.0 +
     (5 - 2 * C + 28 * T - 3 * C2 + 8 * data->edash2 +
      24 * T2) * D5 / 120.0) / coslat;

  *lat = RAD2DEG(phi1 - tlat);
  *lon = RAD2DEG(data->cm + tlon);
}

/* Calculates the "footpoint" latitude, which is the latitude
 * at the central meridian which has the same y coordinate as that
 * of (lat,lon).
 *
 * @param mp   Map projection information.
 * @param M	Meridianal distance at known latitude.
 * @return	Footpoint latitude.
 */
static double tcc_footpoint_latitude(map_proj_t *mp, double M)
{
  trans_mercator_t *data = (trans_mercator_t *)mp->private_data;
  double mu = M / data->mdF0;
  double phi = mu;

  phi += data->fpF0 * sin(2.0 * mu);
  phi += data->fpF1 * sin(4.0 * mu);
  phi += data->fpF2 * sin(6.0 * mu);

  return phi + data->fpF3 * sin(8.0 * mu);
}


/* Lambert Conformal Conic grid projection.
 */
typedef struct {
  double cm;                    /* Centeral meridian in radians. */
  double phi0;                  /* Origin latitude in radians. */
  double phi1;                  /* First standard parallel in radians. */
  double phi2;                  /* Second standard parallel in radians. */
  double e;
  double e2;
  double e4;
  double e6;
  double e8;
  double cosphi1;
  double cosphi2;
  double sinphi0;
  double sinphi1;
  double sinphi2;
  double m1;
  double m2;
  double t0;
  double t1;
  double t2;
  double nn;
  double F;
  double rho0;
  double A, B, C, D;
} lambert_conic_t;


/* Lambert Conformal Conic grid projection.
 *
 * @param mp Pointer to map projection information.
 * @param nargs Number of parameter arguements.
 * @param args Array of parameter arguments as strings. The
 * arguments conform to the proj conventions.
 * @return Pointer to data private to the projection.
 */
static void *lcc_init(map_proj_t *mp, int nargs, char *args[])
{
  double oA, oB, oC, oD;
  lambert_conic_t *data =
    (lambert_conic_t *)malloc(sizeof(lambert_conic_t));
  data->cm = DEG2RAD(CENTRAL_MERIDIAN(nargs, args));
  data->phi1 = DEG2RAD(FIRST_LATITUDE(nargs, args));
  data->phi2 = DEG2RAD(SECOND_LATITUDE(nargs, args));

  data->e2 = (2.0 * mp->ellip_flat) - (mp->ellip_flat * mp->ellip_flat);
  data->e = sqrt(data->e2);
  data->e4 = data->e2 * data->e2;
  data->e6 = data->e2 * data->e2 * data->e2;
  data->e8 = data->e4 * data->e4;

  data->cosphi1 = cos(data->phi1);
  data->cosphi2 = cos(data->phi2);
  data->sinphi1 = sin(data->phi1);
  data->sinphi2 = sin(data->phi2);
  data->m1 =
    data->cosphi1 / sqrt(1 - data->e2 * data->sinphi1 * data->sinphi1);
  data->m2 =
    data->cosphi2 / sqrt(1 - data->e2 * data->sinphi2 * data->sinphi2);
  data->t1 = sqrt(((1 - data->sinphi1) / (1 + data->sinphi1))
                  * pow((1 + data->e * data->sinphi1) /
                        (1 - data->e * data->sinphi1), data->e));
  data->t2 = sqrt(((1 - data->sinphi2) / (1 + data->sinphi2))
                  * pow((1 + data->e * data->sinphi2) /
                        (1 - data->e * data->sinphi2), data->e));
  data->nn = log(data->m1 / data->m2) / log(data->t1 / data->t2);

  /* If lat_0 is not provided, then compute it. */
  if (find_argument("lat_0", nargs, args) < 0) {
    data->sinphi0 = data->nn;
    data->phi0 = asin(data->sinphi0);
  } else {
    data->phi0 = DEG2RAD(CENTRAL_LATITUDE(nargs, args));
    data->sinphi0 = sin(data->phi0);
  }

  data->t0 = sqrt(((1 - data->sinphi0) / (1 + data->sinphi0))
                  * pow((1 + data->e * data->sinphi0) /
                        (1 - data->e * data->sinphi0), data->e));
  data->F = data->m1 / (data->nn * pow(data->t1, data->nn));
  data->rho0 = mp->ellip_major * data->F * pow(data->t0, data->nn);

  oA =
    data->e2 / 2.0 + 5.0 * data->e4 / 24.0 + data->e6 / 12.0 +
    13.0 * data->e8 / 360.0;
  oB =
    7.0 * data->e4 / 48.0 + 29.0 * data->e6 / 240.0 +
    811.0 * data->e8 / 11520.0;
  oC = 7.0 * data->e6 / 120.0 + 81.0 * data->e8 / 1120.0;
  oD = 4279.0 * data->e8 / 161280.0;

  data->A = oA - oC;
  data->B = 2.0 * oB - 4.0 * oD;
  data->C = 4.0 * oC;
  data->D = 8.0 * oD;

  return data;
}

/* Free the private memory associated with a LCC projection.
 */
static void lcc_cleanup(void *data)
{
  if (data != NULL)
    free(data);
}


/* Method to convert latitudes and longitudes
 * to a projected Lambert Conformal Conic eastings and northings.
 *
 * @param mp    Map projection information.
 * @param lat   Latitude in degrees.
 * @param lon   Longitude in degrees.
 * @param east  Pointer to Easting in metres
 * @param north Pointer to Northing in metres
 */
static void lcc_forward(map_proj_t *mp, double lat, double lon,
                        double *east, double *north)
{
  lambert_conic_t *data = (lambert_conic_t *)mp->private_data;
  double rlat = DEG2RAD(lat);
  double rlon = DEG2RAD(lon);
  double theta = data->nn * (rlon - data->cm);
  double sinlat = sin(rlat);
  double esinlat = data->e * sinlat;
  double t = sqrt(((1 - sinlat) / (1 + sinlat))
                  * pow((1 + esinlat) / (1 - esinlat), data->e));
  double rho = mp->ellip_major * data->F * pow(t, data->nn);

  *east = rho * sin(theta) + mp->falsex;
  *north = data->rho0 - rho * cos(theta) + mp->falsey;
}

/* Method to convert projected LCC eastings and northings
 * to latitudes and longitudes.
 *
 * @param mp    Map projection information.
 * @param east  Eastings in metres
 * @param north Northings in metres
 * @param lat   Pointer to Latitude in degrees.
 * @param lon   Pointer to Longitude in degrees.
 */
static void lcc_inverse(map_proj_t *mp, double east, double north,
                        double *lat, double *lon)
{
  lambert_conic_t *data = (lambert_conic_t *)mp->private_data;
  double s = (data->nn >= 0) ? 1 : -1;
  double x = east - mp->falsex;
  double y = data->rho0 - (north - mp->falsey);
  double rho = s * sqrt(x * x + y * y);
  double t = pow(rho / (mp->ellip_major * data->F), 1.0 / data->nn);
  double theta = atan2(s * x, s * y);
  double phi = HALF_PI - 2 * atan(t);
  double sinphi = sin(2 * phi);
  double cosphi = cos(2 * phi);
  double rlat =
    phi + sinphi * (data->A +
                    cosphi * (data->B +
                              cosphi * (data->C + data->D * cosphi)));

  *lat = RAD2DEG(rlat);
  *lon = RAD2DEG(theta / data->nn + data->cm);
}


#define SEROV_X0 -1.55010634e+04
#define SEROV_Y0 7.787528320e+03
#define SEROV_SCALE 1.113231250e+03

/* CMR Remote sensing serov projection. This projection is
 * based on the Lambert Conformal Conic with a offset and
 * scaling centered around a small Russian town called Serov.
 *
 * @param mp Pointer to map projection information.
 * @param nargs Number of parameter arguements.
 * @param args Array of parameter arguments as strings. The
 * arguments conform to the proj conventions.
 * @return Pointer to data private to the projection.
 */
static void *serov_init(map_proj_t *mp, int nargs, char *args[])
{
  int new_nargs;
  char **new_args;

  new_args = add_argument("lon_0=150", nargs, args, &new_nargs);
  new_args = add_argument("lat_1=-40", new_nargs, new_args, &new_nargs);
  new_args = add_argument("lat_2=-10", new_nargs, new_args, &new_nargs);
  return lcc_init(mp, new_nargs, new_args);
}

/* Free the private memory associated with a Serov projection.
 */
static void serov_cleanup(void *data)
{
  lcc_cleanup(data);
}


/* Method to convert latitudes and longitudes
 * to a projected Serov eastings and northings.
 *
 * @param mp    Map projection information.
 * @param lat   Latitude in degrees.
 * @param lon   Longitude in degrees.
 * @param east  Pointer to Easting in metres
 * @param north Pointer to Northing in metres
 */
static void serov_forward(map_proj_t *mp, double lat, double lon,
                          double *east, double *north)
{
  lcc_forward(mp, lat, lon, east, north);

  *east = *east / SEROV_SCALE - SEROV_X0;
  *north = SEROV_Y0 - *north / SEROV_SCALE;
}

/* Method to convert projected Serov eastings and northings
 * to latitudes and longitudes.
 *
 * @param mp    Map projection information.
 * @param east  Eastings in metres
 * @param north Northings in metres
 * @param lat   Pointer to Latitude in degrees.
 * @param lon   Pointer to Longitude in degrees.
 */
static void serov_inverse(map_proj_t *mp, double east, double north,
                          double *lat, double *lon)
{
  east = (east + SEROV_X0) * SEROV_SCALE;
  north = (SEROV_Y0 - north) * SEROV_SCALE;

  lcc_inverse(mp, east, north, lat, lon);
}

#if defined(HAVE_USGS_PROJ)
/* USGS PROJ 4 library.
 *
 * @param mp Pointer to map projection information.
 * @param nargs Number of parameter arguements.
 * @param args Array of parameter arguments as strings. The
 * arguments conform to the proj conventions.
 * @return Pointer to data private to the projection.
 */
static void *usgs_init(map_proj_t *mp, int nargs, char *args[])
{
  return pj_init(nargs, args);
}

/* Free the private memory associated with PROJ 4.
 */
static void usgs_cleanup(void *data)
{
  pj_free((PJ *) data);
}


/* Method to convert latitudes and longitudes
 * to eastings and northings using the PROJ 4 library.
 *
 * @param mp    Map projection information.
 * @param lat   Latitude in degrees.
 * @param lon   Longitude in degrees.
 * @param east  Pointer to Easting in metres
 * @param north Pointer to Northing in metres
 */
static void usgs_forward(map_proj_t *mp, double lat, double lon,
                         double *east, double *north)
{
  LP lp;
  XY xy;
  lp.u = DEG2RAD(lon);
  lp.v = DEG2RAD(lat);
  xy = pj_fwd(lp, (PJ *) mp->private_data);
  *east = xy.u;
  *north = xy.v;
}

/* Method to convert eastings and northings
 * to latitudes and longitudes using the PROJ 4 library.
 *
 * @param mp    Map projection information.
 * @param east  Eastings in metres
 * @param north Northings in metres
 * @param lat   Pointer to Latitude in degrees.
 * @param lon   Pointer to Longitude in degrees.
 */
static void usgs_inverse(map_proj_t *mp, double east, double north,
                         double *lat, double *lon)
{
  XY xy;
  LP lp;
  xy.u = east;
  xy.v = north;
  lp = pj_inv(xy, (PJ *) mp->private_data);
  *lon = RAD2DEG(lp.u);
  *lat = RAD2DEG(lp.v);
}
#endif
