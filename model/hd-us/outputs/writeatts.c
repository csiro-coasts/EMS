/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/outputs/writeatts.c
 *  
 *  Description:
 *  Write the attributes for each variable.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: writeatts.c 7465 2023-12-13 03:52:41Z her127 $
 *
 */

#include <netcdf.h>
#include <string.h>
#include "hd.h"
#include "tracer.h"

static void write_analytic_rect_att(dump_data_t *dumpdata, int cdfid,
                                    int vid, int ilower, int ne1,
                                    int jlower, int ne2, int xylocation);
static void write_analytic_polar_att(dump_data_t *dumpdata, int cdfid,
                                     int vid, int ilower, int ne1,
                                     int jlower, int ne2, int xylocation);
static void write_grid_atts(dump_data_t *dumpdata, int cdfid, int ilower,
                            int jlower);



void write_dump_attributes(dump_data_t *dumpdata, int cdfid,
                           nc_type fptype, int ilower, int nce1, int nfe1,
                           int jlower, int nce2, int nfe2, char *modulo,
			   char *tunit)
{
  /* dimension ids */
  int vid;
  int has_proj = (strlen(projection) > 0);
  int is_geog = has_proj && (strcasecmp(projection, GEOGRAPHIC_TAG) == 0);
  double fill;
  int ns2 = dumpdata->ns2;
  int ns3 = dumpdata->ns3;
  char buf[MAXSTRLEN];

  /* time independent variables */
  vid = ncw_var_id(cdfid, "z_grid");
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name",
                 "Z coordinate at grid layer faces");
  write_text_att(cdfid, vid, "coordinate_type", "Z");

  vid = ncw_var_id(cdfid, "z_centre");
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name",
                 "Z coordinate at grid layer centre");
  write_text_att(cdfid, vid, "coordinate_type", "Z");

  vid = ncw_var_id(cdfid, "x_grid");
  if (is_geog) {
    write_text_att(cdfid, vid, "long_name", "Longitude at grid corners");
    write_text_att(cdfid, vid, "coordinate_type", "longitude");
    write_text_att(cdfid, vid, "units", "degrees_east");
  } else {
    write_text_att(cdfid, vid, "long_name",
                   "X coordinate at grid corners");
    write_text_att(cdfid, vid, "coordinate_type", "X");
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  }
  if (has_proj)
    write_text_att(cdfid, vid, "projection", projection);
  write_analytic_att(dumpdata, cdfid, vid, ilower, nce1,
                     jlower, nce2, CL_GRID);

  vid = ncw_var_id(cdfid, "y_grid");
  if (is_geog) {
    write_text_att(cdfid, vid, "long_name", "Latitude at grid corners");
    write_text_att(cdfid, vid, "coordinate_type", "latitude");
    write_text_att(cdfid, vid, "units", "degrees_north");
  } else {
    write_text_att(cdfid, vid, "long_name",
                   "Y coordinate at grid corners");
    write_text_att(cdfid, vid, "coordinate_type", "Y");
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  }
  if (has_proj)
    write_text_att(cdfid, vid, "projection", projection);
  write_analytic_att(dumpdata, cdfid, vid, ilower, nce1, jlower, nce2,
                     CL_GRID);

  vid = ncw_var_id(cdfid, "x_centre");
  if (is_geog) {
    write_text_att(cdfid, vid, "long_name", "Longitude at cell centre");
    write_text_att(cdfid, vid, "coordinate_type", "longitude");
    write_text_att(cdfid, vid, "units", "degrees_east");
  } else {
    write_text_att(cdfid, vid, "long_name", "X coordinate at cell centre");
    write_text_att(cdfid, vid, "coordinate_type", "X");
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  }
  if (has_proj)
    write_text_att(cdfid, vid, "projection", projection);
  write_analytic_att(dumpdata, cdfid, vid, ilower, nce1, jlower, nce2,
                     CL_CENTRE);

  vid = ncw_var_id(cdfid, "y_centre");
  if (is_geog) {
    write_text_att(cdfid, vid, "long_name", "Latitude at cell centre");
    write_text_att(cdfid, vid, "coordinate_type", "latitude");
    write_text_att(cdfid, vid, "units", "degrees_north");
  } else {
    write_text_att(cdfid, vid, "long_name", "Y coordinate at cell centre");
    write_text_att(cdfid, vid, "coordinate_type", "Y");
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  }
  if (has_proj)
    write_text_att(cdfid, vid, "projection", projection);
  write_analytic_att(dumpdata, cdfid, vid, ilower, nce1, jlower, nce2,
                     CL_CENTRE);

  vid = ncw_var_id(cdfid, "x_left");
  if (is_geog) {
    write_text_att(cdfid, vid, "long_name",
                   "Longitude at centre of left face");
    write_text_att(cdfid, vid, "coordinate_type", "longitude");
    write_text_att(cdfid, vid, "units", "degrees_east");
  } else {
    write_text_att(cdfid, vid, "long_name",
                   "X coordinate at centre of left face");
    write_text_att(cdfid, vid, "coordinate_type", "X");
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  }
  if (has_proj)
    write_text_att(cdfid, vid, "projection", projection);
  write_analytic_att(dumpdata, cdfid, vid, ilower, nce1, jlower, nce2,
                     CL_LEFT);

  vid = ncw_var_id(cdfid, "y_left");
  if (is_geog) {
    write_text_att(cdfid, vid, "long_name",
                   "Latitude at centre of left face");
    write_text_att(cdfid, vid, "coordinate_type", "latitude");
    write_text_att(cdfid, vid, "units", "degrees_north");
  } else {
    write_text_att(cdfid, vid, "long_name",
                   "Y coordinate at centre of left face");
    write_text_att(cdfid, vid, "coordinate_type", "Y");
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  }
  if (has_proj)
    write_text_att(cdfid, vid, "projection", projection);
  write_analytic_att(dumpdata, cdfid, vid, ilower, nce1,
                     jlower, nce2, CL_LEFT);

  vid = ncw_var_id(cdfid, "x_back");
  if (is_geog) {
    write_text_att(cdfid, vid, "long_name",
                   "Longitude at centre of back face");
    write_text_att(cdfid, vid, "coordinate_type", "longitude");
    write_text_att(cdfid, vid, "units", "degrees_east");
  } else {
    write_text_att(cdfid, vid, "long_name",
                   "X coordinate at centre of back face");
    write_text_att(cdfid, vid, "coordinate_type", "X");
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  }
  if (has_proj)
    write_text_att(cdfid, vid, "projection", projection);
  write_analytic_att(dumpdata, cdfid, vid, ilower, nce1, jlower, nce2,
                     CL_BACK);

  vid = ncw_var_id(cdfid, "y_back");
  if (is_geog) {
    write_text_att(cdfid, vid, "long_name",
                   "Latitude at centre of back face");
    write_text_att(cdfid, vid, "coordinate_type", "latitude");
    write_text_att(cdfid, vid, "units", "degrees_north");
  } else {
    write_text_att(cdfid, vid, "long_name",
                   "Y coordinate at centre of back face");
    write_text_att(cdfid, vid, "coordinate_type", "Y");
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  }
  if (has_proj)
    write_text_att(cdfid, vid, "projection", projection);
  write_analytic_att(dumpdata, cdfid, vid, ilower, nce1, jlower, nce2,
                     CL_BACK);

  vid = ncw_var_id(cdfid, "botz");
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name",
                 "Z coordinate at sea-bed at cell centre");
  /*write_text_att(cdfid, vid, "coordinates", "latitude, longitude");&*/

  vid = ncw_var_id(cdfid, "h1au1");
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name",
                 "Cell width at centre of left face");
  write_text_att(cdfid, vid, "coordinates", "x_left, y_left");

  vid = ncw_var_id(cdfid, "h2au1");
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name",
                 "Cell height at centre of left face");
  write_text_att(cdfid, vid, "coordinates", "x_left, y_left");

  vid = ncw_var_id(cdfid, "h1au2");
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name",
                 "Cell width at centre of back face");
  write_text_att(cdfid, vid, "coordinates", "x_back, y_back");

  vid = ncw_var_id(cdfid, "h2au2");
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name",
                 "Cell height at centre of back face");
  write_text_att(cdfid, vid, "coordinates", "x_back, y_back");

  vid = ncw_var_id(cdfid, "h1acell");
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name", "Cell width at cell centre");
  write_text_att(cdfid, vid, "coordinates", "x_centre, y_centre");

  vid = ncw_var_id(cdfid, "h2acell");
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name", "Cell height at cell centre");
  write_text_att(cdfid, vid, "coordinates", "x_centre, y_centre");

  vid = ncw_var_id(cdfid, "h1agrid");
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name", "Cell width at grid corner");
  write_text_att(cdfid, vid, "coordinates", "x_grid, y_grid");

  vid = ncw_var_id(cdfid, "h2agrid");
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name", "Cell height at grid corner");
  write_text_att(cdfid, vid, "coordinates", "x_grid, y_grid");

  vid = ncw_var_id(cdfid, "thetau1");
  write_text_att(cdfid, vid, "units", "radian");
  write_text_att(cdfid, vid, "long_name",
                 "Cell rotation at centre of left face");
  write_text_att(cdfid, vid, "coordinates", "x_left, y_left");

  vid = ncw_var_id(cdfid, "thetau2");
  write_text_att(cdfid, vid, "units", "radian");
  write_text_att(cdfid, vid, "long_name",
                 "Cell rotation at centre of back face");
  write_text_att(cdfid, vid, "coordinates", "x_back, y_back");

  vid = ncw_var_id(cdfid, "coriolis");
  write_text_att(cdfid, vid, "units", " ");
  write_text_att(cdfid, vid, "long_name", "Coriolis parameter");
  write_text_att(cdfid, vid, "coordinates", "x_centre, y_centre");


  /* time dependent variables */
  vid = ncw_var_id(cdfid, "t");
  write_text_att(cdfid, vid, "units", tunit);
  write_text_att(cdfid, vid, "long_name", "Time");
  write_text_att(cdfid, vid, "coordinate_type", "time");
  if (strlen(modulo))
    write_text_att(cdfid, vid, "modulo", modulo);

  if ((vid = ncw_var_id(cdfid, "z_grid_sed")) >= 0) {
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
    write_text_att(cdfid, vid, "long_name",
                   "Z coordinate at sediment layer faces");
    write_text_att(cdfid, vid, "coordinate_type", "Z");
  }

  if ((vid = ncw_var_id(cdfid, "z_centre_sed")) >= 0) {
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
    write_text_att(cdfid, vid, "long_name",
                   "Z coordinate at sediment layer centres");
    write_text_att(cdfid, vid, "coordinate_type", "Z");
  }

  if (dumpdata->large) {
    if ((vid = ncw_var_id(cdfid, "u1av")) >= 0) {
      write_text_att(cdfid, vid, "units", "metre second-1");
      write_text_att(cdfid, vid, "long_name",
		     "I component of depth averaged current at left face");
      write_text_att(cdfid, vid, "coordinates", "t, x_left, y_left");
    }
    
    if ((vid = ncw_var_id(cdfid, "u2av")) >= 0) {
      write_text_att(cdfid, vid, "units", "metre second-1");
      write_text_att(cdfid, vid, "long_name",
		     "J component of depth averaged current at back face");
      write_text_att(cdfid, vid, "coordinates", "t, x_back, y_back");
    }
    
    if ((vid = ncw_var_id(cdfid, "wtop")) >= 0) {
      write_text_att(cdfid, vid, "units", "metre second-1");
      write_text_att(cdfid, vid, "long_name",
		     "Vertical velocity at surface");
      write_text_att(cdfid, vid, "coordinates", "t, x_centre, y_centre");
    }
    
    if ((vid = ncw_var_id(cdfid, "topz")) >= 0) {
      write_text_att(cdfid, vid, "units", dumpdata->lenunit);
      write_text_att(cdfid, vid, "long_name",
		     "Z coordinate for surface cell");
      write_text_att(cdfid, vid, "coordinates", "t, x_centre, y_centre");
    }
  }

  if ((vid = ncw_var_id(cdfid, "eta")) >= 0) {
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
    write_text_att(cdfid, vid, "long_name", "Surface Elevation");
    write_text_att(cdfid, vid, "coordinates", "t, x_centre, y_centre");
  }

  if ((vid = ncw_var_id(cdfid, "thickness_s")) >= 0) {
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
    write_text_att(cdfid, vid, "long_name", "Sediment thickness");
    write_text_att(cdfid, vid, "coordinates", "t, x_centre, y_centre");
  }

  if ((vid = ncw_var_id(cdfid, "resuspension")) >= 0) {
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
    write_text_att(cdfid, vid, "long_name", "Resuspension");
    write_text_att(cdfid, vid, "coordinates", "t, x_centre, y_centre");
  }

  if ((vid = ncw_var_id(cdfid, "resuspension_ac")) >= 0) {
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
    write_text_att(cdfid, vid, "long_name", "Accumulated resuspension");
    write_text_att(cdfid, vid, "coordinates", "t, x_centre, y_centre");
  }

  if ((vid = ncw_var_id(cdfid, "deposition_ac")) >= 0) {
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
    write_text_att(cdfid, vid, "long_name", "Accumulated deposition");
    write_text_att(cdfid, vid, "coordinates", "t, x_centre, y_centre");
  }

  if ((vid = ncw_var_id(cdfid, "wind1")) >= 0) {
    write_text_att(cdfid, vid, "units", "Nm-2");
    write_text_att(cdfid, vid, "long_name",
                   "I component of wind stress at left face");
    write_text_att(cdfid, vid, "coordinates", "t, x_left, y_left");
  }

  if ((vid = ncw_var_id(cdfid, "wind2")) >= 0) {
    write_text_att(cdfid, vid, "units", "Nm-2");
    write_text_att(cdfid, vid, "long_name",
                   "J component of wind stress at back face");
    write_text_att(cdfid, vid, "coordinates", "t, x_back, y_back");
  }

  if ((vid = ncw_var_id(cdfid, "patm")) >= 0) {
    write_text_att(cdfid, vid, "units", "Pa");
    write_text_att(cdfid, vid, "long_name", "Atmospheric pressure");
    write_text_att(cdfid, vid, "coordinates", "t, x_centre, y_centre");
  }

  if (dumpdata->large) {
    if ((vid = ncw_var_id(cdfid, "u1")) >= 0) {
      write_text_att(cdfid, vid, "units", "metre second-1");
      write_text_att(cdfid, vid, "long_name",
		     "I component of current at left face");
      write_text_att(cdfid, vid, "coordinates",
		     "t, x_left, y_left, z_centre");
    }

    if ((vid = ncw_var_id(cdfid, "u2")) >= 0) {
      write_text_att(cdfid, vid, "units", "metre second-1");
      write_text_att(cdfid, vid, "long_name",
		     "J component of current at back face");
      write_text_att(cdfid, vid, "coordinates",
		     "t, x_back, y_back, z_centre");
    }

    if ((vid = ncw_var_id(cdfid, "u")) >= 0) {
      write_text_att(cdfid, vid, "units", "metre second-1");
      write_text_att(cdfid, vid, "long_name",
		     "East component of current at cell centre");
      write_text_att(cdfid, vid, "coordinates",
		     "t, x_centre, y_centre, z_centre");
    }

    if ((vid = ncw_var_id(cdfid, "v")) >= 0) {
      write_text_att(cdfid, vid, "units", "metre second-1");
      write_text_att(cdfid, vid, "long_name",
		     "North component of current at cell centre");
      write_text_att(cdfid, vid, "coordinates",
		     "t, x_centre, y_centre, z_centre");
    }

    if ((vid = ncw_var_id(cdfid, "w")) >= 0) {
      write_text_att(cdfid, vid, "units", "metre second-1");
      write_text_att(cdfid, vid, "long_name",
		     "K component of current at cell centre and Z grid");
      write_text_att(cdfid, vid, "coordinates",
		     "t, x_centre, y_centre, z_centre");
    }

    if ((vid = ncw_var_id(cdfid, "dens")) >= 0) {
      write_text_att(cdfid, vid, "units", "kg metre-3");
      write_text_att(cdfid, vid, "long_name", "Density");
      write_text_att(cdfid, vid, "coordinates",
		     "t, x_centre, y_centre, z_centre");
    }

    if ((vid = ncw_var_id(cdfid, "dens_0")) >= 0) {
      write_text_att(cdfid, vid, "units", "kg metre-3");
      write_text_att(cdfid, vid, "long_name", "Potential density");
      write_text_att(cdfid, vid, "coordinates",
		     "t, x_centre, y_centre, z_centre");
      fill = 1025.0;
      nc_put_att_double(cdfid, vid, "_FillValue", fptype, 1, &fill);
    }

    if ((vid = ncw_var_id(cdfid, "Kz")) >= 0) {
      write_text_att(cdfid, vid, "units", "m2 s-1");
      write_text_att(cdfid, vid, "long_name", "Kz");
      write_text_att(cdfid, vid, "coordinates",
		     "t, x_centre, y_centre, z_centre");
      fill = 0;
      nc_put_att_double(cdfid, vid, "_FillValue", fptype, 1, &fill);
    }

    if ((vid = ncw_var_id(cdfid, "Vz")) >= 0) {
      write_text_att(cdfid, vid, "units", "m2 s-1");
      write_text_att(cdfid, vid, "long_name", "Vz");
      write_text_att(cdfid, vid, "coordinates",
		     "t, x_centre, y_centre, z_centre");
      fill = 0;
      nc_put_att_double(cdfid, vid, "_FillValue", fptype, 1, &fill);
    }

    if ((vid = ncw_var_id(cdfid, "u1bot")) >= 0) {
      write_text_att(cdfid, vid, "units", "metre second-1");
      write_text_att(cdfid, vid, "long_name",
		     "I component of bottom current deviation at left face");
      write_text_att(cdfid, vid, "coordinates", "t, x_left, y_left");
    }

    if ((vid = ncw_var_id(cdfid, "u2bot")) >= 0) {
      write_text_att(cdfid, vid, "units", "metre second-1");
      write_text_att(cdfid, vid, "long_name",
		     "J component of bottom current deviation at back face");
      write_text_att(cdfid, vid, "coordinates", "t, x_back, y_back");
    }
  }

  if ((vid = ncw_var_id(cdfid, "ptconc")) >= 0) {
    write_text_att(cdfid, vid, "units", "kg metre-3");
    write_text_att(cdfid, vid, "long_name", "Particle concentration");
    write_text_att(cdfid, vid, "coordinates",
                   "t, x_centre, y_centre, z_centre");
  }

  if ((vid = ncw_var_id(cdfid, "Cd")) >= 0) {
    write_text_att(cdfid, vid, "units", " ");
    write_text_att(cdfid, vid, "long_name", "Bottom drag coefficient");
    write_text_att(cdfid, vid, "coordinates", "t, x_centre, y_centre");
  }

  if ((vid = ncw_var_id(cdfid, "ustrcw")) >= 0) {
    write_text_att(cdfid, vid, "units", "ms-1");
    write_text_att(cdfid, vid, "long_name", "Bottom friction velocity");
    write_text_att(cdfid, vid, "coordinates", "t, x_centre, y_centre");
  }

  if ((vid = ncw_var_id(cdfid, "u1vh")) >= 0) {
    write_text_att(cdfid, vid, "units", "m2s-1");
    write_text_att(cdfid, vid, "long_name", "Horizontal viscosity e1 direction");
    write_text_att(cdfid, vid, "coordinates", "t, x_centre, y_centre");
  }

   if ((vid = ncw_var_id(cdfid, "u2vh")) >= 0) {
    write_text_att(cdfid, vid, "units", "m2s-1");
    write_text_att(cdfid, vid, "long_name", "Horizontal viscosity e2 direction");
    write_text_att(cdfid, vid, "coordinates", "t, x_centre, y_centre");
  }

  {
    tracer_att_t attr[] = { {"tracer", "true"},
    {"coordinates", "t, x_centre, y_centre, z_centre"}
    };
    tracer_write_nc(cdfid, dumpdata->ntr, dumpdata->trinfo_3d, 2, attr);
  }

  {
    tracer_att_t attr[] = { {"tracer2D", "true"},
    {"coordinates", "t, x_centre, y_centre"}
    };
    tracer_write_nc(cdfid, dumpdata->ntrS, dumpdata->trinfo_2d, 2, attr);
  }

  {
    tracer_att_t attr[] = { {"tracer_sed", "true"},
    {"coordinates", "t, x_centre, y_centre, z_centre_sed"}
    };
    tracer_write_nc(cdfid, dumpdata->nsed, dumpdata->trinfo_sed, 2, attr);
  }

  vid = ncw_var_id(cdfid, "flag");
  if (vid >= 0) {
    write_text_att(cdfid, vid, "long_name", "SHOC masking flags");
    write_text_att(cdfid, vid, "coordinates",
                   "t, x_centre, y_centre, z_centre");
  }

  /* global attributes */
  write_text_att(cdfid, NC_GLOBAL, "title", codeheader);
  write_text_att(cdfid, NC_GLOBAL, "paramhead", parameterheader);
  getcwd(buf, MAXSTRLEN);
  sprintf(buf, "%s/%s", buf, dumpdata->prmname);
  write_text_att(cdfid, NC_GLOBAL, "paramfile", buf);
  write_text_att(cdfid, NC_GLOBAL, "ems_version", version);
  write_text_att(cdfid, NC_GLOBAL, "Conventions", "CMR/Timeseries/SHOC");
  if (dumpdata->runno >= 0)
    write_text_att(cdfid, NC_GLOBAL, "Run_ID", dumpdata->runnoc);
  /*nc_put_att_double(cdfid, NC_GLOBAL, "Run_ID", NC_DOUBLE, 1, &dumpdata->runno);*/
  if (strlen(dumpdata->runcode))
    write_text_att(cdfid, NC_GLOBAL, "Run_code", dumpdata->runcode);
  if (strlen(dumpdata->rev))
    write_text_att(cdfid, NC_GLOBAL, "Parameter_File_Revision", dumpdata->rev);
  if (strlen(dumpdata->trl))
    write_text_att(cdfid, NC_GLOBAL, "Technology_Readiness_Level", dumpdata->trl);
  if (strlen(dumpdata->reference))
    write_text_att(cdfid, NC_GLOBAL, "Output_Reference", dumpdata->reference);
  nc_put_att_int(cdfid, NC_GLOBAL, "nce1", NC_INT, 1, &nce1);
  nc_put_att_int(cdfid, NC_GLOBAL, "nce2", NC_INT, 1, &nce2);
  nc_put_att_int(cdfid, NC_GLOBAL, "nfe1", NC_INT, 1, &nfe1);
  nc_put_att_int(cdfid, NC_GLOBAL, "nfe2", NC_INT, 1, &nfe2);
  nc_put_att_int(cdfid, NC_GLOBAL, "ns3", NC_INT, 1, &ns3);
  nc_put_att_int(cdfid, NC_GLOBAL, "ns2", NC_INT, 1, &ns2);

  write_grid_atts(dumpdata, cdfid, ilower, jlower);
}


void write_analytic_att(dump_data_t *dumpdata, int cdfid, int vid,
                        int ilower, int ne1, int jlower, int ne2,
                        int xylocation)
{
  if (dumpdata->rg)
    write_analytic_rect_att(dumpdata, cdfid, vid, ilower, ne1, jlower, ne2,
                            xylocation);

  else if (dumpdata->pg)
    write_analytic_polar_att(dumpdata, cdfid, vid, ilower, ne1, jlower,
                             ne2, xylocation);
}


static void write_analytic_rect_att(dump_data_t *dumpdata, int cdfid,
                                    int vid, int ilower, int ne1,
                                    int jlower, int ne2, int xylocation)
{
  char line[MAXLINELEN];
  double ioffset = 0.0;
  double joffset = 0.0;
  double xorigin = dumpdata->gridx[jlower][ilower];
  double yorigin = dumpdata->gridy[jlower][ilower];

  switch (xylocation) {
  case CL_CENTRE:
    ioffset = 0.5;
    joffset = 0.5;
    break;

  case CL_GRID:
    ioffset = 0.0;
    joffset = 0.0;
    break;

  case CL_LEFT:
    ioffset = 0.0;
    joffset = 0.5;
    break;

  case CL_BACK:
    ioffset = 0.5;
    joffset = 0.0;
    break;
  }

  sprintf(line, "rectangular %g %g %d %d %g %g %g %g %g",
          ioffset, joffset, ne1 + 1, ne2 + 1,
          xorigin, yorigin, dumpdata->rg->dx,
          dumpdata->rg->dy, dumpdata->rg->th);
  write_text_att(cdfid, vid, "analytic", line);
}


static void write_analytic_polar_att(dump_data_t *dumpdata, int cdfid,
                                     int vid, int ilower, int ne1,
                                     int jlower, int ne2, int xylocation)
{
  char line[MAXLINELEN];
  double ioffset = 0.0;
  double joffset = 0.0;
  double xdiff = (dumpdata->pg->x0 - dumpdata->gridx[jlower][ilower]);
  double ydiff = (dumpdata->pg->y0 - dumpdata->gridy[jlower][ilower]);
  double rmin = sqrt(xdiff * xdiff + ydiff * ydiff);

  switch (xylocation) {
  case CL_CENTRE:
    ioffset = 0.5;
    joffset = 0.5;
    break;

  case CL_GRID:
    ioffset = 0.0;
    joffset = 0.0;
    break;

  case CL_LEFT:
    ioffset = 0.0;
    joffset = 0.5;
    break;

  case CL_BACK:
    ioffset = 0.5;
    joffset = 0.0;
    break;
  }

  sprintf(line, "polar %g %g %d %d %g %g %g %g %g",
          ioffset, joffset, ne1 + 1, ne2 + 1,
          dumpdata->pg->x0, dumpdata->pg->y0, dumpdata->pg->arc,
          rmin, dumpdata->pg->rotation);
  write_text_att(cdfid, vid, "analytic", line);
}


void write_text_att(int cdfid, int varid, const char *name,
                    const char *text)
{
  nc_put_att_text(cdfid, varid, name, strlen(text), text);
}


static void write_grid_atts(dump_data_t *dumpdata, int fid, int ilower,
                            int jlower)
{
  nc_put_att_text(fid, NC_GLOBAL, "gridtype", strlen(dumpdata->gridtype),
                  dumpdata->gridtype);
  if (dumpdata->rg) {
    double xorigin = dumpdata->gridx[jlower][ilower];
    double yorigin = dumpdata->gridy[jlower][ilower];
    nc_put_att_double(fid, NC_GLOBAL, "xorigin", NC_DOUBLE, 1, &xorigin);
    nc_put_att_double(fid, NC_GLOBAL, "yorigin", NC_DOUBLE, 1, &yorigin);
    nc_put_att_double(fid, NC_GLOBAL, "dx",
                      NC_DOUBLE, 1, &dumpdata->rg->dx);
    nc_put_att_double(fid, NC_GLOBAL, "dy",
                      NC_DOUBLE, 1, &dumpdata->rg->dy);
    nc_put_att_double(fid, NC_GLOBAL, "rotation",
                      NC_DOUBLE, 1, &dumpdata->rg->th);
  }

  else if (dumpdata->pg) {
    double xdiff = (dumpdata->pg->x0 - dumpdata->gridx[jlower][ilower]);
    double ydiff = (dumpdata->pg->y0 - dumpdata->gridy[jlower][ilower]);
    double rmin = sqrt(xdiff * xdiff + ydiff * ydiff);
    nc_put_att_double(fid, NC_GLOBAL, "xorigin",
                      NC_DOUBLE, 1, &dumpdata->pg->x0);
    nc_put_att_double(fid, NC_GLOBAL, "yorigin",
                      NC_DOUBLE, 1, &dumpdata->pg->y0);
    nc_put_att_double(fid, NC_GLOBAL, "arc",
                      NC_DOUBLE, 1, &dumpdata->pg->arc);
    nc_put_att_double(fid, NC_GLOBAL, "rmin", NC_DOUBLE, 1, &rmin);
    nc_put_att_double(fid, NC_GLOBAL, "rotation",
                      NC_DOUBLE, 1, &dumpdata->pg->rotation);
  }
}

void read_grid_atts(parameters_t *params, int fid)
{
  char buf[MAXSTRLEN];

  if (nc_get_att_text(fid, NC_GLOBAL, "gridtype", buf) >= 0)
    strcpy(params->gridtype, buf);

  /* Read info for rectangular grid */
  if (strcasecmp(params->gridtype, "rectangular") == 0) {
    if ((params->rg = (grid_rect_t *)malloc(sizeof(grid_rect_t))) == NULL)
      hd_quit("read_grid_atts: No memory for RectGrid structure\n");
    nc_get_att_double(fid, NC_GLOBAL, "xorigin", &params->rg->x0);
    nc_get_att_double(fid, NC_GLOBAL, "yorigin", &params->rg->y0);
    nc_get_att_double(fid, NC_GLOBAL, "dx", &params->rg->dx);
    nc_get_att_double(fid, NC_GLOBAL, "dy", &params->rg->dy);
    nc_get_att_double(fid, NC_GLOBAL, "rotation", &params->rg->th);
    params->rg->sinth = sin(params->rg->th * M_PI / 180.0);
    params->rg->costh = cos(params->rg->th * M_PI / 180.0);

  } else if (strcasecmp(params->gridtype, "polar") == 0) {
    if ((params->pg =
         (grid_polar_t *)malloc(sizeof(grid_polar_t))) == NULL)
      hd_quit("read_grid_atts: No memory for PolarGrid structure\n");
    nc_get_att_double(fid, NC_GLOBAL, "xorigin", &params->pg->x0);
    nc_get_att_double(fid, NC_GLOBAL, "yorigin", &params->pg->y0);
    nc_get_att_double(fid, NC_GLOBAL, "arc", &params->pg->arc);
    nc_get_att_double(fid, NC_GLOBAL, "rmin", &params->pg->rmin);
    nc_get_att_double(fid, NC_GLOBAL, "rotation", &params->pg->rotation);
  }
}
