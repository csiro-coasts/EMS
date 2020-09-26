/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/misc/solar.c
 *
 *  \brief Solar elevation
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: solar.c 5900 2018-08-28 02:09:26Z riz008 $
 */

/*
 * The following routine was embedded within the heatflux forcing of
 * SHOC. Moving it here means the functionality is accessible by others
 */
#include <stdlib.h>
#include <math.h>
#include "ems.h"

/*-------------------------------------------------------------------*/
/* Routine to calculate the year and Julian day given time units     */
/* Moved from heatflux.c                                             */
/*                                                                   */
/* Depricated: Used only for backwards compatibility and for cases   */
/*             without geographic projection                         */
/*-------------------------------------------------------------------*/
static void dtime_ounits(char *output_tunit, char *model_tunit,
			 double time, int *year, double *day)
{
  int yr = 1990;
  double sf = (time > 0.0) ? 1 : -1; /* should be >=, see dtime_adjlon */
  double d1 = 0.0, dc = 0.0;
  char timeunit[MAXSTRLEN];
  char ounit[MAXSTRLEN];
  char *tok, *tzone, *saveptr;
  double ntime = time;

  /* Find the time zone from the output units string, and put in tzone */
  /*
    FR: This is just a hueristic - strictly speaking the output
        timeunit doesn't have to be in local time, although it usally
        is. To be precise the time should be calculated from the longitude
  */
  strcpy(timeunit, output_tunit); 
  tok = strtok_r ( timeunit, " ", &saveptr );
  while (tok != NULL) {
    tzone = tok;
    tok = strtok_r ( NULL, " ", &saveptr);
  }

  /* Add the time zone to the reference date. */
  strcpy(ounit, "days since 1990-01-01 00:00:00 ");
  if (strlen(tzone) < 5) 
       strcat(ounit, tzone);
  
  /* Change the time units. */
  tm_change_time_units(model_tunit, ounit, &ntime, 1);
  while (dc < sf * ntime) { /* should be <=, see dtime_adjlon */
    if (yr % 4 == 0 && (yr % 100 != 0 || yr % 400 == 0))
      d1 = (sf * 366);
    else
      d1 = (sf * 365);
    dc += d1;
    yr += sf;
  }
  dc -= d1;
  yr -= sf;
  *day = ntime - dc + 1.0;
  *year = yr;
}

/* END dtime_ounits()                                                */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Routine to calculate the year and Julian day given time units     */
/* and with longitude adjustment                                     */
/*-------------------------------------------------------------------*/
static void dtime_adjlon(char *output_tunit, char *model_tunit,
			 double time, int *year, double *day, double lon)
{
  int yr = 1990;
  char ounit[MAXSTRLEN];
  double sf = (time >= 0.0) ? 1 : -1;
  double d1 = 0.0, dc = 0.0;
  double ntime = time;

 /* 
  * Change the time units to reference date UTC 
  */
  strcpy(ounit, "days since 1990-01-01 00:00:00 +0");
  tm_change_time_units(model_tunit, ounit, &ntime, 1);

  /* Algorithm below won't work for time before reference year */
  if (ntime < 1)
    quit("dtime: model time before year 1990 in longitude adjustment");
  
  /* Apply longitude correction */
  if (lon > 180)
    lon -= 360;
  ntime += lon/(24*15);
  
  /* Caluclate day and year based on time */
  sf = (ntime > 0.0) ? 1 : -1;
  while (dc <= sf * ntime) {
    if (yr % 4 == 0 && (yr % 100 != 0 || yr % 400 == 0))
      d1 = (sf * 366);
    else
      d1 = (sf * 365);
    dc += d1;
    yr += sf;
  }
  dc -= d1;
  yr -= sf;
  *day = ntime - dc + 1.0;
  *year = yr;
}

/* END dtime_adjlon()                                                */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Routine to calculate the year and Julian day given time units     */
/*  - gateway function                                               */
/*-------------------------------------------------------------------*/
void dtime(char *output_tunit, char *model_tunit,
	   double time, int *year, double *day, double *lon)
{
  if (lon == NULL)
    dtime_ounits(output_tunit, model_tunit, time, year, day);
  else
    dtime_adjlon(output_tunit, model_tunit, time, year, day, *lon);
}

/**
 * Calculates the solar elevation
 * @param ounit output timeunits
 * @param tunit model timeunits
 * @param time time point
 * @param lat latitude
 * @param out_dec Optionally fills out declination, if not null
 * @param lon longitude of cell for dtime
 * @return solar elevation
 */
double calc_solar_elevation(char *ounit, char *tunit, double time, double lat,
			    double *out_dec, double *lon)
{
  int nday;                     /* Day of the year */
  double hrs;                   /* Hour of the day */
  double dec;                   /* Solar declination */
  double h;                     /* The solar elevation */
  double hrang;                 /* Hour angle */
  double d1;                    /* Dummy variable */
  double jday;                  /* Julian day */
  int yr;                       /* Year */

  if (lon == NULL)
    dtime(ounit, tunit, time, &yr, &jday, NULL);
  else
    dtime(NULL, tunit, time, &yr, &jday, lon);
  
  nday = (int)jday;
  hrs  = 24.0 * (jday - (double)nday);

  d1 = nday * 2 * PI / 365.0;
  dec = 0.006918 + 0.070257 * sin(d1) - 0.399912 * cos(d1)
    + 0.000907 * sin(2 * d1) - 0.006758 * cos(2 * d1)
    + 0.00148 * sin(3 * d1) - 0.002697 * cos(3 * d1);

  /*-----------------------------------------------------------------*/
  /* Get the hour angle */
  hrang = (hrs - 12.0) * 180.0 / 12.0;

  /*-----------------------------------------------------------------*/
  /* Get the solar elevation */
  hrang *= PI / 180.0;
  h = asin(sin(lat) * sin(dec) + cos(lat) * cos(dec) * cos(hrang));
  if (h < 0.0)
    h = 0.0;

  /*
   * Optional output
   */
  if (out_dec != NULL)
    *out_dec = dec;

  return(h);
}
