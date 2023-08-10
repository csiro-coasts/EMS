/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/forcings/heatflux.c
 *  
 *  Description:
 *  Update atmospheric variables as required.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: heatflux.c 6953 2021-12-17 03:14:49Z her127 $
 *
 */

#include <math.h>
#include <string.h>
#include "hd.h"

/* Prototypes                                                        */
double jday(double day);
double jyear(double day);
double julian_dif(double j1, double j2);
double julian_sum(double j1, double days);
double lwrad(double sst, double dtd, double tcld);
double lwrado(double sst);
double swrad(int yr, double jday, double tem, double es, double mdc,
             double lat);
double swr_from_mean(double swr, double lat, double t);
void bulkp(double sst, double dtd, double qs, double qq, double wspd,
           double pres, double zref, double fetch, double *he, double *la,
           int bulkf);
int mday(int nday, int year);
double dewe(double tem, double pres);
double estt(double ta, double pr, double ar);

double masagutov(double wsd, double tv, double ts, double zref);
double interp_maov(double a, double b, int mode);
void init_fetch(master_t *master, unsigned long **flag);
double get_fetch(double *fetch, double wdir);
double get_vapour_press(window_t *windat, win_priv_t *wincon, int c, double at, 
			double pres, double sal, double *es, double *esat, double *rh);
double get_albedo(int yr, double jday, double mdc, double lat);
double advt(double ts, double ta, double u, double x, double K,
            double zref);
int sign(double x);

void surf_heat_flux(geometry_t *window, window_t *windat, win_priv_t *wincon);
void calc_heatflux(geometry_t *window, window_t *windat, win_priv_t *wincon);
void surf_salt_flux(geometry_t *window, window_t *windat, win_priv_t *wincon);
void bulk_salt_flux(geometry_t *window, window_t *windat, win_priv_t *wincon);
void comp_heat_mom(geometry_t *window, window_t *windat, win_priv_t *wincon);
void comp_heat_inv(geometry_t *window, window_t *windat, win_priv_t *wincon);
void comp_heat_none(geometry_t *window, window_t *windat, win_priv_t *wincon);

int dayno[2][13] = {
  {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
  {0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}
};

double ad[5] = { 0.0, 0.771, 0.867, 1.2, 0.0 };
double bd[5] = { 1.08, 0.0858, 0.0667, 0.025, 0.073 };
double pd[5] = { -0.15, 1.0, 1.0, 1.0, 1.0 };

double ah[5] = { 0.0, 0.927, 1.15, 1.17, 1.652 };
double bh[5] = { 1.185, 0.0546, 0.01, 0.0075, -0.017 };
double ch[5] = { 0.0, 0.0, 0.0, -0.00045, 0.0 };
double ph[5] = { -0.157, 1.0, 1.0, 1.0, 1.0 };

double ae[5] = { 0.0, 0.969, 1.18, 1.196, 1.68 };
double be[5] = { 1.23, 0.0521, 0.01, 0.008, -0.016 };
double ce[5] = { 0.0, 0.0, 0.0, -0.0004, 0.0 };
double pe[5] = { -0.16, 1.0, 1.0, 1.0, 1.0 };

double ap[9][10] =
  { {0.25, 0.25, 0.13, 0.08, 0.06, 0.05, 0.045, 0.04, 0.04, 0.04},
{0.25, 0.25, 0.13, 0.08, 0.06, 0.05, 0.045, 0.04, 0.04, 0.04},
{0.25, 0.25, 0.13, 0.08, 0.06, 0.05, 0.045, 0.04, 0.04, 0.04},
{0.25, 0.25, 0.13, 0.08, 0.06, 0.05, 0.045, 0.04, 0.04, 0.04},
{0.25, 0.25, 0.13, 0.08, 0.06, 0.05, 0.045, 0.04, 0.04, 0.04},
{0.22, 0.22, 0.12, 0.08, 0.065, 0.055, 0.05, 0.045, 0.045, 0.045},
{0.19, 0.19, 0.11, 0.08, 0.065, 0.055, 0.05, 0.05, 0.05, 0.05},
{0.14, 0.14, 0.09, 0.075, 0.068, 0.063, 0.06, 0.06, 0.06, 0.06},
{0.09, 0.09, 0.075, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07}
};

double bunker_Cd[10][7] =
  { {0.06, 0.06, 0.98, 1.20, 1.32, 1.56, 1.80},
    {0.77, 1.30, 1.43, 1.54, 1.60, 1.78, 1.86},
    {1.47, 1.72, 1.80, 1.87, 1.90, 2.00, 2.10},
    {1.95, 2.04, 2.10, 2.16, 2.22, 2.25, 2.32},
    {2.26, 2.30, 2.35, 2.40, 2.42, 2.44, 2.48},
    {2.52, 2.54, 2.57, 2.60, 2.62, 2.63, 2.64},
    {2.78, 2.79, 2.80, 2.80, 2.80, 2.80, 2.80},
    {3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00},
    {3.20, 3.20, 3.20, 3.20, 3.20, 3.20, 3.20},
    {3.40, 3.40, 3.40, 3.40, 3.40, 3.40, 3.40}
  };

double bunker_Ce[9][7] =
  { {0.07, 0.30, 0.72, 1.32, 1.65, 2.05, 2.52},
    {0.22, 0.67, 1.12, 1.34, 1.45, 1.68, 2.01},
    {0.69, 1.17, 1.36, 1.44, 1.46, 1.58, 1.79},
    {1.06, 1.36, 1.48, 1.53, 1.58, 1.65, 1.79},
    {1.39, 1.58, 1.61, 1.64, 1.68, 1.74, 1.84},
    {1.59, 1.68, 1.75, 1.80, 1.82, 1.86, 1.94},
    {1.74, 1.79, 1.83, 1.86, 1.86, 1.86, 1.93},
    {1.81, 1.84, 1.85, 1.86, 1.87, 1.88, 1.90},
    {1.86, 1.86, 1.86, 1.86, 1.86, 1.86, 1.86}
  };

double maov[9][22] =
  { {0.42, 0.42, 0.42, 0.42, 0.52, 0.72, 0.86, 0.98, 1.06, 1.13, 1.18,
     1.23, 1.28, 1.30, 1.34, 1.37, 1.40, 1.44, 1.46, 1.48, 1.51, 1.54},
{0.44, 0.44, 0.44, 0.58, 0.70, 0.84, 0.93, 1.03, 1.11, 1.18, 1.22, 1.26,
 1.29, 1.32, 1.35, 1.38, 1.42, 1.45, 1.47, 1.49, 1.52, 1.55},
{0.68, 0.68, 0.68, 0.87, 0.95, 1.00, 1.06, 1.11, 1.16, 1.20, 1.25, 1.28,
 1.30, 1.34, 1.37, 1.40, 1.43, 1.46, 1.48, 1.50, 1.53, 1.55},
{0.95, 0.95, 0.95, 1.04, 1.08, 1.10, 1.14, 1.18, 1.22, 1.25, 1.28, 1.30,
 1.32, 1.35, 1.38, 1.41, 1.44, 1.47, 1.49, 1.51, 1.54, 1.56},
{1.34, 1.34, 1.32, 1.29, 1.26, 1.24, 1.23, 1.24, 1.27, 1.29, 1.31, 1.33,
 1.36, 1.38, 1.40, 1.43, 1.44, 1.47, 1.49, 1.51, 1.54, 1.56},
{1.94, 1.71, 1.54, 1.42, 1.36, 1.30, 1.29, 1.29, 1.30, 1.31, 1.32, 1.34,
 1.36, 1.38, 1.40, 1.43, 1.44, 1.47, 1.49, 1.51, 1.54, 1.56},
{1.95, 1.74, 1.57, 1.46, 1.38, 1.33, 1.32, 1.31, 1.32, 1.32, 1.34, 1.36,
 1.37, 1.39, 1.41, 1.44, 1.45, 1.48, 1.50, 1.52, 1.54, 1.56},
{1.95, 1.76, 1.60, 1.49, 1.42, 1.36, 1.34, 1.34, 1.34, 1.34, 1.35, 1.36,
 1.38, 1.40, 1.42, 1.44, 1.46, 1.48, 1.50, 1.52, 1.54, 1.56},
{1.96, 1.79, 1.64, 1.52, 1.44, 1.39, 1.36, 1.36, 1.36, 1.36, 1.37, 1.38,
 1.40, 1.42, 1.44, 1.46, 1.48, 1.50, 1.52, 1.54, 1.56, 1.58}
};
double maod[9][22] =
  { {0.0, 0.0, 0.0, 0.55, 0.66, 0.72, 0.78, 0.84, 0.92, 1.02, 1.10, 1.20,
     1.26, 1.35, 1.42, 1.50, 1.57, 1.64, 1.72, 1.78, 1.86, 1.93},
{0.40, 0.40, 0.40, 0.64, 0.72, 0.78, 0.83, 0.90, 0.98, 1.06, 1.13, 1.22,
 1.29, 1.37, 1.44, 1.51, 1.58, 1.65, 1.73, 1.79, 1.86, 1.93},
{0.54, 0.54, 0.54, 0.76, 0.80, 0.84, 0.88, 0.94, 1.00, 1.08, 1.16, 1.24,
 1.31, 1.39, 1.46, 1.53, 1.59, 1.66, 1.73, 1.79, 1.86, 1.93},
{0.81, 0.81, 0.81, 0.88, 0.90, 0.92, 0.95, 1.00, 1.05, 1.12, 1.19, 1.27,
 1.34, 1.40, 1.48, 1.54, 1.60, 1.66, 1.74, 1.80, 1.86, 1.93},
{1.08, 1.08, 1.08, 1.04, 1.01, 0.99, 1.02, 1.06, 1.10, 1.17, 1.24, 1.30,
 1.36, 1.42, 1.48, 1.54, 1.60, 1.66, 1.74, 1.80, 1.86, 1.93},
{1.66, 1.40, 1.22, 1.12, 1.07, 1.04, 1.06, 1.10, 1.14, 1.20, 1.26, 1.32,
 1.38, 1.44, 1.50, 1.56, 1.62, 1.68, 1.74, 1.80, 1.86, 1.93},
{1.69, 1.45, 1.26, 1.15, 1.09, 1.07, 1.08, 1.11, 1.15, 1.21, 1.27, 1.33,
 1.38, 1.45, 1.51, 1.57, 1.63, 1.68, 1.74, 1.80, 1.87, 1.93},
{1.72, 1.50, 1.30, 1.17, 1.12, 1.09, 1.10, 1.13, 1.17, 1.22, 1.28, 1.34,
 1.39, 1.45, 1.52, 1.58, 1.63, 1.69, 1.75, 1.81, 1.87, 1.93},
{1.74, 1.54, 1.34, 1.20, 1.14, 1.11, 1.12, 1.14, 1.18, 1.24, 1.29, 1.34,
 1.40, 1.46, 1.52, 1.58, 1.64, 1.69, 1.76, 1.81, 1.87, 1.93}
};

typedef struct {
  double dt;                    /* Time step */
  timeseries_t *ts;             /* Timeseries file */
  int id;                       /* TS id for the variable */
} heatf_data_t;

typedef struct {
  double dt;                    /* Time step */
  timeseries_t *ts;             /* Timeseries file */
  int id;                       /* TS id for the variable */
} latent_data_t;
typedef struct {
  double dt;                    /* Time step */
  timeseries_t *ts;             /* Timeseries file */
  int id;                       /* TS id for the variable */
} sensible_data_t;
typedef struct {
  double dt;                    /* Time step */
  timeseries_t *ts;             /* Timeseries file */
  int id;                       /* TS id for the variable */
} longin_data_t;
typedef struct {
  double dt;                    /* Time step */
  timeseries_t *ts;             /* Timeseries file */
  int id;                       /* TS id for the variable */
} longout_data_t;


/*-------------------------------------------------------------------*/
/* Functions for reading the schedule the heatflux forcings.         */
/*-------------------------------------------------------------------*/
int heatf_init(sched_event_t *event)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  parameters_t *params = master->params;
  heatf_data_t *data = NULL;

  /* Read parameters */
  prm_set_errfn(hd_silent_warn);

  data = (heatf_data_t *)malloc(sizeof(heatf_data_t));
  schedSetPrivateData(event, data);
  data->ts = frcw_read_cell_ts(master, params->hf, params->hf_dt,
                              "heatflux", "W m-2", &data->dt, &data->id,
                              &master->heatf);

  if (data->ts == NULL) {
    free(data);
    return 0;
  }

  return 1;
}

/* END heatf_init()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read the net heatflux from file.                       */
/*-------------------------------------------------------------------*/
double heatf_event(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  heatf_data_t *data = (heatf_data_t *)schedGetPrivateData(event);

  if (t >= (event->next_event - SEPS)) {
    /* Convert surface heat flux ms-1K */
    frc_ts_eval_grid(master, t, data->ts, data->id, master->heatf,
                     -1.0 / (4.0e3 * 1025.0));
    event->next_event += data->dt;
  }

  return event->next_event;
}

/* END heatf_event()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Clean up the net heat flux                                        */
/*-------------------------------------------------------------------*/
void heatf_cleanup(sched_event_t *event, double t)
{
  heatf_data_t *data = (heatf_data_t *)schedGetPrivateData(event);

  if (data != NULL) {
    if (master->heatf) free(master->heatf);
    hd_ts_free(master, data->ts);
    free(data);
  }
}

/* END heatf_cleanup()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Functions for reading the schedule the heatflux forcings where    */
/* the net heatflux is assembled from heatflux components.           */
/* Longwave output input is initialised here.                        */
/*-------------------------------------------------------------------*/
int longout_init(sched_event_t *event)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  parameters_t *params = master->params;
  longout_data_t *data = NULL;

  /* Read parameters */
  prm_set_errfn(hd_silent_warn);

  data = (longout_data_t *)malloc(sizeof(longout_data_t));
  schedSetPrivateData(event, data);

  if ((data->ts = frcw_read_cell_ts(master, params->hf, params->hf_dt,
				   "lwr", "W m-2", &data->dt, &data->id,
				   &master->lwro)) == NULL )
    data->ts = frc_read_cell_ts(master, params->hf, params->hf_dt,
				"lwr", "W m-2", &data->dt, &data->id,
				&master->lwro);

  if (data->ts == NULL) {
    free(data);
    return 0;
  } else
    master->heatflux |= COMP_LWO;

  return 1;
}

/* END longout_init()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read the longwave output from file.                    */
/*-------------------------------------------------------------------*/
double longout_event(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  /*UR unused
   * geometry_t *geom = master->geom; */
  longout_data_t *data = (longout_data_t *)schedGetPrivateData(event);

  if (t >= (event->next_event - SEPS)) {
    /* Convert surface heat flux ms-1K */
    frc_ts_eval_grid(master, t, data->ts, data->id, master->lwro,
                     -1.0 / (4.0e3 * 1025.0));
    event->next_event += data->dt;
  }

  return event->next_event;
}

/* END longout_event()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Clean up the longwave output data                                 */
/*-------------------------------------------------------------------*/
void longout_cleanup(sched_event_t *event, double t)
{
  longout_data_t *data = (longout_data_t *)schedGetPrivateData(event);

  if (data != NULL) {
    free(master->lwro);
    hd_ts_free(master, data->ts);
    free(data);
  }
}

/* END longout_cleanup()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Functions for reading the schedule the heatflux forcings where    */
/* the net heatflux is assembled from heatflux components.           */
/* Longwave input input is initialised here.                         */
/*-------------------------------------------------------------------*/
int longin_init(sched_event_t *event)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  parameters_t *params = master->params;
  longin_data_t *data = NULL;

  /* Read parameters */
  prm_set_errfn(hd_silent_warn);

  data = (longin_data_t *)malloc(sizeof(longin_data_t));
  schedSetPrivateData(event, data);
  if (params->heatflux & ADVANCED) 
    data->ts = frc_read_cell_ts(master, params->lwri, params->lwri_dt,
				"lwr_in", "W m-2", &data->dt, &data->id,
				&master->lwri);
  else
    data->ts = frc_read_cell_ts(master, params->hf, params->hf_dt,
				"lwr_in", "W m-2", &data->dt, &data->id,
				&master->lwrd);

  if (data->ts == NULL) {
    free(data);
    return 0;
  } else
    master->heatflux |= COMP_LWI;

  return 1;
}

/* END longin_init()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read the longwave input from file. This is then added  */
/* to the net heat array.                                            */
/*-------------------------------------------------------------------*/
double longin_event(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  longin_data_t *data = (longin_data_t *)schedGetPrivateData(event);

  if (t >= (event->next_event - 2 * SEPS)) {
    /* Convert surface heat flux ms-1K */
    if (params->heatflux & ADVANCED)
      frc_ts_eval_grid(master, t, data->ts, data->id, master->lwri, 1.0);		       
    else
      frc_ts_eval_grid(master, t, data->ts, data->id, master->lwrd,
		       -1.0 / (4.0e3 * 1025.0));

    event->next_event += data->dt;
  }

  return event->next_event;
}

/* END longin_event()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Clean up the longwave input data                                  */
/*-------------------------------------------------------------------*/
void longin_cleanup(sched_event_t *event, double t)
{
  longin_data_t *data = (longin_data_t *)schedGetPrivateData(event);

  if (data != NULL) {
    if (params->heatflux & ADVANCED)
      free(master->lwri);
    else
      free(master->lwrd);
    hd_ts_free(master, data->ts);
    free(data);
  }
}

/* END longin_cleanup()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Functions for reading the schedule the heatflux forcings where    */
/* the net heatflux is assembled from heatflux components.           */
/* Sensible heat input is initialised here.                          */
/*-------------------------------------------------------------------*/
int sheat_init(sched_event_t *event)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  parameters_t *params = master->params;
  sensible_data_t *data = NULL;

  /* Read parameters */
  prm_set_errfn(hd_silent_warn);

  data = (sensible_data_t *)malloc(sizeof(sensible_data_t));
  schedSetPrivateData(event, data);
  if ((data->ts = frcw_read_cell_ts(master, params->hf, params->hf_dt,
				"sensible", "W m-2", &data->dt, &data->id,
				   &master->shfd)) == NULL)
    data->ts = frc_read_cell_ts(master, params->hf, params->hf_dt,
				"sensible", "Wm-2", &data->dt, &data->id,
				&master->shfd);

  if (data->ts == NULL) {
    free(data);
    return 0;
  }

  return 1;
}

/* END sheat_init()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read the sensible heat input from file. This is then   */
/* added to the net heat array.                                      */
/*-------------------------------------------------------------------*/
double sheat_event(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  sensible_data_t *data = (sensible_data_t *)schedGetPrivateData(event);

  if (t >= (event->next_event - 3 * SEPS)) {
    /* Convert surface heat flux ms-1K */
    frc_ts_eval_grid(master, t, data->ts, data->id, master->shfd,
                     -1.0 / (4.0e3 * 1025.0));
    event->next_event += data->dt;
  }

  return event->next_event;
}

/* END sheat_event()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Clean up the sensible heat flux input data                        */
/*-------------------------------------------------------------------*/
void sheat_cleanup(sched_event_t *event, double t)
{
  sensible_data_t *data = (sensible_data_t *)schedGetPrivateData(event);

  if (data != NULL) {
    free(master->shfd);
    hd_ts_free(master, data->ts);
    free(data);
  }
}

/* END sheat_cleanup()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Functions for reading the schedule the heatflux forcings where    */
/* the net heatflux is assembled from heatflux components.           */
/* Latent heat input is initialised here.                            */
/*-------------------------------------------------------------------*/
int lheat_init(sched_event_t *event)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  parameters_t *params = master->params;
  latent_data_t *data = NULL;

  /* Read parameters */
  prm_set_errfn(hd_silent_warn);

  data = (latent_data_t *)malloc(sizeof(latent_data_t));
  schedSetPrivateData(event, data);
  if (params->heatflux & COMP_HEAT_MOM)
    data->ts = frc_read_cell_ts(master, params->hf, params->hf_dt,
				"latent", "kg/m2/sec", &data->dt, &data->id,
				&master->lhfd);
  else
    data->ts = frcw_read_cell_ts(master, params->hf, params->hf_dt,
				"latent", "W m-2", &data->dt, &data->id,
				&master->lhfd);

  if (data->ts == NULL) {
    free(data);
    return 0;
  }

  return 1;
}

/* END lheat_init()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read the latent heat input from file. This is then     */
/* truncated to a maximum of zero to remove condensation effects and */
/* added to the net heat array.                                      */
/*-------------------------------------------------------------------*/
double lheat_event(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  latent_data_t *data = (latent_data_t *)schedGetPrivateData(event);

  if (t >= (event->next_event - 4 * SEPS)) {
    /* Convert surface heat flux ms-1K */
    frc_ts_eval_grid(master, t, data->ts, data->id, master->lhfd,
                     -1.0 / (4.0e3 * 1025.0));
    event->next_event += data->dt;
  }

  return event->next_event;
}

/* END lheat_event()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Clean up the latent heat flux input data                          */
/*-------------------------------------------------------------------*/
void lheat_cleanup(sched_event_t *event, double t)
{
  latent_data_t *data = (latent_data_t *)schedGetPrivateData(event);

  if (data != NULL) {
    free(master->lhfd);
    hd_ts_free(master, data->ts);
    free(data);
  }
}

/* END lheat_cleanup()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Functions for reading the schedule the heatflux forcings where    */
/* the net heatflux is assembled from MOM compatibel heatflux        */
/* components. The shortwave is an integrated value over swr_dt.     */
/*-------------------------------------------------------------------*/
int swr_mean_init(sched_event_t *event)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  parameters_t *params = master->params;
  longout_data_t *data = NULL;

  /* Read parameters */
  prm_set_errfn(hd_silent_warn);

  data = (longout_data_t *)malloc(sizeof(longout_data_t));
  schedSetPrivateData(event, data);
  if ((data->ts = frcw_read_cell_ts(master, params->hf, params->hf_dt,
				    "swr", "Wm-2", &data->dt, &data->id,
				    &master->swrd)) == NULL)
    data->ts = frc_read_cell_ts(master, params->hf, params->hf_dt,
				"swr", "W m-2", &data->dt, &data->id,
				&master->swrd);

  master->albedo = max(0.0, master->albedo);
  if (data->ts == NULL) {
    free(data);
    return 0;
  }

  return 1;
}

/* END swr_mean_init()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read the mean swr output from file.                    */
/*-------------------------------------------------------------------*/
double swr_mean_event(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  geometry_t *geom = master->geom;
  longout_data_t *data = (longout_data_t *)schedGetPrivateData(event);
  double day_dt = 86400.0;
  double start = floor(t) + (10 * 86400);
  int na = (int)(86400 / data->dt);
  int c, cc, n;

  if (t >= (event->next_event - SEPS)) {
    if (params->heatflux & COMP_HEAT_NONE) {
      frc_ts_eval_grid(master, t, data->ts, data->id, master->swr, 1.0);
      for (cc = 1; cc <= geom->b2_t; cc++) {
	c = geom->w2_t[cc];
	master->swrd[c] = master->swr[c];
      }
      event->next_event += data->dt;
    } else {
      memset(master->swrd, 0, geom->sgsizS * sizeof(double));
      for (n = 0; n < na; n++) {
	/* Convert surface heat flux ms-1K */
	t = start + (double)n * data->dt;
	frc_ts_eval_grid(master, t, data->ts, data->id, master->swr, 1.0);
	for (cc = 1; cc <= geom->b2_t; cc++) {
	  c = geom->w2_t[cc];
	  master->swrd[c] += (master->swr[c] / (double)na);
	}
      }
      event->next_event += day_dt;
    }
  }
  return event->next_event;
}

/* END swr_mean_event()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Clean up the mean swr output data                                 */
/*-------------------------------------------------------------------*/
void swr_mean_cleanup(sched_event_t *event, double t)
{
  longout_data_t *data = (longout_data_t *)schedGetPrivateData(event);

  if (data != NULL) {
    free(master->swrd);
    hd_ts_free(master, data->ts);
    free(data);
  }
}

/* END swr_mean_cleanup()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to perform the selected heatflux                          */
/*-------------------------------------------------------------------*/
void calc_heatf(geometry_t *window,
		window_t   *windat, 
		win_priv_t *wincon
  )
{
  double fact = -4.0e3 * 1025.0;  /* Conversion Wm-2 to ms-1K */
  int dof = 0;

  /* Surface relaxation is done in the windows */
  if (!(wincon->heatflux & (ADVANCED | INVERSE | NET_HEAT | COMP_HEAT | COMP_HEAT_MOM | COMP_HEAT_NONE)))
    return;

  /* No heatflux before the ramp time */
  if (windat->t < wincon->hf_ramp) {
    if (windat->heatf)
      memset(windat->heatf, 0, window->sgsizS * sizeof(double));
    if (windat->swr)
      memset(windat->swr, 0, window->sgsizS * sizeof(double));
    return;
  }

  /* Get the heat flux for ADVANCED and INVERSE methods. Time series */
  /* of heatflux is read from file via the scheduler if a file is */
  /* present in the input file.  */
  if (wincon->heatflux & ADVANCED) {
    surf_heat_flux(window, windat, wincon);
  } else if (wincon->heatflux & INVERSE) {
    calc_heatflux(window, windat, wincon);
  } else if (wincon->heatflux & NET_HEAT) {
    /* Nothing to do : heatflux read from file via the scheduler */
    /* Copy heatfluxs to the diagnostic tracers if required      */
    if (windat->heatf) {
      int c, cc;
      for (cc = 1; cc <= window->b2_t; cc++) {
	c = window->w2_t[cc];
	windat->nhfd[c] = windat->heatf[c] * fact;
      }
    }
    if (windat->swr) {
      int c, cc;
      for (cc = 1; cc <= window->b2_t; cc++) {
	c = window->w2_t[cc];
	windat->swrd[c] = windat->swr[c];
      }
    }
  } else if (wincon->heatflux & COMP_HEAT) {
    comp_heat_inv(window, windat, wincon);
  } else if (wincon->heatflux & COMP_HEAT_MOM) {
    comp_heat_mom(window, windat, wincon);
  } else if (wincon->heatflux & COMP_HEAT_NONE) {
    comp_heat_none(window, windat, wincon);
  }
  if (dof) {
  switch (wincon->heatflux) {
  case ADVANCED:
    surf_heat_flux(window, windat, wincon);
    break;
  case INVERSE:
    calc_heatflux(window, windat, wincon);
    break;
  case NET_HEAT:
    /* Nothing to do : heatflux read from file via the scheduler */
    /* Copy heatfluxs to the diagnostic tracers if required      */
    if (windat->heatf) {
      int c, cc;
      for (cc = 1; cc <= window->b2_t; cc++) {
	c = window->w2_t[cc];
	windat->nhfd[c] = windat->heatf[c] * fact;
      }
    }
    if (windat->swr) {
      int c, cc;
      for (cc = 1; cc <= window->b2_t; cc++) {
	c = window->w2_t[cc];
	windat->swrd[c] = windat->swr[c];
      }
    }
    break;
  case COMP_HEAT:
    comp_heat_inv(window, windat, wincon);
    break;
  case COMP_HEAT_MOM:
    comp_heat_mom(window, windat, wincon);
    break;
  case COMP_HEAT_NONE:
    comp_heat_none(window, windat, wincon);
    break;
  }
  }
}

/* END calc_heatf()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to perform the selected saltflux                          */
/*-------------------------------------------------------------------*/
void calc_saltf(geometry_t *window,
		window_t   *windat, 
		win_priv_t *wincon
  )
{

  /* Surface relaxation is done in the windows */
  if (wincon->saltflux & NONE)
    return;

  /* Get the salt flux for ADVANCED & ORIGINAL methods. Time series */
  /* of saltflux is read from file via the scheduler if a file is */
  /* present in the input file.  */
  switch (wincon->saltflux) {
  case ADVANCED:
    surf_salt_flux(window, windat, wincon);
    break;
  case BULK:
    bulk_salt_flux(window, windat, wincon);
    break;
  case ORIGINAL:
    surf_salt_flux(window, windat, wincon);
    break;
  }
}

/* END calc_saltf()                                                  */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Routine to transfer heat components as is                         */
/*-------------------------------------------------------------------*/
void comp_heat_none(geometry_t *window,
		    window_t   *windat, 
		    win_priv_t *wincon)
{
  int c, cc;
  double factH = -4.0e3 * 1025.0;   /* Conversion Wm-2 to ms-1K */
  double factS = 1000 * 86400;      /* Conversion m/s to mm/day */
  
  /*
   * Note: For single window some of the memory locations will be
   *       the same on both sides
   */

  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->w2_t[cc];
    if (windat->swrd)
      windat->swr[c] = windat->swrd[c];
    if (windat->lwro) {
      /* See conversion factor in longout_event */
      windat->tr_wcS[windat->lwrn][c] = factH * windat->lwro[c];
    }
    if (windat->shfd) {
      windat->tr_wcS[windat->shfn][c] = factH * windat->shfd[c];
    }
    if (windat->lhfd) {
      windat->tr_wcS[windat->lhfn][c] = factH * windat->lhfd[c];
    }
    if (windat->evap) {
      /* evap_event converts mm/day to m/s - here we reverse that */
      windat->tr_wcS[windat->evapn][c] = factS * windat->evap[c];
    }
    if (windat->precip) {
      /* precip_event converts mm/day to m/s - here we reverse that */
      windat->tr_wcS[windat->precipn][c] = factS * windat->precip[c];
    }
  }
}

/*-------------------------------------------------------------------*/
/* Routine to add the components of the net heat flux.               */
/*-------------------------------------------------------------------*/
void comp_heat_mom(geometry_t *window,
		   window_t   *windat, 
		   win_priv_t *wincon
  )
{
  double fact = -4.0e3 * 1025.0;         /* Conversion Wm-2 to ms-1K */
  double hlv;                            /* Latent heat of evaporation (J/kg) */
  double ang = 7.29e-5;                  /* Earth's angular velocity (s-1) */
  double sc = 1.0;                       /* Sign convention scale */
  double lhf;                            /* Latent heat factor */
  double lat;                            /* Latitude (deg) */
  int c, cc;                             /* Counters                 */

  /* hlv = 4.1868 * (597.31 - 0.56525 * at) * 1e3                    */
  /* MOM : set in shared/constants/constants.F90. Corresponds to     */
  /* ~.35 deg C.                                                     */
  hlv = -2.5e6;
  /* Corresponds to ~25 deg C                                        */
  hlv = -2.44e6;

  lhf = fact * hlv;

  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->w2_t[cc];
 
    lat = asin(wincon->coriolis[c] / (2.0 * ang));

    /*
     * The handling of swr/swrd is a bit confusing here
     * - The scheduler controls the allocation of swrd
     *   BUT fills in swr and then computes and stores into swrd
     * - swr is a tracer which was temporarily used by the
     *   swr_mean_event to fill swrd
     * - Finally the tracer, swr, is filled on the line below
     */
    windat->swr[c] = (1.0 - wincon->albedo) *
      swr_from_mean(windat->swrd[c], lat, windat->t);

    windat->heatf[c] = (windat->lwro[c] + 
			windat->shfd[c] + windat->lhfd[c] * hlv);
    windat->nhfd[c] = windat->heatf[c] * fact + windat->swr[c];

    /* Note : the diagnostic tracers for lwr, shf and lhf were not   */
    /* pointed to by master->lwr etc. in load_tracer() since the     */
    /* master arrays are used to read in the data from file (i.e.    */
    /* they contain different data to the diagnostics. Hence the use */
    /* of the tracer numbers master->lwrn etc, set in load_tracer(). */
    windat->tr_wcS[windat->lwrn][c] = windat->lwro[c] * fact;
    windat->tr_wcS[windat->lhfn][c] = windat->lhfd[c] * lhf;
    windat->tr_wcS[windat->shfn][c] = windat->shfd[c] * fact;
  }
}

/* END comp_heat_mom()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to add the components of the net heat flux, with longwave */
/* input computed assuming clear skies.                              */
/*-------------------------------------------------------------------*/
void comp_heat_inv(geometry_t *window,
		   window_t   *windat, 
		   win_priv_t *wincon
  )
{
  int c, cc;                    /* Counters */
  double at;                    /* Dry bulb temperature (deg C) */
  double dtw;                   /* Wet bulb temperature (deg C) */
  double tcld;                  /* Fractional cloud cover */
  double wt;                    /* Water temperature (deg C) */
  double sal;                   /* Salinity (psu) */
  double wspd;                  /* Windspeed (ms-1) */
  double wdir;                  /* Wind direction (deg T) */
  double pres;                  /* Air pressure (HPa) */
  double es;                    /* Vapour pressure at the air temp (HPa) */
  double esat;                  /* Saturation vapour pressure (HPa) */
  double ew;                    /* Vapour pressure at water temp. (Hpa) */
  double rh;                    /* Relative humidity (%) */
  double lg, lwr;               /* Long wave radiation flux (Wm-2) */
  double sg;                    /* Short wave radiation flux (Wm-2) */
  double lat;                   /* Latitude (deg) */
  int yr;                       /* Year */
  double day;                   /* Julian day */
  double dt;                    /* Time step (s) */
  double lhf;                   /* Latent heat factor */
  double sc = 1.0;              /* Sign convention scale */
  double Cv = 4e3;              /* Specific heat at constant volume */
  double ang = 7.29e-5;         /* Earth's angular velocity (s-1) */
  double fact = -4.0e3*1025.0;  /* Conversion Wm-2 to ms-1K */

  /* Get the year and Julian day */
  /*dtime(master->t / 86400, &yr, &day);*/
  dt = windat->dt;
  es = sg = 0.0;
  lhf = fact;

  /* Loop over the water cells */
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->w2_t[cc];

    wt = windat->temp[c];
    sal = windat->sal[c];
    at = (windat->airtemp) ? windat->airtemp[c] : 20.0;
    pres = windat->patm[c] / 100.0; /* Convert pressure to HPa */
    wspd = windat->windspeed[c];
    wdir = windat->winddir[c];

    /* Get the wet bulb temperature (not currently used). */
    if (windat->rh) {
      rh = windat->rh[c];
      dtw = estt(at, pres, rh);
      /* Saturation vapour pressure in the air.  */
      esat =
	dtw * (dtw *
	       (dtw *
		(dtw *
		 (dtw * (dtw * 6.136821e-11 + 2.034081e-8) +
		  3.03124e-6) + 2.650648e-4) + 1.428946e-2) +
	       0.4436519) + 6.1078;
      esat *= (1.0 - 0.000537 * sal); /* esat at salinity sal */
      /* Actual vapour pressure in the air.  */
      es = esat - (6.6e-4 * (1 + 1.5e-3 * dtw) * (at - dtw) * pres);
    } /* else {
      hd_quit_and_dump
	("Realtive humidity must be supplied as input.\n");
    }
      */

    /* Get the clear sky long wave radiation */
    /* Computed incoming longwave assuming clear skies (Wm-2) */
    lg = lwrad(wt, at, 0.0) - lwrado(windat->temp[c]);
    /* Incoming longwave read from file, taking clouds into account (Wm-2) */
    lwr = windat->lwrd[c] * fact;
    /* Comment this to use computed incoming longwave */
    /*lg = lwr;*/
    /* Total cloud inferred from the above (not currently used) */
    /* tcld = min((1.0 - lwr / lg) / 0.63, 1.0); */

    /* Get the short wave radiation (not currently used) */
    /*
    lat = asin(master->coriolis[c] / (2.0 * ang));
    sg = swrad(yr, day, at, es, tcld, lat);
    */
    
    /* Get the net heat flux */
    windat->heatf[c] = (lg / fact + windat->lwro[c] +
			windat->shfd[c] + max(windat->lhfd[c], 0.0));
    windat->nhfd[c] = windat->heatf[c] * fact;
    windat->swrd[c] = windat->swr[c];
    /* Note : the diagnostic tracers for lwr, shf and lhf were not   */
    /* pointed to by master->lwr etc. in load_tracer() since the     */
    /* master arrays are used to read in the data from file (i.e.    */
    /* they contain different data to the diagnostics. Hence the use */
    /* of the tracer numbers master->lwrn etc, set in load_tracer(). */
    windat->tr_wcS[windat->lwrn][c] = lg + windat->lwro[c] * fact;
    windat->tr_wcS[windat->lhfn][c] = (max(windat->lhfd[c], 0.0)) * lhf;
    windat->tr_wcS[windat->shfn][c] = windat->shfd[c] * fact * sc;
  }
}

/* END comp_heat_inv()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to implement surface heat flux terms                      */
/*-------------------------------------------------------------------*/
void surf_heat_flux(geometry_t *window,
		    window_t   *windat, 
		    win_priv_t *wincon
  )
{
  int c, cc;                    /* Counters */
  double at;                    /* Dry bulb temperature (deg C) */
  double dtw;                   /* Wet bulb temperature (deg C) */
  double tcld;                  /* Fractional cloud cover */
  double wt;                    /* Water temperature (deg C) */
  double sal;                   /* Salinity (psu) */
  double wspd;                  /* Windspeed (ms-1) */
  double wdir;                  /* Wind direction (deg T) */
  double pres;                  /* Air pressure (HPa) */
  double es;                    /* Vapour pressure at the air temp (HPa) */
  double esat;                  /* Saturation vapour pressure (HPa) */
  double ew;                    /* Vapour pressure at water temp. (Hpa) */
  double rh;                    /* Relative humidity (%) */
  double q;                     /* Specific humidity at pressure es
                                   (kg/kg) */
  double qs;                    /* Specific humidity at the surface
                                   (kg/kg) */
  double bow;                   /* Bowen ratio */
  double lg;                    /* Long wave radiation flux (Wm-2) */
  double sg;                    /* Short wave radiation flux (Wm-2) */
  double he;                    /* Sensible heat flux (Wm-2) */
  double la;                    /* Latent heat flux (Wm-2) */
  double lat;                   /* Latitude (deg) */
  double fetch;                 /* Wind fetch (m) */
  int yr;                       /* Year */
  double day;                   /* Julian day */
  double dt;                    /* Time step (s) */
  double Cv = 4e3;              /* Specific heat at constant volume */
  double ang = 7.29e-5;         /* Earth's angular velocity (s-1) */

  dt = windat->dt;
  es = sg = 0.0;

  /* Loop over the water cells */
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->w2_t[cc];

    /* Get the year and Julian day */
    if (window->is_geog)
      dtime(NULL, master->timeunit, windat->t, &yr, &day, &window->cellx[c]);
    else
      dtime(master->params->output_tunit, master->timeunit, windat->t, &yr, &day, NULL);
    
    wt = (windat->hftemp) ? windat->hftemp[c] : windat->temp[c];
    sal = windat->sal[c];
    at = windat->airtemp[c];
    tcld = min(max((windat->cloud) ? windat->cloud[c] : 0.0, 0.0), 1.0);
    pres = windat->patm[c] / 100.0; /* Convert pressure to HPa */
    wspd = windat->windspeed[c];
    wdir = windat->winddir[c];

    /* Get the vapour pressures */
    dtw = get_vapour_press(windat, wincon, c, at, pres, sal, &es, &esat, &rh);

    /* Saturation vapour pressure over the water. Since the wet bulb */
    /* temp. over water is not available, ew represents the maximum */
    /* amount of water the air can accomodate ie. assume that wet */
    /* temp. = SST and rel humitity = 100% */
    ew = wt * (wt * (wt * (wt * (wt * (wt * 6.136821e-11 + 2.034081e-8) +
                                 3.03124e-6) + 2.650648e-4) +
                     1.428946e-2) + 0.4436519) + 6.1078;
    ew *= (1.0 - 0.000537 * sal); /* ew at salinity sal */

    /* Get the specific humidities (kg/kg) */
    qs = 0.622 * ew / (pres - (0.378 * ew));
    q = 0.622 * es / (pres - (0.378 * es));

    /* Get the Bowen ratio */
    bow = 0.0;
    if (ew != es)
      bow = 0.62 * (wt - at) / (ew - es);

    /* Get the long wave radiation */
    if (wincon->heatflux & COMP_LWI) {
      lg = lwrado(wt) + windat->lwri[c];
    } else
      lg = lwrad(wt, at, tcld);

    /* Get the short wave radiation */
    if (fabs(wincon->albedo) <= 1.0) {
      if (wincon->albedo >= 0.0)
	sg = windat->swr[c];
      else {
	lat = asin(wincon->coriolis[c] / (2.0 * ang));
	sg = windat->swr[c] * (1.0 - get_albedo(yr, day, tcld, lat));
      }
    } else {
      /* Get the latitude from the Coriolis parameter */
      lat = asin(wincon->coriolis[c] / (2.0 * ang));
      sg = windat->swr[c] = swrad(yr, day, at, es, tcld, lat);
    }

    /* Get the wind fetch if required */
    fetch = (wincon->fetch) ? get_fetch(wincon->fetch[c], wdir) : 0.0;

    /* Get the bulk transfer coefficients, L and HE */
    bulkp(wt, at, qs, q, wspd, pres, wincon->zref, fetch, &he, &la,
          wincon->bulkf);

    /* Get the net heat flux (Wm-2) */
    q = lg + he + la;
    windat->heatf[c] = (windat->swr_attn) ? q : sg + q;

    /* Save the heatflux to diagnostic tracers */
    windat->nhfd[c] = sg + q;
    windat->swrd[c] = sg;
    windat->lwrd[c] = lg;
    windat->shfd[c] = he;
    windat->lhfd[c] = la;

    /* Scale to ms-1K */
    windat->heatf[c] /= (-Cv * windat->dens[c]);
  }
}

/* END surf_heat_flux()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to implement surface salt flux based on evaporation from  */
/* the bulk latent heat of evaporation.                              */
/*-------------------------------------------------------------------*/
void bulk_salt_flux(geometry_t *window, /* Processing window */
		    window_t   *windat, /* Window data structure */
		    win_priv_t *wincon  /* Window geometry / constants */
  )
{
  int c, cc;                    /* Counters */
  double sal;                   /* Salinity (psu) */
  double at;                    /* Dry bulb temperature (deg C) */
  double dtw;                   /* Wet bulb temperature (deg C) */
  double q;                     /* Specific humidity at pressure es
                                   (kg/kg) */
  double es;                    /* Vapour pressure at the air temp (HPa) */
  double esat;                  /* Saturation vapour pressure (HPa) */
  double pres;                  /* Air pressure (HPa) */
  double rh;                    /* Relative humidity (%) */
  double lv;                    /* Latent heat of vaporization (J/kg) */
  double rho;                   /* Density of moist air (kg/m3) */
  double tv;                    /* Virtual temperature at ref height (deg
                                   C) */

  /* Evaporation is added to the bulk flux here based on the latent */
  /* heat of evaporation calculated in surf_heat_flux(). If this */
  /* heat flux component is not computed set to zero and return.  */
  /* Precipitation must be added.  */
  if (!windat->lhfd) {
    memset(windat->nsfd, 0, window->sgsizS);
    return;
  }

  /* Loop over the water cells */
  es = 0.0;
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->w2_t[cc];

    sal = windat->sal[c];
    at = (windat->airtemp) ? windat->airtemp[c] : 20.0;
    pres = windat->patm[c];     /* Pressure in Pa */

    /* Get the vapour pressures */
    dtw = get_vapour_press(windat, wincon, c, at, pres, sal, &es, &esat, &rh);

    /* Get the specific humidities (kg/kg) */
    q = 0.622 * es / (pres - (0.378 * es));

    /* Get the virtual temperature at the reference height (deg C) */
    tv = ((at + 273.16) * (1.0 + (0.608 * q))) - 273.16;

    /* Get the density of moist air (kg/m3) */
    rho = (3.4838e-3) * pres / (tv + 273.16);
    rho = 999.0;

    /* Get the latent heat of vaporization (J/kg) */
    lv = 4.1868 * (597.31 - 0.56525 * at) * 1e3;

    /* Calculate the evaporation in ms-1 (note lhfd is usually */
    /* negative).  */
    windat->nsfd[c] = windat->lhfd[c] / (lv * rho);
  }

  /* Add the precipitation (ms-1) */
  if (windat->precip) {
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      windat->nsfd[c] += windat->precip[c];
    }
  }
}

/* END bulk_salt_flux()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to get the vapour pressure and saturation vapour pressure */
/*-------------------------------------------------------------------*/
double get_vapour_press(window_t   *windat, 
			win_priv_t *wincon,
			int c,              /* Sparse coordinate     */
			double at,          /* Air temperature       */ 
			double pres,        /* Air pressure          */
			double sal,         /* Salinity              */
			double *es,         /* Vapour pressure       */
			double *esat,       /* Saturation pressure   */
			double *rh          /* Relative humidity     */
			)
{
  double dtw;   /* Wet bulb temperature                              */

  /* Get the vapour pressures. Options are:                          */
  /* 1. Wet bulb supplied directly from file                         */
  /* 2. Dew point supplied directly from file                        */
  /* 2. Relative humidity only supplied. Calculate dtw inversely     */
  /* using Newton's method.                                          */
  if (wincon->sh_f & WETBULB) {
    dtw = windat->wetb[c];
    /* Saturation vapour pressure in the air.                        */
    *esat =
      dtw * (dtw *
	     (dtw *
	      (dtw *
	       (dtw * (dtw * 6.136821e-11 + 2.034081e-8) + 3.03124e-6) +
	       2.650648e-4) + 1.428946e-2) + 0.4436519) + 6.1078;
    *esat *= (1.0 - 0.000537 * sal); /* esat at salinity sal */
    /* Actual vapour pressure in the air.                            */
    *es = *esat - (6.6e-4 * (1 + 1.5e-3 * dtw) * (at - dtw) * pres);
    *rh = 100.0 * (*es / *esat);
  } else if (wincon->sh_f & DEWPOINT) {
    dtw = windat->wetb[c];
    *es = dewe(dtw, pres);
    *esat = dewe(at, pres);
    *rh = 100.0 * (*es / *esat);
    dtw = estt(at, pres, *rh);
  } else if (wincon->sh_f & SPECHUM) {
    double q = windat->rh[c];
    /* Back out es from q = 0.622 * es / (pres - (0.378 * es)) */
    *es = (q * pres) / (0.622 + q * 0.378);
    *esat = dewe(at, pres);
    *rh = 100.0 * (*es / *esat);
    dtw = estt(at, pres, *rh);
  } else if (wincon->sh_f & RELHUM) {
    if (windat->rh) {
      *rh = windat->rh[c];
      dtw = estt(at, pres, *rh);
      /* Saturation vapour pressure in the air.                      */
      *esat =
	dtw * (dtw *
	       (dtw *
		(dtw *
		 (dtw * (dtw * 6.136821e-11 + 2.034081e-8) +
		  3.03124e-6) + 2.650648e-4) + 1.428946e-2) +
	       0.4436519) + 6.1078;
      *esat *= (1.0 - 0.000537 * sal); /* esat at salinity sal */
      /* Actual vapour pressure in the air.                          */
      *es = *esat - (6.6e-4 * (1 + 1.5e-3 * dtw) * (at - dtw) * pres);
    } else {
      hd_quit_and_dump
	("Insufficient information to calculate latent heat flux\n");
    }
  }
  return(dtw);
}

/* END get_vapour_pres()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to solve es/esat=ar for wet bulb temperature using        */
/* Newtons method, with ar given.                                    */
/*-------------------------------------------------------------------*/
double estt(double ta,          /* Air temperature (deg C) */
            double pr,          /* Air pressure (HPa) */
            double ar           /* Relative humidity (%) */
  )
{
  int nmax = 100;               /* Maximum number of iterations */
  double to = 2.0;              /* Tolerance = 2% */
  double tw;                    /* Wet bulb temperature (deg C) */
  double fes;                   /* Derivative of esat with respect to tw */
  double fe;                    /* Derivative of es with respect to tw */
  double esat;                  /* Saturation vapour pressure (HPa) */
  double es;                    /* Vapour pressure at the air temp (HPa) */
  int n = 0;                    /* Iterative counter */

  /* Scale */
  ar /= 100.0;
  to /= 100.0;

  /* First approximation of tw */
  tw = 0.6 * ta;

  /* Do the iterations until relh is ar% +/- to% */
  do {
    /* Get saturation vapour pressure */
    esat = tw * (tw * (tw * (tw * (tw * (tw * 6.136821e-11 + 2.034081e-8) +
                                   3.03124e-6) + 2.650648e-4) +
                       1.428946e-2) + 0.4436519) + 6.1078;
    esat *= (1.0 - 0.000537 * 35.0);  /* esat at salinity of 35.00 ppt */
    /* Get the vapur pressure */
    es = esat - (6.6e-4 * (1 + 1.5e-3 * tw) * (ta - tw) * pr);
    /* Get the derivative of esat at tw */
    fes =
      tw * (tw *
            (tw *
             (tw * (tw * 6.0 * 6.136821e-11 + 5.0 * 2.034081e-8) +
              4.0 * 3.03124e-6) + 3.0 * 2.650648e-4) + 2.0 * 1.428946e-2) +
      0.4436519;
    fes *= (1.0 - 0.000537 * 35.0); /* esat at salinity of 35.00 ppt */
    /* Get the derivative of es at tw */
    fe = fes - (6.6e-4 * pr * (1.5e-3 * ta - 1.0 - 2.0 * 1.5e-3 * tw));
    /* Find the next iteration of tw */
    tw = tw - (es - ar * esat) / (fe - ar * fes);
    n++;
  }
  while (n < nmax && (es / esat >= ar + to || es / esat <= ar - to));
  return (tw);
}

/* END estt()                                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to get the vapour pressure given the dew point temp. or   */
/* the saturation vapour pressure given the dry bulb temp.	     */
/*-------------------------------------------------------------------*/
double dewe(double tem, double pres /* Air pressure (HPa) */
  )
{
  double a1 = 373.16 / (tem + 273.16);
  double b1 = -3.49149;
  double b2 = 11.344;
  double b3 = 5.02808;
  double b4 = -7.90298;
  double b5 = 8.1328e-3;
  double b6 = -1.3816e-7;
  double a2 = a1 - 1.0;
  double a3 = 1.0 - 1.0 / a1;
  double a4 = pow(10.0, a2 * b1) - 1.0;
  double a5 = pow(10.0, a3 * b2) - 1.0;
  double e;
  e = pres * pow(a1, b3) * pow(10.0, a2 * b4 + a4 * b5 + a5 * b6);
  return (e);
}

/* END dewe()                                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the net long wave radiation at the surface.  */
/*-------------------------------------------------------------------*/
double lwrad(double sst,        /* Sea surface temperature (deg C) */
             double dtd,        /* Air temperature (deg C) */
             double tcld        /* Fractional cloud cover */
  )
{
  double em = 0.97;             /* Total hemispherical emissivity of water
                                   surf. */
  double s = 5.67e-8;           /* Stefan-Boltzmann constant (Wm-2degK-4) */
  double k = 0.92e-5;           /* Constant to evaluate emissivity of
                                   atmosphere */
  double td, d1;                /* Dummy variable */
  double lg;                    /* Long wave flux (Wm-2) */

  sst += 273.0;                 /* SST in deg K */
  td = dtd + 273.0;             /* Air temp in deg K */
  d1 = td * td;
  lg = em * s * d1 * d1 * (k * d1 - 1.0);
  /* Include the temperature jump term */
  lg = lg - 4.0 * em * s * d1 * td * (sst - td);
  lg *= (1.0 - 0.63 * tcld);
  return (lg);
}

/* END lwrad()                                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate outgoing long wave radiation at the surface. */
/*-------------------------------------------------------------------*/
double lwrado(double sst          /* Sea surface temperature (deg C) */
  )
{
  double em = 0.97;             /* Total hemispherical emissivity of water 
                                   surf. */
  double s = 5.67e-8;           /* Stefan-Boltzmann constant (Wm-2degK-4) */
  double lg;                    /* Long wave flux (Wm-2) */

  sst += 273.0;                 /* SST in deg K */
  lg = - em * s * sst * sst * sst * sst;
  return (lg);
}

/* END lwrado()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the incoming short wave radiation for clear  */
/* skies.							     */
/*-------------------------------------------------------------------*/
double swrad(int yr,            /* Year */
             double jday,       /* Day of the year */
             double tem,        /* Air temperature */
             double es,         /* Vapour pressure (HPa) */
             double mdc,        /* Fractional cloud cover */
             double lat         /* Latitude */
  )
{
  int nday;                     /* Day of the year */
  double hrs;                      /* Hour of the day */
  double se;                    /* Solar beam irradiance */
  double dec;                   /* Solar declination */
  double h;                     /* The solar elevation */
  double hrang;                 /* Hour angle */
  double pi = 3.14159;          /* Value of pi */
  double sg;                    /* Short wave radiation */
  double d1, d2;                /* Dummies */
  int iw;                       /* Counters */
  int tcld;                     /* Cloud cover (oktas) */
  double alp;                   /* Albedo of the sea surface */
  int day;                      /* Day of the month */

  nday = (int)jday;
  /*hrs = (int)(24 * (jday - (double)nday));*/
  hrs = 24 * (jday - (double)nday);
  tcld = (int)(mdc * 8);

  /*-----------------------------------------------------------------*/
  /* Get the declination and solar beam irradiance */
  d1 = nday * 2 * pi / 365.0;
  dec = 0.006918 + 0.070257 * sin(d1) - 0.399912 * cos(d1)
    + 0.000907 * sin(2 * d1) - 0.006758 * cos(2 * d1)
    + 0.00148 * sin(3 * d1) - 0.002697 * cos(3 * d1);

  /* Solar beam irradiance from Zillman (1972). Maximum irradiance 
     assumed to occur on 3rd January. */
  se = 1380.0 + 46.67 * cos((nday - 3) * pi / 180.0);

  /*-----------------------------------------------------------------*/
  /* Get the hour angle */
  hrang = (hrs - 12.0) * 180.0 / 12.0;

  /*-----------------------------------------------------------------*/
  /* Get the solar elevation */
  hrang *= pi / 180.0;
  /* dec*=pi/180.0; */
  /* lat*=pi/180.0; */
  h = sin(lat) * sin(dec) + cos(lat) * cos(dec) * cos(hrang);
  if (h < 0.0)
    h = 0.0;

  /*-----------------------------------------------------------------*/
  /* Get the incoming short wave radiation for clear skies */
  sg = se * h * h / ((h + 2.7) * 1e-3 * es + 1.085 * h + 0.1);

  /*-----------------------------------------------------------------*/
  /* Get the incoming short wave radiation for cloudy skies
     1. Parkinson & Washington (1979):
	sg *= 1.0 - 0.6 * (mdc*mdc*mdc);
     2. Gill (1982)
        sg *= 1.0 - 0.7 * mdc;
     3. Reed (1977)
       sg *= 1.0 - 0.62 * mdc + 0.0019 * asin(d2); */
  d2 = sin(lat) * sin(dec) + cos(lat) * cos(dec);
  if (tcld > 2)
       sg *= 1.0 - 0.62 * mdc + 0.0019 * asin(d2);

  /*-----------------------------------------------------------------*/
  /* Get the albedo of the sea surface */
  h = asin(h) * 180.0 / pi;
  iw = (int)(h / 10.0);
  d1 = (double)iw *10.0;
  if (iw == 0)
    alp = ap[tcld][0];
  else
    alp =
      ((ap[tcld][iw + 1] - ap[tcld][iw]) / 10.0) * (h - d1) + ap[tcld][iw];

  /*-----------------------------------------------------------------*/
  /* Get the net short wave radiation */
  sg *= (1.0 - alp);
  if (sg < 0.0 || h < 0.0)
    sg = 0.0;

  return (sg);
}

/* END swrad()                                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the albedo as a function of hour angle and   */
/* cloud amount.                                                     */
/*-------------------------------------------------------------------*/
double get_albedo(int yr,            /* Year */
		  double jday,       /* Day of the year */
		  double mdc,        /* Fractional cloud cover */
		  double lat         /* Latitude */
  )
{
  int nday;                     /* Day of the year */
  double hrs;                      /* Hour of the day */
  double dec;                   /* Solar declination */
  double h;                     /* The solar elevation */
  double hrang;                 /* Hour angle */
  double pi = 3.14159;          /* Value of pi */
  double d1, d2;                /* Dummies */
  int iw;                       /* Counters */
  int tcld;                     /* Cloud cover (oktas) */
  double alp;                   /* Albedo of the sea surface */
  int day;                      /* Day of the month */

  nday = (int)jday;
  /*hrs = (int)(24 * (jday - (double)nday));*/
  hrs = 24 * (jday - (double)nday);
  tcld = (int)(mdc * 8);

  /*-----------------------------------------------------------------*/
  /* Get the declination */
  d1 = nday * 2 * pi / 365.0;
  dec = 0.006918 + 0.070257 * sin(d1) - 0.399912 * cos(d1)
    + 0.000907 * sin(2 * d1) - 0.006758 * cos(2 * d1)
    + 0.00148 * sin(3 * d1) - 0.002697 * cos(3 * d1);

  /*-----------------------------------------------------------------*/
  /* Get the hour angle */
  hrang = (hrs - 12.0) * 180.0 / 12.0;

  /*-----------------------------------------------------------------*/
  /* Get the solar elevation */
  hrang *= pi / 180.0;
  /* dec*=pi/180.0; */
  /* lat*=pi/180.0; */
  h = sin(lat) * sin(dec) + cos(lat) * cos(dec) * cos(hrang);
  if (h < 0.0)
    h = 0.0;

  /*-----------------------------------------------------------------*/
  /* Get the albedo of the sea surface */
  h = asin(h) * 180.0 / pi;
  iw = (int)(h / 10.0);
  d1 = (double)iw *10.0;
  if (iw == 0)
    alp = ap[tcld][0];
  else
    alp =
      ((ap[tcld][iw + 1] - ap[tcld][iw]) / 10.0) * (h - d1) + ap[tcld][iw];

  return (alp);
}

/* END get_albedo()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the sensible heat and latent heat fluxes     */
/* using bulk transfer coefficients. Also calculated are the bulk    */
/* coefficients for neutral conditions and the atmospheric vertical  */
/* eddy diffusivities.                                               */
/*-------------------------------------------------------------------*/
void bulkp(double sst,          /* Sea surface temperature (deg C) */
           double dtd,          /* Air temperature at the reference height
                                   (deg C) */
           double qs,           /* Specific humidity at sea surface
                                   (kg/kg) */
           double qq,           /* Specific humidity at reference height
                                   (kg/kg) */
           double wspd,         /* Wind speed at the reference height
                                   (ms-1) */
           double pres,         /* Atmospheric pressure (HPa) */
           double zref,         /* Height to scale bulk parameters to */
           double fetch,        /* Fetch (m) */
           double *he,          /* Sensible heat flux (Wm-2) */
           double *la,          /* Latent heat flux (Wm-2) */
           int bulkf            /* Bulk scheme type */
           /* bulkf=KONDO     : Kondo (1975) */
           /* bulkf=LARGEPOND : Large and Pond (1982) */
           /* bulkf=MASAG     : Masagutov (1981) */
           /* bulkf=KITIAG    : Kitiagorodskii et al (1973) */
           /* bulkf=BUNKER    : Bunker (1976) */
  )
{
  double ta;                    /* Air temperature (deg C) */
  double pr;                    /* Air pressure (Pa) */
  double lv;                    /* Latent heat of vaporiz ation (J/kg) */
  double cp;                    /* Specific heat of moist air at const.
                                   pres.  */
  double tv;                    /* Virtual temperature at reference height
                                   (deg C) */
  double ts;                    /* Virtual temperature at the sea surface
                                   (deg C) */
  double rho;                   /* Density of moist air (kg/m3) */
  double c_h;                   /* Bulk parameter for sensible heat */
  double c_e;                   /* Bulk parameter for latent heat */

  /* Get the air temp in deg C */
  ta = dtd;

  /* Get the air pressure in Pa */
  pr = pres * 100.0;

  /* Get the latent heat of vaporization (J/kg) */
  lv = 4.1868 * (597.31 - 0.56525 * ta) * 1e3;

  /* Get the virtual temperature at the reference height (deg C) */
  tv = ((ta + 273.16) * (1.0 + (0.608 * qq))) - 273.16;

  /* Get the specific heat of moist air at constant pressure (J/kgK) */
  cp = 1.004 * (1.0 + 0.9 * qq) * 1e3;

  /* Get the density of moist air (kg/m3) */
  rho = (3.4838e-3) * pr / (tv + 273.16);

  if (bulkf == KONDO) {
    /* Kondo (1975) formulation */
    kondo(sst, ta, wspd, qs, qq, zref, &c_e, &c_h);
  } else if (bulkf == LARGEPOND) {
    /* Large and Pond (1982) formulation */
    large_pond(sst, ta, wspd, zref, &c_e, &c_h);
  } else if (bulkf == MASAG) {
    /* Masagutov (1981) formulation */
    ts = ((sst + 273.16) * (1.0 + (0.608 * qs))) - 273.16;
    c_h = c_e = masagutov(wspd, tv, ts, zref);
  } else if (bulkf == KITIAG) {
    /* Kitiagorodskii et al (1973) formulation.  */
    kitiagorodskii(wspd, zref, &c_e, &c_h);
  } else if (bulkf == BUNKER) {
    /* Bunker (1976) formulation.  */
    bunker(sst, ta, wspd, zref, &c_e, &c_h);
  }

  /* Correct for advection */
  if (fetch > 0.0) {
    double k_h = c_h * wspd * zref;
    double k_e = c_e * wspd * zref;
    ta = advt(sst, ta, wspd, fetch, k_h, zref);
    qq = advt(qs, qq, wspd, fetch, k_e, zref);
  }

  /* Get the sensible heat flux */
  *he = rho * cp * c_h * wspd * (ta - sst);

  /* Get the latent heat flux */
  *la = rho * lv * c_e * wspd * (qq - qs);

  /* Limit */
  *la = min(*la, 0.0);
}

/* END bulkp()                                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to return the bulk parameters using the scheme of         */
/* Kondo (1975).                                                     */
/*-------------------------------------------------------------------*/
double kondo(double sst,  /* Sea surface temperature (deg C)         */
	     double ta,   /* Air temperature (deg C)                 */
	     double wspd, /* Wind speed at reference height (ms-1)   */
	     double qs,   /* Specific humidity at sea surface (kg/kg)*/
	     double qq,   /* Specific humidity at 10m (kg/kg)        */
	     double zref, /* Height to scale bulk parameters to      */
	     double *c_e, /* Bulk parameter for latent heat          */
	     double *c_h  /* Bulk parameter for sensible heat        */
  )
{
  double z10 = 10.0;            /* 10m reference height */
  int jw = 0;                   /* Wind speed catagory factor */
  double k = 0.4;               /* Von Karman constant */
  double ws;                    /* Wind speed at 10m height (ms-1) */
  double sb;                    /* Stability parameter */
  double c_d;                   /* Drag bulk parameter */
  double c_d10;                 /* Drag bulk parameter at 10m height */
  double d;                     /* Dummy */

  /* Get the wind speed index assuming ws(10)=ws(ref) */
  if (wspd < 0.3)
    wspd = 0.3;
  if (wspd <= 2.2)
    jw = 0;
  if (wspd > 2.2 && wspd <= 5.0)
    jw = 1;
  if (wspd > 5.0 && wspd <= 8.0)
    jw = 2;
  if (wspd > 8.0 && wspd <= 25.0)
    jw = 3;
  if (wspd > 25.0)
    jw = 4;

  /* Get the bulk drag coefficient for neutral conditions */
  c_d10 = ad[jw] + bd[jw] * pow(wspd, pd[jw]);
  c_d10 *= 1e-3;

  /* Transform the wind speed to 10m height */
  d = exp(log(z10) - k / sqrt(c_d10));
  ws = wspd * log(z10 / d) / log(zref / d);

  /* Get the wind speed index using ws(10) */
  if (ws <= 2.2)
    jw = 0;
  if (ws > 2.2 && ws <= 5.0)
    jw = 1;
  if (ws > 5.0 && ws <= 8.0)
    jw = 2;
  if (ws > 8.0 && ws <= 25.0)
    jw = 3;
  if (ws > 25.0)
    jw = 4;

  /* Get the bulk drag coefficient for neutral conditions */
  *c_h =
    ah[jw] + bh[jw] * pow(ws, ph[jw]) + ch[jw] * (ws - 8.0) * (ws - 8.0);
  *c_e =
    ae[jw] + be[jw] * pow(ws, pe[jw]) + ce[jw] * (ws - 8.0) * (ws - 8.0);
  *c_h *= 1e-3;
  *c_e *= 1e-3;

  /* Transform to the referece height */
  d = k / sqrt(c_d10) - log(z10 / zref);
  c_d = k * k / (d * d);
  *c_h = k * sqrt(c_d) / (k * sqrt(c_d10) / (*c_h) + log(zref / z10));
  *c_e = k * sqrt(c_d) / (k * sqrt(c_d10) / (*c_e) + log(zref / z10));

  /* Sensible and latent heat flux under neutral conditions */
  /* hen=rho*cp*c_h*ws*(ta-sst); */
  /* lan=rho*lv*c_e*ws*(qq-qs); */

  if (qs == NOTVALID || qq == NOTVALID)
    return(c_d);

  /* Get the stability parameter */
  sb = (sst - ta) + 0.61 * (ta + 0.0098 * zref) * (qs - qq);
  sb =
    sb / (wspd * wspd * (1.0 + log10(z10 / zref)) *
          (1.0 + log10(z10 / zref)));
  sb = sb * (fabs(sb) / (fabs(sb) + 0.01));

  /* Stable conditions */
  if (sst < ta && sb < 0.0) {
    if (sb > -3.3 && sb < 0) {
      d = 0.1 + 0.03 * sb + 0.9 * exp(4.8 * sb);
      *c_h *= d;
      *c_e *= d;
      c_d *= d;
    } else if (sb <= -3.3) {
      *c_h = 0.0;
      *c_e = 0.0;
      c_d = 0;
    }
  }

  /* Unstable conditions */
  if (sst > ta && sb > 0.0) {
    *c_h *= 1.0 + 0.63 * sqrt(sb);
    *c_e *= 1.0 + 0.63 * sqrt(sb);
    c_d *= 1.0 + 0.47 * sqrt(sb);
  }

  /* Get the vertical eddy diffusivities (m2s-1) */
  /*
   *c_h=*c_h*ws*zref;
   *c_e=*c_e*ws*zref;
   */
  return(c_d);
}

/* END kondo()                                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to return the bulk parameters using the scheme of         */
/* Large and Pond (1982).                                            */
/*-------------------------------------------------------------------*/
double large_pond(double sst,  /* Sea surface temperature (deg C)    */
		  double ta,   /* Air temperature (deg C)            */
		  double ws,   /* Wind speed at reference height (ms-1) */
		  double zref, /* Height to scale bulk parameters to */
		  double *c_e, /* Bulk parameter for sensible heat   */
		  double *c_h  /* Bulk parameter for latent heat     */
  )
{
  double z10 = 10.0;            /* 10m reference height */
  double k = 0.4;               /* Von Karman constant */
  double ws_min = 4.0;          /* Minimum windspeed */
  double ws_max = 25.0;         /* Maximum windspeed */
  double ws_thr = 11.0;         /* Windspeed threshold */
  double c_d;                   /* Bulk drag coefficient */
  double c_d10;                 /* Drag bulk parameter at 10m height */
  double d;                     /* Dummy */

  /* Get the drag coefficient assuming ws(10)=ws(ref) */
  ws = max(ws, ws_min);
  ws = min(ws, ws_max);
  if (ws < ws_thr)
    c_d10 = 1.2;
  else
    c_d10 = 0.49 + 0.065 * ws;
  c_d10 *= 1e-3;
  
  /* Get the bulk parameters at 10m height */
  *c_e = 1.15e-3;
  *c_h = 0.66e-3;
  if (sst >= ta)
    *c_h = 1.13e-3;

  /* Transform to the referece height */
  d = k / sqrt(c_d10) - log(z10 / zref);
  c_d = k * k / (d * d);
  *c_h = k * sqrt(c_d) / (k * sqrt(c_d10) / (*c_h) + log(zref / z10));
  *c_e = k * sqrt(c_d) / (k * sqrt(c_d10) / (*c_e) + log(zref / z10));
  return(c_d);
}

/* END large_pond()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to return the bulk parameters using the scheme of         */
/* Masagutov (1981).                                                 */
/*-------------------------------------------------------------------*/
double masagutov(double wsd,    /* Wind speed at reference height (ms-1) */
                 double tv,     /* Virtual temperature at reference height
                                 */
                 double ts,     /* Virtual temperature at tthe sea surface
                                 */
                 double zref    /* Height to scale bulk parameters to */
  )
{
  double z10 = 10.0;            /* 10m reference height */
  double k = 0.4;               /* Von Karman constant */
  double ws;                    /* Wind index */
  double sb;                    /* Stability index */
  double c_b;                   /* Bulk parameter */
  double c_b10;                 /* Bulk parameter at 10m height */
  double c_d;                   /* Bulk drag coefficient */
  double c_d10;                 /* Drag bulk parameter at 10m height */
  double d;                     /* Dummy */

  if (wsd < 2.0)
    wsd = 2.0;
  sb = ts - tv + 4.0;
  if (sb < 0.0)
    sb = 0.0;
  ws = wsd - 1.0;
  if (ws < 0.0)
    ws = 0.0;
  c_b10 = interp_maov(ws, sb, 1) * 1e-3;
  c_d10 = interp_maov(ws, sb, 0) * 1e-3;

  /* Transform to the referece height */
  d = k / sqrt(c_d10) - log(z10 / zref);
  c_d = k * k / (d * d);
  c_b = k * sqrt(c_d) / (k * sqrt(c_d10) / c_b10 + log(zref / z10));

  return (c_b);
}

/* END masagutov()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to return the bulk parameters using the scheme of         */
/* Kitiagorodskii et al (1973).                                      */
/*-------------------------------------------------------------------*/
double kitiagorodskii(double ws,   /* Wind speed at reference height (ms-1) */
		      double zref, /* Height to scale bulk parameters to */
		      double *c_e, /* Bulk parameter for sensible heat */
		      double *c_h  /* Bulk parameter for latent heat */
  )
{
  double z10 = 10.0;            /* 10m reference height */
  double k = 0.4;               /* Von Karman constant */
  double ws_max = 22.0;         /* Maximum windspeed */
  double c_d;                   /* Bulk drag coefficient */
  double c_d10;                 /* Drag bulk parameter at 10m height */
  double d, r;                  /* Dummy */

  /* Get the drag coefficient assuming ws(10)=ws(ref) */
  ws = min(ws, ws_max);
  if (ws < 1.0)
    ws = 1.0;
  r = pow(10.0, -0.616 + 0.137 * ws);
  c_d10 = 1.2e-3 * pow(r, 0.15);

  /* Get the bulk parameters at 10m height */
  *c_e = 1.0e-3 * pow(r, 0.11);
  *c_h = *c_e;

  /* Transform to the referece height */
  d = k / sqrt(c_d10) - log(z10 / zref);
  c_d = k * k / (d * d);
  *c_h = k * sqrt(c_d) / (k * sqrt(c_d10) / (*c_h) + log(zref / z10));
  *c_e = k * sqrt(c_d) / (k * sqrt(c_d10) / (*c_e) + log(zref / z10));
  return(c_d);
}

/* END kitiagorodskii()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to return the bulk parameters using the scheme of         */
/* Bunker (1976).                                                    */
/*-------------------------------------------------------------------*/
double bunker(double sst,  /* Sea surface temperature (deg C)        */
	      double ta,   /* Air temperature (deg C)                */
	      double ws,   /* Wind speed at reference height (ms-1)  */
	      double zref, /* Height to scale bulk parameters to     */
	      double *c_e, /* Bulk parameter for sensible heat       */
	      double *c_h  /* Bulk parameter for latent heat         */
  )
{
  double z10 = 10.0;            /* 10m reference height */
  double k = 0.4;               /* Von Karman constant */
  double c_d;                   /* Bulk drag coefficient */
  double c_d10;                 /* Drag bulk parameter at 10m height */
  double d;                     /* Dummy */
  int w, t;

  /* Get the drag coefficient assuming ws(10)=ws(ref) */
  /* Get the stability (temp difference) index */
  d = fabs(ta - sst);
  if (d >= 5.0)
    t = 0;
  else if(d < 5 && d >= 1.0)
    t = 1;
  else if(d < 1 && d >= 0.2)
    t = 2;
  else if(d < 0.2 && d >= -0.2)
    t = 3;
  else if(d < -0.2 && d >= -1.0)
    t = 4;
  else if(d < -1.0 && d >= -5.0)
    t = 5;
  else
    t = 6;

  /* Get the wind index */
  w = min((int)(ws / 5.0), 8);

  /* Get the wind index */
  if (ws <= 15.0)
    w = (int)(ws / 3.0);
  else
    w = min((int)(ws / 5.0) + 2, 8);
  c_d10 = bunker_Cd[w][t] * 1e-3;

  /* Get the bulk parameters at 10m height */
  /* Get the wind index */
  if (ws <= 15.0)
    w = (int)(ws / 3.0);
  else
    w = min((int)(ws / 5.0) + 2, 8);
  *c_e = *c_h = bunker_Ce[w][t] * 1e-3;

  /* Transform to the referece height */
  d = k / sqrt(c_d10) - log(z10 / zref);
  c_d = k * k / (d * d);
  *c_h = k * sqrt(c_d) / (k * sqrt(c_d10) / (*c_h) + log(zref / z10));
  *c_e = k * sqrt(c_d) / (k * sqrt(c_d10) / (*c_e) + log(zref / z10));
  return (c_d);
}

/* END bunker()                                                      */
/*-------------------------------------------------------------------*/


/*--------------------------------------------------------------------*/
/* Performs bilinear interpolation at point (a,b,c)                   */
/*--------------------------------------------------------------------*/
double interp_maov(double a, double b,  /* Windspeed and stability */
                   int mode     /* mode=1 : c_e,c_h, mode=0 : c_d */
  )
{
  int ix;                       /* Counter */
  int l, m;                     /* Integer parts of a,b */
  double vals[4];               /* Mesh points surrounding A[b][a] */
  double wgt[4];                /* Weights for bilinear interpolation */
  double nval;                  /* Interpolated value */
  double p, q;                  /* Decimal parts of a,b */
  int l_nce1 = 21;              /* Limit of grid in the x direction */
  int l_nce2 = 8;               /* Limit of grid in the y direction */
  int lp1;
  int mp1;

  /* Get the fractional parts of a,b and c */
  p = a - floor(a);
  q = b - floor(b);

  /* Get the weights */
  wgt[0] = (1.0 - p) * (1.0 - q);
  wgt[1] = (1.0 - p) * q;
  wgt[2] = p * (1.0 - q);
  wgt[3] = p * q;

  /* Get the integer parts of a and b and surrounding points.  */
  /* If the particle passes outside the domain, then set the position */
  /* equal to the relevant boundary.  */
  /* If the surrounding points are outside the grid, set the relevant */
  /* point to the nearest boundary point.  */
  if (a < 0.0)
    l = lp1 = 0;
  else if (a >= l_nce1)
    l = lp1 = l_nce1;
  else {
    l = (int)a;
    lp1 = l + 1;
  }
  if (b < 0.0)
    m = mp1 = 0;
  else if (b >= l_nce2)
    m = mp1 = l_nce2;
  else {
    m = (int)b;
    mp1 = b + 1;
  }

  /* Get the eight integral neighbours of our non-integral point */
  if (mode) {
    vals[0] = maov[m][l];
    vals[1] = maov[mp1][l];
    vals[2] = maov[m][lp1];
    vals[3] = maov[mp1][lp1];
  } else {
    vals[0] = maod[m][l];
    vals[1] = maod[mp1][l];
    vals[2] = maod[m][lp1];
    vals[3] = maod[mp1][lp1];
  }

  /* Use bilinear interpolation */
  nval = 0.0;
  for (ix = 0; ix < 4; ix++)
    nval += wgt[ix] * vals[ix];

  return (nval);
}

/* END interp_maov()                                                  */
/*--------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to get day no. for a given year, month & day		     */
/*-------------------------------------------------------------------*/
int yrday(int year, int mon, int day)
{
  int i, lp;
  lp = year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
  for (i = 1; i < mon; i++)
    day += dayno[lp][i];
  return (day);
}

/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to get day of the month for a given day number            */
/*-------------------------------------------------------------------*/
int mday(int nday,              /* Julian day */
         int year               /* Year */
  )
{
  int lp;                       /* lp = 1 for leap years, lp = 0 otherwise
                                 */
  int i;                        /* Counter */
  int day;                      /* Last day of the month counter */
  int oday;                     /* Previous last day of the month */

  if (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0))
    lp = 1;
  else
    lp = 0;

  day = oday = dayno[lp][0];
  for (i = 1; i <= 12; i++) {
    if (day >= nday)
      break;
    oday = day;
    day = day + dayno[lp][i];
  }
  return (nday - oday);
}

/* END mday()                                                        */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Returns a year given a fractional Julian day of the form          */
/* YYYYDDD.frac.                                                     */
/*-------------------------------------------------------------------*/
double jyear(double day)
{
  return ((double)floor(day / 1e3));
}

/* END jyear()                                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Returns a Julian day given a fractional Julian day of the form    */
/* YYYYDDD.frac.                                                     */
/*-------------------------------------------------------------------*/
double jday(double day)
{
  if (day / 1e3 < 1.0)
    return (day);
  else
    return (day - jyear(day) * 1e3);
}

/* END jday()                                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Returns the difference between two fractional Julian days.        */
/*-------------------------------------------------------------------*/
double julian_dif(double j1, double j2  /* Fractional Julian days */
  )
{

  if (j1 / 1e3 >= 1.0 && j2 / 1e3 >= 1.0) {
    /* Time format of the form YYYYDDD.frac */
    return (fabs((jyear(j1) - jyear(j2)) * 365.0 + jday(j1) - jday(j2)));
  } else {
    /* Time format of the form DDD.frac */
    return (fabs(j1 - j2));
  }
}

/* END julian_dif()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Returns the sum of two fractional Julian days, j1+days.           */
/*-------------------------------------------------------------------*/
double julian_sum(double j1,    /* Fractional Julian day.  */
                  double days   /* Number of days to add.  */
  )
{
  int year;
  double yearlen;
  double jsum;

  if (j1 / 1e3 >= 1.0) {
    /* Time format of the form YYYYDDD.frac.  */
    /* Get the year and total Julian days.  */
    year = (int)jyear(j1);
    jsum = jday(j1) + days;

    /* Get the length of the year in days.  */
    if (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0))
      yearlen = 366.0;
    else
      yearlen = 365.0;

    /* Find the new year and residual number of days.  */
    while (jsum > yearlen) {
      year = year + 1;
      jsum = jsum - yearlen;
      if (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0))
        yearlen = 366.0;
      else
        yearlen = 365.0;
    }

    /* Get the new fractional julian day in YYYYDDD.frac.  */
    return ((double)year * 1e3 + jsum);
  } else {
    /* Time format of the form DDD.frac */
    return (j1 + days);
  }
}

/* END julian_sum()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the value of the error function erf(x)	     */
/*-------------------------------------------------------------------*/
double erfn(double x)
{
  double a1 = 0.254829592;
  double a2 = -0.284496736;
  double a3 = 1.421413741;
  double a4 = -1.453152027;
  double a5 = 1.061405429;
  double p = 0.3275911;
  double t;

  t = 1.0 / (1.0 + p * x);
  x = 1 - t * (a1 + t * (a2 + t * (a3 + t * (a4 + t * a5)))) * exp(-x * x);
  return (x);
}

/* END erfn()                                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to find the air temperaure at 2m height with wind u       */
/* blowing in the offshore direction.				     */
/*-------------------------------------------------------------------*/
double advt(double ts,          /* Sea temperature (boundary condition at
                                   z=0, x>0 */
            double ta,          /* Initial air temperature profile x=0,
                                   z>0 */
            double u,           /* Wind speed */
            double x,           /* Distance offshore */
            double K,           /* Eddy conductivity */
            double zref         /* Height to calculate air temperature */
  ) {

  double t;                     /* Air temperature at height z, distance x
                                 */
  if (K == 0.0 || x == 0.0)
    t = ta;
  else
    t = ts + (ta - ts) * erfn((zref * sqrt(u)) / sqrt(4.0 * K * x));
  return (t);
}

/* END advt()                                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the heatflux required to raise the temp      */
/* to that observed.                                                 */
/*-------------------------------------------------------------------*/
void calc_heatflux(geometry_t *window,
		   window_t   *windat, 
		   win_priv_t *wincon
		   )
{
  double Cv = 4e3;              /* Specific heat at constant volume */
  int c, cc;                    /* Counters */
  double obs;                   /* Observed sea surface temperature */

  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->w2_t[cc];
    obs = windat->hftemp[c];
    windat->heatf[c] = windat->nhfd[c] = Cv * windat->mixl[c] *
      windat->dens[c] * (windat->temp[c] - obs) / wincon->hftc;
    /* Scale to ms-1K */
    windat->heatf[c] /= (Cv * windat->dens[c]);
  }
}

/* END calc_heatflux()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to distribute the short wave radiation with depth         */
/*-------------------------------------------------------------------*/
void set_sbc(geometry_t *window,  /* Processing window               */
             window_t *windat,    /* Window data structure           */
             win_priv_t *wincon   /* Window geometry / constants     */
  )
{
  double i_top;                 /* Short wave at layer top           */
  double i_bot;                 /* Short wave at layer bottom        */
  double ctop;                  /* Depth of the top of a layer       */
  double dz;                    /* Layer thickness                   */
  int c, cc, cs, cb, c3;        /* Counters                          */
  double Cv = 4e3;              /* Specific heat at constant volume  */
  double attn;                  /* SWR attenuation                   */

  /*-----------------------------------------------------------------*/
  /* If the heat flux is calculated and an attenuation coefficient */
  /* is supplied the distribute the short wave radiation throughout */
  /* the water column.  */
  if (wincon->heatflux & (ADVANCED|NET_HEAT|COMP_HEAT|COMP_HEAT_MOM)
      && windat->swr_attn) {
    if (!(windat->swr_attn1)) {
      for (cc = 1; cc <= wincon->vcs1; cc++) {
	c = wincon->s1[cc];
	cs = window->m2d[c];  /* 2D coordinate corresponding to c */
	cb = wincon->i1[cc];  /* Bottom coordinate for cells to process */

	ctop = windat->eta[cs];

	/* No heatflux for columns one cell deep */
	if ((c == cb) || (windat->eta[cs] - window->gridz[c] <= wincon->hmin)) {
	  windat->heatf[cs] = 0.0;
	  continue;
	}
	windat->heatf[cs] -= ((1.0 - windat->swr_tran[cs]) * windat->swr[cs] /
			      (Cv * windat->dens[c]));
	i_top = windat->swr_tran[cs] * windat->swr[cs];

	c3 = c;
	while (c3 != cb) {
	  dz = (ctop - window->gridz[c3]) * wincon->Ds[cs];
	  attn = (wincon->swr_type & SWR_2D) ? windat->swr_attn[cs] : windat->swr_attn[c3];
	  i_bot = decay_exact(i_top, attn, dz);

	  windat->temp[c3] += windat->dt * (i_top - i_bot) /
	    (Cv * dz * windat->dens[c3]);
	  i_top = i_bot;
	  ctop = window->gridz[c3];
	  c3 = window->zm1[c3];
	}
	dz = (ctop - window->botz[cs]) * wincon->Ds[cs];
	attn = (wincon->swr_type & SWR_2D) ? windat->swr_attn[cs] : windat->swr_attn[cb];
	i_bot = decay_exact(i_top, attn, dz);
	windat->temp[cb] += windat->dt * (i_top - windat->swr_babs[cs] * i_bot) /
	  (Cv * dz * windat->dens[cb]);
	/*
	windat->temp[cb] += windat->dt * i_top /
	  (Cv * (ctop - window->botz[cs]) * windat->dens[cb]);
	*/
      }
    } else {
      for (cc = 1; cc <= wincon->vcs1; cc++) {
	c = wincon->s1[cc];
	cs = window->m2d[c]; /* 2D coordinate corresponding to c */
	cb = wincon->i1[cc]; /* Bottom coordinate for cells to process */

	ctop = windat->topz[cs];

	/* No heatflux for columns one cell deep */
	if ((c == cb) || (windat->eta[cs] - 
			  window->gridz[c] * wincon->Ds[cs] <= wincon->hmin)) {
	  windat->heatf[cs] = 0.0;
	  continue;
	}

	i_top = windat->swr[cs];

	c3 = c;
	while (c3 != cb) {
	  dz = (ctop - window->gridz[c3]) * wincon->Ds[cs];
	  i_bot = i_top *
	    (decay_exact(windat->swr_tran[cs], windat->swr_attn[cs], dz) +
	     decay_exact((1.0 - windat->swr_tran[cs]), windat->swr_attn1[cs], dz));

	  windat->temp[c3] += windat->dt * (i_top - i_bot) /
	    (Cv * dz * windat->dens[c3]);
	  i_top = i_bot;
	  ctop = window->gridz[c3];
	  c3 = window->zm1[c3];
	}
	dz = (ctop - window->botz[cs]) * wincon->Ds[cs];
	i_bot = i_top * 
	  (decay_exact(windat->swr_tran[cs], windat->swr_attn[cs], dz) +
	   decay_exact((1.0 - windat->swr_tran[cs]), windat->swr_attn1[cs], dz));
	windat->temp[cb] += windat->dt * (i_top - windat->swr_babs[cs] * i_bot) /
	  (Cv * dz * windat->dens[cb]);
	/*
	windat->temp[cb] += windat->dt * i_top /
	  (Cv * (ctop - window->botz[cs]) * wincon->Ds[cs] * windat->dens[cb]);
	*/
      }
    }
  }
}

/* END set_sbc()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to distribute the short wave radiation with depth and     */
/* return the temperature source for vertical diffusion.             */
/*-------------------------------------------------------------------*/
void calc_swr(geometry_t *window,  /* Processing window              */
	      window_t *windat,    /* Window data structure          */
	      win_priv_t *wincon,  /* Window geometry / constants    */
	      double *So,          /* Temperature source term        */
	      int cci              /* Sparse index                   */
  )
{
  double i_top;                 /* Short wave at layer top           */
  double i_bot;                 /* Short wave at layer bottom        */
  double ctop;                  /* Depth of the top of a layer       */
  double dz;                    /* Layer thickness                   */
  int c, cc, cs, cb, c3;        /* Counters                          */
  double Cv = 4e3;              /* Specific heat at constant volume  */
  double attn;                  /* SWR attenuation                   */
  int ccs = 1;
  int cce = wincon->vcs;

  if (cci) {
    ccs = cci;
    cce = cci;
  }

  /*-----------------------------------------------------------------*/
  /* If the heat flux is calculated and an attenuation coefficient */
  /* is supplied the distribute the short wave radiation throughout */
  /* the water column.  */
  if (wincon->heatflux & (ADVANCED|NET_HEAT|COMP_HEAT|COMP_HEAT_MOM)
      && windat->swr_attn) {
    if (!(windat->swr_attn1)) {
      for (cc = ccs; cc <= cce; cc++) {
	c = wincon->s1[cc];
	cs = window->m2d[c];  /* 2D coordinate corresponding to c */
	cb = wincon->i1[cc];  /* Bottom coordinate for cells to process */

	ctop = windat->eta[cs];

	/* No heatflux for columns one cell deep */
	if ((c == cb) || (windat->eta[cs] - window->gridz[c] <= wincon->hmin)) {
	  windat->heatf[cs] = 0.0;
	  continue;
	}
	windat->heatf[cs] -= ((1.0 - windat->swr_tran[cs]) * windat->swr[cs] /
			      (Cv * windat->dens[c]));
	i_top = windat->swr_tran[cs] * windat->swr[cs];

	c3 = c;
	while (c3 != cb) {
	  dz = (ctop - window->gridz[c3]) * wincon->Ds[cs];
	  attn = (wincon->swr_type & SWR_2D) ? windat->swr_attn[cs] : windat->swr_attn[c3];
	  i_bot = decay_exact(i_top, attn, dz);

	  So[c3] = (i_top - i_bot) / (Cv * dz * windat->dens[c3]);
	  i_top = i_bot;
	  ctop = window->gridz[c3];
	  c3 = window->zm1[c3];
	}
	dz = (ctop - window->botz[cs]) * wincon->Ds[cs];
	attn = (wincon->swr_type & SWR_2D) ? windat->swr_attn[cs] : windat->swr_attn[cb];
	i_bot = decay_exact(i_top, attn, dz);
	So[cb] = (i_top - windat->swr_babs[cs] * i_bot) /
	  (Cv * dz * windat->dens[cb]);
      }
    } else {
      for (cc = 1; cc <= cce; cc++) {
	c = wincon->s1[cc];
	cs = window->m2d[c]; /* 2D coordinate corresponding to c */
	cb = wincon->i1[cc]; /* Bottom coordinate for cells to process */

	ctop = windat->topz[cs];

	/* No heatflux for columns one cell deep */
	if ((c == cb) || (windat->eta[cs] - 
			  window->gridz[c] * wincon->Ds[cs] <= wincon->hmin)) {
	  windat->heatf[cs] = 0.0;
	  continue;
	}
	i_top = windat->swr[cs];

	c3 = c;
	while (c3 != cb) {
	  dz = (ctop - window->gridz[c3]) * wincon->Ds[cs];
	  i_bot = i_top *
	    (decay_exact(windat->swr_tran[cs], windat->swr_attn[cs], dz) +
	     decay_exact((1.0 - windat->swr_tran[cs]), windat->swr_attn1[cs], dz));
	  So[c3] = (i_top - i_bot) / (Cv * dz * windat->dens[c3]);
	    
	  i_top = i_bot;
	  ctop = window->gridz[c3];
	  c3 = window->zm1[c3];
	}
	dz = (ctop - window->botz[cs]) * wincon->Ds[cs];
	i_bot = i_top * 
	  (decay_exact(windat->swr_tran[cs], windat->swr_attn[cs], dz) +
	   decay_exact((1.0 - windat->swr_tran[cs]), windat->swr_attn1[cs], dz));
        So[cb] = (i_top - windat->swr_babs[cs] * i_bot) /
	  (Cv * dz * windat->dens[cb]);
      }
    }
  }
}


/* END calc_swr()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to relax the sea surface to observations                  */
/*-------------------------------------------------------------------*/
void surf_relax(geometry_t *window, /* Processing window */
                window_t *windat, /* Window data structure */
                win_priv_t *wincon  /* Window geometry / constants */
  )
{
  int c, cc, cs;            /* Counters */
  double *sst;

  if (wincon->heatflux & SURF_RELAX)
    sst = windat->hftemp;
  else if(wincon->heatflux & AVHRR)
    sst = windat->avhrr;
  else if(wincon->heatflux & GHRSST)
    sst = windat->ghrsst;
  else
    return;

  for (cc = 1; cc <= wincon->vcs1; cc++) {
    c = wincon->s1[cc];
    cs = window->m2d[c];        /* 2D coordinate corresponding to c */

    windat->temp[c] = sst[cs] +
      decay_exact(windat->temp[c] - sst[cs], 1.0 / wincon->hftc, windat->dt);

  }
}

/* END surf_relax()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the surface salt flux term                   */
/*-------------------------------------------------------------------*/
void surf_salt_flux(geometry_t *window, /* Processing window */
		    window_t   *windat, /* Window data structure */
		    win_priv_t *wincon  /* Window geometry / constants */
  )
{
  int c, cc; /* Counters */

  memset(windat->nsfd, 0, window->sgsizS * sizeof(double));

  if (windat->evap) {
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      windat->nsfd[c] -= windat->evap[c];
    }
  }
  if (windat->tr_wcS[windat->evapn]) {
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      windat->tr_wcS[windat->evapn][c] = windat->evap[c];
    }
  }
  if (windat->precip) {
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      windat->nsfd[c] += windat->precip[c];
      windat->tr_wcS[windat->precipn][c] = windat->precip[c];
    }
  }
  if (windat->tr_wcS[windat->precipn]) {
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      windat->tr_wcS[windat->precipn][c] = windat->precip[c];
    }
  }
}

/* END surf_salt_flux()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to get the swr constructed from a sinusoidal profile      */
/* given a daily mean. Uses Eq. 2 Schiller and Godfrey (2003),       */
/* Journal of Climate, 16, 21-39.                                    */
/*-------------------------------------------------------------------*/
double swr_from_mean(double swr, double lat, double t)
{
  double es = 26.0;     /* Default vapour pressure                     */
  double at = 20.0;     /* Default air temperature                     */
  double cld = 0.0;     /* Default cloud amount                        */
  double inc = 600.0;   /* Precision (s)                               */
  double dawn, dusk;    /* Times of dawn and dusk to nearest inc       */
  double st;            /* Start of current day                        */
  int yr;               /* Year                                        */
  double day;           /* Julian day                                  */
  double s, ps;         /* swr and previous swr                        */
  double n;             /* Counter                                     */

  st = (floor(t / 86400.0)) * 86400.0;
  /* Get time of dawn to the nearest inc                             */
  ps = 0.0;
  for (n = st; n < st + 86400.0; n += inc) {
    dtime(master->params->output_tunit, master->timeunit, n, &yr, &day, NULL);
    s = swrad(yr, day, at, es, cld, lat);    
    if (ps <= 0.0 && s > 0) {
      dawn = n;
      break;
    }
    ps = s;
  }
  /* Get time of dusk to the nearest inc                             */
  ps = 0.0;
  for (n = st + 86400.0; n > st; n -= inc) {
    dtime(master->params->output_tunit, master->timeunit, n, &yr, &day, NULL);
    s = swrad(yr, day, at, es, cld, lat);    
    if (ps <= 0.0 && s > 0) {
      dusk = n;
      break;
    }
    ps = s;
  }
  if (t >= dawn && t <= dusk)
    return(max(0.0, swr * PI * sin(2 * PI * (t - dawn) / 86400.0)));
  else
    return(0.0);
}

/* END swr_from_mean()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Returns 1 if it's night, 0 if it's day                            */
/*-------------------------------------------------------------------*/
double is_night(double lat, double t)
{
  return (swr_from_mean(600.0, lat, t)) ? 0 : 1;
}
/* END is_night()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
double get_cloud(master_t *master,  /* Master data structure */
		 double rh,         /* Relative humidity     */
		 int c
		 )
{

  double wt = master->temp[c];
  double sal = master->sal[c];
  double at = (master->airtemp) ? master->airtemp[c] : 20.0;
  double pres = master->patm[c] / 100.0; /* Convert pressure to HPa */
  double wspd = master->windspeed[c];
  double wdir = master->winddir[c];
  double ang = 7.29e-5;         /* Earth's angular velocity (s-1) */
  double lat = asin(master->coriolis[c] / (2.0 * ang));
  double dtw;
  double esat, es;

  dtw = estt(at, pres, rh);
  /* Saturation vapour pressure in the air.  */
  esat =
    dtw * (dtw *
	   (dtw *
	    (dtw *
	     (dtw * (dtw * 6.136821e-11 + 2.034081e-8) +
	      3.03124e-6) + 2.650648e-4) + 1.428946e-2) +
	   0.4436519) + 6.1078;
  esat *= (1.0 - 0.000537 * sal); /* esat at salinity sal */
  /* Actual vapour pressure in the air.  */
  es = esat - (6.6e-4 * (1 + 1.5e-3 * dtw) * (at - dtw) * pres);

  /*sg = swrad(yr, day, at, es, 0.0, lat);*/
  return 0.0;
}

/* END get_cloud()                                                   */
/*-------------------------------------------------------------------*/
