/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/control/dp_none.c
 *  
 *  Description:
 *  Manage distributed processing of windows in
 *  a series. That's don't distribute the processing
 *  acrossM multiple CPUs.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: dp_none.c 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <time.h>
#include "hd.h"

void mode2d_step_window_p1(master_t *master,
			   geometry_t *window,
			   window_t *windat, win_priv_t *wincon);
void mode2d_step_window_p2(master_t *master,
			   geometry_t *window,
			   window_t *windat, win_priv_t *wincon);
void mode2d_step_window_p3(master_t *master,
			   geometry_t *window,
			   window_t *windat, win_priv_t *wincon);
void tracer_step_window(master_t *master,
                        geometry_t *window,
                        window_t *windat, win_priv_t *wincon);
void mode3d_step_window_p1(master_t *master,
                           geometry_t *window,
                           window_t *windat, win_priv_t *wincon);
void mode3d_step_window_p2(master_t *master,
                           geometry_t *window,
                           window_t *windat, win_priv_t *wincon);
void mode3d_post_window_p1(master_t *master,
                           geometry_t *window,
                           window_t *windat, win_priv_t *wincon);
void mode3d_post_window_p2(master_t *master,
                           geometry_t *window,
                           window_t *windat, win_priv_t *wincon);

/* Perform the DP in serial mode. */
void dp_none_init(dp_window_t *dpw)
{
}

void dp_none_cleanup(dp_window_t *dpw)
{
}

void dp_none_vel2d_step_p1(dp_window_t *dpw)
{
  TIMING_SET;
  mode2d_step_window_p1(dpw->master, dpw->geom, dpw->windata, dpw->wincon);
  TIMING_DUMP_WIN(3, "   vel2d_s1_win",dpw->window_id);
}

void dp_none_vel2d_gather_step_p1(dp_window_t *dpw)
{
}

void dp_none_vel2d_step_p2(dp_window_t *dpw)
{
  TIMING_SET;
  mode2d_step_window_p2(dpw->master, dpw->geom, dpw->windata, dpw->wincon);
  TIMING_DUMP_WIN(3, "   vel2d_s2_win",dpw->window_id);
}

void dp_none_vel2d_gather_step_p2(dp_window_t *dpw)
{
}

void dp_none_vel2d_step_p3(dp_window_t *dpw)
{
  TIMING_SET;
  mode2d_step_window_p3(dpw->master, dpw->geom, dpw->windata, dpw->wincon);
  TIMING_DUMP_WIN(3, "   vel2d_s3_win",dpw->window_id);
}

void dp_none_vel2d_gather_step_p3(dp_window_t *dpw)
{
}

void dp_none_vel3d_step_p1(dp_window_t *dpw)
{
  TIMING_SET;
  mode3d_step_window_p1(dpw->master, dpw->geom, dpw->windata, dpw->wincon);
  TIMING_DUMP_WIN(3, "   vel3d_s1_win",dpw->window_id);
}

void dp_none_vel3d_gather_step_p1(dp_window_t *dpw)
{
}

void dp_none_vel3d_step_p2(dp_window_t *dpw)
{
  TIMING_SET;
  mode3d_step_window_p2(dpw->master, dpw->geom, dpw->windata, dpw->wincon);
  TIMING_DUMP_WIN(3, "   vel3d_s2_win", dpw->window_id);
}

void dp_none_vel3d_gather_step_p2(dp_window_t *dpw)
{
}

void dp_none_vel3d_post_p1(dp_window_t *dpw)
{
  mode3d_post_window_p1(dpw->master, dpw->geom, dpw->windata, dpw->wincon);
}

void dp_none_vel3d_gather_post_p1(dp_window_t *dpw)
{
}

void dp_none_vel3d_post_p2(dp_window_t *dpw)
{
  mode3d_post_window_p2(dpw->master, dpw->geom, dpw->windata, dpw->wincon);
}

void dp_none_vel3d_gather_post_p2(dp_window_t *dpw)
{
}


void dp_none_tracer_step(dp_window_t *dpw)
{
  tracer_step_window(dpw->master, dpw->geom, dpw->windata, dpw->wincon);
}

void dp_none_tracer_gather_step(dp_window_t *dpw)
{
}
