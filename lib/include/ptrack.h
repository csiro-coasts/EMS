/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/ptrack.h
 *
 *  \brief Include file for particle tracking routines
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: ptrack.h 7419 2023-10-05 02:14:39Z her127 $
 */


#ifndef _PTRACK_H
#define _PTRACK_H 1

typedef struct {
  double e1;
  double e2;
  double e3;
  int c;
  short flag;
  int dumpf;
  double age;
  unsigned char out_age;
  double size;
  unsigned char out_size;
  double svel;
} particle_t;

#define PT_ACTIVE 0x001
#define PT_LOST   0x002
#define PT_AGE    0x004
#define PT_SIZE   0x008
#define PT_IN     0x020
#define PT_OUT    0x040
#define PT_WIND   0x080
#define PT_FATT   0x100

int pt_create(char *name, long np, char *t_units, int dumpf);
void pt_read(char *name, int rec, long *np, particle_t **p, double *t,
             char *t_units, int *ndump);
void pt_write(int fid, int rec, double t, long np, particle_t *p);
void pt_write_a(int fid, int rec, double t, long np, particle_t *p);

#endif                          /* _PTRACK_H */
