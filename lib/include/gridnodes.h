/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  File: model/lib/grid/include/gridnodes.h
 *
 *  Description: Header file for handling grid node arrays
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: gridnodes.h 6595 2020-09-03 03:36:52Z riz008 $
 *
 */

#if !defined(_GRIDNODES_H)
#define GRIDNODES_H

typedef enum {
    NT_NONE = 0,                /* not specified */
    NT_DD = 1,                  /* double density */
    NT_CEN = 2,                 /* cell centers */
    NT_COR = 3                  /* cell corners */
} NODETYPE;

typedef enum {
    CT_X = 0,
    CT_Y = 1,
    CT_XY = 2
} COORDTYPE;

struct gridnodes;
typedef struct gridnodes gridnodes;

extern char* nodetype2str[];

gridnodes* gridnodes_read(char* fname, NODETYPE type);
gridnodes* gridnodes_create(int nx, int ny, NODETYPE type);
void gridnodes_readnextpoint(gridnodes* gn, double x, double y);
void gridnodes_destroy(gridnodes* gn);

void gridnodes_applymask(gridnodes* gn, int** mask);
void gridnodes_calcstats(gridnodes* gn);
gridnodes* gridnodes_copy(gridnodes* gn);
gridnodes* gridnodes_subgrid(gridnodes* gn, int imin, int imax, int jmin, int jmax);
gridnodes* gridnodes_transform(gridnodes* gn, NODETYPE newtype);
void gridnodes_validate(gridnodes* gn);
void gridnodes_write(gridnodes* gn, char* fname, COORDTYPE ctype);

int gridnodes_getnx(gridnodes* gn);
int gridnodes_getny(gridnodes* gn);
double** gridnodes_getx(gridnodes* gn);
double** gridnodes_gety(gridnodes* gn);
int gridnodes_getnce1(gridnodes* gn);
int gridnodes_getnce2(gridnodes* gn);

#endif
