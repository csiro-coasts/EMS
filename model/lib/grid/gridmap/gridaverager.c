/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  File: model/lib/grid/gridmap/gridaverager.c
 *  
 *  Description:
 *  Calculate average value of a 2D field a cell of a grid
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: gridaverager.c 5839 2018-06-28 03:38:07Z riz008 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include "nan.h"
#include "gridmap.h"
#include "gridaverager.h"
#include "grid_utils.h"
#include "ems.h"

struct gridaverager {
    gridmap* gm;
    double** v;
    int** n;
};

gridaverager* ga_create(gridmap* gm)
{
    gridaverager* ga = malloc(sizeof(gridaverager));

    ga->gm = gm;
    ga->v = d_alloc_2d(gridmap_getnce1(gm), gridmap_getnce2(gm));
    ga->n = i_alloc_2d(gridmap_getnce1(gm), gridmap_getnce2(gm));

    return ga;
}

void ga_destroy(gridaverager* ga)
{
    d_free_2d(ga->v);
    i_free_2d(ga->n);
    free(ga);
}

void ga_addpoints(gridaverager* ga, int n, point points[])
{
    gridmap* gm = ga->gm;
    int ii;

    for (ii = 0; ii < n; ++ii) {
        point* p = &points[ii];
        int i, j;

        if (gridmap_xy2ij(gm, p->x, p->y, &i, &j) != 0) {
            ga->v[j][i] += p->z;
            ga->n[j][i]++;
        }
    }
}

void ga_getvalue(gridaverager* ga, point * p)
{
    int i, j, n;

    if (gridmap_xy2ij(ga->gm, p->x, p->y, &i, &j) != 0 && (n = ga->n[j][i]) != 0)
        p->z = ga->v[j][i] / (double) n;
    else
        p->z = NaN;
}
