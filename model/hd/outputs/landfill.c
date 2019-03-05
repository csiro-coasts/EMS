/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/outputs/landfill.c
 *  
 *  Description:
 *  Extrapolate the tracer values onto land.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: landfill.c 5988 2018-10-17 04:43:39Z her127 $
 *
 */

#include <math.h>
#include <string.h>
#include "hd.h"

int is_layer_solid(dump_data_t *dumpdata, int k);
int find_closest(dump_data_t *dumpdata, int oi, int oj, int k, int *ci,
                 int *cj);
void def_fill_land(dump_data_t *dumpdata);
void zero_fill_land(dump_data_t *dumpdata);
void no_fill_land(dump_data_t *dumpdata);
void missing_land(dump_data_t *dumpdata);
void cs_fill_land(dump_data_t *dumpdata);
void cs_fill_land_orig(dump_data_t *dumpdata);
void cs_landfill_map(dump_data_t *dumpdata);

landfillfn_t locate_landfill_function(const char *tag)
{
  int i = 0;
  static landfill_map_t standard[] = {
    {"default", def_fill_land},
    {"cascade_search", cs_fill_land},
    {"cascade_search_orig", cs_fill_land_orig},
    {"zero_fill", zero_fill_land},
    {"no_fill", no_fill_land},
    {"missing_value", missing_land},
    {NULL, NULL}
  };
  i = 0;

  while (standard[i].name != NULL) {
    if (strcmp(standard[i].name, tag) == 0)
      return standard[i].landfill;
    ++i;
  }

  hd_quit("Unable to locate the land filling function: %s\n", tag);
  return NULL;
}

void dump_bathy_mask(dump_data_t *dumpdata, double bathyf)
{
  int i, j, k, n;

  /* Filling the surface elevation array. */
  if (bathyf != 9999) {
    k = dumpdata->nz - 1;
    for (j = 0; j < dumpdata->nce2; ++j) {
      for (i = 0; i < dumpdata->nce1; ++i) {
      if (dumpdata->flag[k][j][i] & (SOLID | OUTSIDE))
        dumpdata->eta[j][i] = 0.0;
      }
    }
  }

  /* Filling the 3D water column tracer arrays. */
  if (bathyf != 9999) {
    for (n = 0; n < dumpdata->ntr; ++n) {
      double ***tr = dumpdata->tr_wc[n];
      for (k = 0; k < dumpdata->nz; ++k) {
        for (j = 0; j < dumpdata->nce2; ++j) {
          for (i = 0; i < dumpdata->nce1; ++i) {
	    if (!(dumpdata->flag[k][j][i] & (SOLID | OUTSIDE))) {
	      if (bathyf > 0.0 && fabs(dumpdata->botz[j][i]) > fabs(bathyf))
                tr[k][j][i] = NaN;
	      if (bathyf < 0.0 && fabs(dumpdata->botz[j][i]) < fabs(bathyf))
                tr[k][j][i] = NaN;
            }
          }
        }
      }
    }
    for (n = 0; n < dumpdata->ntrS; ++n) {
      double **tr = dumpdata->tr_wcS[n];
      k = dumpdata->nz - 1;
      for (j = 0; j < dumpdata->nce2; ++j) {
        for (i = 0; i < dumpdata->nce1; ++i) {
	  if (!(dumpdata->flag[k][j][i] & (SOLID | OUTSIDE))) {
	    if (bathyf > 0.0 && fabs(dumpdata->botz[j][i]) > fabs(bathyf))
              tr[j][i] = NaN;
	    if (bathyf < 0.0 && fabs(dumpdata->botz[j][i]) < fabs(bathyf))
              tr[j][i] = NaN;
          }
        }
      }
    }
  }
}

void def_fill_land(dump_data_t *dumpdata)
{
  int i, j, k, n;

  /* Filling the surface elevation array. */
  k = dumpdata->nz - 1;
  for (j = 0; j < dumpdata->nce2; ++j) {
    for (i = 0; i < dumpdata->nce1; ++i) {
      if (dumpdata->flag[k][j][i] & (SOLID | OUTSIDE))
        dumpdata->eta[j][i] = 0.0;
    }
  }

  /* Filling the 3D water column tracer arrays. */
  for (n = 0; n < dumpdata->ntr; ++n) {
    double ***tr = dumpdata->tr_wc[n];
    for (k = 0; k < dumpdata->nz; ++k) {
      for (j = 0; j < dumpdata->nce2; ++j) {
        for (i = 0; i < dumpdata->nce1; ++i) {
          if (dumpdata->flag[k][j][i] & (SOLID | OUTSIDE)) {
            /*UR reverted change made by Jason for the restart capability 
             *  tr[k][j][i] = dumpdata->trinfo_3d[n].valid_range_wc[0]-1;*/
            tr[k][j][i] = NaN;
	    dumpdata->dens_0[k][j][i] = NaN;
          }
        }
      }
    }
  }

}

void no_fill_land(dump_data_t *dumpdata)
{
  return;
}

void zero_fill_land(dump_data_t *dumpdata)
{
  int i, j, k, kk, n;
  /* 
   * Fill the surface cell centred arrays
   */
  k = dumpdata->nz - 1;
  for (j = 0; j < dumpdata->nce2; ++j) {
    for (i = 0; i < dumpdata->nce1; ++i) {
      if (dumpdata->flag[k][j][i] & (SOLID | OUTSIDE) ||
	  isnan(dumpdata->eta[j][i])) {
	dumpdata->eta[j][i] = 0.0;
	dumpdata->patm[j][i] = 0.0;
	dumpdata->wind1[j][i] = 0.0;
	dumpdata->wind2[j][i] = 0.0;
      }
    }
  }
  /* 
   * u1 faces
   */
  for (j = 0; j < dumpdata->nce2; ++j) {
    for (i = 0; i < dumpdata->nfe1; ++i) {
      if (dumpdata->flag[k][j][i] & (SOLID | OUTSIDE) ||
   	  isnan(dumpdata->wind1[j][i])) {
	dumpdata->wind1[j][i] = 0.0;
	dumpdata->u1av[j][i]  = 0.0;
      }
    }
  }
  /* 
   * u2 faces
   */
  for (j = 0; j < dumpdata->nfe2; ++j) {
    for (i = 0; i < dumpdata->nce1; ++i) {
      if (dumpdata->flag[k][j][i] & (SOLID | OUTSIDE) ||
	  isnan(dumpdata->wind2[j][i])) {
	dumpdata->wind2[j][i] = 0.0;
	dumpdata->u2av[j][i]  = 0.0;
      }
    }
  }

  /*
   * Fill the surface (2D) tracers
   */
  for (n = 0; n < dumpdata->ntrS; ++n) {
    double **tr = dumpdata->tr_wcS[n];
    for (j = 0; j < dumpdata->nce2; ++j) {
      for (i = 0; i < dumpdata->nce1; ++i) {
	if (dumpdata->flag[k][j][i] & (SOLID | OUTSIDE) || isnan(tr[j][i]))
	  tr[j][i] = 0.0;
      }
    }
  }
  /*
   * Fill all 3D tracers
   */
  for (n = 0; n < dumpdata->ntr; ++n) {
    double ***tr = dumpdata->tr_wc[n];
    for (k = 0; k < dumpdata->nz; ++k) {
      for (j = 0; j < dumpdata->nce2; ++j) {
        for (i = 0; i < dumpdata->nce1; ++i) {
	  if (dumpdata->flag[k][j][i] & (SOLID | OUTSIDE) || isnan(tr[k][j][i]))
	    tr[k][j][i] = 0.0;
	}
      }
    }
  }

  /* Sediment tracers */
  k = dumpdata->nz - 1;
  for (n = 0; n < dumpdata->nsed; ++n) {
    double ***tr = dumpdata->tr_sed[n];
    for (kk = 0; kk < dumpdata->sednz; ++kk) {
      for (j = 0; j < dumpdata->nce2; ++j) {
        for (i = 0; i < dumpdata->nce1; ++i) {
	  if (dumpdata->flag[k][j][i] & (SOLID | OUTSIDE) || isnan(tr[kk][j][i]))
	    tr[kk][j][i] = 0.0;
	}
      }
    }
  }
}


void missing_land(dump_data_t *dumpdata)
{
  int i, j, k, kk, n;

  /* Filling the surface elevation array. */
  k = dumpdata->nz - 1;
  for (j = 0; j < dumpdata->nce2; ++j) {
    for (i = 0; i < dumpdata->nce1; ++i) {
      if (dumpdata->flag[k][j][i] & (SOLID | OUTSIDE) ||
	  isnan(dumpdata->eta[j][i])) {
	dumpdata->eta[j][i] = 0.0;
	dumpdata->patm[j][i] = 0.0;
	dumpdata->wind1[j][i] = 0.0;
	dumpdata->wind2[j][i] = 0.0;
      }
    }
  }
  for (j = 0; j < dumpdata->nce2; ++j) {
    for (i = 0; i < dumpdata->nfe1; ++i) {
      if (dumpdata->flag[k][j][i] & (SOLID | OUTSIDE) ||
	  isnan(dumpdata->eta[j][i])) {
	dumpdata->wind1[j][i] = 0.0;
      }
    }
  }
  for (j = 0; j < dumpdata->nfe2; ++j) {
    for (i = 0; i < dumpdata->nce1; ++i) {
      if (dumpdata->flag[k][j][i] & (SOLID | OUTSIDE) ||
	  isnan(dumpdata->eta[j][i])) {
	dumpdata->wind2[j][i] = 0.0;
      }
    }
  }

  for (n = 0; n < dumpdata->ntrS; ++n) {
    double **tr = dumpdata->tr_wcS[n];
    for (j = 0; j < dumpdata->nce2; ++j) {
      for (i = 0; i < dumpdata->nce1; ++i) {
	if (dumpdata->flag[k][j][i] & (SOLID | OUTSIDE) || isnan(tr[j][i]))
	  tr[j][i] = dumpdata->trinfo_2d[n].missing_value;
      }
    }
  }

  for (n = 0; n < dumpdata->ntr; ++n) {
    double ***tr = dumpdata->tr_wc[n];
    for (k = 0; k < dumpdata->nz; ++k) {
      for (j = 0; j < dumpdata->nce2; ++j) {
        for (i = 0; i < dumpdata->nce1; ++i) {
	  if (dumpdata->flag[k][j][i] & (SOLID | OUTSIDE) || isnan(tr[k][j][i]))
	    tr[k][j][i] = dumpdata->trinfo_3d[n].missing_value;
	}
      }
    }
  }
  /* Sediments */
  k = dumpdata->nz - 1;
  for (n = 0; n < dumpdata->nsed; ++n) {
    double ***tr = dumpdata->tr_sed[n];
    for (kk = 0; kk < dumpdata->sednz; ++kk) {
      for (j = 0; j < dumpdata->nce2; ++j) {
        for (i = 0; i < dumpdata->nce1; ++i) {
	  if (dumpdata->flag[k][j][i] & (SOLID | OUTSIDE) || isnan(tr[kk][j][i]))
	    tr[kk][j][i] = dumpdata->trinfo_sed[n].missing_value;
	}
      }
    }
  }
}

void cs_landfill_map_custom(dump_data_t *dumpdata)
{
  int i, j, k;
  int ci;
  int cj;
  int last_layer_solid = 1;

  for (k = 0; k < dumpdata->nz; ++k) {
    if (last_layer_solid && is_layer_solid(dumpdata, k)) {
      continue;
    } else {
      last_layer_solid = 0;
      for (j = 0; j < dumpdata->nce2; ++j) {
        for (i = 0; i < dumpdata->nce1; ++i) {
          if (dumpdata->flag[k][j][i] & (SOLID | OUTSIDE) ||
	      isnan(dumpdata->tr_wcS[0][j][i])) {
            if (find_closest_nonnan(dumpdata, dumpdata->tr_wcS[0], i, j, k, &ci, &cj)) {
	      dumpdata->fmap_i[k][j][i] = ci;
	      dumpdata->fmap_j[k][j][i] = cj;	
	    } else {
	      dumpdata->fmap_i[k][j][i] = -1;
	      dumpdata->fmap_j[k][j][i] = -1;	
	    }	      
	  }
	}
      }
    }
  }
}


void cs_landfill_map(dump_data_t *dumpdata)
{
  int i, j, k;
  int ci;
  int cj;
  int last_layer_solid = 1;

  for (k = 0; k < dumpdata->nz; ++k) {
    if (last_layer_solid && is_layer_solid(dumpdata, k)) {
      continue;
    } else {
      last_layer_solid = 0;
#if defined(HAVE_OMP)
#pragma omp parallel for private(i,ci,cj)
#endif
      for (j = 0; j < dumpdata->nce2; ++j) {
        for (i = 0; i < dumpdata->nce1; ++i) {
          if (dumpdata->flag[k][j][i] & (SOLID | OUTSIDE)) {
            if (find_closest(dumpdata, i, j, k, &ci, &cj)) {
	      dumpdata->fmap_i[k][j][i] = ci;
	      dumpdata->fmap_j[k][j][i] = cj;	
	    } else {
	      dumpdata->fmap_i[k][j][i] = -1;
	      dumpdata->fmap_j[k][j][i] = -1;	
	    }	      
	  }
	}
      }
    }
  }
}




void cs_fill_land(dump_data_t *dumpdata)
{
  int i, j, k, kk, n;
  int ci;
  int cj;
  int last_layer_solid = 1;

  /* Get the landfill map                                           */
  if (dumpdata->fmap_i == NULL && dumpdata->fmap_j == NULL) {
    dumpdata->fmap_i = i_alloc_3d(dumpdata->nce1,dumpdata->nce2,dumpdata->nz);
    dumpdata->fmap_j = i_alloc_3d(dumpdata->nce1,dumpdata->nce2,dumpdata->nz);
    cs_landfill_map(dumpdata);
  }

  /* Filling the 2D arrays */
  k = dumpdata->nz - 1;
  for (j = 0; j < dumpdata->nce2; ++j) {
    for (i = 0; i < dumpdata->nce1; ++i) {
      if (dumpdata->flag[k][j][i] & (SOLID | OUTSIDE)) {
	ci = dumpdata->fmap_i[k][j][i];
	cj = dumpdata->fmap_j[k][j][i];
        if (ci >= 0 && cj >= 0) {
          dumpdata->eta[j][i] = dumpdata->eta[cj][ci];
          dumpdata->u1av[j][i] = dumpdata->u1av[cj][ci];
          dumpdata->u2av[j][i] = dumpdata->u2av[cj][ci];
          dumpdata->patm[j][i] = dumpdata->patm[cj][ci];
          dumpdata->wind1[j][i] = dumpdata->wind1[cj][ci];
          dumpdata->wind2[j][i] = dumpdata->wind2[cj][ci];
	  for (n = 0; n < dumpdata->ntrS; ++n)
	    dumpdata->tr_wcS[n][j][i] = dumpdata->tr_wcS[n][cj][ci];
	}
        else {
	  /* Should this be an error? */
          dumpdata->eta[j][i] = 0;
          dumpdata->u1av[j][i] = 0.0;
          dumpdata->u2av[j][i] = 0.0;
	  dumpdata->patm[j][i] = 0.0;
	  dumpdata->wind1[j][i] = 0.0;
	  dumpdata->wind2[j][i] = 0.0;
	  for (n = 0; n < dumpdata->ntrS; ++n)
	    dumpdata->tr_wcS[n][j][i] = 0.0;
	}
	if (isnan(dumpdata->eta[j][i])) dumpdata->eta[j][i] = 0.0;
      }

      // last e1 face value
      if (i == dumpdata->nce1-1 && dumpdata->flag[k][j][i] & (SOLID | OUTSIDE))
	dumpdata->u1av[j][i+1] = dumpdata->u1av[j][i];

      // last e2 face value
      if (j == dumpdata->nce2-1 && dumpdata->flag[k][j][i] & (SOLID | OUTSIDE))
	dumpdata->u2av[j+1][i] = dumpdata->u2av[j][i];
    }
  }

  for (j = 0; j < dumpdata->nce2; ++j)
    dumpdata->wind1[j][dumpdata->nce1] = dumpdata->wind1[j][dumpdata->nce1-1];
  for (i = 0; i < dumpdata->nce1; ++i)
    dumpdata->wind2[dumpdata->nce2][i] = dumpdata->wind1[dumpdata->nce2-1][i];

  /* Filling the 3D water column tracer arrays. */
  for (k = 0; k < dumpdata->nz; ++k) {
    if (last_layer_solid && is_layer_solid(dumpdata, k)) {
      for (j = 0; j < dumpdata->nce2; ++j)
        for (i = 0; i < dumpdata->nce1; ++i)
          for (n = 0; n < dumpdata->ntr; ++n) {
            dumpdata->tr_wc[n][k][j][i]
              = dumpdata->trinfo_3d[n].fill_value_wc;
	    // Reset velocities
	    dumpdata->u1[k][j][i] = 0.0;
	    dumpdata->u2[k][j][i] = 0.0;
	  }

    } else {
      last_layer_solid = 0;

      for (j = 0; j < dumpdata->nce2; ++j) {
        for (i = 0; i < dumpdata->nce1; ++i) {
          if (dumpdata->flag[k][j][i] & (SOLID | OUTSIDE)) {
	    ci = dumpdata->fmap_i[k][j][i];
	    cj = dumpdata->fmap_j[k][j][i];
            for (n = 0; n < dumpdata->ntr; ++n) {
	      dumpdata->tr_wc[n][k][j][i] = dumpdata->tr_wc[n][k][cj][ci];
	    }
	    dumpdata->u1[k][j][i] = dumpdata->u1[k][cj][ci];
	    dumpdata->u2[k][j][i] = dumpdata->u2[k][cj][ci];
          }
	  for (n = 0; n < dumpdata->ntr; ++n) {
	    if (isnan(dumpdata->tr_wc[n][k][j][i]))
	      dumpdata->tr_wc[n][k][j][i] = dumpdata->trinfo_3d[n].fill_value_wc;
	  }
	  
	  if (isnan(dumpdata->u1[k][j][i]))
	    dumpdata->u1[k][j][i] = 0.0;
	  if (isnan(dumpdata->u2[k][j][i]))
	    dumpdata->u2[k][j][i] = 0.0;

	  // last e1 face value
	  if (i == dumpdata->nce1-1)
	    dumpdata->u1[k][j][i+1] = dumpdata->u2[k][j][i];
	  
	  // last e2 face value
	  if (j == dumpdata->nce2-1)
	    dumpdata->u2[k][j+1][i] = dumpdata->u1[k][j][i];

        }
      }
    }
  }
  /* Fill sediment tracers */
  k = dumpdata->nz - 1;
  for (j = 0; j < dumpdata->nce2; ++j) {
    for (i = 0; i < dumpdata->nce1; ++i) {
      if (dumpdata->flag[k][j][i] & (SOLID | OUTSIDE)) {
	ci = dumpdata->fmap_i[k][j][i];
	cj = dumpdata->fmap_j[k][j][i];
	if (ci >=0 && cj >= 0) {
	  for (n = 0; n < dumpdata->nsed; ++n) {
	    for (kk = 0; kk < dumpdata->sednz; ++kk)
	      dumpdata->tr_sed[n][kk][j][i] = dumpdata->tr_sed[n][kk][cj][ci];
	  }
	} else
	  for (n = 0; n < dumpdata->nsed; ++n) {
	    for (kk = 0; kk < dumpdata->sednz; ++kk)
	      dumpdata->tr_sed[n][kk][j][i] = 
		               dumpdata->trinfo_sed[n].fill_value_sed;
	  }
      }
    }
  }
}

void cs_fill_land_orig(dump_data_t *dumpdata)
{
  int i, j, k, n;
  int ci;
  int cj;
  int last_layer_solid = 1;


  /* Filling the surface elevation array. */
  k = dumpdata->nz - 1;
  for (j = 0; j < dumpdata->nce2; ++j) {
    for (i = 0; i < dumpdata->nce1; ++i) {
      if (dumpdata->flag[k][j][i] & (SOLID | OUTSIDE)) {
        if (find_closest(dumpdata, i, j, k, &ci, &cj))
          dumpdata->eta[j][i] = dumpdata->eta[cj][ci];
        else
          dumpdata->eta[j][i] = 0;
      }
    }
  }

  /* Filling the water column tracer arrays. */
  for (k = 0; k < dumpdata->nz; ++k) {
    if (last_layer_solid && is_layer_solid(dumpdata, k)) {
      for (j = 0; j < dumpdata->nce2; ++j)
        for (i = 0; i < dumpdata->nce1; ++i)
          for (n = 0; n < dumpdata->ntr; ++n)
            dumpdata->tr_wc[n][k][j][i]
              = dumpdata->trinfo_3d[n].fill_value_wc;
    } else {
      last_layer_solid = 0;

      for (j = 0; j < dumpdata->nce2; ++j) {
        for (i = 0; i < dumpdata->nce1; ++i) {
          if (dumpdata->flag[k][j][i] & (SOLID | OUTSIDE)) {
            find_closest(dumpdata, i, j, k, &ci, &cj);
            for (n = 0; n < dumpdata->ntr; ++n)
              dumpdata->tr_wc[n][k][j][i] = dumpdata->tr_wc[n][k][cj][ci];
          }
        }
      }
    }
  }
}

int is_layer_solid(dump_data_t *dumpdata, int k)
{
  int i, j;
  for (j = 0; j < dumpdata->nce2; ++j)
    for (i = 0; i < dumpdata->nce1; ++i)
      if (!(dumpdata->flag[k][j][i] & (SOLID | OUTSIDE)))
        return 0;

  return 1;
}


/* Locate the closest cell to oi, oj that contains valid data.
 * This is not very efficient but so what, we don't run this
 * program very often.
 */
int find_closest(dump_data_t *dumpdata, int oi, int oj, int k, int *ci,
                 int *cj)
{
  int i, j;
  double mindist;
  int level;

  mindist = 1e38;
  level = 1;
  *ci = -1;
  *cj = -1;
  while (*ci < 0) {
    int finishedLevel = 0;
    int edge = 0;

    /* Scanned the whole grid, must be solid every where. */
    if (level > dumpdata->nce1 && level > dumpdata->nce2)
      return 0;

    while (!finishedLevel) {

      int jfrom = oj - level;
      int jto = oj + level;
      int ifrom = oi - level;
      int ito = oi + level;

      switch (edge) {
      case 0:                  /* Left edge */
        ito = ifrom;
        break;

      case 1:                  /* Bottom edge */
        ++ifrom;
        jto = jfrom;
        break;

      case 2:                  /* Right edge */
        ifrom = ito;
        ++jfrom;
        break;

      case 3:                  /* Top edge */
        ++ifrom;
        jfrom = jto;
        finishedLevel = 1;
        break;
      }
      ++edge;

      for (j = jfrom; j <= jto; ++j) {
        for (i = ifrom; i <= ito; ++i) {
          if ((i < 0) || (j < 0) || (i >= dumpdata->nce1) ||
              (j >= dumpdata->nce2))
            continue;
          if (!(dumpdata->flag[k][j][i] & (SOLID | OUTSIDE))) {
            double dx = dumpdata->cellx[oj][oi] - dumpdata->cellx[j][i];
            double dy = dumpdata->celly[oj][oi] - dumpdata->celly[j][i];
            double dist = dx * dx + dy * dy;
	    /* Curvilinear grids have cellx & celly = nan at locations */
	    /* not included in curvilinear index space - use the */
	    /* minimum distance based on (i,j) location in this case. */
	    if (isnan(dumpdata->cellx[oj][oi]) || 
		isnan(dumpdata->celly[oj][oi]))
	      dist = (oi - i) * (oi - i) + (oj - j) * (oj - j);
            if (dist < mindist) {
              *ci = i;
              *cj = j;
              mindist = dist;
            }
          }
        }
      }
    }
    ++level;
  }

  return 1;
}



int find_closest_nonnan(dump_data_t *dumpdata, double **a, 
			int oi, int oj, int k, int *ci,
			int *cj)
{
  int i, j;
  double mindist;
  int level;

  mindist = 1e38;
  level = 1;
  *ci = -1;
  *cj = -1;
  while (*ci < 0) {
    int finishedLevel = 0;
    int edge = 0;

    /* Scanned the whole grid, must be solid every where. */
    if (level > dumpdata->nce1 && level > dumpdata->nce2)
      return 0;

    while (!finishedLevel) {

      int jfrom = oj - level;
      int jto = oj + level;
      int ifrom = oi - level;
      int ito = oi + level;

      switch (edge) {
      case 0:                  /* Left edge */
        ito = ifrom;
        break;

      case 1:                  /* Bottom edge */
        ++ifrom;
        jto = jfrom;
        break;

      case 2:                  /* Right edge */
        ifrom = ito;
        ++jfrom;
        break;

      case 3:                  /* Top edge */
        ++ifrom;
        jfrom = jto;
        finishedLevel = 1;
        break;
      }
      ++edge;

      for (j = jfrom; j <= jto; ++j) {
        for (i = ifrom; i <= ito; ++i) {
          if ((i < 0) || (j < 0) || (i >= dumpdata->nce1) ||
              (j >= dumpdata->nce2))
            continue;

          if (!(dumpdata->flag[k][j][i] & (SOLID | OUTSIDE)) &&
	      !(isnan(a[j][i]))) {
            double dx = dumpdata->cellx[oj][oi] - dumpdata->cellx[j][i];
            double dy = dumpdata->celly[oj][oi] - dumpdata->celly[j][i];
            double dist = dx * dx + dy * dy;
	    /* Curvilinear grids have cellx & celly = nan at locations */
	    /* not included in curvilinear index space - use the */
	    /* minimum distance based on (i,j) location in this case. */
	    if (isnan(dumpdata->cellx[oj][oi]) || 
		isnan(dumpdata->celly[oj][oi]))
	      dist = (oi - i) * (oi - i) + (oj - j) * (oj - j);
            if (dist < mindist) {
              *ci = i;
              *cj = j;
              mindist = dist;
            }
          }
        }
      }
    }
    ++level;
  }

  return 1;
}


/* Locate the closest cell to oi, oj that contains valid data.
 * This is not very efficient but so what, we don't run this
 * program very often.
 */
int find_closest_w(geometry_t *geom, int oi, int oj, int *ci, int *cj)
                 
{
  int i, j, k = geom->nz-1;
  int c, co = geom->map[k][oj][oi];
  double mindist;
  int level;

  mindist = 1e38;
  level = 1;
  *ci = -1;
  *cj = -1;
  while (*ci < 0) {
    int finishedLevel = 0;
    int edge = 0;

    /* Scanned the whole grid, must be solid every where. */
    if (level > geom->nce1 && level > geom->nce2)
      return 0;

    while (!finishedLevel) {

      int jfrom = oj - level;
      int jto = oj + level;
      int ifrom = oi - level;
      int ito = oi + level;

      switch (edge) {
      case 0:                  /* Left edge */
        ito = ifrom;
        break;

      case 1:                  /* Bottom edge */
        ++ifrom;
        jto = jfrom;
        break;

      case 2:                  /* Right edge */
        ifrom = ito;
        ++jfrom;
        break;

      case 3:                  /* Top edge */
        ++ifrom;
        jfrom = jto;
        finishedLevel = 1;
        break;
      }
      ++edge;

      for (j = jfrom; j <= jto; ++j) {
        for (i = ifrom; i <= ito; ++i) {
          if ((i < 0) || (j < 0) || (i >= geom->nce1) ||
              (j >= geom->nce2))
            continue;
	  c = geom->map[k][j][i];
          if (c > 0 && c <= geom->ewetS) {
            double dx = geom->cellx[co] - geom->cellx[c];
            double dy = geom->celly[co] - geom->celly[c];
            double dist = dx * dx + dy * dy;
	    /* Curvilinear grids have cellx & celly = nan at locations */
	    /* not included in curvilinear index space - use the */
	    /* minimum distance based on (i,j) location in this case. */
	    if (!co || isnan(geom->cellx[co]) || 
		isnan(geom->celly[co]))
	      dist = (oi - i) * (oi - i) + (oj - j) * (oj - j);
            if (dist < mindist) {
              *ci = i;
              *cj = j;
              mindist = dist;
            }
          }
        }
      }
    }
    ++level;
  }

  return(geom->map[k][*cj][*ci]);
}


/* Locate the closest cell to lat, lon that contains valid data.
 * This is not very efficient but so what, we don't run this
 * program very often.
 */
int find_closest_b(parameters_t *params, double **bathy, double lat, double lon, int *ci, int *cj, int mode)
{
  int i, j, ii, jj, io, jo, nz = params->nz;
  double mindist = 1e36;
  double minlat, minlon, maxlat, maxlon;
  double dx, dy, dist;

  /* Get the grid bounds */
  minlon = minlat = 360.0;
  maxlon = maxlat = -360.0;
  for (i = 0 ; i < params->nce1; i++)
    for (j = 0 ; j < params->nce2; j++) {
      if (!isnan(params->x[2*j+1][2*i+1])) {
	minlon = min(minlon, params->x[2*j+1][2*i+1]);
	maxlon = max(maxlon, params->x[2*j+1][2*i+1]);
      }
      if (!isnan(params->y[2*j+1][2*i+1])) {
	minlat = min(minlat, params->y[2*j+1][2*i+1]);
	maxlat = max(maxlat, params->y[2*j+1][2*i+1]);
      }
    }

  *ci = -1; *cj = -1;
  if (lon < minlon || lon > maxlon || lat < minlat || lat > maxlat)
    return(1);

  /* Get the (i,j) of (lon,lat) */
  for (i = 0 ; i < params->nce1; i++) {
    for (j = 0 ; j < params->nce2; j++) {
      ii = 2 * i + 1;
      jj = 2 * j + 1;
      if (isnan(params->x[jj][ii]) || isnan(params->y[jj][ii]))
	continue;
      dx = lon - params->x[jj][ii];
      dy = lat - params->y[jj][ii];
      dist = dx * dx + dy * dy;
      if (dist < mindist) {
	io = i;
	jo = j;
	mindist = dist;
      }
    }
  }
  if (mode) {
    *ci = io;
    *cj = jo;
    return(0);
  }
  /* Find the nearest land cell if (lon,lat) is a wet cell */
  if (bathy[jo][io] < params->layers[nz] && bathy[jo][io] < params->etamax &&
      bathy[jo][io] != NOTVALID) {
    mindist = 1e36;
    for (i = 0 ; i < params->nce1; i++) {
      for (j = 0 ; j < params->nce2; j++) {
	ii = 2 * i + 1;
	jj = 2 * j + 1;
	if (isnan(params->x[jj][ii]) || isnan(params->y[jj][ii]))
	  continue;
	if (bathy[j][i] < params->layers[nz] && bathy[j][i] < params->etamax &&
	    bathy[j][i] != NOTVALID)
	  continue;
	dx = lon - params->x[jj][ii];
	dy = lat - params->y[jj][ii];
	dist = dx * dx + dy * dy;
	if (dist < mindist) {
	  *ci = i;
	  *cj = j;
	  mindist = dist;
	}
      }
    }
  } else {
    double nlat, nlon;
    /* Find the nearest wet cell */
    mindist = 1e36;
    for (i = 0 ; i < params->nce1; i++) {
      for (j = 0 ; j < params->nce2; j++) {
	ii = 2 * i + 1;
	jj = 2 * j + 1;
	if (isnan(params->x[jj][ii]) || isnan(params->y[jj][ii]))
	  continue;
	if (bathy[j][i] >= params->layers[nz] || bathy[j][i] >= params->etamax ||
	    bathy[j][i] == NOTVALID)
	  continue;
	dx = lon - params->x[jj][ii];
	dy = lat - params->y[jj][ii];
	dist = dx * dx + dy * dy;
	if (dist < mindist) {
	  nlon = params->x[jj][ii];
	  nlat = params->y[jj][ii];
	  mindist = dist;
	}
      }
    }
    lat = nlat;
    lon = nlon;
    /* Find the nearest land cell */
    mindist = 1e36;
    for (i = 0 ; i < params->nce1; i++) {
      for (j = 0 ; j < params->nce2; j++) {
	ii = 2 * i + 1;
	jj = 2 * j + 1;
	if (isnan(params->x[jj][ii]) || isnan(params->y[jj][ii]))
	  continue;
	if (bathy[j][i] < params->layers[nz] && bathy[j][i] < params->etamax &&
	    bathy[j][i] != NOTVALID)
	  continue;
	dx = lon - params->x[jj][ii];
	dy = lat - params->y[jj][ii];
	dist = dx * dx + dy * dy;
	if (dist < mindist) {
	  *ci = i;
	  *cj = j;
	  mindist = dist;
	}
      }
    }
  }
  if (*ci == -1 || *cj == -1)
    return(1);
  return(0);
}
