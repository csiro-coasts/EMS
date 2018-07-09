/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file misc/colourtable.c
 *  
 *  \brief Colour table routines
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: colourtable.c 5831 2018-06-26 23:48:06Z riz008 $
 */

#include <stdio.h>
#include <stdlib.h>
#include "ems.h"


/** Opens and reads a colour table file, and returns an
  * instance of a colour_table_t structure.
  * @param fname colour_table_t file name.
  * @return pointer to a colour_table_t structure.
  */
colour_table_t *ct_read(char *fname)
{
  FILE *fp;
  int i;
  int ntvals;
  char line[MAXLINELEN];
  colour_table_t *p;

  if ((fp = fopen(fname, "r")) == NULL)
    return (NULL);

  /* Count the number of table entries */
  ntvals = 0;
  while (fgets(line, MAXLINELEN, fp) != NULL)
    if (line[0] != '#' && line[0] != 0)
      ntvals++;
  if (ntvals < 1) {
    fclose(fp);
    return (NULL);
  }

  /* Allocate memory for table */
  if ((p = (colour_table_t *)malloc(sizeof(colour_table_t))) == NULL) {
    fclose(fp);
    return (NULL);
  }
  p->v = d_alloc_1d(ntvals);
  p->r = d_alloc_1d(ntvals);
  p->g = d_alloc_1d(ntvals);
  p->b = d_alloc_1d(ntvals);

  p->n = ntvals;

  /* Read table */
  fseek(fp, 0L, 0);
  for (i = 0; fgets(line, MAXLINELEN, fp) != NULL;) {
    double v, r, g, b;

    if (line[0] != '#' && line[0] != 0) {
      if (sscanf(line, "%lf %lf %lf %lf", &v, &r, &g, &b) != 4)
        quit("ct_read: Can't read table values\n");
      p->v[i] = v;
      p->r[i] = r;
      p->g[i] = g;
      p->b[i] = b;
      i++;
    }
  }

  /* Check values are monotonic increasing */
  for (i = 1; i < ntvals; i++)
    if (p->v[i] <= p->v[i - 1])
      quit("ct_read: Table out of order\n");

  /* Close file and return colour table pointer */
  fclose(fp);
  return (p);
}

void ct_write(void)
{
  /* NOT IMPLEMENTED YET */
}


/** Get an RGB colour from the colour_table_t.
  *
  * Get the colour (r, g, b values) corresponding to
  * a particular value from the colour table.
  *
  * @param v A value within the colour_table_t range.
  * @param ct Colour table.
  * @param r Red corresponding to specified value.
  * @param g Green corresponding to specified value.
  * @param b Blue corresponding to specified value.
  */
void
ct_get_RGB(double v, colour_table_t ct, double *r, double *g, double *b)
{
  double frac;
  int ilow;
  int imid;
  int ihigh;

  /* first check whether v is within the table range */
  if (v <= ct.v[0]) {
    *r = ct.r[0];
    *g = ct.g[0];
    *b = ct.b[0];
    return;
  }
  if (v >= ct.v[ct.n - 1]) {
    *r = ct.r[ct.n - 1];
    *g = ct.g[ct.n - 1];
    *b = ct.b[ct.n - 1];
    return;
  }

  /* perform binary chop to determine values either side of v */
  ilow = 0;
  ihigh = ct.n - 1;
  while (ihigh - ilow > 1) {
    imid = (ilow + ihigh) / 2;
    if (v >= ct.v[imid])
      ilow = imid;
    else
      ihigh = imid;
  }

  /* Calculate fractional position */
  frac = (v - ct.v[ilow]) / (ct.v[ihigh] - ct.v[ilow]);

  if (ct.n >= 100) {
    /* High-resolution colour table - return the nearest colour. This
       allows an application to use a known set of colours */
    int index = (frac <= 0.5) ? ilow : ihigh;
    *r = ct.r[index];
    *g = ct.g[index];
    *b = ct.b[index];
  } else {
    /* Interpolate red, green and blue values */
    *r = ct.r[ilow] * (1.0 - frac) + ct.r[ihigh] * frac;
    *g = ct.g[ilow] * (1.0 - frac) + ct.g[ihigh] * frac;
    *b = ct.b[ilow] * (1.0 - frac) + ct.b[ihigh] * frac;
  }

  return;
}


#if TEST
int main(int argc, char *argv[])
{
  colour_table_t *ct;
  double start;
  double stop;
  double step;
  double s;
  double r = 0.0;
  double g = 0.0;
  double b = 0.0;
  double minv;
  double maxv;

  if (argc != 5) {
    fprintf(stderr, "Usage: %s colour-table start stop step\n", argv[0]);
    exit(-1);
  }

  ct = ct_read(argv[1]);
  start = atof(argv[2]);
  stop = atof(argv[3]);
  step = atof(argv[4]);
  minv = ct->v[0];
  maxv = ct->v[ct->n - 1];

  s = start;
  while (s <= stop) {
    ct_get_RGB(((s - start) / (stop - start)) * (maxv - minv) + minv, *ct,
               &r, &g, &b);
    printf("%g %g %g %g\n", s, r, g, b);
//      printf("%g %d %d %d\n", s, ((int)r*255), ((int)g*255), ((int)b*255));
    s += step;
  }

  exit(0);
}
#endif
