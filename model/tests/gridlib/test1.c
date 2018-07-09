/*
 * Simple test of the grid library
 *
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include "gridlib.h"

int main(int argc, char *argv[])
{
  int num = 16, i;
  /*
   * Define a regular grid
   */
  double x[] = {0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3};
  double y[] = {0,1.2,2,3,0,1.1,2,3.3,0,1.44,2.3,3.8,0,1.9,2,3};
  double z[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  GRID_SPECS *gs = NULL;
  double xy[] = {-1,-1};
  double pinterp = -99;
  FILE *fp;

  // Make truly random
  //  srand(now);

  for (i=0; i<num; i++) {
    // Get a value between 0 and 10
    z[i] = rand()/(double)RAND_MAX * 10;
  }

  // Initialise the grid
  gs = grid_interp_init(x, y, z, num, GRID_NN_NONSIBSONIAN);

  // Get the x & y from the command line
  xy[0] = atof(argv[1]);
  xy[1] = atof(argv[2]);

  pinterp = grid_interp_on_point(gs, xy[0], xy[1]);

  // Print out
  fp = fopen("foo.txt", "w");
  for (i=0; i<num; i++) {
    fprintf(fp,"%f %f %f\n", x[i],y[i],z[i]);
  };
  fprintf(fp,"%f %f %f\n", xy[0],xy[1],pinterp);
  fprintf(fp,"\n");
  fclose(fp);

  // destroy to make sure the function works okay
  grid_specs_destroy(gs);

  return(0);
}


// EOF
