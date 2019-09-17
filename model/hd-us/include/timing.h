/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/include/timing.h
 *  
 *  Description:
 *  Macros to do profile timing
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: timing.h 6187 2019-03-25 06:11:22Z riz008 $
 *
 */

#ifndef _TIMING_H
#define _TIMING_H

#include <sys/time.h>

/*
 * Macros and function prototypes are only defined when DO_TIMING is
 * set
 */
#ifdef DO_TIMING
/*
 * Initialises some timing globals
 */
void init_timing(void);
void print_timing_header(char *s, void *m);

#define INIT_TIMING init_timing();
#define END_TIMING fclose(tfp);
#define FLUSH_TIMING if (tfp != NULL) fflush(tfp);
#define TIMING_LEVEL_DEF 3
/*
 * Consider making this an environment variable
 */
#define TIMING_NUM_STEPS 25

#define TIMING_PRINT(x) \
         fprintf(tfp, x);

#define PRINT_TIMING_HEADER(s,m) print_timing_header(s,m);

/*
 * The opening brace is closed out in _TIMING_DUMP
 */
#define TIMING_SET \
         { \
           struct timeval tm1, tm2; \
           gettimeofday(&tm1, NULL);

#define TIMING_PAD0 " \t"
#define TIMING_PAD1 " \t\t"
#define TIMING_PAD2 " \t\t\t"
#define TIMING_PAD3 " \t\t\t\t"
#define TIMING_PAD4 " \t\t\t\t\t"

/*
 * The final brace closes out the one in TIMING_SET
 */
#define TIMING_DUMP(pad,x) \
           if (pad <= timing_level) {   \
             gettimeofday(&tm2, NULL);  \
             fprintf(tfp, " %-15s\t: %s%.5f\n", x , TIMING_PAD##pad, \
		       ((tm2.tv_sec  - tm1.tv_sec) +       \
		       (tm2.tv_usec - tm1.tv_usec)*1e-6)); \
             } \
         }

#define TIMING_DUMP_WIN(pad,x,win) \
           if (pad <= timing_level) {   \
             gettimeofday(&tm2, NULL);  \
             fprintf(tfp, " %-15s\t: %s%.5f(%i)\n", x , TIMING_PAD##pad, \
		       ((tm2.tv_sec  - tm1.tv_sec) +            \
		       (tm2.tv_usec - tm1.tv_usec)*1e-6), win); \
             } \
         }

#define TIMING_COUNTER \
         if ((timing_counter++ % TIMING_NUM_STEPS) == 0) { \
              (void)fseek(tfp, tfp_pos, SEEK_SET); \
         }
#else
// define empty
#define PRINT_TIMING_HEADER
#define INIT_TIMING
#define END_TIMING
#define TIMING_SET
#define TIMING_DUMP(x,y)
#define TIMING_PRINT(x,...)
#define TIMING_COUNTER
#define TIMING_DUMP_WIN(x,y,z)
#define FLUSH_TIMING
#endif


#endif                          /* _TIMING_H */

