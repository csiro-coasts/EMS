/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/misc/standalone.c
 *
 *  \brief Routines for non EMS standalone applications
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: standalone.c 5831 2018-06-26 23:48:06Z riz008 $
 */


/**
 * Generic quit function
 */
void quit(char* format, ...)
{
    va_list args;

    fflush(stdout);
    fprintf(stderr, "\nerror: grid library: ");
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);

    exit(1);
}

// EOF

