/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: lib/ecology/ediagn.c
 *  
 *  Description:
 *  A simple diagnostic utility for ecology.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: ediagn.c 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include "ecology.h"
#include "utils.h"
#include "einterface.h"

static int verbose = 0;
static int printtracers = 0;
static char* tunits = "";

static void usage()
{
    printf("Usage: ediagn <prm file> [-v]\n");
    printf("       ediagn -v\n");
    printf("Options:\n");
    printf("  -t -- print list of tracers\n");
    printf("  -v -- verbose / version\n");

    exit(0);
}

static void version()
{
    printf("ediagn/libecology version %s\n", ecology_version);
    exit(0);
}

static void parse_commandline(int argc, char* argv[], char** fname)
{
    int i;

    if (argc == 2 && strcmp(argv[1], "-v") == 0)
        version();

    i = 1;
    while (i < argc) {
        if (argv[i][0] != '-') {
            if (*fname == NULL)
                *fname = argv[i];
            else
                usage();
            i++;
        } else
            switch (argv[i][1]) {
            case 't':
                i++;
                printtracers = 1;
                break;
            case 'v':
                i++;
                verbose = 1;
                break;
            default:
                usage();
                break;
            }
    }

    if (*fname == NULL)
        usage();
}

int einterface_getverbosity(void* model)
{
    return verbose;
}

char* einterface_gettimeunits(void* model)
{
    return tunits;
}

int einterface_getntracers(void* model)
{
    e_quit("must have \"internal_tracers 1\" in the parameter file for ediagn to run\n");
    return 0;
}

char* einterface_gettracername(void* model, int i)
{
    e_quit("programming error: not supposed to get here\n");
    return NULL;
}

int einterface_gettracerdiagnflag(void* model, char* name)
{
    return ecology_getdiagnflag(name);
}

int einterface_getnepis(void* model)
{
    e_quit("programming error: not supposed to get here\n");
    return 0;
}

char* einterface_getepiname(void* model, int i)
{
    e_quit("programming error: not supposed to get here\n");
    return NULL;
}

int einterface_getepidiagnflag(void* model, char* name)
{
    return ecology_getdiagnflag(name);
}

int einterface_getnumbercolumns(void* model)
{
    return 0;
}

int einterface_getnumberwclayers(void* model)
{
    return 0;
}

int einterface_getnumbersedlayers(void* model)
{
    return 0;
}

int einterface_isboundarycolumn(void* model, int b)
{
    e_quit("programming error: not supposed to get here\n");
    return -1;
}

int einterface_getwctopk(void* model, int b)
{
    e_quit("programming error: not supposed to get here\n");
    return -1;
}

int einterface_getwcbotk(void* model, int b)
{
    e_quit("programming error: not supposed to get here\n");
    return -1;
}

int einterface_getsedtopk(void* model, int b)
{
    e_quit("programming error: not supposed to get here\n");
    return -1;
}

int einterface_getsedbotk(void* model, int b)
{
    e_quit("programming error: not supposed to get here\n");
    return -1;
}

double* einterface_getwccellthicknesses(void* model, int b)
{
    e_quit("programming error: not supposed to get here\n");
    return NULL;
}

double* einterface_getsedcellthicknesses(void* model, int b)
{
    e_quit("programming error: not supposed to get here\n");
    return NULL;
}

double** einterface_getwctracers(void* model, int b)
{
    e_quit("programming error: not supposed to get here\n");
    return NULL;
}

double** einterface_getsedtracers(void* model, int b)
{
    e_quit("programming error: not supposed to get here\n");
    return NULL;
}

double** einterface_getepivars(void* model, int b)
{
    e_quit("programming error: not supposed to get here\n");
    return NULL;
}

double* einterface_getporosity(void* model, int b)
{
    e_quit("programming error: not supposed to get here\n");
    return NULL;
}

void* einterface_gete_quitfn()
{
    return NULL;
}

void einterface_ecologyinit(void* model, void* ecology)
{
}

double einterface_getmodeltime(void* model)
{
    e_quit("programming error: not supposed to get here\n");
    return 0.0 / 0.0;
}

double einterface_getustrcw(void* model, int b)
{
    e_quit("programming error: not supposed to get here\n");
    return 0.0 / 0.0;
}

double einterface_getlighttop(void* model, int b)
{
    e_quit("programming error: not supposed to get here\n");
    return 0.0 / 0.0;
}

quitfntype einterface_getquitfn()
{
    return NULL;
}

int einterface_gettracerparticflag(void* model, char* name)
{
    return 0;
}

int main(int argc, char* argv[])
{
    char* fname = NULL;
    ecology* e = NULL;

    parse_commandline(argc, argv, &fname);

    e = ecology_build(NULL, fname);

    if (printtracers) {
        int n = ecology_getntracers(e);
        int i;

        emslog(LDEBUG, "tracers = ");
        if (!verbose)
            for (i = 0; i < n; ++i)
                emslog(LDEBUG, "%s ", ecology_gettracername(e, i));
        else
            for (i = 0; i < n; ++i) {
                char* name = ecology_gettracername(e, i);

                emslog(LDEBUG, "\n%4d %-9s %d", i, name, ecology_getdiagnflag(name));
            }
        emslog(LDEBUG, "\n");

        n = ecology_getnepis(e);
        emslog(LDEBUG, "epis = ");
        if (!verbose)
            for (i = 0; i < n; ++i)
                emslog(LDEBUG, "%s ", ecology_getepiname(e, i));
        else
            for (i = 0; i < n; ++i) {
                char* name = ecology_getepiname(e, i);

                emslog(LDEBUG, "\n%4d %-9s %d", i, name, ecology_getdiagnflag(name));
            }
        emslog(LDEBUG, "\n");
    }

    ecology_destroy(e);

    return 0;
}
