/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: lib/ecology/eprocess.c
 *  
 *  Description:
 *  Process servicing code -- implementation
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: eprocess.c 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "ecology_internal.h"
#include "stringtable.h"
#include "utils.h"
#include "eprocess.h"
#include "parameter_info.h"

#define STRBUFSIZE 2048

/* these tags identify process groups in the process parameter file
 */
char* eprocess_type_tags[] = {
    "generic",
    "water",
    "sediment",
    "epibenthos"
};

/** Removes heading and trailing garbage.
 */
static char* clear_string(char* s)
{
    int len;

    while (isspace((int)s[0]))
        s++;
    while ((len = strlen(s)) > 0 && isspace((int)s[len - 1]))
        s[len - 1] = 0;
    return s;
}

/** Prints recreated process entry string to a string buffer.
 * @param p Process
 * @param buf Buffer to print in
 */
static void eprocess_sprint(eprocess* p, char* buf)
{
    stringtable* st = p->prms;
    int i;

    sprintf(buf, "%s(", p->name);
    if (st->n > 0)
        strcat(buf, st->se[0]->s);
    for (i = 1; i < st->n; ++i) {
        strcat(buf, ",");
        strcat(buf, st->se[i]->s);
    }
    strcat(buf, ")");
}

/** Creates and initialises process defined by the process entry string.
 * @param e Pointer to ecology
 * @param type Process type
 * @param entry Process entry string
 * @return Pointer to eprocess
 */
eprocess* eprocess_create(ecology* e, EPROCESSTYPE type, char* entry)
{
    char* grouptag = eprocess_type_tags[type];
    eprocess* p = malloc(sizeof(eprocess));
    char seps[] = "()\n";
    char buf[STRBUFSIZE];
    char* token;
    int i;

    p->prms = stringtable_create("process parameters");
    p->prms->unique = 0;        /* duplicating entries allowed */
    p->ecology = e;

    p->type = type;

    if ((token = strtok(entry, seps)) == NULL)
        e->quitfn("ecology: error: \"%s\": %s: empty process name\n", e->processfname, eprocess_type_tags[type]);
    p->name = strdup(clear_string(token));

    if (type == PT_GEN)
        e->quitfn("ecology: error: \"%s\": \"%s\": this process may be entered in \"%s\", \"%s\" or \"%s\" groups only\n", e->processfname, p->name, eprocess_type_tags[PT_WC], eprocess_type_tags[PT_SED], eprocess_type_tags[PT_EPI]);

    if ((token = strtok(NULL, seps)) != NULL) {
        char* prmstr = strdup(clear_string(token));
        char* prmseps = " ,";

        if ((token = strtok(prmstr, prmseps)) != NULL) {
            stringtable_add(p->prms, token, -1);
            while ((token = strtok(NULL, prmseps)) != NULL)
                stringtable_add(p->prms, token, -1);
        }

        free(prmstr);
    }
    /*
     * set the standard procedures for the process 
     */
    for (i = 0; i < NEPROCESSES; ++i) {
        eprocess_entry* pe = &eprocesslist[i];
        char* tag = eprocess_type_tags[pe->type];

        if (strcasecmp(pe->name, p->name) == 0) {
            /*
             * check if the process type is right 
             */
            if (strcasecmp(tag, grouptag) != 0) {
                if (strcasecmp(tag, "generic") != 0)
                    e->quitfn("ecology: error: \"%s\": \"%s\": \"%s\": process group = \"%s\"; %s expected\n", e->processfname, grouptag, p->name, grouptag, tag);
                else if (strcasecmp(grouptag, "water") != 0 && strcasecmp(grouptag, "sediment") != 0)
                    e->quitfn("ecology: error: \"%s\": \"%s\": \"%s\": process group = \"%s\"; \"water\" or \"sediment\" expected\n", e->processfname, grouptag, p->name, grouptag);
            }

            if (pe->nparams >= 0 && pe->nparams != p->prms->n)
                e->quitfn("ecology: error: \"%s\": %s: \"%s\": %d parameters entered; %d parameters expected\n", e->processfname, eprocess_type_tags[type], p->name, p->prms->n, pe->nparams);

            p->delayed = pe->delayed;
            p->init = pe->init;
            p->postinit = pe->postinit;
            p->destroy = pe->destroy;
            p->precalc = pe->precalc;
            p->calc = pe->calc;
            p->postcalc = pe->postcalc;
            break;
        }
    }
    if (i == NEPROCESSES)
        e->quitfn("ecology: error: \"%s\": %s: process \"%s\" not found in the process library\n", e->processfname, eprocess_type_tags[type], p->name);

    /*
     * store the full entry 
     */
    eprocess_sprint(p, buf);
    p->fullname = strdup(buf);

    emstag(LDEBUG, "ecology:eprocess:eprocess_create:","init: %s - %s \n", eprocess_type_tags[type], p->fullname);

    /*
     * initialise the process 
     */
    if (p->init != NULL) {
      eco_write_setup(e, "Calling %s_init\n", p->name);
      p->init(p);
      eco_write_setup(e, "\n");
    }
    return p;
}

/** Process destructor.
 * @param p Process
 */
void eprocess_destroy(eprocess* p)
{
    p->destroy(p);
    free(p->name);
    free(p->fullname);
    stringtable_destroy(p->prms);
    free(p);
}
