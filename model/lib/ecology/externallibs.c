/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: lib/ecology/externallibs.c
 *  
 *  Description:
 *  Provides an interface to external libraries. Present only for
 *  backward compatibility with SJWlib.
 *  
 *  Shall be eventually removed when SJWlib is retired.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: externallibs.c 5846 2018-06-29 04:14:26Z riz008 $
 *
 */
#include <stdio.h>
#include "utils.h"
#include "externallibs.h"

/*
 * prm_seterrorfn()
 */

prm_seterrfn_fn prm_seterrorfn = prm_set_errfn;


/*
 * prm_geterrorfn()
 */

extern error_fn keyprm_errfn;
error_fn prm_geterrorfn()
{
    return keyprm_errfn;
}


/*
 * prm_readstring()
 */
prm_readstring_fn prm_readstring = prm_read_char;


/*
 * prm_readdouble()
 */
prm_readdouble_fn prm_readdouble = prm_read_double;


/*
 * prm_readint()
 */
prm_readint_fn prm_readint = prm_read_int;


/*
 * prm_getkey()
 */
prm_getkey_fn prm_getkey = prm_get_key;

/*
 * prm_parseline
 */
prm_parseline_fn prm_parseline = parseline;



