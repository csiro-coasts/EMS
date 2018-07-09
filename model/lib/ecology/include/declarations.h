/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: lib/ecology/include/declarations.c
 *  
 *  Description:
 *  C stuff; used for handling structures in the ecology code
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: declarations.h 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#if !defined(_DECLARATIONS_H)

#if !defined(_ECOLOGY_H)
struct ecology;
typedef struct ecology ecology;
#endif

#if !defined(_EPROCESS_H)
struct eprocess;
typedef struct eprocess eprocess;
#endif

#if !defined(STRINGTABLE_STRUCT)
struct stringtable;
typedef struct stringtable stringtable;

#define STRINGTABLE_STRUCT
#endif

#if !defined(_PARAMETER_INFO_H)
struct parameter_info;
typedef struct parameter_info parameter_info;
#endif

#if !defined(_COLUMN_H)
struct column;
typedef struct column column;
#endif

#if !defined(_COLUMN_H)
struct intargs;
typedef struct intargs intargs;
#endif

#if !defined(_CELL_H)
struct cell;
typedef struct cell cell;
#endif

#define _DECLARATIONS_H
#endif
