/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/carbon_chemistry_wc.h
 *  
 *  Description:
 *  Header file for process carbon_chemistry_wc
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: carbon_chemistry_wc.h 5846 2018-06-29 04:14:26Z riz008 $
 *
 */
#if !defined(_carbon_chemistry_wc_H)

void carbon_chemistry_wc_init(eprocess* p);
void carbon_chemistry_wc_postinit(eprocess* p);
void carbon_chemistry_wc_destroy(eprocess* p);
void carbon_chemistry_wc_precalc(eprocess* p, void* pp);
void carbon_chemistry_wc_calc(eprocess* p, void* pp);
void carbon_chemistry_wc_postcalc(eprocess* p, void* pp);
#define _carbon_chemistry_wc_H
#endif


/* Adding a process to the library!
 * 
 * Once you have decided over a unique process name copy this file and 
 * process_template_content.c to new files with your process name, but maintain 
 * the suffix of each file. 
 * In your copy replace all occcuences of REPLACE_ME_WITH_PROCESS_NAME with 
 * your process name. Make sure that only above string is replaced and not what 
 * comes before or after.
 * Read the comments in your 'process name'.c file and modify as needed.
 * Once you are satisfied with the functions you will need to add the process to
 * the list of processes. In the file 'allprocesses.c' is a list of all available processes.
 * It requires an entry with the new process of the type
 * 
 * <process name as string>,<process type>,<process arguments>,<delayed flag>,
 * <REPLACE_ME_WITH_PROCESS_NAME_init>,<REPLACE_ME_WITH_PROCESS_NAME_postinit>,
 * <REPLACE_ME_WITH_PROCESS_NAME_destroy>,<REPLACE_ME_WITH_PROCESS_NAME_precalc>,
 * <REPLACE_ME_WITH_PROCESS_NAME_calc>,<REPLACE_ME_WITH_PROCESS_NAME_postcalc>
 * 
 * on one line enclosed by curly brackets.
 * 
 * <process name> as string - the name as a quoted string with which you like to identify 
 *                          the process in the 'processes.prm' file, generally the same 
 *                          as your process name
 * 
 * <process type>           - the type as defined in include/eprocesses.h:EPROCESSTYPE
 *                          valid values are
 *                            PT_GEN - processes which can be run in water or sediment
 *                            PT_SED - process which can be run only in the sediment
 *                            PT_WC  - process which can be run only in the water column
 *                            PT_EPI - process which can be run only in the epi benthos
 * 
 * 
 * <process arguments>      - the option exists to pass parameters to the function in the 'processes.prm'
 *                          file in the way of <process name >(argument-list).
 *                          The number of arguments as comma separated list 
 *                          in 'argument-list' must be given as integer here or 0
 *                          if none are expected.
 * 
 * 
 * <delayed flag>           - processes are executed in the sequence they are entered in the 
 *                          'processes.prm' file. If a process must be executed after all otheres have finished
 *                          set this to 1, 0 otherwise.
 * 
 * 
 * <REPLACE_ME_WITH_PROCESS_NAME_init>      - your process init function unquoted
 * <REPLACE_ME_WITH_PROCESS_NAME_postinit>  - your process postinit function unquoted
 * <REPLACE_ME_WITH_PROCESS_NAME_destroy>   - your process destroy function unquoted 
 * <REPLACE_ME_WITH_PROCESS_NAME_precalc>   - your process precalc function unquoted
 * <REPLACE_ME_WITH_PROCESS_NAME_calc>      - your process calc function unquoted
 * <REPLACE_ME_WITH_PROCESS_NAME_postcalc>  - your process postcalc function unquoted
 * 
 * 
 * 
 * If you have introduced a new tracer than that should also be added to the 
 * allprocesses.c:diagnflags array in form of
 * 
 * <tracer name>,<diagnostic flag>
 * 
 * this entry is used to verify the correct setting in the host model
 * 
 * <tracer name>      - quoted string of the tracer name
 * <diagnostic flag>  - integer indicating the type of diagnostic this tracer represents
 *                    valid values are 
 *                      0 - not a diagnostic
 *                      1 - flux diagnostic
 *                      2 - value diagnostic
 * 
 * all diagnostic tracers are set to 0.0 before the time step, flux diagnostics 
 * are divided by the time step to provide a rate instead of the total value 
 * the time step.
 * Please not that the <TARCERNAME>.diagn entry must reflect the entry in the 
 * allprocesses.c:diagnflags entry.
 * 
 * 
 * Furthermore add the your copy of this file as an #include "REPLACE_ME_WITH_PROCESS_NAME.h"
 * to the allprocesses.c file
 *  
 * 
 * Finally add an entry <location of process file>/REPLACE_ME_WITH_PROCESS_NAME.o 
 * to the make file you are using. The location is indicated in the example make file.
 * 
 * 
 * Last but not least don't forget to add your process to the 'processes.prm' you
 * are using.
 * 
 * 
 * Have fun 
 *    Uwe Rosebrock
 */



