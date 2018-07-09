/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/string_utils.h
 *
 *  \brief Header for string utilities
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: string_utils.h 5834 2018-06-27 00:55:37Z riz008 $
 */


#ifndef	_STRING_UTILS_H
#define _STRING_UTILS_H


#define WILDCARD '*'
 
/* For Window 32 if required */
#if defined(_WIN32)  &&  !defined(__MINGW32__)
int strcasecmp(const char *s1, const char *s2);
int strncasecmp(const char *s1, const char *s2, int n);
#endif

int is_blank(int c);

/* Prototypes */
char *contains_token(const char *haystack, const char *needle);
void strip(char *str, const char *seps);
void stripend(char *str);
int parseline(char *line, char **str, int max);
int find_token(const char* haystack, const char* needle, 
	       char *buf, const char term);

/* test if the string contains the common descriptios for true or false
 * including 1 and 0
 */
int is_true(const char *tag);

/* test if the string is not null, only space and length greater 0 */
int is_valid(const char* s);

/* test if the string contains relevant characters */
int only_space(const char* s);
/**
 * extract the int ranges from a mixed string of individual ints comma separated
 * and ranges separated by an underscore
 *
 * @param arg - the string containing the content
 * @param nv - reference to an int which contain the number of entries in the result array
 * @return  an one dimensional array of ints extracted from the string argument
 */
int* split_int_list(char* arg, int* nv );

/**
 * Establish if char c is one of the characters in the string elements
 * this isthe same as calling 'contains_char(c,elements,0,strlen(elements))'
 * @param c - the character to test on
 * @param elements the string to interogate
 * @param nelm the length of the string elements
 * @return 1 if c is an element of elements, 0 otherwise
 */
int contains_char(char c, char* elements);


/**
 * Establish if char c is one of the characters in the string elements 
 * between the start index and the end index
 * @param c - the character to test on
 * @param elements the string to interogate
 * @param start the first of the string elements tobe tested
 * @param end the last element+1 of the string elements tobe tested
 * @return 1 if c is an element of elements, 0 otherwise
 */
int containsn_char(char c, char* elements,int start, int end);



/** Parse a string and split it into fields as delimited
  * by elements of del repescting quotes.
  * Whitespaces can be delimiters, whitespace at the
  * start of the string is removed.
  *
  * @param line pointer to storage for line read
  * @param str returned array of string fields.
  * @param del the string which chars are serving as delimiter .
  *
  * @return number of fields read.
  */
int split(char *line, char **str, char* del);

/**
 * Match if any of the pattern of the list of patterns in 'patternlist' are matched in c.
 * The function provides a wrapper around the underlying fnmatch of the standart C library
 * @param c  - the string to interogate
 * @param pattern - the pattern list
 * @return 1 if a match can be found, 0 otherwise 
 *   
 */
int contains_pattern(const char* c, const char* patternlist);


/** Finds if the substring `needle' is equal to the end of a string `haystack'.
  *
  * @param needle substring
  * @param haystack string
  * @return 1 if above is true, 0 otherwise
  */
int endswith(const char *haystack, const char *needle);


/**
 * Convert string (i.e. char *) to double
 * @param token - the string
 * @param value - pointer to the output double value
 * @return 1 on success, 0 with *value = NaN for empty or null string
 */
int str2double(char* token, double* value);

/* .. and to int */
int str2int(char* token, int* value);


/** Finds if the substring `needle' is equal to the start of a string `haystack'.
  *
  * @param needle substring
  * @param haystack string
  * @return 1 if above is true, 0 otherwise
  */
int startswith(const char *haystack, const char *needle);


/** Finds if the substring `needle' up to 'needle_length' is equal to the start of a string `haystack'.
  *
  * @param needle substring
  * @param haystack string
  * @param needle_length length to interrogate needle for
  * @return 1 if above is true, 0 otherwise
  */
int startsnwith(const char *haystack, const char *needle, int needle_length);


/**
 * Expand a list of space delimited files containing the possible
 * wildcard '*' to a list of filenames
 * 
 * @param line - the line to parse
 * @param str - an unallocated 2d array of chars
 * @param max - the maximum number of entries
 * @return the number of entries found     
 */
char ** expand_files(char *line, int* nelem, int max);
/*int expand_files(char *line, char **str, int max); */
#endif                          /* _STRING_UTILS_H */
