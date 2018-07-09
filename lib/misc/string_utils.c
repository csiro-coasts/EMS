/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/misc/string_utils.c
 *
 *  \brief String utilities
 *
 *  String parsing utilities
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: string_utils.c 5831 2018-06-26 23:48:06Z riz008 $
 */


#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <fnmatch.h>
#include <glob.h>
#include <unistd.h>
#include "ems.h"

/*  #define DOTEST 1 */


/** Parse a line and split it into fields as delimited
  * by whitespace respecting quotes.
  *
  * @param line pointer to storage for line read
  * @param str returned array of string fields.
  * @param max mamimum number of fields.
  *
  * @return number of fields read.
  */
int parseline(char *line, char **str, int max)
{
  int quote;
  int n;

  if (line == NULL)
    return(-1);

   /* Skip whitespace */
  while (*line && isspace((int)*line))
    line++;

  /* Loop to find string separators */
  for (n = 0; *line && n < max; n++) {
    /* Check for start of quote */
    if (*line == '\'' || *line == '"')
      quote = *line++;
    else
      quote = 0;
    /* Store string start */
    str[n] = line;
    /* Loop to end of string */
    while (*line &&
           ((quote && (*line != quote)) || (!quote && !isspace((int)*line)
            )))
      line++;
    if (*line)
      *line++ = 0;
    while (*line && isspace((int)*line))
      line++;
  }
  if(n >= max &&  *line)
  	emstag(LFATAL,"lib:util:string_utils","string outside allowable number of elements > %u, string = '%s'\n",max, line);
  return (*line ? n + 1 : n);
}


/** Finds if the substring `needle' is a token of a string `haystack'.
  * Possible delimiters are: \<space\>.
  *
  * @param needle substring
  * @param haystack string
  * @return pointer to he beginning of substring, or NULL if the substring
  * is not found or is not a separate token.
  */
char *contains_token(const char *haystack, const char *needle)
{
  char *position;
  int length = strlen(needle);

  while ((position = strstr(haystack, needle)) != NULL) {
    if (position == haystack || position[-1] == ' ')
      if (position[length] == ' ' || position[length] == 0)
        break;
    haystack = position;
    haystack++;
  }

  return position;
}




/** Finds if the substring `needle' is equal to the start of a string `haystack'.
  * (same as startsnwith(haystack, needle, strlen(needle)))
  * @param needle substring
  * @param haystack string
  * @return 1 if above is true, 0 otherwise
  */
int startswith(const char *haystack, const char *needle)
{
	return startsnwith(haystack, needle, strlen(needle));
}


/** Finds if the substring `needle' up to 'needle_length' is equal to the start of a string `haystack'.
  *
  * @param needle substring
  * @param haystack string
  * @param needle_length length to interrogate needle for
  * @return 1 if above is true, 0 otherwise
  */
int startsnwith(const char *haystack, const char *needle, int needle_length)
{
  char *position;
  int i,l1= strlen(haystack);
  /*,l2 = strlen(needle); */
	l1= strlen(haystack);
	if(needle_length > l1)
		return 0;
	
	for(i = 0;i < needle_length; i++ )
	{
		if(haystack[i] != needle[i])
			return 0;
	}
	return 1;
}




/** Finds if the substring `needle' is equal to the end of a string `haystack'.
  *
  * @param needle substring
  * @param haystack string
  * @return 1 if above is true, 0 otherwise
  */
int endswith(const char *haystack, const char *needle)
{
  char *position;
  int i,j,l1,l2 = strlen(needle);
	l1= strlen(haystack);
	if(l2 > l1)
		return 0;
	
	for(i = l1-1, j = l2-1;i >= 0 && j >= 0 ; i--,j-- )
	{
		if(haystack[i] != needle[j])
			return 0;
	}
	return 1;
}

/** Finds token and returns string up to but not including the
 *  'term' character
 * @param haystack string to search
 * @param needle   string to find in haystack
 * @param out      output string buffer
 * @param term     terminating character
 *
 */
int find_token(const char* haystack, const char* needle, 
	       char *out, const char term)
{
  char *str, *src, *dst;
  char buf[MAXSTRLEN];
  char buf1[MAXSTRLEN];
  out[0] = '\0';

  /* See if the token exists at all */
  if ( (str = strstr(haystack, needle)) == NULL )
    return(0);

  /* Fast forward past the token */
  sprintf(buf, "%s%s%c", needle, "%s", term);
  if (!sscanf(str, buf, buf1))
    return(0);
  
  /* Copy in characters up to term */
  src = buf1;
  dst = out;
  while (*src != term)
    *(dst++) = *(src++);

  /* Terminate string */
  *dst = '\0';

  return(1);
}

/** Strips a string from specified characters.
 * @param str input / output
 * @param seps array of characters to be stripped from `str'
 */
void strip(char *str, const char *seps)
{
  char *token;
  char buf[MAXLINELEN];
  char buf2[MAXLINELEN];

  if (str == NULL)
    return;
  strcpy(buf, str);

  if ((token = strtok(buf, seps)) == NULL)
    return;
  strcpy(buf2, token);

  while ((token = strtok(NULL, seps)) != NULL)
    strcpy(buf2 + strlen(buf2), token);
  strcpy(str, buf2);
}

/** Removes all characters after '.'  in a string of specified characters.
 * @param str input / output
 */
void stripend(char *str)
{
  char *token;
  char buf[MAXLINELEN];
  int n;

  if (str == NULL)
    return;

  for (n = 0; n < strlen(str); n++) {
    if (str[n] == '.') 
      break;
    else
      buf[n] = str[n];
  }
  buf[n] = '\0';
  strcpy(str, buf);
}

/* Check whether the parameter tag, matchs one of the
 * valid boolean definitions.
 */
int is_true(const char *tag)
{
  char ltag[MAXLINELEN];
  int i;
  int n = strlen(tag);

  for (i = 0; i < n; ++i)
    ltag[i] = tolower(tag[i]);
  ltag[n] = 0;

  if (strcmp(ltag, "true") == 0)
    return 1;
  else if (strcmp(ltag, "false") == 0)
    return 0;
  else if (strcmp(ltag, "yes") == 0)
    return 1;
  else if (strcmp(ltag, "no") == 0)
    return 0;
  else if (strcmp(ltag, "1") == 0)
    return 1;
  else if (strcmp(ltag, "0") == 0)
    return 0;
  else if (strcmp(ltag, "include") == 0)
    return 1;
  else if (strcmp(ltag, "exclude") == 0)
    return 0;
  else if (strcmp(ltag, "positive") == 0)
    return 1;
  else if (strcmp(ltag, "negative") == 0)
    return 0;

  quit("The boolean tag '%s' is unknown.\n", ltag);

  return 0;
}


/*UR why aren't we using isspace?
 */
int is_blank(int c)
{
  return isspace(c);
/*  return ((c == ' ') || (c == '\t')); */
}


#if defined(_WIN32)  &&  !defined(__MINGW32__)

int strcasecmp(const char *s1, const char *s2)
{
  unsigned int i, n1, n2;

  n1 = strlen(s1);
  n2 = strlen(s2);

  if (n1 != n2)
    return (-1);

  for (i = 0; (s1[i] != '\0'); i++)
    if (tolower(s1[i]) != tolower(s2[i]))
      return (-1);

  return (0);
}


int strncasecmp(const char *s1, const char *s2, int n)
{

  unsigned int i;

  if (n > min(strlen(s1), strlen(s2)))
    return (-1);

  for (i = 0; i <= (n - 1); i++)
    if (tolower(s1[i]) != tolower(s2[i]))
      return (-1);

  return (0);
}

/** Strip newline character
 * @param str - the string to process
 */
void stripnl(char *str) {
  while(strlen(str) && ( (str[strlen(str) - 1] == 13) ||
       ( str[strlen(str) - 1] == 10 ))) {
    str[strlen(str) - 1] = 0;
  }
}

#endif

/** Test if the string contains charactrs other
 * than space, tab or newline
 * @param s - the string to interogate
 * @return 0 if there are other characters than space or tab, 1 otherwise
 */
int only_space(const char* s)
{
  int i,len;
  if(s ==NULL)
   return 1;
  len = strlen (s);
  for(i = 0; i < len ;i++ )
  {
    if(s[i] != ' ' ||
        s[i] != '\t' ||
        s[i] != 10 ||
        s[i] != 13)
      return 0;
  }
  return 1;
}


/** Test if this string or char is a valid string
 * which mustbe not space, not null, of length greater than 0 and not only whitespace
 * @param s - the string to test
 * @return 1 if not space, not null, of length greater than 0 and not only whitespace, 0 otherwise
 */
int is_valid(const char* s)
{
  return (s != NULL && strlen(s) > 0 && !only_space(s));
}


/**
 * extract the int ranges from a mixed string of individual ints comma separated
 * and ranges separated by an underscore
 *
 * @param arg - the string containing the content
 * @param nv - reference to an int which will contain the number of entries in the result array
 * @return  an one dimensional array of ints extracted from the string argument
 */
int* split_int_list(char* arg, int* nv )
{
  char* tok;
  char* tok2;
  char tbuf[MAXSTRLEN];
  char ttbuf[MAXSTRLEN][MAXSTRLEN];
  int lng = strlen(arg);
  int sz=0;
  int nt;
  int idash=0;
  int iind = 0;
  int d,n,nn,nnn;
  int** dashranges;
  int* indvals;
  int* vals;

  /*count how many we have */
  for(n=0;n<lng;n++)
  {
    if(arg[n] == '_')
      idash++;
    else if(arg[n] == ',')
      iind++;
  }
  emstag(LTRACE,"lib:misc:string_utils:split_int_list","Found individual vals: %d - ranges: %d from %s",iind,idash,arg);

  /* do we have a single value */
  if(iind == 0 && idash == 0)
  {
    vals = i_alloc_1d(1);
    vals[0] = atoi(arg);
    *nv = 1;
    return vals;
  }else if(iind == 0)/* only one range */
  {
    tok = strchr(arg,'_');
    n = atoi(strncpy(tbuf,arg,tok-arg));
    nn = atoi(tok+1);
    if(n==nn) /* silly*/
    {
      vals = i_alloc_1d(1);
      vals[0] = nn;
      *nv = 1;
      return vals;
    }
    if(nn > n)
    {
        vals = i_alloc_1d(nn-n);
        for(nnn=0;n <= nn;n++,nnn++)
        vals[nnn] = n;
    }else
    {
      vals = i_alloc_1d(n-nn);
        for(nnn=0;nn <= n;nn++,nnn++)
        vals[nnn] = nn;
    }
    *nv = nnn;
    return vals;
  }else if(idash == 0)/* no range */
  {
      vals = i_alloc_1d(iind+1);
      tok = strtok(arg,",");
      n=0;
      while(tok!= NULL)
      {
        vals[n] = atoi(tok);
        tok = strtok(NULL,",");
        n++;
      }
      *nv = n;
      return vals;
  }
  /*now we need to do some serious parsing */

  sz = iind + 1 - idash;/* shouldn't be here if sz < 0 */
  dashranges = i_alloc_2d(idash,2);
  indvals = i_alloc_1d(sz);

  n = 0;
  tok2 = strtok(arg,",");
  strcpy(ttbuf[n],tok2);
  n++;

  while(tok2 != NULL)
  {
    tok2 = strtok(NULL,",");
    if(tok2 != NULL)
    {
      strcpy(ttbuf[n],tok2);
      n++;
    }else
      break;
  }

  iind = sz;
  d = 0;
  nnn = 0;
  for(nn=0 ;nn<n;nn++)
  {
    tok = strchr(ttbuf[nn],'_');
    if(tok == NULL)/* single value */
    {
      indvals[nnn]= atoi(ttbuf[nn]);
      nnn++;
    }else
    {
      dashranges[d][1] = atoi(tok+1);
      nt = strlen(tok) - 1;
      ttbuf[nn][nt]='\0';
      dashranges[d][0] = atoi(ttbuf[nn]);
      if(dashranges[d][0] > dashranges[d][1])
      {
        nt = dashranges[d][0];
        dashranges[d][0] = dashranges[d][1];
        dashranges[d][1] = nt;
      }
      sz = sz + abs(dashranges[d][1]-dashranges[d][0])+1;
      d++;
    }
  }
  emstag(LTRACE,"lib:misc:string_utils:split_int_list","Found total vals: %d ",sz);

  vals = i_alloc_1d(sz);
  nnn=0;
  for(n =0;n < d;n++)
  {
    for(; dashranges[n][0] <= dashranges[n][1] ;)
    {
       vals[nnn] =  dashranges[n][0];
       dashranges[n][0]++;
       nnn++;
    }
  }

  for(n =0; n < iind && nnn < sz ;n++)
  {
    vals[nnn] =   indvals[n];
    nnn++;
  }

  tbuf[0] = '\0';
  if(is_log_enabled(LTRACE))
  {
    for(n= 0; n < sz ;n++)
    {
        sprintf(tbuf,"%s %d",tbuf,vals[n]);
    }
    emstag(LTRACE,"lib:misc:string_utils:split_int_list","Values %s ",tbuf);
  }
  free_2d(dashranges);
  free_1d(indvals);

  *nv = sz;
  return vals;
}
/* END split_ints ------------------------------------------------*/


/*
 * find string tname in the array of strings names
 * @param names - 2d char to interogate
 * @param ntr -length of names
 * @param tname - the string to match
 * @return the index of name in names if the string name exists in names, -1 otherwise
 */
int contains_name(char** names, int ntr, char* name)
{
  int i;
  for(i = 0; i < ntr; i++)
  {
    if(strcmp(names[i],name) == 0)
      return i;
  }
  return -1;
}


/**
 * Establish if char c is one of the characters in the string elements
 * @param c the character to test on
 * @param elements the string to interogate
 * @return 1 if c is an element of elements, 0 otherwise
 */
int contains_char(char c, char* elements)
{
	return containsn_char(c,elements,0,strlen(elements));
	/*return (strchr(elements,c) != NULL); */
}

/**
 * Establish if char c is one of the characters in the string elements
 * within range of indicies
 * @param c the character to test on
 * @param elements - the string to interogate
 * @param start start index to search from
 * @param end end index to search to
 * @return 1 if c is an element of elements, 0 otherwise
 */
int containsn_char(char c, char* elements, int start, int end)
{
  int i;
  if(elements == NULL )
    return 0;
    
  for(i = start;i < end;i++)
  {
    if(c == elements[i])
      return 1;
  }
  return 0;
}


/**
 * Convert string (i.e. char *) to double
 * @param token - the string
 * @param value - pointer to the output double value
 * @return 1 on success, 0 with *value = NaN for empty or null string
 */
int str2double(char* token, double* value)
{
    char* end = NULL;

    if (token == NULL) {
        *value = NaN;
        return 0;
    }

    *value = strtod(token, &end);

    if (end == token) {
        *value = NaN;
        return 0;
    }

    return 1;
}

/**
 * Convert string (i.e. char *) to int
 * @param token - the string
 * @param value - pointer to the output integer value
 * @return 1 on success, 0 with *value = NaN for empty or null string
 */
int str2int(char* token, int* value)
{
    char* end = NULL;

    if (token == NULL) {
        *value = INT_MIN;
        return 0;
    }

    *value = strtol(token, &end, 10);

    if (end == token) {
        *value = INT_MIN;
        return 0;
    }

    return 1;
}

/**
 * Match if any of the pattern of the list of patterns in 'patternlist' are matched in c.
 * The function provides a wrapper around the underlying fnmatch of the standart C library
 * @param c  the string to interogate
 * @param patternlist the pattern list
 * @return 1 if a match can be found, 0 otherwise 
 *   
 */
int contains_pattern(const char* c, const char* patternlist)
{
	/* char[] s = {"+",c}; */ 
	
	int val = fnmatch(patternlist,c,FNM_PERIOD);/* FNM_NOESCAPE);/ * (casesensitive?FNM_EXTMATCH:FNM_EXTMATCH | FNM_CASEFOLD)); */
	if(val == 0)
		return 1;
	if(val != FNM_NOMATCH)
	  emstag(LWARN,"An error occured testing for string match mask '%s' - to '%s', return val %d",(char*)patternlist,c,val);
		
	return 0;
}


/* int expand_files(char *line, char **str, int max) */
char ** expand_files(char *line, int* nelem, int max)
{
  int nfiles = 0;
  char files[max][MAXSTRLEN]; 
  char** str = c_alloc_2d(MAXSTRLEN,max);
  
  int ff,j,i, nentries = parseline(line,((char**)files),max);
  
  
  glob_t* matches = malloc(sizeof(glob_t));
  
  for(i = 0 ; i < nentries;i++)
    {
      if(contains_char(WILDCARD,((char**)files)[i]))
	{
	  ff = glob(((char**)files)[i],GLOB_NOSORT,NULL,matches);
	  if(ff == 0)
	    {
	      for(j =0;j< matches->gl_pathc;j++)
		{
		  strcpy(str[nfiles],matches->gl_pathv[j]);	
		  nfiles++;
		}
	    }
	}else
	  {
	    strcpy(str[nfiles],((char**)files)[i]);	
	    nfiles++;
	  }
    }
  
  
  
  globfree(matches);
  *nelem = nfiles;
  
  return str;
}




/** Parse a string and split it into fields as delimited
  * by elements of del representing quotes.
  * Whitespaces can be delimiters, whitespace at the
  * start of the string are removed.
  *
  * @param line - pointer to storage for line read
  * @param str - returned array of string fields.
  * @param del - the string which chars are serving as delimiter .
  *
  * @return number of fields read.
  */
int split(char *line, char **str, char* del)
{
  int quote;
  int n;
  int ndel= strlen(del);

  /* Skip whitespace at the start*/
  while (*line && isspace((int)*line))
    line++;

  /* Loop to find string separators */
  for (n = 0; *line ; n++) {
    /* Check for start of quote */
    if (*line == '\'' || *line == '"')
      quote = *line++;
    else
      quote = 0;
    /* Store string start */
    str[n] = line;
    /* Loop to end of string */
    while (*line &&
           ((quote && (*line != quote)) || (!quote && !startsnwith(line,del,ndel))))
      line++;
    if (*line)
      *line++ = 0;
    while (*line && containsn_char((int)*line,del,0,ndel))
      line++;
  }
  return (*line ? n + 1 : n);
}

/**
 * Trim the argument string of all surrounding spaces
 * @param s - the string to strip the spaces
 * @return a the string excluding the surrounding spaces 
 * 					or a valid string of length if param s is NULL 0
 */
char* trim(char* s)
{
  int len;
  if(s == NULL)
    return "\000";
  while (isspace((int)s[0]))
      s++;
  while ((len = strlen(s)) > 0 && isspace((int)s[len - 1]))
      s[len - 1] = 0;
  return s;
    /*UR old, too complicated
  int n,l=0,i1=0,i2 = 0;
  char* buf;
  if(s == NULL)
    return "\000";
  l = strlen(s);
  i2 = l -1;
  for(n= 0 ;n<l;n++)
  {
    if(isspace(s[n]))
      i1++;
    else
      break;
  }

  for(n = i2 ;n>= 0;n--)
  {
    if(isspace(s[n]))
      i2--;
    else
      break;
  }

  if(i1 >= i2)
    return "\000";

  buf = malloc((i2-i1+1) * sizeof(char));

  for(n = 0 ;n <= i2;n++)
  {
    buf[n] = s[i1 + n];
  }

  return buf;
  */
}





#ifdef DOTEST

static int test_expand_files()
{
  /* should result in 3,1,1,1 files*/
  char line[] = {"./misc/*_utils.c  ./misc/*y.c ./math/*y.c ./misc/string_utils.c"};
  int nelem;
  // As BJR notes, the function signature has changed. -FR
  (void)expand_files(((char*)line), &nelem, 20);
  if(nelem == 6)
    return 1;
  
  return 0;
}



static int test_parseline()
{
	return 1;
}



static int test_contains_token()
{
	return 1;
}

int main(int argc, char *argv[])
{
	int t,ts= 0,s = 0,i;
	for(i = 1;i< argc;i++)
	{
		if(strcmp(argv[i],"expand_files") == 0)
		{
			ts++;
			t= test_expand_files();
			s += t;
			fprintf(stderr, "expand files test: %s \n",(t?"Success":"FAILED"));
		}else if(strcmp(argv[i],"parseline")== 0)
		{
			ts++;
			t= test_parseline();
			s += t;
			fprintf(stderr, "parseline test: %s \n",(t?"Success":"FAILED"));
		}else if(strcmp(argv[i],"contains_token")== 0)
		{
			ts++;
			t= test_contains_token();
			s += t;
			fprintf(stderr, "contains token test: %s \n",(t?"Success":"FAILED"));
		}else
			fprintf(stderr, "Option not recognized: %s \n",argv[i]);
	}
	
	fprintf(stderr, "Runtests, %u of %u successful \n\n",s,ts);
	
	
	
	exit(0);
}

#endif


// EOF
