#!/bin/sh
#
# Figures out the subversion revision string and overwrites the header
# file with a macro
#
# Internal CSIRO development only
#
BASE=$1
FNAME=$BASE/lib/include/svn_rev.h
REV=unknown
if which svnversion > /dev/null 2>&1
then
    # -n prevents newline character
    REV=`svnversion -n`
    if [ $REV = exported ]
    then
    REV=unknown
    fi
fi

# Now test to see if the header already exists
doPrint=true
if [ -f $FNAME ]
then
    # Don't update if we don't have a valid revision string
    if [ $REV = unknown ]
    then
 	doPrint=false
    fi
fi

# Write the file, if needed
if ( $doPrint )
then
    \rm -f $FNAME
    printf "#define SVN_REV_STR \"$REV\"\n" > $FNAME
fi

# Always
exit 0
