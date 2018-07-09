#!/bin/sh
#
# Test for hdf5 threadsafety
#
# This seems a bit convoluted but I have not figured out how to do
# this with a runtime call. Its further complicated by the fact that
# the HDF5 library is implicitly linked in via the NETCDF library
#
# Farhan Rizwi 2014-11-11
#
# NETCDF should be available by a previous check
NC=`which ncdump`
H1=`ldd $NC | grep libhdf5.so. | cut -d '>' -f 2 | cut -d ' ' -f 2`
# Bail out if hdf5 not found
if [ -z "$H1" ] ; then
    exit 0
fi
H2=`dirname $H1`
HF="$H2/libhdf5.settings"
# Bail out if settings file not found
if [ ! -e ${HF} ] ; then
    exit 0
fi
# Grep for feature
H3=`grep -i threadsafety $HF`
# Final grep
echo $H3 | grep yes > /dev/null
if [ $? -eq 0 ] ; then
    # All good
    echo yes
fi

exit 0
