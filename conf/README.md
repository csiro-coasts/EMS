# EMS configure and build

These are the configure scripts used to build EMS

This directory is essential to configure and build any EMS application

For details on the configure options type:
```
conf/configure -h
```

*Note that the configure script must always be run from the root EMS directory*

### Supported platforms
EMS is fully tested and maintained on most modern 32 & 64 bit Linux
distribuitions. This is the same architecture type that is available
at most of the high-performance computing centres. With user
contributions, it may be possible to support other platforms.

## Note on NetCDF
If the NetCDF library is not already installed on your system they
may be downloaded from the UNIDATA website
[HERE](https://www.unidata.ucar.edu/software/netcdf/)

If the NetCDF library is not installed in common system areas, such
as */usr* or */usr/local* then the `NETCDF_ROOT` environment variable
must be set to the root dir of the NetCDF installation prior to
running `configure`. This is common in an HPC environment where the
"module" command is used to load/unload the library.

## Most common options
### OpenMP support
OpenMP support allows parallel processing via shared memory
To use:
```
conf/configure --enable-omp
```

### MPI support
The Message Passing Interface allows parallelisation across compute
nodes in a distributed environment such HPC clusters
```
conf/configure --with-mpi=<full_path_to_mpi_library>
```

## EMS build structure
EMS requires certain core libraries, these are:
```
conf
lib
model/lib/grid
```
EMS then needs atleast one of the 2 hydrodynamic models, which are also the "driver" applications. The top-level directories for these are:
##### SHOC
```
model/hd
```
##### COMPAS
```
model/hd-us
ext/jigsaw
```

Depending on the application required, it is not necessary to clone the entire EMS repository but just the directories needed documented above. Alternatively, the entire repository may be cloned but unneeded directories may be deleted prior to running the `configure` script

#### Other EMS libraries
There are several optional EMS libraries available that work with both SHOC and COMPAS
##### Tracerstats
```
model/lib/tracerstats
```
The [tracerstats](../model/lib/tracerstats/README.md) library allows certain statistics to be calculated inline during model runtime.  See the SHOC manual for more details.
##### Waves
```
model/lib/waves
```
The [waves](../model/lib/waves/README.md) library may be used to calculate the wave variables of period,
amplitude, direction, orbital velocity, bottom stress, enhanced bottom drag and tangential
radiation stresses. See the SHOC manual for more details.
##### Sediments
```
model/lib/sediments
model/lib/waves     # dependency
```
The [sediments](../model/lib/sediments/README.md) library allows for sediment layers below the water column and is responsible for the deposition & resuspension of particulate matter. See the SEDIMENT manuals for more details.
##### Ecology
```
model/lib/ecology
model/lib/sediments # dependency
model/lib/waves     # dependency
```
The [ecology](../model/lib/ecology/README.md) library provides nutrient cycling and spectral light propagation through the water/sediment columns plus benthic processes. See the ECOLOGY manual for more details.
##### DA
```
model/lib/da
```
The [da](../model/lib/ecology/README.md) library is an implementation
of the Ensemble Optimal Interpolator (EnOI) scheme. See the SHOC
manual for more details

This also has a dependency on the [GNU Scientific library](https://www.gnu.org/software/gsl/)