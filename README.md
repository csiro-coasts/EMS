# EMS-pre-release

## Features
* Hydrodynamic models for [structured](model/hd) and [unstructured](model/hd-us) grids
* [Sediment transport, resuspension and deposition](model/lib/sediments)
* [Spectrally-resolved semi-empirical optical model](model/lib/ecology)
* [Biogeochemical (nutrients, plankton, carbon chemistry) and benthic (seagrass, corals) processes](model/lib/ecology)

## Background
The CSIRO Coastal Environmental Modelling (CEM) team develops, maintains and uses the EMS software that allows investigation of the physical, sediment and biogeochemical processes in marine environments. This is achieved by a ‘driver’ hydrodynamic code into which are linked various libraries to perform sediment transport and biogeochemistry, all supported by a core library. The ‘driver’ may be any model that manages the tracers required for sediments and biogeochemistry. The sediment and biogeochemical libraries are stand-alone modules that are linked to the driver via an interface, and in principle may be linked to any hydrodynamic code. Currently the ‘drivers’ available are a full hydrodynamic mode, a transport model that uses offline data to advect and diffuse sediment /  biogeochemical variables, and a box model. The hydrodynamic code may further operate in reduced dimensions of 1D, 2D vertically averaged or 2D laterally averaged. A waves and tracer statistic library also exist; the latter allowing various operations to be performed during run-time on any tracers supported by the driver (e.g. means, fluxes, vertical integrals).

Additional software exists to generate the complex orthogonal curvilinear grids that are typically used for case studies. These grids allow variable resolution over the domain, useful for representing areas of interest with high resolution and less critical regions with coarser resolution. The curvilinear grid may also allow a dimensionality to be reduced from 3-D to 2-D within the same grid. This is useful when representing rivers or narrow estuaries, since the cross-river coordinate becomes very small in these areas and therefore becomes the defining grid size for setting the model time-step. Eliminating these small grid cells by making rivers or estuaries 2-D laterally averaged allows larger time-steps, hence a faster model. The curvilinear grids require dedicated software for visualizisation of model output, and the CEM supports several visualisation platforms to archive this. These software packages allow publication quality images and animations to be produced, and allow exploration of the data in 4 dimensions for analysis purposes.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

Things you need to build EMS applications

```
C/C++ compiler
make
NetCDF library
```

### Building
To get started, type the following commands on a Linux command line to build SHOC:

```
conf/configure
make
# The SHOC executable will be model/hd/shoc
# Type the following to get version information
model/hd/shoc -v
```
See [this page](conf/README.md) for full explaination of the build procedure

## Running the tests

There are many tests available for each module. See [model/tests](model/tests) for a complete list.

## Contributing
Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.

## Authors

EMS is wrtten and maintained by the [Coastal Environmental Modelling team](https://research.csiro.au/cem/people/)

See also the list of [contributors](https://github.com/csiro-coasts/EMS-pre-release/contributors) who participated in this project.

## License

This project is licensed under the CSIRO open-source License - see the [LICENSE.md](LICENSE.md) file for details

## Website

Further information about the CSIRO CEM Team may be found by [visiting our website](https://research.csiro.au/cem).
