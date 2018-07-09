### Tracerstats library test: basic

This test consists of a single water and 2 sediment cells and makes sure that
the statistics are basically doing the right thing.

Salinity is setup to be reset every timestep (in the water column)
from the values in the timeseries file.  Sediment tracers remain
constant.

To run the test type `./run_test`

type `out` in Matlab to see the results
