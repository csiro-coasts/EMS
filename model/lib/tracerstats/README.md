# EMS Tracerstats library

This library allows for computing statistics on tracer values during
model simulation at the model timesteps. This can greatly reduce the
storage requirements needed to compute the statistics.

#### Current list of supported stats:
* fluxe1 : Flux of 3D tracer in the e1 direction
* fluxe2 : Flux of 3D tracer in the e2 direction
* fluxw  : Flux of 3D tracer in the vertical
* meanfluxe1 : Mean flux of 3D tracer in the e1 direction
* meanfluxe2 : Mean flux of 3D tracer in the e2 direction
* meanfluxw  : Mean flux of 3D tracer in the vertical
* mean : Mean of 2D or 3D tracer
* variance : Variance of 2D or 3D tracer
* stdev : Standard deviation of 2D or 3D tracer
* corr : Correlation coefficient of two 2D or 3D tracers
* cov  : Covariance of two 2D or 3D tracers
* sum : Sum of 2D or 3D tracers
* diff : Difference of two 2D or 3D tracers
* max : Maximum tracer value over a given time period
* min : Minimum tracer value over a given time period
* vmax : Maximum water column value of a 3D tracer
* vmin : Minimum water column value of a 3D tracer
* vint : Vertical integral of 3D tracer
* vmean : Vertical mean of 3D tracer
* vdiff : The ratio of a series of layers of a 3d tracer
* copy : Copy of a tracer
* sectionflux : Integrate the flux over a defined area and time
* rmse : Compute the RMS error between two tracers

##### See section 13 of the [SHOC User guide](https://research.csiro.au/cem/?ddownload=124) for full details
      