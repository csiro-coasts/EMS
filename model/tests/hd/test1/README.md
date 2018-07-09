### SHOC : Test1 : See Section 11.1 of the User's Manual

Test 1 is an extremely simple, null case, test, where no forcing is 
applied to a closed model domain. The purpose is to demonstrate that 
no model variables deviate from their initial values. A rectangular 
grid is used with a horizontal grid of 5 by 10 cells, and 5 layers in 
the vertical having 1m vertical spacing. The bathymetry varies and the 
water is initially vertically stratified. Vertical diffusion of salt 
and heat is turned off in this case to avoid diffusive changes. Each 
run consists of a 1000 second integration with no externally applied 
forcing. The initial variable values (zero elevation, zero velocity, 
salinity and temperature) should remain unchanged for the duration of 
the integration.

Testing is performed for both 'z' and sigma vertical coordinates.
The input parameter file for the 'z' model is : test1.prm
The input parameter file for the sigma model is : test1_s.prm

The 'z' output is stored in the file : out1_z.nc
the sigma output is stored in the file : out1_s.nc

To run type :
```
run_test1
```
To view output execute 'out' in Matlab, e.g.

```
matlab
> out
```

A plot is saved in test1.eps