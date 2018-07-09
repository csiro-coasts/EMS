### SHOC : Test3 : See Section 11.3 of the User's Manual

This test examines the set-up due to a steady wind applied to a 
1-layer (depth-averaged) model domain. If a constant wind is 
applied to a homogeneous closed basin of constant depth then 
depth averaged velocities are equal to zero in the steady state. 
For a linear model and constant wind stress in the x direction 
the surface slope should balance the applied wind stress, and 
the equations of motion reduce to an expression for the sea level 
gradient in the x direction:

                d eta/dx = -ts / (rho.g.D)

where D is the water depth. If a wind stress of 1 Nm-2 is applied 
to a homogeneous ocean of temperature 20oC and salinity 35 psu so 
rho = 1024.76 kgm-3, then the slope in a 10m deep basin is equal 
to 9.958x10-6, and the depth averaged velocities should be near zero.

Testing is performed for both 'z' and sigma vertical coordinates.
The input parameter file for the 'z' model is : test3.prm
The input parameter file for the sigma model is : test3_s.prm

The 'z' output is stored in the file : out3_z.nc
the sigma output is stored in the file : out3_s.nc

To run type :
```
run_test3
```

To view output execute 'out' in Matlab, e.g.
```
matlab
>out
```
This script provides surface slope for each model.
A plot is saved in test3.eps
