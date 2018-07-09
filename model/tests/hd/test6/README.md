### SHOC : Test6 : See Section 11.6 of the User's Manual

This test simulates a wetting bore propagating along an 
initially dry channel. The model domain represents a 
channel 2km wide and 100km long with uniform (flat) 
bathymetry. A constant velocity of 1 ms-1 is applied at one 
end of the initially dry channel. A bore propagates along 
the channel, with a parabolic shape (surface elevation profile) 
determined by a balance between the quadratic bottom friction 
and surface slope. The length of the bore is related to the 
depth at the inflow via:

	L = g.D^2 / 2.Cd.U^2

Since only the 'z' model is a wetting and drying model, this 
test is performed for the 'z'  model only.
The input parameter file for the 'z' model is : test6.prm

The 'z' output is stored in the file : out6_z.nc

To run type :
```
run_test6
```

