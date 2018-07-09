### SHOC : Test5 : See Section 11.5 of the User's Manual

A wind applied perpendicular to an infinitely long coastline 
(on the southern boundary in this case) will result in an 
elevation setup against the coast with zero depth averaged 
currents everywhere. The on-shore wind stress drives depth 
averaged flow to the west (east) in the northern (southern) 
hemisphere, and the sea level gradient resulting from setup 
at the coast drives this flow to the east (west). In a perfect 
situation the sea level gradient and wind stress forces balance 
resulting in no flow. Boundary effects and numerical error may 
make one of these forces dominate, leading to non-zero flow in 
the east (west) direction in the northern (southern) hemisphere 
if the sea level pressure gradient dominates, and vice versa if 
wind stress dominates. Assuming a linear model with linear bottom 
friction, the analytical solution for sea level profile is given 
by (Chapman,  1985, eqn 4.13a):

           d^2eta/dy^2 = -(ts/rho.g.D^2)(d D/dy)	(1)

Using a domain with two cyclic cross-shelf open boundaries and 
one offshore boundary with elevation clamped to zero, and the 
linear depth profile used by  Chapman (1985, eqn 4.1), then the 
sea level profile is given by Chapman (1985, eqn 4.14) and shown 
in Table 1. The boundary conditions of eqn. (1) for this domain are:

d D/dy=-Do at y=0 : d eta/dy=ts/rho.g.Do and eta=0 at y=L (2)

where L is the distance to the offshore boundary.

Table 12.5.1 : Sea level profile for onshore wind stress.

Offshore distance (km)	Surface elevation (m)
	5			0.0301
	15			0.0189
	25			0.0137
	35			0.0103
	45			0.0077
	55			0.0057
	65			0.0039
	75			0.0025
	85			0.0012
	95			0.0001

Again, a non-linear 3-D model with quadratic bottom friction is 
expected to give different results, and bottom friction is 
generally required to be increased for solutions to match theory. 
Specifically, adequate solutions were obtained using the CONSTANT 
mixing scheme with background vertical viscosity VZ0 = 0.0507 
(see Test2) and minimum bottom drag coefficient of QBFC = 0.003 
(UF = Z0 =1e-8). Horizontal viscosity of U1VH = U2VH = 800 is 
also required for stability.

Testing is performed for both 'z' and sigma vertical coordinates.
The input parameter file for the 'z' model is : test5.prm
The input parameter file for the sigma model is : test5_s.prm

The 'z' output is stored in the file : out5_z.nc
the sigma output is stored in the file : out5_s.nc

To run type :
```
run_test5
```

To view output execute 'out' in Matlab, e.g.
```
matlab
>out
```
This script provides the velocity components mid-domain for each model.
A plot is saved in test5.eps
