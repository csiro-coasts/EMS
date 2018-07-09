### SHOC : Test4 : See Section 11.4 of the User's Manual

An analytical solution exists for a linear model of constant 
wind stress applied in a longshore direction along an infinitely 
long coast. Assuming cross-shelf transport and alongshore sea 
level gradient are small, then along shelf transport, U, is 
given by (Chapman,  1985, eqn 4.5):

            U = (D.ts / rho.r)(1 - exp(-t.r/D))		(1)

with a steady state velocity (t=infinity) given by:

            U = (D.ts / rho.r)				(2)

Note U is the transport, hence velocity u = U/D.
The sea surface slope is given by (Chapman,  1985, eqn 4.6):

d eta/dy = -f.U/g.D = (f.ts / rho.r.g)(1 - exp(-t.r/D))	(3)

with steady state sea level given by:

            eta(y) = (f.ts / rho.r.g)(L/2 - y)		(4)

where L is the width of the channel and it is assumed eta = 0  
at y = -L/2.
Cyclic open boundaries are used to represent an infinite coastline. 
Using a wind stress of 0.1 Nm-2 in a channel 500km wide 
(the dimensions of this domain are the same as the test domain 
used by Palma and Matano (1998) except the Southern hemisphere 
is considered) with linear resistance coefficient r = 0.0005, 
Coriolis = -1.028e-4 and rho = 1024.76 kgm-3, the along-shore 
depth averaged velocity is 0.195 ms-1, the cross-shore depth 
averaged velocity is zero and elevation at the coast is -0.49 m 
(slope of 2.05x10-6. Note the first elevation cell center is 
found at y = 10km). This result assumes a linear depth averaged 
model is used, and a non-linear 3-D model with quadratic bottom 
friction is expected to give different results.

Testing is performed for both 'z' and sigma vertical coordinates.
The input parameter file for the 'z' model is : test4.prm
The input parameter file for the sigma model is : test4_s.prm

The 'z' output is stored in the file : out4_z.nc
the sigma output is stored in the file : out4_s.nc

To run type :
```
run_test4
```

To view output execute 'out' in Matlab, e.g.
```
matlab
>out
```
This script provides elevation and velocity components at
the coast for each model.
A plot is saved in test4.eps
