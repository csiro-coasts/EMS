### SHOC : Test2 : See Section 11.2 of the User's Manual

A wind of constant stress and direction is blown over a 
homogeneous open ocean of constant depth. The model uses 
cyclic open boundaries reflect the open ocean condition, 
and utilizes constant vertical viscosity and linear bottom 
friction for simplicity. According to Kowalik and Murty (1993, p27) 
the linear resistance coefficient is related to the bottom 
drag coefficient via:

                  r = rho.Cd.|v|			(1)

where Cd is the drag coefficient, v is the bottom current speed 
and rho  is the density. Given a constant eddy viscosity, the Ekman 
depth, DE, is given by (Pond and Pickard, 1983, p108):

                 DE = pi.sqrt(2.Vz/|f|)			(2)

where Vz is the eddy viscosity and f is the Coriolis parameter. 
The surface current speed, Vo, is then given by 
(Pond and Pickard, 1983, eqn 9.10):

                 Vo = (sqrt(2).pi.ts) / (DE.rho.|f|)	(3)

where ts is the wind stress, and the current speed, Vb, at the 
Ekman layer depth is (Pond and Pickard, 1983, p108):

                 Vb = Vo.exp(-pi) ~ 0.04Vo		(4)

Therefore, using a wind stress of 0.01 Nm-2 on and f-plane 
with f = 1e-4, and Vz = 0.0507 m2s-1, the Ekman layer 
depth DE = 100m. Furthermore, using (3) with rho = 1025 kgm-3 
gives Vo = 4.33x10-3 ms-1 and Vb = 1.87x10-4 ms-1. 
Using the bottom velocity in (1) and a nominal drag 
coefficient of Cd = 0.003 gives a resistance coefficient of 
r = 0.00058. 
Linear friction is achieved by setting the parameter UF 
to a large value and Z0 to a low value. 
Using the value of QBFC quoted above, Z0 < 7.7x10-8. Using these 
values with UF = 1.0 provides linear bottom friction with the 
required resistance coefficient. Using the above configuration, 
model results should show an Ekman spiral with velocities rotating 
clockwise with depth, surface current speed ~ 0.0043 ms-1 and 
bottom current speed ~ 0.00019 ms-1. Surface elevation should 
be equal to zero.

Testing is performed for both 'z' and sigma vertical coordinates.
The input parameter file for the 'z' model is : test2.prm
The input parameter file for the sigma model is : test2_s.prm

The 'z' output is stored in the file : out2_z.nc
the sigma output is stored in the file : out2_s.nc

To run type :
```
run_test2
```

To view output execute 'out' in Matlab, e.g.
```
matlab
>out
```
This script provides surface and bottom velocities for each model.
A plot is saved in test2.eps
