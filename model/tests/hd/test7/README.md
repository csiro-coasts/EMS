### SHOC : Test7 : See Section 12.7 of the User's Manual

Wind stress possessing curl applied to a closed basin with 
a gradient of f/D results in the formation of a gyre due to 
conservation of potential vorticity which is biased to the 
east if f/D < 0 and biased to the west if f/D > 0 
(e.g. Herzfeld and Tomczak, 1999). The gradient of f/D may 
result from a gradient of f (beta effect) or a change in 
topography. This test consists of a closed basin in the 
southern hemisphere with constant depth in the east - west 
direction, 50m depth at the southern coast and 100m depth 
at the northern coast. Wind stress in the e1 direction is 
applied, with 0.1 Nm-2 at the southern boundary and 0 at the 
northern boundary, hence this stress possesses negative curl. 
The gradient of f/D is positive in this case, thus an 
anticyclonic gyre biased to the west is expected, generated 
by topographically induced conservation of potential vorticity. 
Theory predicts that:
<pre>
	A negative gradient of f/D (i.e. CORIOLIS = 1.0e-4) results 
	  in an eastward biased gyre.
	A flat bottom (BATHYMAX = 50) results in an unbiased gyre.
	Wind stress with positive curl (WIND_SPEED SCALE = -1) 
	  results in a cyclonic gyre with unaltered bias.
</pre>

Testing is performed for both 'z' and sigma vertical coordinates.
The 'z' model is also run with all the model diagnostics tested.
This test invokes the source/sink and particle tracking code. The
file test7_diag.prm may be altered to use an initial release 
distribution or point source for particle release.
The 'z' model is also run linked to waves and sediment transport
libraries and run with 3 sediment fractions under the influence
of wave induced bottom stress.

* The input parameter file for the 'z' model is : test7.prm
* The input parameter file for diagnostic tests is : test7_diag.prm
* The input parameter file for explicit mapping test is : test7_em.prm
* The input parameter file for grid refinement test is : test7_zoom.prm
* The input parameter file for the sigma model is : test7_s.prm
* The input parameter file for the sediment/wave model is : test7_sed.prm

Particle output is stored in out_pt.nc and may be viewed in 
matlab using 'plum'.
The polygon seeded input particle file was generated in matlab using:
seedpt('plum.poly','in7.nc',10000,0,-100,0,1,'part_poly.nc');
The polygon was creaded in 'plum' using extras->polygons.

* The 'z' output is stored in the file : out7_z.nc
* The diagnostic test output is stored in the file : out7_d.nc
Additional output is generated for this test in files:
<pre>
	alert.txt      	: ALERT mode summary.
	flushing.ts	: flushing tracer timeseries output.
	totals.ts	: total heat, mass, temp, salt and sediment output.
</pre>
* The explicit mapping test output is stored in the file : out7_e.nc
* The grid refinement test output is stored in the file : out7_r.nc
* The sigma output is stored in the file : out7_s.nc
* The sediment/wave output is stored in the file : out7_sed.nc

Multiple windows are tested in out7_1w.nc and out7_4w.nc. To check if
variables are the same, use (e.g. for eta):
```
ncdiff -v eta out7_1w.nc out7_4w.nc out7_4w.nc (type 'a'ppend on prompt).
```
To run type :
```
run_test7
```

To view output execute 'out' in Matlab, e.g.
```
matlab
>out
```
This script provides the velocity components mid-domain for each model.
A plot is saved in test7.eps
Totals plot saved in total.eps
Total sediment mass in sediment.eps
Flushing tracer plot saved in flush.eps
