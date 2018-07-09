### SHOC : Test8

This test exercises the wetting and drying capability of the
'z' model. A 50 by 27 by 24 closed basin is set up with sloping 
bottom (4m offshore and 1.2m inshore) and one offshore open 
boundary is forced with oscillating elevation at the offshore
boundary. The is initially dry throughout, fills with water as
elevation rises, then empties. 

Since only the 'z' model is a wetting and drying model, this 
test is performed for the 'z'  model only.
The input parameter file for the 'z' model is : test8.prm

The 'z' output is stored in the file : out8_z.nc

To run type :
```
run_test8
```
To view output execute 'out' in Matlab, e.g.
```
matlab
>out
```
This script provides the velocity components mid-domain for each model.
A plot is saved in test8.eps
