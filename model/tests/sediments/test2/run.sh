#!/bin/bash

if [ $# == 0 ] 
then
   echo " "
   echo " $# "
   echo "Purpose: Test mass balance of 3D sediment transport model"
   echo "Usage:"
   echo "./run.sh test_flat       (runs hd first and then sediment transport)"
   echo "./run.sh test_flat sedi  (runs only sediment transport; hd must be precomputed)"
   echo "Date: August 2018"
   echo "Author: Nugzar Margvelashvili"
   echo ""
   exit
else
   echo "run './run.sh' for info"
fi

if [ $# == 1 ] 
then
   echo "Cleaning up"
   rm outputs/*
   rm tran_outputs/*
   rm inputs/in_test.nc
   rm inputs/test_trans_1990-01.nc
   rm diag.txt runlog mass.txt sedlog.txt setup.txt trans.mnc sed_mass_*.txt
   echo "Running $1: hd first and then sediment transport"
  ./shoc -g inputs/test.prm inputs/in_test.nc
  ./shoc -p inputs/test.prm
  ./shoc -t inputs/$1.tran
fi

if [ $# == 2 ] 
then
   echo "Running $1: only sediment transport"
   rm tran_outputs/*
   rm mass.txt sedlog.txt sed_mass_*.txt
  ./shoc -t inputs/$1.tran
fi

# post-processing
awk 'BEGIN{a=0;} {if ($NR eq 3) a=$4;} END{if((a-$4)<0.1) \
print("OK: mass conservation of dissolved tracer"); \
else print("ERROR: mass conservation of dissolved tracer");}' sed_mass_end.txt
awk 'BEGIN{a=0;} {if ($NR eq 3) a=$5;} END{if((a-$5)<0.1) \
print("OK: mass conservation of gravel");\
else print("ERROR: mass conservation of gravel");}' sed_mass_end.txt
awk 'BEGIN{a=0;} {if ($NR eq 3) a=$6;} END{if((a-$6)<0.1) \
print("OK: mass conservation of sand");\
else print("ERROR: mass conservation of sand");}' sed_mass_end.txt
awk 'BEGIN{a=0;} {if ($NR eq 3) a=$7;} END{if((a-$7)<0.1) \
print("OK: mass conservation of mud"); \
else print("ERROR: mass conservation of mud");}' sed_mass_end.txt

rm sed_mass_*.txt
echo "Done"

