#!/bin/bash -e

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

SHOC="${SHOC:-./shoc}"
mkdir -p outputs tran_outputs

if [ $# == 1 ] 
then
   echo "Cleaning up"
   rm -f outputs/*
   rm -f tran_outputs/*
   rm -f inputs/in_test.nc
   rm -f inputs/test_trans_1990-01.nc
   rm -f diag.txt runlog mass.txt sedlog.txt setup.txt trans.mnc sed_mass_*.txt
   echo "Running $1: hd first and then sediment transport"
   $SHOC -g inputs/test.prm inputs/in_test.nc
   $SHOC -p inputs/test.prm
   $SHOC -t inputs/$1.tran
fi

if [ $# == 2 ] 
then
   echo "Running $1: only sediment transport"
   rm -f tran_outputs/*
   rm -f mass.txt sedlog.txt sed_mass_*.txt
   $SHOC -t inputs/$1.tran
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

rm -f sed_mass_*.txt
echo "Done"
