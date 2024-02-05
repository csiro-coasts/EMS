#!/bin/bash -e

if [ $# == 0 ] 
then
   echo " "
   echo " $# "
   echo "Purpose: Test mass balance of 1D vertical sediment transport model"
   echo "Usage:"
   echo "./run.sh prm1       (runs hd first and then sediment transport)"
   echo "./run.sh prm1 prm2  (runs only sediment transport; hd must be precomputed)"
   echo "The first parameter (prm1) must be set to:"
   echo "test_r   (to run resuspension test)"
   echo "test_c   (to run compaction test)"
   echo "test_d   (to run diffusion test)"
   echo "test_p   (to run pollutant test)"
   echo "test_rcd (to run resuspension, diffusion and compaction processes together)"
   echo "When the second parameter (prm2) is set to sedi,"
   echo "the model simulates only sediment transport (hd must be precomputed)"   
   echo "Date: August 2018"
   echo "Author: Nugzar Margvelashvili"
   echo ""
   exit;
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
cat ./vars sed_mass_end.txt > mass.txt
rm sed_mass_*.txt
#python plot.py
echo "Done"
