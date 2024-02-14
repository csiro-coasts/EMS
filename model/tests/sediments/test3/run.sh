#!/bin/bash -e

if [ $# == 0 ] 
then
   echo " "
   echo " $# "
   echo "Purpose: 3D sediment transport model in an idealised estuarine channel"
   echo "Usage:"
   echo "./run.sh prm1       (runs hd first and then sediment transport)"
   echo "./run.sh prm1 prm2  (runs only sediment transport; hd must be precomputed)"
   echo "Takes about 10 min to complete"
   echo "First parameter (prm1) must be set to either"
   echo "test_basic to run basic configuration or"
   echo "test_advanced to run advanced configuration"
   echo "Second parameter (prm2) can be set to"
   echo "sedi to simulate only sediment transport (hd must be precomputed)"   
   echo "Date: September 2018"
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
   rm -f diag.txt runlog sedlog.txt setup.txt trans.mnc
   echo "Running $1: hd first and then sediment transport"
  $SHOC -g inputs/test.prm inputs/in_test.nc
  $SHOC -p inputs/test.prm
  $SHOC -t inputs/$1.tran
fi

if [ $# == 2 ] 
then
   echo "Running $1: only sediment transport"
   rm -f tran_outputs/*
   rm -f sedlog.txt
   $SHOC -t inputs/$1.tran
fi

echo "Done"
