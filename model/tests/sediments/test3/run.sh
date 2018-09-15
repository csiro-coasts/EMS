#!/bin/bash

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

if [ $# == 1 ] 
then
   echo "Cleaning up"
   rm outputs/*
   rm tran_outputs/*
   rm inputs/in_test.nc
   rm inputs/test_trans_1990-01.nc
   rm diag.txt runlog sedlog.txt setup.txt trans.mnc
   echo "Running $1: hd first and then sediment transport"
  ./shoc -g inputs/test.prm inputs/in_test.nc
  ./shoc -p inputs/test.prm
  ./shoc -t inputs/$1.tran
fi

if [ $# == 2 ] 
then
   echo "Running $1: only sediment transport"
   rm tran_outputs/*
   rm sedlog.txt
  ./shoc -t inputs/$1.tran
fi

echo "Done"

