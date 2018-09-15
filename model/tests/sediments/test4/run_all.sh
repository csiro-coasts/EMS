#!/bin/bash

if [ $# == 1 ] 
then
   echo " "
   echo "Purpose: Test mass balance of a coupled 3D hd-sediment transport model"
   echo "Usage: ./run_all.sh"
   echo "Takes about 7 min to complete the run"
   echo "Date: August 2018"
   echo "Author: Nugzar Margvelashvili"
   echo ""
   exit;
else
  echo ""
fi
  
echo "Cleaning up"
rm outputs/*
rm inputs/in_test.nc
rm diag.txt runlog sedlog.txt setup.txt sed_mass_*.txt
echo ""

echo "Simulating coupled hd and sediment transport in a wind driven, closed basin"
./shoc -g inputs/test.prm inputs/in_test.nc
./shoc -p inputs/test.prm
echo "OK: Run complete"
echo ""

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

echo ""
rm  sed_mass_*.txt
echo "3D TEST COMPLETE"
