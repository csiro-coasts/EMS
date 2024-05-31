#!/bin/bash -e

if [ $# == 1 ] 
then
   echo " "
   echo "Purpose: Test mass balance of a 3D sediment transport model"
   echo "Usage: ./run_all.sh"
   echo "Takes about 5 min to run"
   echo "Date: August 2018"
   echo "Author: Nugzar Margvelashvili"
   echo ""
   exit;
else
  echo ""
fi
  
SHOC="${SHOC:-./shoc}"
mkdir -p outputs tran_outputs

echo "Cleaning up"
rm -f outputs/*
rm -f tran_outputs/*
rm -f inputs/in_test.nc
rm -f inputs/test_trans_1990-01.nc
rm -f diag.txt runlog mass.txt sedlog.txt setup.txt trans.mnc sed_mass_*.txt
echo ""

echo "Simulating hd in a wind driven, closed basin"
$SHOC -g inputs/test.prm inputs/in_test.nc
$SHOC -p inputs/test.prm
echo "OK: HD forcing files complete"
echo ""


echo "Simulating transport of particulate and dissolver tracers"
$SHOC -t inputs/test_flat.tran
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
rm -f sed_mass_*.txt
echo "3D TEST COMPLETE"
