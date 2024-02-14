#!/bin/bash -e

if [ $# == 1 ] 
then
   echo " "
   echo "Purpose: Test mass balance of 1D vertical sediment transport models"
   echo "Usage: ./run_all.sh"
   echo "Should take less than 5 min to complete"
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

echo "Running hd"
$SHOC -g inputs/test.prm inputs/in_test.nc
$SHOC -p inputs/test.prm
echo "OK: HD forcing files complete"
echo ""

echo "Running Diffusion test"
$SHOC -t inputs/test_d.tran
awk 'BEGIN{a=0;} {if ($NR eq 3) a=$4;} END{if((a-$4)<0.1) \
print("OK: mass conservation of dissolved tracer"); \
else print("ERROR: mass conservation of dissolved tracer");}' sed_mass_end.txt
echo ""

echo "Running Resuspension test"
rm -f tran_outputs/*   
rm -f sedlog.txt sed_mass_*.txt
$SHOC -t inputs/test_r.tran
awk 'BEGIN{a=0;} {if ($NR eq 3) a=$5;} END{if((a-$5)<0.1) \
print("OK: mass conservation of gravel"); \
else print("ERROR: mass conservation of gravel");}' sed_mass_end.txt
awk 'BEGIN{a=0;} {if ($NR eq 3) a=$6;} END{if((a-$6)<0.1) \
print("OK: mass conservation of sand"); \
else print("ERROR: mass conservation of sand");}' sed_mass_end.txt
awk 'BEGIN{a=0;} {if ($NR eq 3) a=$7;} END{if((a-$7)<0.1)\
print("OK: mass conservation of mud"); \
else print("ERROR: mass conservationof mud");}' sed_mass_end.txt
echo ""

echo "Running Compaction test"
rm -f tran_outputs/*   
rm -f sedlog.txt sed_mass_*.txt
$SHOC -t inputs/test_c.tran
awk 'BEGIN{a=0;} {if ($NR eq 3) a=$4;} END{if((a-$4)<0.1) \
print("OK: mass conservation of dissolved tracer"); \
else print("ERROR: mass conservation of dissolved tracer");}' sed_mass_end.txt
awk 'BEGIN{a=0;} {if ($NR eq 3) a=$7;} END{if((a-$7)<0.1) \
print("OK: mass conservation of mud"); \
else print("ERROR: mass conservation of mud");}' sed_mass_end.txt
echo ""

echo "Running Desorption of Pollutant test"
rm -f tran_outputs/*   
rm -f sedlog.txt sed_mass_*.txt
$SHOC -t inputs/test_p.tran
awk 'BEGIN{a=0;} {if ($NR eq 3) a=$4+$8;} END{if((a-$4-$8)<0.1) \
print("OK: mass conservation of pollutant"); \
else print("ERROR: mass conservation of pollutant");}' sed_mass_end.txt
echo ""

echo "Running Resuspension, Diffusion and Compaction together"
rm -f tran_outputs/*
rm -f sedlog.txt sed_mass_*.txt
$SHOC -t inputs/test_rcd.tran
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
echo "1D TESTS COMPLETE"
