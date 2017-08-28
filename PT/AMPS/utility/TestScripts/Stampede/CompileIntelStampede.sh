#!/bin/csh

set WorkDir = $SCRATCH 

module load intel/17.0.4 

echo -n "Compiling Intel....."                 
cd $WorkDir/Tmp_AMPS_test/Intel/AMPS           


echo $WorkDir/Tmp_AMPS_test/Intel/AMPS
echo "Start"


make test_compile >>& test_amps.log          
echo " done."

cd $WorkDir/Tmp_AMPS_test
echo Done > AmpsCompilingIntelComplete



