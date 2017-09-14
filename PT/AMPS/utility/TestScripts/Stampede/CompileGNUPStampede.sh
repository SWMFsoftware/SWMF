#!/bin/csh

set WorkDir = $SCRATCH

module load gcc 

echo -n "Compiling GNU....."                 

echo $WorkDir/Tmp_AMPS_test/GNU/AMPS 
echo "Begin"

cd $WorkDir/Tmp_AMPS_test/GNU/AMPS           
make test_compile >>& test_amps.log          
echo " done."

cd $WorkDir/Tmp_AMPS_test
echo Done > AmpsCompilingGNUComplete

