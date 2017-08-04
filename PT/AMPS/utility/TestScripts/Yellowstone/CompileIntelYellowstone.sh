#!/bin/csh

source $HOME/.cshrc
set WorkDir = /glade2/scratch2/`whoami`

echo -n "Compiling Intel....."                 
cd $WorkDir/Tmp_AMPS_test/Intel/AMPS           
make test_compile >>& test_amps.log          
echo " done."

cd $WorkDir/Tmp_AMPS_test
echo Done > AmpsCompilingIntelComplete
