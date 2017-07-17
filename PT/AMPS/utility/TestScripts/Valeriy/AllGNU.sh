#!/bin/csh

set WorkDir = $HOME

source $WorkDir/Tmp_AMPS_test/AMPS/utility/TestScripts/CompilerSetup/set_mpi_gnu.valeriy
echo -n "Compiling GNU....."

cd $WorkDir/Tmp_AMPS_test/GNU/AMPS  
./Config.pl -link-option=-lmpi_cxx
make test_compile >>& test_amps.log
echo " done."

echo -n "Executing tests GNU....."
make test_run >>& test_amps.log
echo " done."
