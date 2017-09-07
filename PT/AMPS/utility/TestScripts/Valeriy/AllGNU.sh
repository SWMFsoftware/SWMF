#!/bin/csh

set WorkDir = $HOME

source $WorkDir/Tmp_AMPS_test/AMPS/utility/TestScripts/CompilerSetup/set_mpi_gnu.valeriy
echo -n "Compiling GNU....."

cd $WorkDir/Tmp_AMPS_test/GNU/AMPS  
./Config.pl -link-option=-lmpi_cxx
make test_compile >>& test_amps.log

cd $WorkDir/Tmp_AMPS_test/GNU/AMPS/

echo " done."

echo -n "Executing tests GNU....."
make TESTMPIRUN4="mpirun -np 4"  MPIRUN="mpirun -np 4" TESTMPIRUN1="export DYLD_LIBRARY_PATH=/Users/ccmc/boost/lib:/opt/intel/compilers_and_libraries/mac/lib;mpirun -np 1" test_run >>& test_amps.log
echo " done."
