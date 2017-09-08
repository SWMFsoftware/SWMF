#!/bin/csh

set WorkDir = $HOME

source $WorkDir/Tmp_AMPS_test/AMPS/utility/TestScripts/CompilerSetup/set_mpi_pgi.valeriy
echo -n "Compiling PGI....."

cd $WorkDir/Tmp_AMPS_test/PGI/AMPS  
make test_compile >>& test_amps.log


echo " done."

echo -n "Executing tests PGI....."
make TESTMPIRUN4="mpirun -np 4" MPIRUN="mpirun -np 4" test_run >>& test_amps.log
echo " done."
