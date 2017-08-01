#!/bin/csh

set WorkDir = $HOME

source $WorkDir/Tmp_AMPS_test/AMPS/utility/TestScripts/CompilerSetup/set_mpi_gnu.valeriy
echo -n "Compiling GNU....."

cd $WorkDir/Tmp_AMPS_test/GNU/AMPS  
./Config.pl -link-option=-lmpi_cxx
make test_compile >>& test_amps.log

#making of test_CG-PostProcess--Read-Trajectories--off is not standard - used a custom build
cd $WorkDir/Tmp_AMPS_test/GNU/AMPS/srcCG-PostProcess; ./Config.pl -extraflag=-fopenmp 
cd $WorkDir/Tmp_AMPS_test/GNU/AMPS/utility/PostProcess; ./Config.pl -extraflag=-fopenmp

cd $WorkDir/Tmp_AMPS_test/GNU/AMPS/
make EXTRALINKEROPTIONS="./utility/PostProcess/pproc.a -fopenmp" test_CG-PostProcess--Read-Trajectories--off_compile -j
make test_CG-PostProcess--Read-Trajectories--off_rundir

echo " done."

echo -n "Executing tests GNU....."
make TESTMPIRUN1="mpirun -np 1" test_run >>& test_amps.log
echo " done."
