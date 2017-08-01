#!/bin/csh

set WorkDir = $HOME

source $WorkDir/Tmp_AMPS_test/AMPS/utility/TestScripts/CompilerSetup/set_mpi_pgi.valeriy
echo -n "Compiling PGI....."

cd $WorkDir/Tmp_AMPS_test/PGI/AMPS  
make test_compile >>& test_amps.log

#making of test_CG-PostProcess--Read-Trajectories--off is not standard - used a custom build
cd $WorkDir/Tmp_AMPS_test/PGI/AMPS/srcCG-PostProcess; ./Config.pl -extraflag=-openmp
cd $WorkDir/Tmp_AMPS_test/PGI/AMPS/utility/PostProcess; ./Config.pl -extraflag=-openmp

cd $WorkDir/Tmp_AMPS_test/PGI/AMPS/
make EXTRALINKEROPTIONS="./utility/PostProcess/pproc.a -openmp" test_CG-PostProcess--Read-Trajectories--off_compile -j
mpicxx -o amps srcTemp/main/main.a srcTemp/libAMPS.a -lstdc++ /Users/vtenishe/SPICE/Toolkit/cspice/lib/cspice.a ./utility/PostProcess/pproc.a
make test_CG-PostProcess--Read-Trajectories--off_rundir

echo " done."

echo -n "Executing tests PGI....."
make TESTMPIRUN1="mpirun -np 1" test_run >>& test_amps.log
echo " done."
