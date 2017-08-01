#!/bin/csh

set WorkDir = $HOME

source $WorkDir/Tmp_AMPS_test/AMPS/utility/TestScripts/CompilerSetup/set_mpi_intel.valeriy
echo -n "Compiling Intel....."

cd $WorkDir/Tmp_AMPS_test/Intel/AMPS; 
make test_compile >>& test_amps.log

#making of test_CG-PostProcess--Read-Trajectories--off is not standard - used a custom build
cd $WorkDir/Tmp_AMPS_test/Intel/AMPS/srcCG-PostProcess; ./Config.pl -extraflag=-fopenmp
cd $WorkDir/Tmp_AMPS_test/Intel/AMPS/utility/PostProcess; ./Config.pl -extraflag=-fopenmp

cd $WorkDir/Tmp_AMPS_test/Intel/AMPS/
make EXTRALINKEROPTIONS="./utility/PostProcess/pproc.a -fopenmp" test_CG-PostProcess--Read-Trajectories--off_compile -j
mpicxx -o amps srcTemp/main/main.a srcTemp/libAMPS.a -lstdc++ -lmpi_cxx    /Users/vtenishe/SPICE/Toolkit/cspice/lib/cspice.a ./utility/PostProcess/pproc.a /Users/vtenishe/Debugger/eclipse-workspace/MERCURYAMPS/AMPS/libiomp5.a
make test_CG-PostProcess--Read-Trajectories--off_rundir

echo " done."

echo -n "Executing tests Intel....."
make TESTMPIRUN1="mpirun -np 1" test_run >>& test_amps.log
echo " done."
