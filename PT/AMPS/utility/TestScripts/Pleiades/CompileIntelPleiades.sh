#!/bin/csh

source /usr/share/modules/init/csh
source $HOME/.cshrc

set WorkDir = /nobackup/`whoami`

module purge;
module load comp-intel/2016.2.181;           
module load mpi-sgi/mpt;

echo -n "Compiling Intel....."                 
cd $WorkDir/Tmp_AMPS_test/Intel/AMPS           
make test_compile >>& test_amps.log          
echo " done."

cd $WorkDir/Tmp_AMPS_test
echo Done > AmpsCompilingIntelComplete
