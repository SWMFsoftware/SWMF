#!/bin/csh

source /usr/share/modules/init/csh
source $HOME/.cshrc

set WorkDir = /nobackup/`whoami`

module purge;
module load comp-pgi/17.1           
module load mpi-sgi/mpt

echo -n "Compiling PGI....."                 
cd $WorkDir/Tmp_AMPS_test/PGI/AMPS           
make test_compile >>& test_amps.log          
echo " done."

cd $WorkDir/Tmp_AMPS_test
echo Done > AmpsCompilingPGIComplete
