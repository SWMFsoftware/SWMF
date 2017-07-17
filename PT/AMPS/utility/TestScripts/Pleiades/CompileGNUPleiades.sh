#!/bin/csh

source /usr/share/modules/init/csh
source $HOME/.cshrc

set WorkDir = /nobackup/`whoami`

module purge;
module load gcc/6.2                         
module load mpi-sgi;

echo -n "Compiling GNU....."                 
cd $WorkDir/Tmp_AMPS_test/GNU/AMPS           
make test_compile >>& test_amps.log          
echo " done."

cd $WorkDir/Tmp_AMPS_test
echo Done > AmpsCompilingGNUComplete
