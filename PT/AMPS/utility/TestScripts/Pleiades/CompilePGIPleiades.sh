#!/bin/csh

source /usr/share/modules/init/csh
source $HOME/.cshrc

set WorkDir = /nobackup/`whoami`

module unload comp-intel
module load comp-pgi 

echo -n "Compiling PGI....."                 
cd $WorkDir/Tmp_AMPS_test/PGI/AMPS           
make test_compile >>& test_amps.log          
echo " done."

cd $WorkDir/Tmp_AMPS_test
echo Done > AmpsCompilingPGIComplete
