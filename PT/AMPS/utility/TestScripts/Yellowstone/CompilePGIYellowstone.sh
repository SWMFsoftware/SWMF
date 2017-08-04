#!/bin/csh

source $HOME/.cshrc
set WorkDir = /glade2/scratch2/`whoami`

module swap intel pgi 

echo -n "Compiling PGI....."                 
cd $WorkDir/Tmp_AMPS_test/PGI/AMPS           
make test_compile >>& test_amps.log          
echo " done."

cd $WorkDir/Tmp_AMPS_test
echo Done > AmpsCompilingPGIComplete
