#!/bin/csh

source $HOME/.cshrc
set WorkDir = /glade2/scratch2/`whoami`

module swap intel gnu

echo -n "Compiling GNU....."                 
cd $WorkDir/Tmp_AMPS_test/GNU/AMPS           
make test_compile >>& test_amps.log          
echo " done."

cd $WorkDir/Tmp_AMPS_test
echo Done > AmpsCompilingGNUComplete
