#!/bin/csh
# This script checks out the latest version of the AMPS,
# executes the AMPS tests
# This script is meant to be used by the AMPS developers.

# The script can be executed by simply typing
#
# ./run_test_amps.sh
#
# To run the AMPS tests perodically, you can use the crontab facility.
# Type 'crontab -e' to add a new entry to your crontab.
# Here is an example entry for nightly runs at 12:30 am:
#
# 30 0 * * * $HOME/bin/run_test_amps.sh

# NOTE:
# the body of the script is divided into blocks of the format
#  #>BlockName ##############
#  #command_1               #
#  # ...                    #
#  #command_last           <#
# certain blocks will be uncommented at the test installation,

# source shell run commands (ONLY csh AND tcsh SHELLS ARE USED FOR NOW)
# to set CVSROOT or CVS_RSH variables
# note: it is better to have these variables set in the beginning of rc file
source $HOME/.cshrc

#Set the working directory
set WorkDir = $HOME  

#Go to your home directory
cd $WorkDir

#Create a temporary directory for the tests
mkdir -p Tmp_AMPS_test
cd Tmp_AMPS_test

#Remove the previous test directory if necessary
rm -rf AMPS */AMPS

#Checkout the latest code version
cvs co -D "`date +%m/%d/%Y` 23:20" AMPS 
cd AMPS
cvs co -D "`date +%m/%d/%Y` 23:20" AMPS_data 
cd ..

#Create separate folders for different compilers
rm -rf GNU
mkdir -p GNU;   cp -r AMPS GNU/;

rm -rf Intel
mkdir -p Intel; cp -r AMPS Intel/;

rm -rf PGI
mkdir -p PGI;   cp -r AMPS PGI/;

#Install AMPS
cd $WorkDir/Tmp_AMPS_test/GNU/AMPS                                         
./Config.pl -install -compiler=gfortran,gcc_mpicc -f-link-option=-lmpi_cxx   >& test_amps.log    

cd $WorkDir/Tmp_AMPS_test/Intel/AMPS                                       
./Config.pl -f-link-option=-lc++ -install -compiler=ifort,iccmpicxx -link-option=-cxxlib >& test_amps.log

cd $WorkDir/Tmp_AMPS_test/PGI/AMPS                                         
./Config.pl -f-link-option=-lc++ -install -compiler=pgf90,pgccmpicxx  -link-option=-lc++    >& test_amps.log    

#Execute the tests
$WorkDir/Tmp_AMPS_test/AMPS/utility/TestScripts/Valeriy/AllGNU.sh & 
$WorkDir/Tmp_AMPS_test/AMPS/utility/TestScripts/Valeriy/AllIntel.sh &   
$WorkDir/Tmp_AMPS_test/AMPS/utility/TestScripts/Valeriy/AllPGI.sh &


