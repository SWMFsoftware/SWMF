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

# init command to load modules      
#>Pleiades ############################
#source /usr/share/modules/init/csh  <#

# source shell run commands (ONLY csh AND tcsh SHELLS ARE USED FOR NOW)
# to set CVSROOT or CVS_RSH variables
# note: it is better to have these variables set in the beginning of rc file
source $HOME/.cshrc

# set the working directory
set WorkDir = $HOME  

#>Pleiades ############################
#set WorkDir = /nobackup/`whoami`    <#

#>Yellowstone ###########
#set WorkDir =         <#

#>Stampede #############
 set WorkDir = $SCRATCH

# Go to your home directory
cd $WorkDir

# Create a temporary directory for the tests
mkdir -p Tmp_AMPS_test
cd Tmp_AMPS_test

# Remove the previous test directory if necessary
rm -rf AMPS */AMPS

# Checkout the latest code version
cvs co -D "`date +%m/%d/%Y` 23:20" AMPS 
cd AMPS
cvs co -D "`date +%m/%d/%Y` 23:20" AMPS_data 
cd ..

# Update data files for test at supercomputers
#>Pleiades>Yellowstone>Stampede #########################
rsync -r amps@tower-left.engin.umich.edu:/Volumes/Data01/AMPS_DATA_TEST/ $WorkDir/AMPS_DATA_TEST 

# create separate folders for different compilers
#>GNUAll ############################
mkdir -p GNU;   cp -r AMPS GNU/;  
#>IntelAll ##########################
mkdir -p Intel; cp -r AMPS Intel/;
#>PGIAll ############################
#mkdir -p PGI;   cp -r AMPS PGI/;  <#

# install AMPS
#>GNUAll ###################################################################
cd $WorkDir/Tmp_AMPS_test/GNU/AMPS                                        #
./Config.pl -install -compiler=gfortran,gcc_mpicc    >& test_amps.log    
#>IntelAll #################################################################
cd $WorkDir/Tmp_AMPS_test/Intel/AMPS                                      #
./Config.pl -install -compiler=ifortmpif90,iccmpicxx >& test_amps.log    
#>Valeriy ##################################################################
#cd $WorkDir/Tmp_AMPS_test/Intel/AMPS                                      #
#./Config.pl -link-option=-lc++ -install -compiler=ifort,iccmpicxx -link-option=-cxxlib >>& test_amps.log<#
#>PGIAll ###################################################################
#cd $WorkDir/Tmp_AMPS_test/PGI/AMPS                                        #
#./Config.pl -install -compiler=pgf90,pgccmpicxx      >& test_amps.log    <#


# copy job files to the AMPS directory on supercomputers
# #>Pleiades ###############################################
# #cp AMPS/utility/TestScripts/test_amps.pleiades.*.job . <#
# #>Stampede ###############################################
cp utility/TestScripts/test_amps.stampede.*.job ../.. 

# Compile AMPS tests
cd $WorkDir/Tmp_AMPS_test
rm -f AmpsCompilingIntelComplete
rm -f AmpsCompilingGNUComplete

AMPS/utility/TestScripts/Stampede/CompileGNUPStampede.sh &
AMPS/utility/TestScripts/Stampede/CompileIntelStampede.sh &

#waite untill all compilation is complete
while ((! -f AmpsCompilingIntelComplete) || (! -f AmpsCompilingGNUComplete) )
  sleep 60
end

rm -f AmpsCompilingIntelComplete
rm -f AmpsCompilingGNUComplete

echo Compiling of AMPS is completed

# Run test
#>Stampede ####################################
set submit = '/usr/bin/sbatch'              

cd $WorkDir/Tmp_AMPS_test                    #
rm -f AmpsTestDone

foreach job (test_amps.stampede.all*.job)                #
  sbatch $job 

  while (! -f AmpsTestDone)
    sleep 60
  end

  sleep 180
  rm -f AmpsTestDone
end                                         

#>Yellowstone #################################
#/usr/bin/bsub < test_amps.job               <#

# Other machines

# GNU compiled tests

#>Valeriy ########################################################################
#source $WorkDir/Tmp_AMPS_test/AMPS/utility/TestScripts/CompilerSetup/set_mpi_gnu.valeriy <#

#>GNUOther ####################################
#cd $WorkDir/Tmp_AMPS_test/GNU/AMPS           #
#make test_run >>& test_amps.log             <#

# Intel compiled tests

#>Valeriy ########################################################################
#source $WorkDir/Tmp_AMPS_test/AMPS/utility/TestScripts/CompilerSetup/set_mpi_intel.valeriy <#

#>IntelOther ##################################
#cd $WorkDir/Tmp_AMPS_test/Intel/AMPS         #
#make test_run >>& test_amps.log             <#

# PGI compiled tests

#>PGIOther ####################################
#cd $WorkDir/Tmp_AMPS_test/PGI/AMPS           #
#make test_run >>& test_amps.log             <#

