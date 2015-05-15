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

# create separate folders for different compilers
#>GNUAll ############################
#mkdir -p GNU;   cp -r AMPS GNU/;  <#
#>IntelAll ##########################
#mkdir -p Intel; cp -r AMPS Intel/;<#
#>PGIAll ############################
#mkdir -p PGI;   cp -r AMPS PGI/;  <#

# copy job files to the AMPS directory on supercomputers
#>Pleiades #############################################
#cp AMPS/utility/TestScripts/test_amps.pleiades.job test_amps.job <#

# install AMPS
#>GNUAll ######################################
#cd $WorkDir/Tmp_AMPS_test/GNU/AMPS           #
#./Config.pl -install >& test_amps.log       <#
#>IntelAll ####################################
#cd $WorkDir/Tmp_AMPS_test/Intel/AMPS         #
#./Config.pl -install >& test_amps.log       <#
#>PGIAll ######################################
#cd $WorkDir/Tmp_AMPS_test/PGI/AMPS           #
#./Config.pl -install >& test_amps.log       <#

# compile AMPS tests

# GNU compiler

#>Pleiades>Yellowstone ########################
#module purge;                               <#

#>Pleiades ####################################
#module load mpi-openmpi/1.6.5-gcc;          <#

#>GNUAll ###################################################
#cd $WorkDir/Tmp_AMPS_test/GNU/AMPS                        #
#make test_compile COMPILE.mpicxx=mpicxx >>& test_amps.log<#

# Intel compiler 

#>Pleiades>Yellowstone ########################
#module purge;                               <#

#>Pleiades ####################################
#module load comp-intel/2015.0.090;           #
#module load mpi-openmpi/1.6.5-intel;        <#


#>IntelAll #################################################
#cd $WorkDir/Tmp_AMPS_test/Intel/AMPS                      #
#make test_compile COMPILE.mpicxx=mpicxx >>& test_amps.log<#

# PGI compiler

#>Pleiades>Yellowstone ########################
#module purge;                               <#

#>Pleiades ####################################
#module load comp-pgi/15.3;                   #
#module load mpi-mvapich2/1.8/pgi12.4;       <#


#>PGIAll ###################################################
#cd $WorkDir/Tmp_AMPS_test/PGI/AMPS                        #
#make test_compile COMPILE.mpicxx=mpicxx >>& test_amps.log<#


# Run test

# Super computers

#>Pleiades ####################################
#cd $WorkDir/Tmp_AMPS_test                    #
#/PBS/bin/qsub test_amps.job                 <#
#>Yellowstone #################################
#/usr/bin/bsub < test_amps.job               <#

# Other machines

# GNU compiled tests

#>GNUOther ####################################
#cd $WorkDir/Tmp_AMPS_test/GNU/AMPS           #
#make test_run >>& test_amps.log             <#

# Intel compiled tests

#>IntelOther ##################################
#cd $WorkDir/Tmp_AMPS_test/Intel/AMPS         #
#make test_run >>& test_amps.log             <#

# PGI compiled tests

#>PGIOther ####################################
#cd $WorkDir/Tmp_AMPS_test/PGI/AMPS           #
#make test_run >>& test_amps.log             <#

