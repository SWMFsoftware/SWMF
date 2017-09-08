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
source /usr/share/modules/init/csh  

# source shell run commands (ONLY csh AND tcsh SHELLS ARE USED FOR NOW)
# to set CVSROOT or CVS_RSH variables
# note: it is better to have these variables set in the beginning of rc file
source $HOME/.cshrc

module load boost/1.62
# set the working directory
set WorkDir = $HOME  

#>Pleiades ############################
set WorkDir = /nobackup/`whoami`    

#>Yellowstone ###########
#set WorkDir =         <#

#>Stampede #############
# set WorkDir = $WORK <#

# Go to your home directory
cd $WorkDir

# Create a temporary directory for the tests
mkdir -p Tmp_AMPS_test
cd Tmp_AMPS_test

# Remove previous job scripts
rm -f *.job

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

#cp /home3/vtenishe/Table /nobackupp8/vtenishe/Tmp_AMPS_test/AMPS/MakefileTest

# create separate folders for different compilers
#>GNUAll ############################
rm -rf GNU
mkdir -p GNU;   cp -r AMPS GNU/;  
#>IntelAll ##########################
rm -rf Intel
mkdir -p Intel; cp -r AMPS Intel/;
#>PGIAll ############################
rm -rf PGI
mkdir -p PGI;   cp -r AMPS PGI/;  

# install AMPS
#>GNUAll ###################################################################
cd $WorkDir/Tmp_AMPS_test/GNU/AMPS                                        #
./Config.pl -install -compiler=gfortran,gcc_mpicc    >& test_amps.log    
#>IntelAll #################################################################
cd $WorkDir/Tmp_AMPS_test/Intel/AMPS                                      #
./Config.pl -install -compiler=ifortmpif90,iccmpicxx >& test_amps.log    
#>PGIAll ###################################################################
cd $WorkDir/Tmp_AMPS_test/PGI/AMPS                                        #
./Config.pl -install -compiler=pgf90,pgccmpicxx -cpp-compiler=pgc++ -link-option=-L/nasa/sgi/mpt/2.14r19/lib,-lmpi++,-lmpi       >& test_amps.log    

# copy job files to the AMPS directory on supercomputers
#>Pleiades ###############################################
cp utility/TestScripts/test_amps.pleiades.*.job ../..
#>Stampede ###############################################
#cp AMPS/utility/TestScripts/test_amps.stampede.*.job . <#


# Exeptions
echo -n "Set Exeptions....."

#>Valeriy ##################################################################
#cd $WorkDir/Tmp_AMPS_test/Intel/AMPS                                                                      #
#./Config.pl -link-option=-lc++ -install -compiler=ifort,iccmpicxx -link-option=-cxxlib >>& test_amps.log <# 

#>Pleiades ##############################################
cd $WorkDir/Tmp_AMPS_test/Intel/AMPS                   #
./Config.pl -install -compiler=ifort,icc
./Config.pl -cpp-link-option=-lmpi 
./Config.pl -f-link-option=-lmpi                     #
./Config.pl -cpplib-rm=-lmpi_cxx
./Config.pl -noopenmp

cd $WorkDir/Tmp_AMPS_test/GNU/AMPS                     #
./Config.pl -cpplib-rm=-lmpi_cxx
./Config.pl -f-link-option=-lmpi++                  

cd $WorkDir/Tmp_AMPS_test/PGI/AMPS
./Config.pl -f-link-option=-pgf90libs,-pgc++libs

cd $WorkDir/Tmp_AMPS_test
rm -f AmpsCompilingIntelComplete
rm -f AmpsCompilingGNUComplete
rm -f AmpsCompilingPGIComplete
rm -f AmpsTestComplete

echo " done."

# compile AMPS tests

# GNU compiler

#>Valeriy ######################################################################
#source $WorkDir/Tmp_AMPS_test/AMPS/utility/TestScripts/CompilerSetup/set_mpi_gnu.valeriy <#

#>Pleiades ####################################
module load gcc/6.2                         #
module load mpi-sgi;                        
$HOME/bin/CompileGNUPleiades.sh &

#>Stampede ####################################
#module load gcc/4.9.1                        #
#module load mvapich2/2.1                    <#

#>GNUAll ######################################
#jecho -n "Compiling GNU....."                 # 
#cd $WorkDir/Tmp_AMPS_test/GNU/AMPS           #
#make test_compile >>& test_amps.log          #
#echo " done."                                

# Intel compiler 

#>Valeriy ########################################################################
#source $WorkDir/Tmp_AMPS_test/AMPS/utility/TestScripts/CompilerSetup/set_mpi_intel.valeriy <#

#>Pleiades>Yellowstone>Stampede ###############
module purge;                               

#>Pleiades ####################################
$HOME/bin/CompileIntelPleiades.sh &

#>Stampede ####################################
#module load intel/15.0.2;                    #
#module load mvapich2/2.1;                   <#


#>IntelAll ####################################
#echo -n "Compiling Intel....."               #
#cd $WorkDir/Tmp_AMPS_test/Intel/AMPS         #
#make test_compile >>& test_amps.log          #
#echo " done."                                

# PGI compiler

#>Pleiades>Yellowstone ########################
module purge;                               

$HOME/bin/CompilePGIPleiades.sh &

#>Pleiades ####################################


#>PGIAll ######################################
#echo -n "Compiling PGI....."                 # 
#cd $WorkDir/Tmp_AMPS_test/PGI/AMPS           #
#make test_compile >>& test_amps.log          #
#echo " done."                                


# Run test
# Super computers
#>Pleiades ####################################
cd $WorkDir/Tmp_AMPS_test

#waite untill all compilation is complete
while ((! -f AmpsCompilingIntelComplete) || (! -f AmpsCompilingGNUComplete) || (! -f AmpsCompilingPGIComplete)) 
  sleep 60
end

rm -f AmpsCompilingIntelComplete
rm -f AmpsCompilingGNUComplete
rm -f AmpsCompilingPGIComplete

echo Compiling of AMPS is completed

#########################################

rm -f AmpsTestComplete

foreach job (test_amps.pleiades.all.*job)                #
  #execute the next part of the tests
  /PBS/bin/qsub $job

  #waite while the set of test is finished
  #and then proceed with the new part of 
  #the tests 
  while (! -f AmpsTestComplete) 
    sleep 60
  end

  sleep 300 
  rm -f AmpsTestComplete
end                                         


#>Stampede ####################################
#set submit = '/usr/bin/sbatch'              <#
#>Stampede ####################################
#cd $WorkDir/Tmp_AMPS_test                    #
#@ delay = 2                                  #
#foreach job (test_amps.*.job)                #
#  echo "$submit $job"|at now+$delay minutes  #
#  @ delay = $delay + 121                     #
#end                                         <#
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

