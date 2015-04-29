#!/bin/csh
# This script checks out the latest version of the AMPS,
# executes the AMPS tests and sends the results to the server:
# herot.engin.umich.edu for further processing
#******************************************************************************
# the script is specific to PLEIADES supercomputer;
# it performs test with 2 different C++ compilers: GCC and Intel
#******************************************************************************
# The script can be executed by simply typing
#
# ./run_test_amps
#
# To run the AMPS tests perodically, you can use the crontab facility.
# Type 'crontab -e' to add a new entry to your crontab.
# Here is an example entry for nightly runs at 12:30 am:
#
# 30 0 * * * $HOME/bin/run_test_amps.pleiades.sh


# init command to load modules
source /usr/share/modules/init/csh

# Go to your home directory
cd $HOME

# Create a temporary directory for the tests
mkdir -p Tmp_AMPS_test
cd Tmp_AMPS_test

#------------------------------------------------------------------------------
# load Intel compiler
module load mpi-openmpi/1.6.5-intel
module load comp-intel/2015.0.090

# Remove the previous test directory if necessary
rm -rf AMPS

# Checkout and install the latest code
cvs co -D "`date +%m/%d/%Y` 0:30" AMPS 
cd AMPS
cvs co -D "`date +%m/%d/%Y` 0:30" AMPS_data 

# Install the AMPS
./Config.pl -install > test_amps.log

# Run test
make test >> test_amps.log
          
# Create file with results' summary                                                  
ls -ltr  *diff > test_amps.res
echo '=============================================================='\
                >> test_amps.res
head -100 *diff >> test_amps.res

# Copy test_amps.log and test_amps.res to the server
 scp test_amps.res test_amps.log \
    dborovik@herot.engin.umich.edu:Sites/Current/`hostname -s`/ 
#------------------------------------------------------------------------------
# load GCC compiler
module load mpi-openmpi/1.6.5-gcc

# Remove the previous srcTemp directory
cd ..
rm -rf AMPS

# Checkout and install the latest code
cvs co -D "`date +%m/%d/%Y` 0:30" AMPS 
cd AMPS
cvs co -D "`date +%m/%d/%Y` 0:30" AMPS_data 

# Install the AMPS
./Config.pl -install > test_amps.log

# Run test
make test >> test_amps.log
          
# Create file with results' summary
ls -ltr  *diff > test_amps.res
echo '=============================================================='\
                >> test_amps.res
head -100 *diff >> test_amps.res

# Copy test_amps.log and test_amps.res to the server
scp test_amps.res test_amps.log \
   dborovik@herot.engin.umich.edu:Sites/Current/`hostname -s`/ 


