#!/bin/csh
# This script checks out the latest version of the AMPS,
# executes the AMPS tests and sends the results to the server:
# herot.engin.umich.edu for further processing
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

# Go to your home directory
cd $HOME

# Create a temporary directory for the tests
mkdir -p Tmp_AMPS_test
cd Tmp_AMPS_test

# Remove the previous test directory if necessary
rm -rf AMPS

# Checkout and install the latest code
cvs co -D "`date +%m/%d/%Y` 0:30" AMPS 
cd AMPS
cvs co -D "`date +%m/%d/%Y` 0:30" AMPS_data 

# Install the AMPS
./Config.pl -install > test_amps.log

# Run test
make -j4 test  >> test_amps.log
          
# Create file with results' summary
ls -ltr  *diff > test_amps.res
echo '=============================================================='\
                >> test_amps.res
head -100 *diff >> test_amps.res

# Copy test_amps.log and test_amps.res to the server
 scp test_amps.res test_amps.log \
    dborovik@herot.engin.umich.edu:Sites/Current/`hostname -s`/ 


