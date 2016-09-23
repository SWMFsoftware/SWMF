#!/bin/csh
# This sends the results of AMPS nightly tests to the server:
#   herot.engin.umich.edu 
# for further processing
# This script is meant to be used by the AMPS developers.

# The script can be executed by simply typing
#
# ./scp_run_amps.sh
#
# To scp the AMPS tests' resuults perodically, you can use the crontab facility.
# Type 'crontab -e' to add a new entry to your crontab.
# Here is an example entry for for scp'ing nightly runs results at 8:00 am:
#
# 30 8 * * * $HOME/bin/scp_test_amps.sh

# set the working directory
set WorkDir = $HOME  
set Server  = dborovik@herot.engin.umich.edu

#>Pleiades ############################
#set WorkDir = /nobackup/`whoami`    <#

#>Stampede ############################
#set WorkDir = $WORK                 <#

#>Yellowstone ###########
#set WorkDir =         <#

# Go to your home directory
cd $WorkDir/Tmp_AMPS_test

# Create file with results' summary and scp to server

# GNU compiled test
#>GNUAll #####################################################################
#cd GNU/AMPS                                                                 #
#rm -rf test_amps.res
#ls -ltr  *diff > test_amps.res                                              #
#echo '===============================================' >> test_amps.res     # 
#head -100 *diff >> test_amps.res                                            #
#scp test_amps.res test_amps.log ${Server}:Sites/Current/`hostname -s`_gnu/  #
#cd ../..                                                                   <#


# Intel compiled test
#>IntelAll ###################################################################
#cd Intel/AMPS                                                               #
#rm -rf test_amps.res
#ls -ltr  *diff > test_amps.res                                              #
#echo '===============================================' >> test_amps.res     # 
#head -100 *diff >> test_amps.res                                            #
#scp test_amps.res test_amps.log ${Server}:Sites/Current/`hostname -s`_intel/#
#cd ../..                                                                   <#

# PGI compiled test
#>PGIAll #####################################################################
#cd PGI/AMPS                                                                 #
#rm -rf test_amps.res
#ls -ltr  *diff > test_amps.res                                              #
#echo '===============================================' >> test_amps.res     # 
#head -100 *diff >> test_amps.res                                            #
#scp test_amps.res test_amps.log ${Server}:Sites/Current/`hostname -s`_pgi/  #
#cd ../..                                                                   <#

