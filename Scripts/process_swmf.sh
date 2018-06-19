#!/bin/csh
#  Copyright (C) 2002 Regents of the University of Michigan, 
#  portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
cd $HOME/Sites

# Create directory for new tests
mkdir -p SWMF_TEST_RESULTS/`date +%Y/%m/%d`

# Remove code from yesterday as well as various logs
rm -rf SWMF_yesterday code.diff manual.log manual.err

if(-d SWMF)then
    # Move the previous SWMF repository out of the way
    cd SWMF
    # uninstall now to remove the manuals created yesterday
    ./Config.pl -uninstall
    cd ..
    mv SWMF SWMF_yesterday
endif

# Get current version of SWMF and make it unreadable for the web
gitclone SWMF
cd SWMF; Config.pl -install -compiler=gfortran; Config.pl -uninstall
cd $HOME/Sites
chmod -R go-r SWMF

touch code.diff
if(-d SWMF_yesterday)then
    # Compare the current SWMF with yesterday's version
    diff -r SWMF SWMF_yesterday >& code.diff
endif

# Create manuals
cd ~/Sites/SWMF
Config.pl -install -compiler=gfortran >& ~/Sites/manual.log
make PDF           >>& ~/Sites/manual.log
#make HTML          >>& ~/Sites/manual.log
chmod -R go+r doc/ Copyrights/
chmod -R go-r doc/Tex

# Commented out due to LATEX style/class files not found
cd ~/Sites/SWMF/util/CRASH/doc/Tex
make PDF           >>& ~/Sites/manual.log
chmod go+r ../*.pdf ../index.html ../RELEASENOTES

cd ~/Sites/SWMF/GM/BATSRUS
make PDF           >>& ~/Sites/manual.log
#make HTML          >>& ~/Sites/manual.log
chmod -R go+r Doc/
chmod -R go-r Doc/Tex

cd ~/Sites/SWMF/PW/PWOM
make PDF           >>& ~/Sites/manual.log
chmod go+r doc/PWOM.pdf

cd ~/Sites/SWMF
make clean         >>& ~/Sites/manual.log

# Check if the manual creation worked
cd ~/Sites
grep -v logo_small.png manual.log | grep -C10 Error > manual.err

