#!/bin/csh
#  Copyright (C) 2002 Regents of the University of Michigan, 
#  portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
cd $HOME/Sites

# Create directory for new tests
mkdir -p SWMF_TEST_RESULTS/`date +%Y/%m/%d`

# Remove code from yesterday as well as various logs
rm -rf SWMF_yesterday code.diff manual.log manual.err fortran.err

if(-d SWMF)then
    # Move the previous SWMF repository out of the way
    cd SWMF
    # uninstall now to remove the manuals created yesterday
    ./Config.pl -uninstall
    # There were created after clone during install and normally kept after uninstall
    rm -rf gitinfo.txt PC/ALTOR/srcBATL_orig UA/GITM2/srcData GM/BATSRUS/src/ModUser.f90.safe \
	    share/Python/swmfpy
    cd ..
    mv SWMF SWMF_yesterday
endif

# Get current version of SWMF+BATL and make it unreadable for the web
gitclone SWMF
cd SWMF; Config.pl -clone
gitclone BATL
cd BATL; Config.pl -clone
cd $HOME/Sites
chmod -R go-r SWMF

# Create diff file excluding .git directories
touch code.diff
if(-d SWMF_yesterday)then
    # Compare the current SWMF with yesterday's version
    diff -r -x .git SWMF SWMF_yesterday >& code.diff
endif

# Create manuals
cd ~/Sites/SWMF
Config.pl -install -compiler=gfortran >& ~/Sites/manual.log
make MANUAL >>& ~/Sites/manual.log

chmod -R go+r doc/ Copyrights/
chmod -R go-r doc/Tex
chmod -R go+r GM/BATSRUS/Doc/
chmod -R go-r GM/BATSRUS/Doc/Tex
chmod go+r PW/PWOM/doc/PWOM.pdf
chmod go+r share/IDL/doc/idl.pdf
chmod go+r util/CRASH/doc/*

cd ~/Sites/SWMF
make clean         >>& ~/Sites/manual.log

# Check if the manual creation worked
cd ~/Sites
grep -v logo_small.png manual.log | grep -C10 Error > manual.err

# Check Fortran formatting
cd ~/Sites/SWMF
make FORMATF90 > ~/Sites/fortran.err
# return to original code (we could also commit it!)
gitall checkout .
