#!/bin/csh
#  Copyright (C) 2002 Regents of the University of Michigan, 
#  portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
cd $HOME/Sites

# Remove code from yesterday as well as various logs
rm -rf SWMF_yesterday code.diff manual.log manual.err fortran.err \
    name.err process_swmf.log

# Create directory for new tests
mkdir -p SWMF_TEST_RESULTS/`date +%Y/%m/%d` >& ~/Sites/process_swmf.log

if(-d SWMF)then
    # Move the previous SWMF repository out of the way
    cd SWMF
    # uninstall now to remove the manuals created yesterday
    ./Config.pl -uninstall >>& ~/Sites/process_swmf.log
    # These were created after clone during install and normally kept after uninstall
    rm -rf gitinfo.txt PC/ALTOR/srcBATL_orig UA/GITM2/srcData GM/BATSRUS/src/ModUser.f90.safe
    cd ..
    mv SWMF SWMF_yesterday
endif

# Get current version of SWMF+BATL and make it unreadable for the web
gitclone -s SWMF >>& ~/Sites/process_swmf.log
cd SWMF; Config.pl -clone >>& ~/Sites/process_swmf.log
gitclone -s BATL >>& ~/Sites/process_swmf.log
cd BATL; Config.pl -clone >>& ~/Sites/process_swmf.log
cd $HOME/Sites
chmod -R go-r SWMF >>& ~/Sites/process_swmf.log

# Create diff file excluding .git directories
if(-d SWMF_yesterday)then
    # Compare the current SWMF with yesterday's version
    diff -r -x .git SWMF SWMF_yesterday >& code.diff
endif
touch code.diff

# Create manuals
cd ~/Sites/SWMF
Config.pl -install -compiler=gfortran >& ~/Sites/manual.log
rm -rf UA/GITM/share UA/GITM/util ### clean does not work correctly
share/Scripts/gitall pull --depth 100
make MANUAL >>& ~/Sites/manual.log

chmod -R go+r doc/ Copyrights/
chmod -R go-r doc/Tex
chmod -R go+r GM/BATSRUS/Doc/
chmod -R go-r GM/BATSRUS/Doc/Tex
chmod go+r PW/PWOM/doc/PWOM.pdf
chmod go+r share/IDL/doc/idl.pdf
chmod go+r util/CRASH/doc/*

cd ~/Sites/SWMF
make clean >>& ~/Sites/manual.log

# Check if the manual creation worked
cd ~/Sites
grep -v logo_small.png manual.log | grep -C10 Error > manual.err

# Check Fortran formatting
cd ~/Sites/SWMF
make FORMATF90 > ~/Sites/fortran.err
# commit and push formatted code
gitall "commit -m FormatFortran .; git push" >>& ~/Sites/process_swmf.log

# Check data names
cd ~/Sites/SWMF/GM/BATSRUS
~/Sites/SWMF/share/Scripts/CheckDataName.pl src/*.f90 srcEquation/*.f90 \
    srcBATL/*.f90 srcInterface/*.f90 >& ~/Sites/name.err

