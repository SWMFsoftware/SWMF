#!/bin/csh
cd $HOME/Sites

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
cvs co SWMF
chmod -R go-r SWMF

touch code.diff

if(-d SWMF_yesterday)then
    # Compare the current SWMF with yesterday's version
    rm -rf SWMF_yesterday/GM/BATSRUS/src/ModUser.f90.safe
    diff -r -x 'CVS' SWMF SWMF_yesterday >& code.diff
endif

# Create manuals
cd ~/Sites/SWMF
Config.pl -install  >& ~/Sites/manual.log
make PDF           >>& ~/Sites/manual.log
make HTML          >>& ~/Sites/manual.log
chmod -R go+r doc/ Copyrights/
chmod -R go-r doc/Tex

cd ~/Sites/SWMF/util/CRASH/doc/Tex
make PDF           >>& ~/Sites/manual.log
chmod go+r ../*.pdf ../index.html

cd ~/Sites/SWMF/GM/BATSRUS
make PDF           >>& ~/Sites/manual.log
make HTML          >>& ~/Sites/manual.log
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

