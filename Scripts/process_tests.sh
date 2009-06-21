#!/bin/csh
cd $HOME/Sites

# Remove test results and code from yesterday
rm -rf SWMF_yesterday test.diff code.diff manual.log manual.err

# Move the previous SWMF repository out of the way
cd SWMF
# uninstall now to remove the manuals created yesterday
./Config.pl -uninstall
cd ..
mv SWMF SWMF_yesterday

# Create directory for new test results
setenv NEWTESTDIR SWMF_TEST_RESULTS/`date -v-1d +%Y/%m/%d`
mkdir -p ${NEWTESTDIR}

# Copy over test results but preserve the subdirectories in Current
cp -r Current/* ${NEWTESTDIR}/
rm -f Current/*/*

# Get current version of SWMF and make it unreadable for the web
cvs co -D "19:00EDT" SWMF
chmod -R go-r SWMF

# Compare the current SWMF with yesterday's version
rm -rf SWMF_yesterday/GM/BATSRUS/src/ModUser.f90.safe
diff -r -x 'CVS' SWMF SWMF_yesterday >& code.diff

# Create new index.html and list changes in test.diff
./process_tests.pl

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

