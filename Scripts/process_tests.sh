#!/bin/csh
cd $HOME/Sites

# Remove test results and code from yesterday
rm -rf 7days SWMF_yesterday test.diff code.diff manual.log manual.err

# Shift Today to Yesterday and Current to Today
cd SWMF
./Config.pl -uninstall
cd ..
mv SWMF SWMF_yesterday
mv 6days 7days
mv 5days 6days
mv 4days 5days
mv 3days 4days
mv 2days 3days
mv Yesterday 2days
mv Today Yesterday
mv Current Today

# Recreat Current directory for new results
mkdir Current
cd Current
mkdir grendel grendel_ompi grid mesh nyx nyx_pgf90 columbia xena xena_xlf
cd ..

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

## Send a notification if the test results have changed
#if(-s test.diff)then
#    /usr/bin/mail -s 'SWMF test results have changed!' gtoth@umich.edu <<EOF
#    Test results have changed relative to the previous day: see listing below.
#    Check the table at http://herot.engin.umich.edu/~gtoth/ 
#    for more information. Capitalized items indicate changes.
#
#`cat test.diff`
#
#First 20 lines of changes in the SWMF code (if any):
#`head -20 code.diff`
#
#EOF
#endif
