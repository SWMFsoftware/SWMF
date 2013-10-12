#!/bin/sh
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

#
# Make and install the current SWMF Manual
#

echo ''
echo 'Making SWMF and REFERENCE manuals ...'
echo '  ... this may take a minutes or two (output redirected) ...'
echo ''
mkdir tmp
rsync -avq --exclude gui ../. tmp/
cd tmp
./Config.pl -install >& log1
make HTML >& log2
if ls ../htdocs | grep -q Manuals
then
  rm -rf ../htdocs/Manuals
fi
mkdir ../htdocs/Manuals
mv doc/HTML/SWMF ../htdocs/Manuals/
mv doc/HTML/REFERENCE ../htdocs/Manuals/
cd ..
rm -rf tmp
echo ''
echo 'Finished.  SWMF and REFERENCE manuals created and installed.'
echo ''
