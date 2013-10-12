#!/bin/csh
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
cd $HOME/Sites

# Remove code from yesterday as well as various logs
rm -rf test.diff

# Create directory for new test results
setenv NEWTESTDIR SWMF_TEST_RESULTS/`date -v-1d +%Y/%m/%d`
mkdir -p ${NEWTESTDIR}

# Copy over test results but preserve the subdirectories in Current
cp -r Current/* ${NEWTESTDIR}/
rm -f Current/*/*

# Create new index.html and list changes in test.diff
./process_tests.pl
