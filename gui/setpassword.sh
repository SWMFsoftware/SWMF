#!/bin/sh
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

#
# Set the web password
#

DIR=`dirname $0`
echo $DIR
if [ "$DIR" != "." ] ; then
  echo 'ERROR: Needs to be run from gui directory.'
  exit
fi

echo 'Resetting the password for user: gui'
apache/bin/htpasswd -m security/passwords gui

