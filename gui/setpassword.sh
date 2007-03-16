#!/bin/sh

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

