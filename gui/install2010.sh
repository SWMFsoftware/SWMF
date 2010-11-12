#!/bin/sh
if [ "$1" = "" ] ; then
  echo 'Usage: '$0' port#'
  echo '  The port number should be > 1500 and preferably in the 8000s'
  echo 'Example:'
  echo '  '$0' 8080'
  exit
fi

# Set PORT
port=$1
rm -f PORT*
touch PORT=$port

# Set GUI HOMEDIR
hdir=`pwd`
echo '<?php  //Base directory of GUI install'  > htdocs/config.php
echo '  $HOMEDIR = "'${hdir}'/htdocs";'       >> htdocs/config.php
echo '?>'                                     >> htdocs/config.php


adir='apache'
echo 'Making directories ...'
if ls | grep -q $adir
then
  rm -rf $adir
fi
mkdir $adir
cd $adir

idir=`pwd`
echo 'Install directory = '$idir

mkdir src
cd src


#
# apache
#

Vapache='httpd-2.2.17'
echo '--> Downloading apache ...'
if ls ../../src_external | grep -q httpd
then
  cp ../../src_external/${Vapache}.tar.gz .
else
  wget http://apache.mirrors.versehost.com/httpd/${Vapache}.tar.gz
  cp -p ${Vapache}.tar.gz ../../src_external/ 
fi

echo '--> Unpacking apache ...'
echo '--> Unpacking apache ...' >> log_apache
gunzip < ${Vapache}.tar.gz | tar -xvf - >> log_apache

echo '--> Configuring apache ...'
echo '--> Configuring apache ...' >> log_apache
cd $Vapache
./configure --prefix=${idir} --with-port=${port} >> ../log_apache

echo '--> Compiling apache ...'
echo '--> Compiling apache ...' >> ../log_apache
make >> ../log_apache

echo '--> Installing apache ...'
echo '--> Installing apache ...' >> ../log_apache
make install >> ../log_apache
cd ..


#
# php
#

Vphp='php-5.3.3'
echo '--> Downloading php ...'
if ls ../../src_external | grep -q php
then
  cp ../../src_external/${Vphp}.tar.gz .
else
  wget http://www.php.net/distributions/${Vphp}.tar.gz
  cp -p ${Vphp}.tar.gz ../../src_external/
fi

echo '--> Unpacking php ...'
echo '--> Unpacking php ...' >> log_php
gunzip < ${Vphp}.tar.gz | tar -xvf - >> log_php

echo '--> Configuring php ...'
echo '--> Configuring php ...' >> log_php
cd $Vphp
./configure --prefix=${idir}/php --with-config-file-path=${idir}/php --with-apxs2=${idir}/bin/apxs >> ../log_php

echo '--> Compiling php ...'
echo '--> Compiling php ...' >> ../log_php
make >> ../log_php

echo '--> Installing php ...'
echo '--> Installing php ...' >> ../log_php
make install >> ../log_php

echo '--> Extras for php ...'
echo '--> Extras for php ...' >> ../log_php
#cp -p php.ini-recommended ${idir}/php/php.ini
cat php.ini-production | sed s/post_max_size\ =\ 8M/post_max_size\ =\ 500M/ | sed s/upload_max_filesize\ =\ 2M/upload_max_filesize\ =\ 500M/> ${idir}/php/php.ini
cd ..


# Add to httpd.conf file
cd ${idir}/conf
mv httpd.conf httpd.conf_ORIG
cp httpd.conf_ORIG httpd.conf
echo ''                                                              >> httpd.conf
echo '###'                                                           >> httpd.conf
echo '### Customizations'                                            >> httpd.conf
echo ''                                                              >> httpd.conf
echo '# Add index.php to your DirectoryIndex line:'                  >> httpd.conf
echo 'DirectoryIndex index.html index.php'                           >> httpd.conf
echo ''                                                              >> httpd.conf
echo 'AddType application/x-httpd-php .php'                          >> httpd.conf
echo ''                                                              >> httpd.conf
echo '# PHP Syntax Coloring'                                         >> httpd.conf
echo '# (optional but useful for reading PHP source for debugging):' >> httpd.conf
echo 'AddType application/x-httpd-php-source phps'                   >> httpd.conf
echo ''                                                              >> httpd.conf
cd ..

# Point to new documents directory
mv htdocs htdocs_ORIG
ln -s ../htdocs .

# Cleanup
cd ${idir}/src
rm -rf $Vapache
rm -rf $Vphp

