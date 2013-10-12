#!/bin/sh
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
if [ "$1" = "" ] ; then
  echo 'Usage: '$0' port#'
  echo '  The port number should be > 1500 and preferably in the 8000s'
  echo 'Example:'
  echo '  '$0' 8080'
  exit
fi

port=$1
rm -f PORT*
touch PORT=$port

sdir=`pwd`
echo 'Starting directory = '$sdir
echo 'Port # = '$port

echo 'Making directories ...'
if ls | grep -q apache
then
  rm -rf apache
fi
mkdir apache
mkdir apache/src
cd apache/src

#
# apache
#

Vapache='httpd-2.2.3'
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
./configure --prefix=${sdir}/apache --with-port=${port} --enable-so >> ../log_apache

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

##php5
#Vphp='php-5.1.4'
##php4
Vphp='php-4.4.4'
##
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
./configure --with-apxs2=${sdir}/apache/bin/apxs --prefix=${sdir}/apache/php --with-config-file-path=${sdir}/apache/php --disable-cgi --with-zlib  >> ../log_php

echo '--> Compiling php ...'
echo '--> Compiling php ...' >> ../log_php
make >> ../log_php

echo '--> Installing php ...'
echo '--> Installing php ...' >> ../log_php
make install >> ../log_php

echo '--> Extras for php ...'
echo '--> Extras for php ...' >> ../log_php
##php5
#cp -p libs/libphp5.so ${sdir}/apache/modules/
##php4
##
#cp -p php.ini-recommended ${sdir}/apache/php/php.ini
cat php.ini-recommended | sed s/post_max_size\ =\ 8M/post_max_size\ =\ 500M/ | sed s/upload_max_filesize\ =\ 2M/upload_max_filesize\ =\ 500M/> ${sdir}/apache/php/php.ini
cd ${sdir}/apache/conf
mv httpd.conf httpd.conf_ORIG
cp httpd.conf_ORIG httpd.conf
echo ''                                                              >> httpd.conf
echo '###'                                                           >> httpd.conf
echo '### Customizations'                                            >> httpd.conf
echo ''                                                              >> httpd.conf
echo '<Directory "'${sdir}'/apache/htdocs/*_*">'                     >> httpd.conf
echo '    Options Indexes FollowSymLinks'                            >> httpd.conf
echo '    AllowOverride None'                                        >> httpd.conf
echo '    Order allow,deny'                                          >> httpd.conf
echo '    Allow from all'                                            >> httpd.conf
echo ''                                                              >> httpd.conf
echo '    AuthType Basic'                                            >> httpd.conf
echo '    AuthName "GUI Restricted"'                                 >> httpd.conf
echo '    AuthUserFile '${sdir}'/security/passwords'                 >> httpd.conf
echo '    Require valid-user'                                        >> httpd.conf
echo '</Directory>'                                                  >> httpd.conf
echo ''                                                              >> httpd.conf
##php5
#echo '# Load PHP 5.x:'                                               >> httpd.conf
#echo 'LoadModule php5_module        modules/libphp5.so'              >> httpd.conf
##php4
echo '# Load PHP 4.x:'                                               >> httpd.conf
echo 'LoadModule php4_module        modules/libphp4.so'              >> httpd.conf
##
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


#
# Cleanup
#
rm -rf src/$Vapache
rm -rf src/$Vphp

mv htdocs htdocs_ORIG
ln -s ../htdocs .

cd ..


#
# Set default password
#

echo 'Setting default site user/password to gui/go'
touch security/passwords
apache/bin/htpasswd -bm security/passwords gui go


#
# Set GUI HOMEDIR
#

echo '<?php  //Base directory of GUI install'  > htdocs/config.php
echo '  $HOMEDIR = "'${sdir}'/htdocs";'       >> htdocs/config.php
echo '?>'                                     >> htdocs/config.php
