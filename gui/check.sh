#!/bin/sh

#
# Check the httpd.conf configuration file
#

`dirname $0`/apache/bin/apachectl configtest

