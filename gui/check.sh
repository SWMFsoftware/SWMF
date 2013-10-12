#!/bin/sh
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

#
# Check the httpd.conf configuration file
#

`dirname $0`/apache/bin/apachectl configtest

