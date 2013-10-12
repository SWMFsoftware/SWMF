#!/usr/bin/perl -s
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

if (!$type) { $type = $ENV{'HOSTTYPE'} }
if (!$host) { $host = $ENV{'HOST'} }

if ($type =~ /powermac/) { 
    $compiler = ''; 
    $mpiversion = '';
}

if ($type =~ /linux/ && $host =~ /cfe/) { 
    $compiler = 'ifort'; 
    $mpiversion = 'Altix'; 
}


print $type,"\n" if $t;
print $host,"\n" if $h;
print $compiler,"\n" if $c;
print $mpiversion,"\n" if $m;

exit(1);



