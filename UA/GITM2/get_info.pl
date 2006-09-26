#!/usr/bin/perl -s

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



