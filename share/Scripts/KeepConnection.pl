#!/usr/bin/perl

# This script can keep connection open on machines that would 
# otherwise close it. Simply start it on the X window as
#
# share/Scripts/KeepConnection.pl

$t = 0;
{
    print "time from start=$t\n";
    sleep 10;
    $t += 10;
    redo;
}
