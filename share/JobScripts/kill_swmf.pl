#!/usr/bin/perl
#^CFG COPYRIGHT UM

open(TEST, "ps xw |");

while (<TEST>) {
    if (/SWMF.exe/ or /mpirun/) {
	print "Killing processes : \n";
	print $_;
	split;
	system "kill -9 ".$_[0];
    }
}

close(TEST);


