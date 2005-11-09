#!/usr/bin/perl -s
# Preplot multiple files

my $Help = $h;
my $Keep = $k;
my $Gzip = $g;
use strict;

if($Help or not @ARGV){
    print "
Purpose:

   Preplot (and compress) multiple Tecplot data files.
   If the data file is compressed, use gunzip before preplot.

Usage: 

   preplot.pl [-h] [-k] [-g] FILE1 [FILE2 FILE3 ...] 

   -h    Print this help message
   -k    Keep original files
   -g    Compress .plt files

Examples:

   Preplot .dat files and keep the originals:

preplot.pl -k *.dat

   Preplot compressed .dat files and compress results too:

preplot.pl -g *.dat.gz

";
    exit;
}
my $file;
foreach $file (@ARGV){
    my $datfile = $file;
    `gunzip -c $file > $datfile` if $datfile =~ s/\.dat\.gz$/.dat/;

    if($datfile !~ /\.dat$/){
	warn "WARNING in preplot.pl: extension should be .dat or .dat.gz: ".
	    "$file\n";
	next;
    }
    `preplot $datfile`;
    if($Gzip){
	my $pltfile = $datfile;
	$pltfile =~ s/.dat/.plt/;
	`gzip $pltfile`;
    }
    # Remove uncompressed data file if any
    unlink $datfile if $file =~ /\.dat\.gz$/;
    # Remove original file
    unlink $file unless $Keep;
}
