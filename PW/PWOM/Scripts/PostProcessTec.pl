#!/usr/bin/perl -s
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

my $Help = ($h or $help);
my $Debug = ($d or $debug);
my $Verbose = ($v or $verbose);

use strict;
use FileHandle;
&print_help if $Help;

my $plot;
my $file;
my %plot_handle;
my @test;
my $i=1;
my $iLine;
my $nHeader=5;
my $pi=3.14159;
my $rPlanet=1.0;

print STDOUT "Enter Alt Slice For Output: ";
my $AltSlice = <STDIN>;
chop($AltSlice);

my @plot_file = sort glob "plots_iline????.out";

for $file (@plot_file) {
    my $fh;
    open($fh,"<",$file) or next;
    my $a=$fh;
    $plot_handle{$file} = $a;
}

my $nfile = @plot_file;

my $nData; # number of grid points in radial direction



my $iPlot;
my $iPlotString = "0000";
SNAPSHOT:{

    my $fh;
    $iPlot++;
    $iPlotString++;
    #print "snapshot $iPlot\n" if $Verbose;
    my $variables;
    # Read header and discard
    foreach $plot (@plot_file) {
	$fh = $plot_handle{$plot};
	for $iLine (1..$nHeader){
	    my $line;
	    $line = <$fh>;
	    print "header $iLine: $line\n" if $Debug;
	    #before discarding get nData and variable list
	    if ($iLine eq 5) {
		$variables =$line;
	    }elsif ($iLine eq 3){
		$line =~ /\s*(\d+)/ 
		    or die "Could not match nData in $line";
		$nData = $1;
	    }
	}
    }
    
    # loop over altitude and write out each fieldlines alt slice
    my $iAlt;
    for $iAlt (1..$nData){
		    
	# open a file for output
	if($iAlt eq $AltSlice){
	    open OUT, ">plot_snapshot${iPlotString}_iAlt${iAlt}.dat";
	}
	print "Save into plot_snapshot${iPlot}_iAlt${iAlt}.dat\n" if $Verbose;
	
	# write Tecplot header into file
	$variables=~s/\s+g//;
	$variables=~s/r/X/;
	$variables=~s/Lat/Y/;
	$variables=~s/Lon/Z/;
	$variables=~s/\s+(\w+)/ "$1"/g;
	
	
	my $header = "Variables = $variables";
	if($iAlt eq $AltSlice){
	    print OUT $header;
	}
	# write data into files
	foreach $plot (@plot_file) {
	    print "file = $plot\n" if $Verbose;
	    
	    # open get filehandle from hash for reading
	    $fh = $plot_handle{$plot};
	    print "fh = $fh\n" if $Debug;
	    
	    # read line from current alt
	    my $line = <$fh>;
	    
	    # parse line and convert the coords to XYZ
	    my @linesplit= split/\s+/, $line ;
	    my $r = $rPlanet;
	    my $theta = ( 90.0 - $linesplit[2] ) * $pi/180.0;
	    my $phi = $linesplit[3] * $pi/180.0;
	    my $x = $r*sin($theta)*cos($phi);
	    my $y = $r*sin($theta)*sin($phi);
	    my $z = $r*cos($theta);
	    
	    my $nLineElements=scalar@linesplit;
	    $linesplit[1]=$x;
	    $linesplit[2]=$y;
	    $linesplit[3]=$z;
	    $line=" ";
	    
	    #reform the line with XYZ inplace of r lat lon
	    my $iLineElement;
	    for $iLineElement (1..$nLineElements-1){
		$line=$line . $linesplit[$iLineElement] . " ";
	    }
	    $line=$line."\n";
	    
	    # write out line
	    if($iAlt eq $AltSlice){
		print OUT $line;
	    }
	}
	
    }
    redo SNAPSHOT unless eof($fh);
}


exit;
##############################################################################

sub print_help{
    print "
Purpose: extract an altitude slices for ploting in tecplot

Usage: PostProcess.pl [-v] [-d] [-h]

-h -help    - print help message and exit
-v -verbose - print verbose information
-d -debug   - print debugging information

";
    exit;
}
