#!/usr/bin/perl -s

my $Help    = ($h or $help);
my $Debug   = ($d or $debug);
my $Verbose = ($v or $verbose);
my $Append  = ($M);
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
my $nLine = @plot_file;
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
my @Head;
SNAPSHOT:{

    my $fh;
    $iPlot++;
    $iPlotString++;
    #print "snapshot $iPlot\n" if $Verbose;
    my $variables;
    # Read header and discard
    my $nStep;
    my $Time;
    foreach $plot (@plot_file) {
	$fh = $plot_handle{$plot};
	for $iLine (1..$nHeader){
	    my $line;
	    $line = <$fh>;
	    print "header $iLine: $line\n" if $Debug;
	    $Head[$iLine]=$line;
	    #before discarding get nData and variable list
	    if ($iLine eq 2) {
		my @linesplit= split/\s+/, $line ;
		$nStep = $linesplit[1];
		$Time  = $linesplit[2];
	    }
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
	    if($Append){
		open OUT, ">>plot_movie_iAlt${iAlt}.dat";
	    }else{
		open OUT, ">plot_snapshot${iPlotString}_iAlt${iAlt}.dat";
	    }
	}
	print "Save into plot_snapshot${iPlot}_iAlt${iAlt}.dat\n" if $Verbose;
	
	# write Header into file
	$variables=~s/Lat/X/;
	$variables=~s/r/Y/;
	$variables=~s/Lon/Z/;
	
	my $header = "$variables";
	if($iAlt eq $AltSlice){
	    print OUT "PWOM output_hd23\n";
	    print OUT "$nStep  $Time  -2  1 19\n";
	    print OUT "$nLine 1\n";
	    print OUT $Head[4];
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
	    my $y = -1.0*$r*sin($theta)*sin($phi); #mult y by -1 for later rotation
	    my $z = $r*cos($theta);
	    
	    $x=sprintf("%.10f",$x);
	    $y=sprintf("%.10f",$y);
	    $z=sprintf("%.10f",$z);
	    
	    $x=$x.'E+00';
	    $y=$y.'E+00';
	    $z=$z.'E+00';

	    my $nLineElements=scalar@linesplit;
	    #Interchange x and y so as to rotate resulting idl plot to have 
	    #noon point up.
	    $linesplit[2]=$x;
	    $linesplit[1]=$y;
	    $linesplit[3]=$z;

	    # remove extra "-" in front of zeros
	    my $iLineElement;
	    for $iLineElement (1..$nLineElements-2){
		if($linesplit[$iLineElement] == 0){
		    $linesplit[$iLineElement]=~ s/-//;
		}
	    }

	    if ($linesplit[1]>=0){
		$line="  ";
	    }else{
		$line=" ";
	    }
	    
	    #reform the line with XYZ inplace of r lat lon
	    my $iLineElement;
	    for $iLineElement (1..$nLineElements-2){
		if ($linesplit[$iLineElement+1] >= 0){
		    $line=$line . $linesplit[$iLineElement] . "  ";
		}else{
		    $line=$line . $linesplit[$iLineElement] . " ";
		}
	    }
	    $line=$line . $linesplit[$nLineElements-1] . "\n";
	    
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
Purpose: Extract 2D slices in altitude for ploting with idl. This routine 
         reads all of the individual field-line output and extracts the 
         values at a given altitude. It is assumed that the individual lines 
         are in the same directory as this script. Examples of use are given 
         in the PWOM manual.

Usage: PostProcess.pl [-v] [-d] [-h]

-h -help    - print help message and exit
-v -verbose - print verbose information
-d -debug   - print debugging information
-M          - Create a single file containing all the snapshots at a 
              given altitude for animation in idl. If this flag is not 
              specified, then output a different file for each snapshot.


";
    exit;
}
