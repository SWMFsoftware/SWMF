#!/usr/bin/perl -s

my $Help = ($h or $help);
use strict;

&print_help if $Help;

my $nCircle = $ARGV[0];

if($nCircle < 6){ $nCircle=6}

my$InputFile1 = 'restart_iline0001.dat';



# Create latitude-longitude pairs
my $nLine= $nCircle * 6.0 * (1.0+$nCircle)/2.0;
my $dLat = 40.0/$nCircle;
my $Lat = 90.0;
my $Lon = 0.0;
my $i = 0;
my @LatLon;
for my $iCircle (1..$nCircle){
    $Lat=$Lat-$dLat;
    my $dLon = 360.0/($iCircle * 6.0);
    for my $iPoint (0..($iCircle * 6.0-1.0)){
	$Lon = $iPoint*$dLon;
	$LatLon[$i] = "$Lat $Lon\n";
	$i++;
    }
}

open (FILE, $InputFile1)
    or die "Could not open first input file $InputFile1\n";

print "Reading in first input file $InputFile1 ...\n";
my $Line;
my $iLine;
my @Lines;
while ($Line = <FILE>){
    $Lines[$iLine++] = $Line;
}
close FILE;
print "Read $iLine lines from first input file $InputFile1\n";

print "Writing out input files ...\n";

my $iFile;
my $File;
for $iFile (0..($nLine-1)){
    my $iFileString = "0" x (4 - length($iFile+1)) . ($iFile + 1);
    $File = $InputFile1;
    $File =~ s/0001/$iFileString/
	or die "First input file name $InputFile1 does not contain 0001!\n";

    print "File = $File\n";

    open(FILE, ">$File")
	or die "Could not open file $File for writing!\n";
    $Lines[1] = $LatLon[$iFile];
    print FILE @Lines;
    close FILE;
}

exit;

#############################################################################

sub print_help{
    print "
Usage: CreateRestart.pl nCircle

Where nCircle is the number of circles of points filling the space between 
90 and 50 degrees latitude. nCircle default and minimum is 6. 
The points are chosen to ensure uniform layout. The total number of lines 
can be determined from nCircles according to:

3*nCircle(1+nCircle)
 
The script must be run in a restart directory with a seed line 
restart_iline0001.out

where the Latitude Longitude (degrees) pairs define the footpoint 
for the field lines.
"
}
