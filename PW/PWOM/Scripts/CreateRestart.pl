#!/usr/bin/perl -s

my $Help = ($h or $help);
use strict;

&print_help if $Help;

my $ParamFile = $ARGV[0];
&print_help unless $ParamFile;

open(PARAMFILE, $ParamFile) or
    die "Could not open parameter file $ParamFile!\n";

my $InputFile1 = <PARAMFILE> 
    or die "Could not read first input file name from $ParamFile\n";
chop($InputFile1);

my $nInputFile = <PARAMFILE>
    or die "ould not read number of input files from $ParamFile\n";
chop($nInputFile);

# Read in latitude-longitude pairs
my $i;
my @LatLon;
for $i (0..$nInputFile-1){
    $LatLon[$i] = <PARAMFILE>
	or die "Could not read lat-lon line $i from $ParamFile\n";
}
close PARAMFILE;

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
for $iFile (0..($nInputFile-1)){
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
Usage: CreateRestart.pl PARAMFILE

where PARAMFILE contains the following lines

FirstInputFileName
NumberOfInputFiles
Lat Lon
Lat Lon
...

where the Latitude Longitude (degrees) pairs define the footpoint 
for the field lines.
"
}
