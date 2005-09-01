#!/usr/bin/perl -s -i

# Put switches into properly named variables
my $Help         = $h or $H or $help or $Help;
my $GridSize     = $g;
my $ShowHuman    = $s;
my $ShowCompact  = not($g or $s);

use strict;

# Set the names of the files that contains the grid size
my $NameFile="src/ModSize.f90";

# Print help if required
&print_help if $Help;

# Set error string for error messages
my $ERROR = 'IE_GridSize_ERROR:';

# Set grid size variables based on the arguments
my ($nTheta, $nPsi);

if($GridSize=~/^\d+,\d+$/){
    ($nTheta, $nPsi)= split(',', $GridSize);
}elsif($GridSize){
    die "$ERROR -g=$GridSize should be 2 integers separated with a comma\n";
}

# Read previous grid size and set IsNewGrid to true if there is change
my $IsNewGrid;
&read_grid_size;

# Write the new grid into $NameFile if necessary
if($IsNewGrid){
    &set_grid_size;
    # Check if it worked by rereading info. Now IsNewGrid should be false.
    &read_grid_size;
    if($IsNewGrid){
	&show_grid_size;
	die "$ERROR incorrect grid information was saved";
    }
}

&show_grid_size if $ShowHuman;

print "GridSize.pl -g=$nTheta,$nPsi\n" if $ShowCompact;

exit 0;

#############################################################################

sub read_grid_size{

    # Local variables to be read
    my ($nTheta_, $nPsi_);

    # Read size of the grid from $NameFile
    open(FILE,$NameFile) or die "$ERROR could not open $NameFile\n";
    while(<FILE>){
	next if /^\s*!/; # skip commented out lines
	$nTheta_=$1 if /\bIONO_nTheta\s*=\s*(\d+)/i;
        $nPsi_  =$1 if /\bIONO_nPsi\s*=\s*(\d+)/i;
	last if $nTheta_ and $nPsi_;
    }
    close FILE;

    die "$ERROR could not read nTheta from $NameFile\n" 
	unless length($nTheta_);
    die "$ERROR could not read nPsi from $NameFile\n" 
	unless length($nPsi_);

    $nTheta      = $nTheta_      unless $nTheta;
    $nPsi        = $nPsi_        unless $nPsi;

    $IsNewGrid = "$nTheta,$nPsi" ne "$nTheta_,$nPsi_";
}

#############################################################################

sub set_grid_size{

    print "Writing new grid size into $NameFile...\n";

    @ARGV = ($NameFile);

    while(<>){
	if(/^\s*!/){print; next} # Skip commented out lines
	s/\b(IONO_nTheta\s*=[^\d]*)(\d+)/$1$nTheta/i;
	s/\b(IONO_nPsi\s*=[^\d]*)(\d+)/$1$nPsi/i;
	print;
    }

}

##############################################################################

sub show_grid_size{

    print "Number of colatitudes/hemisphere  : nTheta = $nTheta\n";
    print "Number of longitudes              : nPsi   = $nPsi\n";
}

#############################################################################

sub print_help{

    print "
Purpose:

    Get and set grid size information from and into 
    $NameFile and check if settings are correct.
    If called without any switch the current settings are shown in
    a compact format which is suitable for processing by other codes.

Usage:

    GridSize.pl [-h] [-g=nTheta,nPsi] [-s]

  -h              print help message and stop

  -g=nTheta,nPsi  set the grid size

  -s              show current settings in human readable form

Examples:

  Show current settings in a compact form:

GridSize.pl

  Show current settings in human readable form:

GridSize.pl -s

  Set grid size to nTheta=91, nPsi=361 (1 deg by 1 deg):

GridSize.pl -g=91,361

";
    exit 0;
}
