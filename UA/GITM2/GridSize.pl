#!/usr/bin/perl -s -i

# Put switches into properly named variables
my $Help         = $h or $H or $help or $Help;
my $nCells       = $c;
my $MaxBlock     = $b;
my $GridSize     = $g;
my $ShowHuman    = $s;
my $ShowCompact  = not($c or $b or $g or $s);

use strict;

# Set the names of the files that contains the grid size
my $NameFile="src/ModSize.f90";

# Print help if required
&print_help if $Help;

# Set error string for error messages
my $ERROR = 'GridSize_ERROR:';

# Check parameters
die "$ERROR Do not use -g together with -b or -c\n"
    if $GridSize and ($MaxBlock or $nCells);

# Set grid size variables based on the arguments
my ($nLon, $nLat, $nAlt);

if($GridSize=~/^\d+,\d+,\d+,\d+$/){
    ($nLon, $nLat, $nAlt, $MaxBlock)= split(',', $GridSize);
}elsif($GridSize){
    die "$ERROR -g=$GridSize should be 4 integers separated with commas\n";
}

if($nCells=~/^\d+,\d+,\d+/){
    ($nLon, $nLat, $nAlt)= split(',',$nCells);
}elsif($nCells){
    die "$ERROR -c=$nCells should be 3 integers separated with commas\n";
}

if($MaxBlock and $MaxBlock !~ /^\d+$/){
    die "$ERROR -b=$MaxBlock should be a positive integer\n";
}

# Read previous grid size and set IsNewGrid to true if there is change
my $IsNewGrid;
&read_grid_size;

# Check the grid size (to be set)
die "$ERROR MaxBlock=$MaxBlock must be positive integers\n" 
    if $MaxBlock<1 or $MaxBlock != int($MaxBlock);

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

print "GridSize.pl -g=$nLon,$nLat,$nAlt,$MaxBlock\n" if $ShowCompact;

exit 0;

#############################################################################

sub read_grid_size{

    # Local variables to be read
    my ($nLon_, $nLat_, $nAlt_, $MaxBlock_);

    # Read size of the grid from $NameFile
    open(FILE,$NameFile) or die "$ERROR could not open $NameFile\n";
    while(<FILE>){
	next if /^\s*!/; # skip commented out lines
        $nLon_=$1         if /\bnLons\s*=\s*(\d+)/i;
	$nLat_=$1         if /\bnLats\s*=\s*(\d+)/i;
	$nAlt_=$1         if /\bnAlts\s*=\s*(\d+)/i;
	$MaxBlock_=$1     if /\bnBlocksMax\s*=\s*(\d+)/i;
	last if $nLon_ and $nLat_ and $nAlt_ and $MaxBlock_;
    }
    close FILE;

    die "$ERROR could not read nLon from $NameFile\n" unless length($nLon_);
    die "$ERROR could not read nLat from $NameFile\n" unless length($nLat_);
    die "$ERROR could not read nAlt from $NameFile\n" unless length($nAlt_);
    die "$ERROR could not read MaxBlock from $NameFile\n" 
	unless length($MaxBlock_);

    $nLon        = $nLon_        unless $nLon;
    $nLat        = $nLat_        unless $nLat;
    $nAlt        = $nAlt_        unless $nAlt;
    $MaxBlock    = $MaxBlock_    unless $MaxBlock;

    $IsNewGrid = "$nLon,$nLat,$nAlt,$MaxBlock" ne 
	"$nLon_,$nLat_,$nAlt_,$MaxBlock_";
}

#############################################################################

sub set_grid_size{

    print "Writing new grid size into $NameFile...\n";

    @ARGV = ($NameFile);

    while(<>){
	if(/^\s*!/){print; next} # Skip commented out lines
	s/\b(nLons\s*=[^\d]*)(\d+)/$1$nLon/i;
	s/\b(nLats\s*=[^\d]*)(\d+)/$1$nLat/i;
	s/\b(nAlts\s*=[^\d]*)(\d+)/$1$nAlt/i;
	s/\b(nBlocksMax\s*=[^\d]*)(\d+)/$1$MaxBlock/i;
	print;
    }

}

##############################################################################

sub show_grid_size{

    print "Number of cells in a block: nLon=$nLon, nLat=$nLat, nAlt=$nAlt\n";
    print "Max. number of blocks     : MaxBlock=$MaxBlock\n";
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

    GridSize.pl [-h] [-c=nLon,nLat,nAlt] [-b=MaxBlock] 
                [-g=nLon,nLat,nAlt,MaxBlock] [-s]

  -h              print help message and stop

  -c=nLon,nLat,nAlt 
                  set the block size to be nLon*nLat*nAlt

  -b=MaxBlock
                  set the maximum number of blocks to MaxBlock

  -g=nLon,nLat... set the block size and the maximum number of blocks together

  -s              show current settings in human readable form

Examples:

  Show current settings in a compact form:

GridSize.pl

  Set block size to 18x18x25 and the number of blocks to 4
  and show current settings in a human readable form:

GridSize.pl -c=18,18,25 -b=4 -s

  Set blocks size to 36x36x50 and the number of blocks to 16:

GridSize.pl -g=36,36,50,16


";
    exit 0;
}
