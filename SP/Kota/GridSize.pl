#!/usr/bin/perl -s -i

# Put switches into properly named variables
my $Help         = $h or $H or $help or $Help;
my $GridSize     = $g;
my $ShowHuman    = $s;
my $ShowCompact  = not($g or $s);

use strict;

# Print help if required
&print_help if $Help;

# Set error string for error messages
my $ERROR = 'GridSize_ERROR:';

# Set the name of the file that contains the grid size
my $NameFile="src/param.h";

# Set grid size variables based on the arguments
my ($nRMax, $nMuMax, $nPMax);
if($GridSize=~/^\d+,\d+,\d+$/x){
    ($nRMax, $nMuMax, $nPMax)= split(',', $GridSize);
}elsif($GridSize){
    die "$ERROR -g=$GridSize should be 3".
	" integers separated with commas\n";
}

# Read previous grid size and set IsNewGrid to true if there is change
my $IsNewGrid;
&read_grid_size;

# Write the new grid into $NameFile if necessary
if($IsNewGrid){
    &set_grid_size;
    # Check if it worked be rereading info. Now IsNewGrid should be false.
    &read_grid_size;
    if($IsNewGrid){
	&show_grid_size;
	die "$ERROR incorrect grid information was saved";
    }
}

&show_grid_size if $ShowHuman;

print "GridSize.pl -g=$nRMax,$nMuMax,$nPMax\n" if $ShowCompact;

exit 0;

#############################################################################

sub read_grid_size{

    # Local variables to be read
    my ($nRMax_, $nMuMax_, $nPMax_);

    # Read size of the grid from $NameFile
    open(MODSIZE,$NameFile) or die "$ERROR could not open $NameFile\n";
    while(<MODSIZE>){
	next if /^\s*!/; # skip commented out lines
        $nRMax_=$1           if /\bnRMax\s*=\s*(\d+)/i;
	$nMuMax_=$1           if /\bnMuMax\s*=\s*(\d+)/i;
	$nPMax_=$1           if /\bnPMax\s*=\s*(\d+)/i;
    }
    close MODSIZE;

    die "$ERROR could not read nRMax from $NameFile\n" unless length($nRMax_);
    die "$ERROR could not read nMuMax from $NameFile\n" unless length($nMuMax_);
    die "$ERROR could not read nPMax from $NameFile\n" unless length($nPMax_);
 
    $nRMax           = $nRMax_            unless $nRMax;
    $nMuMax          = $nMuMax_           unless $nMuMax;
    $nPMax           = $nPMax_            unless $nPMax;
   
    $IsNewGrid = "$nRMax,$nMuMax,$nPMax" ne 
	"$nRMax_,$nMuMax_,$nPMax_";
}

#############################################################################

sub set_grid_size{

    print "Writing new grid size into $NameFile...\n";

    @ARGV = ($NameFile);

    while(<>){
	if(/^\s*!/){print; next} # Skip commented out lines
	s/\b(nRMax\s*=[^\d]*)(\d+)/$1$nRMax/i;
	s/\b(nMuMax\s*=[^\d]*)(\d+)/$1$nMuMax/i;
	s/\b(nPMax\s*=[^\d]*)(\d+)/$1$nPMax/i;
	print;
    }

}

##############################################################################

sub show_grid_size{

    print "nRMax=$nRMax, nMuMax=$nMuMax, nPMax=$nPMax\n";
}

#############################################################################

sub print_help{

    print "
Purpose:

    Get and set grid size information from and into src/param.h.
    If called without any switch the current settings are shown in
    a compact format which is suitable for processing by other codes.

Usage:

  GridSize.pl [-h] [-g=nR,nMu,nP] [-s]

  -h              print help message and stop

  -g=nR,nMu,nP    set the grid size to nR x nMu x nP

  -s              show current settings in human readable form

Examples:

  Show the current gridsize in compact form:

GridSize.pl

  Show the current grid size in human readable form:

GridSize.pl -s

  Set the grid size to 1500x60x100:

GridSize.pl -g=1500,60,100
\n\n";
    exit 0;
}
