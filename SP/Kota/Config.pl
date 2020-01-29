#!/usr/bin/perl

# Allow in-place editing                                                        
$^I = "";

# Add local directory to search                                                 
push @INC, ".";

use strict;

our @Arguments       = @ARGV;
my $config     = "share/Scripts/Config.pl";
require "../../$config";

# These are inherited from $config
our %Remaining;   # Arguments not handled by share/Scripts/Config.pl
our $Show;
our $Help;
our $ERROR;
our $WARNING;
our $NewGridSize;
our $ShowGridSize;

&print_help if $Help;

# Grid size variables
my $NameFile="src/param.h";
my $GridSize;
my ($nRMax, $nMuMax, $nPMax);

&get_settings;

foreach (@Arguments){
    if(/^-s$/)                {$Show=1;  next};
    warn "WARNING: Unknown flag $_\n" if $Remaining{$_};
}

&set_grid_size if $NewGridSize and $NewGridSize ne $GridSize;

&show_settings if $Show;

print "Config.pl -g=$GridSize\n" if $ShowGridSize and not $Show;

exit 0;

#############################################################################

sub get_settings{

    # Read size of the grid from $NameFile
    open(MODSIZE,$NameFile) or die "$ERROR could not open $NameFile\n";
    while(<MODSIZE>){
	next if /^\s*!/; # skip commented out lines
        $nRMax=$1           if /\bnRMax\s*=\s*(\d+)/i;
	$nMuMax=$1           if /\bnMuMax\s*=\s*(\d+)/i;
	$nPMax=$1           if /\bnPMax\s*=\s*(\d+)/i;
    }
    close MODSIZE;

    die "$ERROR could not read nRMax from $NameFile\n" unless length($nRMax);
    die "$ERROR could not read nMuMax from $NameFile\n" unless length($nMuMax);
    die "$ERROR could not read nPMax from $NameFile\n" unless length($nPMax);
 
    $GridSize = "$nRMax,$nMuMax,$nPMax";
}

#############################################################################

sub set_grid_size{

    $GridSize = $NewGridSize;

    if($GridSize=~/^\d+,\d+,\d+$/x){
	($nRMax, $nMuMax, $nPMax)= split(',', $GridSize);
    }else{
	die "$ERROR -g=$GridSize should be 3 integers separated with commas\n";
    }

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

sub show_settings{

    print "Spatial dimension:     nR =$nRMax\n";
    print "Pitch angle dimension: nMu=$nMuMax\n";
    print "Momentum dimension:    nP =$nPMax\n"

}

#############################################################################

sub print_help{

print "Set the grid size to nR=1500, nMu=60, nP=100:

    Config.pl -g=1500,60,100

Show settings for SP/Kota:

    Config.pl -s
\n\n";
    exit 0;
}
