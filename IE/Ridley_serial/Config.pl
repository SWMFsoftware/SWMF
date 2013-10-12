#!/usr/bin/perl -i
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
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

foreach (@Arguments){
    if(/^-s$/)                {$Show=1;  next};
    warn "WARNING: Unknown flag $_\n" if $Remaining{$_};
}

# Grid size variables
my $NameFile="src/ModSize.f90";
my $GridSize;
my ($nTheta, $nPsi);

&get_settings;

&set_grid_size if $NewGridSize and $NewGridSize ne $GridSize;

&show_settings if $Show;

print "Config.pl -g=$nTheta,$nPsi\n" if $ShowGridSize and not $Show;

exit 0;

#############################################################################

sub get_settings{

    # Read size of the grid from $NameFile
    open(FILE,$NameFile) or die "$ERROR could not open $NameFile\n";
    while(<FILE>){
	next if /^\s*!/; # skip commented out lines
	$nTheta=$1 if /\bIONO_nTheta\s*=\s*(\d+)/i;
        $nPsi  =$1 if /\bIONO_nPsi\s*=\s*(\d+)/i;
	last if $nTheta and $nPsi;
    }
    close FILE;

    die "$ERROR could not read nTheta from $NameFile\n" 
	unless length($nTheta);
    die "$ERROR could not read nPsi from $NameFile\n" 
	unless length($nPsi);

    $GridSize = "$nTheta,$nPsi";
}

#############################################################################

sub set_grid_size{

    $GridSize = $NewGridSize;

    if($GridSize=~/^\d+,\d+$/){
	($nTheta, $nPsi)= split(',', $GridSize);
    }else{
	die "$ERROR -g=$GridSize should be 2 integers separated with a comma\n"
    }

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

sub show_settings{

    print "Number of colatitudes/hemisphere  : nTheta = $nTheta\n";
    print "Number of longitudes              : nPsi   = $nPsi\n";
}

#############################################################################

sub print_help{

print "Set grid size to nTheta=91, nPsi=361 (1 deg by 1 deg):

    Config.pl -g=91,361

Show settings specific to IE/Ridley_serial:

    Config.pl -s

";
    exit 0;
}
