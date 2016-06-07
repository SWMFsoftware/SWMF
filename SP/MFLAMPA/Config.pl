#!/usr/bin/perl -i
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

# Allow in-place editing
$^I = "";

use strict;

our $Component = "SP";
our $Code = "MFLAMPA";
our @Arguments= @ARGV;

require "../../share/Scripts/Config.pl";

# Variables inherited from share/Scripts/Config.pl
our %Remaining; # Unprocessed arguments
our $ERROR;
our $WARNING;
our $Help;
our $Verbose;
our $Show;
our $ShowGridSize;
our $NewGridSize;
our $NewGhostCell;


&print_help if $Help;

my $Src = 'src';

# Grid size variables
my $NameSizeFile = "$Src/ModSize.f90";
my $GridSize;
my ($iPMin, $iPMax, $nLat, $nLon);
my $nP;

# Read previous grid size
&get_settings;

foreach (@Arguments){
    if(/^-s$/)                {$Show=1;                        next};
    warn "WARNING: Unknown flag $_\n" if $Remaining{$_};
}


# Set new grid size
&set_grid_size if ($NewGridSize and $NewGridSize ne $GridSize);

# Show current settings
my $Settings = &current_settings; print $Settings if $Show;

exit 0;

#############################################################################

sub get_settings{

    # Read size of the grid from $NameSizeFile
    open(FILE, $NameSizeFile) or die "$ERROR could not open $NameSizeFile\n";
    while(<FILE>){
	next if /^\s*!/;
        $iPMin=$1 if /\biParticleMin\s*=\s*(-?\d+)/i;
	$nP   =$1 if /\bnParticle\s*=[^0-9]*(\d+)/i;
	$nLat =$1 if /\bnLat\s*=[^0-9]*(\d+)/i;
	$nLon =$1 if /\bnLon\s*=[^0-9]*(\d+)/i;
    }
    close FILE;

    die "$ERROR could not read iParticleMin from $NameSizeFile\n" 
	unless length($iPMin);

    die "$ERROR could not read nParticle from $NameSizeFile\n" 
	unless length($nP);                         

    die "$ERROR could not read nLat from $NameSizeFile\n" 
	unless length($nLat);                         

    die "$ERROR could not read nLon from $NameSizeFile\n" 
	unless length($nLon);                         

    $iPMax = '' . (-1 + $nP + $iPMin);
    $GridSize = "$nP,$iPMax,$nLat,$nLon";

}

#############################################################################

sub set_grid_size{

    $GridSize = $NewGridSize if $NewGridSize;

    if($GridSize =~ /^[0-9]\d*,[0-9]\d*,[1-9]\d*,[1-9]\d*$/){
	($nP,$iPMax,$nLon,$nLat) = split(',', $GridSize);
    }elsif($GridSize){
	die "$ERROR -g=$GridSize must be 4 integers\n";
    }

    $iPMin = '' . (1 + $iPMax - $nP);

    # Check the grid size (to be set)
    die "$ERROR nParticle=$nP must be positive\n" if $nP<=0;
    die "$ERROR iParticleMin=$iPMin must not be positive\n" if $iPMin>0;
    die "$ERROR iParticleMax=$iPMax must not be negative\n" if $iPMax<0;
    die "$ERROR nLat=$nLat must be positive\n" if $nLat<=0;
    die "$ERROR nLon=$nLon must be positive\n" if $nLon<=0;

    print "Writing new grid size $GridSize into ".
	"$NameSizeFile ...\n";

    @ARGV = ($NameSizeFile);
    while(<>){
	if(/^\s*!/){print; next} # Skip commented out lines
	s/\b(iParticleMin\s*=[^-0-9]*)(-?\d+)/$1$iPMin/i;
	s/\b(nParticle\s*=[^0-9]*)(\d+)/$1$nP/i;
	s/\b(nLat\s*=[^0-9]*)(\d+)/$1$nLat/i;
	s/\b(nLon\s*=[^0-9]*)(\d+)/$1$nLon/i;
	print;
    }
}

#############################################################################

sub print_help{

    print "
Additional options for MFLAMPA/Config.pl:

-g=nP,iPMax,nLon,nLat 
                Set grid size. 
                nP is maximum number of particles per field line,
                iPMax is the maximum allowed particle index,
                nLon, nLat are the grid size at the origin surface.
\n";
    exit 0;
}

#############################################################################

sub current_settings{

    $Settings .= 
	"Number of particles per line   : nParticle=$nP\n";
    $Settings .= 
	"Min. particle index            : iParticleMin=$iPMin\n";
    $Settings .=
	"Size of grid on source surface : nLon=$nLon, nLat=$nLat\n";

}

