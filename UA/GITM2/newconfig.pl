#!/usr/bin/perl -i
use strict;

our $Component       = 'UA';
our $Code            = 'GITM2';
our $MakefileDefOrig = 'srcMake/Makefile.def';
our @Arguments       = @ARGV;

my $config     = "share/Scripts/config.pl";
if(-f $config){
    require $config;
}else{
    require "../../$config";
}

# These are inherited from $config
our %Remaining;   # Arguments not handled by share/Scripts/config.pl
our $Show;
our $Help;
our $ERROR;
our $WARNING;
our $NewGridSize;
our $ShowGridSize;

&print_help if $Help;

# Grid size variables
my $NameGridFile="src/ModSize.f90"; # File that contains the grid size
my $GridSize;
my ($nLon, $nLat, $nAlt, $MaxBlock);

# Planet variables
my $ModPlanet = "ModPlanet.f90";
my $Planet;
my $PlanetOrig;

&get_settings;

foreach (@Arguments){
    if(/^-install/)           {$Planet="Earth" unless $Planet; next};
    if(/^-mars$/i)            {$Planet="Mars";                 next};
    if(/^-earth$/i)           {$Planet="Earth";                next};
    if(/^-s$/)                {$Show=1;                        next};

    warn "WARNING: Unknown flag $_\n" if $Remaining{$_};
}


&set_grid_size if $NewGridSize and $NewGridSize ne $GridSize;

&set_planet if $Planet and $Planet ne $PlanetOrig;

&show_settings if $Show; 

print "config.pl -g=$GridSize\n" if $ShowGridSize and not $Show;

exit 0;

############################################################################
sub get_settings{

    # Read size of the grid from $NameGridFile
    open(FILE,$NameGridFile) or die "$ERROR could not open $NameGridFile\n";
    while(<FILE>){
	next if /^\s*!/; # skip commented out lines
        $nLon=$1         if /\bnLons\s*=\s*(\d+)/i;
	$nLat=$1         if /\bnLats\s*=\s*(\d+)/i;
	$nAlt=$1         if /\bnAlts\s*=\s*(\d+)/i;
	$MaxBlock=$1     if /\bnBlocksMax\s*=\s*(\d+)/i;
	last if $nLon and $nLat and $nAlt and $MaxBlock;
    }
    close FILE;

    die "$ERROR could not read nLon from $NameGridFile\n" unless length($nLon);
    die "$ERROR could not read nLat from $NameGridFile\n" unless length($nLat);
    die "$ERROR could not read nAlt from $NameGridFile\n" unless length($nAlt);
    die "$ERROR could not read MaxBlock from $NameGridFile\n" 
	unless length($MaxBlock);

    $GridSize = "$nLon,$nLat,$nAlt,$MaxBlock";

    # Get the current planet
    if(-l "src/$ModPlanet"){
        my $link = `ls -l src/$ModPlanet`;
	$link =~ /Mod(\w+)\.f90$/ or 
	    warn "GITM2/config: Could not find planet in $link";
        $PlanetOrig = $1;
    }
}

#############################################################################

sub set_grid_size{

    print "Writing new grid size $NewGridSize into $NameGridFile...\n";

    $GridSize = $NewGridSize;

    if($GridSize=~/^\d+,\d+,\d+,\d+$/){
	($nLon, $nLat, $nAlt, $MaxBlock)= split(',', $GridSize);
    }elsif($GridSize){
	die "$ERROR -g=$GridSize should be 4 integers separated with commas\n";
    }

    if($MaxBlock and $MaxBlock !~ /^\d+$/){
	die "$ERROR -b=$MaxBlock should be a positive integer\n";
    }

    @ARGV = ($NameGridFile);
    while(<>){
	if(/^\s*!/){print; next} # Skip commented out lines
	s/\b(nLons\s*=[^\d]*)(\d+)/$1$nLon/i;
	s/\b(nLats\s*=[^\d]*)(\d+)/$1$nLat/i;
	s/\b(nAlts\s*=[^\d]*)(\d+)/$1$nAlt/i;
	s/\b(nBlocksMax\s*=[^\d]*)(\d+)/$1$MaxBlock/i;
	print;
    }

}

############################################################################

sub set_planet{
    $PlanetOrig = $Planet;

    chdir "src" or die "Could not change directory to src\n";

    print "Configuring GITM for $Planet!!\n"; 

    &shell_command("rm -f ModPlanet.f90 ModPlanet.o planet.f90");
    &shell_command("ln -s Mod$Planet.f90 ModPlanet.f90");
    &shell_command("ln -s $Planet.f90 planet.f90");

    my $file;
    foreach $file (glob("*.$Planet.f90")) {
	$file =~ /(.*)\.$Planet\.f90/;
	my $name = $1;
	unlink "$name.f90";
	&shell_command("ln -s $name.$Planet.f90 $name.f90");
    }

    chdir "..";

    &shell_command("cd srcData ; cp UAM.in.$Planet UAM.in");
}

############################################################################

sub show_settings{

    print "Number of cells in a block: nLon=$nLon, nLat=$nLat, nAlt=$nAlt\n";
    print "Max. number of blocks     : MaxBlock=$MaxBlock\n";
    print "Planet=$PlanetOrig\n";
}

############################################################################

sub print_help{
    print "Additional switches for GITM2/config.pl:

-Mars       Configure GITM2 for Mars. This flag is case insensitive.

-Earth      Configure GITM2 for Earth. This flag is case insensitive.

-s          Show current planet.

Additional examples for GITM2/config.pl:

Install for Mars:

    config.pl -install -Mars

Reconfigure GITM for Earth:

    config.pl -Earth

Set grid to nLon=36, nLat=36, nAlt=50 and the number of blocks to 16:

    config.pl -g=36,36,50,16

 ";
    exit 0;
}
