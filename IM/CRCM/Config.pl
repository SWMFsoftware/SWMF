#!/usr/bin/perl -i
use strict;
our @Arguments       = @ARGV;
our $MakefileDefOrig = "src/Makefile.def";
our $Component = "IM";
our $Code      = "CRCM";
our $Show;
our $Help;
our $ERROR;
our $WARNING;
our %Remaining;   # Arguments not handled by share/Scripts/Config.pl

# Planet variables
my $ConfigLog = "config.log";
my $Planet;
my $NewPlanet;
my $Grid;
my $NewGrid;
# Make sure that config.log exists
`touch $ConfigLog`;

my $config = "share/Scripts/Config.pl";
if(-f $config){
    require $config;
}else{
    require "../../$config";
}

&print_help if $Help;

# Select planet: Current choices are Earth with  H+,O+,He+,e or H+,O+,e, or H+,e
# Select Grid:   Current choices are default and expanded
# More planet configurations will be added in future
foreach (@Arguments){
    if(/^-earthhohe$/i)         {$NewPlanet="EarthHOHe";           next};
    if(/^-earthho$/i)           {$NewPlanet="EarthHO";             next};
    if(/^-earthh$/i)            {$NewPlanet="EarthH";              next};
    if(/^-griddefault$/i)       {$NewGrid="GridDefault";           next};
    if(/^-gridexpanded$/i)      {$NewGrid="GridExpanded";          next};
    if(/^-s$/)                  {$Show=1;                          next};


    warn "WARNING: Unknown flag $_\n" if $Remaining{$_};
}
&get_settings;

&set_planet if $NewPlanet and $NewPlanet ne $Planet;

&set_grid   if $NewGrid and $NewGrid ne $Grid;

&show_settings if $Show; 

exit 0;
############################################################################

sub get_settings{

    #Define defaults
    $Planet='EarthHO';
    $Grid='GridDefault';
    
    open(FILE, $ConfigLog) or 
	die "$ERROR could not open $ConfigLog\n";
    
    while(<FILE>){
        if(/\s*PLANET\s*=\s*(\w+)/){
	    $Planet=$1;
	    last;
	}
    }
    close FILE;
    
    open(FILE, $ConfigLog) or 
	die "$ERROR could not open $ConfigLog\n";
    
    while(<FILE>){
	if(/\s*Grid\s*=\s*(\w+)/){
	    $Grid=$1;
	    last;
	}
    }
    close FILE;

    #die "$ERROR could not find PLANET name in $ConfigLog\n" 
    #	unless $Planet;
}

#############################################################################

sub set_planet{
    my $Files;
    $Planet = $NewPlanet;

    my $Dir = "src";
    die "Directory $Dir is missing\n" unless -d $Dir;

    $Files .= " $Dir/ModEarthHOHe.f90" if $Planet eq "EarthHOHe";
    $Files .= " $Dir/ModEarthHO.f90"   if $Planet eq "EarthHO";
    $Files .= " $Dir/ModEarthH.f90"    if $Planet eq "EarthH";

    &shell_command("cp $Files src/ModPlanet.f90");

    &shell_command("echo PLANET=$NewPlanet > config.log");
    &shell_command("echo Grid=$NewGrid >> config.log");
}
#############################################################################

sub set_grid{
    my $Files;
    $Grid = $NewGrid;

    my $Dir = "src";
    die "Directory $Dir is missing\n" unless -d $Dir;

    $Files .= " $Dir/ModGrid_default.f90"   if $Grid eq "GridDefault";
    $Files .= " $Dir/ModGrid_expanded.f90"  if $Grid eq "GridExpanded";

    &shell_command("cp $Files src/ModGrid.f90");

    &shell_command("echo PLANET=$NewPlanet > config.log");
    &shell_command("echo Grid=$NewGrid >> config.log");
}
############################################################################

sub show_settings{
    print "Planet=$Planet\n";
    print "Grid=$Grid\n";
}

############################################################################

sub print_help{
    print "Additional options for CRCM/Config.pl:

-EarthHOHe       Configure CRCM for Earth with H+,O+,He+,e. 

-EarthHO         Configure CRCM for Earth with H+,O+,e. 

-EarthH          Configure CRCM for Earth with H+,e. 

-GridDefault     Configure CRCM with the default grid

-GridExpanded    Configure CRCM with a grid expanded to higher lat

These flags are case insensitive.

In the future, other ions could be added, or other planets such as Jupiter or 
Saturn could be included.

EXAMPLE: Reconfigure CRCM for EarthH and show selected planet:

    Config.pl -EarthH -s
";
    exit 0;
}
