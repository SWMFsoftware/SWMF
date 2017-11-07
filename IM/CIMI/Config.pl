#!/usr/bin/perl -i
use strict;
our @Arguments       = @ARGV;
our $MakefileDefOrig = "src/Makefile.def";
our $Component = "IM";
our $Code      = "CIMI";
our $Show;
our $Help;
our $ERROR;
our $WARNING;
our $Install;
our $Compiler;

our %Remaining;   # Arguments not handled by share/Scripts/Config.pl

# Planet variables
my $ConfigLog = "config.log";
my $Planet;
my $NewPlanet;
my $Grid;
my $NewGrid;

$NewPlanet="EarthHO" if $Install;

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
    if(/^-compiler=(.*)$/i)     {$Compiler=$1;                     next};        
    warn "WARNING: Unknown flag $_\n" if $Remaining{$_};
}

&get_settings;

&set_planet if $NewPlanet and $NewPlanet ne $Planet;

&set_grid   if $NewGrid and $NewGrid ne $Grid;

&set_compiler  if $Install;

&show_settings if $Show; 

exit 0;
############################################################################

sub get_settings{

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
    my $File;
    $Planet = $NewPlanet;

    my $Dir = "src";
    die "$ERROR Directory $Dir is missing\n" unless -d $Dir;

    $File = " $Dir/ModEarthHOHe.f90" if $Planet eq "EarthHOHe";
    $File = " $Dir/ModEarthHO.f90"   if $Planet eq "EarthHO";
    $File = " $Dir/ModEarthH.f90"    if $Planet eq "EarthH";

    &shell_command("cp $File src/ModPlanet.f90");

    &shell_command("echo PLANET=$NewPlanet > config.log");
    &shell_command("echo Grid=$NewGrid >> config.log");
}
#############################################################################

sub set_compiler{

    my $Dir = "srcSAMI3";
    die "$ERROR Directory $Dir is missing\n" unless -d $Dir;
    my $File = "$Dir/Makefile.suff";

    if(-e "$File\.$Compiler"){
	&shell_command("cp $File\.$Compiler $File");
    }else{
	&shell_command("rm -f $File; touch $File");
    }
}
#############################################################################


sub set_grid{
    my $File;
    $Grid = $NewGrid;

    my $Dir = "src";
    die "$ERROR Directory $Dir is missing\n" unless -d $Dir;

    $File = " $Dir/ModGrid_default.f90"   if $Grid eq "GridDefault";
    $File = " $Dir/ModGrid_expanded.f90"  if $Grid eq "GridExpanded";

    &shell_command("cp $File src/ModGrid.f90");

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
    print "Additional options for CIMI/Config.pl:

-EarthHOHe       Configure CIMI for Earth with H+,O+,He+,e. 

-EarthHO         Configure CIMI for Earth with H+,O+,e. 

-EarthH          Configure CIMI for Earth with H+,e. 

-GridDefault     Configure CIMI with the default grid

-GridExpanded    Configure CIMI with a grid expanded to higher lat

These flags are case insensitive.

In the future, other ions could be added, or other planets such as Jupiter or 
Saturn could be included.

EXAMPLE: Reconfigure CIMI for EarthH and show selected planet:

    Config.pl -EarthH -s
";
    exit 0;
}
