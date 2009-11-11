#!/usr/bin/perl -i
use strict;

our $Component       = 'PW';
our $Code            = 'PWOM';
our $MakefileDefOrig = 'src/Makefile.def';
our @Arguments       = @ARGV;

# Planet variables
my $MakefilePlanet = "Makefile.planet";
my $Planet;
my $NewPlanet;

# Make sure that Makefile.planet exists
`touch $MakefilePlanet`;

my $config     = "share/Scripts/Config.pl";
if(-f $config){
    require $config;
}else{
    require "../../$config";
}

# These are inherited from $config
our %Remaining;   # Arguments not handled by share/Scripts/Config.pl
our $Show;
our $Help;
our $ERROR;
our $WARNING;
our $NewGridSize;
our $ShowGridSize;
our $Install;

&print_help if $Help;

$NewPlanet = "Earth" if $Install;

foreach (@Arguments){
    if(/^-saturn$/i)          {$NewPlanet="Saturn";            next};
    if(/^-earth$/i)           {$NewPlanet="Earth";             next};
    if(/^-s$/)                {$Show=1;                        next};

    warn "WARNING: Unknown flag $_\n" if $Remaining{$_};
}

&get_settings;

&set_planet if $NewPlanet and $NewPlanet ne $Planet;

&show_settings if $Show; 

exit 0;

############################################################################
sub get_settings{

    open(FILE, $MakefilePlanet) or 
	die "$ERROR could not open $MakefilePlanet\n";

    while(<FILE>){
        if(/\s*PLANET\s*=\s*(\w+)/){
	    $Planet=$1;
	    last;
	}
    }
    close FILE;

    #die "$ERROR could not find PLANET name in $MakefilePlanet\n" 
    #	unless $Planet;
}

#############################################################################

sub set_planet{
    $Planet = $NewPlanet;

    my $Dir = "src$Planet";
    die "Directory $Dir is missing\n" unless -d $Dir;
    my @file = glob("$Dir/*.f");
    my $file;
    for $file (@file){
	my $outfile = $file;
	$outfile =~ s/^$Dir/src/;
	$outfile =~ s/\.f$/_planet\.f/;
	&shell_command("cp $file $outfile");
    }
    my $Files = "$Dir/PLANET $Dir/Makefile.planet $Dir/ModCommonPlanet.f90".
	" $Dir/upper_heat_conduction.f90";
    $Files .= " $Dir/get_rate.f90" if $Planet eq "Saturn";
    $Files .= " $Dir/ModGlow.f90"  if $Planet eq "Earth";
    &shell_command("cp $Files src/");

    &shell_command("echo PLANET=$Planet > Makefile.planet");
}

############################################################################

sub show_settings{
    print "Planet=$Planet\n";
}

############################################################################

sub print_help{
    print "Additional options for PWOM/Config.pl:

-Earth      Configure PWOM for Earth. This flag is case insensitive.

-Saturn     Configure PWOM for Saturn. This flag is case insensitive.

-s          Show current planet.

Additional examples for PWOM/Config.pl:

Install for Saturn:

    Config.pl -install -Saturn

Reconfigure PWOM for Earth and show selected planet:

    Config.pl -Earth -s
";
    exit 0;
}
