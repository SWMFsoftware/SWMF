#!/usr/bin/perl -i
use strict;

our @Arguments = ('srcMake/Makefile.def', 'UA/GITM2', @ARGV);
our %Remaining;   # Arguments not handled by share/Scripts/config.pl
my $config     = "share/Scripts/config.pl";

if(-f $config){
    require $config;
}else{
    require "../../$config";
}

my $ModPlanet = "ModPlanet.f90";

my $Planet;
my $PlanetOrig;
my $Show;
my $Help;

$Show = 1 if not @ARGV;

&get_settings;

foreach (@ARGV){
    if(/^-install/)           {$Planet="Earth" unless $Planet; next};
    if(/^-mars$/i)            {$Planet="Mars";                 next};
    if(/^-earth$/i)           {$Planet="Earth";                next};
    if(/^-s(how)?$/)          {$Show  =1;                      next};
    if(/^-h(elp)?$/i)         {&print_help;                    exit};

    warn "WARNING: Unknown flag $_\n" if $Remaining{$_};
}


if ($Planet ne $PlanetOrig) { 
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

if($Show){
    print "GITM2 is configured for planet=$PlanetOrig\n";
}

exit 0;

############################################################################
sub get_settings{

    if(-l "src/$ModPlanet"){
        my $link = `ls -l src/$ModPlanet`;
	$link =~ /Mod(\w+)\.f90$/ or 
	    warn "GITM2/config: Could not find planet in $link";
        $PlanetOrig = $1;
    }

    print "PlanetOrig = $PlanetOrig\n";
}

############################################################################

sub print_help{
    print "Additional switches for GITM2/config.pl:

   -p=PLANET   Configure for PLANET. Possible values are 'Mars' or 'Earth'
               Default planet is Earth.

   -Mars       Configure GITM2 for Mars. This flag is case insensitive.

   -Earth      Configure GITM2 for Earth. This flag is case insensitive.

Examples:

Install for Mars:

   config.pl -install -planet=Mars

Configure GITM for Earth:

   config.pl -Earth

 ";
    exit 0;
}
