#!/usr/bin/perl -s
my $Planet = ($p or $planet);
my $Help   = ($h or $help);
use strict;
die "Usage: config.pl -p=PLANETNAME\n" if $Help or not $Planet;

my $Dir = "src$Planet";
die "Directory $Dir is missing\n" unless -d $Dir;
my @file = glob("$Dir/*.f");
my $file;
for $file (@file){
    my $outfile = $file;
    $outfile =~ s/^$Dir/src/;
    $outfile =~ s/\.f$/_planet\.f/;
    `cp $file $outfile`;
}
`cp $Dir/PLANET $Dir/Makefile.planet src/`;
`cp $Dir/PLANET $Dir/ModCommonPlanet.f90 src/`;

`echo "PLANET=$Planet" > Makefile.planet`;
