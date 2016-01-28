#!/usr/bin/perl -i
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
use strict;

our $Component       = 'PC';
our $Code            = 'IPIC3D';
our $MakefileDefOrig = 'Makefile.def.ipic3d';
our @Arguments       = @ARGV;

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

my @mycpp;
my $ng;
my $pseudrand; # TRUE
my $h5hut; #FALSE

my $makefilelocal = "src/Makefile.local";

#print "IPIC3D argv = @Arguments \n";

foreach (@Arguments){
    if(/^-s$/)                 {$Show=1;  next};
    if(/^-ng=(.*)/i)           {$ng="$1";  next};
    if(/^-pseudrand=(.*)/i)    {$pseudrand="$1";  next};
    if(/^-h5hut=(.*)/i)        {$h5hut="$1";  next};
    warn "WARNING: Unknown flag $_\n" if $Remaining{$_};
}

&get_settings;

# set local Makefile options
&set_cpp_options;

&print_help if $Help;
&show_settings if $Show;

exit 0;
#############################################################################

sub get_settings{
    my $mypseudorand = "false";
    my $myh5hut      = "false";
    my $myngp;

    # Read size of the grid from topology file
    open(FILE, $makefilelocal) or die "$ERROR could not open $makefilelocal\n";
    while(<FILE>){
	next if /^\s*!/; # skip commented out lines
	$mypseudorand="true" if/^.*-DPSEUDRAND(\D+)/i;
	$myh5hut="true" if/^.*-DUSEH5HUT(\D+)/i;
	$myngp=$1 if/^.*-DNG_P=(\d+)/i;
    }
    close $makefilelocal;

    @mycpp =  ($mypseudorand, $myh5hut, $myngp);

}

################################################################################
sub set_cpp_options{

    die "$ERROR File $makefilelocal does not exist!\n" unless -f $makefilelocal;

    if($h5hut eq "true"){
	print "\nh5hut is not included yet! \n";
	print "http://www-vis.lbl.gov/Research/H5hut/ \n";
	$h5hut = "false";
    }

    # Check if it already is set
    open(FILE,$makefilelocal);
    while (<FILE>) {
	if($_ =~ m/^.*-DPSEUDRAND.*/i and $pseudrand eq "true"){
	    $pseudrand = "";
	}
	if($_ =~ m/^.*-DUSEH5HUT.*/i and $h5hut eq "true"){
	    $h5hut = "";
	}
    }
    close FILE;

    die "$ERROR File $makefilelocal does not exist!\n" unless -f $makefilelocal;

    @ARGV = ($makefilelocal);
    while(<>){
	if($ng ne "") {s/DNG_P=[0-9*]/DNG_P=$ng/;};
	if($pseudrand eq "false") {s/\s*-DPSEUDRAND \\\n//;}; # always try to remove 
	if($h5hut     eq "false") {s/\s*-DUSEH5HUT \\\n//;}; # always try to remove 
	if($pseudrand eq "true") {s/FLAGC_LOCAL=/FLAGC_LOCAL=-DPSEUDRAND \\\n\t/;};
	if($h5hut     eq "true") {s/FLAGC_LOCAL=/FLAGC_LOCAL=-DUSEH5HUT \\\n\t/;};
	print;
    }
}

################################################################################
sub show_settings{
    print "Number of ghost cell layers for the particles: $mycpp[2]\n";
    print "Use pseudo-random numbers (instead of rand()): $mycpp[0]\n";
    print "Use the H5hut library                        : $mycpp[1]\n";
}
################################################################################

sub print_help{

    print "
Additional options for IPIC3D/Config.pl:

-ng=NG               Use NG ghost cell layers for particles. Default is 3.

-pseudrand=BOOLEAN   Use pseudo-random numbers or the built-in random numbers.
                     Default is true.

-h5hut=BOOLEAN       Use the H5HUT library? Default is false.
-s                   Show IPIC3D settings

Examples:
Set number of ghostcells to 1, switch on pseudo-random numbers and H5HUT library

   Config.pl -ng=1 -pseudorand=true -h5hut=true

Show settings

   Config.pl -s

";
    exit -0;
}
