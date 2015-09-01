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

my $npx;
my $npy;
my $npz;

my $mytopology;
my $mybc;
my @mycpp;

my $periodicX;
my $periodicY;
my $periodicZ;

my $ng;
my $pseudrand; # TRUE
my $h5hut; #FALSE

my $makefilelocal = "src/Makefile.local";
my $VCtopologyfile = "src/include/VCtopology3D.h";


#print "IPIC3D argv = @Arguments \n";

foreach (@Arguments){
    if(/^-s$/)                 {$Show=1;  next};
    if(/^-npx=(.*)/i)          {$npx.="$1";  next};
    if(/^-npy=(.*)/i)          {$npy.="$1";  next};
    if(/^-npz=(.*)/i)          {$npz.="$1";  next};
    if(/^-periodicx=(.*)/i)    {$periodicX.="$1";  next};
    if(/^-periodicy=(.*)/i)    {$periodicY.="$1";  next};
    if(/^-periodicz=(.*)/i)    {$periodicZ.="$1";  next};
    if(/^-ng=(.*)/i)           {$ng="$1";  next};
    if(/^-pseudrand=(.*)/i)    {$pseudrand="$1";  next};
    if(/^-h5hut=(.*)/i)        {$h5hut="$1";  next};
    warn "WARNING: Unknown flag $_\n" if $Remaining{$_};
}


&get_settings;

# Set the number of proc in each direction
&set_topolegy;
# Set periodic boundary or not per direction
&set_boundary;

# set local Makefile options
&set_cpp_options;

&print_help if $Help;
&show_settings if $Show;

exit 0;
#############################################################################

sub get_settings{

    my $mynpx;
    my $mynpy;
    my $mynpz;
    my $myperiodicx;
    my $myperiodicy;
    my $myperiodicz;

    # Read size of the grid from topology file
    open(FILE, $VCtopologyfile) or die "$ERROR could not open $VCtopologyfile\n";
    while(<FILE>){
	next if /^\s*!/; # skip commented out lines
        $mynpx=$1     if /\bXLEN\s*=\s*(\d+)/i;
        $mynpy=$1     if /\bYLEN\s*=\s*(\d+)/i;
        $mynpz=$1     if /\bZLEN\s*=\s*(\d+)/i;
	$myperiodicx=$1 if/\bPERIODICX\s*=\s*(true|false)/i;
	$myperiodicy=$1 if/\bPERIODICY\s*=\s*(true|false)/i;
	$myperiodicz=$1 if/\bPERIODICZ\s*=\s*(true|false)/i;
    }
    close($VCtopologyfile);

    $mytopology = "$mynpx, $mynpy, $mynpz";
    $mybc       = "$myperiodicx, $myperiodicy, $myperiodicz";
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
sub set_topolegy{

    if($mytopology ne "$npx, $npy, $npz"){  

	die "$ERROR File $VCtopologyfile  does not exist!\n" unless -f $VCtopologyfile ;

	@ARGV = ($VCtopologyfile);
	while(<>){
	    # should test for that its a number
	    if($npx ne ""){s/\b(XLEN\s*=[^0-9]*)(\d+)/$1$npx/i; }
	    if($npy ne ""){s/\b(YLEN\s*=[^0-9]*)(\d+)/$1$npy/i; }
	    if($npz ne ""){s/\b(ZLEN\s*=[^0-9]*)(\d+)/$1$npz/i; }
	    print;
	}
    }
}
################################################################################
sub set_boundary{


    if($mybc ne "$periodicX, $periodicY, $periodicZ"){

	die "$ERROR File $VCtopologyfile  does not exist!\n" unless -f $VCtopologyfile ;

	@ARGV = ($VCtopologyfile);
	while(<>){
	    # instead of testing for not empty we should test equal "true" or "false"
	    if($periodicX ne ""){s/\b(PERIODICX\s*=[^.]*;)/PERIODICX = $periodicX;/i;}
	    if($periodicX ne ""){s/\b(PERIODICX_P\s*=[^.]*;)/PERIODICX_P = $periodicX;/i;}
	    if($periodicY ne ""){s/\b(PERIODICY\s*=[^.]*;)/PERIODICY = $periodicY;/i;}
	    if($periodicY ne ""){s/\b(PERIODICY_P\s*=[^.]*;)/PERIODICY_P = $periodicY;/i;}
	    if($periodicZ ne ""){s/\b(PERIODICZ\s*=[^.]*;)/PERIODICZ = $periodicZ;/i;}
	    if($periodicZ ne ""){s/\b(PERIODICZ_P\s*=[^.]*;)/PERIODICZ_P = $periodicZ;/i;}
	    print;
	}
    }
}
################################################################################
sub set_cpp_options{

    die "$ERROR File $makefilelocal does not exist!\n" unless -f $makefilelocal;

    if($h5hut eq "true"){
	print "\nh5hut is not included yeat! \n";
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

    print "Number of subdomains in X, Y and Z directions: $mytopology\n";
    print "Periodicity          in X, Y and Z directions: $mybc\n";
    print "Number of ghost cell layers for the particles: $mycpp[2]\n";
    print "Use pseudo-random numbers (instead of rand()): $mycpp[0]\n";
    print "Use the H5hut library                        : $mycpp[1]\n";
}
################################################################################

sub print_help{

    print "
Additional options for IPIC3D/Config.pl:

-npx=NPX             NPX is the number of subdomains in the X direction
-npy=NPY             NPY is the number of subdomains in the X direction
-npz=NPZ             NPZ is the number of subdomains in the X direction

-periodicx=BOOLEAN   Periodic in X direction: BOOLEAN is true or false
-periodicy=BOOLEAN   Periodic in Y direction: BOOLEAN is true or false
-periodicz=BOOLEAN   Periodic in Z direction: BOOLEAN is true or false

-ng=NG               Use NG ghost cell layers for particles. Default is 3.

-pseudrand=BOOLEAN   Use pseudo-random numbers or the built-in random numbers.
                     Default is true.

-h5hut=BOOLEAN       Use the H5HUT library? Default is false.
-s                   Show IPIC3D settings

Examples:

Set 2 x 2 subdomains in the X and Y directions and periodicity in the Z direction.

   Config.pl -npx=2 -npy=2 -npz=1 -periodicx=false -periodicy=false -periodicz=true

Set number of ghostcells to 1, switch on pseudo-random numbers and H5HUT library

   Config.pl -ng=1 -pseudorand=true -h5hut=true

Show settings

   Config.pl -s

";
    exit -0;
}
