#!/usr/bin/perl -s

my $Verbose = $v; undef $v;
my $Help    = ($h or $help or $H); undef $h; undef $help; undef $H;

use strict;

if($Help){
    print "
Purpose:

    Check the correctness of the input parameter file using the XML description
    in PARAM.XML. The grid size returned by GridSize.pl, and the PRECISION 
    defined in Makefile.conf are also used during the check, if available.
    This simple script calls the general script in share/Scripts/CheckParam.pl.

Usage:

    TestParam.pl [-h] [-v] [PARAMFILE]

  -h            print help message and stop

  -v            print verbose information

  PARAMFILE     check parameters in PARAMFILE. Default value is 'run/PARAM.in'

";
    exit 0;
}

my $CheckParamScript  = 'share/Scripts/CheckParam.pl';
my $XmlFile           = 'PARAM.XML';
my $GridSizeScript    = 'GridSize.pl';
my $MakefileConf      = 'Makefile.conf';
my $GridSize;
my $Precision;

if(open(MAKEFILE,$MakefileConf)){
    $Precision = 'unknown';
    while(<MAKEFILE>){
	$Precision = 'double' if /^\s*PRECISION\s*=\s*(\-r8|\-real_size\s*64)/;
	$Precision = 'single' if /^\s*PRECISION\s*=\s*(\-r4|\-real_size\s*32)/;
    }
    if($Precision eq 'unknown'){
	warn "WARNING Could not find PRECISION setting in $MakefileConf,".
	    " assuming double precision\n";
	$Precision = 'double';
    }
}else{
    warn "WARNING could not open $MakefileConf, assuming double precision\n";
    $Precision = 'double';
}
close(MAKEFILE);

if(-x $GridSizeScript){
    $GridSize = `$GridSizeScript`; chop($GridSize);
    $GridSize =~ s/$GridSizeScript(\s*-g=)?//;
}else{
    warn "WARNING could not execute $GridSizeScript\n";
}

my @command = (
	       $CheckParamScript,
	       "-S",
	       "-v=$Verbose",
	       "-g=$GridSize",
	       "-p=$Precision",
	       "-x=$XmlFile",
	       @ARGV);

print "@command\n" if $Verbose;
system(@command);
