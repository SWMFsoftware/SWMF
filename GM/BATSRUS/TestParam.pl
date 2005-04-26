#!/usr/bin/perl -s

my $Verbose     = $v; undef $v;
my $Help        = $h; undef $h;
my $HelpXmlParam= $H; undef $H;
my $HelpXml     = $X; undef $X;
my $Save        = $s; undef $s;

use strict;

# This script uses the CheckParam.pl script 
my $CheckParamScript  = 'share/Scripts/CheckParam.pl';

# The -H and -X flags are transferred to CheckParam.pl
exec("$CheckParamScript -X") if $HelpXml;
exec("$CheckParamScript -H") if $HelpXmlParam;

if($Help){
    print "
Purpose:

    Check the correctness of the input parameter file using the XML description
    in PARAM.XML. The grid size returned by GridSize.pl, and the PRECISION 
    defined in Makefile.conf are also used during the check, if available.
    This simple script calls the general script in share/Scripts/CheckParam.pl.

Usage:

    TestParam.pl [-h] [-H] [-X] [-v] [-s] [PARAMFILE]

  -h            print help message and stop

  -H            print help about the XML tags used in PARAM.XML files and stop

  -X            print a short introduction to the XML language and stop

  -v            print verbose information

  -s            save PARAM.pl based on PARAM.XML 
                (requires XML-PARSER::EasyTree Perl package)

  PARAMFILE     check parameters in PARAMFILE. Default value is 'run/PARAM.in'

Examples:

  Check the default parameter file run/PARAM.in:

      TestParam.pl

  Check another parameter file:

      TestParam.pl run/test.000/PARAM.expand

  Convert the XML file PARAM.XML into the Perl tree file PARAM.pl:

      TestParam.pl -s

";
    exit 0;
}

my $XmlFile           = 'PARAM.XML';
my $GridSizeScript    = 'GridSize.pl';
my $MakefileConf      = 'Makefile.conf';
my $StandAloneCode    = 'src/stand_alone.f90';
my $NameComp;
my $Precision;
my $GridSize;

if($Save){
    system($CheckParamScript,"-s","-x=$XmlFile");
    exit 0;
}

if(open(CODE,$StandAloneCode)){
    while(<CODE>){
	$NameComp = $1 if /^\s*NameThisComp\s*=\s*[\'\"]([A-Z][A-Z])/;
    }
    close(CODE);
    if(not $NameComp){
       warn "WARNING Could not find NameThisComp setting in $StandAloneCode\n"
	   ."Assuming component=GM\n";
       $NameComp = 'GM';
   }
}else{
    warn "WARNING Could not open $StandAloneCode,".
	" assuming component=GM\n";
    $NameComp = 'GM';
}

if(open(MAKEFILE,$MakefileConf)){
    $Precision = 'unknown';
    while(<MAKEFILE>){
	$Precision = 'double' if 
	    /^\s*PRECISION\s*=\s*(\-r8|\-real_size\s*64)/;
	$Precision = 'single' if 
	    /^\s*PRECISION\s*=\s*(\-r4|\-real_size\s*32)?\s*$/;
    }
    close(MAKEFILE);
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
	       "-c=$NameComp",
	       "-v=$Verbose",
	       "-g=$GridSize",
	       "-p=$Precision",
	       "-x=$XmlFile",
	       @ARGV);

print "@command\n" if $Verbose;
system(@command);
