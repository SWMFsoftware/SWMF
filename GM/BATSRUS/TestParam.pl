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
    in PARAM.XML. The grid size returned by Config.pl -g, and the PRECISION 
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
my $ConfigScript      = 'Config.pl';
my $NameComp='GM';
my $Precision;
my $GridSize;

if($Save){
    system($CheckParamScript,"-s","-x=$XmlFile");
    exit 0;
}

if(-x $ConfigScript){
    $GridSize = `$ConfigScript -g`; chop($GridSize);
    $GridSize =~ s/$ConfigScript(\s*-g=)?//;

    $Precision = `Config.pl -show`;
    $Precision = $1 if $Precision =~ /\b(single|double)\b/i;

}else{
    warn "WARNING could not execute $ConfigScript\n";
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
