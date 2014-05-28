#!/usr/bin/perl -s
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

# Read command line options
my $Debug       = $D; undef $D;
my $Help        = $h; undef $h;
my $HelpXmlParam= $H; undef $H;
my $HelpXml     = $X; undef $X;
my $LayoutFile  = $l; undef $l;
my $nProc       = $n; undef $n;
my $Verbose     = $v; undef $v;

use strict;

# Print help message and exit if -h switch was used
&print_help if $Help;

# The script and the XML file names to check the parameters
my $CheckParamScript  = 'share/Scripts/CheckParam.pl';
my $XMLFile           = 'PARAM.XML';
my $ConfigPl           = 'Config.pl';

# The -H, -X and -s flags are transferred to CheckParam.pl
exec("$CheckParamScript -X") if $HelpXml;
exec("$CheckParamScript -H") if $HelpXmlParam;

# Error and warning strings
my $ERROR = 'TestParam_ERROR:';
my $WARNING='TestParam_WARNING:';

# Error flag to be returned to Unix
my $IsError;

# Set default values
my $ParamFileDefault  = 'run/PARAM.in';
my $LayoutFileDefault = 'run/LAYOUT.in';

# Get the name of the PARAM file 
my $ParamFile = ($ARGV[0] or $ParamFileDefault);
die "$ERROR could not find $ParamFile\n" unless -f $ParamFile;

# Try to guess the name of the LAYOUT file if not given
if(not $LayoutFile){
    if($ParamFile =~ /PARAM/){
	$LayoutFile = $ParamFile;
	# Replace PARAM with LAYOUT
	$LayoutFile =~ s/PARAM/LAYOUT/;
	# Replace .expand with .in
	$LayoutFile =~ s/\.expand/.in/;
	# Delete .bak
	$LayoutFile =~ s/\.bak//;
	# If the LAYOUT file based on the PARAM file name is not found 
	# try the default
	$LayoutFile = $LayoutFileDefault 
	    if (not -f $LayoutFile and -f $LayoutFileDefault);
    }else{
	$LayoutFile = $LayoutFileDefault;
    }
}
warn "$WARNING No layout file $LayoutFile was found\n" unless -f $LayoutFile;

die "$ERROR could not find $ConfigPl" unless -f $ConfigPl;

# Global variables containing the SWMF settings
my $Precision;
my %Version;
my %GridSize;

&get_settings;

# Read layout info from LAYOUT file and check against the SWMF settings
my %Layout;  # Root,Last,Stride for a component ID
my %nProc;   # Number of PE-s for a component ID
if(-f $LayoutFile){
    &get_layout;
    &check_layout if $nProc;
}

# Check parameters against the XML descriptions
die "$ERROR could not find executable $CheckParamScript\n"
    unless -x $CheckParamScript;

# Check CON parameters
my $Registered = join(',',sort keys %Layout);
my $check = 
    "$CheckParamScript -C=$Registered -n=$nProc -p=$Precision $ParamFile";
print "$check\n" if $Verbose;
my $Error = `$check`;
if($Error){
    $IsError = 1;
    warn "$ERROR parameter errors for CON in $ParamFile:\n\n$Error\n";
}

# Check component parameters
my $Comp;
foreach $Comp (sort keys %Layout){
    my $xml = "$Comp/$Version{$Comp}/$XMLFile";
    if(not -f $xml){
	warn "$WARNING could not find XML description $xml\n";
	next;
    }
    my $check = "$CheckParamScript -c=$Comp -n=$nProc{$Comp} -p=$Precision ".
	"-x=$xml -g=$GridSize{$Comp} $ParamFile";
    print "$check\n" if $Verbose;
    $Error = `$check`;
    if($Error){
	$IsError = 1;
	warn "$ERROR parameter errors for $Comp:\n\n$Error\n";
    }
}

exit $IsError;

###############################################################################

sub get_settings{

    my $Settings = `./$ConfigPl -show | head -9; ./$ConfigPl -s`;

    die "$ERROR $ConfigPl -s did not provide the settings\n" unless $Settings;

    my $Installed;
    foreach (split(/\n/,$Settings)){
	$Installed = 1  if /is installed/;
	$Precision = $1 if /(single|double) precision/;
	if(/^\t([A-Z][A-Z])\/(\w+)/){
	    my $Comp=$1; my $Version=$2;
	    $Version{$Comp}  = $Version;
	    $GridSize{$Comp} = $1 if /grid:\s*(.*)/i;
	}
    }

    die "$ERROR SWMF is not installed based on $Settings" unless $Installed;
    die "$ERROR precision could not be found in $Settings\n" unless $Precision;
    die "$ERROR component versions could not be found in $Settings\n" 
	unless %Version;

    print 
       "Installed = $Installed\n".
       "Precision = $Precision\n".
       "Versions  = ",join('; ',%Version),"\n".
       "GridSize  = ",join('; ',%GridSize),"\n" if $Verbose;

}

###############################################################################
sub get_layout{

    open(LAYOUT,$LayoutFile) or 
	die "$ERROR could not open layout file $LayoutFile\n";

    my $start;
    while(<LAYOUT>){

	if(/^\#COMPONENTMAP/){$start=1; next}

	next unless $start; # Skip lines before #COMPONENTMAP
	last if /^\#END/;   # Ignore lines after #END

	# extract layout information from one line in the component map
	/^([A-Z][A-Z])\s+(\d+)\s+(\d+)\s+(\d+)/ or
	    die "$ERROR incorrect syntax at line $. in $LayoutFile:\n$_";

	die "$ERROR invalid component ID $1 at line $. in $LayoutFile:\n$_"
	    unless $Version{$1};

	die "$ERROR $1 component has Empty version ".
	    "at line $. in $LayoutFile:\n$_"
	    if $Version{$1} eq 'Empty';

	die "$ERROR root PE rank=$3 should not exceed last PE rank=$2\n".
	    "\tat line $. in $LayoutFile:\n$_"
	    if $2 > $3;

	die "$ERROR stride=$4 must be positive at line $. in $LayoutFile:\n$_"
	    if $4 < 1;

	$Layout{$1}="$2,$3,$4";
	
    }
    close(LAYOUT);
    die "$ERROR #COMPONENTMAP was not found in $LayoutFile\n" unless $start;

    print "Layout    = ",join('; ',%Layout),"\n" if $Verbose;

}
###############################################################################
sub check_layout{

    my $iProc;
    for $iProc (0..$nProc-1){
	my $Comps; # the components that use PE $iProc
	my $Comp;
	foreach $Comp (sort keys %Layout){
	    my ($RootProc,$LastProc,$StrideProc)=split(/,/,$Layout{$Comp});

	    if($iProc >= $RootProc and $iProc <=$LastProc 
	       and (($iProc-$RootProc) % $StrideProc)==0){
		$Comps .= "$Comp,";
		$nProc{$Comp}++;
	    }
	}
	die "$ERROR processor $iProc is not used by any component\n".
	    "\tfor the layout in $LayoutFile and $nProc processors.\n" 
	    unless $Comps;

    }

    my $Comp;
    foreach $Comp (sort keys %Layout){
	die "$ERROR component $Comp would run on 0 processors\n".
	    "\tfor the layout in $LayoutFile and $nProc processors.\n"
	    unless $nProc{$Comp}
    }

    print "nProc     = ",join('; ',%nProc),"\n" if $Verbose;
}
###############################################################################
#!QUOTE: \clearpage
#BOP
#!QUOTE: \subsection{Checking Layout and Parameters with Scripts/TestParam.pl}
#!ROUTINE: TestParam.pl - check layout and parameter file before running SWMF
#!DESCRIPTION:
# A lot of time and frustration can be saved if the consistency of the
# SWMF configuration with the layout and input parameters are checked before
# running the code. The TestParam.pl script does extensive checks and provides
# detailed warnings and error messages. 
# The configuration is obtained from Config.pl and it is checked against the LAYOUT file.
# The parameters of CON and the physics components are checked against the XML
# descriptions in the PARAM.XML files using the share/Scripts/CheckParam.pl script.
#
# If there are no errors, no output is produced.
#
#!REVISION HISTORY:
# 09/01/2003 G. Toth - initial version
#                      many improvements and extensions
#EOP

sub print_help{

    print 
#BOC
"Purpose:

    Based on the settings shown by Config.pl and the registration and layout 
    information contained in the layout file, check the consistency of 
    settings and the correctness of the input parameter file.

Usage:

  Scripts/TestParam.pl [-h] [-H] [-X] [-v]
                       [-l=LAYOUTFILE] [-n=NPROC] [PARAMFILE]

  -h            print help message and stop

  -H            print help about the XML tags used in PARAM.XML files and stop

  -X            print a short introduction to the XML language and stop

  -v            print verbose information

  -l=LAYOUTFILE obtain layout from LAYOUTFILE. 
                Default name of the LAYTOUFILE is obtained from the name
                of the PARAMFILE: 
                'PARAM' is replaced with 'LAYOUT' and '.expand' with '.in' 
                and if that file is not found, run/LAYOUT.in is used 
                if it exists

  -n=NPROC      assume that SWMF will run on NPROC processors

  PARAMFILE     check parameters in PARAMFILE. Default value is 'run/PARAM.in'


Examples:

  Check the default parameter file run/PARAM.in without checking the layout:

Scripts/TestParam.pl

  Check another parameter and layout file for a 16 processor execution:

Scripts/TestParam.pl -n=16 run/test.000/PARAM.expand"
#EOC
    ,"\n\n";
    exit 0;
}
