#!/usr/bin/perl -s
#  Copyright (C) 2002 Regents of the University of Michigan, 
#  portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

# Read command line options
my $Debug       = $D; undef $D;
my $Help        = $h; undef $h;
my $HelpXmlParam= $H; undef $H;
my $HelpXml     = $X; undef $X;
my $nProc       = $n; undef $n;
my $nThread     = ($t or 1); undef $t;
my $Verbose     = $v; undef $v;
my $Format      = $F; undef $F;

use strict;

# Print help message and exit if -h switch was used
&print_help if $Help;

# The script and the XML file names to check the parameters
my $CheckParamScriptOrig  = 'share/Scripts/CheckParam.pl';
my $CheckParamScript = $CheckParamScriptOrig;
$CheckParamScript .= " -D" if $Debug;
$CheckParamScript .= " -F" if $Format;

my $XMLFile  = 'PARAM.XML';
my $ConfigPl = 'Config.pl';

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

# Get the name of the PARAM file 
my $ParamFile = ($ARGV[0] or $ParamFileDefault);
die "$ERROR could not find $ParamFile\n" unless -f $ParamFile;

die "$ERROR could not find $ConfigPl" unless -f $ConfigPl;

# Global variables containing the SWMF settings
my $Precision;
my %Version;
my %GridSize;

&get_settings;

# Read layout info from #COMPONENTMAP/#LAYOUT and check against the SWMF settings
my %Layout;  # Root,Last,Stride for a component ID
my %nProc;   # Number of PE-s for a component ID

&get_layout;
&check_layout if $nProc;

# Check parameters against the XML descriptions
die "$ERROR could not find executable $CheckParamScriptOrig\n"
    unless -x $CheckParamScriptOrig;

# Check CON parameters
my $Registered = join(',',sort keys %Layout);
my $check = 
    "$CheckParamScript -C=$Registered -n=$nProc -p=$Precision $ParamFile";
print "$check\n" if $Verbose;
my $Error = `$check`;

# Extract list of reformatted files from the output: "CheckParam.pl: diff FILE FILE_orig_\n"
my %Formatted;
$Formatted{$1} = 1 while $Error =~ s/CheckParam.pl: diff (\S+) \1_orig_\n//;

# If there is any error left, show it
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
    my $check = 
	"$CheckParamScript -c=$Comp -n=$nProc{$Comp} -p=$Precision ".
	"-x=$xml -g=$GridSize{$Comp} $ParamFile";
    print "$check\n" if $Verbose;
    $Error = `$check`;
    # Hash of reformatted files
    $Formatted{$1} = 1 while $Error =~ s/CheckParam.pl: diff (\S+) \1_orig_\n//;
    # Show error if anything is left
    if($Error){
	$IsError = 1;
	warn "$ERROR parameter errors for $Comp:\n\n$Error\n";
    }
}

foreach my $File (sort keys %Formatted){
    my $Orig = $File."_orig_";
    if(`diff $File $Orig`){
	print "TestParam formatting: diff $File $Orig\n";
    }else{
	unlink $Orig;
    }
}

exit $IsError;

###############################################################################

sub get_settings{

    my $Settings = `./$ConfigPl -show | head -10; ./$ConfigPl -s`;

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

    open(LAYOUT, $ParamFile) or
	die "$ERROR could not open param file $ParamFile\n";

    my $start;
    while(<LAYOUT>){

	if(/^\#(COMPONENTMAP|LAYOUT)/){$start=1; next}

	next unless $start; # Skip lines before #COMPONENTMAP
	last if /^\#END/ or /^$/;   # Ignore lines after #END or empty line

	# extract layout information from one line in the component map
	/^([A-Z][A-Z])\s+(\-?\d+)\s+(\-?\d+)\s+(\-?\d+)\s*(\-?\d*)/ or
	    die "$ERROR incorrect syntax at line $. in $ParamFile:\n$_";

	my $id       = $1;
	my $proc0    = $2;
	my $proclast = $3;
	my $stride   = ($4 or 1);
	my $thread   = ($5 or 1);

	# Evaluate negative values
	$proc0    += $nProc if $proc0 < 0;    $proc0    = 0 if $proc0    < 0;
	$proclast += $nProc if $proclast < 0; $proclast = 0 if $proclast < 0;
	$stride    = $nThread/abs($stride) if $stride < 0;
	$thread    = $nThread/abs($thread) if $thread < 0;
	
	die "$ERROR invalid component ID $id at line $. in $ParamFile:\n$_"
	    unless $Version{$id};

	die "$ERROR $id component has Empty version ".
	    "at line $. in $ParamFile:\n$_"
	    if $Version{$id} eq 'Empty';

	# If root PE exceeds nProc, it will be pushed down. Not an error.
        #die "$ERROR root PE rank=$proc0 should not exceed last PE rank=".
	#    "$proclast\n\tat line $. in $ParamFile:\n$_"
	#    if $proc0 > $proclast;

	$Layout{$1}="$proc0,$proclast,$stride,$thread";
    }
    close(LAYOUT);
    die "$ERROR #COMPONENTMAP was not found in $ParamFile\n" unless $start;

    print "Layout    = ",join('; ',%Layout),"\n" if $Verbose;

}
###############################################################################
sub check_layout{

    my $iProc;
    for $iProc (0..$nProc-1){
	my $Comps; # the components that use PE $iProc
	my $Comp;
	foreach $Comp (sort keys %Layout){
	    my ($RootProc,$LastProc,$StrideProc,$nThread)
		=split(/,/,$Layout{$Comp});

	    if($iProc >= $RootProc and $iProc <=$LastProc 
	       and (($iProc-$RootProc) % $StrideProc)<=$nThread-1){
		$Comps .= "$Comp,";
		$nProc{$Comp}++;
	    }
	}
	die "$ERROR processor $iProc is not used by any component\n".
	    "\tfor the layout in $ParamFile and $nProc processors.\n" 
	    unless $Comps;

    }

    my $Comp;
    foreach $Comp (sort keys %Layout){
	die "$ERROR component $Comp would run on 0 processors\n".
	    "\tfor the layout in $ParamFile and $nProc processors.\n"
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
# The configuration is obtained from Config.pl and it is checked against the 
# component map/layout defined in the PARAM.in file.
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

    Based on the settings shown by Config.pl, check the consistency of 
    settings and the correctness of the input parameter file.

Usage:

  Scripts/TestParam.pl [-h] [-H] [-X] [-v] [-F]
                       [-n=NPROC] [-t=NTHREAD] [PARAMFILE]

  -h            print help message and stop

  -H            print help about the XML tags used in PARAM.XML files and stop

  -X            print a short introduction to the XML language and stop

  -v            print verbose information

  -F            format the PARAMFILE and the included files. 
                The original files are copied into FILENAME_orig_.

  -n=NPROC      assume that SWMF will run on NPROC processors

  -t=NTHREAD    assume that SWMF will run with NTHREAD maximum OpenMP threads

  PARAMFILE     check parameters in PARAMFILE. Default value is 'run/PARAM.in'


Examples:

  Check and format the default parameter file run/PARAM.in:

Scripts/TestParam.pl -F

  Check another parameter file for a 16-processor and 8-thread execution:

Scripts/TestParam.pl -n=16 -t=8 run/test.000/PARAM.expand"
#EOC
    ,"\n\n";
    exit 0;
}
