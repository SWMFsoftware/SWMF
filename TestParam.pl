#!/usr/bin/perl -s

# Read command line options
my $Debug       = $D; undef $D;
my $Help        = $h; undef $h;
my $LayoutFile  = $l; undef $l;
my $nProc       = $n; undef $n;
my $Verbose     = $v; undef $v;

use strict;

# Print help message and exit if -h switch was used
&print_help if $Help;

# Error and warning strings
my $ERROR = 'TestParam_ERROR:';
my $WARNING='TestParam_WARNING:';

# Set default values
my $ParamFileDefault  = 'run/PARAM.in';
my $LayoutFileDefault = 'run/LAYOUT.in';

# The script and the XML file names to check the parameters
my $CheckParamScript  = 'Common/CheckParam.pl';
my $XMLFile           = 'PARAM.XML';

my $SetSWMF = 'SetSWMF.pl';

# Overwrite default is necessary
my $ParamFile = ($ARGV[0] or $ParamFileDefault);
$LayoutFile = $LayoutFileDefault unless $LayoutFile;

die "$ERROR could not find $SetSWMF" unless -f $SetSWMF;

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
}else{
    print "$WARNING No layout file $LayoutFile was found\n";
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
die "$ERROR parameter errors for CON in $ParamFile:\n\n$Error\n" if $Error;

# Check component parameters
my $Comp;
foreach $Comp (sort keys %Layout){
    my $xml = "$Comp/$Version{$Comp}/$XMLFile";
    if(not -f $xml){
	print "$WARNING could not find XML description $xml\n";
	next;
    }
    my $check = "$CheckParamScript -c=$Comp -n=$nProc{$Comp} -p=$Precision ".
	"-x=$xml -g=$GridSize{$Comp} $ParamFile";
    print "$check\n" if $Verbose;
    $Error = `$check`;
    die "$ERROR parameter errors for $Comp:\n\n$Error\n" if $Error;
}

exit 0;

###############################################################################

sub get_settings{

    my $Settings = `./$SetSWMF -s`;

    die "$ERROR $SetSWMF -s did not provide the settings\n" unless $Settings;

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
    die "$ERROR code versions could not be found in $Settings\n" 
	unless %Version;

    print 
       "Installed = $Installed\n".
       "Precision = $Precision\n".
       "Versions  = ",join('; ',%Version),"\n".
       "GridSize  = ",join('; ',%GridSize),"\n" if $Verbose;

    die "$ERROR IH/UofM_share requires GM/UofM and not GM/$Version{GM}\n"
	if $Version{IH} eq 'UofM_share' and $Version{GM} ne 'UofM';

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
	if($iProc==10000){
	    print "Hi $Comps=$Comps\n";
	}
	die "$ERROR processor $iProc is not used by any component\n".
	    "\tfor the layout in $LayoutFile and $nProc processors.\n" 
	    unless $Comps;

	die "$ERROR IH/UofM_share overlaps with GM\n".
	    "\tfor the layout in $LayoutFile and $nProc processors.\n"
	    if $Version{IH} eq "UofM_share" and $Comps =~ /GM,.*IH,/;
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

sub print_help{

    print "
Purpose:

    Based on the settings shown be SetSWMF.pl and the registration and layout 
    information contained in the layout file, check the consistency of 
    settings and the correctness of the input parameter file.

Usage:

    TestParam.pl [-h] [-v] [-l=LAYOUTFILE] [-n=NPROC] [PARAMFILE]

  -h            print help message and stop

  -v            print verbose information

  -l=LAYOUTFILE obtain layout from LAYOUTFILE. Default value is 'run/LAYOUT.in'

  -n=NPROC      assume that SWMF will run on NPROC processors

  PARAMFILE     check parameters in PARAMFILE. Default value is 'run/PARAM.in'

";
    exit 0;
}



